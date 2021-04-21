#include "RooUnfoldBlobel.h"

#include <iostream>
#include <iomanip>
// #include <Eigen/Core>
// #include <Eigen/Eigenvalues>
#include "TClass.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TDecompSVD.h"

#include "RooUnfoldHelpers.h"
#include "RooUnfoldResponse.h"
// #include "RooUnfoldTH1Helpers.h"
// #include "RooUnfoldFitHelpers.h"

using namespace RooUnfolding;

using std::cout;
using std::cerr;
using std::endl;

ClassImp(RooUnfoldBlobel)

RooUnfoldBlobel::RooUnfoldBlobel (const RooUnfoldBlobel& rhs)
  : RooUnfold (rhs)
{
  //duplicate constructor
  Init();
  CopyData (rhs);
}

RooUnfoldBlobel::RooUnfoldBlobel (const RooUnfoldResponse* res, const TH1* meas, Int_t kreg,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  //! Constructor with response matrix object and measured unfolding input histogram.
  //! The regularisation parameter is kreg.
  Init();
}

RooUnfoldBlobel::RooUnfoldBlobel (const RooUnfoldResponse* res, const TH1* meas, Int_t kreg, Int_t ntoyssvd,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  //! Constructor with old ntoys argument. No longer required, left over from svd, kept for posterity
  Init();
  _NToys = ntoyssvd;
}

RooUnfoldBlobel*
RooUnfoldBlobel::Clone (const char* newname) const
{
  RooUnfoldBlobel* unfold= new RooUnfoldBlobel(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

void
RooUnfoldBlobel::Reset()
{
  // destroy and re-initialize this object
  Destroy();
  Init();
  RooUnfold::Reset();
}

void
RooUnfoldBlobel::Destroy()
{
  //! delete all members of this object
  delete _meas1d;
  delete _train1d;
  delete _truth1d;
  delete _reshist;
}

void
RooUnfoldBlobel::Init()
{
  //! initialize this object with zero values
  _meas1d= _train1d= _truth1d= 0;
  _reshist= _meascov= 0;
  Hessian = 0;
  GetSettings();
}

void
RooUnfoldBlobel::Assign (const RooUnfoldBlobel& rhs)
{
  //! assign data from another instance
  RooUnfold::Assign (rhs);
  CopyData (rhs);
}

void
RooUnfoldBlobel::CopyData (const RooUnfoldBlobel& rhs)
{
  //! copy data from another instance
  _kreg= rhs._kreg;
}

void
RooUnfoldBlobel::Unfold()
{
  //! perform the unfolding
  if (_res->GetDimensionTruth() != 1 || _res->GetDimensionMeasured() != 1) {
    std::cerr << "RooUnfoldBlobel may not work very well for multi-dimensional distributions" << std::endl;
  }
  if (_res->GetDimensionTruth() != _res->GetDimensionMeasured()){
    std::cerr << "RooUnfoldBlobel may not work very well for distributions of different sizes" << std::endl;
  }
  if (_kreg < 0) {
    std::cerr << "RooUnfoldBlobel invalid kreg: " << _kreg << std::endl;
    return;
  }
  _nb= _nm > _nt ? _nm : _nt;

  if (_kreg > _nb) {
    cerr << "RooUnfoldSvd invalid kreg=" << _kreg << " with " << _nb << " bins" << endl;
    return;
  }
  // this->PrepareHistograms();

  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _meas1d=  HistNoOverflow (_meas,             _overflow);
  _train1d= HistNoOverflow (_res->Hmeasured(), _overflow);
  _truth1d= HistNoOverflow (_res->Htruth(),    _overflow);
  _reshist= _res->HresponseNoOverflow();
  Resize (_meas1d,  _nb);
  Resize (_train1d, _nb);
  Resize (_truth1d, _nb);
  Resize (_reshist, _nb, _nb);
  if (_res->FakeEntries()) {
    TVectorD fakes= _res->Vfakes();
    Double_t fac= _res->Vmeasured().Sum();
    if (fac!=0.0) fac=  Vmeasured().Sum() / fac;
    if (_verbose>=1) cout << "Subtract " << fac*fakes.Sum() << " fakes from measured distribution" << endl;
    for (Int_t i= 1; i<=_nm; i++)
      _meas1d->SetBinContent (i, _meas1d->GetBinContent(i)-(fac*fakes[i-1]));
  }
  // if (_verbose>=1) std::cout << "RUN(Blobel) init " << nBins(_reshist) << " x " << nBins(_reshist) << " bins, kreg=" << _kreg << std::endl;
  if(!_meas1d) throw std::runtime_error("no meas1d given!");
  if(!_train1d) throw std::runtime_error("no train1d given!");
  if(!_truth1d) throw std::runtime_error("no truth1d given!");
  if(!_reshist) throw std::runtime_error("no reshist given!");

  _nb= _nm > _nt ? _nm : _nt;
  \
  // cout << "Likelihood is: " << Likelihood << endl;
  //first find objective function of response and observed
  //Normalize response hist
  TMatrixD reshistmatrix(_nb, _nb);
  TVectorD reshistnorm(_nb);
  for(int i = 0; i < _nb; i++){
    for(int j = 0; j < _nb; j++){
      reshistnorm(i) += _reshist->GetBinContent(i+1,j+1);
      reshistmatrix(i, j) = _reshist->GetBinContent(i+1,j+1);
    }
  }
  reshistmatrix.Print();
  bool tryNorm = true;
  // tryNorm = false;
  // reshistmatrix(0,0)=0.23101425;
  // reshistmatrix(0,1)=0.00154646;
  // reshistmatrix(0,0)=0.0;
  // reshistmatrix(0,0)=0.0;
  // reshistmatrix(1,0)=0.00771165;
  // reshistmatrix(1,1)=0.5302493;
  // reshistmatrix(1,2)=0.02543888;
  // reshistmatrix(1,3)=0.0;
  // reshistmatrix(2,0)=0.0;
  // reshistmatrix(2,1)=0.02940941;
  // reshistmatrix(2,2)=0.6366643;
  // reshistmatrix(2,3)=0.01407598;
  // reshistmatrix(3,0)=0.0;
  // reshistmatrix(3,1)=0.0;;
  // reshistmatrix(3,2)=0.0022497;
  // reshistmatrix(3,3)=0.5336764;
  // reshistmatrix(0,0)=0.23101425;
  if(tryNorm){
    for(int i = 0; i < _nb; i++){
      for(int j = 0; j < _nb; j++){
        if(reshistnorm(i) != 0.0){
          reshistmatrix(i, j) = reshistmatrix(i,j)/reshistnorm(i);
        }
      }
    }
  }
  TVectorD measured(_nb);
  double flatSpread = 0.0;
  for(int i = 0; i < _nb; i++){
    measured(i) = _meas1d->GetBinContent(i+1);
    flatSpread += measured(i);
  }
  measured.Print();

  //initial estimate is flat estimate using total values of meas1d
  TVectorD _est(_nb);
  double total = 0.0;
  for(int i = 0; i < _nb; i++){
    // total += _meas1d->GetBinContent(i);
    total += measured(i);
  }
  // for(int i = 0; i < _nb; i++){
  //   measured(i) = _est(i);
  // }
  total = total/_nb;
  for(int i = 0; i < _nb; i++){
    _est(i) = total;
  }
  //Get Gradient Vector for initial estimate
  TVectorD _grad(_nb);
  _grad *= 0.0;
  TVectorD _lsqGrad(_nb);
  _lsqGrad *= 0.0;
  for(int i = 0; i < _nb; i++){
    double sumValue = 0.0;
    for (int j = 0; j < _nb; j++) {
      double dot_product = 0.0;
      for (int k = 0; k < _nb; k++) {
        dot_product += reshistmatrix(j,k)*_est(k);
      }
      if(measured(j) != 0){
        sumValue -=
         (reshistmatrix(j,i)*(measured(j)-dot_product)/measured(j));
      }
    }
    _lsqGrad(i) = sumValue;
  }
  //Get Hessian for initial estimated
  TMatrixD _hess(_nb,_nb);
  _rec.ResizeTo (_nb);
  //Get Hessian Matrix for lsq objective
  TMatrixD _lsqHess(_nb, _nb);
  _hess *= 0.0;
  _lsqHess *= 0.0;
  for(int i = 0; i < _nb; i++){
    for (int j = 0; j < _nb; j++) {
      double sumValue = 0.0;
      for (int k = 0; k < _nb; k++) {
        if(measured(k) != 0){
          sumValue += reshistmatrix(k,i)*reshistmatrix(k,j)/measured(k);
        }
      }
      _lsqHess(i,j) = sumValue;
    }
  }
  _lsqHess.Print();
  _lsqGrad.Print();
  TMatrixD Identity(_nb, _nb);
  Identity *= 0.0;
  for(int i = 0; i < _nb; i++){
    Identity(i,i) = 1;
  }
  double regParm = 0;
  TMatrixD curvatureMatrix(_nb, _nb);
  curvatureMatrix *= 0.0;
  TMatrixD _lsqHessInv(_nb, _nb);
  TVectorD _estOld(_nb);
  _estOld = _est;
  _lsqHessInv = _lsqHess.Invert();
  _est -= _lsqHessInv*_lsqGrad;
  double loss = GetLoss(_nb, measured, reshistmatrix, &_est);
  int itNum = 2;
  while(abs(loss) > 1000){
    _hess = GetHess(_nb, measured, reshistmatrix, &_est);
    _grad = GetGrad(_nb, measured, reshistmatrix, &_est);
    TMatrixD _hessInv = _hess.Invert();
    TVectorD mod = _hessInv*_grad;
    mod.Print();
    _est -= _hessInv*_grad;
    TMatrixD eigenVectors(_nb, _nb);
    TMatrixD eigenValues(_nb, _nb);
    TMatrixDEigen eigenHess(_hess);
    eigenValues = eigenHess.GetEigenValues();
    eigenVectors = eigenHess.GetEigenVectors();

    // eigenValues.Print();
    // eigenVectors.Print();
    // TMatrixD eVTrans = eigenVectors.Transpose(eigenVectors);
    // TMatrixD DHalf = eigenValues;
    // TMatrixD DInv = eigenValues.Invert();
    // TMatrixD DHalfInv = DInv;
    // for(int i = 0; i < _nb; i ++){
    //   DHalfInv(i,i) = sqrt(DInv(i,i));
    //   DHalf(i,i) = sqrt(DHalf(i,i));
    // }
    // TMatrixD curvShift(_nb, _nb);
    // curvShift *= 0.0;
    // curvShift = DHalfInv*eVTrans*curvatureMatrix*eigenVectors*DHalfInv;
    // TMatrixD curveEVec(_nb, _nb);
    // TMatrixD curveEVal(_nb, _nb);
    // TMatrixDEigen eigenCurve(curvShift);
    // curveEVec = eigenCurve.GetEigenVectors();
    // curveEVal = eigenCurve.GetEigenValues();
    // regParm = GetRegParm(curveEVal, 0);
    // TMatrixD compOne(_nb, _nb);
    // compOne = (Identity + regParm*curveEVal);
    // compOne = compOne.Invert();
    // TMatrixD compTwo(_nb, _nb);
    // compTwo = (eigenVectors*DHalfInv*curveEVec);
    // compTwo = compTwo.Transpose(compTwo);
    // TVectorD compThree(_nb);
    // compThree = (_hess * _est - _grad);
    // TVectorD estESpace = compOne*compTwo*compThree;
    // _est = (eigenVectors * DHalf * curveEVec) * estESpace;


    loss = GetLoss(_nb, measured, reshistmatrix, &_est);
    cout << "Loss on itteration " << itNum << " is " << loss << endl;
    itNum++;
    if(itNum >= 4){
      break;
    }
  }
  _rec = _est;
  _unfolded= true;
  _haveCov=  false;
}

TMatrixD
RooUnfoldBlobel::GetHess(Int_t _nb, TVectorD measured,  TMatrixD reshistmatrix, TVectorD *_est){
  TVectorD est = *_est;
  TMatrixD _hess(_nb,_nb);
  _hess *= 0.0;
  for(int i = 0; i < _nb; i++){
    for (int j = 0; j < _nb; j++) {
      double sumValue = 0.0;
      for (int k = 0; k < _nb; k++) {
        double dot_product = 0.0;
        for (int l = 0; l < _nb; l++) {
          dot_product += reshistmatrix(k,l)*est(l);
        }
        if(dot_product != 0){
          sumValue += measured(k)*reshistmatrix(k,i)*reshistmatrix(k,j)/(dot_product*dot_product);
        }
      }
      _hess(i,j) = sumValue;
    }
  }
  return _hess;
}

TVectorD
RooUnfoldBlobel::GetGrad(Int_t _nb, TVectorD measured, TMatrixD reshistmatrix, TVectorD *_est){
  TVectorD est = *_est;
  TVectorD _grad(_nb);
  _grad *= 0.0;
  double sumValue;
  double dot_product;
  for (int i = 0; i < _nb; i++) {
    sumValue = 0.0;
    for (int j = 0; j < _nb; j++) {
      dot_product = 0.0;
      for (int k = 0; k < _nb; k++) {
        dot_product += reshistmatrix(j,k)*est(k);
      }
      if(dot_product != 0){
        sumValue += reshistmatrix(j,i)-measured(j)*reshistmatrix(j,i)/dot_product;
      }
    }
    _grad(i) = sumValue;
  }
  return _grad;
}

Double_t
RooUnfoldBlobel::GetLoss(Int_t _nb, TVectorD measured, TMatrixD reshistmatrix, TVectorD *_est){
  double loss = 0.0;
  TVectorD est = *_est;
  for(int i = 0; i < _nb; i++){
    double sum = 0.0;
    for(int j = 0; j < _nb; j++){
      sum += reshistmatrix(i,j) * est(j);
    }
    double value = 0.0;
    if(sum > 0.0){
      value = sum - (measured(i) * log(sum));
    }
    loss += value;
  }
  return loss;
}

Double_t
RooUnfoldBlobel::GetRegParm(TMatrixD diag, int n_df){
 /* Need to work on this next */
 return 0.0;
}

void
RooUnfoldBlobel::GetCov()
{
  return;
//   if (!_svd) return;
//   Bool_t oldstat= TH1::AddDirectoryStatus();
//   TH1::AddDirectory (kFALSE);
//
//   TH2D *unfoldedCov= 0, *adetCov= 0;
//   //Get the covariance matrix for statistical uncertainties on the measured distribution
//   if (_dosys!=2) unfoldedCov= _svd->GetXtau();
//   //Get the covariance matrix for statistical uncertainties on the response matrix
//   //Uses Poisson or Gaussian-distributed toys, depending on response matrix histogram's Sumw2 setting.
//   if (_dosys)        adetCov= _svd->GetAdetCovMatrix (_NToys);
//
//   _cov.ResizeTo (_nt, _nt);
//   for (Int_t i= 0; i<_nt; i++) {
//     for (Int_t j= 0; j<_nt; j++) {
//       Double_t v  = 0;
//       if (unfoldedCov) v  = unfoldedCov->GetBinContent(i+1,j+1);
//       if (adetCov)     v += adetCov    ->GetBinContent(i+1,j+1);
//       _cov(i,j)= v;
//     }
//   }
//
//   delete adetCov;
// #ifdef TSVDUNFOLD_LEAK
//   delete unfoldedCov;
// #endif
//   TH1::AddDirectory (oldstat);
//
//   _haveCov= true;
}

void
RooUnfoldBlobel::GetWgt()
{
  return;
  // //! Get weight matrix
  // if (_dosys) RooUnfold::GetWgt();   // can't add sys errors to weight, so calculate weight from covariance
  // if (!_svd) return;
  // Bool_t oldstat= TH1::AddDirectoryStatus();
  // TH1::AddDirectory (kFALSE);
  //
  // //Get the covariance matrix for statistical uncertainties on the measured distribution
  // TH2D* unfoldedWgt= _svd->GetXinv();
  //
  // _wgt.ResizeTo (_nt, _nt);
  // for (Int_t i= 0; i<_nt; i++) {
  //   for (Int_t j= 0; j<_nt; j++) {
  //     _wgt(i,j)= unfoldedWgt->GetBinContent(i+1,j+1);
  //   }
  // }

// #ifdef TSVDUNFOLD_LEAK
//   delete unfoldedWgt;
// #endif
//   TH1::AddDirectory (oldstat);
//
//   _haveWgt= true;
}

void
RooUnfoldBlobel::GetSettings(){
    _minparm=0;
    _maxparm= _meas ? _meas->GetNbinsX() : 0;
    _stepsizeparm=1;
    _defaultparm=_maxparm/2;
}

void
RooUnfoldBlobel::FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC, Int_t _nb) const
{
   Double_t eps = 0.00001;

   Int_t ndim = tCurv.GetNrows();

   // Init
   tCurv *= 0;
   tC    *= 0;

  for (Int_t i=0; i< _nb; i++) {
     if (i > 0)      tC(i,i-1) = 1.0;
     if (i < _nb-1)  tC(i,i+1) = 1.0;
     tC(i,i) = -2.0;
  }
  tC(0,0) = -1.0;
  tC(_nb-1,_nb-1) = -1.0;


   // Add epsilon to avoid singularities
   for (Int_t i=0; i<ndim; i++) tC(i,i) = tC(i,i) + eps;

   //Get curvature matrix
   for (Int_t i=0; i<ndim; i++) {
      for (Int_t j=0; j<ndim; j++) {
         for (Int_t k=0; k<ndim; k++) {
            tCurv(i,j) = tCurv(i,j) + tC(k,i)*tC(k,j);
         }
      }
   }
}
