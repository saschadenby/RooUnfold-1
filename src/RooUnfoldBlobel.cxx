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

  // cout << "Likelihood is: " << Likelihood << endl;
  //first find objective function of response and observed
  //Normalize response hist
  TMatrixD reshistmatrix(_nb, _nb);
  // TVectorD reshistnorm(_nb);
  // for(int i = 0; i < _nb; i++){
  //   for(int j = 0; j < _nb; j++){
  //     reshistnorm(i) += _reshist->GetBinContent(i+1,j+1);
  //     reshistmatrix(i, j) = _reshist->GetBinContent(i+1,j+1);
  //   }
  // }
  // bool tryNorm = true;
  // tryNorm = false;
  reshistmatrix(0,0)=0.23101425;
  reshistmatrix(0,1)=0.00154646;
  reshistmatrix(0,0)=0.0;
  reshistmatrix(0,0)=0.0;
  reshistmatrix(1,0)=0.00771165;
  reshistmatrix(1,1)=0.5302493;
  reshistmatrix(1,2)=0.02543888;
  reshistmatrix(1,3)=0.0;
  reshistmatrix(2,0)=0.0;
  reshistmatrix(2,1)=0.02940941;
  reshistmatrix(2,2)=0.6366643;
  reshistmatrix(2,3)=0.01407598;
  reshistmatrix(3,0)=0.0;
  reshistmatrix(3,1)=0.0;;
  reshistmatrix(3,2)=0.0022497;
  reshistmatrix(3,3)=0.5336764;
  reshistmatrix(0,0)=0.23101425;
  // if(tryNorm){
  //   for(int i = 0; i < _nb; i++){
  //     for(int j = 0; j < _nb; j++){
  //       if(reshistnorm(i) != 0.0){
  //         reshistmatrix(i, j) = reshistmatrix(i,j)/reshistnorm(i);
  //       }
  //     }
  //   }
  // }
  TVectorD measured(_nb);
  // for(int i = 0; i < _nb; i++){
  //   measured(i) = _meas1d->GetBinContent(i+1);
  // }

  //initial estimate is flat estimate using total values of meas1d
  TVectorD _est(_nb);
  _est(0) = 290;
  _est(1) = 3157;
  _est(2) = 3784;
  _est(3) = 325;
  double total = 0.0;
  for(int i = 0; i < _nb; i++){
    // total += _meas1d->GetBinContent(i);
    total += _est(i);
  }
  for(int i = 0; i < _nb; i++){
    measured(i) = _est(i);
  }
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
      sumValue -=
       (reshistmatrix(j,i)*(measured(j)-dot_product)/measured(j));
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
        sumValue += reshistmatrix(k,i)*reshistmatrix(k,j)/measured(k);
      }
      _lsqHess(i,j) = sumValue;
    }
  }

  TMatrixD _lsqHessInv(_nb, _nb);
  TVectorD _estOld(_nb);
  _estOld = _est;
  _lsqHessInv = _lsqHess.Invert();
  _est -= _lsqHessInv*_lsqGrad;
  double loss = GetLoss(_nb, measured, reshistmatrix, &_est);

  _grad = GetGrad(_nb, measured, reshistmatrix, &_est);
  cout << "Gradient itteration 2" << endl;
  _grad.Print();

  while(loss < -1000){


    break;
    loss = GetLoss(_nb, measured, reshistmatrix, &_est);
  }
  _rec = _est;
  _unfolded= true;
  _haveCov=  false;
  return;


  // TMatrixD eigenVectors(_nb, _nb);
  // TMatrixD eigenValues(_nb, _nb);
  // TMatrixDEigen eigenHess(_lsqHess);
  // eigenValues = eigenHess.GetEigenValues();
  // eigenVectors = eigenHess.GetEigenVectors();
  // // eigenValues.Print();
  // TMatrixD DInv = eigenValues.Invert();
  // // DInv.Print();
  // TMatrixD eVTrans = eigenVectors.Transpose(eigenVectors);
  // // eigenVectors.Print();
  // // eVTrans.Print();
  //
  // //Invert Hessian of lsq objective
  // // _lsqHess.Print();
  // // TDecompSVD svd(_lsqHess);
  // // Bool_t solve;
  // // TMatrixD pseudo  = svd.Invert();
  //
  // _lsqHessInv = eigenVectors * DInv * eVTrans;
  // // _lsqHessInv.Print();
  // // _lsqHessInv.Print();
  // // for (int i = 0; i < _nb; i++) {
  // //   double dot_prod = 0.0;
  // //   for (int j = 0; j < _nb; j++) {
  // //     dot_prod += _lsqHessInv(i,j) * _lsqGrad(j);
  // //   }
  // //   _est(i) -= dot_prod;
  // // }
  // _estOld = _est;
  // TVectorD estAdj(_nb);
  // // estAdj = eigenVectors * (DInv * (eVTrans * _lsqGrad));
  // TMatrixD Dsqrt(_nb,_nb);
  // Dsqrt *= 0.0;
  // for (int i = 0; i < _nb; i++) {
  //   double value = DInv(i,i);
  //   Dsqrt(i,i) = sqrt(value);
  // }
  //
  // estAdj = Dsqrt * (eVTrans * ((_lsqHess*_est)+_lsqGrad));
  // // estAdj.Print();
  // _est = eigenVectors * (Dsqrt * estAdj);
  // // _lsqHessInv.Print();
  //
  // TMatrixD mCurv(_nb, _nb), mC(_nb, _nb);
  // FillCurvatureMatrix(mCurv, mC, _nb);
  //
  // int nItterations = 0;
  // for(int i = 0; i < nItterations; i++){
  //   _grad = GetGrad(_nb, _meas1d, reshistmatrix, &_est, &_estOld);
  //   _hess = GetHess(_nb, _meas1d, reshistmatrix, &_est);
  //   _estOld = _est;
  //   TMatrixD _hessInv(_nb, _nb);
  //   TMatrixD eigenVectorsItt(_nb, _nb);
  //   TMatrixD eigenValuesItt(_nb, _nb);
  //   TMatrixDEigen eigenHessItt(_hess);
  //   eigenValuesItt = eigenHessItt.GetEigenValues();
  //   eigenVectorsItt = eigenHessItt.GetEigenVectors();
  //   // eigenValues.Print();
  //   TMatrixD DInvItt = eigenValuesItt.Invert();
  //   // DInv.Print();
  //   TMatrixD eVTransItt = eigenVectorsItt.Transpose(eigenVectors);
  //   _hessInv = eigenVectors * DInv * eVTrans;
  //   // _lsqHessInv.Print();
  //   for (int i = 0; i < _nb; i++) {
  //     double dot_prod = 0.0;
  //     for (int j = 0; j < _nb; j++) {
  //       dot_prod += _hessInv(i,j) * _grad(j);
  //     }
  //     _est(i) -= dot_prod;
  //   }
  // }
  _rec = _est;
  //Covariance is given by inverse of Hessian
  _unfolded= true;
  _haveCov=  false;
}

TMatrixD
RooUnfoldBlobel::GetHess(Int_t _nb, TH1D *_meas1d,  TMatrixD reshistmatrix, TVectorD *_est){
  TVectorD est = *_est;
  TMatrixD _hess(_nb,_nb);
  _hess *= 0.0;
  for (int j = 0; j < _nb; j++) {
    for (int k = 0; k < _nb; k++) {
      double sum = 0.0;
      for(int i = 0; i < _nb; i++){
        if(est(i) != 0){
          sum += reshistmatrix(i,j) * reshistmatrix(i,k) / (est(i)*est(i));
        }
      }
      _hess(j,k) = sum;
    }
  }
  return _hess;
}

TVectorD
RooUnfoldBlobel::GetGrad(Int_t _nb, TVectorD measured, TMatrixD reshistmatrix, TVectorD *_est){
  TVectorD est = *_est;
  TVectorD _grad(_nb);
  _grad *= 0.0;
  for(int i = 0; i < _nb; i ++){
    double sumValue = 0.0;
    for (int j = 0; j < _nb; j++) {
      double dot_product = 0.0;
      for (int k = 0; k < _nb; k++) {
        dot_product += reshistmatrix(j,k) * est(k);
      }
      sumValue += (reshistmatrix(j,i)-(measured(j)*reshistmatrix(j,i)/dot_product));
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
    if(sum != 0.0){
      value = sum - (measured(i) * log(sum));
    }
    loss += value;
  }
  return loss;
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
