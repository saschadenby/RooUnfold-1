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
#include "TCanvas.h"
#include "TMath.h"

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
  //Measure and check for bin correlation
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
  //Set number of bins equal to size of largest out of measured and truth histogram
  _nb= _nm > _nt ? _nm : _nt;

  //Make sure kreg is not greater than number of bins
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
  _rec.ResizeTo (_nb);
  //Fill response histogram (When saying fill in this section meaning turning histogram into vector)
  TMatrixD reshist(_nb,_nb);
  for(int i = 1; i <= _nb; i++){
    for (int j = 1; j <= _nb; j++) {
      reshist((i-1),(j-1)) = _reshist->GetBinContent(i,j);
    }
  }
  //Fill Measured
  TVectorD meas1d(_nb);
  for (int i = 1; i <= _nb; i++) {
    meas1d((i-1)) = _meas1d->GetBinContent(i);
  }
  //Fill Train
  TVectorD train1d(_nb);
  for (int i = 1; i <= _nb; i++) {
    train1d((i-1)) = _train1d->GetBinContent(i);
  }
  //Fill truth
  TVectorD truth1d(_nb);
  for (int i = 1; i <= _nb; i++) {
    truth1d((i-1)) = _truth1d->GetBinContent(i);
  }
  // initial estimate is flat estimate using total values of meas1d
  TVectorD _est(_nb);
  int total = _meas1d->Integral()/_nb;
  for(int i = 0; i < _nb; i++){
    _est(i) = total;
  }
  //Normalize response
  for (int i = 0; i < _nb; i++) {
    for (int j = 0; j < _nb; j++) {
      reshist(i,j) = reshist(i,j)/train1d(i);
    }
  }
  //Find Loss
  double loss = GetLoss(_nb, meas1d, reshist, &_est);
  cout << "Initial loss is: " << loss << endl;
  TMatrixD resInv = InvertMatrix(reshist, _nb);
  TVectorD grad = GetGrad(_nb, meas1d, reshist, &_est);
  TMatrixD hess = GetHess(_nb, meas1d, reshist, &_est);
  TMatrixD hessInv= InvertMatrix(hess, _nb);
  TMatrixD check(_nb,_nb);
  check = hessInv * hess;
  TVectorD mod(_nb);
  mod = hess * grad;
  _est -= mod;
  loss = GetLoss(_nb, meas1d, reshist, &_est);
  for (int i = 0; i < _kreg; i++) {
    grad = GetGrad(_nb, meas1d, reshist, &_est);
    hess = GetHess(_nb, meas1d, reshist, &_est);
    hessInv = InvertMatrix(hess, _nb);
    mod = hess * grad;
    _est -= mod;
    loss = GetLoss(_nb, meas1d, reshist, &_est);
  }
  double estInt = 0.0;
  double truthInt = 0.0;
  double measInt = 0.0;
  double trainInt = 0.0;
  for (int i = 0; i < _nb; i++) {
    _est(i) -= total;
    measInt += meas1d(i);
    estInt += _est(i);
    truthInt += truth1d(i);
    trainInt += train1d(i);
  }
  double r1 = truthInt/trainInt;
  double resultRatio = measInt/estInt;
  _est *= resultRatio;
  _est *= r1;
  _rec = _est;
  _unfolded= true;
  _haveCov=  false;
  return;
}

TMatrixD
RooUnfoldBlobel::InvertMatrix(TMatrixD a, Int_t _nb){
  TDecompSVD svd(a);
  Bool_t ok = svd.Decompose();
  TMatrixD b(_nb, _nb);
  if (ok)
  b = svd.Invert();
  else {
  cout << "SVD failed, condition: " << svd.Condition() <<endl;
  a.Print();
  }
  return b;
}

TMatrixD
RooUnfoldBlobel::GetHess(Int_t _nb, TVectorD measured,  TMatrixD reshistmatrix, TVectorD *_est){
  TVectorD est = *_est;
  TMatrixD _hess(_nb,_nb);
  _hess *= 0.0;
  for(int i = 0; i < _nb; i++){
    for (int j = 0; j < _nb; j++) {
      double sumVal = 0.0;
      for (int k = 0; k < _nb; k++){
        sumVal += measured(k)*reshistmatrix(k,i)*reshistmatrix(k,j)/(est(k)*est(k));
      }
      _hess(i,j) = sumVal;
    }
  }
  return _hess;
}

TVectorD
RooUnfoldBlobel::GetGrad(Int_t _nb, TVectorD measured, TMatrixD reshistmatrix, TVectorD *_est){
  TVectorD est = *_est;
  TVectorD _grad(_nb);
  _grad *= 0.0;
  for(int i = 0; i < _nb; i++){
    double sumVal = 0.0;
    for (int j = 0; j < _nb; j++) {
      sumVal += (measured(j)/est(j)-1)*reshistmatrix(j,i);
    }
    _grad(i) -= sumVal;
  }
  return _grad;
}

Double_t
RooUnfoldBlobel::GetLoss(Int_t _nb, TVectorD measured, TMatrixD reshistmatrix, TVectorD *_est){
  TVectorD est = *_est;
  double loss = 0.0;
  for (int i = 0; i < _nb; i++) {
    loss += measured(i)*log(est(i))-est(i);
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

TMatrixD
RooUnfoldBlobel::FillCurvatureMatrix(Int_t _nb)
{
   TMatrixD cv(_nb,_nb);
   for (int i = 3; i < _nb; i++) {
     cv(i-3,i) = 1;
     cv(i,i-3) = 1;
   }
   for (int i = 3; i < _nb-3; i++) {
     cv(i,i) = 16;
     cv(i,i-1) = -9;
     cv(i,i+1) = -9;
     cv(i+1,i) = -9;
     cv(i-1,i) = -9;
   }
   cv(0,0) = 2;
   cv(1,1) = 8;
   cv(2,2) = 14;
   cv(_nb-1,_nb-1) = 2;
   cv(_nb-2,_nb-2) = 8;
   cv(_nb-3,_nb-3) = 14;
   cv(0,1) = -3;
   cv(1,0) = -3;
   cv(1,2) = -6;
   cv(2,1) = -6;
   cv(_nb-1,_nb-2) = -3;
   cv(_nb-2,_nb-1) = -3;
   cv(_nb-2,_nb-3) = -6;
   cv(_nb-3,_nb-2) = -6;
   return cv;
}
