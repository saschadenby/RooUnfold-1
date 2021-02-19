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

  //set up expected values given _reshist
  TVectorD expected(_nb);
  for(int i = 0; i < _nb; i++){
    for(int j = 0; j < _nb; j++){
      expected(i) += (_reshist->GetBinContent(i,j) * _meas1d->GetBinContent(j));
    }
  }

  //set up discrete curvature vector using elements of response histogram
  TVectorD curve(_nb);
  for(int i = 0; i < _nb; i++){
    if(i == 0){
      curve(i) = ((-2.0)* _meas1d->GetBinContent(i) + _meas1d->GetBinContent(i+1));
    } else if(i == (_nb - 1)){
      curve(i) = ((-2.0)* _meas1d->GetBinContent(i) + _meas1d->GetBinContent(i-1));
    } else{
      curve(i) = (_meas1d->GetBinContent(i-1) + _meas1d->GetBinContent(i+1) - ((-2.0)*_meas1d->GetBinContent(i)));
    }
  }

  //set up loss function
  // TVectorD LogL(_nb);
  double Likelihood = 0.0;
  for(int i = 0; i < _nb; i++){
    double bin_content = _meas1d->GetBinContent(i);
    double product = 0.0;
    if(bin_content == 0){
      product = 0.0;
    } else{
      product = expected(i) * log(bin_content);
    }
    // cout << "product is: " << product << endl;
    Likelihood += (bin_content - product);
  }
  // cout << "Likelihood is: " << Likelihood << endl;
  //first find objective function of response and observed

  //then find  gradient, right now f represents target histogram
  TVectorD target(_nb);
  for(int i = 0; i < _nb; i++){
    target(i) = 1.0;
  }
  TVectorD _grad(_nb);
  for(int i = 0; i < _nb; i++){
    double sum_grad = 0.0;
    for(int j = 0; j < _nb; j++){
      double dot_product = 0.0;
      for(int k = 0; k < _nb; k++){
        dot_product += _meas1d->GetBinContent(k) * _reshist->GetBinContent(k,j);
      }
      if(dot_product != 0){
        sum_grad += (_reshist->GetBinContent(j,i) - ((_meas1d->GetBinContent(j) * _reshist->GetBinContent(j,i))/dot_product));
      }
    }
    _grad(i) = sum_grad;
  }

  //now hessian
  TMatrixD _hess(_nb,_nb);
  for(int i = 0; i < _nb; i++){
    for(int j = 0; j < _nb; j++){
      _hess(i,j) = 0.0;
    }
  }

  //fill Hessian
  for(int i = 0; i < _nb; i++){
    for(int j = 0; j < _nb; j++){
      double _sum = 0.0;
      for(int k = 0; k < _nb; k++){
        double dot_product = 0.0;
        for(int l = 0; l < _nb; l++){
          dot_product += _meas1d->GetBinContent(l) * _reshist->GetBinContent(l,k);
        }
        if(dot_product != 0){
          _sum += ((_meas1d->GetBinContent(k) * _reshist->GetBinContent(k,j) * _reshist->GetBinContent(k,i)) /
                  (dot_product * dot_product));
        }
      }
      _hess(i,j) = _sum;
    }
  }
  _hess.Print();

  //then create Tikhonov matrix

  //Covariance is given by inverse of Hessian
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
