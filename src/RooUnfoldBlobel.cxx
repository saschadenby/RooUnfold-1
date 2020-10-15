#include "RooUnfoldBlobel.h"

#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "TClass.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include "RooUnfoldHelpers.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldTH1Helpers.h"
#include "RooUnfoldFitHelpers.h"

using namespace RooUnfolding;

template<class Hist,class Hist2D>
RooUnfoldBlobelT<Hist,Hist2D>::RooUnfoldBlobelT (const RooUnfoldBlobelT<Hist,Hist2D>& rhs)
  : RooUnfoldT<Hist,Hist2D> (rhs)
{
  //duplicate constructor
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D>
RooUnfoldBlobelT<Hist,Hist2D>::RooUnfoldBlobelT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg,
						 const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  //! Constructor with response matrix object and measured unfolding input histogram.
  //! The regularisation parameter is kreg.
  Init();
}

template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::Reset()
{
  // destroy and re-initialize this object
  Destroy();
  Init();
  RooUnfoldT<Hist,Hist2D>::Reset();
}

template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::Destroy()
{
  //! delete all members of this object
  delete this->_meas1d;
  delete this->_train1d;
  delete this->_truth1d;
  delete this->_reshist;
}

template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::Init()
{
  //! initialize this object with zero values
  this->_meas1d= this->_train1d= this->_truth1d= 0;
  this->_reshist= this->_meascov= 0;
  GetSettings();
}

template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::Assign (const RooUnfoldBlobelT<Hist,Hist2D>& rhs)
{
  //! assign data from another instance
  RooUnfoldT<Hist,Hist2D>::Assign (rhs);
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::CopyData (const RooUnfoldBlobelT<Hist,Hist2D>& rhs)
{
  //! copy data from another instance
  this->_kreg= rhs._kreg;
}


template<class Hist,class Hist2D> void RooUnfoldBlobelT<Hist,Hist2D>::PrepareHistograms() const{
  //! fill all the histogram members
  this->_meas1d = this->_meas;
  this->_train1d= this->_res->Hmeasured();
  this->_truth1d= this->_res->Htruth();
  this->_reshist= this->_res->Hresponse();
}

namespace{
  TH1* histNoOverflow(const TH1* hist, bool overflow){
    return createHist<TH1>(h2v(hist,overflow),h2ve(hist,overflow),name(hist),title(hist),vars(hist),overflow);
  }
  void subtract(TH1* hist, const TVectorD& vec, double fac){
    for (Int_t i= 1; i<=hist->GetNbinsX()+1; i++){
      hist->SetBinContent (i, hist->GetBinContent(i)-(fac*vec[i-1]));
    }
  }
}
template<> void RooUnfoldBlobelT<TH1,TH2>::PrepareHistograms() const {
  //! fill all the histogram members
  TH1* meas1d = ::histNoOverflow (this->_meas,             this->_overflow);
  TH1* train1d= ::histNoOverflow (this->_res->Hmeasured(), this->_overflow);
  TH1* truth1d= ::histNoOverflow (this->_res->Htruth(),    this->_overflow);
  TH2* reshist= this->_res->HresponseNoOverflow();
  RooUnfolding::resize (meas1d,  this->_nm);
  RooUnfolding::resize (train1d, this->_nm);
  RooUnfolding::resize (truth1d, this->_nt);
  RooUnfolding::resize (reshist, this->_nm, this->_nt);

  // Subtract fakes from measured distribution
  if (this->_res->HasFakes()) {
    TVectorD fakes(this->_res->Vfakes());
    Double_t fac= this->_res->Vmeasured().Sum();
    if (fac!=0.0) fac=  this->Vmeasured().Sum() / fac;
    if (this->_verbose>=1) std::cout << "Subtract " << fac*fakes.Sum() << " fakes from measured distribution" << std::endl;
    ::subtract(meas1d,fakes,fac);
  }

  delete this->_meas1d  ;
  delete this->_train1d ;
  delete this->_truth1d ;
  delete this->_reshist ;

  this->_meas1d  = meas1d  ;
  this->_train1d = train1d;
  this->_truth1d = truth1d;
  this->_reshist = reshist ;
}

template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::Unfold() const
{
  //! perform the unfolding
  if (this->_res->GetDimensionTruth() != 1 || this->_res->GetDimensionMeasured() != 1) {
    std::cerr << "RooUnfoldBlobelT may not work very well for multi-dimensional distributions" << std::endl;
  }
  if (this->_res->GetDimensionTruth() != this->_res->GetDimensionMeasured()){
    std::cerr << "RooUnfoldBlobelT may not work very well for distributions of different sizes" << std::endl;
  }
  if (this->_kreg < 0) {
    std::cerr << "RooUnfoldBlobelT invalid kreg: " << this->_kreg << std::endl;
    return;
  }

  this->PrepareHistograms();

  if (this->_verbose>=1) std::cout << "RUN(Blobel) init " << nBins(this->_reshist,X) << " x " << nBins(this->_reshist,Y) << " bins, kreg=" << this->_kreg << std::endl;
  if(!this->_meas1d) throw std::runtime_error("no meas1d given!");
  if(!this->_train1d) throw std::runtime_error("no train1d given!");
  if(!this->_truth1d) throw std::runtime_error("no truth1d given!");
  if(!this->_reshist) throw std::runtime_error("no reshist given!");

  //Find Hessian of Response Matrix
  int Hess_a = nBins(this->_reshist,X);
  int Hess_b = nBins(this->_reshist,Y);
  double Hessian[Hess_a][Hess_a];
  for(int i = 0; i < Hess_a; i++){
    for(int j = 0; j < Hess_a; j++){
      Hessian[i][j] = 0;
    }
  }
  for(int i = 0; i < Hess_a; i++){
    for(int j = 0; j < Hess_a; j++){
      double value = 0.0;
      double dot = 0.0;
      for(int k = 0; k < Hess_b; k++){
        for(int l = 0; l < Hess_b; l++){
          dot += this->_reshist[k][l] * this->_meas1d[l];
        }
        value += this->_truth1d[i] * this->_reshist[k][i] * this->_reshist[k][j] / (dot * dot);
        dot = 0.0;
      }
      Hessian[i][j] = value;
    }
 }}

  //convert hessian into diagonalization and basis matrixes
  double Basis[Hess_a][Hess_a];
  double Diagonal[Hess_a][Hess_a];
  double Last[Hess_a][Hess_a];
  double Final[Hess_a][Hess_a];
  for(int i = 0; i < Hess_a; i++){
    for(int j = 0; j < Hess_a; j++){
      Basis[i][j] = 0.0;
      Last[i][j] = 0.0;
      Diagonal[i][j] = 0.0;
      Final[i][j] = 0.0;
    }
  }

  //apply hessian decomposition to reco distribution
  for(int i = 0; i < Hess_a; i++){
    double temp_value = 0.0;
    for(int j = 0; j < Hess_a; j++){
      temp_value += Final[i][j] * this->_meas1d[j];
    }
    this->_train1d[i] = temp_value;
  }


//establish eigen library object matrix in order to perform eigen function
Eigen::Matrix<double, Hess_a, Hess_a> Eigenvalues;

//convert Hessian
for(int i = 0; i < Hess_a; i++){
  for(int j = 0; j < Hess_a; j++){
    Eigenvalues(i,j) = Hessian[i][j];
  }
 }

//find eigenvalues
Eigen::SelfAdjointEigenSolver<Matrix<double, Hess_a, Hess_a>> eigensolver(Eigenvalues);
if (eigensolver.info() != Success) abort();
Eigen::VectorXcd eivals = eigensolver.eigenvalues();
for(int i = 0; i < Hess_a; i++){
  Diagonal[i][i] = eivals(i);
 }

//Now find eigenvectors
for(int i = 0; i < Hess_a; i++){
  Eigen::VectorXcd eivals = eigensolver.eigenvectors().col(i);
  for(int j = 0; j < Hess_a; j++){
    Basis_vals[i][j] = eivals(j);
    Last[j][i] = eivals(j);
  }
 }

//calculate final product of these three to apply to the distribution
double Stepping_Stone[Hess_a][Hess_a];
for(int i = 0; i < Hess_a; i++){
  for(int j = 0; j < Hess_a; j++){
    Stepping_Stone[i][j] = 0.0;
  }
 }
for (int i = 0; i < Hess_a; i++) {
  for (int  j= 0; j < Hess_a; j++) {
    for (int k = 0; k < Hess_a; k++) {
      Stepping_Stone[i][j] += Basis_vals[i][k] * Diagonal[k][j];
    }
  }
 }
for (int i = 0; i < Hess_a; i++) {
  for (int  j= 0; j < Hess_a; j++) {
    for (int k = 0; k < Hess_a; k++) {
      Final[i][j] += Stepping_Stone[i][k] * Diagonal[k][j];
    }
  }
 }

double truth_array[Hess_a];
for(int i = 0; i < Hess_a; i++){
  truth_array[i] = 0.0;
 }

//apply hessian decomposition to reco distribution
for(int i = 0; i < Hess_a; i++){
  double temp_value = 0.0;
  for(int j = 0; j < Hess_a; j++){
    temp_value += Final[i][j] * this->_meas1d[j];
  }
  truth_array[i] = temp_value;
 }

TvectorD final_values(truth_array);
this->Vtruth() = final_values;

  if (this->_verbose>=2) {
    printTable (std::cout, h2v(this->_truth1d), h2v(this->_train1d), h2v(this->_meas1d), this->_cache._rec);
    TMatrixD resmat(h2m(this->_reshist));
    printMatrix(resmat,"BlobelUnfold response matrix");
  }

  this->_cache._unfolded= true;
  this->_cache._haveCov=  false;
}

template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::GetCov() const
{
  //! Get covariance matrix, need to figure out still
  this->_cache._haveCov= true;
  return;
}

template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::GetWgt() const
{
  //! Get weight matrix
  if (this->_dosys) RooUnfoldT<Hist,Hist2D>::GetWgt();
  if (!this->_svd) return;

  this->_cache._wgt.ResizeTo(this->_nt, this->_nt);

  this->_cache._wgt = this->_svd->GetXinv();

  this->_cache._haveWgt= true;
}


template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::GetSettings() const {
  this->_cache._minparm=0;
  this->_cache._maxparm= this->_meas ? nBins(this->_meas,X) : 0;
  this->_cache._stepsizeparm=1;
  this->_cache._defaultparm=this->_cache._maxparm/2;
}

template<class Hist,class Hist2D>
RooUnfoldBlobelT<Hist,Hist2D>::RooUnfoldBlobelT()
  : RooUnfoldT<Hist,Hist2D>()
{
  //! Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldBlobelT<Hist,Hist2D>::RooUnfoldBlobelT (const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldBlobelT<Hist,Hist2D>::RooUnfoldBlobelT (const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldBlobelT<Hist,Hist2D>& RooUnfoldBlobelT<Hist,Hist2D>::operator= (const RooUnfoldBlobelT<Hist,Hist2D>& rhs)
{
  //! Assignment operator for copying RooUnfoldSvdT settings.
  Assign(rhs);
  return *this;
}

template<class Hist,class Hist2D>
RooUnfoldBlobelT<Hist,Hist2D>::~RooUnfoldBlobelT()
{
  Destroy();
}


template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::SetKterm (Int_t kreg)
{
  //! Set regularisation parameter
  this->_kreg= kreg;
}


template<class Hist,class Hist2D> Int_t
RooUnfoldBlobelT<Hist,Hist2D>::GetKterm() const
{
  //! Return regularisation parameter
  return this->_kreg;
}

template<class Hist,class Hist2D> void
RooUnfoldBlobelT<Hist,Hist2D>::SetRegParm (Double_t parm)
{
  //! Set regularisation parameter
  SetKterm(parm);
}

template<class Hist,class Hist2D> Double_t
RooUnfoldBlobelT<Hist,Hist2D>::GetRegParm() const
{
  //! Return regularisation parameter
  return GetKterm();
}

template class RooUnfoldBlobelT<TH1,TH2>;
ClassImp (RooUnfoldBlobel)

#ifndef NOROOFIT
typedef RooUnfoldBlobelT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldBlobel;
template class RooUnfoldBlobelT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
ClassImp (RooFitUnfoldBlobel)
#endif
