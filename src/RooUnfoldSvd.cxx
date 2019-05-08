//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      SVD unfolding. Just an interface to TSVDUnfold.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p>Links to TSVDUnfold class which unfolds using Singular Value Decomposition (SVD).</p>
<p>Regularisation parameter defines the level at which values are deemed to be due to statistical fluctuations and are cut out. (Default= number of bins/2)
<p>Returns errors as a full matrix of covariances
<p>Can only handle 1 dimensional distributions
<p>Can account for both smearing and biasing
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldSvd.h"

#include <iostream>
#include <iomanip>

#include "TClass.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include "RooUnfoldHelpers.h"
#include "RooUnfoldResponse.h"

using namespace RooUnfolding;

using std::cout;
using std::cerr;
using std::endl;

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const RooUnfoldSvdT<Hist,Hist2D>& rhs)
  : RooUnfoldT<Hist,Hist2D> (rhs)
{
  // Copy constructor.
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg,
                            const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  // The regularisation parameter is kreg.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg, Int_t ntoyssvd,
                            const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  // Constructor with old ntoyssvd argument. No longer required.
  Init();
  this->_NToys = ntoyssvd;
}

template<class Hist,class Hist2D> RooUnfoldSvdT<Hist,Hist2D>*
RooUnfoldSvdT<Hist,Hist2D>::Clone (const char* newname) const
{
  RooUnfoldSvdT<Hist,Hist2D>* unfold= new RooUnfoldSvdT<Hist,Hist2D>(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Reset()
{
  Destroy();
  Init();
  RooUnfoldT<Hist,Hist2D>::Reset();
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Destroy()
{
  delete this->_svd;
  delete this->_meas1d;
  delete this->_train1d;
  delete this->_truth1d;
  delete this->_reshist;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Init()
{
  this->_svd= 0;
  this->_meas1d= this->_train1d= this->_truth1d= 0;
  this->_reshist= this->_meascov= 0;
  GetSettings();
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Assign (const RooUnfoldSvdT<Hist,Hist2D>& rhs)
{
  RooUnfoldT<Hist,Hist2D>::Assign (rhs);
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::CopyData (const RooUnfoldSvdT<Hist,Hist2D>& rhs)
{
  this->_kreg= rhs._kreg;
}

template<class Hist,class Hist2D> typename RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold*
RooUnfoldSvdT<Hist,Hist2D>::Impl()
{
  return this->_svd;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Unfold()
{
  if (this->_res->GetDimensionTruth() != 1 || this->_res->GetDimensionMeasured() != 1) {
    cerr << "RooUnfoldSvdT may not work very well for multi-dimensional distributions" << endl;
  }
  if (this->_kreg < 0) {
    cerr << "RooUnfoldSvdT invalid kreg: " << this->_kreg << endl;
    return;
  }

  this->_nb= this->_nm > this->_nt ? this->_nm : this->_nt;

  if (this->_kreg > this->_nb) {
    cerr << "RooUnfoldSvdT invalid kreg=" << this->_kreg << " with " << this->_nb << " bins" << endl;
    return;
  }

  Bool_t oldstat= TH1::AddDirectoryStatus();
  this->_meas1d=  this->HistNoOverflow (this->_meas,             this->_overflow);
  this->_train1d= this->HistNoOverflow (this->_res->Hmeasured(), this->_overflow);
  this->_truth1d= this->HistNoOverflow (this->_res->Htruth(),    this->_overflow);
  this->_reshist= this->_res->HresponseNoOverflow();
  RooUnfolding::resize (this->_meas1d,  this->_nb);
  RooUnfolding::resize (this->_train1d, this->_nb);
  RooUnfolding::resize (this->_truth1d, this->_nb);
  RooUnfolding::resize (this->_reshist, this->_nb, this->_nb);

  // Subtract fakes from measured distribution
  if (this->_res->FakeEntries()) {
    TVectorD fakes(this->_res->Vfakes());
    Double_t fac= this->_res->Vmeasured().Sum();
    if (fac!=0.0) fac=  this->Vmeasured().Sum() / fac;
    if (this->_verbose>=1) cout << "Subtract " << fac*fakes.Sum() << " fakes from measured distribution" << endl;
    subtract(this->_meas1d,fakes,fac);
  }

  const TMatrixD& cov(this->GetMeasuredCov());

  if (this->_verbose>=1) cout << "SVD init " << nBins(this->_reshist,X) << " x " << nBins(this->_reshist,Y)
                        << " bins, kreg=" << this->_kreg << endl;
  this->_svd= new SVDUnfold (this->_meas1d, cov, this->_train1d, this->_truth1d, this->_reshist);

  this->_rec.ResizeTo (this->_nt);
  this->_rec = this->_svd->UnfoldV (this->_kreg);

  if (this->_verbose>=2) {
    printTable (cout, h2v(this->_truth1d), h2v(this->_train1d), h2v(this->_meas1d), this->_rec, this->_nb, this->_nb);
    TMatrixD resmat(h2m(this->_reshist));
    printMatrix(resmat,"SVDUnfold response matrix");
  }

  this->_unfolded= true;
  this->_haveCov=  false;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::GetCov()
{
  if (!this->_svd) return;
  this->_cov.ResizeTo(this->_nb,this->_nb);
  //Get the covariance matrix for statistical uncertainties on the measured distribution
  if (this->_dosys!=2) this->_cov = this->_svd->GetXtau();
  //Get the covariance matrix for statistical uncertainties on the response matrix
  //Uses Poisson or Gaussian-distributed toys, depending on response matrix histogram's Sumw2 setting.

  if (this->_dosys){
    add(this->_cov,this->_svd->GetAdetCovMatrix (this->_NToys));
  }
  this->_haveCov= true;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::GetWgt()
{
  // Get weight matrix
  if (this->_dosys) RooUnfoldT<Hist,Hist2D>::GetWgt();   // can't add sys errors to weight, so calculate weight from covariance
  if (!this->_svd) return;

  //Get the covariance matrix for statistical uncertainties on the measured distribution
  this->_wgt = this->_svd->GetXinv();
  
  this->_haveWgt= true;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::GetSettings(){
    this->_minparm=0;
    this->_maxparm= this->_meas ? nBins(this->_meas,X) : 0;
    this->_stepsizeparm=1;
    this->_defaultparm=this->_maxparm/2;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::Streamer (TBuffer &R__b)
{
  // Stream an object of class RooUnfoldSvdT.
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    RooUnfoldSvdT<Hist,Hist2D>::Class()->ReadBuffer  (R__b, this);
    TH1::AddDirectory (oldstat);
  } else {
    RooUnfoldSvdT<Hist,Hist2D>::Class()->WriteBuffer (R__b, this);
  }
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT()
  : RooUnfoldT<Hist,Hist2D>()
{
  // Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>& RooUnfoldSvdT<Hist,Hist2D>::operator= (const RooUnfoldSvdT<Hist,Hist2D>& rhs)
{
  // Assignment operator for copying RooUnfoldSvdT settings.
  Assign(rhs);
  return *this;
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::~RooUnfoldSvdT()
{
  Destroy();
}


template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::SetKterm (Int_t kreg)
{
  // Set regularisation parameter
  this->_kreg= kreg;
}


template<class Hist,class Hist2D> Int_t
RooUnfoldSvdT<Hist,Hist2D>::GetKterm() const
{
  // Return regularisation parameter
  return this->_kreg;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::SetRegParm (Double_t parm)
{
  // Set regularisation parameter
  SetKterm(Int_t(parm+0.5));
}

template<class Hist,class Hist2D> Double_t
RooUnfoldSvdT<Hist,Hist2D>::GetRegParm() const
{
  // Return regularisation parameter
  return GetKterm();
}

template class RooUnfoldSvdT<TH1,TH2>;
ClassImp (RooUnfoldSvd);

template class RooUnfoldSvdT<RooAbsReal,RooAbsReal>;
ClassImp (RooFitUnfoldSvd);

