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

#include "RooUnfoldTH1Helpers.h"
#include "RooUnfoldResponse.h"

using namespace RooUnfolding;

using std::cout;
using std::cerr;
using std::endl;

ClassImp (RooUnfoldSvd);

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const RooUnfoldSvdT<Hist,Hist2D>& rhs)
  : RooUnfold (rhs)
{
  // Copy constructor.
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  // The regularisation parameter is kreg.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg, Int_t ntoyssvd,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _kreg(kreg ? kreg : res->GetNbinsTruth()/2)
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
  RooUnfold::Reset();
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
  RooUnfold::Assign (rhs);
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
  TH1::AddDirectory (kFALSE);
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
    TVectorD fakes= this->_res->Vfakes();
    Double_t fac= this->_res->Vmeasured().Sum();
    if (fac!=0.0) fac=  this->Vmeasured().Sum() / fac;
    if (this->_verbose>=1) cout << "Subtract " << fac*fakes.Sum() << " fakes from measured distribution" << endl;
    for (Int_t i= 1; i<=this->_nm; i++)
      this->_meas1d->SetBinContent (i, this->_meas1d->GetBinContent(i)-(fac*fakes[i-1]));
  }

  this->_meascov= new TH2D ("meascov", "meascov", this->_nb, 0.0, 1.0, this->_nb, 0.0, 1.0);
  const TMatrixD& cov= this->GetMeasuredCov();
  for (Int_t i= 0; i<this->_nm; i++)
    for (Int_t j= 0; j<this->_nm; j++)
      this->_meascov->SetBinContent (i+1, j+1, cov(i,j));

  if (this->_verbose>=1) cout << "SVD init " << this->_reshist->GetNbinsX() << " x " << this->_reshist->GetNbinsY()
                        << " bins, kreg=" << this->_kreg << endl;
  this->_svd= new SVDUnfold (this->_meas1d, this->_meascov, this->_train1d, this->_truth1d, this->_reshist);

  Hist* rechist= this->_svd->Unfold (this->_kreg);

  this->_rec.ResizeTo (this->_nt);
  for (Int_t i= 0; i<this->_nt; i++) {
    this->_rec[i]= rechist->GetBinContent(i+1);
  }

  if (this->_verbose>=2) {
    printTable (cout, this->_truth1d, this->_train1d, 0, this->_meas1d, rechist, this->_nb, this->_nb, kFALSE, kErrors);
    TMatrixD* resmat= RooUnfoldResponseT<Hist,Hist2D>::H2M (this->_reshist, this->_nb, this->_nb);
    RooUnfoldResponseT<Hist,Hist2D>::PrintMatrix(*resmat,"SVDUnfold response matrix");
    delete resmat;
  }

  delete rechist;
  TH1::AddDirectory (oldstat);

  this->_unfolded= true;
  this->_haveCov=  false;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::GetCov()
{
  if (!this->_svd) return;
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);

  TH2 *unfoldedCov= 0, *adetCov= 0;
  //Get the covariance matrix for statistical uncertainties on the measured distribution
  if (this->_dosys!=2) unfoldedCov= this->_svd->GetXtau();
  //Get the covariance matrix for statistical uncertainties on the response matrix
  //Uses Poisson or Gaussian-distributed toys, depending on response matrix histogram's Sumw2 setting.
  if (this->_dosys)        adetCov= this->_svd->GetAdetCovMatrix (this->_NToys);

  this->_cov.ResizeTo (this->_nt, this->_nt);
  for (Int_t i= 0; i<this->_nt; i++) {
    for (Int_t j= 0; j<this->_nt; j++) {
      Double_t v  = 0;
      if (unfoldedCov) v  = unfoldedCov->GetBinContent(i+1,j+1);
      if (adetCov)     v += adetCov    ->GetBinContent(i+1,j+1);
      this->_cov(i,j)= v;
    }
  }

  delete adetCov;
  TH1::AddDirectory (oldstat);

  this->_haveCov= true;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::GetWgt()
{
  // Get weight matrix
  if (this->_dosys) RooUnfold::GetWgt();   // can't add sys errors to weight, so calculate weight from covariance
  if (!this->_svd) return;
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);

  //Get the covariance matrix for statistical uncertainties on the measured distribution
  TH2* unfoldedWgt= this->_svd->GetXinv();

  this->_wgt.ResizeTo (this->_nt, this->_nt);
  for (Int_t i= 0; i<this->_nt; i++) {
    for (Int_t j= 0; j<this->_nt; j++) {
      this->_wgt(i,j)= unfoldedWgt->GetBinContent(i+1,j+1);
    }
  }

  TH1::AddDirectory (oldstat);

  this->_haveWgt= true;
}

template<class Hist,class Hist2D> void
RooUnfoldSvdT<Hist,Hist2D>::GetSettings(){
    this->_minparm=0;
    this->_maxparm= this->_meas ? this->_meas->GetNbinsX() : 0;
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
  : RooUnfold()
{
  // Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const char* name, const char* title)
  : RooUnfold(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::RooUnfoldSvdT (const TString& name, const TString& title)
  : RooUnfold(name,title)
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

