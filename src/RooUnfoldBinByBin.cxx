//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Unfolding class using the bin by bin method of conversion factors. 
//
// Authors: Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p> Uses the correction factor method to unfold the distribution by looking at each bin individually.</p>
<p> This method cannot account for bin migration and as such cannot unfold reliably if a bias/smearing effects are applied.</p>
<p>Can only handle 1 dimensional distributions
<p>True and measured distributions must have the same binning
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldBinByBin.h"

#include <iostream>
#include "TH1.h"
#include "TH2.h"

#include "RooUnfoldResponse.h"

template<class Hist,class Hist2D> 
RooUnfoldBinByBinT<Hist,Hist2D>::RooUnfoldBinByBinT (const RooUnfoldBinByBinT<Hist,Hist2D>& rhs)
  : RooUnfoldT<Hist,Hist2D> (rhs)
{
  // Copy constructor.
  GetSettings();  
}

template<class Hist,class Hist2D>
RooUnfoldBinByBinT<Hist,Hist2D>::RooUnfoldBinByBinT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, 
                            const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  GetSettings();
}

template<class Hist,class Hist2D>
RooUnfoldBinByBinT<Hist,Hist2D>::~RooUnfoldBinByBinT()
{
}

template<class Hist,class Hist2D>
RooUnfoldBinByBinT<Hist,Hist2D>*
RooUnfoldBinByBinT<Hist,Hist2D>::Clone (const char* newname) const
{
    //Clones object
  RooUnfoldBinByBinT<Hist,Hist2D>* unfold= new RooUnfoldBinByBinT<Hist,Hist2D>(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}



template<class Hist,class Hist2D> TVectorD*
RooUnfoldBinByBinT<Hist,Hist2D>::Impl()
{
  return &_factors;
}

template<class Hist,class Hist2D> void
RooUnfoldBinByBinT<Hist,Hist2D>::Unfold()
{
    const TVectorD& vmeas=  this->Vmeasured();
    const TVectorD& vtrain= this->_res->Vmeasured();
    const TVectorD& vtruth= this->_res->Vtruth();

    TVectorD fakes= this->_res->Vfakes();
    Double_t fac= 0.0;
    if (this->_res->FakeEntries()) {
      fac= vtrain.Sum();
      if (fac!=0.0) fac= vmeas.Sum() / fac;
      if (this->_verbose>=1) std::cout << "Subtract " << fac*fakes.Sum() << " fakes from measured distribution" << std::endl;
    }

    this->_rec.ResizeTo(this->_nt);
    this->_factors.ResizeTo(this->_nt);
    Int_t nb= this->_nm < this->_nt ? this->_nm : this->_nt;
    for (int i=0; i<nb; i++) {
      Double_t train= vtrain[i]-fakes[i];
      if (train==0.0) continue;
      Double_t c= vtruth[i]/train;
      this->_factors[i]= c;
      this->_rec[i]= c * (vmeas[i]-fac*fakes[i]);
    }
    this->_unfolded= true;
}

template<class Hist,class Hist2D> void
RooUnfoldBinByBinT<Hist,Hist2D>::GetCov()
{
    const TMatrixD& covmeas= this->GetMeasuredCov();
    this->_cov.ResizeTo(this->_nt,this->_nt);
    Int_t nb= this->_nm < this->_nt ? this->_nm : this->_nt;
    for (int i=0; i<nb; i++)
      for (int j=0; j<nb; j++)
        this->_cov(i,j)= this->_factors[i]*this->_factors[j]*covmeas(i,j);
    this->_haveCov= true;
}

template<class Hist,class Hist2D> void
RooUnfoldBinByBinT<Hist,Hist2D>::GetSettings(){
    this->_minparm=0;
    this->_maxparm=0;
    this->_stepsizeparm=0;
    this->_defaultparm=0;
}


template<class Hist,class Hist2D> 
RooUnfoldBinByBinT<Hist,Hist2D>::RooUnfoldBinByBinT()
  : RooUnfoldT<Hist,Hist2D>()
{
  // Default constructor. Use Setup() to prepare for unfolding.
  GetSettings();
}

template<class Hist,class Hist2D> 
RooUnfoldBinByBinT<Hist,Hist2D>::RooUnfoldBinByBinT (const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  GetSettings();
}

template<class Hist,class Hist2D> 
RooUnfoldBinByBinT<Hist,Hist2D>::RooUnfoldBinByBinT (const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  GetSettings();
}

template<class Hist,class Hist2D> 
RooUnfoldBinByBinT<Hist,Hist2D>& RooUnfoldBinByBinT<Hist,Hist2D>::operator= (const RooUnfoldBinByBinT<Hist,Hist2D>& rhs)
{
  // Assignment operator for copying RooUnfoldBinByBinT settings.
  this->Assign(rhs);
  return *this;
}

template class RooUnfoldBinByBinT<TH1,TH2>;
ClassImp (RooUnfoldBinByBin);

template class RooUnfoldBinByBinT<RooAbsReal,RooAbsReal>;
ClassImp (RooFitUnfoldBinByBin);

