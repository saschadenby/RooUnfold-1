/*! \class RooUnfoldPoissonT
*/

//#define OLDERRS   // restore old (incorrect) error calculation
//#define OLDERRS2  // restore old (incorrect) systematic error calculation
//#define OLDMULT   // restore old (slower) matrix multiplications

#include "RooUnfoldPoisson.h"
#include "RooUnfoldTH1Helpers.h"
#ifndef NOROOFIT
#include "RooUnfoldFitHelpers.h"
#endif

#include <iostream>
#include <iomanip>
#include <math.h>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "RooUnfoldHelpers.h"
#include "RooUnfoldResponse.h"

using namespace RooUnfolding;

template<class Hist,class Hist2D>
RooUnfoldPoissonT<Hist,Hist2D>::RooUnfoldPoissonT (const RooUnfoldPoissonT<Hist,Hist2D>& rhs)
  : RooUnfoldT<Hist,Hist2D> (rhs)
{
  //! Copy constructor.
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D>
RooUnfoldPoissonT<Hist,Hist2D>::RooUnfoldPoissonT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Double_t regparm,
                                const char* name, const char* title)
  : _regparm(regparm), RooUnfoldT<Hist,Hist2D> (res, meas, name, title)
{

  //! Constructor with response matrix object and measured unfolding input histogram.
  //! The regularisation parameter is niter (number of iterations).
  Init();
}

template<class Hist,class Hist2D> void
RooUnfoldPoissonT<Hist,Hist2D>::Init()
{
  GetSettings();
}

template<class Hist,class Hist2D> void 
RooUnfoldPoissonT<Hist,Hist2D>::Reset()
{
  Init();
  RooUnfoldT<Hist,Hist2D>::Reset();
}

template<class Hist,class Hist2D> void 
RooUnfoldPoissonT<Hist,Hist2D>::Assign (const RooUnfoldPoissonT<Hist,Hist2D>& rhs)
{
  RooUnfoldT<Hist,Hist2D>::Assign (rhs);
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldPoissonT<Hist,Hist2D>::CopyData (const RooUnfoldPoissonT<Hist,Hist2D>& rhs)
{
  this->_regparm=    rhs._regparm;
}

template<class Hist,class Hist2D> void
RooUnfoldPoissonT<Hist,Hist2D>::Unfold() const
{
  
  
  this->setup();

  // Set the start values of the truth bins according to some
  // passed truth histogram -> User should define this! Default response truth.

  // Minimize the regularized nllh.
  MinimizeRegLLH();

  this->_cache._rec.ResizeTo(this->_nt);  // drop fakes in final bin
  this->_cache._unfolded= true;
  this->_cache._haveCov=  false;
}

template<class Hist,class Hist2D> void
RooUnfoldPoissonT<Hist,Hist2D>::GetCov() const
{
    // Get covariance.

  this->_cache._cov.ResizeTo (this->_nt, this->_nt);

  this->_cache._haveCov= true;
}

template<class Hist,class Hist2D> void
RooUnfoldPoissonT<Hist,Hist2D>::GetSettings() const
{
 
  this->_cache._minparm=1;
  this->_cache._maxparm=2;
  this->_cache._stepsizeparm=1e-2;
  this->_cache._defaultparm=2;
}

//-------------------------------------------------------------------------
template<class Hist,class Hist2D> void
RooUnfoldPoissonT<Hist,Hist2D>::setup() const
{
  this->_response.ResizeTo(this->_nt,this->_nt);
  this->_response = this->_res->Mresponse(true);
  this->_data.ResizeTo(this->_nm);
  this->_data = this->Vmeasured();
  this->_truth_start.ResizeTo(this->_nt);
  this->_truth_start = this->_res->Vtruth();
}

template<class Hist,class Hist2D> void
RooUnfoldPoissonT<Hist,Hist2D>::Print(Option_t* option) const
{
  
  // Fill this.
}

template<class Hist,class Hist2D> double*
RooUnfoldPoissonT<Hist,Hist2D>::Rmu(const double* truth) const
{
 
 double* Rmu = new double[_response.GetNrows()];

  for (int i = 0; i < _response.GetNrows(); i++){
    double reco_bin = 0;
    for (int j = 0; j < _response.GetNcols(); j++){
      reco_bin += _response[i][j]*truth[j];
    }
    Rmu[i] = reco_bin;
  }

  return Rmu;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldPoissonT<Hist,Hist2D>::NegativeLLH(const double* truth) const
{

  double* nu = Rmu(truth);

  Double_t func_val = 0;

  for (int i = 0; i < _response.GetNrows(); i++){
    func_val += nu[i] - _data[i] * log(nu[i]);
  }

  delete[] nu;

  return func_val;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldPoissonT<Hist,Hist2D>::TikinovReg(const double* truth) const
{

  Double_t func_val = 0;

  for (int i = 0; i < _response.GetNcols() - 2; i++){
    Double_t der = (truth[i+2] - truth[i+1]) - (truth[i+1] - truth[i]);
    func_val += der*der;
  }

  return func_val;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldPoissonT<Hist,Hist2D>::RegLLH(const double* truth) const
{
  return NegativeLLH(truth) + this->_regparm*TikinovReg(truth);
}


template<class Hist,class Hist2D> void
RooUnfoldPoissonT<Hist,Hist2D>::MinimizeRegLLH() const
{
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  min->SetMaxFunctionCalls(10000000);
  min->SetTolerance(0.001);
  min->SetStrategy(2);
  min->SetPrintLevel(1);


  ROOT::Math::Functor f(this, &RooUnfoldPoissonT<Hist,Hist2D>::RegLLH,_response.GetNcols());

  min->SetFunction(f);

  double* step = new double[_response.GetNcols()];
  double* start = new double[_response.GetNcols()];

  for (int i = 0; i < _response.GetNcols(); i++){
    step[i] = 1;
    start[i] = _truth_start[i];

    std::string s = std::to_string(i);
    std::string x("mu");
    x.append(s);

    min->SetVariable(i,x.c_str(),start[i], step[i]);
  }

  // do the minimization
  min->Minimize();

  delete[] step;
  delete[] start;
  delete min;
}



template<class Hist,class Hist2D>
RooUnfoldPoissonT<Hist,Hist2D>::RooUnfoldPoissonT()
  : RooUnfoldT<Hist,Hist2D>()
{

  //! Default constructor. Use Setup() to prepare for unfolding.]
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldPoissonT<Hist,Hist2D>::RooUnfoldPoissonT (const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldPoissonT<Hist,Hist2D>::RooUnfoldPoissonT (const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>RooUnfoldPoissonT<Hist,Hist2D>& 
RooUnfoldPoissonT<Hist,Hist2D>::operator= (const RooUnfoldPoissonT<Hist,Hist2D>& rhs)
{
  //! Assignment operator for copying RooUnfoldPoisson settings.
  Assign(rhs);
  return *this;
}

template<class Hist,class Hist2D>
void  RooUnfoldPoissonT<Hist,Hist2D>::SetRegParm (Double_t parm)
{
  //! Set regularisation parameter (number of iterations)
  this->_regparm = parm;
}

template<class Hist,class Hist2D>
double RooUnfoldPoissonT<Hist,Hist2D>::GetRegParm() const
{
  //! Return regularisation parameter (number of iterations)
  return this->_regparm;
}

template class RooUnfoldPoissonT<TH1,TH2>;
ClassImp (RooUnfoldPoisson)

#ifndef NOROOFIT
template class RooUnfoldPoissonT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
ClassImp (RooFitUnfoldPoisson)
#endif
