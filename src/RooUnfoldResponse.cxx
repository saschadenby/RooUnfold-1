//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Response Matrix
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
 <p> Class to create response object as used in RooUnfold </p>
 <p> Contains measured and truth distributions as TH1s and the response matrix as a TH2. Also contains methods for handling these data</p>
<p> Can handle 1,2 or 3 dimensional histograms and return vectors and matrices of their bin content and error (1 and 2D distributions respectively).
 Conversely can also convert these vectors and matrices into TH1s and TH2Ds. </p>
<p> Can also take a variety of parameters as inputs. This includes maximum and minimum values, distributions and vectors/matrices of values. </p>
<p> This class does the numerical modifications needed to allow unfolding techniques to work in the unfolding routines used in RooUnfold. </p>
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldResponse.h"
#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"
#ifndef NOROOFIT
#include "RooUnfoldFitHelpers.h"
#include "RooDataHist.h"
#endif

#include <iostream>
#include <assert.h>
#include <cmath>

#include "TClass.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TPRegexp.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TCollection.h"

#include "TH1.h"
#include "TH2.h"


using namespace RooUnfolding;

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,18,0)
#define HAVE_RooUnfoldFoldingFunction
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::pow;
using std::sqrt;


#ifdef HAVE_RooUnfoldFoldingFunction
template<class Hist, class Hist2D>
class RooUnfoldFoldingFunction {
public:
  RooUnfoldFoldingFunction<Hist,Hist2D> (const RooUnfoldResponseT<Hist,Hist2D>* res, TF1* func, Double_t eps=1e-12, bool verbose=false)
    : _res(res), _func(func), _eps(eps), _verbose(verbose), _fvals(_res->GetNbinsMeasured()) {
    _ndim= dynamic_cast<TF3*>(_func) ? 3 :
      dynamic_cast<TF2*>(_func) ? 2 : 1;
    if (_ndim>=2 && eps==1e-12) eps= 0.000001;
    FVals();
  }

  double operator() (double* x, double* p) const {
    const Hist* mes= _res->Hmeasured();
    Int_t bin;
    if      (_ndim==1) bin= findBin (mes, x[0]);
    else if (_ndim==2) bin= findBin (mes, x[0], x[1]);
    else               bin= findBin (mes, x[0], x[1], x[2]);
    if (bin<0 || bin>=_res->GetNbinsMeasured()) return 0.0;
    for (Int_t i=0, n=_func->GetNpar(); i<n; i++) {
      if (p[i] == _func->GetParameter(i)) continue;
      _func->SetParameters(p);
      FVals();
      break;
    }
    Double_t fy= _fvals[bin];
    if (_verbose) cout << "x=" << x[0] << ", bin=" << bin << " -> " << fy << endl;
    return fy;
  }

private:
  void FVals() const {
    const Hist* tru= _res->Htruth();
    if (_verbose) {
      cout << "p=";
      for (int i=0, n=_func->GetNpar(); i<n; i++) cout <<_func->GetParameter(i)<<",";
      cout << " f=";
    }
    _fvals.Zero();
    for (Int_t i=0, n=_res->GetNbinsTruth(); i<n; i++) {
      Int_t jx, jy, jz;
      binXYZ(tru,i,jx,jy,jz);
      Double_t fv;
      if (_eps<=0.0) {
        if (_ndim>=2)
          fv= _func->Eval (binCenter(tru,jx,RooUnfolding::X),
                           binCenter(tru,jy,RooUnfolding::Y),
                           binCenter(tru,jz,RooUnfolding::Z));
        else
          fv= _func->Eval (binCenter(tru,jx,RooUnfolding::X));
      } else {
        if        (_ndim==1) {
          Double_t tw= binWidth(tru,jx,RooUnfolding::X), tlo= binLowEdge(tru,jx,RooUnfolding::X), thi= tlo+tw;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
          fv= _func->Integral (tlo, thi, _eps) / tw;
        } else {
          Double_t tlo[3]= { binLowEdge(tru,jx,RooUnfolding::X), binLowEdge(tru,jy,RooUnfolding::Y), binLowEdge(tru,jz,RooUnfolding::Z) };
          Double_t thi[3]= { binHighEdge(tru,jx,RooUnfolding::X), binHighEdge(tru,jy,RooUnfolding::Y), binHighEdge(tru,jz,RooUnfolding::Z) };
          Double_t relerr=0.0;
          fv= _func->IntegralMultiple (_ndim, tlo, thi, _eps, relerr);
          fv /= binWidth(tru,jx,RooUnfolding::X) * binWidth(tru,jy,RooUnfolding::Y);
          if (_ndim>=3) fv /= binWidth(tru,jz,RooUnfolding::Z);
#else
          fv= _func->Integral (tlo, thi, (Double_t*)0, _eps) / tw;
        } else if (_ndim==2) {
          fv= _func->Integral (binLowEdge(tru,jx,RooUnfolding::X), binHighEdge(tru,jx,RooUnfolding::X),
                               binLowEdge(tru,jy,RooUnfolding::Y), binHighEdge(tru,jy,RooUnfolding::Y),
                               _eps);
          fv /=binWidth(tru,jx,RooUnfolding::X) * binWidth(tru,jy,RooUnfolding::Y);
        } else {
          fv= _func->Integral (binLowEdge(tru,jx,RooUnfolding::X), binHighEdge(tru,jx,RooUnfolding::X),
                               binLowEdge(tru,jy,RooUnfolding::Y), binHighEdge(tru,jy,RooUnfolding::Y),
                               binLowEdge(tru,jz,RooUnfolding::Z), binHighEdge(tru,jz,RooUnfolding::Z),
                               _eps);
          fv /=binWidth(tru,jx,RooUnfolding::X) * binWidth(tru,jy,RooUnfolding::Y) * binWidth(tru,jz,RooUnfolding::Z);          
#endif
        }
      }
      if (_verbose) cout << " " << fv;
      for (Int_t bin=0, m=_res->GetNbinsMeasured(); bin<m; bin++) {
        _fvals[bin] += fv * (*_res)(bin,i);
      }
    }
    if (_verbose) cout << endl;
  }

  const RooUnfoldResponseT<Hist,Hist2D>* _res;
  TF1* _func;
  Double_t _eps;
  bool _verbose;
  mutable TVectorD _fvals;
  Int_t _ndim;
};
#endif  

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::setup()
{
  _tru= _mes= _fak= 0;
  _res= 0;
  _nm= _nt= _mdim= _tdim= 0;
  this->ClearCache();
  SetNameTitleDefault ("response", "Response");
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::ClearCache()
{
  if(_vMes){ delete _vMes; _vMes= 0; }
  if(_eMes){ delete _eMes; _eMes= 0; }
  if(_vFak){ delete _vFak; _vFak= 0; }
  if(_vTru){ delete _vTru; _vTru= 0; }
  if(_eTru){ delete _eTru; _eTru= 0; }
  if(_mRes){ delete _mRes; _mRes= 0; }
  if(_eRes){ delete _eRes; _eRes= 0; }
  _cached= false;
}

template <class Hist, class Hist2D> Hist2D*
RooUnfoldResponseT<Hist,Hist2D>::HresponseNoOverflow() const
{
  const Hist2D* res = Hresponse();
  TVectorD vals(h2v<Hist>(res,_overflow,_density));
  TVectorD errs(h2ve<Hist>(res,_overflow,_density));  
  return createHist<Hist2D>(vals,errs,name(res),title(res),vars(res),_overflow);
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::Print (Option_t* /* option */) const
{
  std::cout << "RooUnfoldResponseT @ " << this << std::endl;
  printHistogram(this->Hmeasured());
  printHistogram(this->Hfakes());
  printHistogram(this->Htruth());
  printHistogram(this->Hresponse());
}

template <class Hist, class Hist2D> bool
RooUnfoldResponseT<Hist,Hist2D>::Cached () const
{
  return _cached;
}

template <class Hist, class Hist2D> Hist*
RooUnfoldResponseT<Hist,Hist2D>::ApplyToTruth (const Hist* truth, const char* name) const
{
  // Apply the response matrix to the truth
  // Errors not set, since we assume original truth has no errors
  if (!Htruth()) return 0;  // Needed for checking binning if nothing else

  // If no truth histogram input, use training truth
  // If truth histogram input, make sure its binning is correct
  TVectorD resultvect;;
  if (truth) {
    if (nBins(truth,RooUnfolding::X) != nBins(_tru,RooUnfolding::X) ||
        nBins(truth,RooUnfolding::Y) != nBins(_tru,RooUnfolding::Y) ||
        nBins(truth,RooUnfolding::Z) != nBins(_tru,RooUnfolding::Z))
      cerr << "Warning: RooUnfoldResponseT<Hist,Hist2D>::ApplyToTruth truth histogram is a different size ("
           << (nBins(truth,RooUnfolding::X) * nBins(truth,RooUnfolding::Y) * nBins(truth,RooUnfolding::Z)) << " bins) or shape from response matrix truth ("
           << ( nBins(_tru,RooUnfolding::X) * nBins( _tru,RooUnfolding::Y) * nBins( _tru,RooUnfolding::Z)) << " bins)" << endl;
    resultvect= h2v (truth, _overflow,_density);
  } else {
    resultvect= Vtruth();
  }

  resultvect *= Mresponse();   // v= A*v

  // Turn results vector into properly binned histogram
  const Hist* t = Hmeasured();
  Hist* result= createHist<Hist>(resultvect, RooUnfolding::name(t),name, vars(t), _overflow);
  return result;
}


template <class Hist, class Hist2D> TF1*
RooUnfoldResponseT<Hist,Hist2D>::MakeFoldingFunction (TF1* func, Double_t eps, Bool_t verbose) const
{
  // Creates a function object that applies the response matrix to a user parametric function.
  // This can be fitted to the measured distribution as an alternative to unfolding.
  // The returned object is owned by the caller. The function will be binned.
  // Specify eps=0 to calculate function at bin centers; otherwise integrates over each bin (may be slow).
  // Example:
  //    TF1* func= new TF1 ("func", "gaus", 0, 10);
  //    TF1* fold= respose->MakeFoldingFunction(func);
  //    histMeasured->Fit(fold);
  //    fold->Draw("h"); // draw function fitted to histMeasured
  //    func->Draw();    // draw truth function
#ifdef HAVE_RooUnfoldFoldingFunction
  Int_t np= func->GetNpar();
  RooUnfoldFoldingFunction<Hist,Hist2D> ff (this, func, eps, verbose);
  TString name= func->GetName();
  name += "_folded";
  TF1* f;
  if        (TF3* func3= dynamic_cast<TF3*>(func))
    f= new TF3 (name, ROOT::Math::ParamFunctor(ff),
                func3->GetXmin(), func3->GetXmax(),
                func3->GetYmin(), func3->GetYmax(),
                func3->GetZmin(), func3->GetZmax(), np);
  else if (TF2* func2= dynamic_cast<TF2*>(func))
    f= new TF2 (name, ROOT::Math::ParamFunctor(ff),
                func2->GetXmin(), func2->GetXmax(),
                func2->GetYmin(), func2->GetYmax(), np);
  else
    f= new TF1 (name, ROOT::Math::ParamFunctor(ff),
                func ->GetXmin(), func ->GetXmax(), np);
  f->SetNpx (_nm<=2 ? 4 : _nm==3 ? 6 : _nm);  // TF1 requires Npx>=4
  // Copy parameters in case we set them in func
  f->SetParameters (func->GetParameters());
  f->SetParErrors  (func->GetParErrors());
  for (Int_t i=0; i<np; i++) {
    Double_t plo=0.0, phi=0.0;
    func->GetParLimits (i, plo, phi);
    f   ->SetParLimits (i, plo, phi);
    f->SetParName (i, func->GetParName(i));
  }
  return f;
#else
  cerr << "RooUnfoldResponseT<Hist,Hist2D>::MakeFoldingFunction not supported in this version of ROOT" << endl;
  return 0;
#endif
}


template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>*
RooUnfoldResponseT<Hist,Hist2D>::RunToy() const
{
  // Returns new RooUnfoldResponse object with smeared response matrix elements for use as a toy.
  const Hist2D* hres= this->Hresponse();
  TMatrixD values(h2m(hres));
  TMatrixD errors(h2me(hres));
  for (Int_t i= 0; i<_nm; i++) {
    for (Int_t j= 0; j<_nt; j++) {
      Double_t e= errors(i,j);
      if (e>0.0) {
        Double_t v= values(i,j) + gRandom->Gaus(0.0,e);
        if (v<0.0) v= 0.0;
        values(i,j) = v;
      }
    }
  }

  Hist2D* smeared = createHist<Hist2D>(values,errors,RooUnfolding::name(hres),RooUnfolding::title(hres),var(hres,X),var(hres,Y));
  TString name= GetName();
  name += "_toy";
  RooUnfoldResponseT<Hist,Hist2D>* res= new RooUnfoldResponseT<Hist,Hist2D> (name.Data(),this->GetTitle());

  res->_mes = this->_mes;
  res->_tru = this->_tru;
  res->_res = smeared;

  return res;
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::SetNameTitleDefault (const char* defname, const char* deftitle)
{
  // Set object name and title
  const char* s= GetName();
  if (s[0] == '\0') {
    if (_res) s= name(_res);
    if (s[0] == '\0') {
      if (defname) SetName (defname);
      else if (_mes && _tru) {
        TString n= name(_mes);
        if (n.Length()) n.Append ("_");
        n.Append (name(_tru));
        if (!n.Length()) n= "response";
        SetName (n);
      }
    } else
      SetName (s);
  }
  s= GetTitle();
  if (s[0] == '\0') {
    if (_res) s= title(_res);
    if (s[0] == '\0') {
      if (deftitle) SetTitle (deftitle);
      else if (_mes && _tru) {
        TString n= title(_tru);
        if (n.Length()) n.Append (" #rightarrow ");
        n.Append (title(_mes));
        if (n.Length())
          n.Prepend ("Response ");
        else
          n= "Response";
        SetTitle (n);
      }
    } else
      SetTitle (s);
  }
}

template<class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::Streamer (TBuffer &R__b)
{
  if (R__b.IsReading()) {
    RooUnfoldResponseT<Hist,Hist2D>::Class()->ReadBuffer  (R__b, this);
  } else {
    RooUnfoldResponseT<Hist,Hist2D>::Class()->WriteBuffer (R__b, this);
  }
}

template <> void
RooUnfoldResponseT<TH1,TH2>::Streamer (TBuffer &R__b)
{
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    RooUnfoldResponseT<TH1,TH2>::Class()->ReadBuffer  (R__b, this);
    TH1::AddDirectory (oldstat);
  } else {
    RooUnfoldResponseT<TH1,TH2>::Class()->WriteBuffer (R__b, this);
  }
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT()
  : TNamed()
{
  // RooUnfoldResponseT<Hist,Hist2D> default constructor. Use Setup() to set values.
  this->setup();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const char*    name, const char*    title)
  : TNamed(name,title)
{
  // RooUnfoldResponseT<Hist,Hist2D> default named constructor. Use Setup() to set values.
  this->setup();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const TString& name, const TString& title)
  : TNamed(name,title)
{
  // RooUnfoldResponseT<Hist,Hist2D> default named constructor. Use Setup() to set values.
  this->setup();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const RooUnfoldResponseT<Hist,Hist2D>& rhs) : 
  TNamed(rhs.GetName(),rhs.GetTitle()),
  _mdim(rhs._mdim),     
  _tdim(rhs._tdim),     
  _nm(rhs._nm),       
  _nt(rhs._nt),       
  _mes(clone(rhs._mes)),      
  _fak(clone(rhs._fak)),      
  _tru(clone(rhs._tru)),      
  _res(clone(rhs._res)),      
  _overflow(rhs._overflow),
  _density(rhs._density)
{
  // copy constructor
  this->ClearCache();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const char* name, const char* title, Hist2D* response, Hist* truth, Hist* reco, bool overflow,bool density) :
  TNamed(name,title),
  _mdim(dim(reco)),
  _tdim(dim(truth)),
  _nm(nBins(reco)),
  _nt(nBins(truth)),
  _mes(reco),
  _fak(0),
  _tru(truth),
  _res(response),
  _overflow(overflow),
  _density(density)
{
  // explicit constructor
  this->ClearCache();
}



template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::~RooUnfoldResponseT()
{
  // RooUnfoldResponseT<Hist,Hist2D> destructor
  ClearCache();
}


template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::GetDimensionMeasured() const
{
  // Dimensionality of the measured distribution (1=1D, 2=2D, 3=3D)
  return _mdim;
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::GetDimensionTruth() const
{
  // Dimensionality of the truth distribution (1=1D, 2=2D, 3=3D)
  return _tdim;
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::GetNbinsMeasured() const
{
  // Total number of bins in the measured distribution
  return _nm;
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::GetNbinsTruth() const
{
  // Total number of bins in the truth distribution
  return _nt;
}


template<class Hist, class Hist2D>
const Hist* RooUnfoldResponseT<Hist,Hist2D>::Hmeasured() const
{
  // Measured distribution, including fakes
  return _mes;
}


template<class Hist, class Hist2D>
Hist*         RooUnfoldResponseT<Hist,Hist2D>::Hmeasured()
{
  return _mes;
}


template<class Hist, class Hist2D>
const Hist* RooUnfoldResponseT<Hist,Hist2D>::Hfakes() const
{
  // Fakes distribution
  return _fak;
}


template<class Hist, class Hist2D>
Hist*         RooUnfoldResponseT<Hist,Hist2D>::Hfakes()
{
  return _fak;
}

template<class Hist, class Hist2D>
const Hist*   RooUnfoldResponseT<Hist,Hist2D>::Htruth() const
{
  // Truth distribution, used for normalisation
  return _tru;
}

template<class Hist, class Hist2D>
Hist*         RooUnfoldResponseT<Hist,Hist2D>::Htruth()
{
  return _tru;
}

template<class Hist, class Hist2D>
const Hist2D*   RooUnfoldResponseT<Hist,Hist2D>::Hresponse() const
{
  // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  return _res;
}

template<class Hist, class Hist2D>
Hist2D*         RooUnfoldResponseT<Hist,Hist2D>::Hresponse()
{
  return _res;
}


template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vmeasured() const
{
  // Measured distribution as a TVectorD
  if (!_vMes) _cached= (_vMes= new TVectorD(h2v  (_mes, _overflow,_density)));
  if(_vMes->GetNrows() != this->_nm) throw std::runtime_error("invalid dimensionality in measured vector!");
  return *_vMes;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vfakes() const
{
  // Fakes distribution as a TVectorD
  if(!_fak){
    _cached = (_vFak = new TVectorD(this->_nm) );
  } else {
    if (!_vFak) _cached= (_vFak= new TVectorD(h2v  (_fak, _overflow,_density)));
  }
  if(_vFak->GetNrows() != this->_nm) throw std::runtime_error("invalid dimensionality in fakes vector!");
  return *_vFak;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Emeasured() const
{
  // Measured distribution errors as a TVectorD
  if (!_eMes) _cached= (_eMes= new TVectorD(h2ve (_mes, _overflow,_density)));
  if(_eMes->GetNrows() != this->_nm) throw std::runtime_error("invalid dimensionality in measured uncertainty vector!");
  return *_eMes;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vtruth() const
{
  // Truth distribution as a TVectorD
  if (!_vTru) _cached= (_vTru= new TVectorD(h2v  (_tru, _overflow,_density))); 
  if(_vTru->GetNrows() != this->_nt) throw std::runtime_error("invalid dimensionality in truth vector!");
  return *_vTru;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Etruth() const
{
  // Truth distribution errors as a TVectorD
  if (!_eTru) _cached= (_eTru= new TVectorD(h2ve (_tru, _overflow,_density))); 
  if(_eTru->GetNrows() != this->_nt) throw std::runtime_error("invalid dimensionality in truth uncertainty vector!");
  return *_eTru;
}

template<class Hist, class Hist2D>
const TMatrixD& RooUnfoldResponseT<Hist,Hist2D>::Mresponse(bool norm) const
{
  // Response matrix as a TMatrixD: (row,column)=(measured,truth)

  if (!_mRes){
    if (norm){
      _cached= (_mRes= new TMatrixD(h2mNorm  (_res, _tru, _overflow,_density))); 
    } else {
      _cached = (_mRes = new TMatrixD(h2m (_res, _overflow, _density)));
    }
  }
  return *_mRes;
}

template<class Hist, class Hist2D>
const TMatrixD& RooUnfoldResponseT<Hist,Hist2D>::Eresponse(bool norm) const
{
  // Response matrix errors as a TMatrixD: (row,column)=(measured,truth)
  if (!_eRes){
    if (norm) {
      _cached= (_eRes= new TMatrixD(h2meNorm (_res, _tru, _overflow,_density))); 
    } else {
      _cached= (_eRes= new TMatrixD(h2me (_res, _overflow,_density))); 
    }
  }
  return *_eRes;
}


template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::operator() (Int_t r, Int_t t) const
{
  // Response matrix element (measured,truth)
  return Mresponse()(r,t);
}

template<class Hist, class Hist2D>
void RooUnfoldResponseT<Hist,Hist2D>::UseOverflow (Bool_t set)
{
  // Specify to use overflow bins. Only supported for 1D truth and measured distributions.
  _overflow= (set ? 1 : 0);
}

template<class Hist, class Hist2D>
Bool_t RooUnfoldResponseT<Hist,Hist2D>::UseOverflowStatus() const
{
  // Get UseOverflow setting
  return _overflow;
}

template<class Hist, class Hist2D>
Bool_t RooUnfoldResponseT<Hist,Hist2D>::UseDensityStatus() const
{
  // Get UseDensity setting
  return _density;
}

template<class Hist, class Hist2D>
bool RooUnfoldResponseT<Hist,Hist2D>::HasFakes() const
{
  // Return number of fake entries
  return _fak && !empty(_fak);
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

template class RooUnfoldResponseT<TH1,TH2>;

namespace {
  template<class Hist> int entries(const Hist* hist);

  template<class Hist> int entries(const Hist* hist){
    return hist->GetEntries();
  }

  template int entries<TH1>(TH1 const*);
  template int entries<TH2>(TH2 const*);

  template<class Hist, class Hist2D> void projectY(Hist2D* _res, Hist* _tru, bool overflow);
  template<class Hist, class Hist2D> void projectX(Hist2D* _res, Hist* _mes, bool overflow);
  template<class Hist, class Hist2D> void subtractProjectX(Hist2D* _res, Hist* _mes, Hist* _fak, bool overflow);

  template<> void projectY<TH1>(TH2* _res, TH1* _tru, bool overflow){
    Int_t s= _res->GetSumw2N();
    for (Int_t j= 1-overflow; j<_res->GetNbinsY()+1+overflow; j++) {
      Double_t ntru= 0.0, wtru= 0.0;
      for (Int_t i= 0; i<_res->GetNbinsX()+2; i++) {
        ntru +=      _res->GetBinContent (i, j);
        if (s) wtru += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t b= bin<TH1>(_tru, j, overflow);
      _tru->SetBinContent (b,      ntru);
      if (s) _tru->SetBinError   (b, sqrt(wtru));
    }
  }
  template<> void projectX<TH1>(TH2* _res, TH1* _mes, bool overflow){
    Int_t s= _res->GetSumw2N();
    for (Int_t i= 1-overflow; i<_res->GetNbinsX()+1+overflow; i++) {
      Double_t nmes= 0.0, wmes= 0.0;
      for (Int_t j= 0; j<_res->GetNbinsY(); j++) {
        nmes +=      _res->GetBinContent (i, j);
        if (s) wmes += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t b= RooUnfolding::bin (_mes, i, overflow);
      _mes->SetBinContent (b,      nmes );
      if (s) _mes->SetBinError   (b, sqrt(wmes));
    }
  }
  template<> void subtractProjectX<TH1>(TH2* _res, TH1* _mes, TH1* _fak, bool overflow){
    Int_t s= _res->GetSumw2N();
    Int_t sm= _mes->GetSumw2N(), nfake=0;
    for (Int_t i= 1-overflow; i<_res->GetNbinsX()+1+overflow; i++) {
      Double_t nmes= 0.0, wmes= 0.0;
      for (Int_t j= 0; j<_res->GetNbinsY()+2; j++) {
        nmes +=      _res->GetBinContent (i, j);
        if (s) wmes += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t b= RooUnfolding::bin (_mes, i, overflow);
      Double_t fake= _mes->GetBinContent (b) - nmes;
      if (fake!=0.0) nfake++;
      if (!s) wmes= nmes;
      _fak->SetBinContent (b, fake);
      _fak->SetBinError   (b, sqrt (wmes + (sm ? pow(_mes->GetBinError(b),2) : _mes->GetBinContent(b))));
    }
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,13,0)
    _fak->SetEntries (_fak->GetEffectiveEntries());  // 0 entries if 0 fakes
#else
    _fak->SetEntries (nfake);  // 0 entries if 0 fakes
#endif    
  }


  template<class Hist> int fill(Hist* hist, double x, double w);
  template<class Hist> int fill(Hist* hist, double x, double y, double w);
  template<class Hist> int fill(Hist* hist, double x, double y, double z, double w);
  template<> int fill<TH1>(TH1* hist, double x, double w){
    return hist->Fill (x, w);
  }
  template<> int fill<TH1>(TH1* hist, double x, double y, double w){
    return ((TH2*)hist)->Fill (x, y, w);
  }
  template<> int fill<TH1>(TH1* hist, double x, double y, double z, double w){
    return ((TH3*)hist)->Fill (x, y, z, w);
  }    
  template<> int fill<TH2>(TH2* hist, double x, double y, double w){
    return hist->Fill (x, y, w);
  } 


}




Int_t RooUnfoldResponse::Miss (Double_t xt)
{
  // Fill missed event into 1D Response Matrix
  return Miss1D(xt);
}

Int_t RooUnfoldResponse::Miss (Double_t xt, Double_t w)
{
  // Fill missed event into 1D (with weight) or 2D Response Matrix
  return _tdim==2 ? Miss2D(xt,w) : Miss1D(xt,w);
}

Int_t RooUnfoldResponse::Miss (Double_t xt, Double_t yt, Double_t w)
{
  // Fill missed event into 2D (with weight) or 3D Response Matrix
  return _tdim==3 ? Miss(xt,yt,w,1.0) : Miss2D(xt,yt,w);
}


Int_t RooUnfoldResponse::Fake (Double_t xr)
{
  // Fill fake event into 1D Response Matrix
  return Fake1D(xr);
}

Int_t RooUnfoldResponse::Fake (Double_t xr, Double_t w)
{
  // Fill fake event into 1D (with weight) or 2D Response Matrix
  return _mdim==2 ? Fake2D(xr,w) : Fake1D(xr,w);
}

Int_t RooUnfoldResponse::Fake (Double_t xr, Double_t yr, Double_t w)
{
  // Fill fake event into 2D (with weight) or 3D Response Matrix
  return _mdim==3 ? Fake(xr,yr,w,1.0) : Fake2D(xr,yr,w);
}




RooUnfoldResponse::RooUnfoldResponse(const RooUnfoldResponse& rhs)
  : RooUnfoldResponseT(rhs.GetName(), rhs.GetTitle())
{
  // RooUnfoldResponseT<class Hist, class Hist2D> copy constructor
  Setup(rhs);
}

RooUnfoldResponse&
RooUnfoldResponse::operator=(const RooUnfoldResponse& rhs)
{
  // RooUnfoldResponseT<class Hist, class Hist2D> assignment operator
  if (this == &rhs) return *this;
  Reset();
  SetNameTitle(rhs.GetName(), rhs.GetTitle());
  return Setup(rhs);
}



RooUnfoldResponse&
RooUnfoldResponse::Setup(const RooUnfoldResponse& rhs)
{
  // Copy data from another RooUnfoldResponseT<class Hist, class Hist2D>
  _overflow= rhs._overflow;
  return Setup(clone(rhs.Hmeasured()), clone(rhs.Htruth()), clone(rhs.Hresponse()));
}


RooUnfoldResponse&
RooUnfoldResponse::Reset()
{
  // Resets object to initial state.
  if(_mes) delete  _mes;
  if(_tru) delete  _tru;
  if(_res) delete  _res;
  if(_fak) delete _fak;
  setup();
  return *this;
}


RooUnfoldResponse::RooUnfoldResponse(Int_t nb, Double_t xlo, Double_t xhi,
                                                     const char* name, const char* title)
  : RooUnfoldResponseT(name, title)
{
  // RooUnfoldResponse<class TH1, class TH2> constructor - simple 1D case with same binning, measured vs truth
  Setup(nb, xlo, xhi);
}


RooUnfoldResponse::RooUnfoldResponse(Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi,
                                                     const char* name, const char* title)
  : RooUnfoldResponseT(name, title)
{
  // RooUnfoldResponse<class TH1, class TH2> constructor - simple 1D case
  Setup(nm, mlo, mhi, nt, tlo, thi);
}


RooUnfoldResponse::RooUnfoldResponse(const TH1* measured, const TH1* truth, const TH2* response,
                                                     const char* name, const char* title, bool overflow)
  : RooUnfoldResponseT(name, title)
{
  // RooUnfoldResponse<class TH1, class TH2> constructor - create from already-filled histograms
  // "response" gives the response matrix, measured X truth.
  // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  // in "truth" for unmeasured events (inefficiency).
  // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  // to indicate, respectively, no fakes and/or no inefficiency.
  this->_overflow = overflow;
  Setup(measured, truth, response);
}


RooUnfoldResponse::RooUnfoldResponse(const TH1* measured, const TH1* truth,
                                                     const char* name, const char* title, bool overflow)
  : RooUnfoldResponseT(name, title)
{
  // RooUnfoldResponse<class TH1, class TH2> constructor - measured and truth only used for shape
  this->_overflow = overflow;
  Setup(measured, truth);
}


void RooUnfoldResponse::Add(const RooUnfoldResponse& rhs)
{
  // Add another RooUnfoldResponse<class TH1, class TH2>, accumulating contents
  if (_res == 0) {
    Setup(rhs);
    return;
  }
  assert (_mdim==rhs._mdim);
  assert (_tdim==rhs._tdim);
  assert (_mes != 0 && rhs._mes != 0);
  assert (_fak != 0 && rhs._fak != 0);
  assert (_tru != 0 && rhs._tru != 0);
  assert (_res != 0 && rhs._res != 0);
  if (Cached()) ClearCache();
  _mes->Add(rhs._mes);
  _fak->Add(rhs._fak);
  _tru->Add(rhs._tru);
  _res->Add(rhs._res);
}


 Long64_t
RooUnfoldResponse::Merge (TCollection* others)
{
  // Add all RooUnfoldResponse<class TH1, class TH2> objects in the collection to this one.
  // This allows merging with hadd and TFileMerger.
  for (TIter it= others; TObject* o= it();) {
    if (RooUnfoldResponse* other= dynamic_cast<RooUnfoldResponse*>(o))
      Add (*other);
  }
  return Long64_t(::entries(_res));
}

RooUnfoldResponse&
RooUnfoldResponse::Setup(Int_t nb, Double_t xlo, Double_t xhi)
{
  // constructor -  simple 1D case with same binning, measured vs truth
  return Setup(nb, xlo, xhi, nb, xlo, xhi);
}




 RooUnfoldResponse&
RooUnfoldResponse::Setup(Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi)
{
  // set up simple 1D case
  Reset();
  _mdim= _tdim= 1;
  _nm= nm;
  _nt= nt;
  _mes= createHist<TH1>("measured", "Measured",   Variable<TH1>(nm, mlo, mhi,"xm"));
  _fak= createHist<TH1>("fakes",    "Fakes",      Variable<TH1>(nm, mlo, mhi,"xm"));
  _tru= createHist<TH1>("truth",    "Truth",      Variable<TH1>(nt, tlo, thi,"xt"));
  _res= createHist<TH2>("response", "Response", Variable<TH2>(nm, mlo, mhi, "xm"), Variable<TH2>(nt, tlo, thi, "xt"));
  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup(const TH1* measured, const TH1* truth)
{
  // set up - measured and truth only used for shape
  Reset();
  _mes= createHist<TH1>(measured->GetName(),measured->GetTitle(), vars(measured));
  _fak= createHist<TH1>("fakes","Fakes",vars(measured));
  _tru= createHist<TH1>("truth",truth->GetTitle(), vars(truth));
  _mdim= dim(_mes);
  _tdim= dim(_tru);
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  SetNameTitleDefault();
  _nm= nBins(_mes);
  _nt= nBins(_tru);
  _res=createHist<TH2>(GetName(), GetTitle(), Variable<TH2>(_nm, 0.0, _nm, "xm"), Variable<TH2>(_nt, 0.0, _nt, "xt"));
  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup(const TH1* measured, const TH1* truth, const TH2* response)
{
  // Set up from already-filled histograms.
  // "response" gives the response matrix, measured X truth.
  // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  // in "truth" for unmeasured events (inefficiency).
  // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  // to indicate, respectively, no fakes and/or no inefficiency.
  Reset();
  _res= clone(response);
  if (measured) {
    _mes= clone(measured);
    _mdim= dim(_mes);
  } else {
    _mes= createHist<TH1>("measured", "Measured", Variable<TH1>(nBins(response,RooUnfolding::X), 0.0, 1.0, "xm"));
    _mdim= 1;
  }
  _fak= createHist<TH1>("fakes","Fakes",vars(measured));
  if (truth) {
    _tru= clone(truth);
    _tdim= dim(_tru);
  } else {
    _tru= createHist<TH1>("truth",    "Truth",    Variable<TH1>(nBins(response,RooUnfolding::Y), 0.0, 1.0, "xt"));
    _tdim= 1;
  }
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  _nm= nBins(_mes);
  _nt= nBins(_tru);
  if (_nm != nBins(_res,RooUnfolding::X) || _nt != nBins(_res,RooUnfolding::Y)) {
    cerr << "Warning: RooUnfoldResponse<class TH1, class TH2> measured X truth is " << _nm << " X " << _nt
         << ", but matrix is " << nBins(_res,RooUnfolding::X)<< " X " << nBins(_res,RooUnfolding::Y) << endl;
  }

  if (!measured || ::entries(_mes) == 0.0) {
    // Similar to _res->ProjectionX() but without stupid reset of existing histograms
    // Always include under/overflows in sum of truth.
    projectX(_res,_mes,_overflow);
  } else {
    // Fill fakes from the difference of _mes - _res->ProjectionX()
    // Always include under/overflows in sum of truth.
    subtractProjectX(_res,_mes,_fak,_overflow);
  }

  if (!truth || ::entries(_tru) == 0.0) {
    // similar to _res->ProjectionY() but without stupid reset of existing histograms
    // Always include under/overflows in sum of measurements.
    projectY(_res,_tru,_overflow);
  }

  SetNameTitleDefault();
  return *this;
}


 Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t xt, Double_t w)
{
  // Fill 1D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==1 && _tdim==1);
  if (Cached()) ClearCache();
  fill(_mes,xr,w);
  fill(_tru,xt,w);
  return fill(_res,xr,xt,w);
}

 Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t yr, Double_t xt, Double_t yt, Double_t w)
{
  // Fill 2D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==2 && _tdim==2);
  if (Cached()) ClearCache();
  fill((TH2*)_mes,xr, yr, w);
  fill((TH2*)_tru,xt, yt, w);
  return fill(_res,binCenter(_res,findBin (_mes, xr, yr)+1,RooUnfolding::X),binCenter(_res,findBin (_tru, xt, yt)+1,RooUnfolding::Y), w);
}

 Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t yr, Double_t zr, Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  // Fill 3D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==3 && _tdim==3);
  if (Cached()) ClearCache();
  fill(_mes,xr, yr, zr, w);
  fill(_tru,xt, yt, zt, w);
  return fill(_res,binCenter(_res,findBin (_mes, xr, yr, zr)+1,RooUnfolding::X),binCenter(_res,findBin (_tru, xt, yt, zt)+1,RooUnfolding::Y), w);  
}


 Int_t
RooUnfoldResponse::FindBin(const TH1* h, Double_t x, Double_t y, Double_t z)
{
  // Get vector index (0..nx*ny*nz-1) for bin containing (x,y,z) coordinates
  Int_t nx=   nBins(h,RooUnfolding::X);
  Int_t ny=   nBins(h,RooUnfolding::Y);
  Int_t nz=   nBins(h,RooUnfolding::Z);
  Int_t binx= findBin(h,x,RooUnfolding::X) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny*nz;
  Int_t biny= findBin(h,y,RooUnfolding::Y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny*nz;
  Int_t binz= findBin(h,z,RooUnfolding::Z) - 1;
  if (binz <  0)  return -1;
  if (binz >= nz) return nx*ny*nz;
  return binx + nx*(biny + ny*binz);
}

Int_t
RooUnfoldResponse::Miss1D (Double_t xt, Double_t w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 1D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_tdim==1);
  if (Cached()) ClearCache();
  return fill(_tru, xt, w);
}

 Int_t
RooUnfoldResponse::Miss2D (Double_t xt, Double_t yt, Double_t w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 2D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_tdim==2);
  if (Cached()) ClearCache();
  return fill(_tru, xt, yt, w);
}

 Int_t
RooUnfoldResponse::Miss (Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 3D Response Matrix
  assert (_tru != 0);
  assert (_tdim==3);
  if (Cached()) ClearCache();
  return fill(_tru, xt, yt, zt, w);
}

 Int_t
RooUnfoldResponse::Fake1D (Double_t xr, Double_t w)
{
  // Fill fake event (reconstructed event with no truth) into 1D Response Matrix (with weight)
  assert (_fak != 0 && _mes != 0);
  assert (_mdim==1);
  if (Cached()) ClearCache();
  fill(_mes,xr, w);
  return fill(_fak, xr, w);
}

 Int_t
RooUnfoldResponse::Fake2D (Double_t xr, Double_t yr, Double_t w)
{
  // Fill fake event (reconstructed event with no truth) into 2D Response Matrix (with weight)
  assert (_mes != 0);
  assert (_mdim==2);
  if (Cached()) ClearCache();
  fill(_fak, xr, yr, w);
  return fill(_mes, xr, yr, w);
}

 Int_t
RooUnfoldResponse::Fake (Double_t xr, Double_t yr, Double_t zr, Double_t w)
{
  // Fill fake event (reconstructed event with no truth) into 3D Response Matrix
  assert (_mes != 0);
  assert (_mdim==3);
  if (Cached()) ClearCache();
  fill(_mes, xr, yr, zr, w);
  return fill(_fak, xr, yr, zr, w);
}


Long64_t RooUnfoldResponse::FakeEntries() const
{
  // Return number of fake entries
  return _fak ? entries(_fak) : 0.0;
}

ClassImp (RooUnfoldResponse);


#ifndef NOROOFIT

#include "RooRealVar.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooAddition.h"
#include "RooConstVar.h"
#include "RooHistFunc.h"



namespace {
  template<class T> void makeList(const RooAbsCollection* c,std::vector<T*>& vars,TPRegexp* re = 0){
    RooFIter itr(c->fwdIterator());
    RooRealVar* obj = NULL;
    while((obj = (RooRealVar*)itr.next())){
      if(!obj->InheritsFrom(RooRealVar::Class())) continue;
      if(!re || re->Match(obj->GetName()))
        vars.push_back(obj);
    }
  }
  
  template<class T> std::vector<T*> makeList(const RooAbsCollection* c){
    std::vector<RooAbsArg*> vars;
    makeList<T>(c,vars);
    return vars;
  }
}

template <> RooUnfolding::RooFitHist*
RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::HresponseNoOverflow() const
{
  return const_cast<RooUnfolding::RooFitHist*>(Hresponse());
}



template class RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;

RooFitUnfoldResponse::RooFitUnfoldResponse(const char* name, const char* title, RooUnfolding::RooFitHist* response, RooUnfolding::RooFitHist* truth, RooUnfolding::RooFitHist* reco, bool density) : 
  RooUnfoldResponseT(name,title,response,truth,reco,false,density)
{
  // inherited constructor
}


RooFitUnfoldResponse::RooFitUnfoldResponse(const char* name, const char* title, RooAbsReal* response, RooAbsReal* truth, RooAbsReal* reco, RooAbsReal* fakes, RooRealVar* obs_truth, RooRealVar* obs_reco, bool density){
  if(!truth->dependsOn(*obs_truth)) throw std::runtime_error("truth histogram does not depend on truth observable!");
  if(!reco->dependsOn(*obs_reco)) throw std::runtime_error("reco histogram does not depend on reco observable!");
  if(!response->dependsOn(*obs_truth)) throw std::runtime_error("response histogram does not depend on truth observable!");
  if(!response->dependsOn(*obs_reco)) throw std::runtime_error("response histogram does not depend on reco observable!");
  
  this->_mdim = 1;
  this->_tdim = 1;
  this->_nm = obs_reco->getBins();
  this->_nt = obs_truth->getBins();

  TPRegexp gamma("gamma_stat_.*");

  std::vector<RooRealVar*> allvars;
  RooArgSet* c_reco_vars = reco->getVariables();
  makeList<RooRealVar>(c_reco_vars,allvars,&gamma);
  delete c_reco_vars;
  if(fakes){
    RooArgSet* c_fake_vars = fakes->getVariables();
    makeList<RooRealVar>(c_fake_vars,allvars,&gamma);
    delete c_fake_vars;
  }
  RooArgSet* c_truth_vars = truth->getVariables();
  makeList<RooRealVar>(c_truth_vars,allvars,&gamma);
  delete c_truth_vars;

  this->_mes = new RooFitHist(reco,obs_reco,allvars);
  this->_fak = fakes ? new RooFitHist(fakes,obs_reco,allvars) : 0;
  this->_tru = new RooFitHist(truth,obs_truth,allvars);
  this->_res = new RooFitHist(response,obs_truth,obs_reco,allvars);
  this->_overflow = 0;
  this->_density = density;
}

RooFitUnfoldResponse::RooFitUnfoldResponse(const char* name, const char* title, RooAbsReal* response, RooAbsReal* truth, RooAbsReal* reco, RooAbsReal* fakes, const RooAbsCollection* observables, bool density)
   : RooUnfoldResponseT(name,title) {
  // standard constructor

  RooArgSet obs_in;
  obs_in.add(*observables);
  RooArgSet* obs = response->getObservables(&obs_in);
  RooArgSet* obsset_reco = reco->getObservables(&obs_in);
  RooArgSet* obsset_truth = truth->getObservables(&obs_in);    
 
  if(obs->getSize() != 2) throw std::runtime_error(TString::Format("unsupported dimensionality for response: %d",obs->getSize()).Data());
  if(obsset_reco->getSize() != 1) throw std::runtime_error(TString::Format("unsupported dimensionality for reco: %d",obsset_reco->getSize()).Data());
  if(obsset_truth->getSize() != 1) throw std::runtime_error(TString::Format("unsupported dimensionality for truth: %d",obsset_truth->getSize()).Data());

  RooRealVar* obs_reco  = (RooRealVar*)obsset_reco->first();
  RooRealVar* obs_truth = (RooRealVar*)obsset_truth->first();
  
  this->_mdim = 1;
  this->_tdim = 1;
  this->_nm = obs_reco->getBins();
  this->_nt = obs_truth->getBins();
  this->_mes = new RooFitHist(reco,obs_reco);
  this->_fak = fakes ? new RooFitHist(fakes,obs_reco) : 0;
  this->_tru = new RooFitHist(truth,obs_truth);
  this->_res = new RooFitHist(response,::makeList<RooAbsArg>(obs));
  this->_overflow = 0;
  this->_density = density;
}


RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHistSum(RooAbsReal* a, RooAbsReal* b, double ca, double cb){
  TString name(TString::Format("%s_minus_%s",a->GetName(),b->GetName()));
  RooArgList funcs;
  funcs.add(*a);
  funcs.add(*b);
  RooConstVar* cA = new RooConstVar("cA","cA",ca);
  RooConstVar* cB = new RooConstVar("cB","cB",cb);  
  RooArgList coefs;
  coefs.add(*cA);
  coefs.add(*cB);
  RooAddition* add = new RooAddition(name.Data(),name.Data(),funcs,coefs);
  return this->makeHist(add);
}

RooHistFunc* RooFitUnfoldResponse::makeHistFunc(RooDataHist* dhist){
  if(!dhist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    if(dhist->get()->find(*this->_mes->obs(i))) v.push_back(this->_mes->obs(i));
  }
  for(size_t i=0; i<this->_tru->dim(); ++i){
    if(dhist->get()->find(*this->_tru->obs(i))) v.push_back(this->_tru->obs(i));
  }
  if(v.size() == 0){
    throw std::runtime_error(TString::Format("unable to construct histogram from RooDataHist '%s', does not seem to depend on any known observable!",dhist->GetName()).Data());
  }
  return RooUnfolding::makeHistFunc(dhist,v);
}

RooHistFunc* RooFitUnfoldResponse::makeHistFuncMeasured(const TH1* hist){
  if(!hist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    v.push_back(this->_mes->obs(i));
  }
  RooDataHist* dhist = RooUnfolding::convertTH1(hist,v,this->_overflow,this->_density);
  return RooUnfolding::makeHistFunc(dhist,v);
}
RooHistFunc* RooFitUnfoldResponse::makeHistFuncTruth(const TH1* hist){
  if(!hist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_tru->dim(); ++i){
    v.push_back(this->_tru->obs(i));
  }
  RooDataHist* dhist = RooUnfolding::convertTH1(hist,v,this->_overflow,this->_density);
  return RooUnfolding::makeHistFunc(dhist,v);
}

RooAbsPdf* RooFitUnfoldResponse::makeHistPdf(RooDataHist* dhist){
  if(!dhist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    if(dhist->get()->find(*this->_mes->obs(i))) v.push_back(this->_mes->obs(i));
  }
  for(size_t i=0; i<this->_tru->dim(); ++i){
    if(dhist->get()->find(*this->_tru->obs(i))) v.push_back(this->_tru->obs(i));
  }
  if(v.size() == 0){
    throw std::runtime_error(TString::Format("unable to construct histogram from RooDataHist '%s', does not seem to depend on any known observable!",dhist->GetName()).Data());
  }
  return RooUnfolding::makeHistPdf(dhist,v);
}

RooAbsPdf* RooFitUnfoldResponse::makeHistPdfMeasured(const TH1* hist){
  if(!hist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    v.push_back(this->_mes->obs(i));
  }
  RooDataHist* dhist = RooUnfolding::convertTH1(hist,v,this->_overflow,this->_density);
  return RooUnfolding::makeHistPdf(dhist,v);
}
RooAbsPdf* RooFitUnfoldResponse::makeHistPdfTruth(const TH1* hist){
  if(!hist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_tru->dim(); ++i){
    v.push_back(this->_tru->obs(i));
  }
  RooDataHist* dhist = RooUnfolding::convertTH1(hist,v,this->_overflow,this->_density);
  return RooUnfolding::makeHistPdf(dhist,v);
}

RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHistMeasured(const TH1* hist){
  return this->makeHist(this->makeHistFuncMeasured(hist));
}
RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHistTruth(const TH1* hist){
  return this->makeHist(this->makeHistFuncTruth(hist));  
}


RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHist(RooAbsReal* object){
  if(!object) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    if(object->dependsOn(*this->_mes->obs(i))) v.push_back(this->_mes->obs(i));
  }
  for(size_t i=0; i<this->_tru->dim(); ++i){
    if(object->dependsOn(*this->_tru->obs(i))) v.push_back(this->_tru->obs(i));
  }
  if(v.size() == 0){
    throw std::runtime_error(TString::Format("unable to construct histogram from RooAbsReal '%s', does not seem to depend on any known observable!",object->GetName()).Data());
  }
  return new RooUnfolding::RooFitHist(object,v);
}

RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHist(RooDataHist* object){
  if(!object) return NULL;
  std::vector<RooAbsArg*> v;
  TPRegexp obs("obs_.*");
  makeList<RooAbsArg>(object->get(),v,&obs);
  return new RooUnfolding::RooFitHist(object,v);
}

ClassImp (RooFitUnfoldResponse);

#endif
