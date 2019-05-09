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
#include "RooUnfoldFitHelpers.h"

#include <iostream>
#include <assert.h>
#include <cmath>

#include "TClass.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TCollection.h"

#include "TH1.h"
#include "TH2.h"
#include "RooAbsReal.h"


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


template <class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT (const RooUnfoldResponseT<Hist,Hist2D>& rhs)
  : TNamed (rhs.GetName(), rhs.GetTitle())
{
  // RooUnfoldResponseT<class Hist, class Hist2D> copy constructor
  Init();
  Setup (rhs);
}

template <class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT (Int_t nb, Double_t xlo, Double_t xhi,
                                                     const char* name, const char* title)
  : TNamed (name, title)
{
  // RooUnfoldResponseT<class Hist, class Hist2D> constructor - simple 1D case with same binning, measured vs truth
  Init();
  Setup (nb, xlo, xhi);
}

template <class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi,
                                                     const char* name, const char* title)
  : TNamed (name, title)
{
  // RooUnfoldResponseT<class Hist, class Hist2D> constructor - simple 1D case
  Init();
  Setup (nm, mlo, mhi, nt, tlo, thi);
}

template <class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT (const Hist* measured, const Hist* truth, const Hist2D* response,
                                                     const char* name, const char* title, bool overflow)
  : TNamed (name, title)
{
  // RooUnfoldResponseT<class Hist, class Hist2D> constructor - create from already-filled histograms
  // "response" gives the response matrix, measured X truth.
  // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  // in "truth" for unmeasured events (inefficiency).
  // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  // to indicate, respectively, no fakes and/or no inefficiency.
  Init();
  this->_overflow = overflow;
  Setup (measured, truth, response);
}

template <class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT (const Hist* measured, const Hist* truth,
                                                     const char* name, const char* title, bool overflow)
  : TNamed (name, title)
{
  // RooUnfoldResponseT<class Hist, class Hist2D> constructor - measured and truth only used for shape
  Init();
  this->_overflow = overflow;
  Setup (measured, truth);
}

template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::operator= (const RooUnfoldResponseT<Hist,Hist2D>& rhs)
{
  // RooUnfoldResponseT<class Hist, class Hist2D> assignment operator
  if (this == &rhs) return *this;
  Reset();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  return Setup (rhs);
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::Add (const RooUnfoldResponseT<Hist,Hist2D>& rhs){
  throw std::runtime_error(TString::Format("adding not supported for %s",this->ClassName()).Data());
}

template <> void
RooUnfoldResponseT<TH1,TH2>::Add (const RooUnfoldResponseT<TH1,TH2>& rhs)
{
  // Add another RooUnfoldResponseT<class Hist, class Hist2D>, accumulating contents
  if (_res == 0) {
    Setup (rhs);
    return;
  }
  assert (_mdim==rhs._mdim);
  assert (_tdim==rhs._tdim);
  assert (_mes != 0 && rhs._mes != 0);
  assert (_fak != 0 && rhs._fak != 0);
  assert (_tru != 0 && rhs._tru != 0);
  assert (_res != 0 && rhs._res != 0);
  if (_cached) ClearCache();
  _mes->Add(rhs._mes);
  _fak->Add(rhs._fak);
  _tru->Add(rhs._tru);
  _res->Add(rhs._res);
}


template <class Hist, class Hist2D> Long64_t
RooUnfoldResponseT<Hist,Hist2D>::Merge (TCollection* others)
{
  // Add all RooUnfoldResponseT<class Hist, class Hist2D> objects in the collection to this one.
  // This allows merging with hadd and TFileMerger.
  for (TIter it= others; TObject* o= it();) {
    if (RooUnfoldResponseT<Hist,Hist2D>* other= dynamic_cast<RooUnfoldResponseT<Hist,Hist2D>*>(o))
      Add (*other);
  }
  return Long64_t(entries(_res));
}


template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Reset()
{
  // Resets object to initial state.
  ClearCache();
  if(_mes) maybeDelete( _mes);
  if(_fak) delete _fak;
  if(_tru) maybeDelete( _tru);
  if(_res) maybeDelete( _res);
  return Setup();
}

template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Init()
{
  _overflow= 0;
  return Setup();
}

template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Setup()
{
  _tru= _mes= _fak= 0;
  _res= 0;
  _vMes= _eMes= _vFak= _vTru= _eTru= 0;
  _mRes= _eRes= 0;
  _nm= _nt= _mdim= _tdim= 0;
  SetNameTitleDefault ("response", "Response");
  _cached= false;
  return *this;
}

template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Setup (const RooUnfoldResponseT<Hist,Hist2D>& rhs)
{
  // Copy data from another RooUnfoldResponseT<class Hist, class Hist2D>
  _overflow= rhs._overflow;
  return Setup (rhs.Hmeasured(), rhs.Htruth(), rhs.Hresponse());
}

template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Setup (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi)
{
  // set up simple 1D case
  Reset();
  _mdim= _tdim= 1;
  _nm= nm;
  _nt= nt;
  _mes= createHist<Hist>("measured", "Measured",   Variable<Hist>(nm, mlo, mhi,"xm"));
  _fak= createHist<Hist>("fakes",    "Fakes",      Variable<Hist>(nm, mlo, mhi,"xm"));
  _tru= createHist<Hist>("truth",    "Truth",      Variable<Hist>(nt, tlo, thi,"xt"));
  _res= createHist<Hist2D>("response", "Response", Variable<Hist2D>(nm, mlo, mhi, "xm"), Variable<Hist2D>(nt, tlo, thi, "xt"));
  return *this;
}

template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Setup (const Hist* measured, const Hist* truth)
{
  // set up - measured and truth only used for shape
  Reset();
  _mes= createHist<Hist>(measured->GetName(),measured->GetTitle(), vars(measured));
  _fak= createHist<Hist>("fakes","Fakes",vars(measured));
  _tru= createHist<Hist>("truth",truth->GetTitle(), vars(truth));
  _mdim= dim(_mes);
  _tdim= dim(_tru);
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  SetNameTitleDefault();
  _nm= nBins(_mes);
  _nt= nBins(_tru);
  _res=createHist<Hist2D>(GetName(), GetTitle(), Variable<Hist2D>(_nm, 0.0, _nm, "xm"), Variable<Hist2D>(_nt, 0.0, _nt, "xt"));
  return *this;
}

template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Setup (const Hist* measured, const Hist* truth, const Hist2D* response)
{
  // Set up from already-filled histograms.
  // "response" gives the response matrix, measured X truth.
  // "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  // but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  // in "truth" for unmeasured events (inefficiency).
  // "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  // to indicate, respectively, no fakes and/or no inefficiency.
  Reset();
  _res= maybeCopy(response);
  if (measured) {
    _mes= maybeCopy(measured);
    _mdim= dim(_mes);
  } else {
    _mes= createHist<Hist>("measured", "Measured", Variable<Hist>(nBins(response,RooUnfolding::X), 0.0, 1.0, "xm"));
    _mdim= 1;
  }
  _fak= createHist<Hist>("fakes","Fakes",vars(measured));
  if (truth) {
    _tru= maybeCopy(truth);
    _tdim= dim(_tru);
  } else {
    _tru= createHist<Hist>("truth",    "Truth",    Variable<Hist>(nBins(response,RooUnfolding::Y), 0.0, 1.0, "xt"));
    _tdim= 1;
  }
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  _nm= nBins(_mes);
  _nt= nBins(_tru);
  if (_nm != nBins(_res,RooUnfolding::X) || _nt != nBins(_res,RooUnfolding::Y)) {
    cerr << "Warning: RooUnfoldResponseT<class Hist, class Hist2D> measured X truth is " << _nm << " X " << _nt
         << ", but matrix is " << nBins(_res,RooUnfolding::X)<< " X " << nBins(_res,RooUnfolding::Y) << endl;
  }

  if (!measured || entries(_mes) == 0.0) {
    // Similar to _res->ProjectionX() but without stupid reset of existing histograms
    // Always include under/overflows in sum of truth.
    projectX(_res,_mes,_overflow);
  } else {
    // Fill fakes from the difference of _mes - _res->ProjectionX()
    // Always include under/overflows in sum of truth.
    subtractProjectX(_res,_mes,_fak,_overflow);
  }

  if (!truth || entries(_tru) == 0.0) {
    // similar to _res->ProjectionY() but without stupid reset of existing histograms
    // Always include under/overflows in sum of measurements.
    projectY(_res,_tru,_overflow);
  }

  SetNameTitleDefault();
  return *this;
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::ClearCache()
{
  if(_vMes) delete _vMes; _vMes= 0;
  if(_eMes) delete _eMes; _eMes= 0;
  if(_vFak) delete _vFak; _vFak= 0;
  if(_vTru) delete _vTru; _vTru= 0;
  if(_eTru) delete _eTru; _eTru= 0;
  if(_mRes) delete _mRes; _mRes= 0;
  if(_eRes) delete _eRes; _eRes= 0;
  _cached= false;
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fill (Double_t xr, Double_t xt, Double_t w)
{
  // Fill 1D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==1 && _tdim==1);
  if (_cached) ClearCache();
  fill(_mes,xr,w);
  fill(_tru,xt,w);
  return fill(_res,xr,xt,w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fill (Double_t xr, Double_t yr, Double_t xt, Double_t yt, Double_t w)
{
  // Fill 2D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==2 && _tdim==2);
  if (_cached) ClearCache();
  fill((Hist2D*)_mes,xr, yr, w);
  fill((Hist2D*)_tru,xt, yt, w);
  return fill(_res,binCenter(_res,findBin (_mes, xr, yr)+1,RooUnfolding::X),binCenter(_res,findBin (_tru, xt, yt)+1,RooUnfolding::Y), w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fill (Double_t xr, Double_t yr, Double_t zr, Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  // Fill 3D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==3 && _tdim==3);
  if (_cached) ClearCache();
  fill(_mes,xr, yr, zr, w);
  fill(_tru,xt, yt, zt, w);
  return fill(_res,binCenter(_res,findBin (_mes, xr, yr, zr)+1,RooUnfolding::X),binCenter(_res,findBin (_tru, xt, yt, zt)+1,RooUnfolding::Y), w);  
}


template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::FindBin(const Hist* h, Double_t x, Double_t y, Double_t z)
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

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Miss1D (Double_t xt, Double_t w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 1D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_tdim==1);
  if (_cached) ClearCache();
  return fill(_tru, xt, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Miss2D (Double_t xt, Double_t yt, Double_t w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 2D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_tdim==2);
  if (_cached) ClearCache();
  return fill(_tru, xt, yt, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Miss (Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 3D Response Matrix
  assert (_tru != 0);
  assert (_tdim==3);
  if (_cached) ClearCache();
  return fill(_tru, xt, yt, zt, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fake1D (Double_t xr, Double_t w)
{
  // Fill fake event (reconstructed event with no truth) into 1D Response Matrix (with weight)
  assert (_fak != 0 && _mes != 0);
  assert (_mdim==1);
  if (_cached) ClearCache();
  fill(_mes,xr, w);
  return fill(_fak, xr, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fake2D (Double_t xr, Double_t yr, Double_t w)
{
  // Fill fake event (reconstructed event with no truth) into 2D Response Matrix (with weight)
  assert (_mes != 0);
  assert (_mdim==2);
  if (_cached) ClearCache();
  fill(_fak, xr, yr, w);
  return fill(_mes, xr, yr, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fake (Double_t xr, Double_t yr, Double_t zr, Double_t w)
{
  // Fill fake event (reconstructed event with no truth) into 3D Response Matrix
  assert (_mes != 0);
  assert (_mdim==3);
  if (_cached) ClearCache();
  fill(_mes, xr, yr, zr, w);
  return fill(_fak, xr, yr, zr, w);
}


template <class Hist, class Hist2D> Hist2D*
RooUnfoldResponseT<Hist,Hist2D>::HresponseNoOverflow() const
{
  const Hist2D* res = Hresponse();
  TVectorD vals(h2v<Hist>(res,_overflow));
  TVectorD errs(h2ve<Hist>(res,_overflow));  
  return createHist<Hist2D>(vals,errs,res->GetName(),res->GetTitle(),vars(res),_overflow);
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::Print (Option_t* /* option */) const
{
  printMatrix (Mresponse(), Form("%s response matrix",GetTitle()));
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
    resultvect= h2v (truth, _overflow);
  } else {
    resultvect= Vtruth();
  }

  resultvect *= Mresponse();   // v= A*v

  // Turn results vector into properly binned histogram
  const Hist* t = Hmeasured();
  Hist* result= createHist<Hist>(resultvect, t->GetName(),name, vars(t), _overflow);
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
  TMatrixD values = h2m(hres);
  TMatrixD errors = h2me(hres);
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

  Hist2D* smeared = createHist<Hist2D>(values,errors,hres->GetName(),hres->GetTitle(),var(hres,X),var(hres,Y));
  TString name= GetName();
  name += "_toy";
  RooUnfoldResponseT<Hist,Hist2D>* res= new RooUnfoldResponseT<Hist,Hist2D> (this->_mes,this->_tru,smeared,name.Data(),this->GetTitle());
  delete smeared;
  
  return res;
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::SetNameTitleDefault (const char* defname, const char* deftitle)
{
  // Set object name and title
  const char* s= GetName();
  if (s[0] == '\0') {
    if (_res) s= _res->GetName();
    if (s[0] == '\0') {
      if (defname) SetName (defname);
      else if (_mes && _tru) {
        TString n= _mes->GetName();
        if (n.Length()) n.Append ("_");
        n.Append (_tru->GetName());
        if (!n.Length()) n= "response";
        SetName (n);
      }
    } else
      SetName (s);
  }
  s= GetTitle();
  if (s[0] == '\0') {
    if (_res) s= _res->GetTitle();
    if (s[0] == '\0') {
      if (deftitle) SetTitle (deftitle);
      else if (_mes && _tru) {
        TString n= _tru->GetTitle();
        if (n.Length()) n.Append (" #rightarrow ");
        n.Append (_mes->GetTitle());
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
  Init();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const char*    name, const char*    title)
  : TNamed(name,title)
{
  // RooUnfoldResponseT<Hist,Hist2D> default named constructor. Use Setup() to set values.
  Init();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const TString& name, const TString& title)
  : TNamed(name,title)
{
  // RooUnfoldResponseT<Hist,Hist2D> default named constructor. Use Setup() to set values.
  Init();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::~RooUnfoldResponseT()
{
  // RooUnfoldResponseT<Hist,Hist2D> destructor
  Reset();
}


template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>& RooUnfoldResponseT<Hist,Hist2D>::Setup (Int_t nb, Double_t xlo, Double_t xhi)
{
  // constructor -  simple 1D case with same binning, measured vs truth
  return Setup (nb, xlo, xhi, nb, xlo, xhi);
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
  if (!_vMes) _cached= (_vMes= new TVectorD(h2v  (_mes, _overflow)));
  return *_vMes;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vfakes() const
{
  // Fakes distribution as a TVectorD
  if (!_vFak) _cached= (_vFak= new TVectorD(h2v  (_fak, _overflow)));
  return *_vFak;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Emeasured() const
{
  // Measured distribution errors as a TVectorD
  if (!_eMes) _cached= (_eMes= new TVectorD(h2ve (_mes, _overflow)));
  return *_eMes;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vtruth() const
{
  // Truth distribution as a TVectorD
  if (!_vTru) _cached= (_vTru= new TVectorD(h2v  (_tru, _overflow))); 
  return *_vTru;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Etruth() const
{
  // Truth distribution errors as a TVectorD
  if (!_eTru) _cached= (_eTru= new TVectorD(h2ve (_tru, _overflow))); 
  return *_eTru;
}

template<class Hist, class Hist2D>
const TMatrixD& RooUnfoldResponseT<Hist,Hist2D>::Mresponse() const
{
  // Response matrix as a TMatrixD: (row,column)=(measured,truth)
  if (!_mRes) _cached= (_mRes= new TMatrixD(h2mNorm  (_res, _tru, _overflow))); 
  return *_mRes;
}

template<class Hist, class Hist2D>
const TMatrixD& RooUnfoldResponseT<Hist,Hist2D>::Eresponse() const
{
  // Response matrix errors as a TMatrixD: (row,column)=(measured,truth)
  if (!_eRes) _cached= (_eRes= new TMatrixD(h2meNorm (_res, _tru, _overflow))); 
  return *_eRes;
}


template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::operator() (Int_t r, Int_t t) const
{
  // Response matrix element (measured,truth)
  return Mresponse()(r,t);
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::Miss (Double_t xt)
{
  // Fill missed event into 1D Response Matrix
  return Miss1D(xt);
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::Miss (Double_t xt, Double_t w)
{
  // Fill missed event into 1D (with weight) or 2D Response Matrix
  return _tdim==2 ? Miss2D(xt,w) : Miss1D(xt,w);
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::Miss (Double_t xt, Double_t yt, Double_t w)
{
  // Fill missed event into 2D (with weight) or 3D Response Matrix
  return _tdim==3 ? Miss(xt,yt,w,1.0) : Miss2D(xt,yt,w);
}


template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::Fake (Double_t xr)
{
  // Fill fake event into 1D Response Matrix
  return Fake1D(xr);
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::Fake (Double_t xr, Double_t w)
{
  // Fill fake event into 1D (with weight) or 2D Response Matrix
  return _mdim==2 ? Fake2D(xr,w) : Fake1D(xr,w);
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::Fake (Double_t xr, Double_t yr, Double_t w)
{
  // Fill fake event into 2D (with weight) or 3D Response Matrix
  return _mdim==3 ? Fake(xr,yr,w,1.0) : Fake2D(xr,yr,w);
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
Double_t RooUnfoldResponseT<Hist,Hist2D>::FakeEntries() const
{
  // Return number of fake entries
  return _fak ? entries(_fak) : 0.0;
}

template class RooUnfoldResponseT<TH1,TH2>;
ClassImp (RooUnfoldResponse);

template class RooUnfoldResponseT<RooAbsReal,RooAbsReal>;
ClassImp (RooFitUnfoldResponse);
