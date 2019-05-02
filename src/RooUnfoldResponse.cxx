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

#include <iostream>
#include <assert.h>
#include <cmath>

#include "TClass.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TCollection.h"

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
    if      (_ndim==1) bin= RooUnfoldResponseT<Hist,Hist2D>::FindBin (mes, x[0]);
    else if (_ndim==2) bin= RooUnfoldResponseT<Hist,Hist2D>::FindBin (mes, x[0], x[1]);
    else               bin= RooUnfoldResponseT<Hist,Hist2D>::FindBin (mes, x[0], x[1], x[2]);
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
      Int_t j= RooUnfoldResponseT<Hist,Hist2D>::GetBin(tru, i);
      Int_t jx, jy, jz;
      if (_ndim>=2) tru->GetBinXYZ (j, jx, jy, jz);
      Double_t fv;
      if (_eps<=0.0) {
        if (_ndim>=2)
          fv= _func->Eval (tru->GetXaxis()->GetBinCenter(jx),
                           tru->GetYaxis()->GetBinCenter(jy),
                           tru->GetZaxis()->GetBinCenter(jz));
        else
          fv= _func->Eval (tru->GetBinCenter(j));
      } else {
        if        (_ndim==1) {
          Double_t tw= tru->GetBinWidth(j), tlo= tru->GetBinLowEdge(j), thi= tlo+tw;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
          fv= _func->Integral (tlo, thi, _eps) / tw;
        } else {
          Double_t tlo[3]= { tru->GetXaxis()->GetBinLowEdge(jx), tru->GetYaxis()->GetBinLowEdge(jy), tru->GetZaxis()->GetBinLowEdge(jz) };
          Double_t thi[3]= { tru->GetXaxis()->GetBinUpEdge (jx), tru->GetYaxis()->GetBinUpEdge (jy), tru->GetZaxis()->GetBinUpEdge (jz) };
          Double_t relerr=0.0;
          fv= _func->IntegralMultiple (_ndim, tlo, thi, _eps, relerr);
          fv /= tru->GetXaxis()->GetBinWidth(jx) * tru->GetYaxis()->GetBinWidth(jy);
          if (_ndim>=3) fv /= tru->GetZaxis()->GetBinWidth(jz);
#else
          fv= _func->Integral (tlo, thi, (Double_t*)0, _eps) / tw;
        } else if (_ndim==2) {
          fv= _func->Integral (tru->GetXaxis()->GetBinLowEdge(jx), tru->GetXaxis()->GetBinUpEdge(jx),
                               tru->GetYaxis()->GetBinLowEdge(jy), tru->GetYaxis()->GetBinUpEdge(jy), _eps);
          fv /= tru->GetXaxis()->GetBinWidth(jx) * tru->GetYaxis()->GetBinWidth(jy);
        } else {
          fv= _func->Integral (tru->GetXaxis()->GetBinLowEdge(jx), tru->GetXaxis()->GetBinUpEdge(jx),
                               tru->GetYaxis()->GetBinLowEdge(jy), tru->GetYaxis()->GetBinUpEdge(jy),
                               tru->GetZaxis()->GetBinLowEdge(jz), tru->GetZaxis()->GetBinUpEdge(jz), _eps);
          fv /= tru->GetXaxis()->GetBinWidth(jx) * tru->GetYaxis()->GetBinWidth(jy) * tru->GetZaxis()->GetBinWidth(jz);
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
                                      const char* name, const char* title)
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
  Setup (measured, truth, response);
}

template <class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT (const Hist* measured, const Hist* truth,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  // RooUnfoldResponseT<class Hist, class Hist2D> constructor - measured and truth only used for shape
  Init();
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
RooUnfoldResponseT<Hist,Hist2D>::Add (const RooUnfoldResponseT<Hist,Hist2D>& rhs)
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
  _mes->Add (rhs._mes);
  _fak->Add (rhs._fak);
  _tru->Add (rhs._tru);
  _res->Add (rhs._res);
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
  return Long64_t(_res->GetEntries());
}


template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Reset()
{
  // Resets object to initial state.
  ClearCache();
  delete _mes;
  delete _fak;
  delete _tru;
  delete _res;
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
  Bool_t oldstat= Hist::AddDirectoryStatus();
  Hist::AddDirectory (kFALSE);
  _mes= new TH1F ("measured", "Measured", nm, mlo, mhi);
  _fak= new TH1F ("fakes",    "Fakes",    nm, mlo, mhi);
  _tru= new TH1F ("truth",    "Truth",    nt, tlo, thi);
  _mdim= _tdim= 1;
  _nm= nm;
  _nt= nt;
  SetNameTitleDefault ("response", "Response");
  _res= new TH2D (GetName(), GetTitle(), nm, mlo, mhi, nt, tlo, thi);
  Hist::AddDirectory (oldstat);
  return *this;
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::ReplaceAxis (TAxis* dest, const TAxis* source)
{
  // Replaces an axis with that of a different histogram
  TString  name= dest->GetName();
  TObject* hist= dest->GetParent();
  source->Copy (*dest);
  dest->SetName   (name);
  dest->SetParent (hist);
}

template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Setup (const Hist* measured, const Hist* truth)
{
  // set up - measured and truth only used for shape
  Reset();
  Bool_t oldstat= Hist::AddDirectoryStatus();
  Hist::AddDirectory (kFALSE);
  _mes= (Hist*) measured ->Clone();
  _mes->Reset();
  _fak= (Hist*) _mes     ->Clone("fakes");
  _fak->SetTitle("Fakes");
  _tru= (Hist*) truth    ->Clone();
  _tru->Reset();
  _mdim= _mes->GetDimension();
  _tdim= _tru->GetDimension();
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  SetNameTitleDefault();
  _nm= _mes->GetNbinsX() * _mes->GetNbinsY() * _mes->GetNbinsZ();
  _nt= _tru->GetNbinsX() * _tru->GetNbinsY() * _tru->GetNbinsZ();
  _res= new TH2D (GetName(), GetTitle(), _nm, 0.0, Double_t(_nm), _nt, 0.0, Double_t(_nt));
  if (_mdim==1) ReplaceAxis (_res->GetXaxis(), _mes->GetXaxis());
  if (_tdim==1) ReplaceAxis (_res->GetYaxis(), _tru->GetXaxis());
  Hist::AddDirectory (oldstat);
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
  Bool_t oldstat= Hist::AddDirectoryStatus();
  Hist::AddDirectory (kFALSE);
  _res= (Hist2D*) response->Clone();
  if (measured) {
    _mes= (Hist*) measured->Clone();
    _fak= (Hist*) measured->Clone("fakes");
    _fak->Reset();
    _fak->SetTitle("Fakes");
    _mdim= _mes->GetDimension();
  } else {
    _mes= new TH1F ("measured", "Measured", response->GetNbinsX(), 0.0, 1.0);
    ReplaceAxis (_mes->GetXaxis(), _res->GetXaxis());
    _fak= (Hist*) _mes->Clone("fakes");
    _fak->SetTitle("Fakes");
    _mdim= 1;
  }
  if (truth) {
    _tru= (Hist*)  truth   ->Clone();
    _tdim= _tru->GetDimension();
  } else {
    _tru= new TH1F ("truth",    "Truth",    response->GetNbinsY(), 0.0, 1.0);
    ReplaceAxis (_tru->GetXaxis(), _res->GetYaxis());
    _tdim= 1;
  }
  Hist::AddDirectory (oldstat);
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  _nm= _mes->GetNbinsX() * _mes->GetNbinsY() * _mes->GetNbinsZ();
  _nt= _tru->GetNbinsX() * _tru->GetNbinsY() * _tru->GetNbinsZ();
  if (_nm != _res->GetNbinsX() || _nt != _res->GetNbinsY()) {
    cerr << "Warning: RooUnfoldResponseT<class Hist, class Hist2D> measured X truth is " << _nm << " X " << _nt
         << ", but matrix is " << _res->GetNbinsX()<< " X " << _res->GetNbinsY() << endl;
  }

  Int_t first=1, nm= _nm, nt= _nt, s= _res->GetSumw2N();
  if (_overflow) {
    first= 0;
    nm += 2;
    nt += 2;
  }

  if (!measured || _mes->GetEntries() == 0.0) {
    // Similar to _res->ProjectionX() but without stupid reset of existing histograms
    // Always include under/overflows in sum of truth.
    for (Int_t i= 0; i<nm; i++) {
      Double_t nmes= 0.0, wmes= 0.0;
      for (Int_t j= 0; j<_nt+2; j++) {
               nmes +=      _res->GetBinContent (i+first, j);
        if (s) wmes += pow (_res->GetBinError   (i+first, j), 2);
      }
      Int_t bin= GetBin (_mes, i, _overflow);
             _mes->SetBinContent (bin,      nmes );
      if (s) _mes->SetBinError   (bin, sqrt(wmes));
    }
  } else {
    // Fill fakes from the difference of _mes - _res->ProjectionX()
    // Always include under/overflows in sum of truth.
    Int_t sm= _mes->GetSumw2N(), nfake=0;
    for (Int_t i= 0; i<nm; i++) {
      Double_t nmes= 0.0, wmes= 0.0;
      for (Int_t j= 0; j<_nt+2; j++) {
               nmes +=      _res->GetBinContent (i+first, j);
        if (s) wmes += pow (_res->GetBinError   (i+first, j), 2);
      }
      Int_t bin= GetBin (_mes, i, _overflow);
      Double_t fake= _mes->GetBinContent (bin) - nmes;
      if (fake!=0.0) nfake++;
      if (!s) wmes= nmes;
      _fak->SetBinContent (bin, fake);
      _fak->SetBinError   (bin, sqrt (wmes + (sm ? pow(_mes->GetBinError(bin),2) : _mes->GetBinContent(bin))));
    }
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,13,0)
    _fak->SetEntries (_fak->GetEffectiveEntries());  // 0 entries if 0 fakes
#else
    _fak->SetEntries (nfake);  // 0 entries if 0 fakes
#endif
  }

  if (!truth || _tru->GetEntries() == 0.0) {
    // similar to _res->ProjectionY() but without stupid reset of existing histograms
    // Always include under/overflows in sum of measurements.
    for (Int_t j= 0; j<nt; j++) {
      Double_t ntru= 0.0, wtru= 0.0;
      for (Int_t i= 0; i<_nm+2; i++) {
               ntru +=      _res->GetBinContent (i, j+first);
        if (s) wtru += pow (_res->GetBinError   (i, j+first), 2);
      }
      Int_t bin= GetBin (_tru, j, _overflow);
             _tru->SetBinContent (bin,      ntru);
      if (s) _tru->SetBinError   (bin, sqrt(wtru));
    }
  }

  SetNameTitleDefault();
  return *this;
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::ClearCache()
{
  delete _vMes; _vMes= 0;
  delete _eMes; _eMes= 0;
  delete _vFak; _vFak= 0;
  delete _vTru; _vTru= 0;
  delete _eTru; _eTru= 0;
  delete _mRes; _mRes= 0;
  delete _eRes; _eRes= 0;
  _cached= false;
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fill (Double_t xr, Double_t xt, Double_t w)
{
  // Fill 1D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==1 && _tdim==1);
  if (_cached) ClearCache();
  _mes->Fill (xr, w);
  _tru->Fill (xt, w);
  return _res->Fill (xr, xt, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fill (Double_t xr, Double_t yr, Double_t xt, Double_t yt, Double_t w)
{
  // Fill 2D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==2 && _tdim==2);
  if (_cached) ClearCache();
  ((Hist2D*)_mes)->Fill (xr, yr, w);
  ((Hist2D*)_tru)->Fill (xt, yt, w);
  return _res->Fill (_res->GetXaxis()->GetBinCenter (FindBin (_mes, xr, yr)+1),
                     _res->GetYaxis()->GetBinCenter (FindBin (_tru, xt, yt)+1), w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fill (Double_t xr, Double_t yr, Double_t zr, Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  // Fill 3D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==3 && _tdim==3);
  if (_cached) ClearCache();
  ((TH3*)_mes)->Fill (xr, yr, zr, w);
  ((TH3*)_tru)->Fill (xt, yt, zt, w);
  return _res->Fill (_res->GetXaxis()->GetBinCenter (FindBin (_mes, xr, yr, zr)+1),
                     _res->GetYaxis()->GetBinCenter (FindBin (_tru, xt, yt, zt)+1), w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::FindBin(const Hist* h, Double_t x, Double_t y)
{
  // Get vector index (0..nx*ny-1) for bin containing (x,y) coordinates
  Int_t nx=   h->GetNbinsX();
  Int_t ny=   h->GetNbinsY();
  Int_t binx= h->GetXaxis()->FindBin(x) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny;
  Int_t biny= h->GetYaxis()->FindBin(y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny;
  return binx + nx*biny;
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::FindBin(const Hist* h, Double_t x, Double_t y, Double_t z)
{
  // Get vector index (0..nx*ny*nz-1) for bin containing (x,y,z) coordinates
  Int_t nx=   h->GetNbinsX();
  Int_t ny=   h->GetNbinsY();
  Int_t nz=   h->GetNbinsZ();
  Int_t binx= h->GetXaxis()->FindBin(x) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny*nz;
  Int_t biny= h->GetYaxis()->FindBin(y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny*nz;
  Int_t binz= h->GetZaxis()->FindBin(z) - 1;
  if (binz <  0)  return -1;
  if (binz >= nz) return nx*ny*nz;
  return binx + nx*(biny + ny*binz);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::GetBinDim (const Hist* h, Int_t i)
{
  // Converts from vector index (0..nx*ny-1) or (0..nx*ny*nz-1) to multi-dimensional histogram
  // global bin number (0..(nx+2)*(ny+2)-1) or (0..(nx+2)*(ny+2)*(nz+2)-1), skipping under/overflow bins.
  Int_t ndim= h->GetDimension(), nx= h->GetNbinsX();
  if        (ndim == 2) {
//  cout << i << " -> " << "(" << i%nx+1 << "," << i/nx+1 << ")" << endl;
    return (i%nx+1) + (nx+2)*(i/nx+1);
  } else if (ndim == 3) {
    Int_t ny= h->GetNbinsY();
//  cout << i << " -> " << "(" << i%nx+1 << "," << (i/nx)%ny+1 << "," << i/(nx*ny)+1 << ")" << endl;
    return (i%nx+1) + (nx+2)*((i/nx)%ny+1 + (ny+2)*(i/(nx*ny)+1));
  }
  return i+1;   // not used: 1D handled by inline GetBin() (and handling UseOverflow), don't support >3D.
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Miss1D (Double_t xt, Double_t w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 1D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_tdim==1);
  if (_cached) ClearCache();
  return _tru->Fill (xt, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Miss2D (Double_t xt, Double_t yt, Double_t w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 2D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_tdim==2);
  if (_cached) ClearCache();
  return ((Hist2D*)_tru)->Fill (xt, yt, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Miss (Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  // Fill missed event (not reconstructed due to detection inefficiencies) into 3D Response Matrix
  assert (_tru != 0);
  assert (_tdim==3);
  if (_cached) ClearCache();
  return ((TH3*)_tru)->Fill (xt, yt, zt, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fake1D (Double_t xr, Double_t w)
{
  // Fill fake event (reconstructed event with no truth) into 1D Response Matrix (with weight)
  assert (_fak != 0 && _mes != 0);
  assert (_mdim==1);
  if (_cached) ClearCache();
         _mes->Fill (xr, w);
  return _fak->Fill (xr, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fake2D (Double_t xr, Double_t yr, Double_t w)
{
  // Fill fake event (reconstructed event with no truth) into 2D Response Matrix (with weight)
  assert (_mes != 0);
  assert (_mdim==2);
  if (_cached) ClearCache();
         ((Hist2D*)_fak)->Fill (xr, yr, w);
  return ((Hist2D*)_mes)->Fill (xr, yr, w);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::Fake (Double_t xr, Double_t yr, Double_t zr, Double_t w)
{
  // Fill fake event (reconstructed event with no truth) into 3D Response Matrix
  assert (_mes != 0);
  assert (_mdim==3);
  if (_cached) ClearCache();
         ((TH3*)_mes)->Fill (xr, yr, zr, w);
  return ((TH3*)_fak)->Fill (xr, yr, zr, w);
}

template <class Hist, class Hist2D> Hist*
RooUnfoldResponseT<Hist,Hist2D>::H2H1D(const Hist2D* h, Int_t nb)
{
  Hist* h1d= new TH1F(h->GetName(), h->GetTitle(), nb, 0.0, 1.0);
  Int_t s= h->GetSumw2N();
  for (Int_t i= 0; i < nb; i++) {
    Int_t j= GetBin (h, i);  // don't bother with under/overflow bins (not supported for >1D)
           h1d->SetBinContent (i+1, h->GetBinContent (j));
    if (s) h1d->SetBinError   (i+1, h->GetBinError   (j));
  }
  return h1d;
}

template <class Hist, class Hist2D> Hist*
RooUnfoldResponseT<Hist,Hist2D>::H2H1D(const Hist* h, Int_t nb)
{
  if (dynamic_cast<const Hist*>(h)) return dynamic_cast<Hist*>(h->Clone());
  return H2H1D((Hist2D*)h,nb);
}

template <class Hist, class Hist2D> Hist2D*
RooUnfoldResponseT<Hist,Hist2D>::HresponseNoOverflow() const
{
  const Hist2D* h= Hresponse();
  Int_t nx= h->GetNbinsX(), ny= h->GetNbinsY(), s= h->GetSumw2N();
  if (_overflow) {  // implies truth/measured both 1D
    Double_t xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax(), xb= (xhi-xlo)/nx;
    Double_t ylo= h->GetYaxis()->GetXmin(), yhi= h->GetYaxis()->GetXmax(), yb= (yhi-ylo)/ny;
    nx += 2; ny += 2;
    Hist2D* hx= new TH2D (h->GetName(), h->GetTitle(), nx, xlo-xb, xhi+xb, ny, ylo-yb, yhi+yb);
    for (Int_t i= 0; i < nx; i++) {
      for (Int_t j= 0; j < ny; j++) {
               hx->SetBinContent (i+1, j+1, h->GetBinContent (i, j));
        if (s) hx->SetBinError   (i+1, j+1, h->GetBinError   (i, j));
      }
    }
    return hx;
  } else if (dynamic_cast<const TH2D*>(h)) {
    Hist2D* hx= dynamic_cast<TH2D*>(h->Clone());
    // clear under/overflows
    for (Int_t i= 0; i <= nx+1; i++) {
      hx->SetBinContent (i, 0,    0.0);
      hx->SetBinContent (i, ny+1, 0.0);
    }
    for (Int_t i= 1; i <= ny;   i++) {
      hx->SetBinContent (0,    i, 0.0);
      hx->SetBinContent (nx+1, i, 0.0);
    }
    return hx;
  } else {
    Double_t xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax();
    Double_t ylo= h->GetYaxis()->GetXmin(), yhi= h->GetYaxis()->GetXmax();
    Hist2D* hx= new TH2D (h->GetName(), h->GetTitle(), nx, xlo, xhi, ny, ylo, yhi);
    for (Int_t i= 0; i < nx+2; i++) {
      for (Int_t j= 0; j < ny+2; j++) {
               hx->SetBinContent (i, j, h->GetBinContent (i, j));
        if (s) hx->SetBinError   (i, j, h->GetBinError   (i, j));
      }
    }
    return hx;
  }
}

template <class Hist, class Hist2D> TVectorD*
RooUnfoldResponseT<Hist,Hist2D>::H2V  (const Hist* h, Int_t nb, Bool_t overflow)
{
  // Returns TVectorD of the bin contents of the input histogram
  if (overflow) nb += 2;
  TVectorD* v= new TVectorD (nb);
  if (!h) return v;
  for (Int_t i= 0; i < nb; i++) {
    (*v)(i)= GetBinContent (h, i, overflow);
  }
  return v;
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::V2H (const TVectorD& v, Hist* h, Int_t nb, Bool_t overflow)
{
  // Sets the bin content of the histogram as that element of the input vector
  h->Reset();  // in particular, ensure under/overflows are reset
  if (overflow) nb += 2;
  for (Int_t i= 0; i < nb; i++) {
    Int_t j= GetBin (h, i, overflow);
    h->SetBinContent (j, v(i));
  }
}

template <class Hist, class Hist2D> TVectorD*
RooUnfoldResponseT<Hist,Hist2D>::H2VE (const Hist* h, Int_t nb, Bool_t overflow)
{
  // Returns TVectorD of bin errors for input histogram
  if (overflow) nb += 2;
  TVectorD* v= new TVectorD (nb);
  if (!h) return v;
  for (Int_t i= 0; i < nb; i++) {
    (*v)(i)= GetBinError (h, i, overflow);
  }
  return v;
}

template <class Hist, class Hist2D> TMatrixD*
RooUnfoldResponseT<Hist,Hist2D>::H2M  (const Hist2D* h, Int_t nx, Int_t ny, const Hist* norm, Bool_t overflow)
{
  // Returns Matrix of values of bins in a 2D input histogram
  Int_t first= overflow ? 0 : 1;
  if (overflow) {
    nx += 2;
    ny += 2;
  }
  TMatrixD* m= new TMatrixD (nx, ny);
  if (!h) return m;
  for (Int_t j= 0; j < ny; j++) {
    Double_t fac;
    if (!norm) fac= 1.0;
    else {
      fac= GetBinContent (norm, j, overflow);
      if (fac != 0.0) fac= 1.0/fac;
    }
    for (Int_t i= 0; i < nx; i++) {
      (*m)(i,j)= h->GetBinContent(i+first,j+first) * fac;
    }
  }
  return m;
}

template <class Hist, class Hist2D> TMatrixD*
RooUnfoldResponseT<Hist,Hist2D>::H2ME (const Hist2D* h, Int_t nx, Int_t ny, const Hist* norm, Bool_t overflow)
{
  // Returns matrix of bin errors for a 2D histogram.
  Int_t first= overflow ? 0 : 1;
  if (overflow) {
    nx += 2;
    ny += 2;
  }
  TMatrixD* m= new TMatrixD (nx, ny);
  if (!h) return m;
  for (Int_t j= 0; j < ny; j++) {
    Double_t fac;
    if (!norm) fac= 1.0;
    else {
      fac= GetBinContent (norm, j, overflow);
      if (fac != 0.0) fac= 1.0/fac;
    }
    for (Int_t i= 0; i < nx; i++) {
      // Assume Poisson norm, Multinomial P(mes|tru)
      (*m)(i,j)= h->GetBinError(i+first,j+first) * fac;
    }
  }
  return m;
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::PrintMatrix(const TMatrixD& m, const char* name, const char* format, Int_t cols_per_sheet)
{
   // Print the matrix as a table of elements.
   // Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
   // By default the format "%11.4g" is used to print one element.
   // One can specify an alternative format with eg
   //  format ="%6.2f  "

   if (!m.IsValid()) {
     m.Error("PrintMatrix","%s is invalid",name);
     return;
   }

   const Int_t ncols  = m.GetNcols();
   const Int_t nrows  = m.GetNrows();
   const Int_t collwb = m.GetColLwb();
   const Int_t rowlwb = m.GetRowLwb();

   if (!(format && format[0])) format= "%11.4g ";
   char topbar[1000];
   snprintf(topbar,1000,format,123.456789);
   Int_t nch = strlen(topbar)+1;
   if (nch > 18) nch = 18;
   char ftopbar[20];
   for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
   Int_t nk = 1 + Int_t(log10(ncols));
   snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
   Int_t nch2 = strlen(ftopbar);
   for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
   ftopbar[nch] = '|';
   ftopbar[nch+1] = 0;

   printf("\n%dx%d %s is as follows",nrows,ncols,name);

   if (cols_per_sheet <= 0) {
     cols_per_sheet = 5;
     if (nch <= 8) cols_per_sheet =10;
   }
   nk = 5+nch*(cols_per_sheet<ncols ? cols_per_sheet : ncols);
   for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
   topbar[nk] = 0;
   for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
      printf("\n\n     |");
      for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
         printf(ftopbar,j+collwb-1);
      printf("\n%s\n",topbar);
      if (m.GetNoElements() <= 0) continue;
      for (Int_t i = 1; i <= nrows; i++) {
         printf("%4d |",i+rowlwb-1);
         for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
            printf(format,m(i+rowlwb-1,j+collwb-1));
         printf("\n");
      }
   }
   printf("\n");
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::Print (Option_t* /* option */) const
{
  PrintMatrix (Mresponse(), Form("%s response matrix",GetTitle()));
}


template <class Hist, class Hist2D> Hist*
RooUnfoldResponseT<Hist,Hist2D>::ApplyToTruth (const Hist* truth, const char* name) const
{
  // Apply the response matrix to the truth
  // Errors not set, since we assume original truth has no errors
  if (!Htruth()) return 0;  // Needed for checking binning if nothing else

  // If no truth histogram input, use training truth
  // If truth histogram input, make sure its binning is correct
  TVectorD* resultvect;
  if (truth) {
    if (truth->GetNbinsX() != _tru->GetNbinsX() ||
        truth->GetNbinsY() != _tru->GetNbinsY() ||
        truth->GetNbinsZ() != _tru->GetNbinsZ())
      cerr << "Warning: RooUnfoldResponseT<Hist,Hist2D>::ApplyToTruth truth histogram is a different size ("
           << (truth->GetNbinsX() * truth->GetNbinsY() * truth->GetNbinsZ()) << " bins) or shape from response matrix truth ("
           << ( _tru->GetNbinsX() *  _tru->GetNbinsY() *  _tru->GetNbinsZ()) << " bins)" << endl;
    resultvect= H2V (truth, GetNbinsTruth(), _overflow);
    if (!resultvect) return 0;
  } else {
    resultvect= new TVectorD (Vtruth());
  }

  (*resultvect) *= Mresponse();   // v= A*v

  // Turn results vector into properly binned histogram
  Hist* result= (Hist*) Hmeasured()->Clone (name);
  result->SetTitle (name);
  V2H (*resultvect, result, GetNbinsMeasured(), _overflow);
  delete resultvect;
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
  // Returns new RooUnfoldResponseT<class Hist, class Hist2D> object with smeared response matrix elements for use as a toy.
  TString name= GetName();
  name += "_toy";
  RooUnfoldResponseT<Hist,Hist2D>* res= new RooUnfoldResponseT<Hist,Hist2D> (*this);
  res->SetName(name);
  if (!FakeEntries()) _fak->Reset();
  Hist2D* hres= res->Hresponse();
  for (Int_t i= 1; i<=_nm; i++) {
    for (Int_t j= 1; j<=_nt; j++) {
      Int_t bin= hres->GetBin (i,j);
      Double_t e= hres->GetBinError (bin);
      if (e>0.0) {
        Double_t v= hres->GetBinContent(bin) + gRandom->Gaus(0.0,e);
        if (v<0.0) v= 0.0;
        hres->SetBinContent (bin, v);
      }
    }
  }
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

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::Streamer (TBuffer &R__b)
{
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    Bool_t oldstat= Hist::AddDirectoryStatus();
    Hist::AddDirectory (kFALSE);
    RooUnfoldResponseT<Hist,Hist2D>::Class()->ReadBuffer  (R__b, this);
    Hist::AddDirectory (oldstat);
  } else {
    RooUnfoldResponseT<Hist,Hist2D>::Class()->WriteBuffer (R__b, this);
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
  if (!_vMes) _cached= (_vMes= H2V  (_mes, _nm, _overflow));
  return *_vMes;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vfakes() const
{
  // Fakes distribution as a TVectorD
  if (!_vFak) _cached= (_vFak= H2V  (_fak, _nm, _overflow));
  return *_vFak;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Emeasured() const
{
  // Measured distribution errors as a TVectorD
  if (!_eMes) _cached= (_eMes= H2VE (_mes, _nm, _overflow));
  return *_eMes;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vtruth() const
{
  // Truth distribution as a TVectorD
  if (!_vTru) _cached= (_vTru= H2V  (_tru, _nt, _overflow)); 
  return *_vTru;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Etruth() const
{
  // Truth distribution errors as a TVectorD
  if (!_eTru) _cached= (_eTru= H2VE (_tru, _nt, _overflow)); 
  return *_eTru;
}

template<class Hist, class Hist2D>
const TMatrixD& RooUnfoldResponseT<Hist,Hist2D>::Mresponse() const
{
  // Response matrix as a TMatrixD: (row,column)=(measured,truth)
  if (!_mRes) _cached= (_mRes= H2M  (_res, _nm, _nt, _tru, _overflow)); 
  return *_mRes;
}

template<class Hist, class Hist2D>
const TMatrixD& RooUnfoldResponseT<Hist,Hist2D>::Eresponse() const
{
  // Response matrix errors as a TMatrixD: (row,column)=(measured,truth)
  if (!_eRes) _cached= (_eRes= H2ME (_res, _nm, _nt, _tru, _overflow)); 
  return *_eRes;
}


template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::operator() (Int_t r, Int_t t) const
{
  // Response matrix element (measured,truth)
  return Mresponse()(r,t);
}

template<class Hist, class Hist2D>
Int_t    RooUnfoldResponseT<Hist,Hist2D>::GetBin (const Hist* h, Int_t i, Bool_t overflow)
{
  // vector index (0..nx*ny-1) -> multi-dimensional histogram
  // global bin number (0..(nx+2)*(ny+2)-1) skipping under/overflow bins
  return (h->GetDimension()<2) ? i+(overflow ? 0 : 1) : GetBinDim(h,i);
}

template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::GetBinContent (const Hist* h, Int_t i, Bool_t overflow)
{
  // Bin content by vector index
  return h->GetBinContent (GetBin (h, i, overflow));
}

template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::GetBinError   (const Hist* h, Int_t i, Bool_t overflow)
{
  // Bin error   by vector index
  return h->GetBinError   (GetBin (h, i, overflow));
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
  return _fak ? _fak->GetEntries() : 0.0;
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::FindBin (const Hist* h, Double_t x)
{
  // Return vector index for bin containing (x)
  return h->GetXaxis()->FindBin(x) - 1;
}

template class RooUnfoldResponseT<TH1,TH2>;
ClassImp (RooUnfoldResponse);
