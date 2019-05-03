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

#include "RooHistFunc.h"
#include "RooDataHist.h"


#if ROOT_VERSION_CODE >= ROOT_VERSION(5,18,0)
#define HAVE_RooUnfoldFoldingFunction
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::pow;
using std::sqrt;


namespace {
  int findBinX(const TH1* h, double x){
    return h->GetXaxis()->FindBin(x);
  }
  int findBinY(const TH1* h, double y){
    return h->GetYaxis()->FindBin(y);
  }
  int findBinZ(const TH1* h, double z){
    return h->GetZaxis()->FindBin(z);
  }
  int findBinX(const RooAbsReal* h, double x){
    // TODO
    return 0;
  }
  int findBinY(const RooAbsReal* h, double y){
    // TODO
    return 0;
  }
  int findBinZ(const RooAbsReal* h, double z){
    // TODO
    return 0;
  }
  double xMin(const TH1* hist){
    return hist->GetXaxis()->GetXmin();
  }
  double xMax(const TH1* hist){
    return hist->GetXaxis()->GetXmax();
  }  
  double yMin(const TH1* hist){
    return hist->GetYaxis()->GetXmin();
  }
  double yMax(const TH1* hist){
    return hist->GetYaxis()->GetXmax();
  }  
  double zMin(const TH1* hist){
    return hist->GetZaxis()->GetXmin();
  }
  double zMax(const TH1* hist){
    return hist->GetZaxis()->GetXmax();
  }
  double xMin(const RooAbsReal* hist){
    return 0;
  }
  double xMax(const RooAbsReal* hist){
    return 0;
  }  
  double yMin(const RooAbsReal* hist){
    return 0;
  }
  double yMax(const RooAbsReal* hist){
    return 0;
  }  
  double zMin(const RooAbsReal* hist){
    return 0;
  }
  double zMax(const RooAbsReal* hist){
    return 0;
  }    
  
  int sumW2N(const TH1* hist){
    return hist->GetSumw2N();
  }
  int sumW2N(const RooAbsReal* hist){
    return 0;
  }

  
  void add(TH1* hista, const TH1* histb){
    hista->Add(histb);
  }
  void add(RooAbsReal* hista, RooAbsReal* histb){
    // TODO
  }  
  void projectY(RooAbsReal* _res, RooAbsReal* _tru, bool overflow){
    // TODO
  } 
  void projectY(TH2* _res, TH1* _tru, bool overflow){
    Int_t s= _res->GetSumw2N();
    for (Int_t j= 1-overflow; j<_res->GetNbinsY()+1+overflow; j++) {
      Double_t ntru= 0.0, wtru= 0.0;
      for (Int_t i= 0; i<_res->GetNbinsX()+2; i++) {
               ntru +=      _res->GetBinContent (i, j);
        if (s) wtru += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t bin= RooUnfoldResponseT<TH1,TH2>::GetBin (_tru, j, overflow);
             _tru->SetBinContent (bin,      ntru);
      if (s) _tru->SetBinError   (bin, sqrt(wtru));
    }
  }
  void projectX(RooAbsReal* _res, RooAbsReal* _mes, bool overflow){
    // TODO
  }  
  void projectX(TH2* _res, TH1* _mes, bool overflow){
    Int_t s= _res->GetSumw2N();
    for (Int_t i= 1-overflow; i<_res->GetNbinsX()+1+overflow; i++) {
      Double_t nmes= 0.0, wmes= 0.0;
      for (Int_t j= 0; j<_res->GetNbinsY(); j++) {
               nmes +=      _res->GetBinContent (i, j);
        if (s) wmes += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t bin= RooUnfoldResponseT<TH1,TH2>::GetBin (_mes, i, overflow);
             _mes->SetBinContent (bin,      nmes );
      if (s) _mes->SetBinError   (bin, sqrt(wmes));
    }
  }
  void subtractProjectX(RooAbsReal* _res, RooAbsReal* _mes, RooAbsReal* _fak, bool overflow){
    // TODO
  }
  
  void subtractProjectX(TH2* _res, TH1* _mes, TH1* _fak, bool overflow){
    Int_t s= _res->GetSumw2N();
    Int_t sm= _mes->GetSumw2N(), nfake=0;
    for (Int_t i= 1-overflow; i<_res->GetNbinsX()+1+overflow; i++) {
      Double_t nmes= 0.0, wmes= 0.0;
      for (Int_t j= 0; j<_res->GetNbinsY()+2; j++) {
               nmes +=      _res->GetBinContent (i, j);
        if (s) wmes += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t bin= RooUnfoldResponseT<TH1,TH2>::GetBin (_mes, i, overflow);
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

  int fill(TH1* hist, double x, double w){
    return hist->Fill (x, w);
  }
  int fill(TH1* hist, double x, double y, double w){
    return ((TH2*)hist)->Fill (x, y, w);
  }
  int fill(TH1* hist, double x, double y, double z, double w){
    return ((TH3*)hist)->Fill (x, y, z, w);
  }    
  int fill(TH2* hist, double x, double y, double w){
    return hist->Fill (x, y, w);
  }
  int fill(RooAbsReal* hist, double x, double w){
    return 0;
  }
  int fill(RooAbsReal* hist, double x, double y, double w){
    return 0;
  }  
  int fill(RooAbsReal* hist, double x, double y, double z, double w){
    return 0;
  }  

  
  TH1* copy(const TH1* orighist, bool reset, const char* name = 0, const char* title = 0){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH1* hist = (TH1*)(orighist ->Clone());
    if(name) hist->SetName(name);
    if(title) hist->SetTitle(title);
    if(reset) hist->Reset();
    return hist;
    TH1::AddDirectory (oldstat);
  }
  TH2* copy(const TH2* orighist, bool reset, const char* name = 0, const char* title = 0){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH2* hist = (TH2*)(orighist ->Clone());
    if(name) hist->SetName(name);
    if(title) hist->SetTitle(title);
    if(reset) hist->Reset();
    return hist;
    TH1::AddDirectory (oldstat);
  }  
  RooAbsReal* copy(const RooAbsReal* r, bool reset, const char* name = 0, const char* title = 0){
    RooAbsReal* retval = (RooAbsReal*)(r->clone());
    if(name) retval->SetName(name);
    if(title) retval->SetTitle(title);
    return retval;
  }
  int entries(const TH1* hist){
    return hist->GetEntries();
  }
  int entries(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  int dim(const TH1* hist){
    return hist->GetDimension();
  }
  int dim(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  int nBins(const TH1* hist){
    return hist->GetNbinsX() * hist->GetNbinsY() * hist->GetNbinsZ();
  }
  int nBinsX(const TH1* hist){
    return hist->GetNbinsX();
  }  
  int nBinsY(const TH1* hist){
    return hist->GetNbinsY();
  }
  int nBinsZ(const TH1* hist){
    return hist->GetNbinsZ();
  }  
  int nBins(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  int nBinsX(const RooAbsReal* hist){
    // TODO
    return 0;
  }  
  int nBinsY(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  int nBinsZ(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  template<class Hist2D> Hist2D* createHist(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, int nbinsy, double ymin, double ymax, const char* yname);
  template<>  
  TH2* createHist<TH2>(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, int nbinsy, double ymin, double ymax, const char* yname){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH2* hist = new TH2D (name,title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
    TH1::AddDirectory (oldstat);
    return hist;
  }
  template<>  
  RooAbsReal* createHist<RooAbsReal>(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, int nbinsy, double ymin, double ymax, const char* yname){
    RooRealVar* x = new RooRealVar(xname,xname,nbinsx,xmin,xmax);
    RooRealVar* y = new RooRealVar(yname,yname,nbinsy,ymin,ymax);
    RooArgSet vars(*x,*y);
    RooDataHist* hist = new RooDataHist (name,title,vars);
    return new RooHistFunc(name,title,vars,vars,*hist);
  }
  template<class Hist> Hist* createHist(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname);  
  template<> TH1* createHist<TH1>(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH1* hist = new TH1D (name,title, nbinsx, xmin, xmax);
    TH1::AddDirectory (oldstat);
    return hist;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname){
    RooRealVar* x = new RooRealVar(xname,xname,nbinsx,xmin,xmax);
    RooArgSet vars(*x);    
    RooDataHist* hist = new RooDataHist (name,title,vars);
    return new RooHistFunc(name,title,vars,vars,*hist);
  }
}




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
  add(_mes,rhs._mes);
  add(_fak,rhs._fak);
  add(_tru,rhs._tru);
  add(_res,rhs._res);
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
  if(_mes) delete _mes;
  if(_fak) delete _fak;
  if(_tru) delete _tru;
  if(_res) delete _res;
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
  _mdim= _tdim= 1;
  _nm= nm;
  _nt= nt;
  _mes= createHist<Hist>("measured", "Measured", nm, mlo, mhi,"xm");
  _fak= createHist<Hist>("fakes",    "Fakes",    nm, mlo, mhi,"xm");
  _tru= createHist<Hist>("truth",    "Truth",    nt, tlo, thi,"xt");
  _res= createHist<Hist2D>("response", "Response", nm, mlo, mhi, "xm", nt, tlo, thi, "xt");
  return *this;
}

template <class Hist, class Hist2D> RooUnfoldResponseT<Hist,Hist2D>&
RooUnfoldResponseT<Hist,Hist2D>::Setup (const Hist* measured, const Hist* truth)
{
  // set up - measured and truth only used for shape
  Reset();
  _mes= copy(measured,true);
  _fak= copy(measured,true,"fakes","Fakes");
  _tru= copy(truth,true,"truth");
  _mdim= dim(_mes);
  _tdim= dim(_tru);
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  SetNameTitleDefault();
  _nm= nBins(_mes);
  _nt= nBins(_tru);
  _res=createHist<Hist2D>(GetName(), GetTitle(), _nm, 0.0, Double_t(_nm), "xm", _nt, 0.0, Double_t(_nt), "xt");
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
  _res= copy(response,false);
  if (measured) {
    _mes= copy(measured,false);
    _fak= copy(measured,true,"fakes","Fakes");
    _mdim= dim(_mes);
  } else {
    _mes= createHist<Hist>("measured", "Measured", nBinsX(response), 0.0, 1.0, "xm");
    _fak= copy(_mes,false,"fakes","Fakes");
    _mdim= 1;
  }
  if (truth) {
    _tru= copy(truth,false);
    _tdim= dim(_tru);
  } else {
    _tru= createHist<Hist>("truth",    "Truth",    nBinsY(response), 0.0, 1.0, "xt");
    _tdim= 1;
  }
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  _nm= nBins(_mes);
  _nt= nBins(_tru);
  if (_nm != nBinsX(_res) || _nt != nBinsY(_res)) {
    cerr << "Warning: RooUnfoldResponseT<class Hist, class Hist2D> measured X truth is " << _nm << " X " << _nt
         << ", but matrix is " << nBinsX(_res)<< " X " << nBinsY(_res) << endl;
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
  return fill(_res,GetBinCenterX(_res,FindBin (_mes, xr, yr)+1),GetBinCenterY(_res,FindBin (_tru, xt, yt)+1), w);
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
  return fill(_res,GetBinCenterX(_res,FindBin (_mes, xr, yr, zr)+1),GetBinCenterY(_res,FindBin (_tru, xt, yt, zt)+1), w);  
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::FindBin(const Hist* h, Double_t x, Double_t y)
{
  // Get vector index (0..nx*ny-1) for bin containing (x,y) coordinates
  Int_t nx=   nBinsX(h);
  Int_t ny=   nBinsY(h);
  Int_t binx= findBinX(h,x) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny;
  Int_t biny= findBinY(h,y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny;
  return binx + nx*biny;
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::FindBin(const Hist* h, Double_t x, Double_t y, Double_t z)
{
  // Get vector index (0..nx*ny*nz-1) for bin containing (x,y,z) coordinates
  Int_t nx=   nBinsX(h);
  Int_t ny=   nBinsY(h);
  Int_t nz=   nBinsZ(h);
  Int_t binx= findBinX(h,x) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny*nz;
  Int_t biny= findBinY(h,y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny*nz;
  Int_t binz= findBinZ(h,z) - 1;
  if (binz <  0)  return -1;
  if (binz >= nz) return nx*ny*nz;
  return binx + nx*(biny + ny*binz);
}

template <class Hist, class Hist2D> Int_t
RooUnfoldResponseT<Hist,Hist2D>::GetBinDim (const Hist* h, Int_t i)
{
  // Converts from vector index (0..nx*ny-1) or (0..nx*ny*nz-1) to multi-dimensional histogram
  // global bin number (0..(nx+2)*(ny+2)-1) or (0..(nx+2)*(ny+2)*(nz+2)-1), skipping under/overflow bins.
  Int_t ndim= dim(h), nx= nBinsX(h);
  if        (ndim == 2) {
//  cout << i << " -> " << "(" << i%nx+1 << "," << i/nx+1 << ")" << endl;
    return (i%nx+1) + (nx+2)*(i/nx+1);
  } else if (ndim == 3) {
    Int_t ny= nBinsY(h);
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
  Int_t nx= nBinsX(h), ny= nBinsY(h), s= sumW2N(h);
  if (_overflow) {  // implies truth/measured both 1D
    Double_t xlo= xMin(h), xhi= xMax(h), xb= (xhi-xlo)/nx;
    Double_t ylo= yMin(h), yhi= yMax(h), yb= (yhi-ylo)/ny;
    nx += 2; ny += 2;
    Hist2D* hx= createHist<Hist2D>(h->GetName(), h->GetTitle(), nx, xlo-xb, xhi+xb, "xm", ny, ylo-yb, yhi+yb, "xt");
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
  return (dim(h)<2) ? i+(overflow ? 0 : 1) : GetBinDim(h,i);
}

template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::GetBinContent (const Hist* h, Int_t i, Bool_t overflow)
{
  // Bin content by vector index
  return h->GetBinContent (GetBin (h, i, overflow));
}

template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::GetBinCenterX (const Hist2D* h, Int_t i){
  // Bin center by vector index
  return h->GetXaxis()->GetBinCenter(i);
}
template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::GetBinCenterY (const Hist2D* h, Int_t i){
  // Bin center by vector index
  return h->GetYaxis()->GetBinCenter(i);
}
template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::GetBinLowEdgeX (const Hist2D* h, Int_t i){
  // Bin center by vector index
  return h->GetXaxis()->GetBinLowEdge(i);
}
template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::GetBinLowEdgeY (const Hist2D* h, Int_t i){
  // Bin center by vector index
  return h->GetYaxis()->GetBinLowEdge(i);
}
template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::GetBinHighEdgeX (const Hist2D* h, Int_t i){
  // Bin center by vector index
  return h->GetXaxis()->GetBinLowEdge(i+1);
}
template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::GetBinHighEdgeY (const Hist2D* h, Int_t i){
  // Bin center by vector index
  return h->GetYaxis()->GetBinLowEdge(i+1);
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
  return _fak ? entries(_fak) : 0.0;
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::FindBin (const Hist* h, Double_t x)
{
  // Return vector index for bin containing (x)
  return h->GetXaxis()->FindBin(x) - 1;
}

template class RooUnfoldResponseT<TH1,TH2>;
ClassImp (RooUnfoldResponse);

//template class RooUnfoldResponseT<RooAbsReal,RooAbsReal>;
//ClassImp (RooAbsRealUnfoldResponse);
