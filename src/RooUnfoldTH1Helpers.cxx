#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"
#include "RooUnfoldHelpers.tpp"

#include <iostream>
#include <ostream>
#include <sstream>
#include <iomanip>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace RooUnfolding { 
  template<class Hist> int dim(const Hist* hist){
    return hist->GetDimension();
  }
}


namespace{
  const TAxis* getAxis(const TH1* h, RooUnfolding::Dimension d){
    if(d==RooUnfolding::X) return h->GetXaxis();
    if(d==RooUnfolding::Y) return h->GetYaxis();
    if(d==RooUnfolding::Z) return h->GetZaxis();
    throw std::runtime_error("invalid dimension passed!");
  }
  int binDim (const TH1* h, Int_t i)
  {
    // Converts from vector index (0..nx*ny-1) or (0..nx*ny*nz-1) to multi-dimensional histogram
    // global bin number (0..(nx+2)*(ny+2)-1) or (0..(nx+2)*(ny+2)*(nz+2)-1), skipping under/overflow bins.
    Int_t ndim= RooUnfolding::dim(h), nx= h->GetNbinsX();
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
}

namespace RooUnfolding { 
  template<> void reset<TH1>(TH1* h){
    h->Reset();
  }
  template<> int findBin<TH1>(const TH1* h, double x, RooUnfolding::Dimension d){
    return getAxis(h,d)->FindBin(x);
  }
  template<class Hist> double min(const Hist* hist, RooUnfolding::Dimension d){
    return getAxis(hist,d)->GetXmin();
  }
  template<class Hist> double max(const Hist* hist, RooUnfolding::Dimension d){
    return getAxis(hist,d)->GetXmax();
  }  
  template<> int nBins<TH1>(const TH1* hist, bool overflow){
    int d = dim(hist);
    if(d==1){
      return (hist->GetNbinsX()+2*overflow);
    } else if(d==2){
      return ( (hist->GetNbinsX()+2*overflow) * (hist->GetNbinsY()+2*overflow) );
    } else if(d==3){
      return ( (hist->GetNbinsX()+2*overflow) * (hist->GetNbinsY()+2*overflow) * (hist->GetNbinsZ()+2*overflow) );
    }
    return 0;
  }
  template<> int nBins<TH2>(const TH2* hist, bool overflow){
    return ( (hist->GetNbinsX()+2*overflow) * (hist->GetNbinsY()+2*overflow) );
  }
  template<> int nBins<TH3>(const TH3* hist, bool overflow){
    return ( (hist->GetNbinsX()+2*overflow) * (hist->GetNbinsY()+2*overflow) * (hist->GetNbinsZ()+2*overflow) );
  }
  template<class Hist> int nBins(const Hist* hist, RooUnfolding::Dimension d, bool overflow){
    const TAxis* ax = getAxis(hist,d);
    return ax->GetNbins()+2*overflow;
  }
  template<class Hist> const char* varname(const Hist* h, Dimension d){  
    return "";
  }
  template<class Hist> Variable<Hist> var(const Hist* h, Dimension d){
    return Variable<Hist>(nBins(h,d),min(h,d),max(h,d),"");
  }
  template<class Hist> int sumW2N(const Hist* hist){
    return hist->GetSumw2N();
  }
  template<class Hist> int entries(const Hist* hist){
    return hist->GetEntries();
  }
  template<class Hist> int bin(const Hist* h, Int_t i, Bool_t overflow){
    // vector index (0..nx*ny-1) -> multi-dimensional histogram
    // global bin number (0..(nx+2)*(ny+2)-1) skipping under/overflow bins
    return (dim(h)<2) ? i+(overflow ? 0 : 1) : binDim(h,i);
  }
  template<class Hist> int bin(const Hist* h, int i, int j, Bool_t overflow){
    return h->GetBin(i,j);
  }

  template<class Hist> double binCenter(const Hist*h, int i, RooUnfolding::Dimension d){
    const TAxis* ax = getAxis(h,d);
    return ax->GetBinCenter(i);
  }
  template<class Hist> double binWidth(const Hist*h, int i, RooUnfolding::Dimension d){
    const TAxis* ax = getAxis(h,d);
    return ax->GetBinWidth(i);
  }

  template<> double binHighEdge<TH1>(const TH1*h, int i, RooUnfolding::Dimension d){
    const TAxis* ax = getAxis(h,d);
    return ax->GetBinUpEdge(i);
  }
  template<> double binLowEdge<TH1>(const TH1*h, int i, RooUnfolding::Dimension d){
    const TAxis* ax = getAxis(h,d);
    return ax->GetBinLowEdge(i);
  }
  template<class Hist> void add(Hist* hista, const Hist* histb){
    hista->Add(histb);
  }

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
  template<class Hist> Hist* copy(const Hist* orighist, bool reset, const char* name, const char* title){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    Hist* hist = (Hist*)(orighist ->Clone());
    if(name) hist->SetName(name);
    if(title) hist->SetTitle(title);
    if(reset) hist->Reset();
    return hist;
    TH1::AddDirectory (oldstat);
  }
  template<> void binXYZ<TH1>(const TH1* tru, int i, int& jx, int& jy, int& jz){
    Int_t j= RooUnfolding::bin<TH1>(tru, i, false);
    if (dim(tru)>1) tru->GetBinXYZ (j, jx, jy, jz);
    else {
      jx = j;
      jy = 0;
      jz = 0;
    }
  }
  template<> double binError<TH1>(const TH1* h, Int_t i, Bool_t overflow)
  {
    // Bin error   by vector index
    return h->GetBinError   (bin (h, i, overflow));
  }
  template<> double binContent<TH1> (const TH1* h, Int_t i, Bool_t overflow)
  {
    // Bin content by vector index
    return h->GetBinContent (bin (h, i, overflow));
  }
  template<> double binContent<TH1> (const TH1* h, int i, int j, Bool_t overflow)
  {
    // Bin content by vector index
    return h->GetBinContent (bin (h, i, j, overflow));
  }
  template<> TH1* h2h1d<TH1,TH2>(const TH1* h, int nb){
    return dynamic_cast<TH1*>(h->Clone());
  }  
  template<> TH1* h2h1d<TH1,TH2>(const TH2* h, int nb){
    TH1* h1d= new TH1F(h->GetName(), h->GetTitle(), nb, 0.0, 1.0);
    Int_t s= h->GetSumw2N();
    for (Int_t i= 0; i < nb; i++) {
      Int_t j= bin (h, i, false);  // don't bother with under/overflow bins (not supported for >1D)
      h1d->SetBinContent (i+1, h->GetBinContent (j));
      if (s) h1d->SetBinError   (i+1, h->GetBinError   (j));
    }
    return h1d;
  }

  template<> void h2mNorm<TH1,TH2>  (const TH2* h, TMatrixD& m, const TH1* norm, bool overflow){
    // sets Matrix to values of bins in a 2D input histogram
    m.ResizeTo(h->GetNbinsX()+2*overflow,h->GetNbinsY()+2*overflow);
    for (Int_t j= 0; j < h->GetNbinsY()+2*overflow; ++j) {
      double fac = 1.;
      if (norm){
        fac= norm->GetBinContent(j);
        if (fac != 0.0) fac= 1.0/fac;
      }
      for (Int_t i= 0; i < h->GetNbinsX()+2*overflow; ++i) {
        m(i,j)= h->GetBinContent(i+!overflow,j+!overflow) * fac;
      }
    }
  }
  template<> void h2meNorm<TH1,TH2>  (const TH2* h, TMatrixD& m, const TH1* norm, bool overflow){
    // sets Matrix to values of bins in a 2D input histogram
    m.ResizeTo(h->GetNbinsX()+2*overflow,h->GetNbinsY()+2*overflow);
    for (Int_t j= 0; j < h->GetNbinsY()+2*overflow; ++j) {
      double fac = 1.;
      if (norm){
        fac= norm->GetBinContent(j);
        if (fac != 0.0) fac= 1.0/fac;
      }
      for (Int_t i= 0; i < h->GetNbinsX()+2*overflow; ++i) {
        m(i,j)= h->GetBinContent(i+!overflow,j+!overflow) * fac;
      }
    }
  }
  template<> TMatrixD h2mNorm<TH1,TH2>  (const TH2* h, const TH1* norm, bool overflow){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m(h->GetNbinsX()+2*overflow,h->GetNbinsY()+2*overflow);
    h2mNorm(h,m,norm,overflow);
    return m;
  }
  template<> TMatrixD h2meNorm<TH1,TH2>  (const TH2* h, const TH1* norm, bool overflow){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m(h->GetNbinsX()+2*overflow,h->GetNbinsY()+2*overflow);
    h2meNorm(h,m,norm,overflow);
    return m;
  }
  template<> void h2m  (const TH2* h, TMatrixD& m, bool overflow) { h2mNorm (h,m,(const TH1*)NULL,overflow); }
  template<> void h2me  (const TH2* h, TMatrixD& m, bool overflow){ h2meNorm(h,m,(const TH1*)NULL,overflow); };  

  template<> TMatrixD h2m<TH2>  (const TH2* h,bool overflow){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m(h->GetNbinsX()+2*overflow,h->GetNbinsY()+2*overflow);
    h2m(h,m,overflow);
    return m;
  }
  template<> TMatrixD h2me<TH2>  (const TH2* h,bool overflow){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m(h->GetNbinsX()+2*overflow,h->GetNbinsY()+2*overflow);
    h2me(h,m,overflow);
    return m;
  }

  template<> void h2v<TH1>  (const TH1* h, TVectorD& v, bool overflow){
    // sets Vector to values of bins in an input histogram
    int nbinstotal = nBins(h,true);
    v.ResizeTo(nBins(h,overflow));
    int n = 0;
    for (Int_t i= 0; i < nbinstotal; ++i){
      if(!overflow && (h->IsBinOverflow(i) || h->IsBinUnderflow(i))){
        continue;
      }
      v[n] = h->GetBinContent(i);
      ++n;
    }
  }
  template<> void h2ve<TH1>  (const TH1* h, TVectorD& v, bool overflow){
    // sets Vector to values of bins in an input histogram
    int nbinstotal = nBins(h,true);
    v.ResizeTo(nBins(h,overflow));
    int n = 0;
    for (Int_t i= 0; i < nbinstotal; ++i){
      if(!overflow && (h->IsBinOverflow(i) || h->IsBinUnderflow(i))){
        continue;
      }
      v[n] = h->GetBinError(i);
      ++n;
    }  
  }    
  template<> TVectorD h2v<TH1>  (const TH1* h, bool overflow){
    // Returns Vector of values of bins in an input histogram
    TVectorD v(nBins(h,overflow));
    h2v(h,v,overflow);
    return v;
  }
  template<> TVectorD h2ve<TH1>  (const TH1* h, bool overflow){
    // Returns Vector of values of bins in an input histogram
    TVectorD v(nBins(h,overflow));
    h2ve(h,v,overflow);
    return v;
  }
  
  template<> TH2* copyHistogram<TH2>(const TH2* h, bool includeOverflow){
    Int_t nx= nBins(h,RooUnfolding::X), ny= nBins(h,RooUnfolding::Y), s= h->GetSumw2N();
    if (includeOverflow) {  // implies truth/measured both 1D
      Double_t xlo= min(h,X), xhi= max(h,X), xb= (xhi-xlo)/nx;
      Double_t ylo= min(h,Y), yhi= max(h,Y), yb= (yhi-ylo)/ny;
      nx += 2; ny += 2;
      TH2* hx= new TH2D(h->GetName(), h->GetTitle(), nx, xlo-xb, xhi+xb, ny, ylo-yb, yhi+yb);
      for (Int_t i= 0; i < nx; i++) {
        for (Int_t j= 0; j < ny; j++) {
          hx->SetBinContent (i+1, j+1, h->GetBinContent (i, j));
          if (s) hx->SetBinError   (i+1, j+1, h->GetBinError   (i, j));
        }
      }
      return hx;
    } else if (dynamic_cast<const TH2D*>(h)) {
      TH2* hx= dynamic_cast<TH2*>(h->Clone());
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
      TH2* hx= new TH2D (h->GetName(), h->GetTitle(), nx, xlo, xhi, ny, ylo, yhi);
      for (Int_t i= 0; i < nx+2; i++) {
        for (Int_t j= 0; j < ny+2; j++) {
          hx->SetBinContent (i, j, h->GetBinContent (i, j));
          if (s) hx->SetBinError   (i, j, h->GetBinError   (i, j));
        }
      }
      return hx;
    }
  }
  template<> TH1* createHist<TH1>(const char* name, const char* title, const std::vector<Variable<TH1>>& x){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH1* hist = NULL;
    if(x.size() == 1){
      hist = new TH1D (name,title, x[0]._nBins, x[0]._min, x[0]._max);
    } else if(x.size() == 2){
      hist = new TH2D (name,title, x[0]._nBins, x[0]._min, x[0]._max, x[1]._nBins, x[1]._min, x[1]._max);
    } else if(x.size() == 3){
      hist = new TH3D (name,title, x[0]._nBins, x[0]._min, x[0]._max, x[1]._nBins, x[1]._min, x[1]._max, x[2]._nBins, x[2]._min, x[2]._max);
    } else {
      throw std::runtime_error(TString::Format("invalid dimensionality for ROOT histogram: %d",(int)x.size()).Data());
    }
    TH1::AddDirectory (oldstat);
    return hist;
  }
  template<> TH1* createHist<TH1>(const TVectorD& v, const char* name, const char* title, const std::vector<Variable<TH1>>& x, bool overflow){  
    // Sets the bin content of the histogram as that element of the input vector
    int nb = v.GetNrows();
    TH1* h = createHist<TH1>(name,title,x);
    for (Int_t i= 0; i < nb; i++) {
      Int_t j= RooUnfolding::bin<TH1>(h, i, overflow);
      h->SetBinContent (j, v(i));
    }
    return h;
  }
  template<> TH1* createHist<TH1>(const TVectorD& v, const TVectorD& ve, const char* name, const char* title, const std::vector<Variable<TH1>>& x, bool overflow){  
    // Sets the bin content of the histogram as that element of the input vector
    int nb = v.GetNrows();
    TH1* h = createHist<TH1>(name,title,x);
    if (overflow) nb += 2;
    for (Int_t i= 0; i < nb; i++) {
      Int_t j= RooUnfolding::bin<TH1>(h, i, overflow);
      h->SetBinContent (j, v(i));
      h->SetBinError (j, ve(i));
    }
    return h;
  }
  template<>  
  TH2* createHist<TH2>(const char* name, const char* title, const Variable<TH2>& x, const Variable<TH2>& y){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH2* hist = new TH2D (name,title, x._nBins, x._min, x._max, y._nBins, y._min, y._max);
    TH1::AddDirectory (oldstat);
    return hist;
  }
  template<> TH2* createHist<TH2>(const TMatrixD& m, const char* name, const char* title, const Variable<TH2>& x, const Variable<TH2>& y){  
    // Sets the bin content of the histogram as that element of the input vector
    int nx = m.GetNrows();
    int ny = m.GetNcols();
    TH2* h = createHist<TH2>(name,title,x,y);
    for (Int_t i= 0; i < nx; i++) {
      for (Int_t j= 0; j < nx; j++) {
        Int_t n= RooUnfolding::bin<TH1>(h, i, j, false);
        h->SetBinContent (n, m(i,j));
      }
    }
    return h;
  }
  template<> TH2* createHist<TH2>(const TMatrixD& m, const TMatrixD& me, const char* name, const char* title, const Variable<TH2>& x, const Variable<TH2>& y){  
    // Sets the bin content of the histogram as that element of the input vector
    int nx = m.GetNrows();
    int ny = m.GetNcols();
    TH2* h = createHist<TH2>(name,title,x,y);
    for (Int_t i= 0; i < nx; i++) {
      for (Int_t j= 0; j < nx; j++) {
        Int_t n= RooUnfolding::bin<TH1>(h, i, j, false);
        h->SetBinContent (n, m(i,j));
        h->SetBinError   (n, me(i,j));
      }
    }
    return h;
  }
  template<> TVectorD subtract<TH1,TVectorD>(const TVectorD& orig, const TH1* hist, bool overflow){
    const int n = nBins(hist);
    TVectorD res(n);
    for (Int_t i = 0 ; i < n; i++) {
      Int_t it= RooUnfolding::bin (hist, i, overflow);
      if (hist->GetBinContent(it)!=0.0 || hist->GetBinError(it)>0.0) {
        res[i] = orig[i] - hist->GetBinContent(it);
      }
    }
    return res;
  }

  template<> void printTable<TH1> (std::ostream& o, const TH1* hTrainTrue, const TH1* hTrain,
                   const TH1* hTrue, const TH1* hMeas, const TH1* hReco,
                   Int_t _nm, Int_t _nt, Bool_t _overflow,
                   ErrorTreatment withError, Double_t chi_squ)
  {
    // Prints entries from truth, measured, and reconstructed data for each bin.
    if (withError==kDefault) withError= hReco->GetSumw2N() ? kErrors : kNoError;
    if (_nm<=0) _nm= hTrain    ->GetNbinsX();
    if (_nt<=0) _nt= hTrainTrue->GetNbinsX();
    std::ostringstream fmt;
    fmt.copyfmt (o);
    Int_t dim= hReco->GetDimension(), ntxb= hReco->GetNbinsX()+2, ntyb= hReco->GetNbinsY()+2;
    if (hMeas->GetDimension() != dim || hMeas->GetNbinsX()+2 != ntxb || hMeas->GetNbinsY()+2 != ntyb) dim= 1;
    Int_t iwid= (dim==3) ? 8 : (dim==2) ? 7 : 5;
    const char* xwid= (dim==3) ? "===" : (dim==2) ? "==" : "";
    o << "===============================================================================" << xwid << std::endl
      << std::setw(iwid) << ""      << std::setw(9) << "Train" << std::setw(9) << "Train"    << std::setw(9) << "Test"  << std::setw(9) << "Test"  << std::setw(9) << "Unfolded";
    if (withError)
      o << std::setw(10)<<"Error on"<<std::setw(9) << "Diff" << std::setw(9) << "Pull" << std::endl;
    else
      o << std::setw(9) << "Diff" << std::endl;
    o << std::setw(iwid) << "Bin"   << std::setw(9) << "Truth" << std::setw(9) << "Measured" << std::setw(9) << "Truth" << std::setw(9) << "Input" << std::setw(9) << "Output";
    if (withError)
      o << std::setw(10)<<"Unfolding";
    o << std::endl;
    o << "===============================================================================" << xwid << std::endl;
    Double_t true_train_tot=0;
    Double_t meas_train_tot=0;
    Double_t true_test_tot=0;
    Double_t meas_test_tot=0;
    Double_t unf_tot=0;
    Double_t chi2= 0.0;
    Int_t ndf= 0, first= (_overflow ? 0 : 1);
    Int_t maxbin= _nt < _nm ? _nm : _nt;
    for (Int_t i = 0 ; i < maxbin; i++) {
      Int_t it= RooUnfolding::bin<TH1>(hReco, i, _overflow);
      Int_t im= RooUnfolding::bin<TH1>(hMeas, i, _overflow);

      if (dim==2 || dim==3) {
        Int_t iw= (dim==2) ? 3 : 2;
        Int_t ix= it%ntxb;
        Int_t iy= ((it-ix)/ntxb)%ntyb;
        o << std::setw(iw) << ix << ',' << std::setw(iw) << iy;
        if (dim==3) o << ',' << std::setw(iw) << ((it-ix)/ntxb - iy)/ntyb;
      } else
        o << std::setw(iwid) << i+first;
      o << std::fixed << std::setprecision(0);
      true_train_tot+=hTrainTrue->GetBinContent(it);
      meas_train_tot+=hTrain->GetBinContent(im);
      if (hTrue) true_test_tot+=hTrue->GetBinContent(it);
      meas_test_tot+=hMeas->GetBinContent(im);
      unf_tot+=hReco->GetBinContent(it);
      if (i<_nt){
        o << ' ' << std::setw(8) << hTrainTrue->GetBinContent(it);
      }
      else
        o << std::setw(9) << ' ';
      if (i<_nm)
        o << ' ' << std::setw(8) << hTrain->GetBinContent(im);
      else
        o << std::setw(9) << ' ';
      if (hTrue && i<_nt)
        o << ' ' << std::setw(8) << hTrue->GetBinContent(it);
      else
        o << std::setw(9) << ' ';
      if (i<_nm)
        o << ' ' << std::setw(8) << hMeas->GetBinContent(im);
      else
        o << std::setw(9) << ' ';
      o << std::setprecision(1);
      if (i<_nt) {
        Double_t y= hReco->GetBinContent(it), yerr = hReco->GetBinError(it);
        o << ' ' << std::setw(8) << y;
        if (withError)
          o << ' ' << std::setw(9) << yerr;
        if (hTrue &&
            (                       y!=0.0 || (withError &&                   yerr>0.0)) &&
            (hTrue->GetBinContent(it)!=0.0 || (withError && hTrue->GetBinError(it)>0.0))) {
          Double_t ydiff= y - hTrue->GetBinContent(it);
          o << ' ' << std::setw(8) << ydiff;
          if (withError && yerr>0.0) {
            ndf++;
            Double_t ypull = ydiff/yerr;
            chi2 += ypull*ypull;
            o << ' ' << std::setw(8) << ypull;
          }
        }
      }
      o << std::endl;

    }

    o << "===============================================================================" << xwid << std::endl
      << std::setw(iwid) << "" << std::fixed << std::setprecision(0)
      << ' ' << std::setw(8) << true_train_tot
      << ' ' << std::setw(8) << meas_train_tot;
    if (hTrue)
      o << ' ' << std::setw(8) << true_test_tot;
    else
      o << std::setw(9) << ' ';
    Double_t toterr= 0.0;
    if (meas_test_tot>0.0 && meas_train_tot>0.0) toterr= sqrt(meas_test_tot)*true_train_tot/meas_train_tot;
    o << ' ' << std::setw(8) << meas_test_tot << std::setprecision(1)
      << ' ' << std::setw(8) << unf_tot;
    if (withError) 
      o << ' ' << std::setw(9) << toterr;
    o << ' ' << std::setw(8) << unf_tot-true_test_tot;
    if(withError && toterr>0.0)
      o << ' ' << std::setw(8) <<(unf_tot-true_test_tot)/toterr;
    o << std::endl
      << "===============================================================================" << xwid << std::endl;
    o.copyfmt (fmt);
    if (hTrue) {
      if (chi_squ!=-999.0) {
        o << "Chi^2/NDF=" << chi_squ << "/" << ndf << " (bin-by-bin Chi^2=" << chi2 << ")" << std::endl;
      } else {
        o << "Bin-by-bin Chi^2/NDF=" << chi2 << "/" << ndf << std::endl;
        chi_squ= chi2;
      }
      if (chi_squ<=0.0) std::cerr << "Warning: Invalid Chi^2 Value" << std::endl;
    }
  }
  void printTable (std::ostream& o, const TVectorD& vTrainTrue, const TVectorD& vTrain,
                   const TVectorD& vMeas, const TVectorD& vReco,
                   Int_t nm, Int_t nt)
  {
    if (nm<=0) nm= vTrain    .GetNrows();
    if (nt<=0) nt= vTrainTrue.GetNrows();
    TH1* hTrainTrue =  RooUnfolding::createHist<TH1>(vTrainTrue, "","",Variable<TH1>(nt,0.0,nt,"xt"),false);
    TH1* hTrain     =  RooUnfolding::createHist<TH1>(vTrain,     "","",Variable<TH1>(nm,0.0,nm,"xm"),false);
    TH1* hMeas      =  RooUnfolding::createHist<TH1>(vMeas,      "","",Variable<TH1>(nm,0.0,nm,"xm"),false);
    TH1* hReco      =  RooUnfolding::createHist<TH1>(vReco,      "","",Variable<TH1>(nt,0.0,nt,"xt"),false);
    printTable<TH1> (o, hTrainTrue, hTrain, 0, hMeas, hReco, nm, nt);
    delete hTrainTrue;
    delete hTrain    ;
    delete hMeas     ;
    delete hReco     ;
  }


  template<> TH1* histNoOverflow<TH1> (const TH1* h, Bool_t overflow){
    if (!overflow) {   // also for 2D+
      TH1* hx= h2h1d<TH1,TH2> (h, h->GetNbinsX()*h->GetNbinsY()*h->GetNbinsZ());
      if (!hx) return hx;
      // clear under/overflow bins for cloned Hist
      hx->SetBinContent (0,                 0.0);
      hx->SetBinContent (hx->GetNbinsX()+1, 0.0);
      return hx;
    }
    Int_t nb= h->GetNbinsX(), s= h->GetSumw2N();
    Double_t xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax(), xb= (xhi-xlo)/nb;
    nb += 2;
    TH1* hx= new TH1F (h->GetName(), h->GetTitle(), nb, xlo-xb, xhi+xb);
    for (Int_t i= 0; i < nb; i++) {
      hx->SetBinContent (i+1, h->GetBinContent (i));
      if (s) hx->SetBinError   (i+1, h->GetBinError   (i));
    }
    return hx;
  }

  template<class Hist> Hist* resize (Hist* h, Int_t nx, Int_t ny, Int_t nz)
  {
    // Resize a histogram with a different number of bins.
    // Contents and errors are copied to the same bin numbers (the overflow bin
    // is copied to the new overflow bin) in the new histogram.
    // If the new histogram is larger than the old, the extra bins are zeroed.
    Int_t mx= h->GetNbinsX(), my= h->GetNbinsY(), mz= h->GetNbinsZ();
    Int_t nd= h->GetDimension();
    if (nx<0 || nd<1) nx= mx;
    if (ny<0 || nd<2) ny= my;
    if (nz<0 || nd<3) nz= mz;
    TH1* hc= (TH1*) h->Clone("resize_tmp");

    bool mod= false;
    if (nx!=mx) {
      Double_t xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax();
      h->GetXaxis()->Set (nx, xlo, xlo+((xhi-xlo)/mx)*nx);
      mod= true;
    }
    if (ny!=my) {
      Double_t ylo= h->GetYaxis()->GetXmin(), yhi= h->GetYaxis()->GetXmax();
      h->GetYaxis()->Set (ny, ylo, ylo+((yhi-ylo)/my)*ny);
      mod= true;
    }
    if (nz!=mz) {
      Double_t zlo= h->GetZaxis()->GetXmin(), zhi= h->GetZaxis()->GetXmax();
      h->GetZaxis()->Set (nz, zlo, zlo+((zhi-zlo)/mz)*nz);
      mod= true;
    }

    if (mod) {
      h->SetBinsLength();  // Just copies array, which isn't right for overflows or 2D/3D
      Int_t s= h->GetSumw2N();
      Int_t ox= mx+1, oy= my+1, oz= mz+1;  // old overflow bin
      Int_t px= nx+1, py= ny+1, pz= nz+1;  // new overflow bin

      if        (nd==1) {

        for (Int_t i= 0; i<=nx; i++) {
          h->SetBinContent (i, i>mx ? 0.0 : hc->GetBinContent (i));
          if (s) h->SetBinError   (i, i>mx ? 0.0 : hc->GetBinError   (i));
        }
        h->SetBinContent (px, h->GetBinContent (ox));
        if (s) h->SetBinError   (px, h->GetBinError   (ox));

      } else if (nd==2) {

        for (Int_t i= 0; i<=nx; i++) {
          for (Int_t j= 0; j<=ny; j++) {
            h->SetBinContent (i, j, i>mx||j>my ? 0.0 : hc->GetBinContent (i, j));
            if (s) h->SetBinError   (i, j, i>mx||j>my ? 0.0 : hc->GetBinError   (i, j));
          }
          h->SetBinContent (i, py, i>mx ? 0.0 : hc->GetBinContent (i, oy));
          if (s) h->SetBinError   (i, py, i>mx ? 0.0 : hc->GetBinError   (i, oy));
        }
        for (Int_t j= 0; j<=ny; j++) {
          h->SetBinContent (px, j, j>my ? 0.0 : hc->GetBinContent (ox, j));
          if (s) h->SetBinError   (px, j, j>my ? 0.0 : hc->GetBinError   (ox, j));
        }
        h->SetBinContent (px, py, hc->GetBinContent (ox, oy));
        if (s) h->SetBinError   (px, py, hc->GetBinError   (ox, oy));

      } else if (nd==3) {

        for (Int_t i= 0; i<=nx; i++) {
          for (Int_t j= 0; j<=ny; j++) {
            for (Int_t k= 0; k<=nz; k++) {
              h->SetBinContent (i, j, k, i>mx||j>my||k>mz ? 0.0 : hc->GetBinContent (i, j, k));
              if (s) h->SetBinError   (i, j, k, i>mx||j>my||k>mz ? 0.0 : hc->GetBinError   (i, j, k));
            }
            h->SetBinContent (i, j, pz, i>mx||j>my ? 0.0 : hc->GetBinContent (i, j, oz));
            if (s) h->SetBinError   (i, j, pz, i>mx||j>my ? 0.0 : hc->GetBinError   (i, j, oz));
          }
          h->SetBinContent (i, py, pz, i>mx ? 0.0 : hc->GetBinContent (i, oy, oz));
          if (s) h->SetBinError   (i, py, pz, i>mx ? 0.0 : hc->GetBinError   (i, oy, oz));
        }
        for (Int_t j= 0; j<=ny; j++) {
          for (Int_t k= 0; k<=nz; k++) {
            h->SetBinContent (px, j, k, j>my||k>mz ? 0.0 : hc->GetBinContent (ox, j, k));
            if (s) h->SetBinError   (px, j, k, j>my||k>mz ? 0.0 : hc->GetBinError   (ox, j, k));
          }
          h->SetBinContent (px, j, pz, j>my ? 0.0 : hc->GetBinContent (ox, j, oz));
          if (s) h->SetBinError   (px, j, pz, j>my ? 0.0 : hc->GetBinError   (ox, j, oz));
        }
        for (Int_t k= 0; k<=nz; k++) {
          for (Int_t i= 0; i<=nx; i++) {
            h->SetBinContent (i, py, k, i>mx||k>mz ? 0.0 : hc->GetBinContent (i, oy, k));
            if (s) h->SetBinError   (i, py, k, i>mx||k>mz ? 0.0 : hc->GetBinError   (i, oy, k));
          }
          h->SetBinContent (px, py, k, k>mz ? 0.0 : hc->GetBinContent (ox, oy, k));
          if (s) h->SetBinError   (px, py, k, k>mz ? 0.0 : hc->GetBinError   (ox, oy, k));
        }
        h->SetBinContent (px, py, pz, hc->GetBinContent (ox, oy, oz));
        if (s) h->SetBinError   (px, py, pz, hc->GetBinError   (ox, oy, oz));

      }
    }
    delete hc;
    return h;
  }
  template<> void printHistogram<TH2>(const TH2* hist){
    std::cout << hist->GetName() << " ( " << hist->GetNbinsX() << " x  " << hist->GetNbinsY() << " bins)" << std::endl;
    for(int i=0; i<hist->GetNbinsX()+2; ++i){
      for(int j=0; j<hist->GetNbinsY()+2; ++j){
        std::cout << i << " " << j << " " << hist->GetBinContent(i,j) << "+/-" << hist->GetBinError(i,j) << std::endl;
      }
    }
  }
  template<> void printHistogram<TH3>(const TH3* hist){
    std::cout << hist->GetName() << " ( " << hist->GetNbinsX() << " x  " << hist->GetNbinsY() << " x  " << hist->GetNbinsZ() <<" bins)" << std::endl;
    for(int i=0; i<hist->GetNbinsX()+2; ++i){
      for(int j=0; j<hist->GetNbinsY()+2; ++j){
        for(int k=0; k<hist->GetNbinsZ()+2; ++k){
          std::cout << i << " " << j << " " << k << " " << hist->GetBinContent(i,j,k) << "+/-" << hist->GetBinError(i,j,k) << std::endl;
        }
      }
    }
  }
  template<> void printHistogram<TH1>(const TH1* hist){
    if(hist->InheritsFrom(TH3::Class())){
      printHistogram<TH3>((TH3*)hist);
    } else if(hist->InheritsFrom(TH2::Class())){
      printHistogram<TH2>((TH2*)hist);
    } else {
      std::cout << hist->GetName() << " ( " << hist->GetNbinsX() << " bins)" << std::endl;
      for(int i=0; i<hist->GetNbinsX()+2; ++i){
        std::cout << i << " " << hist->GetBinContent(i) << "+/-" << hist->GetBinError(i) << std::endl;
      }
    }
  }

  template<> void subtract<TH1>(TH1* hist, const TVectorD& vec, double fac){
    for (Int_t i= 1; i<=hist->GetNbinsX()+1; i++){
      hist->SetBinContent (i, hist->GetBinContent(i)-(fac*vec[i-1]));
    }
  }

  template<> int findBin<TH1>(const TH1* h, Double_t x){
    // Get vector index (0..nx*ny-1) for bin containing (x,y) coordinates
    Int_t nx=   nBins(h,RooUnfolding::X);
    Int_t binx= findBin<TH1>(h,x,RooUnfolding::X) - 1;
    if (binx <  0)  return -1;
    if (binx >= nx) return nx;
    return binx;
  }

  template<> int findBin<TH1>(const TH1* h, Double_t x, Double_t y){
    // Get vector index (0..nx*ny-1) for bin containing (x,y) coordinates
    Int_t nx=   nBins(h,RooUnfolding::X);
    Int_t ny=   nBins(h,RooUnfolding::Y);
    Int_t binx= findBin<TH1>(h,x,RooUnfolding::X) - 1;
    if (binx <  0)  return -1;
    if (binx >= nx) return nx*ny;
    Int_t biny= findBin<TH1>(h,y,RooUnfolding::Y) - 1;
    if (biny <  0)  return -1;
    if (biny >= ny) return nx*ny;
    return binx + nx*biny;
  }
  template<> int findBin<TH2>(const TH2* h, Double_t x, Double_t y){
    return findBin<TH1>(h,x,y);
  }

  template<> int findBin<TH1>(const TH1* h, Double_t x, Double_t y, Double_t z){
    // Get vector index (0..nx*ny-1) for bin containing (x,y) coordinates
    Int_t nx=   nBins(h,RooUnfolding::X);
    Int_t ny=   nBins(h,RooUnfolding::Y);
    Int_t nz=   nBins(h,RooUnfolding::Z);
    Int_t binx= findBin<TH1>(h,x,RooUnfolding::X) - 1;
    if (binx <  0)  return -1;
    if (binx >= nx) return nx*ny*nz;
    Int_t biny= findBin<TH1>(h,y,RooUnfolding::Y) - 1;
    if (biny <  0)  return -1;
    if (biny >= ny) return nx*ny*nz;
    Int_t binz= findBin<TH1>(h,z,RooUnfolding::Z) - 1;
    if (binz <  0)  return -1;
    if (binz >= nz) return nx*ny*nz;
    return binx + nx*(biny + ny*binz);
  }
  template<> int findBin<TH3>(const TH3* h, Double_t x, Double_t y, Double_t z){
    return findBin<TH1>(h,x,y,z);
  }

}
  



template void RooUnfolding::add<TH2>(TH2*, TH2 const*);
template int RooUnfolding::sumW2N<TH2>(TH2 const*);
template int RooUnfolding::dim<TH2>(TH2 const*);
template double RooUnfolding::binCenter<TH2>(TH2 const*, int, RooUnfolding::Dimension);
template TH2* RooUnfolding::resize<TH2>(TH2*, int, int, int);
template TH1* RooUnfolding::createHist<TH1>(TVectorT<double> const&, char const*, char const*, RooUnfolding::Variable<TH1> const&, bool);
template double RooUnfolding::binWidth<TH1>(TH1 const*, int, RooUnfolding::Dimension);
template RooUnfolding::Variable<TH2> RooUnfolding::var<TH2>(TH2 const*, RooUnfolding::Dimension);
template int RooUnfolding::entries<TH1>(TH1 const*);
template void RooUnfolding::add<TH1>(TH1*, TH1 const*);
template TH1* RooUnfolding::resize<TH1>(TH1*, int, int, int);
template std::vector<RooUnfolding::Variable<TH1> > RooUnfolding::vars<TH1>(TH1 const*);
template int RooUnfolding::entries<TH2>(TH2 const*);
template RooUnfolding::Variable<TH1> RooUnfolding::var<TH1>(TH1 const*, RooUnfolding::Dimension);
template TH1* RooUnfolding::copy<TH1>(TH1 const*, bool, char const*, char const*);
template double RooUnfolding::binCenter<TH1>(TH1 const*, int, RooUnfolding::Dimension);
template TH2* RooUnfolding::copy<TH2>(TH2 const*, bool, char const*, char const*);

  

