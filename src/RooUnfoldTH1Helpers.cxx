#include "RooUnfoldTH1Helpers.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

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
  void reset(TH1* h){
    h->Reset();
  }
  int findBin(const TH1* h, double x, RooUnfolding::Dimension d){
    return getAxis(h,d)->FindBin(x);
  }
  double min(const TH1* hist, RooUnfolding::Dimension d){
    return getAxis(hist,d)->GetXmin();
  }
  double max(const TH1* hist, RooUnfolding::Dimension d){
    return getAxis(hist,d)->GetXmax();
  }  
  int sumW2N(const TH1* hist){
    return hist->GetSumw2N();
  }
  int entries(const TH1* hist){
    return hist->GetEntries();
  }
  int dim(const TH1* hist){
    return hist->GetDimension();
  }
 
  int nBins(const TH1* hist){
    return hist->GetNbinsX() * hist->GetNbinsY() * hist->GetNbinsZ();
  }
  int nBins(const TH1* hist, RooUnfolding::Dimension d){
    const TAxis* ax = getAxis(hist,d);
    return ax->GetNbins();
  }
  int bin(const TH1* h, Int_t i, Bool_t overflow)
  {
    // vector index (0..nx*ny-1) -> multi-dimensional histogram
    // global bin number (0..(nx+2)*(ny+2)-1) skipping under/overflow bins
    return (dim(h)<2) ? i+(overflow ? 0 : 1) : binDim(h,i);
  }
  int bin(const TH1* h, int i, int j, Bool_t overflow)
  {
    return h->GetBin(i,j);
  }
  double binCenter(const TH1*h, int i, RooUnfolding::Dimension d){
    const TAxis* ax = getAxis(h,d);
    return ax->GetBinCenter(i);
  }
  double binWidth(const TH1*h, int i, RooUnfolding::Dimension d){
    const TAxis* ax = getAxis(h,d);
    return ax->GetBinWidth(i);
  }
  double binHighEdge(const TH1*h, int i, RooUnfolding::Dimension d){
    const TAxis* ax = getAxis(h,d);
    return ax->GetBinUpEdge(i);
  }
  double binLowEdge(const TH1*h, int i, RooUnfolding::Dimension d){
    const TAxis* ax = getAxis(h,d);
    return ax->GetBinLowEdge(i);
  }
  
  void add(TH1* hista, const TH1* histb){
    hista->Add(histb);
  }
  void projectY(TH2* _res, TH1* _tru, bool overflow){
    Int_t s= _res->GetSumw2N();
    for (Int_t j= 1-overflow; j<_res->GetNbinsY()+1+overflow; j++) {
      Double_t ntru= 0.0, wtru= 0.0;
      for (Int_t i= 0; i<_res->GetNbinsX()+2; i++) {
        ntru +=      _res->GetBinContent (i, j);
        if (s) wtru += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t b= bin(_tru, j, overflow);
      _tru->SetBinContent (b,      ntru);
      if (s) _tru->SetBinError   (b, sqrt(wtru));
    }
  }
  void projectX(TH2* _res, TH1* _mes, bool overflow){
    Int_t s= _res->GetSumw2N();
    for (Int_t i= 1-overflow; i<_res->GetNbinsX()+1+overflow; i++) {
      Double_t nmes= 0.0, wmes= 0.0;
      for (Int_t j= 0; j<_res->GetNbinsY(); j++) {
        nmes +=      _res->GetBinContent (i, j);
        if (s) wmes += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t bin= RooUnfolding::bin (_mes, i, overflow);
      _mes->SetBinContent (bin,      nmes );
      if (s) _mes->SetBinError   (bin, sqrt(wmes));
    }
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
      Int_t bin= RooUnfolding::bin (_mes, i, overflow);
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
  TH1* copy(const TH1* orighist, bool reset, const char* name, const char* title){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH1* hist = (TH1*)(orighist ->Clone());
    if(name) hist->SetName(name);
    if(title) hist->SetTitle(title);
    if(reset) hist->Reset();
    return hist;
    TH1::AddDirectory (oldstat);
  }
  TH2* copy(const TH2* orighist, bool reset, const char* name, const char* title){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH2* hist = (TH2*)(orighist ->Clone());
    if(name) hist->SetName(name);
    if(title) hist->SetTitle(title);
    if(reset) hist->Reset();
    return hist;
    TH1::AddDirectory (oldstat);
  }  
  void binXYZ(const TH1* tru, int i, int& jx, int& jy, int& jz){
    Int_t j= RooUnfolding::bin(tru, i, false);
    if (dim(tru)>1) tru->GetBinXYZ (j, jx, jy, jz);
    else {
      jx = j;
      jy = 0;
      jz = 0;
    }
  }
  double binError(const TH1* h, Int_t i, Bool_t overflow)
  {
    // Bin error   by vector index
    return h->GetBinError   (bin (h, i, overflow));
  }
  double binContent (const TH1* h, Int_t i, Bool_t overflow)
  {
    // Bin content by vector index
    return h->GetBinContent (bin (h, i, overflow));
  }
  void setBinContent (TH1* h, int i, double val, Bool_t overflow)
  {
    // Bin content by vector index
    h->SetBinContent (bin (h, i, overflow), val);
  }
  void setBinContent (TH1* h, int i, int j, double val, Bool_t overflow)
  {
    // Bin content by vector index
    h->SetBinContent (bin (h, i, j, overflow), val);
  }

  template<>  
  TH2* createHist<TH2>(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, int nbinsy, double ymin, double ymax, const char* yname){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH2* hist = new TH2D (name,title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
    TH1::AddDirectory (oldstat);
    return hist;
  }
  template<> TH1* createHist<TH1>(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH1* hist = new TH1D (name,title, nbinsx, xmin, xmax);
    TH1::AddDirectory (oldstat);
    return hist;
  }
  TH1* h2h1d(TH1* h, int nb){
    return dynamic_cast<TH1*>(h->Clone());
  }
  
  TH1* h2h1d(TH2* h, int nb){
    TH1* h1d= new TH1F(h->GetName(), h->GetTitle(), nb, 0.0, 1.0);
    Int_t s= h->GetSumw2N();
    for (Int_t i= 0; i < nb; i++) {
      Int_t j= bin (h, i, false);  // don't bother with under/overflow bins (not supported for >1D)
      h1d->SetBinContent (i+1, h->GetBinContent (j));
      if (s) h1d->SetBinError   (i+1, h->GetBinError   (j));
    }
    return h1d;
  }

  TH2* copyHistogram(const TH2* h, bool includeOverflow){
    Int_t nx= nBins(h,RooUnfolding::X), ny= nBins(h,RooUnfolding::Y), s= sumW2N(h);
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
  template<> TH1* createHist<TH1>(const TVectorD& v, const char* name, const char* title, int nb, double min, double max, const char* xname, bool overflow){  
    // Sets the bin content of the histogram as that element of the input vector
    TH1* h = createHist<TH1>(name,title,nb,min,max,xname);
    if (overflow) nb += 2;
    for (Int_t i= 0; i < nb; i++) {
      Int_t j= RooUnfolding::bin(h, i, overflow);
      h->SetBinContent (j, v(i));
    }
    return h;
  }
  const char* varname(const TH1* h, Dimension d){  
    return "";
  }
}




