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
  template<class Hist> bool empty(const Hist* hist){
    return hist->GetEntries() == 0;
  }
  template<class Hist> int sumW2N(const Hist* hist){
    return hist->GetSumw2N();
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

  template<class Hist> void h2v  (const Hist* h, TVectorD& v, bool overflow){
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
  template<class Hist> void h2ve  (const Hist* h, TVectorD& v, bool overflow){
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
  template<class Hist> TVectorD h2v  (const Hist* h, bool overflow){
    // Returns Vector of values of bins in an input histogram
    TVectorD v(nBins(h,overflow));
    h2v(h,v,overflow);
    return v;
  }
  template<class Hist> TVectorD h2ve  (const Hist* h, bool overflow){
    // Returns Vector of values of bins in an input histogram
    TVectorD v(nBins(h,overflow));
    h2ve(h,v,overflow);
    return v;
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
  template<> TH2* createHist<TH2>(const char* name, const char* title, const std::vector<Variable<TH2>>& x){
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH2* hist = NULL;
    if(x.size() == 2){
      hist = new TH2D (name,title, x[0]._nBins, x[0]._min, x[0]._max, x[1]._nBins, x[1]._min, x[1]._max);
    } else {
      throw std::runtime_error(TString::Format("invalid dimensionality for ROOT histogram: %d",(int)x.size()).Data());
    }
    TH1::AddDirectory (oldstat);
    return hist;
  }  
  template<class Hist> Hist* createHist(const TVectorD& v, const char* name, const char* title, const std::vector<Variable<Hist>>& x, bool overflow){  
    // Sets the bin content of the histogram as that element of the input vector
    int nb = v.GetNrows();
    Hist* h = createHist<Hist>(name,title,x);
    for (Int_t i= 0; i < nb; i++) {
      Int_t j= RooUnfolding::bin<TH1>(h, i, overflow);
      h->SetBinContent (j, v(i));
    }
    return h;
  }
  template<class Hist> Hist* createHist(const TVectorD& v, const TVectorD& ve, const char* name, const char* title, const std::vector<Variable<Hist>>& x, bool overflow){  
    // Sets the bin content of the histogram as that element of the input vector
    int nb = v.GetNrows();
    Hist* h = createHist<Hist>(name,title,x);
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
  


template int RooUnfolding::sumW2N<TH1>(TH1 const*);
template int RooUnfolding::sumW2N<TH2>(TH2 const*);
template int RooUnfolding::sumW2N<TH3>(TH3 const*);
template double RooUnfolding::binCenter<TH1>(TH1 const*, int, RooUnfolding::Dimension);
template double RooUnfolding::binCenter<TH2>(TH2 const*, int, RooUnfolding::Dimension);
template double RooUnfolding::binCenter<TH3>(TH3 const*, int, RooUnfolding::Dimension);
template double RooUnfolding::binWidth<TH1>(TH1 const*, int, RooUnfolding::Dimension);
template double RooUnfolding::binWidth<TH2>(TH2 const*, int, RooUnfolding::Dimension);
template double RooUnfolding::binWidth<TH3>(TH3 const*, int, RooUnfolding::Dimension);
template TH1* RooUnfolding::resize<TH1>(TH1*, int, int, int);
template TH2* RooUnfolding::resize<TH2>(TH2*, int, int, int);
template TH3* RooUnfolding::resize<TH3>(TH3*, int, int, int);
template int RooUnfolding::dim<TH1>(TH1 const*);
template int RooUnfolding::dim<TH2>(TH2 const*);
template int RooUnfolding::dim<TH3>(TH3 const*);
template std::vector<RooUnfolding::Variable<TH1> > RooUnfolding::vars<TH1>(TH1 const*);
template std::vector<RooUnfolding::Variable<TH2> > RooUnfolding::vars<TH2>(TH2 const*);
template std::vector<RooUnfolding::Variable<TH3> > RooUnfolding::vars<TH3>(TH3 const*);
template RooUnfolding::Variable<TH1> RooUnfolding::var<TH1>(TH1 const*, RooUnfolding::Dimension);
template RooUnfolding::Variable<TH2> RooUnfolding::var<TH2>(TH2 const*, RooUnfolding::Dimension);
template RooUnfolding::Variable<TH3> RooUnfolding::var<TH3>(TH3 const*, RooUnfolding::Dimension);
template TH1* RooUnfolding::createHist<TH1>(char const*, char const*, std::vector<RooUnfolding::Variable<TH1> > const&);
template TH2* RooUnfolding::createHist<TH2>(char const*, char const*, std::vector<RooUnfolding::Variable<TH2> > const&);
template TH1* RooUnfolding::createHist<TH1>(TVectorT<double> const&, char const*, char const*, RooUnfolding::Variable<TH1> const&, bool);
template TH2* RooUnfolding::createHist<TH2>(TVectorT<double> const&, char const*, char const*, RooUnfolding::Variable<TH2> const&, bool);
template TH1* RooUnfolding::createHist<TH1>(TVectorT<double> const&, TVectorT<double> const&, char const*, char const*, std::vector<RooUnfolding::Variable<TH1> > const&, bool);
template TH2* RooUnfolding::createHist<TH2>(TVectorT<double> const&, TVectorT<double> const&, char const*, char const*, std::vector<RooUnfolding::Variable<TH2> > const&, bool);
template void RooUnfolding::h2v<TH1>(TH1 const*, TVectorT<double>&, bool);
template TVectorT<double> RooUnfolding::h2v<TH1>(TH1 const*, bool);
template TVectorT<double> RooUnfolding::h2ve<TH1>(TH1 const*, bool);
template void RooUnfolding::h2ve<TH1>(TH1 const*, TVectorT<double>&, bool);
template void RooUnfolding::printTable<TH1>(std::ostream&, TH1 const*, TH1 const*, TH1 const*, TH1 const*, TH1 const*, bool, RooUnfolding::ErrorTreatment, double);
template bool RooUnfolding::empty<TH1>(const TH1*);
template bool RooUnfolding::empty<TH2>(const TH2*);
template bool RooUnfolding::empty<TH3>(const TH3*);
  

