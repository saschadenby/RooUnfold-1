#include "RooUnfoldTH1Helpers.h"

#include <iostream>
#include <ostream>
#include <sstream>
#include <iomanip>

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
 
  int nBins(const TH1* hist, bool overflow){
    int d = dim(hist);
    if(d==1){
      return (hist->GetNbinsX()+2*overflow);
    } else if(d==2){
      return ( (hist->GetNbinsX()+2*overflow) * (hist->GetNbinsY()+2*overflow) );
    } else if(d==3){
      return ( (hist->GetNbinsX()+2*overflow) * (hist->GetNbinsY()+2*overflow) * (hist->GetNbinsZ()+2*overflow) );
    }
  }
  int nBins(const TH1* hist, RooUnfolding::Dimension d, bool overflow){
    const TAxis* ax = getAxis(hist,d);
    return ax->GetNbins()+2*overflow;
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
      Int_t b= RooUnfolding::bin (_mes, i, overflow);
      _mes->SetBinContent (b,      nmes );
      if (s) _mes->SetBinError   (b, sqrt(wmes));
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
  double binContent (const TH1* h, int i, int j, Bool_t overflow)
  {
    // Bin content by vector index
    return h->GetBinContent (bin (h, i, j, overflow));
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
  TH1* h2h1d(const TH1* h, int nb){
    return dynamic_cast<TH1*>(h->Clone());
  }
  
  TH1* h2h1d(const TH2* h, int nb){
    TH1* h1d= new TH1F(h->GetName(), h->GetTitle(), nb, 0.0, 1.0);
    Int_t s= h->GetSumw2N();
    for (Int_t i= 0; i < nb; i++) {
      Int_t j= bin (h, i, false);  // don't bother with under/overflow bins (not supported for >1D)
      h1d->SetBinContent (i+1, h->GetBinContent (j));
      if (s) h1d->SetBinError   (i+1, h->GetBinError   (j));
    }
    return h1d;
  }
  void h2m  (const TH2* h, TMatrixD& m){
    // sets Matrix to values of bins in a 2D input histogram
    m.ResizeTo(h->GetNbinsX(),h->GetNbinsY());
    for (Int_t i= 0; i < h->GetNbinsX(); ++i) {
      for (Int_t j= 0; j < h->GetNbinsY(); ++j) {
        m(i,j)= h->GetBinContent(i+1,j+1);
      }
    }
  }
  TVectorD h2v  (const TH1* h, bool overflow){
    // Returns Vector of values of bins in an input histogram
    TVectorD v(nBins(h,overflow));
    h2v(h,v,overflow);
    return v;
  }
  
  void h2v  (const TH1* h, TVectorD& v, bool overflow){
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
  void h2ve  (const TH1* h, TVectorD& v, bool overflow){
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
  void h2me  (const TH2* h, TMatrixD& m){
  // sets Matrix to errors of bins in a 2D input histogram    
    m.ResizeTo(h->GetNbinsX(),h->GetNbinsY());
    for (Int_t i= 0; i < h->GetNbinsX(); ++i) {
      for (Int_t j= 0; j < h->GetNbinsY(); ++j) {
        m(i,j)= h->GetBinError(i+1,j+1);
      }
    }
  }  
  TMatrixD h2m  (const TH2* h){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m(h->GetNbinsX(),h->GetNbinsY());
    h2m(h,m);
    return m;
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

  void setContents(TH1* h,const std::vector<double>& values,const std::vector<double>& errors, bool overflow){
    for (Int_t i= 0; i < values.size(); i++) {
      Int_t b= RooUnfolding::bin (h, i, overflow);
      h->SetBinContent(b,values[i]);
      h->SetBinError(b,errors[i]);
    }
  }

  template<> TVectorD subtract<TVectorD>(const TVectorD& orig, const TH1* hist, bool overflow){
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

  void printTable (std::ostream& o, const TVectorD& vTrainTrue, const TVectorD& vTrain,
                   const TVectorD& vMeas, const TVectorD& vReco,
                   Int_t nm, Int_t nt)
  {
    if (nm<=0) nm= vTrain    .GetNrows();
    if (nt<=0) nt= vTrainTrue.GetNrows();
    TH1* hTrainTrue =  RooUnfolding::createHist<TH1>(vTrainTrue, "","", nt,0.0,nt,"xt",false);
    TH1* hTrain     =  RooUnfolding::createHist<TH1>(vTrain,     "","", nm,0.0,nm,"xm",false);
    TH1* hMeas      =  RooUnfolding::createHist<TH1>(vMeas,      "","", nm,0.0,nm,"xm",false);
    TH1* hReco      =  RooUnfolding::createHist<TH1>(vReco,      "","", nt,0.0,nt,"xt",false);
    printTable (o, hTrainTrue, hTrain, 0, hMeas, hReco, nm, nt);
    delete hTrainTrue;
    delete hTrain    ;
    delete hMeas     ;
    delete hReco     ;
  }

  void printTable (std::ostream& o, const TH1* hTrainTrue, const TH1* hTrain,
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
      Int_t it= RooUnfolding::bin(hReco, i, _overflow);
      Int_t im= RooUnfolding::bin(hMeas, i, _overflow);

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

  TH1* histNoOverflow (const TH1* h, Bool_t overflow){
    if (!overflow) {   // also for 2D+
      TH1* hx= h2h1d (h, h->GetNbinsX()*h->GetNbinsY()*h->GetNbinsZ());
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


  TH1* resize (TH1* h, Int_t nx, Int_t ny, Int_t nz)
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

  void printHistogram(const TH1* hist){
    std::cout << hist->GetName() << " ( " << hist->GetNbinsX() << " bins)" << std::endl;
    for(int i=0; i<hist->GetNbinsX()+2; ++i){
      std::cout << hist->GetBinContent(i) << "+/-" << hist->GetBinError(i) << std::endl;
    }

  }

  void subtract(TH1* hist, const TVectorD& vec, double fac){
    for (Int_t i= 1; i<=hist->GetNbinsX()+1; i++){
      hist->SetBinContent (i, hist->GetBinContent(i)-(fac*vec[i-1]));
    }
  }
}
  




