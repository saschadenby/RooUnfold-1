#include "RooUnfoldFitHelpers.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"

namespace RooUnfolding {
  void reset(RooAbsReal* r){
    // TODO
  }
  int findBin(const RooAbsReal* h, double x, RooUnfolding::Dimension d){
    // TOOD
    return 0;
  }
  double min(const RooAbsReal* hist, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  double max(const RooAbsReal* hist, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  int sumW2N(const RooAbsReal* hist){
    return 0;
  }
  void add(RooAbsReal* hista, RooAbsReal* histb){
    // TODO
  }  
  void projectY(RooAbsReal* _res, RooAbsReal* _tru, bool overflow){
    // TODO
  } 
  void projectX(RooAbsReal* _res, RooAbsReal* _mes, bool overflow){
    // TODO
  }  
  void subtractProjectX(RooAbsReal* _res, RooAbsReal* _mes, RooAbsReal* _fak, bool overflow){
    // TODO
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
  RooAbsReal* copy(const RooAbsReal* r, bool reset, const char* name, const char* title){
    RooAbsReal* retval = (RooAbsReal*)(r->clone());
    if(name) retval->SetName(name);
    if(title) retval->SetTitle(title);
    return retval;
  }
  int entries(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  int dim(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  int nBins(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  int nBins(const RooAbsReal* hist, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  double binCenter(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  double binWidth(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  double binHighEdge(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  double binLowEdge(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  void binXYZ(const RooAbsReal* tru, int i, int& jx, int& jy, int& jz){
    // TODO
  }
  double binError(const RooAbsReal* h, Int_t i, Bool_t overflow)
  {
    // Bin error   by vector index
    // TODO
    return 0;
  }  
  double binContent (const RooAbsReal* h, Int_t i, Bool_t overflow){
    // TODO
    return 0;
  }
  void setBinContent (RooAbsReal* h, int i, double val, Bool_t overflow){
    // TODO
  }
  void setBinContent (RooAbsReal* h, int i, int j, double val, Bool_t overflow){
    // TODO
  }
  template<class Hist2D> Hist2D* createHist(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, int nbinsy, double ymin, double ymax, const char* yname);
  template<> RooAbsReal* createHist<RooAbsReal>(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, int nbinsy, double ymin, double ymax, const char* yname){
    RooRealVar* x = new RooRealVar(xname,xname,nbinsx,xmin,xmax);
    RooRealVar* y = new RooRealVar(yname,yname,nbinsy,ymin,ymax);
    RooArgSet vars(*x,*y);
    RooDataHist* hist = new RooDataHist (name,title,vars);
    return new RooHistFunc(name,title,vars,vars,*hist);
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname){
    RooRealVar* x = new RooRealVar(xname,xname,nbinsx,xmin,xmax);
    RooArgSet vars(*x);    
    RooDataHist* hist = new RooDataHist (name,title,vars);
    return new RooHistFunc(name,title,vars,vars,*hist);
  }
  RooAbsReal* h2h1d(const RooAbsReal* h, int nb){
    // TODO
    return 0;
  }
  RooAbsReal* copyHistogram(const RooAbsReal* h, bool includeOverflow){
    // TODO
    return 0;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const TVectorD& v, const char* name, const char* title, int nb, double min, double max, const char* xname, bool overflow){
    // TOOD
    return 0;
  }
  const char* varname(const RooAbsReal* h, Dimension d){
    // TODO
    return "";
  }
  RooAbsReal* histNoOverflow(const RooAbsReal* h, bool overflow){
    // overflow not supported anyway, do nothing
    return copy(h);
  }

  void setContents(RooAbsReal* h,const std::vector<double>& values,const std::vector<double>& errors, bool overflow){
    // TODO
  }
  template<> TVectorD subtract<TVectorD>(const TVectorD& orig, const RooAbsReal* hist, bool overflow) {
    // TODO
    TVectorD res;
    return res;
  }
  void printTable (std::ostream& o, const RooAbsReal* hTrainTrue, const RooAbsReal* hTrain,
                   const RooAbsReal* hTrue, const RooAbsReal* hMeas, const RooAbsReal* hReco,
                   Int_t _nm, Int_t _nt, Bool_t _overflow,
                   RooUnfolding::ErrorTreatment withError, Double_t chi_squ){
    // TODO
  }
  RooAbsReal* resize (RooAbsReal* h, Int_t nx, Int_t ny, Int_t nz){
    // TOOD
    return h;
  }

  void printHistogram(const RooAbsReal* h){
    // TODO
  }

  void subtract(RooAbsReal* hist, const TVectorD& vec, double fac){
    // TODO
  }


  TVectorD h2v  (const RooAbsReal* h){
    TVectorD v;
    h2v(h,v);
    return v;    
  }
  void h2v  (const RooAbsReal* h, TVectorD& v){}
  void h2ve  (const RooAbsReal* h, TVectorD& v){}    
  void h2me  (const RooAbsReal* h, TMatrixD& m){}
  
  void h2m  (const RooAbsReal* h, TMatrixD& m){
    // sets Matrix to values of bins in a 2D input histogram    
  }
  TMatrixD h2m  (const RooAbsReal* h){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m;
    h2m(h,m);
    return m;
  }
  
  
}  

  


