#include "RooUnfoldHelpers.h"
#include "RooUnfoldFitHelpers.h"
#include <RooUnfoldHelpers.tpp>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"

namespace RooUnfolding {
  Variable<RooAbsReal>::Variable(int nBins,double min,double max,const char* name) : _var(new RooRealVar(name,name,nBins,min,max)) {}
  Variable<RooAbsReal>::Variable(RooRealVar* var) : _var(var) {};

  template<> Variable<RooAbsReal> var(const RooAbsReal* h, Dimension d){
    return Variable<RooAbsReal>(NULL);
  }
  
  template<> int findBin<RooAbsReal>(const RooAbsReal* h, double x, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> int findBin<RooAbsReal>(const RooAbsReal* h, Double_t x){
    // TODO
    return 0;
  }
  template<> int findBin<RooAbsReal>(const RooAbsReal* h, Double_t x, Double_t y){
    // TODO
    return 0;
  }
  template<> int findBin<RooAbsReal>(const RooAbsReal* h, Double_t x, Double_t y, Double_t z){
    // TODO
    return 0;
  }

  template<> double min<RooAbsReal>(const RooAbsReal* hist, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> double max<RooAbsReal>(const RooAbsReal* hist, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> int sumW2N<RooAbsReal>(const RooAbsReal* hist){
    return 0;
  }
  template<> void projectY<RooAbsReal>(RooAbsReal* _res, RooAbsReal* _tru, bool overflow){
    // TODO
  } 
  template<> void projectX<RooAbsReal>(RooAbsReal* _res, RooAbsReal* _mes, bool overflow){
    // TODO
  }  
  template<> void subtractProjectX<RooAbsReal>(RooAbsReal* _res, RooAbsReal* _mes, RooAbsReal* _fak, bool overflow){
    // TODO
  }
  template<> int fill<RooAbsReal>(RooAbsReal* hist, double x, double w){
    return 0;
  }
  template<> int fill<RooAbsReal>(RooAbsReal* hist, double x, double y, double w){
    return 0;
  }  
  template<> int fill<RooAbsReal>(RooAbsReal* hist, double x, double y, double z, double w){
    return 0;
  }  
  template<> RooAbsReal* maybeCopy<RooAbsReal>(const RooAbsReal* r){
    return const_cast<RooAbsReal*>(r);
  }
  template<> bool maybeDelete<RooAbsReal>(RooAbsReal* r){
    return false;
  }  
  template<> int entries<RooAbsReal>(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  template<> int dim<RooAbsReal>(const RooAbsReal* hist){
    // TODO
    return 0;
  }
  template<> int nBins<RooAbsReal>(const RooAbsReal* hist, bool overflow){
    // TODO
    return 0;
  }
  template<> int nBins<RooAbsReal>(const RooAbsReal* hist, RooUnfolding::Dimension d, bool overflow){
    // TODO
    return 0;
  }
  template<> double binCenter<RooAbsReal>(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> double binWidth<RooAbsReal>(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> double binHighEdge<RooAbsReal>(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> double binLowEdge<RooAbsReal>(const RooAbsReal*h, int i, RooUnfolding::Dimension d){
    // TODO
    return 0;
  }
  template<> void binXYZ<RooAbsReal>(const RooAbsReal* tru, int i, int& jx, int& jy, int& jz){
    // TODO
  }
  template<> double binError<RooAbsReal>(const RooAbsReal* h, Int_t i, Bool_t overflow)
  {
    // Bin error   by vector index
    // TODO
    return 0;
  }  
  template<> double binContent<RooAbsReal> (const RooAbsReal* h, Int_t i, Bool_t overflow){
    // TODO
    return 0;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const char* name, const char* title, const Variable<RooAbsReal>& x, const Variable<RooAbsReal>& y){
    RooArgSet vars(*x._var,*y._var);
    RooDataHist* hist = new RooDataHist (name,title,vars);
    return new RooHistFunc(name,title,vars,vars,*hist);
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const char* name, const char* title, const std::vector<Variable<RooAbsReal>>& x){
    return NULL;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const TMatrixD& m, const char* name, const char* title, const Variable<RooAbsReal>& x, const Variable<RooAbsReal>& y){  
    // Sets the bin content of the histogram as that element of the input vector
    // TODO
    return NULL;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const TMatrixD& m, const TMatrixD& me, const char* name, const char* title, const Variable<RooAbsReal>& x, const Variable<RooAbsReal>& y){  
    // Sets the bin content of the histogram as that element of the input vector
    // TODO
    return NULL;
  }

  template<> RooAbsReal* createHist<RooAbsReal>(const TVectorD& v, const char* name, const char* title, const std::vector<Variable<RooAbsReal>>& x, bool overflow){
    // TODO
    return 0;
  }
  template<> RooAbsReal* createHist<RooAbsReal>(const TVectorD& v, const TVectorD& ve, const char* name, const char* title, const std::vector<Variable<RooAbsReal>>& x, bool overflow){
    // TODO
    return 0;
  }
  template<> const char* varname<RooAbsReal>(const RooAbsReal* h, Dimension d){
    // TODO
    return "";
  }

  template<> TVectorD subtract<RooAbsReal,TVectorD>(const TVectorD& orig, const RooAbsReal* hist, bool overflow) {
    // TODO
    TVectorD res;
    return res;
  }
  template<> RooAbsReal* resize<RooAbsReal> (RooAbsReal* h, Int_t nx, Int_t ny, Int_t nz){
    // TOOD
    return h;
  }

  template<> void printHistogram<RooAbsReal>(const RooAbsReal* h){
    // TODO
  }

  template<> void subtract<RooAbsReal>(RooAbsReal* hist, const TVectorD& vec, double fac){
    // TODO
  }

  template<> void h2v<RooAbsReal>  (const RooAbsReal* h, TVectorD& v, bool overflow){}
  template<> void h2ve<RooAbsReal>  (const RooAbsReal* h, TVectorD& v, bool overflow){}    

  template<> TVectorD h2v<RooAbsReal>  (const RooAbsReal* h, bool overflow){
    TVectorD v;
    h2v(h,v);
    return v;    
  }
  template<> TVectorD h2ve<RooAbsReal>  (const RooAbsReal* h, bool overflow){
    TVectorD v;
    h2ve(h,v);
    return v;    
  }
  template<> void h2mNorm<RooAbsReal>  (const RooAbsReal* h, TMatrixD& m, const RooAbsReal* norm, bool overflow){
    // sets Matrix to values of bins in a 2D input histogram    
  }
  template<> void h2meNorm<RooAbsReal>  (const RooAbsReal* h, TMatrixD& m, const RooAbsReal* norm, bool overflow){
    // sets Matrix to errors of bins in a 2D input histogram    
  }
  template<> TMatrixD h2mNorm<RooAbsReal>  (const RooAbsReal* h, const RooAbsReal* norm, bool overflow){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m;
    h2mNorm(h,m,norm, overflow);
    return m;
  }
  template<> TMatrixD h2meNorm<RooAbsReal>  (const RooAbsReal* h, const RooAbsReal* norm, bool overflow){
    // Returns Matrix of errors of bins in a 2D input histogram
    TMatrixD m;
    h2meNorm(h,m,norm, overflow);
    return m;
  }
  template<> void h2m  (const RooAbsReal* h, TMatrixD& m, bool overflow) { h2mNorm (h,m,(const RooAbsReal*)NULL,overflow); }
  template<> void h2me  (const RooAbsReal* h, TMatrixD& m, bool overflow){ h2meNorm(h,m,(const RooAbsReal*)NULL,overflow); };  
  template<> TMatrixD h2m<RooAbsReal>  (const RooAbsReal* h, bool overflow){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m;
    h2m(h,m,overflow);
    return m;
  }
  template<> TMatrixD h2me<RooAbsReal>  (const RooAbsReal* h, bool overflow){
    // Returns Matrix of errors of bins in a 2D input histogram
    TMatrixD m;
    h2me(h,m,overflow);
    return m;
  }
}  

template RooAbsReal* RooUnfolding::createHist<RooAbsReal>(TVectorT<double> const&, char const*, char const*, RooUnfolding::Variable<RooAbsReal> const&, bool);
template std::vector<RooUnfolding::Variable<RooAbsReal> > RooUnfolding::vars<RooAbsReal>(RooAbsReal const*); 
template void RooUnfolding::printTable<RooAbsReal>(std::ostream&, RooAbsReal const*, RooAbsReal const*, RooAbsReal const*, RooAbsReal const*, RooAbsReal const*, int, int, bool, RooUnfolding::ErrorTreatment, double);
