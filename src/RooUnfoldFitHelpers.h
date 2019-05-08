#ifndef ROOUNFOLDFITHELPERS_HH
#define ROOUNFOLDFITHELPERS_HH

class RooAbsReal;

#include "RooUnfoldHelpers.h"

namespace RooUnfolding {
  void reset(RooAbsReal* r);
  int findBin(const RooAbsReal* h, double x, RooUnfolding::Dimension d);
  double min(const RooAbsReal* hist, RooUnfolding::Dimension d);
  double max(const RooAbsReal* hist, RooUnfolding::Dimension d);
  int sumW2N(const RooAbsReal* hist);
  void add(RooAbsReal* hista, RooAbsReal* histb);
  void projectY(RooAbsReal* _res, RooAbsReal* _tru, bool overflow);
  void projectX(RooAbsReal* _res, RooAbsReal* _mes, bool overflow);
  void subtractProjectX(RooAbsReal* _res, RooAbsReal* _mes, RooAbsReal* _fak, bool overflow);
  int fill(RooAbsReal* hist, double x, double w);
  int fill(RooAbsReal* hist, double x, double y, double w);
  int fill(RooAbsReal* hist, double x, double y, double z, double w);
  RooAbsReal* copy(const RooAbsReal* r, bool reset=false, const char* name = 0, const char* title = 0);
  int entries(const RooAbsReal* hist);
  int dim(const RooAbsReal* hist);
  int nBins(const RooAbsReal* hist, bool overflow=false);
  int nBins(const RooAbsReal* hist, RooUnfolding::Dimension d, bool overflow=false);
  double binCenter(const RooAbsReal*h, int i, RooUnfolding::Dimension d);
  double binWidth(const RooAbsReal*h, int i, RooUnfolding::Dimension d);
  double binHighEdge(const RooAbsReal*h, int i, RooUnfolding::Dimension d);
  double binLowEdge(const RooAbsReal*h, int i, RooUnfolding::Dimension d);
  void binXYZ(const RooAbsReal* tru, int i, int& jx, int& jy, int& jz);
  double binError(const RooAbsReal* h, int i, bool overflow);
  double binContent (const RooAbsReal* h, int i, bool overflow);
  double binContent (const RooAbsReal* h, int i, int j, Bool_t overflow);  
  void setBinContent (RooAbsReal* h, int i, double val, bool overflow);
  void setBinContent (RooAbsReal* h, int i, int j, double val, bool overflow);
  RooAbsReal* h2h1d(const RooAbsReal* h, int nb);
  TVectorD h2v  (const RooAbsReal* h, bool oveflow=false);    
  TVectorD h2ve  (const RooAbsReal* h, bool oveflow=false);    
  void h2v  (const RooAbsReal* h, TVectorD& v, bool overflow=false);
  void h2ve  (const RooAbsReal* h, TVectorD& v, bool overflow=false);    
  TMatrixD h2m  (const RooAbsReal* h, bool overflow=false);
  TMatrixD h2me  (const RooAbsReal* h, bool overflow=false);
  void h2m  (const RooAbsReal* h, TMatrixD& m, bool overflow=false);
  void h2me  (const RooAbsReal* h, TMatrixD& m, bool overflow=false);  
  TMatrixD h2mNorm  (const RooAbsReal* h, const RooAbsReal* norm, bool overflow=false);
  TMatrixD h2meNorm  (const RooAbsReal* h, const RooAbsReal* norm, bool overflow=false);
  void h2mNorm  (const RooAbsReal* h, TMatrixD& m, const RooAbsReal* norm, bool overflow=false);
  void h2meNorm  (const RooAbsReal* h, TMatrixD& m, const RooAbsReal* norm, bool overflow=false);  
  RooAbsReal* copyHistogram(const RooAbsReal* h, bool includeOverflow);
  const char* varname(const RooAbsReal* h, Dimension d);
  void setContents(RooAbsReal* h,const std::vector<double>& values,const std::vector<double>& errors, bool overflow);
  template<class V> V subtract(const TVectorD& orig, const RooAbsReal* hist, bool overflow);
  void printHistogram(const RooAbsReal* h);  
  void printTable (std::ostream& o, const RooAbsReal* hTrainTrue, const RooAbsReal* hTrain,
                   const RooAbsReal* hTrue, const RooAbsReal* hMeas, const RooAbsReal* hReco,
                   Int_t _nm=0, Int_t _nt=0, Bool_t _overflow=kFALSE,
                   RooUnfolding::ErrorTreatment withError=RooUnfolding::kDefault, Double_t chi_squ=-999.0);
  RooAbsReal* histNoOverflow(const RooAbsReal* h, bool overflow);
  RooAbsReal* resize (RooAbsReal* h, Int_t nx, Int_t ny = 0, Int_t nz = 0);
  void subtract(RooAbsReal* hist, const TVectorD& vec, double fac);
}  

#endif
