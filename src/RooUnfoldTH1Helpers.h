#ifndef ROOUNFOLDTH1HELPERS_HH
#define ROOUNFOLDTH1HELPERS_HH

class TH1;
class TH2;

#include "RooUnfoldHelpers.h"

namespace RooUnfolding {
  void reset(TH1* h);
  int findBin(const TH1* h, double x, RooUnfolding::Dimension d);
  double min(const TH1* hist, RooUnfolding::Dimension d);
  double max(const TH1* hist, RooUnfolding::Dimension d);
  int sumW2N(const TH1* hist);
  int entries(const TH1* hist);
  int dim(const TH1* hist);
  int nBins(const TH1* hist, bool overflow=false);
  int nBins(const TH1* hist, RooUnfolding::Dimension d, bool overflow=false);
  int bin(const TH1* h, int i, bool overflow);
  int bin(const TH1* h, int i, int j, bool overflow);
  double binCenter(const TH1*h, int i, RooUnfolding::Dimension d);
  double binWidth(const TH1*h, int i, RooUnfolding::Dimension d);
  double binHighEdge(const TH1*h, int i, RooUnfolding::Dimension d);
  double binLowEdge(const TH1*h, int i, RooUnfolding::Dimension d);
  void add(TH1* hista, const TH1* histb);
  void projectY(TH2* _res, TH1* _tru, bool overflow);
  void projectX(TH2* _res, TH1* _mes, bool overflow);
  void subtractProjectX(TH2* _res, TH1* _mes, TH1* _fak, bool overflow);
  int fill(TH1* hist, double x, double w);
  int fill(TH1* hist, double x, double y, double w);
  int fill(TH1* hist, double x, double y, double z, double w);
  int fill(TH2* hist, double x, double y, double w);
  TH1* copy(const TH1* orighist, bool reset = false, const char* name = 0, const char* title = 0);
  TH2* copy(const TH2* orighist, bool reset = false, const char* name = 0, const char* title = 0);
  void binXYZ(const TH1* tru, int i, int& jx, int& jy, int& jz);
  double binError(const TH1* h, int i, bool overflow);
  double binContent (const TH1* h, int i, bool overflow);
  double binContent (const TH1* h, int i, int j, Bool_t overflow);  
  void setBinContent (TH1* h, int i, double val, bool overflow);
  void setBinContent (TH1* h, int i, int j, double val, bool overflow);
  TH1* h2h1d(const TH1* h, int nb);  
  TH1* h2h1d(const TH2* h, int nb);
  TVectorD h2v  (const TH1* h,bool overflow = false);
  TVectorD h2ve  (const TH1* h,bool overflow = false);
  void h2v  (const TH1* h, TVectorD& v,bool overflow = false);
  void h2ve  (const TH1* h, TVectorD& v,bool overflow = false);
  void h2m  (const TH2* h, TMatrixD& m, bool overflow = false);
  void h2me  (const TH2* h, TMatrixD& m, bool overflow = false);  
  TMatrixD h2m  (const TH2* h, bool overflow = false);  
  TMatrixD h2me  (const TH2* h, bool overflow = false);  
  void h2mNorm  (const TH2* h, TMatrixD& m, const TH1* norm, bool overflow = false);
  void h2meNorm  (const TH2* h, TMatrixD& m, const TH1* norm, bool overflow = false);  
  TMatrixD h2mNorm  (const TH2* h, const TH1* norm, bool overflow = false);  
  TMatrixD h2meNorm  (const TH2* h, const TH1* norm, bool overflow = false);  
  TH2* copyHistogram(const TH2* h, bool includeOverflow);
  const char* varname(const TH1* h, Dimension d);
  void setContents(TH1* h,const std::vector<double>& values,const std::vector<double>& errors, bool overflow);
  template<class V> V subtract(const TVectorD& orig, const TH1* hist, bool overflow);
  void subtract(TH1* hist, const TVectorD& vec, double fac);

  void printHistogram(const TH1* h);
  void printTable (std::ostream& o, const TH1* hTrainTrue, const TH1* hTrain,
                   const TH1* hTrue, const TH1* hMeas, const TH1* hReco,
                   Int_t _nm=0, Int_t _nt=0, Bool_t _overflow=kFALSE,
                   RooUnfolding::ErrorTreatment withError=RooUnfolding::kDefault, Double_t chi_squ=-999.0);

  TH1* histNoOverflow (const TH1* h, Bool_t overflow);
  TH1* resize (TH1* h, Int_t nx, Int_t ny = 0, Int_t nz = 0);  
  
}
#endif
