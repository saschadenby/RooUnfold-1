#ifndef ROOUNFOLDHELPERS_HH
#define ROOUNFOLDHELPERS_HH

#include <TVectorD.h>
#include <TMatrixD.h>
#include <iostream>

namespace RooUnfolding {

  enum Algorithm {       // Selection of unfolding algorithm:
    kNone,               //   no unfolding (dummy routines in RooUnfold)
    kBayes,              //   RooUnfoldBayes
    kSVD,                //   RooUnfoldSvd
    kBinByBin,           //   RooUnfoldBinByBin
    kTUnfold,            //   RooUnfoldTUnfold
    kInvert,             //   RooUnfoldInvert
    kDagostini,          //   RooUnfoldDagostini
    kIDS                 //   RooUnfoldIds
  };

  enum ErrorTreatment {  // Error treatment:
    kNoError,            //   no error treatment: returns sqrt(N)
    kErrors,             //   bin-by-bin errors (diagonal covariance matrix)
    kCovariance,         //   covariance matrix from unfolding
    kCovToy,             //   covariance matrix from toy MC
    kDefault=-1          //   not specified
  };
  
  enum Dimension { X, Y, Z };
  template<class Hist> struct Variable;

  template<class Hist> double min(const Hist* hist, RooUnfolding::Dimension d);
  template<class Hist> double max(const Hist* hist, RooUnfolding::Dimension d);
  template<class Hist> int sumW2N(const Hist* hist);
  template<class Hist> int entries(const Hist* hist);
  template<class Hist> int dim(const Hist* hist);
  template<class Hist> int nBins(const Hist* hist, bool overflow=false);
  template<class Hist> int nBins(const Hist* hist, RooUnfolding::Dimension d, bool overflow=false);
  template<class Hist> int bin(const Hist* h, int i, bool overflow);
  template<class Hist> int bin(const Hist* h, int i, int j, bool overflow);
  template<class Hist> double binCenter(const Hist*h, int i, RooUnfolding::Dimension d);
  template<class Hist> double binWidth(const Hist*h, int i, RooUnfolding::Dimension d);
  template<class Hist> double binHighEdge(const Hist*h, int i, RooUnfolding::Dimension d);
  template<class Hist> double binLowEdge(const Hist*h, int i, RooUnfolding::Dimension d);
  template<class Hist> int fill(Hist* hist, double x, double w);
  template<class Hist> int fill(Hist* hist, double x, double y, double w);
  template<class Hist> int fill(Hist* hist, double x, double y, double z, double w);
  template<class Hist> void binXYZ(const Hist* tru, int i, int& jx, int& jy, int& jz);
  template<class Hist> double binError(const Hist* h, int i, bool overflow);
  template<class Hist> double binContent (const Hist* h, int i, bool overflow);
  template<class Hist> double binContent (const Hist* h, int i, int j, Bool_t overflow);  
  template<class Hist> TVectorD h2v  (const Hist* h,bool overflow = false);
  template<class Hist> TVectorD h2ve  (const Hist* h,bool overflow = false);
  template<class Hist> void h2v  (const Hist* h, TVectorD& v,bool overflow = false);
  template<class Hist> void h2ve  (const Hist* h, TVectorD& v,bool overflow = false);
  template<class Hist> const char* varname(const Hist* h, Dimension d);
  template<class Hist> void printHistogram(const Hist* h);
  template<class Hist> void printTable (std::ostream& o, const Hist* hTrainTrue, const Hist* hTrain, const Hist* hTrue, const Hist* hMeas, const Hist* hReco, Int_t _nm=0, Int_t _nt=0, Bool_t _overflow=kFALSE, RooUnfolding::ErrorTreatment withError=RooUnfolding::kDefault, Double_t chi_squ=-999.0);
  template<class Hist> Hist* resize (Hist* h, Int_t nx, Int_t ny = 0, Int_t nz = 0);  
  template<class Hist> void subtract(Hist* hist, const TVectorD& vec, double fac);
  template<class Hist, class V> V subtract(const TVectorD& orig, const Hist* hist, bool overflow);

  template<class Hist2D> int fill(Hist2D* hist, double x, double y, double w);

  template<class Hist> Hist* maybeCopy(const Hist* orighist);
  template<class Hist> bool maybeDelete(Hist* orighist);  
  
  template<class Hist> int findBin(const Hist* h, double x, RooUnfolding::Dimension d);
  template<class Hist1D> int findBin(const Hist1D* h, Double_t x);
  template<class Hist2D> int findBin(const Hist2D* h, Double_t x, Double_t y);
  template<class Hist3D> int findBin(const Hist3D* h, Double_t x, Double_t y, double_t z);

  template<class Hist, class Hist2D> Hist* h2h1d(const Hist* h, int nb);
  template<class Hist, class Hist2D> Hist* h2h1d(const Hist2D* h, int nb);
  template<class Hist, class Hist2D> void h2mNorm  (const Hist2D* h, TMatrixD& m, const Hist* norm, bool overflow = false);
  template<class Hist, class Hist2D> void h2meNorm  (const Hist2D* h, TMatrixD& m, const Hist* norm, bool overflow = false);  
  template<class Hist, class Hist2D> TMatrixD h2mNorm  (const Hist2D* h, const Hist* norm, bool overflow = false);  
  template<class Hist, class Hist2D> TMatrixD h2meNorm  (const Hist2D* h, const Hist* norm, bool overflow = false);  
  template<class Hist, class Hist2D> void projectY(Hist2D* _res, Hist* _tru, bool overflow);
  template<class Hist, class Hist2D> void projectX(Hist2D* _res, Hist* _mes, bool overflow);
  template<class Hist, class Hist2D> void subtractProjectX(Hist2D* _res, Hist* _mes, Hist* _fak, bool overflow);

  template<class Hist2D> void h2m  (const Hist2D* h, TMatrixD& m, bool overflow = false);
  template<class Hist2D> void h2me  (const Hist2D* h, TMatrixD& m, bool overflow = false);
  template<class Hist2D> TMatrixD h2m  (const Hist2D* h, bool overflow = false);  
  template<class Hist2D> TMatrixD h2me  (const Hist2D* h, bool overflow = false);  

  template<class Hist2D> Hist2D* createHist(const char* name, const char* title, const Variable<Hist2D>& x, const Variable<Hist2D>& y);
  template<class Hist2D> Hist2D* createHist(const TMatrixD& m, const char* name, const char* title, const Variable<Hist2D>& x, const Variable<Hist2D>& y);
  template<class Hist2D> Hist2D* createHist(const TMatrixD& m, const TMatrixD& me, const char* name, const char* title, const Variable<Hist2D>& x, const Variable<Hist2D>& y);
  template<class Hist> Hist* createHist(const char* name, const char* title, const std::vector<Variable<Hist>>& x);  
  template<class Hist> Hist* createHist(const char* name, const char* title, const Variable<Hist>& x) { return createHist<Hist>(name,title,std::vector<Variable<Hist>>{x}); };  
  template<class Hist> Hist* createHist(const TVectorD& vec, const char* name, const char* title, const std::vector<Variable<Hist>>& x, bool overflow=false);
  template<class Hist> Hist* createHist(const TVectorD& vec, const char* name, const char* title, const Variable<Hist>& x, bool overflow=false);
  template<class Hist> Hist* createHist(const TVectorD& vec, const TVectorD& errvec, const char* name, const char* title, const std::vector<Variable<Hist>>& x, bool overflow=false);
  template<class Hist> Hist* createHist(const TVectorD& vec, const TVectorD& errvec, const char* name, const char* title, const Variable<Hist>& x, bool overflow=false);
  
  void printTable (std::ostream& o, const TVectorD& vTrainTrue, const TVectorD& vTrain, const TVectorD& vMeas, const TVectorD& vReco, Int_t nm, Int_t nt);
  void printTable (std::ostream& o, int dim, int ntxb, int ntyb,
                   const TVectorD& vTrainTrue, const TVectorD& vTrain, const TVectorD& vTrue,const TVectorD& vMeas, const TVectorD& vReco,
                   ErrorTreatment withError, const TVectorD& eTrue, const TVectorD& eReco, double chi_squ=-999., bool overflow=false);

  void printVector(const char* name, const TVectorD& vec);
  void printMatrix (const TMatrixD& m, const char* name="matrix", const char* format=0, Int_t cols_per_sheet=10);

  void add(TMatrixD& target, const TMatrixD& addition);
  TMatrixD& ABAT (const TMatrixD& a, const TMatrixD& b, TMatrixD& c);
  TMatrixD& ABAT (const TMatrixD& a, const TVectorD& b, TMatrixD& c);
  template<class Hist> RooUnfolding::Variable<Hist> var(const Hist* h, Dimension d);
  template<class Hist> std::vector<RooUnfolding::Variable<Hist>> vars(const Hist* h);
}

#endif
