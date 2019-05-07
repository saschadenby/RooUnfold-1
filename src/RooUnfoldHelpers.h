#ifndef ROOUNFOLDHELPERS_HH
#define ROOUNFOLDHELPERS_HH

#include <TVectorD.h>

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

  template<class Hist> Hist* createHist(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname);  
  template<class Hist2D> Hist2D* createHist(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, int nbinsy, double ymin, double ymax, const char* yname);
  template<class Hist> Hist* createHist(const TVectorD& vec, const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, bool overflow=false);
  void printTable (std::ostream& o, const TVectorD& vTrainTrue, const TVectorD& vTrain,
                   const TVectorD& vMeas, const TVectorD& vReco, Int_t nm, Int_t nt);

  void printVector(const char* name, const TVectorD& vec);
  void printMatrix (const TMatrixD& m, const char* name="matrix", const char* format=0, Int_t cols_per_sheet=10);

  void add(TMatrixD& target, const TMatrixD& addition);
  
}

#endif
