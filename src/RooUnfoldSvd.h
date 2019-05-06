//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      SVD unfolding. Just an interface to RooUnfHistoSvd.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDSVD_HH
#define ROOUNFOLDSVD_HH

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

template<class Hist, class Hist2D>
class RooUnfoldSvdT : public RooUnfoldT<Hist,Hist2D> {
public:
  class SVDUnfold {

public:

   // Constructor
   // Initialisation of unfolding
   // "bdat" - measured data distribution (number of events)
   // "Bcov" - covariance matrix for measured data distribution
   // "bini" - reconstructed MC distribution (number of events)
   // "xini" - truth MC distribution (number of events)
   // "Adet" - detector response matrix (number of events)
   SVDUnfold( const Hist* bdat, const Hist* bini, const Hist* xini, const Hist2D* Adet );
   SVDUnfold( const Hist* bdat, Hist2D* Bcov, const Hist* bini, const Hist* xini, const Hist2D* Adet );
   SVDUnfold( const SVDUnfold& other );

   // Destructor
   virtual ~SVDUnfold(); 

   // Set option to normalize unfolded spectrum to unit area
   // "normalize" - switch 
   void     SetNormalize ( Bool_t normalize ) { fNormalize = normalize; }

   // Do the unfolding
   // "kreg"   - number of singular values used (regularisation)
   Hist*    Unfold       ( Int_t kreg );

   // Determine for given input error matrix covariance matrix of unfolded 
   // spectrum from toy simulation
   // "cov"    - covariance matrix on the measured spectrum, to be propagated
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   Hist2D*    GetUnfoldCovMatrix( const Hist2D* cov, Int_t ntoys, Int_t seed = 1 );

   // Determine covariance matrix of unfolded spectrum from finite statistics in 
   // response matrix
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   // "uncmat" - matrix containing the uncertainty on the detector matrix elements if different from purely statistical without any weights
   Hist2D*    GetAdetCovMatrix( Int_t ntoys, Int_t seed=1, const Hist2D* uncmat=0 );

   // Regularisation parameter
   Int_t    GetKReg() const { return fKReg; }

   // Obtain the distribution of |d| (for determining the regularization)
   Hist*    GetD() const;

   // Obtain the distribution of singular values
   Hist*    GetSV() const;

   // Obtain the computed regularized covariance matrix
   Hist2D*    GetXtau() const;

   // Obtain the computed inverse of the covariance matrix
   Hist2D*    GetXinv() const;
   
   //Obtain the covariance matrix on the data
   Hist2D*    GetBCov() const;

   // Helper functions
   Double_t ComputeChiSquared( const Hist& truspec, const Hist& unfspec );

private: 
   
   // Helper functions for vector and matrix operations
   void            FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC ) const;
   static Double_t GetCurvature       ( const TVectorD& vec, const TMatrixD& curv );

   void            InitHistos  ( );

   // Helper functions
   static void     H2V      ( const Hist* histo, TVectorD& vec   );
   static void     H2Verr   ( const Hist* histo, TVectorD& vec   );
   static void     V2H      ( const TVectorD& vec, Hist& histo   );
   static void     H2M      ( const Hist2D* histo, TMatrixD& mat   );
   static void     M2H      ( const TMatrixD& mat, Hist2D& histo   );
   static TMatrixD MatDivVec( const TMatrixD& mat, const TVectorD& vec, Int_t zero=0 );
   static TVectorD CompProd ( const TVectorD& vec1, const TVectorD& vec2 );

   static TVectorD VecDiv                 ( const TVectorD& vec1, const TVectorD& vec2, Int_t zero = 0 );
   static void     RegularisedSymMatInvert( TMatrixDSym& mat, Double_t eps = 1e-3 );
   
   // Class members
   Int_t       fNdim;        //! Truth and reconstructed dimensions
   Int_t       fDdim;        //! Derivative for curvature matrix
   Bool_t      fNormalize;   //! Normalize unfolded spectrum to 1
   Int_t       fKReg;        //! Regularisation parameter
   Hist*       fDHist;       //! Distribution of d (for checking regularization)
   Hist*       fSVHist;      //! Distribution of singular values
   Hist2D*       fXtau;        //! Computed regularized covariance matrix
   Hist2D*       fXinv;        //! Computed inverse of covariance matrix

   // Input histos
   const Hist* fBdat;        // measured distribution (data)
   Hist2D* fBcov;        // covariance matrix of measured distribution (data)
   const Hist* fBini;        // reconstructed distribution (MC)
   const Hist* fXini;        // truth distribution (MC)
   const Hist2D* fAdet;        // Detector response matrix

   // Evaluation of covariance matrices
   Hist*       fToyhisto;    //! Toy MC histogram
   Hist2D*       fToymat;      //! Toy MC detector response matrix
   Bool_t      fToyMode;     //! Internal switch for covariance matrix propagation
   Bool_t      fMatToyMode;  //! Internal switch for evaluation of statistical uncertainties from response matrix
};

public:

  // Standard methods

  RooUnfoldSvdT(); // default constructor
  RooUnfoldSvdT (const char*    name, const char*    title); // named constructor
  RooUnfoldSvdT (const TString& name, const TString& title); // named constructor
  RooUnfoldSvdT (const RooUnfoldSvdT<Hist,Hist2D>& rhs); // copy constructor
  virtual ~RooUnfoldSvdT(); // destructor
  RooUnfoldSvdT<Hist,Hist2D>& operator= (const RooUnfoldSvdT<Hist,Hist2D>& rhs); // assignment operator
  virtual RooUnfoldSvdT<Hist,Hist2D>* Clone (const char* newname= 0) const;

  // Special constructors

  RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg= 0,
                const char* name= 0, const char* title= 0);
  // compatibility constructor
  RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg, Int_t ntoyssvd,
                const char* name= 0, const char* title= 0);

  void SetKterm (Int_t kreg);
  Int_t GetKterm() const;
  virtual void  SetRegParm (Double_t parm);
  virtual Double_t GetRegParm() const;
  virtual void Reset();
  SVDUnfold* Impl();

protected:
  void Assign (const RooUnfoldSvdT<Hist,Hist2D>& rhs); // implementation of assignment operator
  virtual void Unfold();
  virtual void GetCov();
  virtual void GetWgt();
  virtual void GetSettings();

private:
  void Init();
  void Destroy();
  void CopyData (const RooUnfoldSvdT<Hist,Hist2D>& rhs);

protected:
  // instance variables
  SVDUnfold* _svd;  //! Implementation in TSVDUnfold object (no streamer)
  Int_t _kreg;
  Int_t _nb;

  Hist *_meas1d, *_train1d, *_truth1d;
  Hist2D *_reshist, *_meascov;

public:
  ClassDefT (RooUnfoldSvdT, 0) // SVD Unfolding (interface to TSVDUnfold)
};

typedef RooUnfoldSvdT<TH1,TH2> RooUnfoldSvd;

#endif
