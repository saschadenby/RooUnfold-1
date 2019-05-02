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

class TH1;
class TH2;

class RooUnfoldSvd : public RooUnfold {
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
   SVDUnfold( const TH1* bdat, const TH1* bini, const TH1* xini, const TH2* Adet );
   SVDUnfold( const TH1* bdat, TH2* Bcov, const TH1* bini, const TH1* xini, const TH2* Adet );
   SVDUnfold( const SVDUnfold& other );

   // Destructor
   virtual ~SVDUnfold(); 

   // Set option to normalize unfolded spectrum to unit area
   // "normalize" - switch 
   void     SetNormalize ( Bool_t normalize ) { fNormalize = normalize; }

   // Do the unfolding
   // "kreg"   - number of singular values used (regularisation)
   TH1*    Unfold       ( Int_t kreg );

   // Determine for given input error matrix covariance matrix of unfolded 
   // spectrum from toy simulation
   // "cov"    - covariance matrix on the measured spectrum, to be propagated
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   TH2*    GetUnfoldCovMatrix( const TH2* cov, Int_t ntoys, Int_t seed = 1 );

   // Determine covariance matrix of unfolded spectrum from finite statistics in 
   // response matrix
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   // "uncmat" - matrix containing the uncertainty on the detector matrix elements if different from purely statistical without any weights
   TH2*    GetAdetCovMatrix( Int_t ntoys, Int_t seed=1, const TH2* uncmat=0 );

   // Regularisation parameter
   Int_t    GetKReg() const { return fKReg; }

   // Obtain the distribution of |d| (for determining the regularization)
   TH1*    GetD() const;

   // Obtain the distribution of singular values
   TH1*    GetSV() const;

   // Obtain the computed regularized covariance matrix
   TH2*    GetXtau() const;

   // Obtain the computed inverse of the covariance matrix
   TH2*    GetXinv() const;
   
   //Obtain the covariance matrix on the data
   TH2*    GetBCov() const;

   // Helper functions
   Double_t ComputeChiSquared( const TH1& truspec, const TH1& unfspec );

private: 
   
   // Helper functions for vector and matrix operations
   void            FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC ) const;
   static Double_t GetCurvature       ( const TVectorD& vec, const TMatrixD& curv );

   void            InitHistos  ( );

   // Helper functions
   static void     H2V      ( const TH1* histo, TVectorD& vec   );
   static void     H2Verr   ( const TH1* histo, TVectorD& vec   );
   static void     V2H      ( const TVectorD& vec, TH1& histo   );
   static void     H2M      ( const TH2* histo, TMatrixD& mat   );
   static void     M2H      ( const TMatrixD& mat, TH2& histo   );
   static TMatrixD MatDivVec( const TMatrixD& mat, const TVectorD& vec, Int_t zero=0 );
   static TVectorD CompProd ( const TVectorD& vec1, const TVectorD& vec2 );

   static TVectorD VecDiv                 ( const TVectorD& vec1, const TVectorD& vec2, Int_t zero = 0 );
   static void     RegularisedSymMatInvert( TMatrixDSym& mat, Double_t eps = 1e-3 );
   
   // Class members
   Int_t       fNdim;        //! Truth and reconstructed dimensions
   Int_t       fDdim;        //! Derivative for curvature matrix
   Bool_t      fNormalize;   //! Normalize unfolded spectrum to 1
   Int_t       fKReg;        //! Regularisation parameter
   TH1*       fDHist;       //! Distribution of d (for checking regularization)
   TH1*       fSVHist;      //! Distribution of singular values
   TH2*       fXtau;        //! Computed regularized covariance matrix
   TH2*       fXinv;        //! Computed inverse of covariance matrix

   // Input histos
   const TH1* fBdat;        // measured distribution (data)
   TH2* fBcov;        // covariance matrix of measured distribution (data)
   const TH1* fBini;        // reconstructed distribution (MC)
   const TH1* fXini;        // truth distribution (MC)
   const TH2* fAdet;        // Detector response matrix

   // Evaluation of covariance matrices
   TH1*       fToyhisto;    //! Toy MC histogram
   TH2*       fToymat;      //! Toy MC detector response matrix
   Bool_t      fToyMode;     //! Internal switch for covariance matrix propagation
   Bool_t      fMatToyMode;  //! Internal switch for evaluation of statistical uncertainties from response matrix
};

public:

  // Standard methods

  RooUnfoldSvd(); // default constructor
  RooUnfoldSvd (const char*    name, const char*    title); // named constructor
  RooUnfoldSvd (const TString& name, const TString& title); // named constructor
  RooUnfoldSvd (const RooUnfoldSvd& rhs); // copy constructor
  virtual ~RooUnfoldSvd(); // destructor
  RooUnfoldSvd& operator= (const RooUnfoldSvd& rhs); // assignment operator
  virtual RooUnfoldSvd* Clone (const char* newname= 0) const;

  // Special constructors

  RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t kreg= 0,
                const char* name= 0, const char* title= 0);
  // compatibility constructor
  RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t kreg, Int_t ntoyssvd,
                const char* name= 0, const char* title= 0);

  void SetKterm (Int_t kreg);
  Int_t GetKterm() const;
  virtual void  SetRegParm (Double_t parm);
  virtual Double_t GetRegParm() const;
  virtual void Reset();
  SVDUnfold* Impl();

  void SetNtoysSVD (Int_t ntoyssvd);  // no longer used
  Int_t GetNtoysSVD() const;          // no longer used

protected:
  void Assign (const RooUnfoldSvd& rhs); // implementation of assignment operator
  virtual void Unfold();
  virtual void GetCov();
  virtual void GetWgt();
  virtual void GetSettings();

private:
  void Init();
  void Destroy();
  void CopyData (const RooUnfoldSvd& rhs);

protected:
  // instance variables
  SVDUnfold* _svd;  //! Implementation in TSVDUnfold object (no streamer)
  Int_t _kreg;
  Int_t _nb;

  TH1 *_meas1d, *_train1d, *_truth1d;
  TH2 *_reshist, *_meascov;

public:
  ClassDef (RooUnfoldSvd, 1) // SVD Unfolding (interface to TSVDUnfold)
};

// Inline method definitions

inline
RooUnfoldSvd::RooUnfoldSvd()
  : RooUnfold()
{
  // Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

inline
RooUnfoldSvd::RooUnfoldSvd (const char* name, const char* title)
  : RooUnfold(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

inline
RooUnfoldSvd::RooUnfoldSvd (const TString& name, const TString& title)
  : RooUnfold(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

inline
RooUnfoldSvd& RooUnfoldSvd::operator= (const RooUnfoldSvd& rhs)
{
  // Assignment operator for copying RooUnfoldSvd settings.
  Assign(rhs);
  return *this;
}

inline
RooUnfoldSvd::~RooUnfoldSvd()
{
  Destroy();
}


inline
void RooUnfoldSvd::SetKterm (Int_t kreg)
{
  // Set regularisation parameter
  _kreg= kreg;
}


inline
Int_t RooUnfoldSvd::GetKterm() const
{
  // Return regularisation parameter
  return _kreg;
}

inline void RooUnfoldSvd::SetNtoysSVD (Int_t ntoyssvd) {_NToys=ntoyssvd;}  // no longer used
inline Int_t RooUnfoldSvd::GetNtoysSVD() const { return _NToys; }  // no longer used

inline
void  RooUnfoldSvd::SetRegParm (Double_t parm)
{
  // Set regularisation parameter
  SetKterm(Int_t(parm+0.5));
}

inline
Double_t RooUnfoldSvd::GetRegParm() const
{
  // Return regularisation parameter
  return GetKterm();
}

#endif
