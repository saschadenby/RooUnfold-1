//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Response Matrix
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDRESPONSE_HH
#define ROOUNFOLDRESPONSE_HH

#include "TNamed.h"
#include "TMatrixD.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "TH1.h"
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
#include "TVectorDfwd.h"
#else
class TVectorD;
#endif
class TF1;
class TH2;
class TH2D;
class TAxis;
class TCollection;

#ifdef PrintMatrix
// TMVA in ROOT 6.14/00 added a debugging macro called PrintMatrix in TMVA/DNN/Architectures/Cpu/CpuMatrix.h.
// This interferes with RooUnfoldResponse::PrintMatrix below when the RooUnfoldResponse dictionary is loaded
// (eg. if used in CLING or PyROOT). This happens despite RooUnfold and TMVA being entirely independent,
// but presumably TMVA's macro writes all over the global modulemap.
// The following #undef seems to work round the issue.
// See https://sft.its.cern.ch/jira/browse/ROOT-9799 .
#undef PrintMatrix
#endif

template<class Hist, class Hist2D>
class RooUnfoldResponseT : public TNamed {

public:
  // Standard methods

  RooUnfoldResponseT(); // default constructor
  RooUnfoldResponseT(const char* name, const char* title); // named constructor
  RooUnfoldResponseT(const TString& name, const TString& title); // named constructor
  RooUnfoldResponseT(const RooUnfoldResponseT& rhs); // copy constructor
  virtual ~RooUnfoldResponseT(); // destructor
  virtual RooUnfoldResponseT& operator= (const RooUnfoldResponseT& rhs); // assignment operator

  // Special constructors

  RooUnfoldResponseT(Int_t nb, Double_t xlo, Double_t xhi, const char* name= 0, const char* title= 0);  // constructor -  simple 1D case with same binning, measured vs truth
  RooUnfoldResponseT(Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi, const char* name= 0, const char* title= 0);  // constructor -  simple 1D case
  RooUnfoldResponseT(const Hist* measured, const Hist* truth, const char* name= 0, const char* title= 0);  // constructor - measured and truth only used for shape
  RooUnfoldResponseT(const Hist* measured, const Hist* truth, const Hist2D* response, const char* name= 0, const char* title= 0);  // create from already-filled histograms

  // Set up an existing object

  virtual RooUnfoldResponseT& Reset ();  // clear an existing object
  virtual RooUnfoldResponseT& Setup (const RooUnfoldResponseT& rhs);  // set up based on another instance
  virtual RooUnfoldResponseT& Setup (Int_t nb, Double_t xlo, Double_t xhi);  // set up simple 1D case with same binning, measured vs truth
  virtual RooUnfoldResponseT& Setup (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi);  // set up simple 1D case
  virtual RooUnfoldResponseT& Setup (const Hist* measured, const Hist* truth);  // set up - measured and truth only used for shape
  virtual RooUnfoldResponseT& Setup (const Hist* measured, const Hist* truth, const Hist2D* response);  // set up from already-filled histograms

  // Fill with training data

  virtual Int_t Fill (Double_t xr, Double_t xt, Double_t w= 1.0);  // Fill 1D Response Matrix
  virtual Int_t Fill (Double_t xr, Double_t yr, Double_t xt, Double_t yt, Double_t w= 1.0);  // Fill 2D Response Matrix
  virtual Int_t Fill (Double_t xr, Double_t yr, Double_t zr, Double_t xt, Double_t yt, Double_t zt, Double_t w= 1.0);  // Fill 3D Response Matrix

          Int_t Miss (Double_t xt);  // Fill missed event into 1D Response Matrix
          Int_t Miss (Double_t xt, Double_t w);  // Fill missed event into 1D (with weight) or 2D Response Matrix
          Int_t Miss (Double_t xt, Double_t yt, Double_t w);  // Fill missed event into 2D (with weight) or 3D Response Matrix
  virtual Int_t Miss (Double_t xt, Double_t yt, Double_t zt, Double_t w);  // Fill missed event into 3D Response Matrix

          Int_t Fake (Double_t xr);  // Fill fake event into 1D Response Matrix
          Int_t Fake (Double_t xr, Double_t w);  // Fill fake event into 1D (with weight) or 2D Response Matrix
          Int_t Fake (Double_t xr, Double_t yr, Double_t w);  // Fill fake event into 2D (with weight) or 3D Response Matrix
  virtual Int_t Fake (Double_t xr, Double_t yr, Double_t zr, Double_t w);  // Fill fake event into 3D Response Matrix

  virtual void Add (const RooUnfoldResponseT& rhs);
  virtual Long64_t Merge (TCollection* others);

  // Accessors

  Int_t        GetDimensionMeasured() const;   // Dimensionality of the measured distribution
  Int_t        GetDimensionTruth()    const;   // Dimensionality of the truth distribution
  Int_t        GetNbinsMeasured()     const;   // Total number of bins in the measured distribution
  Int_t        GetNbinsTruth()        const;   // Total number of bins in the truth distribution

  const Hist*   Hmeasured()            const;   // Measured distribution, including fakes
  Hist*         Hmeasured();                    // Measured distribution, including fakes
  const Hist*   Hfakes()               const;   // Fakes distribution
  Hist*         Hfakes();                       // Fakes distribution
  const Hist*   Htruth()               const;   // Truth distribution, used for normalisation
  Hist*         Htruth();                       // Truth distribution, used for normalisation
  const Hist2D* Hresponse()            const;   // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  Hist2D*       Hresponse();                    // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  Hist2D*       HresponseNoOverflow()  const;   // Response matrix with under/overflow bins moved into histogram body

  const TVectorD& Vmeasured()         const;   // Measured distribution as a TVectorD
  const TVectorD& Emeasured()         const;   // Measured distribution errors as a TVectorD
  const TVectorD& Vfakes()            const;   // Fakes distribution as a TVectorD
  const TVectorD& Vtruth()            const;   // Truth distribution as a TVectorD
  const TVectorD& Etruth()            const;   // Truth distribution errors as a TVectorD
  const TMatrixD& Mresponse()         const;   // Response matrix as a TMatrixD: (row,column)=(measured,truth)
  const TMatrixD& Eresponse()         const;   // Response matrix errors as a TMatrixD: (row,column)=(measured,truth)

  Double_t operator() (Int_t r, Int_t t) const;// Response matrix element (measured,truth)

  void   UseOverflow (Bool_t set= kTRUE);      // Specify to use overflow bins
  Bool_t UseOverflowStatus() const;            // Get UseOverflow setting
  Double_t FakeEntries() const;                // Return number of bins with fakes
  virtual void Print (Option_t* option="") const;

  template < typename = typename std::enable_if< !std::is_same<Hist,Hist2D>::value > >
  static Hist*     H2H1D(const Hist2D*  h, Int_t nb);
  static Hist*     H2H1D(const Hist*  h, Int_t nb);
  static TVectorD* H2V  (const Hist*  h, Int_t nb, Bool_t overflow= kFALSE);
  static TVectorD* H2VE (const Hist*  h, Int_t nb, Bool_t overflow= kFALSE);
  static TMatrixD* H2M  (const Hist2D*  h, Int_t nx, Int_t ny, const Hist* norm= 0, Bool_t overflow= kFALSE);
  static TMatrixD* H2ME (const Hist2D*  h, Int_t nx, Int_t ny, const Hist* norm= 0, Bool_t overflow= kFALSE);
  static Int_t   FindBin(const Hist*  h, Double_t x);  // return vector index for bin containing (x)
  static Int_t   FindBin(const Hist*  h, Double_t x, Double_t y);  // return vector index for bin containing (x,y)
  static Int_t   FindBin(const Hist*  h, Double_t x, Double_t y, Double_t z);  // return vector index for bin containing (x,y,z)
  static void PrintMatrix (const TMatrixD& m, const char* name="matrix", const char* format=0, Int_t cols_per_sheet=10);

  Hist* ApplyToTruth (const Hist* truth= 0, const char* name= "AppliedResponse") const; // If argument is 0, applies itself to its own truth
  TF1* MakeFoldingFunction (TF1* func, Double_t eps=1e-12, Bool_t verbose=false) const;

  RooUnfoldResponseT* RunToy() const;

private:

  virtual RooUnfoldResponseT& Init();
  virtual RooUnfoldResponseT& Setup();
  virtual void ClearCache();
  virtual void SetNameTitleDefault (const char* defname= 0, const char* deftitle= 0);
  virtual Int_t Miss1D (Double_t xt, Double_t w= 1.0);  // Fill missed event into 1D Response Matrix (with weight)
  virtual Int_t Miss2D (Double_t xt, Double_t yt, Double_t w= 1.0);  // Fill missed event into 2D Response Matrix (with weight)
  virtual Int_t Fake1D (Double_t xr, Double_t w= 1.0);  // Fill fake event into 1D Response Matrix (with weight)
  virtual Int_t Fake2D (Double_t xr, Double_t yr, Double_t w= 1.0);  // Fill fake event into 2D Response Matrix (with weight)

  static Int_t GetBinDim (const Hist* h, Int_t i);

  // instance variables

  Int_t _mdim;     // Number of measured  dimensions
  Int_t _tdim;     // Number of truth     dimensions
  Int_t _nm;       // Total number of measured  bins (not counting under/overflows)
  Int_t _nt;       // Total number of truth     bins (not counting under/overflows)
  Hist*  _mes;      // Measured histogram
  Hist*  _fak;      // Fakes    histogram
  Hist*  _tru;      // Truth    histogram
  Hist2D*_res;      // Response histogram
  Int_t _overflow; // Use histogram under/overflows if 1

  mutable TVectorD* _vMes;   //! Cached measured vector
  mutable TVectorD* _eMes;   //! Cached measured error
  mutable TVectorD* _vFak;   //! Cached fakes    vector
  mutable TVectorD* _vTru;   //! Cached truth    vector
  mutable TVectorD* _eTru;   //! Cached truth    error
  mutable TMatrixD* _mRes;   //! Cached response matrix
  mutable TMatrixD* _eRes;   //! Cached response error
  mutable Bool_t    _cached; //! We are using cached vectors/matrices

public:

  ClassDefT (RooUnfoldResponseT, 1) // Respose Matrix
};

typedef RooUnfoldResponseT<TH1,TH2> RooUnfoldResponse;
typedef RooUnfoldResponseT<RooAbsReal,RooAbsReal> RooAbsRealUnfoldResponse;

#endif
