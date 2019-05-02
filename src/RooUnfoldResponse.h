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

template<class Hist>
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
  RooUnfoldResponseT(const TH1* measured, const TH1* truth, const char* name= 0, const char* title= 0);  // constructor - measured and truth only used for shape
  RooUnfoldResponseT(const TH1* measured, const TH1* truth, const TH2* response, const char* name= 0, const char* title= 0);  // create from already-filled histograms

  // Set up an existing object

  virtual RooUnfoldResponseT& Reset ();  // clear an existing object
  virtual RooUnfoldResponseT& Setup (const RooUnfoldResponseT& rhs);  // set up based on another instance
  virtual RooUnfoldResponseT& Setup (Int_t nb, Double_t xlo, Double_t xhi);  // set up simple 1D case with same binning, measured vs truth
  virtual RooUnfoldResponseT& Setup (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi);  // set up simple 1D case
  virtual RooUnfoldResponseT& Setup (const TH1* measured, const TH1* truth);  // set up - measured and truth only used for shape
  virtual RooUnfoldResponseT& Setup (const TH1* measured, const TH1* truth, const TH2* response);  // set up from already-filled histograms

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

  const TH1*   Hmeasured()            const;   // Measured distribution, including fakes
  TH1*         Hmeasured();                    // Measured distribution, including fakes
  const TH1*   Hfakes()               const;   // Fakes distribution
  TH1*         Hfakes();                       // Fakes distribution
  const TH1*   Htruth()               const;   // Truth distribution, used for normalisation
  TH1*         Htruth();                       // Truth distribution, used for normalisation
  const TH2*   Hresponse()            const;   // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  TH2*         Hresponse();                    // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  TH2D*        HresponseNoOverflow()  const;   // Response matrix with under/overflow bins moved into histogram body

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

  static TH1D*     H2H1D(const TH1*  h, Int_t nb);
  static TVectorD* H2V  (const TH1*  h, Int_t nb, Bool_t overflow= kFALSE);
  static TVectorD* H2VE (const TH1*  h, Int_t nb, Bool_t overflow= kFALSE);
  static TMatrixD* H2M  (const TH2*  h, Int_t nx, Int_t ny, const TH1* norm= 0, Bool_t overflow= kFALSE);
  static TMatrixD* H2ME (const TH2*  h, Int_t nx, Int_t ny, const TH1* norm= 0, Bool_t overflow= kFALSE);
  static void      V2H  (const TVectorD& v, TH1* h, Int_t nb, Bool_t overflow= kFALSE);
  static Int_t   FindBin(const TH1*  h, Double_t x);  // return vector index for bin containing (x)
  static Int_t   FindBin(const TH1*  h, Double_t x, Double_t y);  // return vector index for bin containing (x,y)
  static Int_t   FindBin(const TH1*  h, Double_t x, Double_t y, Double_t z);  // return vector index for bin containing (x,y,z)
  static Int_t   GetBin (const TH1*  h, Int_t i, Bool_t overflow= kFALSE);  // vector index (0..nx*ny-1) -> multi-dimensional histogram global bin number (0..(nx+2)*(ny+2)-1) skipping under/overflow bins
  static Double_t GetBinContent (const TH1* h, Int_t i, Bool_t overflow= kFALSE); // Bin content by vector index
  static Double_t GetBinError   (const TH1* h, Int_t i, Bool_t overflow= kFALSE); // Bin error   by vector index
  static void PrintMatrix (const TMatrixD& m, const char* name="matrix", const char* format=0, Int_t cols_per_sheet=10);

  TH1* ApplyToTruth (const TH1* truth= 0, const char* name= "AppliedResponse") const; // If argument is 0, applies itself to its own truth
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

  static Int_t GetBinDim (const TH1* h, Int_t i);
  static void ReplaceAxis(TAxis* axis, const TAxis* source);

  // instance variables

  Int_t _mdim;     // Number of measured  dimensions
  Int_t _tdim;     // Number of truth     dimensions
  Int_t _nm;       // Total number of measured  bins (not counting under/overflows)
  Int_t _nt;       // Total number of truth     bins (not counting under/overflows)
  TH1*  _mes;      // Measured histogram
  TH1*  _fak;      // Fakes    histogram
  TH1*  _tru;      // Truth    histogram
  TH2*  _res;      // Response histogram
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

// Inline method definitions

template<class Hist>
RooUnfoldResponseT<Hist>::RooUnfoldResponseT()
  : TNamed()
{
  // RooUnfoldResponseT<Hist> default constructor. Use Setup() to set values.
  Init();
}

template<class Hist>
RooUnfoldResponseT<Hist>::RooUnfoldResponseT(const char*    name, const char*    title)
  : TNamed(name,title)
{
  // RooUnfoldResponseT<Hist> default named constructor. Use Setup() to set values.
  Init();
}

template<class Hist>
RooUnfoldResponseT<Hist>::RooUnfoldResponseT(const TString& name, const TString& title)
  : TNamed(name,title)
{
  // RooUnfoldResponseT<Hist> default named constructor. Use Setup() to set values.
  Init();
}

template<class Hist>
RooUnfoldResponseT<Hist>::~RooUnfoldResponseT()
{
  // RooUnfoldResponseT<Hist> destructor
  Reset();
}


template<class Hist>
RooUnfoldResponseT<Hist>& RooUnfoldResponseT<Hist>::Setup (Int_t nb, Double_t xlo, Double_t xhi)
{
  // constructor -  simple 1D case with same binning, measured vs truth
  return Setup (nb, xlo, xhi, nb, xlo, xhi);
}


template<class Hist>
Int_t RooUnfoldResponseT<Hist>::GetDimensionMeasured() const
{
  // Dimensionality of the measured distribution (1=1D, 2=2D, 3=3D)
  return _mdim;
}

template<class Hist>
Int_t RooUnfoldResponseT<Hist>::GetDimensionTruth() const
{
  // Dimensionality of the truth distribution (1=1D, 2=2D, 3=3D)
  return _tdim;
}

template<class Hist>
Int_t RooUnfoldResponseT<Hist>::GetNbinsMeasured() const
{
  // Total number of bins in the measured distribution
  return _nm;
}

template<class Hist>
Int_t RooUnfoldResponseT<Hist>::GetNbinsTruth() const
{
  // Total number of bins in the truth distribution
  return _nt;
}


template<class Hist>
const TH1* RooUnfoldResponseT<Hist>::Hmeasured() const
{
  // Measured distribution, including fakes
  return _mes;
}


template<class Hist>
TH1*         RooUnfoldResponseT<Hist>::Hmeasured()
{
  return _mes;
}


template<class Hist>
const TH1* RooUnfoldResponseT<Hist>::Hfakes() const
{
  // Fakes distribution
  return _fak;
}


template<class Hist>
TH1*         RooUnfoldResponseT<Hist>::Hfakes()
{
  return _fak;
}

template<class Hist>
const TH1*   RooUnfoldResponseT<Hist>::Htruth() const
{
  // Truth distribution, used for normalisation
  return _tru;
}

template<class Hist>
TH1*         RooUnfoldResponseT<Hist>::Htruth()
{
  return _tru;
}

template<class Hist>
const TH2*   RooUnfoldResponseT<Hist>::Hresponse() const
{
  // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  return _res;
}

template<class Hist>
TH2*         RooUnfoldResponseT<Hist>::Hresponse()
{
  return _res;
}


template<class Hist>
const TVectorD& RooUnfoldResponseT<Hist>::Vmeasured() const
{
  // Measured distribution as a TVectorD
  if (!_vMes) _cached= (_vMes= H2V  (_mes, _nm, _overflow));
  return *_vMes;
}

template<class Hist>
const TVectorD& RooUnfoldResponseT<Hist>::Vfakes() const
{
  // Fakes distribution as a TVectorD
  if (!_vFak) _cached= (_vFak= H2V  (_fak, _nm, _overflow));
  return *_vFak;
}

template<class Hist>
const TVectorD& RooUnfoldResponseT<Hist>::Emeasured() const
{
  // Measured distribution errors as a TVectorD
  if (!_eMes) _cached= (_eMes= H2VE (_mes, _nm, _overflow));
  return *_eMes;
}

template<class Hist>
const TVectorD& RooUnfoldResponseT<Hist>::Vtruth() const
{
  // Truth distribution as a TVectorD
  if (!_vTru) _cached= (_vTru= H2V  (_tru, _nt, _overflow)); 
  return *_vTru;
}

template<class Hist>
const TVectorD& RooUnfoldResponseT<Hist>::Etruth() const
{
  // Truth distribution errors as a TVectorD
  if (!_eTru) _cached= (_eTru= H2VE (_tru, _nt, _overflow)); 
  return *_eTru;
}

template<class Hist>
const TMatrixD& RooUnfoldResponseT<Hist>::Mresponse() const
{
  // Response matrix as a TMatrixD: (row,column)=(measured,truth)
  if (!_mRes) _cached= (_mRes= H2M  (_res, _nm, _nt, _tru, _overflow)); 
  return *_mRes;
}

template<class Hist>
const TMatrixD& RooUnfoldResponseT<Hist>::Eresponse() const
{
  // Response matrix errors as a TMatrixD: (row,column)=(measured,truth)
  if (!_eRes) _cached= (_eRes= H2ME (_res, _nm, _nt, _tru, _overflow)); 
  return *_eRes;
}


template<class Hist>
Double_t RooUnfoldResponseT<Hist>::operator() (Int_t r, Int_t t) const
{
  // Response matrix element (measured,truth)
  return Mresponse()(r,t);
}

template<class Hist>
Int_t    RooUnfoldResponseT<Hist>::GetBin (const TH1* h, Int_t i, Bool_t overflow)
{
  // vector index (0..nx*ny-1) -> multi-dimensional histogram
  // global bin number (0..(nx+2)*(ny+2)-1) skipping under/overflow bins
  return (h->GetDimension()<2) ? i+(overflow ? 0 : 1) : GetBinDim(h,i);
}

template<class Hist>
Double_t RooUnfoldResponseT<Hist>::GetBinContent (const TH1* h, Int_t i, Bool_t overflow)
{
  // Bin content by vector index
  return h->GetBinContent (GetBin (h, i, overflow));
}

template<class Hist>
Double_t RooUnfoldResponseT<Hist>::GetBinError   (const TH1* h, Int_t i, Bool_t overflow)
{
  // Bin error   by vector index
  return h->GetBinError   (GetBin (h, i, overflow));
}


template<class Hist>
Int_t RooUnfoldResponseT<Hist>::Miss (Double_t xt)
{
  // Fill missed event into 1D Response Matrix
  return Miss1D(xt);
}

template<class Hist>
Int_t RooUnfoldResponseT<Hist>::Miss (Double_t xt, Double_t w)
{
  // Fill missed event into 1D (with weight) or 2D Response Matrix
  return _tdim==2 ? Miss2D(xt,w) : Miss1D(xt,w);
}

template<class Hist>
Int_t RooUnfoldResponseT<Hist>::Miss (Double_t xt, Double_t yt, Double_t w)
{
  // Fill missed event into 2D (with weight) or 3D Response Matrix
  return _tdim==3 ? Miss(xt,yt,w,1.0) : Miss2D(xt,yt,w);
}


template<class Hist>
Int_t RooUnfoldResponseT<Hist>::Fake (Double_t xr)
{
  // Fill fake event into 1D Response Matrix
  return Fake1D(xr);
}

template<class Hist>
Int_t RooUnfoldResponseT<Hist>::Fake (Double_t xr, Double_t w)
{
  // Fill fake event into 1D (with weight) or 2D Response Matrix
  return _mdim==2 ? Fake2D(xr,w) : Fake1D(xr,w);
}

template<class Hist>
Int_t RooUnfoldResponseT<Hist>::Fake (Double_t xr, Double_t yr, Double_t w)
{
  // Fill fake event into 2D (with weight) or 3D Response Matrix
  return _mdim==3 ? Fake(xr,yr,w,1.0) : Fake2D(xr,yr,w);
}


template<class Hist>
void RooUnfoldResponseT<Hist>::UseOverflow (Bool_t set)
{
  // Specify to use overflow bins. Only supported for 1D truth and measured distributions.
  _overflow= (set ? 1 : 0);
}

template<class Hist>
Bool_t RooUnfoldResponseT<Hist>::UseOverflowStatus() const
{
  // Get UseOverflow setting
  return _overflow;
}

template<class Hist>
Double_t RooUnfoldResponseT<Hist>::FakeEntries() const
{
  // Return number of fake entries
  return _fak ? _fak->GetEntries() : 0.0;
}

template<class Hist>
Int_t RooUnfoldResponseT<Hist>::FindBin (const TH1* h, Double_t x)
{
  // Return vector index for bin containing (x)
  return h->GetXaxis()->FindBin(x) - 1;
}

typedef RooUnfoldResponseT<TH1> RooUnfoldResponse;

#endif
