#ifndef ROOUNFOLDBLOBEL_HH
#define ROOUNFOLDBLOBEL_HH

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TMatrixD
#include "TMatrixD.h"
#endif
#ifndef ROOT_TVectorD
#include "TVectorD.h"
#endif
#ifndef ROOT_TMatrixDSym
#include "TMatrixDSym.h"
#endif
// #include "RooUnfoldHelpers.h"

class RooUnfoldResponse;
class TH1;
class TH1D;
class TH2D;

class RooUnfoldBlobel : public RooUnfold {


public:

  RooUnfoldBlobel(); // default constructor
  RooUnfoldBlobel (const char*    name, const char*    title); // named constructor
  RooUnfoldBlobel (const TString& name, const TString& title); // named constructor
  RooUnfoldBlobel (const RooUnfoldBlobel& rhs); // copy constructor
  virtual ~RooUnfoldBlobel(); // destructor
  RooUnfoldBlobel& operator= (const RooUnfoldBlobel& rhs); // assignment operator
  virtual RooUnfoldBlobel* Clone (const char* newname= 0) const;



    RooUnfoldBlobel (const RooUnfoldResponse* res, const TH1* meas, Int_t kreg= 0,
                  const char* name= 0, const char* title= 0);
    // compatibility constructor
    RooUnfoldBlobel (const RooUnfoldResponse* res, const TH1* meas, Int_t kreg, Int_t ntoyssvd,
                  const char* name= 0, const char* title= 0);


    void SetKterm (Int_t kreg);
    Int_t GetKterm() const;
    virtual void  SetRegParm (Double_t parm);
    virtual Double_t GetRegParm() const;
    virtual void Reset();




protected:
  void Assign (const RooUnfoldBlobel& rhs);
  // virtual void PrepareHistograms();
   TMatrixD GetHess(Int_t _nb, TH1D *_meas1d, TMatrixD reshistmatrix, TVectorD *_est);
   TVectorD GetGrad(Int_t _nb, TH1D *_meas1d, TMatrixD reshistmatrix, TVectorD *_est);
   Double_t GetLoss(Int_t _nb, TH1D *_meas1d, TMatrixD reshistmatrix, TVectorD *_est);
  virtual void Unfold();
  virtual void GetCov();
  virtual void GetWgt();
  virtual void GetSettings();
  void FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC, Int_t _nb) const;

private:
  void Init();
  void Destroy();
  void CopyData (const RooUnfoldBlobel& rhs);

protected:
  //instance variables
  // TBlobelUnfold* _bru;
  // Class members
  Int_t       _kreg;        //! Reconstructed dimensions
  Int_t       _nb;        //! Truth dimensions

  // Input histos    Bool_t      fMatToyMode;  //! Internal switch for evaluation of statistical uncertainties from response matrix
  TH1D *_meas1d, *_train1d, *_truth1d;
  TH2D *_reshist, *_meascov;

  //Matrix to be used
  TH2D *Hessian;

public:
  ClassDef (RooUnfoldBlobel, 1)

};
inline
RooUnfoldBlobel::RooUnfoldBlobel()
  : RooUnfold()
{
  // Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

inline
RooUnfoldBlobel::RooUnfoldBlobel (const char* name, const char* title)
  : RooUnfold(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

inline
RooUnfoldBlobel::RooUnfoldBlobel (const TString& name, const TString& title)
  : RooUnfold(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

inline
RooUnfoldBlobel& RooUnfoldBlobel::operator= (const RooUnfoldBlobel& rhs)
{
  // Assignment operator for copying RooUnfoldBlobel settings.
  Assign(rhs);
  return *this;
}

inline
RooUnfoldBlobel::~RooUnfoldBlobel()
{
  Destroy();
}


inline
void RooUnfoldBlobel::SetKterm (Int_t kreg)
{
  // Set regularisation parameter
  _kreg= kreg;
}


inline
Int_t RooUnfoldBlobel::GetKterm() const
{
  // Return regularisation parameter
  return _kreg;
}

// inline void RooUnfoldBlobel::SetNtoysSVD (Int_t ntoyssvd) {_NToys=ntoyssvd;}  // no longer used
// inline Int_t RooUnfoldBlobel::GetNtoysSVD() const { return _NToys; }  // no longer used

inline
void  RooUnfoldBlobel::SetRegParm (Double_t parm)
{
  // Set regularisation parameter
  SetKterm(Int_t(parm+0.5));
}

inline
Double_t RooUnfoldBlobel::GetRegParm() const
{
  // Return regularisation parameter
  return GetKterm();
}

#endif
