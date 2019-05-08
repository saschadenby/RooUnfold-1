//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Unfolding framework base class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLD_HH
#define ROOUNFOLD_HH

#include "TNamed.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldHelpers.h"

class TH1;

template<class Hist, class Hist2D>
class RooUnfoldT : public TNamed {
public:
  static RooUnfoldT<Hist,Hist2D>* New (RooUnfolding::Algorithm alg, const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Double_t regparm= -1e30,
                                       const char* name= 0, const char* title= 0);

  // typedefs to ensure compatibility with legacy code
  
  typedef RooUnfolding::Algorithm Algorithm;
  typedef RooUnfolding::ErrorTreatment ErrorTreatment;
  static const Algorithm kNone;
  static const Algorithm kBayes;
  static const Algorithm kSVD;
  static const Algorithm kBinByBin;
  static const Algorithm kTUnfold;
  static const Algorithm kInvert;
  static const Algorithm kDagostini;
  static const Algorithm kIDS;
  static const ErrorTreatment kNoError;
  static const ErrorTreatment kErrors;
  static const ErrorTreatment kCovariance;
  static const ErrorTreatment kCovToy;
  static const ErrorTreatment kDefault;
  
  // Standard methods

  RooUnfoldT(); // default constructor
  RooUnfoldT (const char*    name, const char*    title); // named constructor
  RooUnfoldT (const TString& name, const TString& title); // named constructor
  RooUnfoldT (const RooUnfoldT<Hist,Hist2D>& rhs); // copy constructor
  virtual ~RooUnfoldT(); // destructor
  RooUnfoldT<Hist,Hist2D>& operator= (const RooUnfoldT<Hist,Hist2D>& rhs); // assignment operator
  virtual RooUnfoldT<Hist,Hist2D>* Clone (const char* newname= 0) const;

  // Special constructors

  RooUnfoldT<Hist,Hist2D> (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, const char* name= 0, const char* title= 0);

  // Set up an existing object
  
  virtual RooUnfoldT<Hist,Hist2D>& Setup (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas);
  virtual void SetMeasured (const Hist* meas);
  virtual void SetMeasured (const TVectorD& meas, const TMatrixD& cov);
  virtual void SetMeasured (const TVectorD& meas, const TVectorD& err);
  virtual void SetMeasuredCov (const TMatrixD& cov);
  virtual void SetResponse (const RooUnfoldResponseT<Hist,Hist2D>* res);
  virtual void SetResponse (RooUnfoldResponseT<Hist,Hist2D>* res, Bool_t takeOwnership);

  virtual void Reset ();

  // Accessors

  virtual const RooUnfoldResponseT<Hist,Hist2D>* response() const;
  virtual const Hist* Hmeasured() const;
  virtual       Hist* Hreco (RooUnfolding::ErrorTreatment withError=RooUnfolding::kErrors);
  const    TVectorD& Vmeasured() const;   // Measured distribution as a TVectorD
  const    TVectorD& Emeasured() const;   // Measured distribution errors as a TVectorD
  const    TMatrixD& GetMeasuredCov() const;   // Measured distribution covariance matrix

  virtual TVectorD&  Vreco();
  virtual TMatrixD   Ereco  (RooUnfolding::ErrorTreatment witherror=RooUnfolding::kCovariance);
  virtual TVectorD   ErecoV (RooUnfolding::ErrorTreatment witherror=RooUnfolding::kErrors);
  virtual TMatrixD   Wreco  (RooUnfolding::ErrorTreatment witherror=RooUnfolding::kCovariance);

  virtual Int_t      verbose() const;
  virtual void       SetVerbose (Int_t level);
  virtual void       IncludeSystematics (Int_t dosys= 1);
  virtual Int_t      SystematicsIncluded() const;
  virtual Int_t      NToys() const;         // Number of toys
  virtual void       SetNToys (Int_t toys); // Set number of toys
  virtual Int_t      Overflow() const;
  virtual void       PrintTable (std::ostream& o, const Hist* hTrue= 0, RooUnfolding::ErrorTreatment withError=RooUnfolding::kDefault);
  virtual void       SetRegParm (Double_t parm);
  virtual Double_t   GetRegParm() const; // Get Regularisation Parameter
  Double_t Chi2 (const Hist* hTrue,RooUnfolding::ErrorTreatment DoChi2=RooUnfolding::kCovariance);
  Double_t GetMinParm() const;
  Double_t GetMaxParm() const;
  Double_t GetStepSizeParm() const;
  Double_t GetDefaultParm() const;
  RooUnfoldT<Hist,Hist2D>* RunToy() const;
  void Print(Option_t* opt="") const;

protected:
  void Assign (const RooUnfoldT<Hist,Hist2D>& rhs); // implementation of assignment operator
  virtual void SetNameTitleDefault();
  virtual void Unfold();
  virtual void GetErrors();
  virtual void GetCov(); // Get covariance matrix using errors on measured distribution
  virtual void GetErrMat(); // Get covariance matrix using errors from residuals on reconstructed distribution
  virtual void GetWgt(); // Get weight matrix using errors on measured distribution
  virtual void GetSettings();
  virtual Bool_t UnfoldWithErrors (RooUnfolding::ErrorTreatment withError, bool getWeights=false);

  static TMatrixD CutZeros     (const TMatrixD& ereco);
  static Hist*    HistNoOverflow (const Hist* h, Bool_t overflow);
  static Int_t    InvertMatrix (const TMatrixD& mat, TMatrixD& inv, const char* name="matrix", Int_t verbose=1);

private:
  void Init();
  void Destroy();
  void CopyData (const RooUnfoldT<Hist,Hist2D>& rhs);

protected:
  // instance variables
  Double_t _minparm;       // Minimum value to be used in RooUnfoldParms
  Double_t _maxparm;       // Maximum value to be used in RooUnfoldParms
  Double_t _stepsizeparm;  // StepSize value to be used in RooUnfoldParms
  Double_t _defaultparm;   // Recommended value for regularisation parameter
  Int_t    _verbose;       // Debug print level
  Int_t    _nm;            // Total number of measured bins (including under/overflows if _overflow set)
  Int_t    _nt;            // Total number of truth    bins (including under/overflows if _overflow set)
  Int_t    _overflow;      // Use histogram under/overflows if 1 (set from RooUnfoldResponse)
  Int_t    _NToys;         // Number of toys to be used
  Bool_t   _unfolded;      // unfolding done
  Bool_t   _haveCov;       // have _cov
  Bool_t   _haveWgt;       // have _wgt
  Bool_t   _have_err_mat;  // have _err_mat
  Bool_t   _fail;          // unfolding failed
  Bool_t   _haveErrors;    // have _variances
  Bool_t   _haveCovMes;    // _covMes was set, not just cached
  Int_t    _dosys;         // include systematic errors from response matrix? use _dosys=2 to exclude measurement errors
  const RooUnfoldResponseT<Hist,Hist2D>* _res;   // Response matrix (not owned)
  RooUnfoldResponseT<Hist,Hist2D>* _resmine;     // Owned response matrix
  const Hist*               _meas;  // Measured distribution (not owned)
  Hist*     _measmine;      // Owned measured histogram
  TVectorD _rec;           // Reconstructed distribution
  TMatrixD _cov;           // Reconstructed distribution covariance
  TMatrixD _wgt;           // Reconstructed distribution weights (inverse of _cov)
  TVectorD _variances;     // Error matrix diagonals
  TMatrixD _err_mat;       // Error matrix from toys
  mutable TVectorD* _vMes; //! Cached measured vector
  mutable TVectorD* _eMes; //! Cached measured error
  mutable TMatrixD* _covMes;       // Measurement covariance matrix
  mutable TMatrixD* _covL; //! Cached lower triangular matrix for which _covMes = _covL * _covL^T.
  RooUnfolding::ErrorTreatment _withError; // type of error last calulcated

public:

  ClassDefT (RooUnfoldT, 2) // Unfolding base class: implementations in RooUnfoldBayes, RooUnfoldSvd, RooUnfoldBinByBin, RooUnfoldTUnfold, RooUnfoldInvert, RooUnfoldIds
};

typedef RooUnfoldT<TH1,TH2> RooUnfold;
typedef RooUnfoldT<RooAbsReal,RooAbsReal> RooFitUnfold;
#endif
