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

  // Special constructors

  RooUnfoldT<Hist,Hist2D> (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, const char* name= 0, const char* title= 0);

  // Set up an existing object
  
  virtual RooUnfoldT<Hist,Hist2D>& Setup (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas);
  virtual void SetMeasured (const Hist* meas);
  virtual void SetMeasured (const TVectorD& meas, const TMatrixD& cov);
  virtual void SetMeasured (const TVectorD& meas, const TVectorD& err);
  virtual void SetMeasuredCov (const TMatrixD& cov);
  virtual void SetResponse (const RooUnfoldResponseT<Hist,Hist2D>* res, Bool_t takeOwnership = false);

  virtual void Reset ();

  // Accessors

  virtual const RooUnfoldResponseT<Hist,Hist2D>* response() const;
  virtual RooUnfoldResponseT<Hist,Hist2D>* response();
  virtual const Hist* Hmeasured() const;
  virtual Hist* Hmeasured();
  virtual       Hist* Hreco (RooUnfolding::ErrorTreatment withError=RooUnfolding::kErrors);
  const    TVectorD& Vmeasured() const;   // Measured distribution as a TVectorD
  const    TVectorD& Emeasured() const;   // Measured distribution errors as a TVectorD
  const    TMatrixD& GetMeasuredCov() const;   // Measured distribution covariance matrix

  virtual const TVectorD&  Vreco() const;
  virtual TMatrixD   Ereco  (RooUnfolding::ErrorTreatment witherror=RooUnfolding::kCovariance) const;
  virtual TVectorD   ErecoV (RooUnfolding::ErrorTreatment witherror=RooUnfolding::kErrors) const;
  virtual TMatrixD   Wreco  (RooUnfolding::ErrorTreatment witherror=RooUnfolding::kCovariance) const;

  virtual Int_t      verbose() const;
  virtual void       SetVerbose (Int_t level);
  virtual void       IncludeSystematics (Int_t dosys= 1);
  virtual Int_t      SystematicsIncluded() const;
  virtual Int_t      NToys() const;         // Number of toys
  virtual void       SetNToys (Int_t toys); // Set number of toys
  virtual Int_t      Overflow() const;
  virtual void       PrintTable (std::ostream& o, const Hist* hTrue= 0, RooUnfolding::ErrorTreatment withError=RooUnfolding::kDefault) const;
  virtual void       SetRegParm (Double_t parm);
  virtual Double_t   GetRegParm() const; // Get Regularisation Parameter
  Double_t Chi2 (const Hist* hTrue,RooUnfolding::ErrorTreatment DoChi2=RooUnfolding::kCovariance) const;
  Double_t GetMinParm() const;
  Double_t GetMaxParm() const;
  Double_t GetStepSizeParm() const;
  Double_t GetDefaultParm() const;
  RooUnfoldT<Hist,Hist2D>* RunToy() const;
  void Print(Option_t* opt="") const;
  void ForceRecalculation();

protected:
  void Assign (const RooUnfoldT<Hist,Hist2D>& rhs); // implementation of assignment operator
  virtual void SetNameTitleDefault();
  virtual void Unfold() const;
  virtual Bool_t UnfoldWithErrors (RooUnfolding::ErrorTreatment withError, bool getWeights=false) const;

  static TMatrixD CutZeros     (const TMatrixD& ereco);
  static Int_t    InvertMatrix (const TMatrixD& mat, TMatrixD& inv, const char* name="matrix", Int_t verbose=1);

private:
  void Init();
  void Destroy();  
  void CopyData (const RooUnfoldT<Hist,Hist2D>& rhs);

protected:
  // cache 
  virtual void GetCov() const; // Get covariance matrix using errors on measured distribution
  virtual void GetErrors() const;
  virtual void GetSettings() const;
  virtual void GetErrMat() const; // Get covariance matrix using errors from residuals on reconstructed distribution
  virtual void GetWgt() const; // Get weight matrix using errors on measured distribution

  class Cache {
  public:
    Cache();
    ~Cache();
    typename RooUnfoldT<Hist,Hist2D>::Cache& operator=(const Cache&);
    Double_t _minparm;       // Minimum value to be used in RooUnfoldParms
    Double_t _maxparm;       // Maximum value to be used in RooUnfoldParms
    Double_t _stepsizeparm;  // StepSize value to be used in RooUnfoldParms
    Double_t _defaultparm;   // Recommended value for regularisation parameter
    Bool_t   _unfolded;      // unfolding done
    Bool_t   _haveCov;       // have _cov
    Bool_t   _haveWgt;       // have _wgt
    Bool_t   _have_err_mat;  // have _err_mat
    Bool_t   _fail;          // unfolding failed
    Bool_t   _haveErrors;    // have _variances
    TVectorD _rec;           // Reconstructed distribution
    TMatrixD _cov;           // Reconstructed distribution covariance
    TMatrixD _wgt;           // Reconstructed distribution weights (inverse of _cov)
    TVectorD _variances;     // Error matrix diagonals
    TMatrixD _err_mat;       // Error matrix from toys
    TVectorD* _vMes;         // Cached measured vector
    TVectorD* _eMes;         // Cached measured error
    TMatrixD* _covL;         // Cached lower triangular matrix for which _covMes = _covL * _covL^T.
    TMatrixD* _covMes;       // Measurement covariance matrix    
    RooUnfolding::ErrorTreatment _withError; // type of error last calulcated
  };
  mutable Cache _cache; //!

  TMatrixD* _covMes;                       // Measurement covariance matrix
  Int_t    _verbose;                       // Debug print level
  Int_t    _nm;                            // Total number of measured bins (including under/overflows if _overflow set)
  Int_t    _nt;                            // Total number of truth    bins (including under/overflows if _overflow set)
  Int_t    _overflow;                      // Use histogram under/overflows if 1 (set from RooUnfoldResponse)
  Int_t    _NToys;                         // Number of toys to be used
  Int_t    _dosys;                         // include systematic errors from response matrix? use _dosys=2 to exclude measurement errors
  RooUnfoldResponseT<Hist,Hist2D>* _res;   // Response matrix (not owned)
  Hist*    _meas;                          // Measured distribution (not owned)

public:

  ClassDefT (RooUnfoldT, 2) // Unfolding base class: implementations in RooUnfoldBayes, RooUnfoldSvd, RooUnfoldBinByBin, RooUnfoldTUnfold, RooUnfoldInvert, RooUnfoldIds
};

typedef RooUnfoldT<TH1,TH2> RooUnfold;
#ifndef NOROOFIT
#include <RooAbsPdf.h>
#include <RooAbsReal.h>

namespace RooUnfolding {
  template<class Base> class RooFitWrapper : public Base {
  protected:
    RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* _unfolding;
    double _minVal = 1e-12;
    mutable const RooArgSet* _curNormSet ; //! 
       
  public:

    const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unfolding() const ;

    virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const override;
    virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const override;
    virtual Bool_t isBinnedDistribution(const RooArgSet& obs) const override;
    virtual Double_t evaluate() const override;
    virtual TObject* clone(const char* newname = 0) const override;
    virtual Double_t getValV(const RooArgSet* set=0) const override;
    
    virtual Bool_t checkObservables(const RooArgSet *nset) const override;
    virtual Bool_t forceAnalyticalInt(const RooAbsArg &arg) const override;
    virtual Int_t getAnalyticalIntegralWN(RooArgSet &allVars, RooArgSet &numVars, const RooArgSet *normSet, const char *rangeName = 0) const override;
    virtual Double_t analyticalIntegralWN(Int_t code, const RooArgSet *normSet, const char *rangeName = 0) const override;
    virtual void printMetaArgs(std::ostream &os) const override;
    virtual RooAbsArg::CacheMode canNodeBeCached() const override;
    virtual void setCacheAndTrackHints(RooArgSet &) override;

    virtual Bool_t redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive);

    RooFitWrapper();    
    RooFitWrapper(const char* name, const char* title, const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unf);
    RooFitWrapper(const RooUnfolding::RooFitWrapper<Base>& other);
    RooFitWrapper(const RooUnfolding::RooFitWrapper<Base>* other);    
    virtual ~RooFitWrapper();
    ClassDefT(RooFitWrapper,1)
  };
}

class RooUnfoldFunc : public RooUnfolding::RooFitWrapper<RooAbsReal> {
public:
  RooUnfoldFunc();
  RooUnfoldFunc(const char* name, const char* title, const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unf);
  virtual ~RooUnfoldFunc();  
  virtual TObject* clone(const char* newname = 0) const override;
  ClassDefOverride(RooUnfoldFunc,1)
};
class RooUnfoldPdf : public RooUnfolding::RooFitWrapper<RooAbsPdf> {
public:
  RooUnfoldPdf();
  RooUnfoldPdf(const char* name, const char* title, const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unf);    
  virtual RooAbsPdf::ExtendMode extendMode() const override;
  virtual Double_t expectedEvents(const RooArgSet* nset) const override;
  virtual Double_t expectedEvents(const RooArgSet& nset) const override;
  virtual Bool_t selfNormalized() const override;
  virtual ~RooUnfoldPdf();    
  virtual TObject* clone(const char* newname = 0) const override;
  ClassDefOverride(RooUnfoldPdf,1)
};

#endif
#endif
