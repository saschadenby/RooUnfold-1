//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Unfolding class using inversion of the response matrix. This does not produce
//      good results and is designed to illustrate the need for more sophisticated
//      unfolding techniques
//
// Authors: Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDGP_H_
#define ROOUNFOLDGP_H_

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

class TDecompSVD;

template<class Hist, class Hist2D>
class RooUnfoldGPT : public RooUnfoldT<Hist,Hist2D> {

public:

  enum Kernel 
    {
      kRadial,
      kGibbs
    };

  RooUnfoldGPT(); // default constructor
  RooUnfoldGPT (const char*    name, const char*    title); // named constructor
  RooUnfoldGPT (const TString& name, const TString& title); // named constructor
  RooUnfoldGPT (const RooUnfoldGPT<Hist,Hist2D>& rhs); // copy constructor
  virtual ~RooUnfoldGPT(); // destructor
  RooUnfoldGPT& operator= (const RooUnfoldGPT<Hist,Hist2D>& rhs); // assignment operator
  RooUnfoldGPT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kernel = 1, const char* name=0, const char* title=0);
  virtual void SetRegParm(Double_t parm);
  virtual void Reset() override;
  TDecompSVD* Impl();

protected:
  virtual void Unfold() const override;
  //virtual void GetCov() const override;
  void MLEstimator() const;
  void MLCovariance() const;
  void MAPEstimator() const;
  void MAPCovariance() const;
  Bool_t checkGP(const TVectorD& hist) const;
private:
  void Init();
  void SetBinCenters() const;
  void SetFitSettings() const;
  Bool_t InvertResponse() const;
  TVectorD MultiplyMatrixVector(const TMatrixD& matrix, const TVectorD& vector) const;
  TVectorD MultiplyVectorMatrix(const TVectorD& vector, const TMatrixD& matrix) const;
  

  Double_t MarginalLH(const double *params) const;
  void MinimizeMLH() const;
  void SetMinInit(std::vector<Double_t>& init) const;
  void SetMinStep(std::vector<Double_t>& step) const;

  void EvaluateK(const double *params) const;
  void EvaluateKstar(const double *params) const;
  void EvaluateKstarstar(const double *params) const;

  void EvaluateK(const std::vector<Double_t>&  params) const;
  void EvaluateKstar(const std::vector<Double_t>&  params) const;
  void EvaluateKstarstar(const std::vector<Double_t>&  params) const;

  void GibbsKernel(const double *params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const;
  void RadialKernel(const double *params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const;
  void GibbsKernel(const std::vector<Double_t>& params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const;
  void RadialKernel(const std::vector<Double_t>& params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const;


  void printMatrix(const TMatrixD& matrix) const;
  void printVector(const TVectorD& vector) const;

protected:
  
  // instance variables
  mutable TDecompSVD* _svd;
  mutable TMatrixD*   _resinv;
  mutable Bool_t      _haveMLEst;
  mutable Bool_t      _haveMLCov;
  mutable Bool_t      _MLHConverged;

  mutable TVectorD    _MAPEst;
  mutable TMatrixD    _MAPCov;

  mutable TVectorD    _MLEst;
  mutable TMatrixD    _MLCov;
  mutable TMatrixD    _K;
  mutable TMatrixD    _Kstar;
  mutable TMatrixD    _Kstarstar;

  mutable Double_t    _truMin;
  mutable Double_t    _truMax;
  mutable Double_t    _obsMin;
  mutable Double_t    _obsMax;
  mutable TVectorD    _truBinCenters;
  mutable TVectorD    _obsBinCenters;

  mutable std::vector<Double_t> _kernel_init;
  mutable std::vector<Double_t> _kernel_step;
  mutable std::vector<Double_t> _opt_params;
  

  mutable Int_t _kernel;


public:
  ClassDefT (RooUnfoldGPT, 1)  // Unregularised unfolding
};

typedef RooUnfoldGPT<TH1,TH2> RooUnfoldGP;
#ifndef NOROOFIT
typedef RooUnfoldGPT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldGP;
#endif

#endif /*ROOUNFOLDGP_H_*/
