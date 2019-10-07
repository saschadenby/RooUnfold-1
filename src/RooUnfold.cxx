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

//____________________________________________________________
/*! \class RooUnfoldT
\brief A base class for several unfolding methods.
<p>The unfolding method can either use the constructors for individual unfolding algorithms or the New() method, specifiying the algorithm to be used.
<p>The resultant distribution can be displayed as a plot (Hunfold) or as a bin by bin breakdown of the true, measured and unfolded values (PrintTable)
<p>A covariance matrix can be returned using the Eunfold() method. A vector of its diagonals can be returned with the EunfoldV() method.
<p>A summary of the unfolding algorithms which inherit from this class is below:
<ul>
<li>RooUnfoldBayes: Uses the Bayes method of unfolding based on the method written by D'Agostini (<a href="http://www.slac.stanford.edu/spires/find/hep/www?j=NUIMA,A362,487">NIM A 362 (1995) 487</a>).
<ul>
<li>Works for 1 & 2 dimensional distributions
<li>Returned errors can be either as a diagonal matrix or as a full matrix of covariances
<li>Regularisation parameter sets the number of iterations used in the unfolding (default=4)
<li>Is able to account for bin migration and smearing
<li>Can unfold if test and measured distributions have different binning.
<li>Returns covariance matrices with conditions approximately that of the machine precision. This occasionally leads to very large chi squared values
</ul>
<li> RooUnfoldSVD: Uses the singular value decomposition method of Hocker and Kartvelishvili (<a href="http://arxiv.org/abs/hep-ph/9509307">NIM A 372 (1996) 469</a>)
<ul>
<li>Regularisation parameter defines the level at which values are deemed to be due to statistical fluctuations and are cut out. (Default= number of bins/2)
<li>Returns errors as a full matrix of covariances
<li>Error processing is much the same as with the kCovToy setting with 1000 toys. This is quite slow but can be switched off.
<li>Can only handle 1 dimensional distributions
<li>True and measured distributions must have the same binning
<li>Can account for both smearing and biasing
<li>Returns near singular covariance matrices, again leading to very large chi squared values
</ul>
<li> RooUnfoldIds: Uses the Bayes method of unfolding based on the method written by Malaescu (<a href="http://arxiv.org/abs/1106.3107">CERN-PH-EP-2011-111</a>)
<ul>
<li>Set the number of iterations used to improve the folding matrix
<li>Regularisation parameters define the level at which values are deemed to be due to statistical fluctuations. Used for modifying the folding matrix, as well as unfolding.
<li>Returns errors as a full matrix of covariances
<li>Error processing is much the same as with the kCovToy setting with 1000 toys. This is quite slow but can be switched off.
<li>Can handle 2 dimensional distributions
<li>True and measured distributions must have the same binning
<li>Can account for both smearing and biasing
</ul>
<li> RooUnfoldBinByBin: Unfolds using the method of correction factors.
<ul>
<li>Returns errors as a diagonal matrix.
<li>Is not able to handle bin migration caused by bias/smearing of the distribution
<li>Can only handle 1 dimensional distributions
<li>True and measured distributions must have the same binning
</ul>
<li> RooUnfoldTUnfold: Uses the unfolding method implemented in ROOT's <a href="http://root.cern.ch/root/html/TUnfold.html">TUnfold</a> class
<ul>
<li>Only included in ROOT versions 5.22 and higher
<li>Only able to unfold 1 dimensional distributions
<li>Can account for bin migration and smearing
<li>Errors come as a full covariance matrix.
<li>Will sometimes warn of "unlinked" bins. These are bins with 0 entries and do not effect the results of the unfolding
<li>Regularisation parameter can be either optimised internally by plotting log10(chi2 squared) against log10(tau). The 'kink' in this curve is deemed the optimum tau value. This value can also be set manually (FixTau)
<li>The latest version (TUnfold v15) requires that RooUnfoldResponse::SetOverflow=0. ROOT versions 5.26 or below use v13 and so should be safe to use overflows
</ul>
<li> RooUnfoldInvert: The simplest method of unfolding works by simply inverting the response matrix.
<ul>
<li>For small statistics, this method does not produce useful results.
<li>The inversion method is included largely to illustrate the necessity of a more effective method of unfolding</ul>
</ul>
\class RooUnfold
\brief specialization of RooUnfoldT for TH1/TH2 objects
 */

/////////////////////////////////////////////////////////////

#include "RooUnfold.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>

#include "TClass.h"
#include "TMatrixD.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TDecompChol.h"
#include "TRandom.h"
#include "TMath.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldErrors.h"
// Need subclasses just for RooUnfold::New()
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldInvert.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldGP.h"
#ifndef NOTUNFOLD
#include "RooUnfoldTUnfold.h"
#endif
#ifdef HAVE_DAGOSTINI
#include "RooUnfoldDagostini.h"
#endif
#include "RooUnfoldIds.h"
#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"
#include "RooUnfoldFitHelpers.h"
#include "RooBinning.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;
using std::sqrt;
using std::fabs;

template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kNone = RooUnfolding::kNone;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kBayes = RooUnfolding::kBayes;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kSVD = RooUnfolding::kSVD;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kBinByBin = RooUnfolding::kBinByBin;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kTUnfold = RooUnfolding::kTUnfold;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kInvert = RooUnfolding::kInvert;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kDagostini = RooUnfolding::kDagostini;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kIDS = RooUnfolding::kIDS;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kGP = RooUnfolding::kGP; 
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kNoError = RooUnfolding::kNoError;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kErrors = RooUnfolding::kErrors;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kCovariance = RooUnfolding::kCovariance;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kCovToy = RooUnfolding::kCovToy;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kDefault = RooUnfolding::kDefault;

using namespace RooUnfolding;

template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, const char* name, const char* title)
  : TNamed (name, title)
{

  //! Constructor with response matrix object and measured unfolding input histogram.
  //! Should not normally be used directly - instead, create an instance of one of RooUnfold's subclasses,
  //! or use the New() static constructor.
  Init();
  Setup (res, meas);
}

template<class Hist,class Hist2D> RooUnfoldT<Hist,Hist2D>*
RooUnfoldT<Hist,Hist2D>::New (RooUnfolding::Algorithm alg, const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas,Double_t regparm,
                           const char* name, const char* title)
{

    /*!Unfolds according to the value of the alg enum:
    0 = kNone:     dummy unfolding
    1 = kBayes:    Unfold via iterative application of Bayes theorem
    2 = kSVD:      Unfold using singlar value decomposition (SVD)
    3 = kBinByBin: Unfold bin by bin.
    4 = kTUnfold:  Unfold with TUnfold
    5 = kInvert:   Unfold using inversion of response matrix
    7 = kIDS:      Unfold using iterative dynamically stabilized (IDS) method
    */
  
  RooUnfoldT<Hist,Hist2D>* unfold(NULL);

  switch(alg) {

  case kNone:
    unfold= new RooUnfoldT<Hist,Hist2D>         (res, meas);
    break;

  case kBayes:
    unfold= new RooUnfoldBayesT<Hist,Hist2D>    (res, meas);
    break;

  case kSVD:
    unfold= new RooUnfoldSvdT<Hist,Hist2D>      (res, meas);
    break;

  case kBinByBin:
    unfold= new RooUnfoldBinByBinT<Hist,Hist2D> (res, meas);
    break;

  case kTUnfold:
#ifndef NOTUNFOLD
    unfold= new RooUnfoldTUnfoldT<Hist,Hist2D> (res,meas);
    break;
#else
    cerr << "TUnfold library is not available" << endl;
    return 0;
#endif
    
  case kInvert:
    unfold = new RooUnfoldInvertT<Hist,Hist2D>  (res,meas);
    break;
  case kGP:
    unfold = new RooUnfoldGPT<Hist,Hist2D> (res,meas);
    break;

  case kDagostini:
    cerr << "RooUnfoldDagostini is not available" << endl;
    return 0;
  
  case kIDS:
    unfold= new RooUnfoldIdsT<Hist,Hist2D>      (res, meas,4);
    break;

  default: 
    cerr << "Unknown RooUnfold method " << Int_t(alg) << endl;
    return 0;
  }

  if (name)  unfold->SetName  (name);
  if (title) unfold->SetTitle (title);
  unfold->SetAlgorithm(alg);
  if (regparm != -1e30){
    unfold->SetRegParm(regparm);
  }

  return unfold;
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Destroy()
{
  //! clear the cache
  _cache = Cache();
}

template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::Cache::Cache() :
  _minparm(0),
  _maxparm(0),
  _stepsizeparm(0),
  _defaultparm(0),
  _unfolded(false),
  _fail(false),
  _haveCov(false),
  _haveWgt(false),
  _have_err_mat(false),
  _haveBias(false),
  _haveErrors(false),
  _rec(1),
  _cov(1,1),
  _wgt(1,1),
  _variances(1),
  _err_mat(1,1),
  _bias(1),
  _sigbias(1),
  _vMes(0),
  _eMes(0),
  _covL(0),
  _covMes(0),
  _withError(kDefault)
{
  //! default constructor
}

template<class Hist,class Hist2D> 
typename RooUnfoldT<Hist,Hist2D>::Cache& RooUnfoldT<Hist,Hist2D>::Cache::operator= ( const RooUnfoldT<Hist,Hist2D>::Cache & other ){
  //! assignment operator
  _minparm = other._minparm;
  _maxparm = other._maxparm;
  _stepsizeparm = other._stepsizeparm;
  _defaultparm = other._defaultparm;
  _unfolded = other._unfolded;
  _haveCov = other._haveCov;
  _fail = other._fail;
  _have_err_mat = other._have_err_mat;
  _haveBias =  other._haveBias;
  _haveErrors = other._haveErrors;
  _haveWgt = other._haveWgt;
  _withError = other._withError;
  _rec.ResizeTo(other._rec);
  _cov.ResizeTo(other._cov);
  _wgt.ResizeTo(other._wgt);
  _variances.ResizeTo(other._variances);
  _err_mat.ResizeTo(other._err_mat);
  _rec = other._rec;
  _cov = other._cov;
  _wgt = other._wgt;
  _variances = other._variances;
  _err_mat = other._err_mat;
  _bias = other._bias;
  _sigbias = other._sigbias;
  _vMes = other._vMes;
  _eMes = other._eMes;
  _covL = other._covL;
  _covMes = other._covMes;
  return *this;
}


template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::Cache::~Cache(){
  //! destructor
  delete this->_vMes;
  delete this->_eMes;
  delete this->_covMes;
  delete this->_covL;
}

template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const RooUnfoldT<Hist,Hist2D>& rhs)
  : TNamed (rhs.GetName(), rhs.GetTitle())
{
  //! Copy constructor.
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Assign (const RooUnfoldT<Hist,Hist2D>& rhs)
{
  //! assign data from another unfolding object
  if (this == &rhs) return;
  Reset();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::CopyData (const RooUnfoldT<Hist,Hist2D>& rhs)
{
  //! copy data from another unfolding object
  Setup (new RooUnfoldResponseT<Hist,Hist2D>(*(rhs.response())), clone(rhs.Hmeasured()));
  SetVerbose (rhs.verbose());
  SetNToys   (rhs.NToys());
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetAlgorithm (RooUnfolding::Algorithm alg)
{
  //! set the unfolding algorithm to be used
  _alg = alg;
}

template<class Hist,class Hist2D> RooUnfolding::Algorithm
RooUnfoldT<Hist,Hist2D>::GetAlgorithm () const
{
  //! return the unfolding algorithm used
  return _alg;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Reset()
{
  //! clear and reinitialize
  Destroy();
  Init();
}

namespace { 
  void setBinning(RooRealVar* obs, const TAxis* ax, bool includeUnderflowOverflow){
    int n = ax->GetNbins()+(includeUnderflowOverflow?2:0);
    if(ax->IsVariableBinSize()){
      std::vector<double> bounds;
      for(int i=0; i<ax->GetNbins()+1; ++i){
        bounds.push_back(ax->GetBinLowEdge(i+1));
      }
      RooBinning bins(n,&((bounds[0])));
      obs->setBinning(bins);
    } else {
      obs->setBins(n);
    }
  }
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Init()
{
  //! initialize an object with zero
  _res= 0;
  _meas= 0;
  _nm= _nt= 0;
  _verbose= 1;
  _overflow= 0;
  _dosys= false;
  _covMes= 0;
  _NToys=50;
  GetSettings();
}

template<class Hist,class Hist2D> RooUnfoldT<Hist,Hist2D>&
RooUnfoldT<Hist,Hist2D>::Setup (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas)
{
  //! setup object from a response
  Reset();
  SetResponse (res);
  SetMeasured (meas);
  return *this;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasured (const Hist* meas)
{
  
  //! Set measured distribution and errors. RooUnfold does not own the histogram.
  _meas= clone(meas);
  _cache = Cache();
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasured (const TVectorD& meas, const TVectorD& err)
{
  //! Set measured distribution and errors. Should be called after setting response matrix.
  const Hist* orig = _res->Hmeasured();
  _meas = RooUnfolding::createHist<Hist>(meas,GetName(),GetTitle(),var(orig,X));
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasured (const TVectorD& meas, const TMatrixD& cov)
{
  //! Set measured distribution and its covariance matrix. Should be called after setting response matrix.
  SetMeasuredCov (cov);
  SetMeasured (meas, Emeasured());
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasuredCov (const TMatrixD& cov)
{
  //! Set covariance matrix on measured distribution.
  _cache = Cache();
  _covMes= new TMatrixD (cov);
}

template<class Hist,class Hist2D> const TMatrixD&
RooUnfoldT<Hist,Hist2D>::GetMeasuredCov() const
{
  //! Get covariance matrix on measured distribution.
  if (_covMes) return *_covMes;
  const TVectorD& err(Emeasured());
  _cache._covMes= new TMatrixD (_nm,_nm);
  for (Int_t i= 0 ; i<_nm; i++) {
    Double_t e= err[i];
    (*_cache._covMes)(i,i)= e*e;
  }
  return *_cache._covMes;
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::ForceRecalculation () {
  //! clear and rebuild the cache
  this->_cache = Cache();
  this->_res->ClearCache();
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetResponse (const RooUnfoldResponseT<Hist,Hist2D>* res, Bool_t takeOwnership){
  //! Set response matrix for unfolding, optionally taking ownership of the RooUnfoldResponseT<Hist,Hist2D> object
  if(!res) throw std::runtime_error("cannot set response to invalid value!");
  if(takeOwnership) _res= const_cast<RooUnfoldResponseT<Hist,Hist2D>*>(res);
  else _res = new RooUnfoldResponseT<Hist,Hist2D>(*res);
  _overflow= _res->UseOverflowStatus() ? 1 : 0;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  
  SetNameTitleDefault();
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Unfold() const
{

  //! Dummy unfolding - just copies input
  cout << "********************** " << ClassName() << ": dummy unfolding - just copy input **********************" << endl;

  _cache._rec.ResizeTo (_nt);
  Int_t nb= _nm < _nt ? _nm : _nt;
  TVectorD vmeas(Vmeasured());
  for (Int_t i= 0; i < nb; i++) {
    _cache._rec(i)= vmeas(i);
  }
  _cache._unfolded= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetErrors() const
{
    //!Creates vector of diagonals of covariance matrices.
    //!This may be overridden if it can be computed more quickly without the covariance matrix.
    if (!_cache._haveCov) GetCov();
    if (!_cache._haveCov) return;
    _cache._variances.ResizeTo(_nt);
    for (Int_t i= 0; i < _nt; i++) {
      _cache._variances(i)= _cache._cov(i,i);
    }
    _cache._haveErrors= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetCov() const
{
  //!Dummy routine to get covariance matrix. It should be overridden by derived classes.
  const TMatrixD& covmeas(GetMeasuredCov());
  Int_t nb= std::min(_nm,_nt);
  _cache._cov.ResizeTo (nb, nb);
  for (int i=0; i<nb; i++)
    for (int j=0; j<nb; j++)
      _cache._cov(i,j)= covmeas(i,j);
  _cache._haveCov= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetWgt() const
{
  //! Creates weight matrix
  //! This may be overridden if it can be computed directly without the need for inverting the matrix
  if (!_cache._haveCov) GetCov();
  if (!_cache._haveCov) return;
  if (!InvertMatrix (_cache._cov, _cache._wgt, "covariance matrix", _verbose)) return;
  _cache._haveWgt= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetErrMat() const
{
  //! Get covariance matrix from the variation of the results in toy MC tests
  if (_NToys<=1) return;
  _cache._err_mat.ResizeTo(_nt,_nt);
  TVectorD xisum (_nt);
  TMatrixD xijsum(_nt,_nt);
  for (Int_t k=0; k<_NToys; k++){
    RooUnfoldT<TH1,TH2>* unfold= RunToy();
    const TVectorD& x(unfold->Vunfold());
    for (Int_t i=0; i<_nt;i++){
      Double_t xi= x[i];
      xisum[i] += xi;
      for (Int_t j=0; j<_nt; j++) xijsum(i,j) += xi * x[j];
    }
    delete unfold;
  }
  for (Int_t i=0; i<_nt; i++){
    for (Int_t j=0; j<_nt; j++){
      _cache._err_mat(i,j)= (xijsum(i,j) - (xisum[i]*xisum[j])/_NToys) / (_NToys-1);
    }
  }
  _cache._have_err_mat=true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::CalculateBias(Int_t ntoys, const Hist* hTrue) const
{
  //! TODO: document
  if (ntoys<=1) return;
  _cache._bias.ResizeTo(_nt);
  _cache._sigbias.ResizeTo(_nt);

  TMatrixD pull_results(ntoys,_nt);

  TVectorD truth(_nt);

  if (!hTrue){
    truth = _res->Vtruth();
  } else {
    truth = h2v(hTrue,false);
  }


  for ( int i = 0; i < ntoys; i++){
    RooUnfoldT<TH1,TH2>* unfold = RunToy();
    const TVectorD& toy_meas(unfold->Vmeasured());
    const TVectorD& toy_unfold(unfold->Vunfold());
    const TVectorD& toy_error(unfold->EunfoldV());
    const TVectorD& truth(_res->Vtruth());

    for (int j = 0; j < toy_unfold.GetNrows(); j++){
      if (toy_error(j) != 0){
	pull_results(i, j) = (toy_unfold(j) - truth(j)) / toy_error(j);
	_cache._bias(j) += pull_results(i, j);
      } 
    }  
    delete unfold;
  }
  for (int i = 0; i < _nt; i++){
    _cache._bias(i) = _cache._bias(i) / ntoys;
  }

  for (int j = 0; j < _nt; j++){
    for (int i = 0; i < ntoys; i++){  
      _cache._sigbias(j) += sqrt( (pull_results(i, j) - _cache._bias(j)) * (pull_results(i, j) - _cache._bias(j)) );
    }
    _cache._sigbias(j) = _cache._sigbias(j) / ntoys;
  }

  _cache._haveBias=true;
}

template<class Hist,class Hist2D> Bool_t
RooUnfoldT<Hist,Hist2D>::UnfoldWithErrors (ErrorTreatment withError, bool getWeights) const
{
  //! This method initializes the unfolding.

  if (!_cache._unfolded) {

    if (_cache._fail) return false;

    const Hist* rmeas= _res->Hmeasured();
    if (dim(_meas) != dim(rmeas) ||
        nBins(_meas,X)    != nBins(rmeas,X)    ||
        nBins(_meas,Y)    != nBins(rmeas,Y)    ||
        nBins(_meas,Z)    != nBins(rmeas,Z)) {
      cerr << "Warning: measured "              << nBins(_meas,X);
      if (dim(_meas)>=2) cerr << "x" << nBins(_meas,Y);
      if (dim(_meas)>=3) cerr << "x" << nBins(_meas,Z);
      cerr << "-bin histogram does not match "  << nBins(rmeas,X);
      if (dim(rmeas)>=2) cerr << "x" << nBins(rmeas,Y);
      if (dim(rmeas)>=3) cerr << "x" << nBins(rmeas,Z);
      cerr << "-bin measured histogram from RooUnfoldResponse" << endl;
    }

    this->Unfold();

    if (!_cache._unfolded) {
      _cache._fail= true;
      return false;
    }
  }

  Bool_t ok;
  _cache._withError= withError;
  if (getWeights && (withError==kErrors || withError==kCovariance)) {
      if   (!_cache._haveWgt)      GetWgt();
      ok= _cache._haveWgt;
  } else {
    switch (withError) {
    case kErrors:
      if   (!_cache._haveErrors)   GetErrors();
      ok= _cache._haveErrors;
      break;
    case kCovariance:
      if   (!_cache._haveCov)      GetCov();
      ok= _cache._haveCov;
      break;
    case kCovToy:
      if   (!_cache._have_err_mat) GetErrMat();
      ok= _cache._have_err_mat;
      break;
    default:
      ok= true;
    }
  }

  if (!ok) _cache._fail= true;
  
  return ok;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::Chi2(const Hist* hTrue,ErrorTreatment DoChi2) const {
    /*!Calculates Chi squared. Method depends on value of DoChi2
    0: sum of (residuals/error)squared
    1: use errors propagated through the unfolding
    2: use covariance matrix returned from unfolding
    3: use covariance matrix from the variation of the results in toy MC tests
    Returns warnings for small determinants of covariance matrices and if the condition is very large.
    If a matrix has to be inverted also removes rows/cols with all their elements equal to 0*/

    if (!UnfoldWithErrors (DoChi2)) return -1.0;

    TVectorD res(subtract<Hist,TVectorD>(_cache._rec,hTrue,_overflow));

    Double_t chi2= 0.0;
    if (DoChi2==kCovariance || DoChi2==kCovToy) {
      TMatrixD wgt(Wunfold(DoChi2));
      if (_cache._fail) return -1.0;
      TMatrixD resmat(1,_nt), chi2mat(1,1);
      TMatrixDRow(resmat,0)= res;
      ABAT (resmat, wgt, chi2mat);
      chi2= chi2mat(0,0);
    } else {
      TVectorD eunfold(EunfoldV(DoChi2));
      if (_cache._fail) return -1.0;
      for (Int_t i = 0 ; i < _nt; i++) {
        Double_t e= eunfold[i];
        if (e<=0.0) continue;
        Double_t ypull = res[i] / e;
        chi2 += ypull*ypull;
      }
    }
    return chi2;
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::PrintTable (std::ostream& o, const Hist* hTrue, ErrorTreatment withError) const
{

  //! Prints entries from truth, measured, and unfolded data for each bin.
  if (withError==kDefault) withError= _cache._withError;
  if (withError==kDefault) withError= kErrors;

  if (!UnfoldWithErrors (withError)) withError= kNoError;

  if (!hTrue){
    hTrue = response()->Htruth();
  }

  const Hist* hTrainTrue = response()->Htruth();
  const Hist* hTrain = response()->Hmeasured();
  const Hist* hMeas = Hmeasured();

  int ntxb= nBins(_res->Htruth(),X)+2*this->_overflow;
  int ntyb= nBins(_res->Htruth(),Y)+2*this->_overflow;

  int d = dim(_res->Htruth());
  if (!_cache._unfolded) return;
  Double_t chi_squ= -999.0;
  if (hTrue && (withError==kCovariance || withError==kCovToy)) chi_squ = Chi2(hTrue,withError);

  printTable(o,d,
             ntxb,ntyb,
             h2v(hTrainTrue,this->_overflow, this->response()->UseDensityStatus()),
             h2v(hTrain,this->_overflow, this->response()->UseDensityStatus()),
             hTrue ? h2v(hTrue,this->_overflow, this->response()->UseDensityStatus()) : TVectorD(this->_nt) ,             
             h2v(hMeas,this->_overflow),
             this->Vunfold(),
             withError,
             hTrue ? h2ve(hTrue,this->_overflow, this->response()->UseDensityStatus()) : TVectorD(this->_nt) ,
             this->EunfoldV(withError),
             chi_squ);
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetNameTitleDefault()
{
  if (!_res) return;
  const char* s= GetName();
  if (s[0] == '\0') SetName (_res->GetName());
  s= GetTitle();
  if (s[0] == '\0') {
    TString title= "Unfold ";
    title += _res->GetTitle();
    SetTitle (title);
  }
}

template<class Hist,class Hist2D> Hist*
RooUnfoldT<Hist,Hist2D>::Hunfold (ErrorTreatment withError)
{
    /*!Creates unfolded distribution. Error calculation varies by withError:
    0: No errors
    1: Errors from the square root of the diagonals of the covariance matrix given by the unfolding
    2: Errors from the square root of of the covariance matrix given by the unfolding
    3: Errors from the square root of the covariance matrix from the variation of the results in toy MC tests
    */

  if (!UnfoldWithErrors (withError)) withError= kNoError;
  const Hist* t = _res->Htruth();
  if (!_cache._unfolded){
    return RooUnfolding::createHist(name(t),title(t),vars(t));
  } else {
    TVectorD rec(this->Vunfold());
    TVectorD errors(this->EunfoldV());
    return RooUnfolding::createHist(rec,errors,name(t),title(t),vars(t),_overflow);
  }
}

template<class Hist,class Hist2D> TH1*
RooUnfoldT<Hist,Hist2D>::TH1unfold() {
  //! TODO: document
  Double_t min = RooUnfolding::min(this->response()->Htruth(),RooUnfolding::X);
  Double_t max = RooUnfolding::max(this->response()->Htruth(),RooUnfolding::X);
  Int_t nbins;
  Int_t i_st;

  const TVectorD& unfolded(this->Vunfold());
  const TVectorD& err(this->EunfoldV());
  if (this->_overflow){
    nbins = unfolded.GetNrows() - 2;
    i_st = 0;
    min = min + (max - min)/(nbins + 2);
    max = max - (max - min)/(nbins + 2);
  } else {
    i_st = 1;
    nbins = unfolded.GetNrows();
  }

  TH1D* unf_th1 = new TH1D("test","test",nbins,min,max);
  for (int i = i_st; i < unfolded.GetNrows(); i++){
    unf_th1->SetBinContent(i, unfolded[i - i_st]);
    unf_th1->SetBinError(i, err[i - i_st]);
  }

  return unf_th1;
}

template<class Hist, class Hist2D> TH1*
RooUnfoldT<Hist,Hist2D>::TH1bias() {
  //! TODO: document
  Double_t min = RooUnfolding::min(this->response()->Hmeasured(),RooUnfolding::X);
  Double_t max = RooUnfolding::max(this->response()->Hmeasured(),RooUnfolding::X);

  TH1D* bias = new TH1D("bias","bias",_cache._bias.GetNrows(),min,max);

  for (int i = 0; i < bias->GetNbinsX(); i++){
    bias->SetBinContent(i + 1, _cache._bias(i));
    bias->SetBinError(i + 1, _cache._sigbias(i));
  }

  return bias;
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetSettings() const
{
    //!Gets maximum and minimum parameters and step size
    _cache._minparm=0;
    _cache._maxparm=0;
    _cache._stepsizeparm=0;
    _cache._defaultparm=0;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetMinParm() const
{
    //!Get minimum regularisation parameter for unfolding method
    return _cache._minparm;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetMaxParm() const
{
    //!Get maximum regularisation parameter for unfolding method
    return _cache._maxparm;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetStepSizeParm() const
{
    //!Get suggested step size for unfolding distribution
    return _cache._stepsizeparm;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetDefaultParm() const
{
    //!Get suggested regularisation parameter.
    return _cache._defaultparm;
}

template<class Hist,class Hist2D> RooUnfoldT<TH1,TH2>*
RooUnfoldT<Hist,Hist2D>::RunToy() const
{
  //! Returns new RooUnfold object with smeared measurements and
  //! (if IncludeSystematics) response matrix for use as a toy.
  //! Use multiple toys to find spread of unfolding results.
  TString name= GetName();
  name += "_toy";

  const TVectorD& measv(Vmeasured());
  const TVectorD& errv(Emeasured());
  const TVectorD& truv(_res->Vtruth());
  const TVectorD& recov(_res->Vmeasured());
  TMatrixD respm(_res->Mresponse(false));
  
  Double_t min_reco = RooUnfolding::min(this->response()->Hmeasured(),RooUnfolding::X);
  Double_t max_reco = RooUnfolding::max(this->response()->Hmeasured(),RooUnfolding::X);
  Double_t min_truth = RooUnfolding::min(this->response()->Htruth(),RooUnfolding::X);
  Double_t max_truth = RooUnfolding::max(this->response()->Htruth(),RooUnfolding::X);

  
  TH1D *newmeas = new TH1D("newhist","newhist",measv.GetNrows(),min_reco,max_reco);
  TH1D *tru = new TH1D("tru","tru",truv.GetNrows(),min_truth,max_truth);
  TH1D *reco = new TH1D("reco","reco",recov.GetNrows(),min_reco,max_reco);
  TH2D *newresp = new TH2D("newresp","newresp",respm.GetNrows(),min_reco,max_reco,truv.GetNrows(),min_truth,max_truth);
  
  // Smear the response matrix?
			    
  // Change the creation of toy measured by folding truth hist
  // and add gaussian/poisson noise to that?
  
  for (Int_t i= 0; i<_nm; i++) {
    
    reco->SetBinContent(i + 1, recov(i));

    Double_t e= errv(i);
    gRandom->SetSeed(0);
    Double_t g= gRandom->Gaus(0,e);
    if (e>0.0) newmeas->SetBinContent(i+1, measv(i) + g);
    
    for (Int_t j = 0; j < _nt; j++){
      tru->SetBinContent(j + 1, truv(j));
      newresp->SetBinContent(i + 1, j + 1, respm(i,j));
    }
  }

  RooUnfoldResponseT<TH1,TH2>* res= new RooUnfoldResponseT<TH1,TH2> (_res->GetName(),_res->GetTitle(), newresp, tru, reco, false, false);
  
  RooUnfolding::Algorithm alg = this->GetAlgorithm();
  Double_t regparm = GetRegParm();

  RooUnfoldT<TH1,TH2>* unfold = RooUnfoldT<TH1,TH2>::New(alg,res,newmeas,regparm);

  delete newmeas;
  delete tru;
  delete reco;
  delete newresp;

  return unfold;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Print(Option_t* /*opt*/) const
{
  //! print a summary of the configuration
  cout << ClassName() << "::" << GetName() << " \"" << GetTitle()
       << "\"," << " regularisation parameter=" << GetRegParm() << ", ";
  if (_covMes) cout << "with measurement covariance, ";
  if (_dosys)      cout << "calculate systematic errors, ";
  if (dim(_meas)==1) cout << _nm;
  else {
    cout <<        nBins(_meas,X)
         << "x" << nBins(_meas,Y);
    if (dim(_meas)>=3)
      cout << "x" << nBins(_meas,Z);
    cout << " (" << _nm << ")";
  }
  cout << " bins measured, ";
  const Hist* rtrue= _res->Htruth();
  if (dim(rtrue)==1) cout << _nt;
  else {
    cout <<        nBins(rtrue,X)
         << "x" << nBins(rtrue,Y);
    if (dim(rtrue)>=3)
      cout << "x" << nBins(rtrue,Z);
    cout << " (" << _nt << ")";
  }
  cout << " bins truth";
  if (_overflow) cout << " including overflows";
  cout << endl;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Dump() const {
  //! dump the contents of the member variables
  std::cout << "covMes=" <<  _covMes << std::endl;
  std::cout << "verbose=" <<  _verbose << std::endl;
  std::cout << "nm=" <<  _nm << std::endl;
  std::cout << "nt=" <<  _nt << std::endl;
  std::cout << "overflow=" <<  _overflow << std::endl;
  std::cout << "NToys=" <<  _NToys << std::endl;
  std::cout << "dosys=" <<  _dosys << std::endl;
  std::cout << "res=" <<  _res << std::endl;
  std::cout << "meas=" <<  _meas << std::endl;
  _res->Print();
  _meas->Print();
}

template<class Hist,class Hist2D> TMatrixD
RooUnfoldT<Hist,Hist2D>::CutZeros(const TMatrixD& ereco)
{
    //!Removes row & column if all their elements are 0.
    vector<int> diags;
        int missed=0;
        for (int i=0; i<ereco.GetNrows(); i++){
            double coltot=0;
            for (int j=0;j<ereco.GetNrows();j++){
                coltot+=ereco(i,j);
            }
            if (coltot==0){
                diags.push_back(i);
                missed++;
            }
        }
        int x=ereco.GetNrows()-missed;
        int y=ereco.GetNcols()-missed;
        TMatrixD ereco_cut(x,y);
        unsigned int v=0;
        for (int i=0;i<ereco.GetNrows();i++){
            if(v<diags.size() && diags[v]==i){
                v++;
            }
            else{
                for (int j=0; j<ereco_cut.GetNcols();j++){
                    ereco_cut(i-v,j)=ereco(i,j+v);
                    }
                }
        }
    return ereco_cut;
}

template<class Hist,class Hist2D> TMatrixD
RooUnfoldT<Hist,Hist2D>::Eunfold(ErrorTreatment withError) const
{
    /*!Returns covariance matrices for error calculation of type withError
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests
    */
    if (!UnfoldWithErrors (withError)) return TMatrixD(_nt,_nt);

    switch(withError){
    case kNoError: {
      TMatrixD Eunfold_m(_nt,_nt);
      for (int i=0; i<_nt; i++){
        Eunfold_m(i,i)=_cache._rec(i);
      }
      return Eunfold_m;
      break; }
    case kErrors: {
      TMatrixD Eunfold_m(_nt,_nt);
      for (int i=0; i<_nt;i++){
        Eunfold_m(i,i)=_cache._variances(i);
      }
      return Eunfold_m;
      break;
    }
    case kCovariance:
      return _cache._cov;
      break;
    case kCovToy:
      return _cache._err_mat;
      break;
    default:
      throw std::runtime_error("Error, unrecognised error method");
    }
}

template<class Hist,class Hist2D> TVectorD
RooUnfoldT<Hist,Hist2D>::EunfoldV(ErrorTreatment withError) const
{
    /*!Returns vector of unfolding errors computed according to the withError flag:
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests
    */

    TVectorD Eunfold_v(_nt);
    if (!UnfoldWithErrors (withError)) return Eunfold_v;

    switch(withError){
      case kNoError:
        for (int i=0; i<_nt; i++){
          Eunfold_v(i)=sqrt (fabs (_cache._rec(i)));
        }
        break;
      case kErrors:
        for (int i=0; i<_nt; i++){
          Eunfold_v(i)=sqrt (fabs (_cache._variances(i)));
        }
        break;
      case kCovariance:
        for (int i=0; i<_nt; i++){
          Eunfold_v(i)=sqrt (fabs (_cache._cov(i,i)));
        }
        break;
      case kCovToy:
        for (int i=0; i<_nt; i++){
          Eunfold_v(i)=sqrt (fabs (_cache._err_mat(i,i)));
        }
        break;
      default:
        throw std::runtime_error("Error, unrecognised error method");
    }
    return Eunfold_v;
}

template<class Hist,class Hist2D> TMatrixD
RooUnfoldT<Hist,Hist2D>::Wunfold(ErrorTreatment withError) const
{
  //! TODO: document
    TMatrixD Wunfold_m(_nt,_nt);
    if (!UnfoldWithErrors (withError, true)) return Wunfold_m;

    switch(withError){
      case kNoError:
        for (int i=0; i<_nt; i++){
          if (_cache._rec(i)!=0.0) Wunfold_m(i,i)=1.0/_cache._rec(i);
        }
        break;
      case kErrors:
        for (int i=0; i<_nt;i++){
          Wunfold_m(i,i)=_cache._wgt(i,i);
        }
        break;
      case kCovariance:
        Wunfold_m=_cache._wgt;
        break;
      case kCovToy:
        InvertMatrix (_cache._err_mat, Wunfold_m, "covariance matrix from toys", _verbose);
        break;
      default:
        cerr<<"Error, unrecognised error method= "<<withError<<endl;
    }
    return Wunfold_m;
}

template<class Hist,class Hist2D> Int_t
RooUnfoldT<Hist,Hist2D>::InvertMatrix(const TMatrixD& mat, TMatrixD& inv, const char* name, Int_t verbose)
{
  //! Invert a matrix using Single Value Decomposition: inv = mat^-1.
  //! Can use InvertMatrix(mat,mat) to invert in-place.
  Int_t ok= 1;
  TDecompSVD svd (mat);
  const Double_t cond_max= 1e17;
  Double_t cond= svd.Condition();
  if (verbose >= 1) {
    Double_t d1=0,d2=0;
    svd.Det(d1,d2);
    Double_t det= d1*TMath::Power(2.,d2);
    cout << name << " condition="<<cond<<", determinant="<<det;
    if (d2!=0.0) cout <<" ("<<d1<<"*2^"<<d2<<")";
    cout <<", tolerance="<<svd.GetTol()<<endl;
  }
  if        (cond<0.0){
    cerr <<"Warning: bad "<<name<<" condition ("<<cond<<")"<<endl;
    ok= 2;
  } else if (cond>cond_max) {
    cerr << "Warning: poorly conditioned "<<name<<" - inverse may be inaccurate (condition="<<cond<<")"<<endl;
    ok= 3;
  }
  inv.ResizeTo (mat.GetNcols(), mat.GetNrows());  // pseudo-inverse of A(r,c) is B(c,r)
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,13,4)  /* TDecompSVD::Invert() didn't have ok status before 5.13/04. */
  Bool_t okinv= false;
  inv= svd.Invert(okinv);
  if (!okinv) {
    cerr << name << " inversion failed" << endl;
    return 0;
  }
#else
  inv= svd.Invert();
#endif
  if (verbose>=1) {
    TMatrixD I (mat, TMatrixD::kMult, inv);
    if (verbose>=3) printMatrix(I,"V*V^-1");
    Double_t m= 0.0;
    for (Int_t i= 0; i<I.GetNrows(); i++) {
      Double_t d= fabs(I(i,i)-1.0);
      if (d>m) m= d;
      for (Int_t j= 0; j<i; j++) {
        d= fabs(I(i,j)); if (d>m) m= d;
        d= fabs(I(j,i)); if (d>m) m= d;
      }
    }
    cout << "Inverse " << name << " " << 100.0*m << "% maximum error" << endl;
  }
  return ok;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Streamer (TBuffer &R__b)
{
  //! Stream an object of class RooUnfold.
  if (R__b.IsReading()) {
    RooUnfoldT<Hist,Hist2D>::Class()->ReadBuffer  (R__b, this);
  } else {
    RooUnfoldT<Hist,Hist2D>::Class()->WriteBuffer (R__b, this);
  }
}

template<> void
RooUnfoldT<TH1,TH2>::Streamer (TBuffer &R__b)
{
  //! Stream an object of class RooUnfold.
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    RooUnfoldT<TH1,TH2>::Class()->ReadBuffer  (R__b, this);
    TH1::AddDirectory (oldstat);
  } else {
    RooUnfoldT<TH1,TH2>::Class()->WriteBuffer (R__b, this);
  }
}

template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::RooUnfoldT()
  : TNamed()
{
  //! Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const char*    name, const char*    title)
  : TNamed(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const TString& name, const TString& title)
  : TNamed(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>::~RooUnfoldT()
{
  // destructor
  Destroy();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>& RooUnfoldT<Hist,Hist2D>::operator= (const RooUnfoldT<Hist,Hist2D>& rhs)
{
  //! Assignment operator for copying RooUnfold settings.
  Assign(rhs);
  return *this;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::verbose() const
{
  //! Get verbosity setting which controls amount of information to be printed
  return _verbose;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::NToys()     const
{
  //! Get number of toys used in kCovToy error calculation.
  return _NToys;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::Overflow()  const
{
  //! Histogram under/overflow bins are used?
  return _overflow;
}

template<class Hist,class Hist2D> 
const RooUnfoldResponseT<Hist,Hist2D>* RooUnfoldT<Hist,Hist2D>::response()  const
{
   //! Response matrix object
  return _res;
}

template<class Hist,class Hist2D> 
RooUnfoldResponseT<Hist,Hist2D>* RooUnfoldT<Hist,Hist2D>::response()
{  
   //! Response matrix object
  return _res;
}

template<class Hist,class Hist2D> 
const Hist*               RooUnfoldT<Hist,Hist2D>::Hmeasured() const
{
  //! Measured Distribution as a histogram
  return _meas;
}

template<class Hist,class Hist2D> 
Hist*               RooUnfoldT<Hist,Hist2D>::Hmeasured()
{
  //! Measured Distribution as a histogram
  return _meas;
}

template<class Hist,class Hist2D> 
const TVectorD&                RooUnfoldT<Hist,Hist2D>::Vunfold() const
{
  //! Unfolded (reconstructed) distribution as a vector
  if (!_cache._unfolded) {
    if (!_cache._fail){
      this->Unfold();
    }
    if (!_cache._unfolded) {
      _cache._fail= true;
      if (_nt > 0 && _cache._rec.GetNrows() == 0) _cache._rec.ResizeTo(_nt);   // need something
    }
  }

  return _cache._rec;
}

template<class Hist,class Hist2D> 
const TVectorD&          RooUnfoldT<Hist,Hist2D>::Vmeasured() const
{
  //! Measured distribution as a vector.
  if (!_cache._vMes){
    _cache._vMes = new TVectorD(h2v (_meas, _overflow, this->response()->UseDensityStatus()));
  }
  return *_cache._vMes;
}

template<class Hist,class Hist2D> 
const TVectorD&          RooUnfoldT<Hist,Hist2D>::Emeasured() const
{
  //! Measured errors as a vector.
  if (!_cache._eMes){
    if(_covMes){
      _cache._eMes= new TVectorD(_nm);
      for (Int_t i= 0; i<_nm; i++) {
        Double_t e= this->_cache._cov(i,i);
        if (e>0.0) (*_cache._eMes)[i]= sqrt(e);
      }
    } else {
      _cache._eMes = new TVectorD(h2ve (_meas, _overflow, this->response()->UseDensityStatus()));
    }
  }
  return *_cache._eMes;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetVerbose (Int_t level)
{
  //! Set verbosity level which controls amount of information to be printed
  _verbose= level;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetOverflow (Int_t overflow)
{
  //! set the usage of the overflow bin
  _overflow= overflow;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetNToys (Int_t toys)
{
  //! Set number of toys used in kCovToy error calculation.
  _NToys= toys;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetRegParm (Double_t regparm)
{
  //! Set Regularisation parameter
}

template<class Hist,class Hist2D> 
Double_t RooUnfoldT<Hist,Hist2D>::GetRegParm() const
{
  //! Get regularisation parameter.
  return -1;
}

template<class Hist,class Hist2D> 
void RooUnfoldT<Hist,Hist2D>::IncludeSystematics (Int_t dosys)
{
  //! Include systematic errors from response matrix?
  //! Use dosys=2 to exclude measurement errors.
  if (dosys!=_dosys){
    _cache = Cache();
    _dosys= dosys;
  }
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::SystematicsIncluded() const
{
  //! return setting for whether to include systematic errors from response matrix
  return _dosys;
}

template class RooUnfoldT<TH1,TH2>;
ClassImp (RooUnfold)

#ifndef NOROOFIT
#include "RooHistFunc.h"
#include "RooRealSumFunc.h"
#include "RooPrintable.h"
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooProduct.h"


template<> void RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::SetResponse (const RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* res, Bool_t takeOwnership){
  //! Set response matrix for unfolding, optionally taking ownership of the RooUnfoldResponseT<Hist,Hist2D> object
  if(!res) throw std::runtime_error("cannot set response to invalid value!");
  _res = new RooFitUnfoldResponse(res);
  _overflow= _res->UseOverflowStatus() ? 1 : 0;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  SetNameTitleDefault();
}

template class RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
typedef RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooUnfoldT_RooFitHist;
ClassImp (RooUnfoldT_RooFitHist)

template<class Base>RooUnfolding::RooFitWrapper<Base>::RooFitWrapper(const char* name, const char* title, const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unf) : Base(name,title){
  auto unfolding = RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::New(unf->GetMethod(),unf->response(),unf->Hmeasured(),unf->GetRegParm(),unf->GetName(),unf->GetTitle());
  this->_unfolding = dynamic_cast<RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>*>(unfolding);
  this->_unfolding->SetVerbose(0);
  const RooUnfoldResponseT<RooFitHist,RooFitHist>* res = this->_unfolding->response();
  if(res){
    const RooFitHist* htruth = res->Htruth();
    if(htruth){
      this->addServer(*(htruth->func()));
      for(int i=0; i<dim(htruth); ++i){
        this->addServer(*htruth->obs(i));
      }
    }
    const RooFitHist* hfakes = res->Hfakes();
    if(hfakes){
      this->addServer(*(hfakes->func()));
      for(int i=0; i<dim(hfakes); ++i){
        this->addServer(*hfakes->obs(i));
      }
    }
    const RooFitHist* hresponse = res->Hresponse();
    if(hresponse){
      this->addServer(*(hresponse->func()));
      for(int i=0; i<dim(hresponse); ++i){
        this->addServer(*hresponse->obs(i));
      }
    }
    const RooFitHist* hmeasured = res->Hmeasured();
    if(hmeasured){
      this->addServer(*(hmeasured->func()));    
      for(int i=0; i<dim(hmeasured); ++i){
        this->addServer(*hmeasured->obs(i));
      }
    }
  }
  const RooFitHist* hmeasured = this->_unfolding->Hmeasured();

  if(hmeasured){
    this->addServer(*(hmeasured->func()));    
    for(int i=0; i<dim(hmeasured); ++i){
      this->addServer(*hmeasured->obs(i));
    }
  }
}
RooUnfoldFunc::RooUnfoldFunc(const char* name, const char* title, const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unf) : RooFitWrapper(name,title,unf) {}
RooUnfoldPdf::RooUnfoldPdf(const char* name, const char* title, const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unf) : RooFitWrapper(name,title,unf) {}
template<class Base>RooUnfolding::RooFitWrapper<Base>::RooFitWrapper() : _unfolding(NULL) {}
template<class Base>RooUnfolding::RooFitWrapper<Base>::~RooFitWrapper(){
  
  delete _unfolding;
}
template<class Base> Bool_t RooUnfolding::RooFitWrapper<Base>::redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive){
  RooUnfoldResponseT<RooFitHist,RooFitHist>* res = this->_unfolding->response();
  if(res){
    RooFitHist* htruth = res->Htruth();
    if(htruth){
      htruth->replace(newServerList);
    }
    RooFitHist* hfakes = res->Hfakes();
    if(hfakes){
      hfakes->replace(newServerList);
    }
    RooFitHist* hresponse = res->Hresponse();
    if(hresponse){
      hresponse->replace(newServerList);
    }
    RooFitHist* hmeasured = res->Hmeasured();
    if(hmeasured){
      hmeasured->replace(newServerList);
    }
  }
  RooFitHist* hmeasured = this->_unfolding->Hmeasured();
  if(hmeasured){
    hmeasured->replace(newServerList);
  }
  return Base::redirectServersHook(newServerList,mustReplaceAll,nameChange,isRecursive);
}
 




template<class Base>const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* RooUnfolding::RooFitWrapper<Base>::unfolding() const { 
  return this->_unfolding;
}

template <class Base>
std::list<Double_t>* RooUnfolding::RooFitWrapper<Base>::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
  //! retrieve the list of bin boundaries
  return this->_unfolding->response()->Htruth()->func()->binBoundaries(obs,xlo,xhi);
}

template <class Base>
std::list<Double_t>* RooUnfolding::RooFitWrapper<Base>::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
  //! retrieve the sampling hint
  return this->_unfolding->response()->Htruth()->func()->plotSamplingHint(obs,xlo,xhi);
}

template <class Base>
Double_t RooUnfolding::RooFitWrapper<Base>::getValV(const RooArgSet* set) const
{
  this->_curNormSet = set ;
  return Base::getValV(set) ;
}

template <class Base>
Double_t RooUnfolding::RooFitWrapper<Base>::evaluate() const {
  //! call getVal on the internal function
  std::map<std::string,double> snapshot;
  this->_unfolding->response()->Hresponse()->saveSnapshot(snapshot);
  int bin = this->_unfolding->response()->Htruth()->bin();
  this->_unfolding->ForceRecalculation();
  this->_unfolding->response()->Htruth()->checkValidity();
  double v = std::max(this->_unfolding->Vunfold()[bin],this->_minVal);
  if(this->_unfolding->response()->UseDensityStatus()){
    v /= binVolume(this->_unfolding->response()->Htruth(),bin,false);
  }
  this->_unfolding->response()->Hresponse()->loadSnapshot(snapshot);
  return v;  
}

template <class Base>
Bool_t  RooUnfolding::RooFitWrapper<Base>::isBinnedDistribution(const RooArgSet& obs) const {
  //! check if this PDF is a binned distribution in the given observable
  return this->_unfolding->response()->Hresponse()->func()->isBinnedDistribution(obs);
}

template<class Base>
Bool_t RooUnfolding::RooFitWrapper<Base>::checkObservables(const RooArgSet *nset) const {
  return this->_unfolding->response()->Hresponse()->func()->checkObservables(nset);
}

template<class Base>
Bool_t RooUnfolding::RooFitWrapper<Base>::forceAnalyticalInt(const RooAbsArg &arg) const {
  return this->_unfolding->response()->Htruth()->func()->forceAnalyticalInt(arg);
}

template<class Base>
Int_t RooUnfolding::RooFitWrapper<Base>::getAnalyticalIntegralWN(RooArgSet &allVars, RooArgSet &numVars, const RooArgSet *normSet, const char *rangeName) const {
  return this->_unfolding->response()->Htruth()->func()->getAnalyticalIntegralWN(allVars,numVars,normSet,rangeName);
}

template<class Base>
Double_t RooUnfolding::RooFitWrapper<Base>::analyticalIntegralWN(Int_t code, const RooArgSet *normSet, const char *rangeName) const {
  double val = 0;
  auto vec = this->_unfolding->Vunfold();
  for(int i=0; i<vec.GetNrows(); ++i){
    // assuming that density correction has been applied already
    val += vec[i];
  }
  return val;
}

template<class Base>
void RooUnfolding::RooFitWrapper<Base>::printMetaArgs(std::ostream &os) const {
  return this->_unfolding->response()->Htruth()->func()->printMetaArgs(os);
}

template<class Base>
RooAbsArg::CacheMode RooUnfolding::RooFitWrapper<Base>::canNodeBeCached() const {
  return this->_unfolding->response()->Htruth()->func()->canNodeBeCached();
}
template<class Base>
void RooUnfolding::RooFitWrapper<Base>::setCacheAndTrackHints(RooArgSet& arg) {
  this->_unfolding->response()->Htruth()->func()->setCacheAndTrackHints(arg);
}
template<class Base> TObject* RooUnfolding::RooFitWrapper<Base>::clone(const char* newname) const {
  return new RooUnfolding::RooFitWrapper<Base>(newname ? newname : this->GetName(),this->GetTitle(),this->_unfolding);
}



template class RooUnfolding::RooFitWrapper<RooAbsReal>;
template class RooUnfolding::RooFitWrapper<RooAbsPdf>;

RooUnfoldFunc::RooUnfoldFunc() : RooFitWrapper() {
}
RooUnfoldFunc::~RooUnfoldFunc() {
}
TObject* RooUnfoldFunc::clone(const char* newname) const {
  return new RooUnfoldFunc(newname ? newname : this->GetName(),this->GetTitle(),this->_unfolding);
}
ClassImp (RooUnfoldFunc)
RooUnfoldPdf::RooUnfoldPdf() : RooFitWrapper() {
}
RooUnfoldPdf::~RooUnfoldPdf( ) {}
TObject* RooUnfoldPdf::clone(const char* newname) const {

  RooUnfoldPdf* retval = new RooUnfoldPdf(newname ? newname : this->GetName(),this->GetTitle(),this->_unfolding);
  
  return retval;
}
ClassImp (RooUnfoldPdf)

RooAbsPdf::ExtendMode RooUnfoldPdf::extendMode() const {
    //! Return extended mode capabilities
  return RooAbsPdf::MustBeExtended;
}
Double_t RooUnfoldPdf::expectedEvents(const RooArgSet* nset) const {
  //! Return expected number of events for extended likelihood calculation
  //! which is the sum of all coefficients
  std::map<std::string,double> snapshot;
  this->_unfolding->response()->Hresponse()->saveSnapshot(snapshot);
  this->_unfolding->ForceRecalculation();
  this->_unfolding->response()->Htruth()->checkValidity();
  double events = 0;
  const auto vec(this->_unfolding->Vunfold());
  for(int bin = 0; bin<vec.GetNrows(); ++bin){
    events += vec[bin];
  }
  this->_unfolding->response()->Hresponse()->loadSnapshot(snapshot);
  return events;
}
Bool_t RooUnfoldPdf::selfNormalized() const {
  //! P.d.f is self normalized
  return kTRUE ;
}

namespace {
  bool readToken(TString& instr, std::vector<TString>& tokens){
    int pos = instr.First(",");
    if(pos == -1){
      tokens.push_back(instr);
      return false;
    } else {
      tokens.push_back(instr(0,pos));
      instr.Remove(0,pos+1);
      return true;
    }
  }
}



#include <RooAbsCollection.h>

#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
namespace {
  
  struct IteratorHelper {
    RooFIter itr;
    RooAbsArg* nextItem;
    IteratorHelper(const RooAbsCollection& c);
    IteratorHelper();
    RooAbsArg* operator++();
    bool operator!=(const IteratorHelper& other);
    bool operator!=(const RooAbsArg* other);  
    RooAbsArg* operator*();
  };
  
  IteratorHelper::IteratorHelper(const RooAbsCollection& c) : itr(c.fwdIterator()), nextItem(itr.next()) {}
  IteratorHelper::IteratorHelper() : itr(), nextItem(NULL) {}  
  RooAbsArg* IteratorHelper::operator++(){
    nextItem = itr.next();
    return nextItem;
  }
  bool IteratorHelper::operator!=(const IteratorHelper& other){
    return this->nextItem != other.nextItem;
  }
  bool IteratorHelper::operator!=(const RooAbsArg* other){
    return this->nextItem != other;
  }
  
  RooAbsArg* IteratorHelper::operator*(){
    return nextItem;
  }
}

::IteratorHelper begin(const RooAbsCollection& c){
  return ::IteratorHelper(c);
}

::IteratorHelper end(const RooAbsCollection&){
  return ::IteratorHelper();
}
#endif

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth, const char* obs_truth, const TH1* reco, const char* obs_reco, const TH2* response, const TH1* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  RooUnfoldSpec(name,title,truth,obs_truth,reco,obs_reco,response,0,data,includeUnderflowOverflow,errorThreshold,useDensity)
{
  //! constructor forwarding
}


RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth, const char* obs_truth, const TH1* reco, const char* obs_reco, const TH2* response, const TH1* bkg, const TH1* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  int d =dim(truth);
  if(d!=dim(reco)){
    throw std::runtime_error("inconsistent dimensionality between truth and reco histograms!");
  }
  
  TString obs_truth_s(obs_truth);
  TString obs_reco_s(obs_reco);
  std::vector<TString> obs_truth_v;
  std::vector<TString> obs_reco_v;
  bool more_truth = false;
  bool more_reco = false;
  for(int i=0; i<d; ++i){
    more_truth = ::readToken(obs_truth_s,obs_truth_v);
    more_reco = ::readToken(obs_reco_s,obs_reco_v);
    if(!more_truth || !more_reco) break;
  }
  if(more_truth) throw std::runtime_error(TString::Format("encountered additional characters on truth observable list: '%s'",obs_truth_s.Data()).Data());
  if(more_reco) throw std::runtime_error(TString::Format("encountered additional characters on reco observable list: '%s'",obs_reco_s.Data()).Data());
  if((int)obs_truth_v.size() != d) throw std::runtime_error(TString::Format("truth observable list is too short for %d dimensions: '%s'",d,obs_truth).Data());
  if((int)obs_reco_v.size() != d) throw std::runtime_error(TString::Format("reco observable list is too short for %d dimensions: '%s'",d,obs_truth).Data());
    
  RooArgList truth_vars;
  for(int i=0; i<d; ++i){
    const TAxis* ax = RooUnfolding::getAxis(truth,(RooUnfolding::Dimension)i);
    double min = ax->GetBinLowEdge(!includeUnderflowOverflow);
    double max = ax->GetBinUpEdge(ax->GetNbins()+includeUnderflowOverflow);
    RooRealVar* obs = new RooRealVar(obs_truth_v[i],ax->GetTitle() ? ax->GetTitle() : obs_truth_v[i].Data(),min,min,max);
    setBinning(obs,ax,includeUnderflowOverflow);
    obs->setConstant(true);
    truth_vars.add(*obs);
  }
  RooArgList reco_vars;
  for(int i=0; i<d; ++i){
    const TAxis* ax = RooUnfolding::getAxis(reco,(RooUnfolding::Dimension)i);
    double min = ax->GetBinLowEdge(!includeUnderflowOverflow);
    double max = ax->GetBinUpEdge(ax->GetNbins()+includeUnderflowOverflow);
    RooRealVar* obs = new RooRealVar(obs_reco_v[i],ax->GetTitle() ? ax->GetTitle() : obs_reco_v[i].Data(),min,min,max);
    setBinning(obs,ax,includeUnderflowOverflow);
    obs->setConstant(true);
    reco_vars.add(*obs);
  }

  this->setup(truth,truth_vars,reco,reco_vars,response,bkg,data,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, const RooArgList& obs_truth, const TH1* reco_th1, const RooArgList& obs_reco, const TH2* response_th1, const TH1* bkg, const TH1* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  this->setup(truth_th1,obs_truth,reco_th1,obs_reco,response_th1,bkg,data,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, const RooArgList& obs_truth, const TH1* reco_th1, const RooArgList& obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) :
  TNamed(name,title)
{
  this->_bkg.setNominal(bkg);
  this->_data.setNominal(RooUnfolding::makeHistFunc(data,obs_reco));
  this->setup(truth_th1,obs_truth,reco_th1,obs_reco,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, const RooArgList& obs_truth, RooAbsReal* reco, const RooArgList& obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) :
  TNamed(name,title)
{
  this->_reco.setNominal(reco);
  this->_bkg.setNominal(bkg);
  this->_data.setNominal(RooUnfolding::makeHistFunc(data,obs_reco));
  this->setup(truth_th1,obs_truth,NULL,obs_reco,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, RooAbsReal* reco, RooAbsArg* obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : RooUnfoldSpec(name,title,truth_th1,RooArgList(*obs_truth),reco,RooArgList(*obs_reco),response_th1,bkg,data,includeUnderflowOverflow,errorThreshold,useDensity) {}

namespace {
  RooAbsReal* makeParamHistFunc(const char* name, const char* title, const RooArgList& observables, const RooArgList& parameters, bool multiplyDensity){
    ParamHistFunc* phf = new ParamHistFunc(name,title,observables,parameters);
    if(!multiplyDensity) return phf;
    RooArgList inverseBinWidths;
    RooRealVar* obs = (RooRealVar*)(observables.first());
    for(int i=0; i<obs->getBinning().numBins(); ++i){
      double width = obs->getBinWidth(i);
      RooRealVar* bw = new RooRealVar(TString::Format("%s_bin_%d_widthCorr",obs->GetName(),i),TString::Format("density correction of %s in bin %d",obs->GetTitle(),i),1./width);
      bw->setConstant(true);
      inverseBinWidths.add(*bw);
    }
    ParamHistFunc* bwcorr = new ParamHistFunc(TString::Format("%s_widthCorr",obs->GetName()),TString::Format("density corrections for %s",obs->GetTitle()),observables,inverseBinWidths);
    RooArgList elems;
    elems.add(*phf);
    elems.add(*bwcorr);
    RooProduct* prod = new RooProduct(TString::Format("%s_x_%s",phf->GetName(),bwcorr->GetName()),title,elems);
    return prod;
  }
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const RooArgList& reco_bins, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgList& bkg_bins, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_reco.setNominal(::makeParamHistFunc(TString::Format("signal_reco_%s_differential",obs_reco->GetName()).Data(),obs_reco->GetTitle(),obs_reco_list,reco_bins,useDensity));
  this->_bkg.setNominal(::makeParamHistFunc(TString::Format("bkg_reco_%s_differential",obs_reco->GetName()).Data(),obs_reco->GetTitle(),obs_reco_list,bkg_bins,useDensity));
  this->_data.setNominal(RooUnfolding::makeHistFunc(data,obs_reco_list));
  this->setup(truth_th1,obs_truth_list,NULL,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco_th1, RooAbsArg* obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : RooUnfoldSpec(name,title,truth_th1,RooArgList(*obs_truth),reco_th1,RooArgList(*obs_reco),response_th1,bkg,data,includeUnderflowOverflow,errorThreshold,useDensity) {}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco_th1, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgList& bkg_bins, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);  
  this->_bkg.setNominal(::makeParamHistFunc(TString::Format("bkg_reco_%s_differential",obs_reco->GetName()),obs_reco->GetTitle(),obs_reco_list,bkg_bins,useDensity));
  this->_data.setNominal(RooUnfolding::makeHistFunc(data,obs_reco_list));
  this->setup(truth_th1,obs_truth_list,reco_th1,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco_th1, RooAbsArg* obs_reco, const TH2* response_th1, RooAbsReal* measured, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : TNamed(name,title)
{
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_data.setNominal(measured);
  this->setup(truth_th1,obs_truth_list,reco_th1,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco_th1, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgList& measured_bins, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_data.setNominal(::makeParamHistFunc(TString::Format("measured_reco_%s_differential",obs_reco->GetName()),obs_reco->GetTitle(),obs_reco_list,measured_bins,useDensity));
  this->setup(truth_th1,obs_truth_list,reco_th1,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

namespace {
  RooRealSumFunc* makeRooRealSumFunc(const char* name, const char* title, const RooArgSet& contributions, double densityCorrection){
    RooRealVar* binWidth = new RooRealVar(TString::Format("%s_binWidthCorrection",name),"bin width correction",densityCorrection);
    binWidth->setConstant(true);
    RooArgList functions;
    RooArgList coefs;
    RooFIter itr(contributions.fwdIterator());
    RooAbsArg* obj;
    while((obj = itr.next())){
      functions.add(*obj);
      coefs.add(*binWidth);
    }
    return new RooRealSumFunc(name,title,functions,coefs);
  }
}


RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, RooAbsReal* reco, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgSet& bkg_contributions, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) :
  TNamed(name,title)
{
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_reco.setNominal(reco);
  this->_data.setNominal(RooUnfolding::makeHistFunc(data,obs_reco_list));
  double densityCorr = 1;
  if(useDensity && obs_reco->InheritsFrom(RooRealVar::Class())){
    densityCorr = 1./((RooRealVar*)(obs_reco))->getBinning().averageBinWidth();
  }
  this->_bkg.setNominal(::makeRooRealSumFunc(TString::Format("bkg_reco_%s_differential",obs_reco->GetName()),obs_reco->GetTitle(),bkg_contributions,densityCorr));
  this->setup(truth_th1,obs_truth_list,NULL,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, RooAbsReal* truth, RooAbsArg* obs_truth,  RooAbsReal* reco, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgSet& bkg_contributions, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) :
  TNamed(name,title)
{
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_truth.setNominal(truth);
  this->_reco.setNominal(reco);
  this->_data.setNominal(RooUnfolding::makeHistFunc(data,obs_reco_list));
  double densityCorr = 1;
  if(useDensity && obs_reco->InheritsFrom(RooRealVar::Class())){
    densityCorr = 1./((RooRealVar*)(obs_reco))->getBinning().averageBinWidth();
  }
  this->_bkg.setNominal(::makeRooRealSumFunc(TString::Format("bkg_reco_%s_differential",obs_reco->GetName()),obs_reco->GetTitle(),bkg_contributions,densityCorr));
  this->setup(NULL,obs_truth_list,NULL,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}


void RooUnfoldSpec::setup(const TH1* truth_th1, const RooArgList& obs_truth, const TH1* reco_th1, const RooArgList& obs_reco, const TH2* response_th1, const TH1* bkg_th1, const TH1* data_th1, bool includeUnderflowOverflow, double errorThreshold, bool useDensity){
  this->_includeUnderflowOverflow = includeUnderflowOverflow;
  this->_useDensity = useDensity;
  this->_errorThreshold = errorThreshold;
  if(truth_th1) this->_truth.setNominal(RooUnfolding::makeHistFunc(truth_th1,obs_truth,includeUnderflowOverflow,this->_useDensity));
  if(reco_th1)  this->_reco.setNominal(RooUnfolding::makeHistFunc(reco_th1,obs_reco,includeUnderflowOverflow,this->_useDensity));
  this->_obs_truth.add(obs_truth);  
  this->_obs_all.add(obs_truth);
  this->_obs_reco.add(obs_reco);
  this->_obs_all.add(obs_reco);  
  if(response_th1) this->_res.setNominal(RooUnfolding::makeHistFunc(response_th1,this->_obs_all,includeUnderflowOverflow,this->_useDensity));
  if(bkg_th1)  this->_bkg.setNominal(RooUnfolding::makeHistFunc(bkg_th1,obs_reco,includeUnderflowOverflow,this->_useDensity));
  
  if(data_th1) this->_data.setNominal(RooUnfolding::makeHistFunc(data_th1,obs_reco,includeUnderflowOverflow,this->_useDensity));
}

RooUnfoldSpec::~RooUnfoldSpec(){
}


RooUnfolding::RooFitHist* RooUnfoldSpec::makeHistogram(const TH1* hist){
    
  return new RooUnfolding::RooFitHist(hist, this->_obs_truth, this->_includeUnderflowOverflow, this->_errorThreshold, this->_useDensity);
}


RooUnfolding::RooFitHist* RooUnfoldSpec::makeHistogram(const HistContainer& source, double errorThreshold){
  RooAbsReal* hf = source._nom;
  RooAbsReal* func = hf;
  if(!hf) return 0;
  std::vector<RooAbsArg*> obs;
  RooArgList obslist;
  if(this->_obs_all.getSize() < 1){
    throw std::runtime_error("in RooUnfoldSpec::makeHistogram: no observables known!");
  }
  RooFIter itr(this->_obs_all.fwdIterator());
  RooAbsArg* arg = NULL;
  while((arg = itr.next())){
    if(!arg) continue;
    if(!hf->dependsOn(*arg)) continue;
    obs.push_back(arg);
    obslist.add(*arg);
  }
  if(obslist.getSize() < 1){
    std::stringstream ss;
    ss << "in RooUnfoldSpec::makeHistogram: function '";
    hf->printStream(ss,0,RooPrintable::kStandard,"");
    ss << "' does not depend on any of the known observables '";
    this->_obs_all.printStream(ss,0,RooPrintable::kStandard,"");
    ss << "'!";
    throw std::runtime_error(ss.str());
  }
  if(source._shapes.size() > 0){
    RooArgList up, dn;
    RooArgList params;
    for(auto var:source._shapes){
      TString sysname(var.first);
      if(var.second.size() != 2){
        throw std::runtime_error(TString::Format("unable to process systematics '%s' with size %d != 2",var.first.c_str(),(int)(var.second.size())).Data());
      }
      up.add(*var.second[0]);
      dn.add(*var.second[1]);
      TString name = TString::Format("alpha_%s",var.first.c_str());
      RooRealVar* p = (RooRealVar*)(this->_alphas.find(name));
      if(!p){
        p = new RooRealVar(name,name,0,-5,5);
        p->setError(1);
        this->addGaussNP(p);
      }
      params.add(*p);
    }
    TString name = TString::Format("%s_HistoSystematics",hf->GetName());
    hf->SetName(TString::Format("%s_Nominal",hf->GetName()));
    func = new PiecewiseInterpolation(name.Data(),name.Data(),*hf,up,dn,params);
  }
  RooArgList components; 
  if(source._norms.size() > 0){
    std::vector<double> up,dn;
    RooArgList params;
    for(auto var:source._norms){
      TString sysname(var.first);
      up.push_back(var.second.first);
      dn.push_back(var.second.second);
      TString name = TString::Format("alpha_%s",var.first.c_str());
      RooRealVar* p = (RooRealVar*)(this->_alphas.find(name));
      if(!p){
        p = new RooRealVar(name,name,0,-5,5);
        p->setError(1);
        this->addGaussNP(p);
      }
      params.add(*p);
    }
    TString name = TString::Format("%s_OverallSystematics",hf->GetName());
    components.add(*(new RooStats::HistFactory::FlexibleInterpVar(name.Data(),name.Data(),params,1.,up,dn)));
  }  
  std::vector<RooRealVar*> gammas;
  if(errorThreshold >= 0 && hf->InheritsFrom(RooHistFunc::Class())){
    gammas = RooUnfolding::createGammas(&(((RooHistFunc*)hf)->dataHist()),obslist,errorThreshold);
    for(auto g:gammas){
      this->addPoissonNP(g);
    }
    RooAbsReal* phf = RooUnfolding::makeParamHistFunc(hf->GetName(),hf->GetTitle(),obslist,gammas);
    if(phf) components.add(*phf);
  }
  if(components.getSize() > 0){
    components.add(*func);
    TString name(hf->GetName());
    hf->SetName(TString::Format("%s_hist",hf->GetName()));
    func = new RooProduct(name.Data(),hf->GetTitle(),components);  
  }
  return new RooUnfolding::RooFitHist(func,obs,gammas);
}

RooHistFunc* RooUnfoldSpec::makeHistFuncTruth(const TH1* hist){
  if(!hist) return NULL;
  return RooUnfolding::makeHistFunc(hist,this->_obs_truth, this->_includeUnderflowOverflow, this->_useDensity);
}

RooHistFunc* RooUnfoldSpec::makeHistFuncMeasured(const TH1* hist){
  if(!hist) return NULL;
  return RooUnfolding::makeHistFunc(hist,this->_obs_reco, this->_includeUnderflowOverflow, this->_useDensity);
}

RooProdPdf* RooUnfoldSpec::makeConstraints(){
  RooArgList constraints;
  for(auto a:this->_alphas){
    RooRealVar* p = (RooRealVar*)(a);
    RooRealVar* mean = new RooRealVar(TString::Format("%s_nom",p->GetName()),TString::Format("Mean value for %s",p->GetName()),p->getVal());
    RooRealVar* sigma = new RooRealVar(TString::Format("%s_sigma",p->GetName()),TString::Format("Sigma for %s",p->GetName()),p->getError());
    mean->setConstant(true);
    sigma->setConstant(true);
    RooGaussian* gaus = new RooGaussian(TString::Format("%s_constraint",p->GetName()),TString::Format("Gaussian constraint term for %s",p->GetName()),*p,*mean,*sigma);
    constraints.add(*gaus);
  }
  for(auto g:this->_gammas){
    RooRealVar* p = (RooRealVar*)(g);
    double n = pow(p->getError(),-2);
    RooRealVar* tau = new RooRealVar(TString::Format("%s_tau",p->GetName()),TString::Format("Tau parameter value for %s",p->GetName()),n);
    tau->setConstant(true);
    RooArgList params(*p,*tau);
    RooProduct* prod = new RooProduct(TString::Format("%s_nEvents",p->GetName()),TString::Format("Number of events for %s",p->GetName()),params);
    RooPoisson* pois = new RooPoisson(TString::Format("%s_constraint",p->GetName()),TString::Format("Poisson constraint term for %s",p->GetName()),*prod,*tau);
    constraints.add(*pois);
  }
  return new RooProdPdf(TString::Format("%s_constraints",this->GetName()),"Unfolding constraint terms",constraints);
}

void RooUnfoldSpec::addGaussNP(RooRealVar* v){
  if(v) this->_alphas.add(*v);
}
void RooUnfoldSpec::addPoissonNP(RooRealVar* v){
  if(v) this->_gammas.add(*v);
}

void RooUnfoldSpec::makeBackground(){
  this->_locked = true;
  if(!this->_cache._bkg){
    this->_cache._bkg = this->makeHistogram(this->_bkg,this->_errorThreshold);
  }
}
void RooUnfoldSpec::makeData(){
  this->_locked = true;
  if(!this->_cache._data){
    this->_cache._data = this->makeHistogram(this->_data,0);
  }
}
void RooUnfoldSpec::makeResponse(){
  this->_locked = true;
  if(!this->_cache._res){
    this->makeReco();
    this->makeTruth();
    if(!this->_res._nom) throw std::runtime_error("no response input given!");        
    this->_cache._res = this->makeHistogram(this->_res,this->_errorThreshold);
    this->_cache._response = new RooFitUnfoldResponse(this->GetName(),this->GetTitle(),this->_cache._res,this->_cache._truth,this->_cache._reco,this->_useDensity);
  }
}
void RooUnfoldSpec::makeTruth(){
  this->_locked = true;
  if(!this->_cache._truth){
    if(!this->_truth._nom) throw std::runtime_error("no truth input given!");    
    this->_cache._truth = this->makeHistogram(this->_truth,this->_errorThreshold);
  }
}
void RooUnfoldSpec::makeReco(){
  this->_locked = true;
  if(!this->_cache._reco){
    if(!this->_reco._nom) throw std::runtime_error("no measured input given!");
    this->_cache._reco = this->makeHistogram(this->_reco,this->_errorThreshold);
  }
}
void RooUnfoldSpec::makeDataMinusBackground(){
  this->_locked = true;
  this->makeData();
  if(!this->_cache._data_minus_bkg){
    if(this->_bkg._nom){
      this->makeResponse();
      this->makeBackground();
      this->_cache._data_minus_bkg = this->_cache._response->makeHistSum(this->_cache._data->func(),this->_cache._bkg->func(),1.,-1.);
    } else {
      this->_cache._data_minus_bkg = this->_cache._data;
    }
    TString name(TString::Format("%s_data_minus_bkg",this->GetName()));
    this->_cache._data_minus_bkg->func()->SetName(name);
    this->_cache._data_minus_bkg->func()->SetTitle(name);
  }
}

RooAbsReal* RooUnfoldSpec::getBackground(){ this->makeBackground(); if(this->_cache._bkg) return this->_cache._bkg->func(); else return NULL; }
RooAbsReal* RooUnfoldSpec::getData(){ this->makeData(); return this->_cache._data->func(); }
RooAbsReal* RooUnfoldSpec::getResponse(){ this->makeResponse(); return this->_cache._res->func(); }
RooAbsReal* RooUnfoldSpec::getTruth(){ this->makeTruth(); return this->_cache._truth->func(); }
RooAbsReal* RooUnfoldSpec::getReco(){ this->makeReco(); return this->_cache._reco->func(); }
RooAbsReal* RooUnfoldSpec::getDataMinusBackground(){ this->makeDataMinusBackground(); return this->_cache._data_minus_bkg->func(); }

RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* RooUnfoldSpec::unfold(Algorithm alg, Double_t regparam){
  this->makeResponse();
  this->makeDataMinusBackground();

  RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unfolding = RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::New(alg,this->_cache._response,this->_cache._data_minus_bkg,regparam);

  unfolding->SetOverflow(this->_includeUnderflowOverflow);

  return unfolding;
}


void RooUnfoldSpec::HistContainer::setNominal(RooAbsReal* nom){
  this->_nom = nom;
}

void RooUnfoldSpec::HistContainer::addShape(const char* name, RooAbsReal* up, RooAbsReal* dn){
  this->_shapes[name] = {up,dn};
}

void RooUnfoldSpec::HistContainer::addNorm(const char* name, double up, double dn){
  this->_norms[name] = {up,dn};
}

RooUnfoldSpec::HistContainer::~HistContainer(){}

void RooUnfoldSpec::lockCheck(){
  if(this->_locked){
    throw std::runtime_error("this instance of RooUnfoldSpec is locked - it has already been used to produce results and can no longer be modified. please create a new instance for modifications!");
  }
}

void RooUnfoldSpec::registerSystematic(Contribution c, const char* name, const TH1* up, const TH1* down){
  this->lockCheck();
  switch(c){
  case kTruth:
    this->_truth.addShape(name,
                          RooUnfolding::makeHistFunc(TString::Format("truth_%s_%s_up",up->GetName(),name),up,this->_obs_truth,this->_includeUnderflowOverflow,this->_useDensity),
                          RooUnfolding::makeHistFunc(TString::Format("truth_%s_%s_dn",down->GetName(),name),down,this->_obs_truth,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  case kMeasured:
    this->_reco.addShape(name,
                         RooUnfolding::makeHistFunc(TString::Format("meas_%s_%s_up",up->GetName(),name),up,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity),
                         RooUnfolding::makeHistFunc(TString::Format("meas_%s_%s_dn",down->GetName(),name),down,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  case kData:
    this->_data.addShape(name,
                         RooUnfolding::makeHistFunc(TString::Format("data_%s_%s_up",up->GetName(),name),up,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity),
                         RooUnfolding::makeHistFunc(TString::Format("data_%s_%s_dn",down->GetName(),name),down,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  case kResponse:
    this->_res.addShape(name,
                        RooUnfolding::makeHistFunc(TString::Format("resp_%s_%s_up",up->GetName(),name),up,this->_obs_all,this->_includeUnderflowOverflow,this->_useDensity),
                        RooUnfolding::makeHistFunc(TString::Format("resp_%s_%s_dn",down->GetName(),name),down,this->_obs_all,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  case kBackground:
    this->_bkg.addShape(name,
                        RooUnfolding::makeHistFunc(TString::Format("bkg_%s_%s_up",up->GetName(),name),up,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity),
                        RooUnfolding::makeHistFunc(TString::Format("bkg_%s_%s_dn",down->GetName(),name),down,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  }
}

void RooUnfoldSpec::registerSystematic(Contribution c, const char* name, double up, double down){
  this->lockCheck();
  switch(c){
  case kTruth:
    this->_truth.addNorm(name,up,down);
    break;
  case kMeasured:
    this->_reco.addNorm(name,up,down);
    break;
  case kData:
    this->_data.addNorm(name,up,down);
    break;
  case kResponse:
    this->_res.addNorm(name,up,down);
    break;
  case kBackground:
    this->_bkg.addNorm(name,up,down);
    break;
  }
}




RooAbsPdf* RooUnfoldSpec::makePdf(Algorithm alg, Double_t regparam){

  RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unfold = this->unfold(alg, regparam);

  RooUnfoldPdf* pdf = new RooUnfoldPdf(this->GetName(),this->GetTitle(),unfold);
  RooProdPdf* constraints = this->makeConstraints();
  RooArgList comps(*pdf,*constraints);
  RooProdPdf* prod = new RooProdPdf(TString::Format("%s_x_constraints",this->GetName()),"Unfolding pdf, including constraints",comps);
  
  delete unfold;
  return prod;
}


RooAbsReal* RooUnfoldSpec::makeFunc(Algorithm alg, Double_t regparam){
  RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unfold = this->unfold(alg, regparam);

  RooUnfoldFunc* func = new RooUnfoldFunc(this->GetName(),this->GetTitle(),this->unfold(alg, regparam));

  delete unfold;
  return func;
}


ClassImp(RooUnfoldSpec)

#endif


