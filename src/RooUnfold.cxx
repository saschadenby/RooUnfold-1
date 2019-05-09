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
/* BEGIN_HTML
<p>A base class for several unfolding methods.
<p>The unfolding method can either use the constructors for individual unfolding algorithms or the New() method, specifiying the algorithm to be used.
<p>The resultant distribution can be displayed as a plot (Hreco) or as a bin by bin breakdown of the true, measured and reconstructed values (PrintTable)
<p>A covariance matrix can be returned using the Ereco() method. A vector of its diagonals can be returned with the ErecoV() method.
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
<li>Only able to reconstruct 1 dimensional distributions
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
END_HTML */

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
  // Constructor with response matrix object and measured unfolding input histogram.
  // Should not normally be used directly - instead, create an instance of one of RooUnfold's subclasses,
  // or use the New() static constructor.
  Init();
  Setup (res, meas);
}

template<class Hist,class Hist2D> RooUnfoldT<Hist,Hist2D>*
RooUnfoldT<Hist,Hist2D>::New (RooUnfolding::Algorithm alg, const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas,Double_t regparm,
                           const char* name, const char* title)
{
    /*Unfolds according to the value of the alg enum:
    0 = kNone:     dummy unfolding
    1 = kBayes:    Unfold via iterative application of Bayes theorem
    2 = kSVD:      Unfold using singlar value decomposition (SVD)
    3 = kBinByBin: Unfold bin by bin.
    4 = kTUnfold:  Unfold with TUnfold
    5 = kInvert:   Unfold using inversion of response matrix
    7 = kIDS:      Unfold using iterative dynamically stabilized (IDS) method
    */
  RooUnfoldT<Hist,Hist2D>* unfold;
  switch (alg) {
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
//    case kTUnfold:
//#ifndef NOTUNFOLD
//      unfold= new RooUnfoldTUnfold  (res,meas);
//      break;
//#else
//      cerr << "TUnfold library is not available" << endl;
//      return 0;
//#endif
//    case kInvert:
//      unfold = new RooUnfoldInvert  (res,meas);
//      break;
//    case kDagostini:
//#ifdef HAVE_DAGOSTINI
//      unfold = new RooUnfoldDagostini (res,meas);
//      break;
//#else
//      cerr << "RooUnfoldDagostini is not available" << endl;
//      return 0;
//#endif
//    case kIDS:
//      unfold= new RooUnfoldIds      (res, meas);
//      break;
//    default:
      cerr << "Unknown RooUnfold method " << Int_t(alg) << endl;
      return 0;
  }
  if (name)  unfold->SetName  (name);
  if (title) unfold->SetTitle (title);
  if (regparm != -1e30){
    unfold->SetRegParm(regparm);
  }
  return unfold;
}


template<> RooUnfoldT<TH1,TH2>*
RooUnfoldT<TH1,TH2>::New (RooUnfolding::Algorithm alg, const RooUnfoldResponseT<TH1,TH2>* res, const TH1* meas,Double_t regparm,
                           const char* name, const char* title)
{
    /*Unfolds according to the value of the alg enum:
    0 = kNone:     dummy unfolding
    1 = kBayes:    Unfold via iterative application of Bayes theorem
    2 = kSVD:      Unfold using singlar value decomposition (SVD)
    3 = kBinByBin: Unfold bin by bin.
    4 = kTUnfold:  Unfold with TUnfold
    5 = kInvert:   Unfold using inversion of response matrix
    7 = kIDS:      Unfold using iterative dynamically stabilized (IDS) method
    */
  RooUnfold* unfold;
  switch (alg) {
    case kNone:
      unfold= new RooUnfold         (res, meas);
      break;
    case kBayes:
      unfold= new RooUnfoldBayes    (res, meas);
      break;
    case kSVD:
      unfold= new RooUnfoldSvd      (res, meas);
      break;
    case kBinByBin:
      unfold= new RooUnfoldBinByBin (res, meas);
      break;
    case kTUnfold:
#ifndef NOTUNFOLD
      unfold= new RooUnfoldTUnfold  (res,meas);
      break;
#else
      cerr << "TUnfold library is not available" << endl;
      return 0;
#endif
    case kInvert:
      unfold = new RooUnfoldInvert  (res,meas);
      break;
    case kDagostini:
#ifdef HAVE_DAGOSTINI
      unfold = new RooUnfoldDagostini (res,meas);
      break;
#else
      cerr << "RooUnfoldDagostini is not available" << endl;
      return 0;
#endif
    case kIDS:
      unfold= new RooUnfoldIds      (res, meas);
      break;
    default:
      cerr << "Unknown RooUnfold method " << Int_t(alg) << endl;
      return 0;
  }
  if (name)  unfold->SetName  (name);
  if (title) unfold->SetTitle (title);
  if (regparm != -1e30){
    unfold->SetRegParm(regparm);
  }
  return unfold;
}


template<class Hist,class Hist2D> RooUnfoldT<Hist,Hist2D>*
RooUnfoldT<Hist,Hist2D>::Clone (const char* newname) const
{
  // Creates a copy of the unfold object
  RooUnfoldT<Hist,Hist2D>* unfold= new RooUnfoldT<Hist,Hist2D>(*this);
  if (newname) unfold->SetName(newname);
  return unfold;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Destroy()
{
  delete _measmine;
  delete _vMes;
  delete _eMes;
  delete _covMes;
  delete _covL;
  delete _resmine;
}

template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const RooUnfoldT<Hist,Hist2D>& rhs)
  : TNamed (rhs.GetName(), rhs.GetTitle())
{
  // Copy constructor.
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Assign (const RooUnfoldT<Hist,Hist2D>& rhs)
{
  if (this == &rhs) return;
  Reset();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::CopyData (const RooUnfoldT<Hist,Hist2D>& rhs)
{
  Setup (rhs.response(), rhs.Hmeasured());
  SetVerbose (rhs.verbose());
  SetNToys   (rhs.NToys());
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Reset()
{
  Destroy();
  Init();
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Init()
{
  _res= _resmine= 0;
  _vMes= _eMes= 0;
  _covMes= _covL= 0;
  _meas= _measmine= 0;
  _nm= _nt= 0;
  _verbose= 1;
  _overflow= 0;
  _dosys= _unfolded= _haveCov= _haveCovMes= _fail= _have_err_mat= _haveErrors= _haveWgt= false;
  _withError= kDefault;
  _NToys=50;
  GetSettings();
}

template<class Hist,class Hist2D> RooUnfoldT<Hist,Hist2D>&
RooUnfoldT<Hist,Hist2D>::Setup (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas)
{
  Reset();
  SetResponse (res);
  SetMeasured (meas);
  return *this;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasured (const Hist* meas)
{
  // Set measured distribution and errors. RooUnfold does not own the histogram.
  _meas= meas;
  delete _vMes; _vMes= 0;
  delete _eMes; _eMes= 0;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasured (const TVectorD& meas, const TVectorD& err)
{
  // Set measured distribution and errors. Should be called after setting response matrix.
  
  const Hist* orig = _res->Hmeasured();
  if (_measmine) {
    delete _measmine;
  }
  _measmine = RooUnfolding::createHist<Hist>(meas,GetName(),GetTitle(),var(orig,X));
  SetMeasured (_measmine);
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasured (const TVectorD& meas, const TMatrixD& cov)
{
  // Set measured distribution and its covariance matrix. Should be called after setting response matrix.
  SetMeasuredCov (cov);
  SetMeasured (meas, Emeasured());
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasuredCov (const TMatrixD& cov)
{
  // Set covariance matrix on measured distribution.
  delete _covL; _covL= 0;
  delete _eMes;
  delete _covMes;
  _eMes= new TVectorD(_nm);
  for (Int_t i= 0; i<_nm; i++) {
    Double_t e= cov(i,i);
    if (e>0.0) (*_eMes)[i]= sqrt(e);
  }
  _covMes= new TMatrixD (cov);
  _haveCovMes= true;
}

template<class Hist,class Hist2D> const TMatrixD&
RooUnfoldT<Hist,Hist2D>::GetMeasuredCov() const
{
  // Get covariance matrix on measured distribution.
  if (_covMes) return *_covMes;
  const TVectorD& err(Emeasured());
  _covMes= new TMatrixD (_nm,_nm);
  for (Int_t i= 0 ; i<_nm; i++) {
    Double_t e= err[i];
    (*_covMes)(i,i)= e*e;
  }
  return *_covMes;
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetResponse (const RooUnfoldResponseT<Hist,Hist2D>* res)
{
  // Set response matrix for unfolding.
  delete _resmine; _resmine= 0;
  _res= res;
  _overflow= _res->UseOverflowStatus() ? 1 : 0;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  if (_overflow) {
    _nm += 2;
    _nt += 2;
  }
  SetNameTitleDefault();
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetResponse (RooUnfoldResponseT<Hist,Hist2D>* res, Bool_t takeOwnership)
{
  // Set response matrix for unfolding, optionally taking ownership of the RooUnfoldResponseT<Hist,Hist2D> object
  SetResponse (res);
  if (takeOwnership) _resmine= res;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Unfold()
{
  // Dummy unfolding - just copies input
  cout << "********************** " << ClassName() << ": dummy unfolding - just copy input **********************" << endl;
  _rec.ResizeTo (_nt);
  Int_t nb= _nm < _nt ? _nm : _nt;
  for (Int_t i= 0; i < nb; i++) {
    _rec(i)= Vmeasured()(i);
  }
  _unfolded= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetErrors()
{
    //Creates vector of diagonals of covariance matrices.
    //This may be overridden if it can be computed more quickly without the covariance matrix.
    if (!_haveCov) GetCov();
    if (!_haveCov) return;
    _variances.ResizeTo(_nt);
    for (Int_t i= 0; i < _nt; i++) {
      _variances(i)= _cov(i,i);
    }
    _haveErrors= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetCov()
{
  //Dummy routine to get covariance matrix. It should be overridden by derived classes.
  const TMatrixD& covmeas(GetMeasuredCov());
  Int_t nb= _nm < _nt ? _nm : _nt;
  _cov.ResizeTo (_nt, _nt);
  for (int i=0; i<nb; i++)
    for (int j=0; j<nb; j++)
      _cov(i,j)= covmeas(i,j);
  _haveCov= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetWgt()
{
  // Creates weight matrix
  // This may be overridden if it can be computed directly without the need for inverting the matrix
  if (!_haveCov) GetCov();
  if (!_haveCov) return;
  if (!InvertMatrix (_cov, _wgt, "covariance matrix", _verbose)) return;
  _haveWgt= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetErrMat()
{
  // Get covariance matrix from the variation of the results in toy MC tests
  if (_NToys<=1) return;
  _err_mat.ResizeTo(_nt,_nt);
  TVectorD xisum (_nt);
  TMatrixD xijsum(_nt,_nt);
  for (Int_t k=0; k<_NToys; k++){
    RooUnfoldT<Hist,Hist2D>* unfold= RunToy();
    const TVectorD& x= unfold->Vreco();
    for (Int_t i=0; i<_nt;i++){
      Double_t xi= x[i];
      xisum[i] += xi;
      for (Int_t j=0; j<_nt; j++) xijsum(i,j) += xi * x[j];
    }
    delete unfold;
  }
  for (Int_t i=0; i<_nt; i++){
    for (Int_t j=0; j<_nt; j++){
      _err_mat(i,j)= (xijsum(i,j) - (xisum[i]*xisum[j])/_NToys) / (_NToys-1);
    }
  }
  _have_err_mat=true;
}

template<class Hist,class Hist2D> Bool_t
RooUnfoldT<Hist,Hist2D>::UnfoldWithErrors (ErrorTreatment withError, bool getWeights)
{
  if (!_unfolded) {
    if (_fail) return false;
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
    if (!_unfolded) {
      _fail= true;
      return false;
    }
  }
  Bool_t ok;
  _withError= withError;
  if (getWeights && (withError==kErrors || withError==kCovariance)) {
      if   (!_haveWgt)      GetWgt();
      ok= _haveWgt;
  } else {
    switch (withError) {
    case kErrors:
      if   (!_haveErrors)   GetErrors();
      ok= _haveErrors;
      break;
    case kCovariance:
      if   (!_haveCov)      GetCov();
      ok= _haveCov;
      break;
    case kCovToy:
      if   (!_have_err_mat) GetErrMat();
      ok= _have_err_mat;
      break;
    default:
      ok= true;
    }
  }
  if (!ok) _fail= true;
  return ok;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::Chi2(const Hist* hTrue,ErrorTreatment DoChi2)
{
    /*Calculates Chi squared. Method depends on value of DoChi2
    0: sum of (residuals/error)squared
    1: use errors propagated through the unfolding
    2: use covariance matrix returned from unfolding
    3: use covariance matrix from the variation of the results in toy MC tests
    Returns warnings for small determinants of covariance matrices and if the condition is very large.
    If a matrix has to be inverted also removes rows/cols with all their elements equal to 0*/

    if (!UnfoldWithErrors (DoChi2)) return -1.0;

    TVectorD res = subtract<Hist,TVectorD>(_rec,hTrue,_overflow);

    Double_t chi2= 0.0;
    if (DoChi2==kCovariance || DoChi2==kCovToy) {
        TMatrixD wgt= Wreco(DoChi2);
        if (_fail) return -1.0;
        TMatrixD resmat(1,_nt), chi2mat(1,1);
        TMatrixDRow(resmat,0)= res;
        ABAT (resmat, wgt, chi2mat);
        chi2= chi2mat(0,0);
    } else {
        TVectorD ereco= ErecoV(DoChi2);
        if (_fail) return -1.0;
        for (Int_t i = 0 ; i < _nt; i++) {
          Double_t e= ereco[i];
          if (e<=0.0) continue;
          Double_t ypull = res[i] / e;
          chi2 += ypull*ypull;
        }
    }
    return chi2;
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::PrintTable (std::ostream& o, const Hist* hTrue, ErrorTreatment withError)
{
  // Prints entries from truth, measured, and reconstructed data for each bin.
  if (withError==kDefault) withError= _withError;
  if (withError==kDefault) withError= kErrors;
  const Hist* hReco= Hreco (withError);
  if (!_unfolded) return;
  Double_t chi_squ= -999.0;
  if (hTrue && (withError==kCovariance || withError==kCovToy)) chi_squ = Chi2(hTrue,withError);
  printTable (o, response()->Htruth(), response()->Hmeasured(), hTrue, Hmeasured(), hReco,
              _nm, _nt, _overflow, withError, chi_squ);
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
RooUnfoldT<Hist,Hist2D>::Hreco (ErrorTreatment withError)
{
    /*Creates reconstructed distribution. Error calculation varies by withError:
    0: No errors
    1: Errors from the square root of the diagonals of the covariance matrix given by the unfolding
    2: Errors from the square root of of the covariance matrix given by the unfolding
    3: Errors from the square root of the covariance matrix from the variation of the results in toy MC tests
    */

  if (!UnfoldWithErrors (withError)) withError= kNoError;
  if (!_unfolded){
    Hist* reco= copy(_res->Htruth(),true,GetName(),GetTitle());
    return reco;
  } else {

    TVectorD rec(this->_nt);
    TVectorD errors(this->_nt);
    for (Int_t i= 0; i < _nt; i++) {
      rec[i] = _rec[i];
      if        (withError==kErrors){
        errors[i] = (sqrt( fabs (_variances(i))));
      } else if (withError==kCovariance){
        errors[i] = (sqrt (fabs (_cov(i,i))));
      } else if (withError==kCovToy){
        errors[i] = (sqrt (fabs (_err_mat(i,i))));
      }
    }

    const Hist* t = _res->Htruth();
    Hist* reco = createHist<Hist>(rec,errors,t->GetName(),t->GetTitle(),vars(t),_overflow);
    return reco;
  }
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetSettings()
{
    //Gets maximum and minimum parameters and step size
    _minparm=0;
    _maxparm=0;
    _stepsizeparm=0;
    _defaultparm=0;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetMinParm() const
{
    //Get minimum regularisation parameter for unfolding method
    return _minparm;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetMaxParm() const
{
    //Get maximum regularisation parameter for unfolding method
    return _maxparm;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetStepSizeParm() const
{
    //Get suggested step size for unfolding distribution
    return _stepsizeparm;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetDefaultParm() const
{
    //Get suggested regularisation parameter.
    return _defaultparm;
}

template<class Hist,class Hist2D> RooUnfoldT<Hist,Hist2D>*
RooUnfoldT<Hist,Hist2D>::RunToy() const
{
  // Returns new RooUnfold object with smeared measurements and
  // (if IncludeSystematics) response matrix for use as a toy.
  // Use multiple toys to find spread of unfolding results.
  TString name= GetName();
  name += "_toy";
  RooUnfoldT<Hist,Hist2D>* unfold = Clone(name);

  // Make new smeared response matrix
  if (_dosys) unfold->SetResponse (_res->RunToy(), kTRUE);
  if (_dosys==2) return unfold;

  if (_haveCovMes) {

    // _covL is a lower triangular matrix for which the covariance matrix, V = _covL * _covL^T.
    if (!_covL) {
      TDecompChol c(*_covMes);
      c.Decompose();
      TMatrixD U(c.GetU());
      _covL= new TMatrixD (TMatrixD::kTransposed, U);
      if (_verbose>=2) printMatrix(*_covL,"decomposed measurement covariance matrix");
    }
    TVectorD newmeas(_nm);
    for (Int_t i= 0; i<_nm; i++) newmeas[i]= gRandom->Gaus(0.0,1.0);
    newmeas *= *_covL;
    newmeas += Vmeasured();
    unfold->SetMeasured(newmeas,*_covMes);

  } else {

    TVectorD newmeas= Vmeasured();
    const TVectorD& err= Emeasured();
    for (Int_t i= 0; i<_nm; i++) {
      Double_t e= err[i];
      if (e>0.0) newmeas[i] += gRandom->Gaus(0,e);
    }
    unfold->SetMeasured(newmeas,err);

  }
  return unfold;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Print(Option_t* /*opt*/) const
{
  cout << ClassName() << "::" << GetName() << " \"" << GetTitle()
       << "\", regularisation parameter=" << GetRegParm() << ", ";
  if (_haveCovMes) cout << "with measurement covariance, ";
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

template<class Hist,class Hist2D> TMatrixD
RooUnfoldT<Hist,Hist2D>::CutZeros(const TMatrixD& ereco)
{
    //Removes row & column if all their elements are 0.
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
RooUnfoldT<Hist,Hist2D>::Ereco(ErrorTreatment withError)
{
    /*Returns covariance matrices for error calculation of type withError
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests
    */
    TMatrixD Ereco_m(_nt,_nt);
    if (!UnfoldWithErrors (withError)) return Ereco_m;

    switch(withError){
      case kNoError:
        for (int i=0; i<_nt; i++){
          Ereco_m(i,i)=_rec(i);
        }
        break;
      case kErrors:
        for (int i=0; i<_nt;i++){
          Ereco_m(i,i)=_variances(i);
        }
        break;
      case kCovariance:
        Ereco_m=_cov;
        break;
      case kCovToy:
        Ereco_m=_err_mat;
        break;
      default:
        cerr<<"Error, unrecognised error method= "<<withError<<endl;
    }
    return Ereco_m;
}

template<class Hist,class Hist2D> TVectorD
RooUnfoldT<Hist,Hist2D>::ErecoV(ErrorTreatment withError)
{
    /*Returns vector of unfolding errors computed according to the withError flag:
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests
    */
    TVectorD Ereco_v(_nt);
    if (!UnfoldWithErrors (withError)) return Ereco_v;

    switch(withError){
      case kNoError:
        for (int i=0; i<_nt; i++){
          Ereco_v(i)=sqrt (fabs (_rec(i)));
        }
        break;
      case kErrors:
        for (int i=0; i<_nt; i++){
          Ereco_v(i)=sqrt (fabs (_variances(i)));
        }
        break;
      case kCovariance:
        for (int i=0; i<_nt; i++){
          Ereco_v(i)=sqrt (fabs (_cov(i,i)));
        }
        break;
      case kCovToy:
        for (int i=0; i<_nt; i++){
          Ereco_v(i)=sqrt (fabs (_err_mat(i,i)));
        }
        break;
      default:
        cerr<<"Error, unrecognised error method= "<<withError<<endl;
    }
    return Ereco_v;
}

template<class Hist,class Hist2D> TMatrixD
RooUnfoldT<Hist,Hist2D>::Wreco(ErrorTreatment withError)
{
    TMatrixD Wreco_m(_nt,_nt);
    if (!UnfoldWithErrors (withError, true)) return Wreco_m;

    switch(withError){
      case kNoError:
        for (int i=0; i<_nt; i++){
          if (_rec(i)!=0.0) Wreco_m(i,i)=1.0/_rec(i);
        }
        break;
      case kErrors:
        for (int i=0; i<_nt;i++){
          Wreco_m(i,i)=_wgt(i,i);
        }
        break;
      case kCovariance:
        Wreco_m=_wgt;
        break;
      case kCovToy:
        InvertMatrix (_err_mat, Wreco_m, "covariance matrix from toys", _verbose);
        break;
      default:
        cerr<<"Error, unrecognised error method= "<<withError<<endl;
    }
    return Wreco_m;
}

template<class Hist,class Hist2D> Hist*
RooUnfoldT<Hist,Hist2D>::HistNoOverflow (const Hist* h, Bool_t overflow)
{
  return histNoOverflow(h,overflow);
}

template<class Hist,class Hist2D> Int_t
RooUnfoldT<Hist,Hist2D>::InvertMatrix(const TMatrixD& mat, TMatrixD& inv, const char* name, Int_t verbose)
{
  // Invert a matrix using Single Value Decomposition: inv = mat^-1.
  // Can use InvertMatrix(mat,mat) to invert in-place.
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
  // Stream an object of class RooUnfold.
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    RooUnfoldT<Hist,Hist2D>::Class()->ReadBuffer  (R__b, this);
    TH1::AddDirectory (oldstat);
  } else {
    RooUnfoldT<Hist,Hist2D>::Class()->WriteBuffer (R__b, this);
  }
}

template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::RooUnfoldT()
  : TNamed()
{
  // Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const char*    name, const char*    title)
  : TNamed(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const TString& name, const TString& title)
  : TNamed(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>::~RooUnfoldT()
{
  Destroy();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>& RooUnfoldT<Hist,Hist2D>::operator= (const RooUnfoldT<Hist,Hist2D>& rhs)
{
  // Assignment operator for copying RooUnfold settings.
  Assign(rhs);
  return *this;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::verbose() const
{
  // Get verbosity setting which controls amount of information to be printed
  return _verbose;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::NToys()     const
{
  // Get number of toys used in kCovToy error calculation.
  return _NToys;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::Overflow()  const
{
  // Histogram under/overflow bins are used?
  return _overflow;
}

template<class Hist,class Hist2D> 
const RooUnfoldResponseT<Hist,Hist2D>* RooUnfoldT<Hist,Hist2D>::response()  const
{
   // Response matrix object
  return _res;
}

template<class Hist,class Hist2D> 
const Hist*               RooUnfoldT<Hist,Hist2D>::Hmeasured() const
{
  // Measured Distribution as a histogram
  return _meas;
}

template<class Hist,class Hist2D> 
TVectorD&                RooUnfoldT<Hist,Hist2D>::Vreco()
{
  // Unfolded (reconstructed) distribution as a vector
  if (!_unfolded) {
    if (!_fail) Unfold();
    if (!_unfolded) {
      _fail= true;
      if (_nt > 0 && _rec.GetNrows() == 0) _rec.ResizeTo(_nt);   // need something
    }
  }
  return _rec;
}

template<class Hist,class Hist2D> 
const TVectorD&          RooUnfoldT<Hist,Hist2D>::Vmeasured() const
{
  // Measured distribution as a vector.
  if (!_vMes){
    _vMes = new TVectorD(h2v (_meas, _overflow));
  }
  return *_vMes;
}

template<class Hist,class Hist2D> 
const TVectorD&          RooUnfoldT<Hist,Hist2D>::Emeasured() const
{
  // Measured errors as a vector.
  if (!_eMes){
    _eMes = new TVectorD(h2ve (_meas, _overflow));
  }
  return *_eMes;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetVerbose (Int_t level)
{
  // Set verbosity level which controls amount of information to be printed
  _verbose= level;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetNToys (Int_t toys)
{
  // Set number of toys used in kCovToy error calculation.
  _NToys= toys;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetRegParm (Double_t)
{
  // Set Regularisation parameter
}

template<class Hist,class Hist2D> 
Double_t RooUnfoldT<Hist,Hist2D>::GetRegParm() const
{
  // Get regularisation parameter.
  return -1;
}

template<class Hist,class Hist2D> 
void RooUnfoldT<Hist,Hist2D>::IncludeSystematics (Int_t dosys)
{
  // Include systematic errors from response matrix?
  // Use dosys=2 to exclude measurement errors.
  if (dosys!=_dosys) _haveWgt= _haveErrors= _haveCov= _have_err_mat= kFALSE;
  _dosys= dosys;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::SystematicsIncluded() const
{
  // return setting for whether to include systematic errors from response matrix
  return _dosys;
}

template class RooUnfoldT<TH1,TH2>;
ClassImp (RooUnfold);

template class RooUnfoldT<RooAbsReal,RooAbsReal>;
ClassImp (RooFitUnfold);


