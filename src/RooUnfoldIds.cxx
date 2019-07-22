// Author: Bogdan Malaescu <bogdan.malaescu@cern.ch>
// Author: Christopher Meyer <chris.meyer@cern.ch>
//
// Inspired by Tim Adye code for RooUnfoldSvd
// For support, contact: chris.meyer@cern.ch

#include "RooUnfoldIds.h"
#include "RooUnfoldTH1Helpers.h"
#ifndef NOROOFIT
#include "RooUnfoldFitHelpers.h"
#endif


#include <iostream>

#include "TClass.h"
#include "TBuffer.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH2D.h"

#include "RooUnfoldHelpers.h"
#include "RooUnfoldResponse.h"

using namespace RooUnfolding;

ClassImp(RooUnfoldIds)

//______________________________________________________________________________
template<class Hist,class Hist2D>
RooUnfoldIdsT<Hist,Hist2D>::RooUnfoldIdsT(const RooUnfoldIdsT &rhs)
: RooUnfoldT<Hist,Hist2D>(rhs)
{
   // Copy constructor
   Init();
   CopyData(rhs);
}

//______________________________________________________________________________
template<class Hist,class Hist2D>
RooUnfoldIdsT<Hist,Hist2D>::RooUnfoldIdsT(const RooUnfoldResponseT<Hist,Hist2D> *res, const Hist *meas, Int_t niter, const char* name, const char* title)
: RooUnfoldT<Hist,Hist2D>(res, meas)
, _niter(niter)
, _lambdaL(0.)
, _lambdaUmin(0.5)
, _lambdaMmin(0.5)
, _lambdaS(0.)
{
   // Constructor with response matrix object and measured unfolding input histogram.
   Init();
}

// //______________________________________________________________________________
// RooUnfoldIds*
// RooUnfoldIdsT<Hist,Hist2D>::Clone(const char *newname) const
// {
//    RooUnfoldIds *unfold = new RooUnfoldIdsT(*this);
//    if (newname && strlen(newname)) unfold->SetName(newname);
//    return unfold;
// }

// //______________________________________________________________________________

template<class Hist,class Hist2D> void
RooUnfoldIdsT<Hist,Hist2D>::Reset()
{
   Destroy();
   Init();
   RooUnfoldT<Hist,Hist2D>::Reset();
}

//______________________________________________________________________________
template<class Hist,class Hist2D> void
RooUnfoldIdsT<Hist,Hist2D>::Destroy()
{
   delete _meas1d;
   delete _train1d;
   delete _truth1d;
   delete _reshist;
}

//______________________________________________________________________________
template<class Hist,class Hist2D> void
RooUnfoldIdsT<Hist,Hist2D>::Init()
{
   _meas1d = _train1d = _truth1d = 0;
   _reshist = 0;
   this->GetSettings();
}

//______________________________________________________________________________
template<class Hist,class Hist2D> void
RooUnfoldIdsT<Hist,Hist2D>::Assign(const RooUnfoldIdsT &rhs)
{
  RooUnfoldT<Hist,Hist2D>::Assign(rhs);
  CopyData(rhs);
}

//______________________________________________________________________________
template<class Hist,class Hist2D>void
RooUnfoldIdsT<Hist,Hist2D>::CopyData(const RooUnfoldIdsT&rhs)
{
   _niter = rhs._niter;
   _lambdaL = rhs._lambdaL;
   _lambdaUmin = rhs._lambdaUmin;
   _lambdaMmin = rhs._lambdaMmin;
   _lambdaS = rhs._lambdaS;
}

namespace{
  TH1* histNoOverflow(const TH1* hist, bool overflow){
    return createHist<TH1>(h2v(hist,overflow),h2ve(hist,overflow),name(hist),title(hist),vars(hist),overflow);
  }
}

  
//______________________________________________________________________________

template<class Hist,class Hist2D> void
RooUnfoldIdsT<Hist,Hist2D>::Unfold() const
{
   // Data and MC reco/truth must have the same number of bins
   if (this->_res->HasFakes()) {
      _nb = this->_nt+1;
      if (this->_nm>_nb) _nb = this->_nm;
   } else {
      _nb = this->_nm > this->_nt ? this->_nm : this->_nt;
   }

   Bool_t oldstat= TH1::AddDirectoryStatus();
   TH1::AddDirectory (kFALSE);

   TH1* _meas1d  = ::histNoOverflow(this->_meas            , this->_overflow); // data
   TH1* _train1d = ::histNoOverflow(this->_res->Hmeasured(), this->_overflow); // reco
   TH1* _truth1d = ::histNoOverflow(this->_res->Htruth()   , this->_overflow); // true
   TH2* _reshist = this->_res->HresponseNoOverflow();

   RooUnfolding::resize(_meas1d,  _nb);
   RooUnfolding::resize(_train1d, _nb);
   RooUnfolding::resize(_truth1d, _nb);
   RooUnfolding::resize(_reshist, _nb, _nb);

   // Something about fakes here?
   if (this->_res->HasFakes()) {
      TVectorD fakes = this->_res->Vfakes();
      Double_t nfakes = fakes.Sum();
      if (this->_verbose >= 1) std::cout << "Add truth bin for " << nfakes << " fakes" << std::endl;
      for (Int_t i = 0; i < this->_nm; ++i) _reshist->SetBinContent(i+1, this->_nt+1, fakes[i]);
      _truth1d->SetBinContent(this->_nt+1, nfakes);
   }

   if (this->_verbose >= 1) std::cout << "IDS init " << _reshist->GetNbinsX() << " x " << _reshist->GetNbinsY() << std::endl;

   // Perform IDS unfolding
   TH1 *rechist = dynamic_cast<TH1*>(GetIDSUnfoldedSpectrum(_train1d, _truth1d, _reshist, _meas1d, _niter));

   this->_cache._rec.ResizeTo(this->_nt);
   for (Int_t i = 0; i < this->_nt; ++i) {
     this->_cache._rec[i] = rechist->GetBinContent(i+1);
   }

   delete rechist;
   TH1::AddDirectory(oldstat);

   this->_cache._unfolded = kTRUE;
   this->_cache._haveCov = kFALSE;
}

//______________________________________________________________________________
template<class Hist,class Hist2D>void
RooUnfoldIdsT<Hist,Hist2D>::GetCov() const
{
   if (!_meas1d) return;

   Bool_t oldstat = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   TH2 *meascov = new TH2D("meascov", "meascov", _nb, 0.0, 1.0, _nb, 0.0, 1.0);
   const TMatrixD& cov = this->GetMeasuredCov();
   for (Int_t i = 0; i < this->_nm; ++i)
     for (Int_t j = 0; j < this->_nm; ++j)
       meascov->SetBinContent(i+1, j+1, cov(i,j));

   // Need to fill _cov with unfolded result
   TH2 *unfoldedCov = GetUnfoldCovMatrix(meascov, this->_NToys);
   TH2 *adetCov     = GetAdetCovMatrix(this->_NToys);


   this->_cache._cov.ResizeTo(this->_nt, this->_nt);
   for (Int_t i = 0; i < this->_nt; i++) {
      for (Int_t j = 0; j < this->_nt; ++j) {
         this->_cache._cov(i,j) = unfoldedCov->GetBinContent(i+1, j+1) + adetCov->GetBinContent(i+1, j+1);
      }
   }

   delete adetCov;
   delete unfoldedCov;
   delete meascov;
   TH1::AddDirectory(oldstat);

   this->_cache._haveCov = kTRUE;
}

//______________________________________________________________________________
template<>TH2*
RooUnfoldIdsT<TH1,TH2>::GetUnfoldCovMatrix(const TH2 *cov, Int_t ntoys, Int_t seed) const
{
   // Determine for given input error matrix covariance matrix of unfolded
   // spectrum from toy simulation given the passed covariance matrix on measured spectrum
   // "cov"    - covariance matrix on the measured spectrum, to be propagated
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   // Note that this covariance matrix will contain effects of forced normalisation if spectrum is normalised to unit area.

   TH1* _meas1d  = ::histNoOverflow(this->_meas            , this->_overflow); // data
   TH1* _train1d = ::histNoOverflow(this->_res->Hmeasured(), this->_overflow); // reco
   TH1* _truth1d = ::histNoOverflow(this->_res->Htruth()   , this->_overflow); // true
   TH2* _reshist = this->_res->HresponseNoOverflow();

   TH1* unfres = 0;
   TH2* unfcov = (TH2*)_reshist->Clone("unfcovmat");
   unfcov->SetTitle("Toy covariance matrix");
   for (Int_t i = 1; i <= _nb; ++i)
      for(Int_t j = 1; j <= _nb; ++j)
         unfcov->SetBinContent(i, j, 0.);

   // Code for generation of toys (taken from TSVDUnfold [took from RooResult] and modified)
   // Calculate the elements of the upper-triangular matrix L that
   // gives Lt*L = C, where Lt is the transpose of L (the "square-root method")
   TMatrixD L(_nb, _nb); L *= 0;

   for (Int_t iPar = 0; iPar < _nb; ++iPar) {

      // Calculate the diagonal term first
      L(iPar, iPar) = cov->GetBinContent(iPar+1, iPar+1);
      for (Int_t k = 0; k < iPar; ++k) L(iPar, iPar) -= TMath::Power(L(k, iPar), 2);
      if (L(iPar, iPar) > 0.0) L(iPar, iPar) = TMath::Sqrt(L(iPar,iPar));
      else                     L(iPar, iPar) = 0.0;

      // ...then the off-diagonal terms
      for (Int_t jPar = iPar+1; jPar < _nb; ++jPar) {
         L(iPar, jPar) = cov->GetBinContent(iPar+1, jPar+1);
         for (Int_t k = 0; k < iPar; k++) L(iPar, jPar) -= L(k, iPar)*L(k, jPar);
         if (L(iPar,iPar) != 0.) L(iPar, jPar) /= L(iPar, iPar);
         else                    L(iPar, jPar) = 0;
      }
   }

   // Remember it
   TMatrixD *Lt = new TMatrixD(TMatrixD::kTransposed, L);
   TRandom3 random(seed);

   // Needed to build covariance matrix
   TVectorD avgtoy(_nb);
   TMatrixD toys(ntoys, _nb);
   for (Int_t j = 0; j < _nb; ++j) {
      avgtoy[j] = 0.0;
      for (Int_t i = 0; i < ntoys; ++i) {
         toys[i][j] = 0.0;
      }
   }

   // Get the mean of the toys first
   TH1 *toyhist = (TH1*)_meas1d->Clone("toyhisto");
   for (Int_t i = 0; i < ntoys; i++) {

      // create a vector of unit Gaussian variables
      TVectorD g(_nb);
      for (Int_t k = 0; k < _nb; ++k) g(k) = random.Gaus(0.,1.);

      // Multiply this vector by Lt to introduce the appropriate correlations
      g *= (*Lt);

      // Add the mean value offsets and store the results
      for (Int_t j = 1; j <= _nb; ++j) {
         toyhist->SetBinContent(j, _meas1d->GetBinContent(j) + g(j-1));
         toyhist->SetBinError(j, _meas1d->GetBinError(j));
      }

      // Perform IDS unfolding
      unfres = dynamic_cast<TH1*>(GetIDSUnfoldedSpectrum(_train1d, _truth1d, _reshist, toyhist, _niter));

      for (Int_t j = 0; j < _nb; ++j) {
         toys[i][j] = unfres->GetBinContent(j+1);
         avgtoy[j] += unfres->GetBinContent(j+1)/ntoys;
      }

      delete unfres;

   }
   delete toyhist;
   delete Lt;

   for (Int_t i = 0; i < ntoys; ++i) {
      for (Int_t j = 0; j < _nb; ++j) {
         for (Int_t k = 0; k < _nb; ++k) {
            unfcov->SetBinContent(j+1, k+1, unfcov->GetBinContent(j+1, k+1) + (toys[i][j] - avgtoy[j])*(toys[i][k] - avgtoy[k])/ntoys);
         }
      }
   }

   return unfcov;
}

//______________________________________________________________________________
template<>TH2*
RooUnfoldIdsT<TH1,TH2>::GetAdetCovMatrix(Int_t ntoys, Int_t seed) const
{
   // Determine covariance matrix of unfolded spectrum from finite statistics in
   // response matrix using pseudo experiments
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments

   TH1* _meas1d  = ::histNoOverflow(this->_meas            , this->_overflow); // data
   TH1* _train1d = ::histNoOverflow(this->_res->Hmeasured(), this->_overflow); // reco
   TH1* _truth1d = ::histNoOverflow(this->_res->Htruth()   , this->_overflow); // true
   TH2* _reshist = this->_res->HresponseNoOverflow();

   TH1 *unfres = 0;
   TH2 *unfcov = (TH2*)_reshist->Clone("unfcovmat");
   unfcov->SetTitle("Toy covariance matrix");
   for(Int_t i = 1; i <= _nb; ++i)
      for(Int_t j = 1; j <= _nb; ++j)
         unfcov->SetBinContent(i, j, 0.0);

   //Now the toys for the detector response matrix
   TRandom3 random(seed);

   // Needed to build covariance matrix
   TVectorD avgtoy(_nb);
   TMatrixD toys(ntoys, _nb);
   for (Int_t j = 0; j < _nb; ++j) {
      avgtoy[j] = 0.0;
      for (Int_t i = 0; i < ntoys; ++i) {
         toys[i][j] = 0.0;
      }
   }

   TH2 *toymat = (TH2*)_reshist->Clone("toymat");
   Double_t fluc = -1.0;
   for (Int_t i = 0; i < ntoys; ++i) {
      for (Int_t k = 1; k <= _nb; ++k) {
         for (Int_t m = 1; m <= _nb; ++m) {
            if (_reshist->GetBinContent(k, m)) {
               // fToymat->SetBinContent(k, m, random.Poisson(fAdet->GetBinContent(k,m)));
               fluc = -1.0;
               while (fluc < 0.0) {
                  fluc = random.Gaus(_reshist->GetBinContent(k, m), _reshist->GetBinError(k, m));
               }
               
               toymat->SetBinContent(k, m, fluc);
            }
         }
      }

      // Perform IDS unfolding
      unfres = dynamic_cast<TH1*>(GetIDSUnfoldedSpectrum(_train1d, _truth1d, toymat, _meas1d, _niter));

      for (Int_t j = 0; j < _nb; ++j) {
         toys[i][j] = unfres->GetBinContent(j+1);
         avgtoy[j] += unfres->GetBinContent(j+1)/ntoys;
      }

      delete unfres;
   }

   delete toymat;

   for (Int_t i = 0; i < ntoys; ++i) {
      for (Int_t j = 0; j < _nb; ++j) {
         for (Int_t k = 0; k < _nb; ++k) {
            unfcov->SetBinContent(j+1, k+1, unfcov->GetBinContent(j+1, k+1) + (toys[i][j] - avgtoy[j])*(toys[i][k] - avgtoy[k])/ntoys);
         }
      }
   }

   return unfcov;
}

//______________________________________________________________________________
template<class Hist,class Hist2D>Hist*
RooUnfoldIdsT<Hist,Hist2D>::GetIDSUnfoldedSpectrum(const Hist *h_RecoMC, const Hist *h_TruthMC, const Hist2D *h_2DSmear, const Hist *h_RecoData, Int_t iter) const
{

  Int_t nbinsx = h_RecoData->GetNbinsX();
  Int_t nbinsy = h_RecoData->GetNbinsY();
  Int_t nbins  = nbinsx*nbinsy;

   // Sanity checks
   if (h_TruthMC->GetNbinsX() != nbinsx || h_TruthMC->GetNbinsY() != nbinsy ||
       h_RecoMC->GetNbinsX()  != nbinsx || h_RecoMC->GetNbinsY()  != nbinsy ||
       h_2DSmear->GetNbinsX() != nbins) {
      std::cout << "Bins of input histograms don't all match, exiting IDS unfolding and returning NULL." << std::endl;
      return NULL;
   }

   // Put inputs into vectors, and if necessary turn 2-D inputs into 1-D inputs
   TVectorD reco(nbins), truth(nbins);
   TVectorD data(nbins), dataerror(nbins);
   Int_t i = 0;
   for (Int_t by = 1; by <= nbinsy; ++by) { // loop over pt
      for (Int_t bx = 1; bx <= nbinsx; ++bx) { // loop over gap_size, for each pt value
         reco[i]  = h_RecoMC->GetBinContent(bx, by);
         truth[i] = h_TruthMC->GetBinContent(bx, by);
         data[i]  = h_RecoData->GetBinContent(bx, by);
         if (data[i] > 0.0) {
            dataerror[i] = h_RecoData->GetBinError(bx, by);
         } else {
            dataerror[i] = 1.0;
         }
         i++;
      }
   }

   // Make transfer matrix and project matched MC spectra
   TMatrixD migmatrix(nbins, nbins);
   TVectorD recomatch(nbins);
   TVectorD truthmatch(nbins);
   for (Int_t i = 0; i < nbins; ++i) {
      recomatch[i] = 0.0;
      truthmatch[i] = 0.0;
      for (Int_t j = 0; j < nbins; ++j) {
         recomatch[i]   += h_2DSmear->GetBinContent(i+1, j+1);
         truthmatch[i]  += h_2DSmear->GetBinContent(j+1, i+1);
         migmatrix[i][j] = h_2DSmear->GetBinContent(i+1, j+1);
      }
   }

   // Apply matching inefficiency from reco MC to data
   for (Int_t i = 0; i < nbins; ++i) {
      if (reco[i] > 0.0) {
         data[i] *= recomatch[i]/reco[i];
         dataerror[i] *= recomatch[i]/reco[i];
      } else {
         data[i] = 0.0;
         dataerror[i] = 0.0;
      }
   }

   // Matrix unfolding
   // Double_t lambdaL = 0.;
   // Double_t lambdaUmin = 0.0;
   // Double_t lambdaMmin = 0.0;
   // Double_t lambdaS = 0.;

   TVectorD result0(nbins), result(nbins);
   PerformIterations(data, dataerror, migmatrix, nbins,
                     _lambdaL, iter, _lambdaUmin, _lambdaMmin, _lambdaS,
                     &result0, &result);

   // Apply matching efficiency from truth MC to unfolded matched data
   for (Int_t i = 0; i < nbins; ++i) {
      if (truthmatch[i] > 0.0) {
         result[i] *= truth[i]/truthmatch[i];
      } else {
         result[i] = 0.0;
      }
   }

   // Make 1-D or 2-D histogram
   TH1 *h_DataUnfolded = (TH1*)h_RecoData->Clone("unfolded");
   h_DataUnfolded->SetTitle("unfolded");
   h_DataUnfolded->Reset();

   i = 0;
   for (Int_t by = 1; by <= nbinsy; ++by) {
      for (Int_t bx = 1; bx <= nbinsx; ++bx) {
         h_DataUnfolded->SetBinContent(bx, by, result[i++]);
      }
   }

   // Return result
   return (Hist*)h_DataUnfolded;
}

//______________________________________________________________________________
template<class Hist,class Hist2D>Double_t
RooUnfoldIdsT<Hist,Hist2D>::Probability( Double_t deviation, Double_t sigma, Double_t lambda ) const
{
   if( lambda < 0.00001 ){ return 1.;}
   if( lambda > 1000. ){ return 0.;}
   return 1-exp(-pow(deviation/(sigma*lambda),2) );
}

//______________________________________________________________________________
template<class Hist,class Hist2D>Double_t
RooUnfoldIdsT<Hist,Hist2D>::MCnormalizationCoeff(const TVectorD *vd, const TVectorD *errvd, const TVectorD *vRecmc, const Int_t dim, const Double_t estNknownd, const Double_t Nmc, const Double_t lambda, const TVectorD *soustr_ ) const
{
   Double_t Nknownd = estNknownd;
   for(Int_t i=0; i<dim; i++){
      if( ((*vd)[i] - (*soustr_)[i] >= 0.) && ((*vRecmc)[i]!=0) ){
         // The first test shall also be done when computing Nmc & estNknownd
         Double_t ef = Probability( fabs((*vd)[i]- (*soustr_)[i] - estNknownd/Nmc*(*vRecmc)[i]), sqrt( pow((*errvd)[i],2) /* + pow(estNknownd/Nmc,2)*fabs((*vRecmc)[i]) */ ),lambda );
         Nknownd += (1-ef)*( (*vd)[i] - (*soustr_)[i] - estNknownd/Nmc*(*vRecmc)[i] );
      }
   }

   return Nknownd;
}

//______________________________________________________________________________
template<class Hist,class Hist2D>Double_t
RooUnfoldIdsT<Hist,Hist2D>::MCnormalizationCoeffIter(const TVectorD *vd, const TVectorD *errvd, const TVectorD *vRecmc, const Int_t dim, const Double_t estNknownd, const Double_t Nmc, const TVectorD *soustr_, Double_t lambdaN, Int_t NiterMax, Int_t messAct) const
{
   Double_t Nkd = 0., estNknownd_ = estNknownd;
   for(Int_t i=0; i<NiterMax; i++ ){
      Nkd = MCnormalizationCoeff( vd, errvd, vRecmc, dim, estNknownd_, Nmc, lambdaN, soustr_ );
      if( fabs((estNknownd_ -  Nkd)/Nkd) < 1e-6 ){ break; }
      if( (i >= NiterMax-1) && (messAct != 0)  ){ std::cout << "WARNING: a number of " << i+1 << " steps were performed for the normalization " << (estNknownd_ -  Nkd)/Nkd << "  " << (estNknownd -  Nkd)/Nkd << std::endl;}

      estNknownd_ = Nkd;

   }
   return Nkd;
}

//______________________________________________________________________________
template<class Hist,class Hist2D>void
RooUnfoldIdsT<Hist,Hist2D>::IdsUnfold( const TVectorD &b, const TVectorD &errb, const TMatrixD &A, const Int_t dim, const Double_t lambda, TVectorD *soustr_, TVectorD *unf) const
{
   // compute the mc true and reco spectra and normalize them
   TVectorD reco_mcN(dim), true_mcN(dim);
   Double_t estNkd = 0., Nkd = 0. , Nmc = 0.;
   for(Int_t i=0; i<dim; i++ ){
      reco_mcN[i] = 0.;
      true_mcN[i] = 0.;
      for(Int_t j=0; j<dim; j++ ){
         reco_mcN[i] += A[i][j];
         true_mcN[i] += A[j][i];
      }
      if(b[i] - (*soustr_)[i] >= 0.){
         Nmc += reco_mcN[i];
         estNkd += b[i] - (*soustr_)[i];
      }
   }

   // # known data
   Nkd = MCnormalizationCoeffIter(&b, &errb, &reco_mcN, dim, estNkd, Nmc, soustr_ );

   for(Int_t i=0; i<dim; i++ ){
      reco_mcN[i] *= Nkd/Nmc;
      true_mcN[i] *= Nkd/Nmc;
      (*unf)[i] = true_mcN[i] + (*soustr_)[i];
   }

   // compute the prob(j|i) matrix
   TMatrixD prob(dim, dim);
   for(Int_t i=0; i<dim; i++){
      Double_t si = 0.;
      for(Int_t j=0; j<dim; j++){
         si += A[i][j];
      }
      for(Int_t j=0; j<dim; j++){
         if( si!=0. )
            prob[i][j] = A[i][j]/si;
         else
            prob[i][j] = 0.;
      }

      if (reco_mcN[i] != 0.0 && (b[i] - (*soustr_)[i]) > 0.0 /*&&((*b)[i]>0.)*/) {
         Double_t ef = Probability( fabs(b[i]-(*soustr_)[i]-reco_mcN[i]), sqrt( pow(errb[i],2) /* + Nkd/Nmc*fabs((*reco_mcN)[i]) */ ), lambda );
         for(Int_t j=0; j<dim; j++){
            (*unf)[j] += ef * prob[i][j]*(b[i]-(*soustr_)[i]-reco_mcN[i]);
         }
         (*unf)[i] += (1-ef) * (b[i]-(*soustr_)[i]-reco_mcN[i]);
      } else {
         (*unf)[i] += b[i]-(*soustr_)[i]-reco_mcN[i];
      }
   }

}

//______________________________________________________________________________
template<class Hist,class Hist2D>void
RooUnfoldIdsT<Hist,Hist2D>::ComputeSoustrTrue( const TMatrixD *A, const TVectorD *unfres, const TVectorD *unfresErr, Int_t N, TVectorD *soustr_, Double_t lambdaS ) const
{

   TVectorD *true_mcT = new TVectorD(N), *active = new TVectorD(N);
   Double_t estNkd = 0., Nmc=0., NkUR=0.;
   for(Int_t j=0; j<N; j++){

      (*active)[j] = 1.;

      (*true_mcT)[j] = 0.;
      for(Int_t i=0; i<N; i++){
         (*true_mcT)[j] += ((*A)[i][j]);
      }
      if((*unfres)[j] - (*soustr_)[j] >= 0.){
         Nmc += (*true_mcT)[j];
         estNkd += (*unfres)[j] - (*soustr_)[j];
      }

      if( ((*true_mcT)[j]) < 0.){
         std::cout << "found problematic (*true_mcT)[j] " << (*true_mcT)[j] << " " << j << " " << (*unfres)[j] << std::endl;
         exit(1);
      }
   }

   for(Int_t k=0; k<2; k++){
      NkUR = MCnormalizationCoeffIter( unfres, unfresErr, true_mcT, N, estNkd, Nmc, soustr_ );

      for(Int_t j=0; j<N; j++){
         if( ((*unfres)[j] - (*soustr_)[j]>0.) && ((*true_mcT)[j]!=0.) && ((*active)[j]>0.) ){
            Double_t ef = Probability(fabs((*unfres)[j] /*- (*soustr_)[j]*/ -NkUR/Nmc*(*true_mcT)[j]), sqrt(pow((*unfresErr)[j],2) /* +pow(NkUR/Nmc,2)*fabs((*true_mcT)[j]) */ ), lambdaS);
            ((*soustr_)[j]) /*+*/= (1-ef) * ( (*unfres)[j] /*- (*soustr_)[j]*/ - (*true_mcT)[j]/(Nmc/NkUR));
         }
         else{
            ((*soustr_)[j]) /*+*/= (*unfres)[j] /*- (*soustr_)[j]*/ - (*true_mcT)[j]/(Nmc/NkUR);
            (*active)[j] = -1.;
         }
      }
      estNkd = NkUR;
   }

   delete true_mcT; delete active;
}

//______________________________________________________________________________
template<class Hist,class Hist2D>void
RooUnfoldIdsT<Hist,Hist2D>::ModifyMatrix( TMatrixD *Am, const TMatrixD *A, const TVectorD *unfres, const TVectorD *unfresErr, Int_t N, const Double_t lambdaM_, TVectorD *soustr_, const Double_t lambdaS_ ) const
{
   ComputeSoustrTrue( A, unfres, unfresErr, N, soustr_, lambdaS_ );

   TVectorD *true_mcT = new TVectorD(N);
   Double_t estNkd = 0., Nmc=0., NkUR=0.;
   for(Int_t j=0; j<N; j++){
      (*true_mcT)[j] = 0.;
      for(Int_t i=0; i<N; i++){
         (*true_mcT)[j] += ((*A)[i][j]);
         ((*Am)[i][j]) = ((*A)[i][j]);
      }
      if((*unfres)[j] - (*soustr_)[j] >= 0.){
         Nmc += (*true_mcT)[j];
         estNkd += (*unfres)[j] - (*soustr_)[j];
      }

      if( ((*true_mcT)[j]) < 0.){
         std::cout << "found problematic (*true_mcT)[j] " << (*true_mcT)[j] << " " << j << " " << (*unfres)[j] << std::endl;
         exit(1);
      }
   }

   NkUR = MCnormalizationCoeffIter( unfres, unfresErr, true_mcT, N, estNkd, Nmc, soustr_ );

   for(Int_t j=0; j<N; j++){
      if( (*unfres)[j] - (*soustr_)[j]>0. && (*true_mcT)[j]!=0. ){
         Double_t ef = Probability(fabs((*unfres)[j] - (*soustr_)[j] -NkUR/Nmc*(*true_mcT)[j]), sqrt(pow((*unfresErr)[j],2) /* +pow(NkUR/Nmc,2)*fabs((*true_mcT)[j]) */ ), lambdaM_);
         for(Int_t i=0; i<N; i++){
            ((*Am)[i][j]) += ef * ( ((*unfres)[j] - (*soustr_)[j])*(Nmc/NkUR) - (*true_mcT)[j]) * ((*A)[i][j])/((*true_mcT)[j]);
         }
      }
   }

   delete true_mcT;
}

//______________________________________________________________________________
template<class Hist,class Hist2D>void
RooUnfoldIdsT<Hist,Hist2D>::PerformIterations(const TVectorD &data, const TVectorD &dataErr, const TMatrixD &A_, const Int_t &N_, const Double_t lambdaL_, const Int_t NstepsOptMin_, const Double_t lambdaU_, const Double_t lambdaM_, const Double_t lambdaS_, TVectorD* unfres1IDS_, TVectorD* unfres2IDS_) const
{
   TVectorD soustr(N_);
   for (Int_t i = 0; i < N_; i++) soustr[i] = 0.;

   IdsUnfold(data, dataErr, A_, N_, lambdaL_, &soustr, unfres1IDS_); // 1 step
   
   for (Int_t i = 0; i < N_; i++) (*unfres2IDS_)[i] = (*unfres1IDS_)[i];

   TMatrixD Am_(N_, N_);
   for (Int_t k = 0; k < NstepsOptMin_; k++) {
      ModifyMatrix(&Am_, &A_, unfres2IDS_, &dataErr, N_, lambdaM_, &soustr, lambdaS_);

      // UNFOLDING
      IdsUnfold(data, dataErr, Am_, N_, lambdaU_, &soustr, unfres2IDS_); // full iterations
   }
}

//______________________________________________________________________________
template<class Hist,class Hist2D>TMatrixD*
RooUnfoldIdsT<Hist,Hist2D>::GetSqrtMatrix( const TMatrixD& covMat )
{
   // production of square-root matrix (required for correlated random number generation)
   Double_t sum = 0;
   const Int_t size = covMat.GetNrows();;
   TMatrixD* sqrtMat = new TMatrixD( size, size );

   for (Int_t i=0; i<size; i++) {

      sum = 0;
      for (Int_t j=0;j< i; j++) sum += (*sqrtMat)(i,j) * (*sqrtMat)(i,j);

      (*sqrtMat)(i,i) = TMath::Sqrt(TMath::Abs(covMat(i,i) - sum));
      if ((*sqrtMat)(i,i) <= 0) {
         std::cout << "<GetSqrtMatrix> Covariance matrix has incomplete rang" << std::endl;
         exit(1);
      }

      for (Int_t k=i+1 ;k<size; k++) {

         sum = 0;
         for (Int_t l=0; l<i; l++) sum += (*sqrtMat)(k,l) * (*sqrtMat)(i,l);

         (*sqrtMat)(k,i) = (covMat(k,i) - sum) / (*sqrtMat)(i,i);

      }
   }
   return sqrtMat;
}

//______________________________________________________________________________
template<class Hist,class Hist2D>void
RooUnfoldIdsT<Hist,Hist2D>::GenGaussRnd( TArrayD& v, const TMatrixD& sqrtMat, TRandom3& R ) const
{
   // generate vector of correlated Gaussian-distributed random numbers
   // sanity check
   const Int_t size = sqrtMat.GetNrows();
   if (size != v.GetSize()) {
      std::cout << "<GenGaussRnd1> Too short input vector: " << size << " " << v.GetSize() << std::endl;
      exit(1);
   }

   Double_t* tmpVec = new Double_t[size];

   for (Int_t i=0; i<size; i++) {
      Double_t x, y, z;
      y = R.Rndm();
      z = R.Rndm();
      x = 2*TMath::Pi()*z;
      tmpVec[i] = TMath::Sin(x) * TMath::Sqrt(-2.0*TMath::Log(y));
   }

   for (Int_t i=0; i<size; i++) {
      v[i] = 0;
      for (Int_t j=0; j<=i; j++) v[i] += sqrtMat(i,j) * tmpVec[j];
   }

   delete [] tmpVec;
}

//______________________________________________________________________________
// template<class Hist,class Hist2D>void
// RooUnfoldIdsT<Hist,Hist2D>::Streamer(TBuffer &R__b)
// {
//    // Stream an object of class RooUnfoldIds.
//    if (R__b.IsReading()) {
//       // Don't add our histograms to the currect directory.
//       // We own them and we don't want them to disappear when the file is closed.
//       Bool_t oldstat = TH1::AddDirectoryStatus();
//       TH1::AddDirectory(kFALSE);
//       RooUnfoldIdsT<Hist,Hist2D>::Class()->ReadBuffer  (R__b, this);
//       TH1::AddDirectory(oldstat);
//    } else {
//       RooUnfoldIdsT<Hist,Hist2D>::Class()->WriteBuffer (R__b, this);
//    }
// }



template<class Hist,class Hist2D>
RooUnfoldIdsT<Hist,Hist2D>::RooUnfoldIdsT()
: RooUnfoldT<Hist,Hist2D>()
{
   // Default constructor. Use Setup() to prepare for unfolding.
   Init();
}

template<class Hist,class Hist2D>
RooUnfoldIdsT<Hist,Hist2D>::RooUnfoldIdsT(const char *name, const char *title)
: RooUnfoldT<Hist,Hist2D>(name, title)
{
   // Basic named constructor. Use Setup() to prepare for unfolding.
   Init();
}

template<class Hist,class Hist2D>
RooUnfoldIdsT<Hist,Hist2D>::RooUnfoldIdsT(const TString &name, const TString &title)
: RooUnfoldT<Hist,Hist2D>(name, title)
{
   // Basic named constructor. Use Setup() to prepare for unfolding.
   Init();
}

template<class Hist,class Hist2D>
RooUnfoldIdsT<Hist,Hist2D>& RooUnfoldIdsT<Hist,Hist2D>::operator=(const RooUnfoldIdsT<Hist,Hist2D>& rhs)
{
   // Assignment operator for copying RooUnfoldIdsTsettings.
   Assign(rhs);
   return *this;
}

template<class Hist,class Hist2D>
RooUnfoldIdsT<Hist,Hist2D>::~RooUnfoldIdsT()
{
   Destroy();
}

template<class Hist,class Hist2D>
void  RooUnfoldIdsT<Hist,Hist2D>::SetRegParm (Double_t parm)
{
  // Set regularisation parameter (number of iterations)
  SetNIter(Int_t(parm+0.5));
}

template<class Hist,class Hist2D>
Double_t RooUnfoldIdsT<Hist,Hist2D>::GetRegParm() const
{
  // Return regularisation parameter (number of iterations)
  return GetNIter();
}

template<class Hist,class Hist2D>
void RooUnfoldIdsT<Hist,Hist2D>::SetNIter(Int_t niter)
{
   // Set number of iterations
   _niter = niter;
}

template<class Hist,class Hist2D>
Int_t RooUnfoldIdsT<Hist,Hist2D>::GetNIter() const
{
   // Return number of iterations
   return _niter;
}

template<class Hist,class Hist2D>
void RooUnfoldIdsT<Hist,Hist2D>::SetLambdaM(Double_t lambdaM)
{
   // Set number of iterations
   _lambdaMmin = lambdaM;
}

template<class Hist,class Hist2D>
Double_t RooUnfoldIdsT<Hist,Hist2D>::GetLambdaM() const
{
   // Return number of iterations
   return _lambdaMmin;
}

template<class Hist,class Hist2D>
void RooUnfoldIdsT<Hist,Hist2D>::SetLambdaU(Double_t lambdaU)
{
   // Set number of iterations
   _lambdaUmin = lambdaU;
}

template<class Hist,class Hist2D>
Double_t RooUnfoldIdsT<Hist,Hist2D>::GetLambdaU() const
{
   // Return number of iterations
   return _lambdaUmin;
}

template<class Hist,class Hist2D>
void RooUnfoldIdsT<Hist,Hist2D>::SetLambdaL(Double_t lambdaL)
{
   // Set number of iterations
   _lambdaL = lambdaL;
}

template<class Hist,class Hist2D>
Double_t RooUnfoldIdsT<Hist,Hist2D>::GetLambdaL() const
{
   // Return number of iterations
  return _lambdaL;
}

template<class Hist,class Hist2D>
void RooUnfoldIdsT<Hist,Hist2D>::SetLambdaS(Double_t lambdaS)
{
   // Set number of iterations
   _lambdaS = lambdaS;
}

template<class Hist,class Hist2D>
Double_t RooUnfoldIdsT<Hist,Hist2D>::GetLambdaS() const
{
   // Return number of iterations
   return _lambdaS;
}


template class RooUnfoldIdsT<TH1,TH2>;
ClassImp (RooUnfoldIds);

#ifndef NOROOFIT
template class RooUnfoldIdsT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
ClassImp (RooFitUnfoldIds);
#endif
