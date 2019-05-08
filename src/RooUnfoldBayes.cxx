//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Bayesian unfolding.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p>Links to the RooUnfoldBayesImpl class which uses Bayesian unfolding to reconstruct the truth distribution.</p>
<p>Works for 2 and 3 dimensional distributions
<p>Returned errors can be either as a diagonal matrix or as a full matrix of covariances
<p>Regularisation parameter sets the number of iterations used in the unfolding (default=4)
<p>Is able to account for bin migration and smearing
<p>Can unfold if test and measured distributions have different binning.
<p>Returns covariance matrices with conditions approximately that of the machine precision. This occasionally leads to very large chi squared values
END_HTML */

/////////////////////////////////////////////////////////////

//#define OLDERRS   // restore old (incorrect) error calculation
//#define OLDERRS2  // restore old (incorrect) systematic error calculation
//#define OLDMULT   // restore old (slower) matrix multiplications

#include "RooUnfoldBayes.h"

#include <iostream>
#include <iomanip>
#include <math.h>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"

#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"
#include "RooUnfoldFitHelpers.h"
#include "RooUnfoldResponse.h"

using namespace RooUnfolding;

using std::min;
using std::cerr;
using std::endl;
using std::cout;
using std::setw;
using std::left;
using std::right;

ClassImp (RooUnfoldBayes);

template<class Hist,class Hist2D>
RooUnfoldBayesT<Hist,Hist2D>::RooUnfoldBayesT (const RooUnfoldBayesT<Hist,Hist2D>& rhs)
  : RooUnfoldT<Hist,Hist2D> (rhs)
{
  // Copy constructor.
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D>
RooUnfoldBayesT<Hist,Hist2D>::RooUnfoldBayesT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, int niter, Bool_t smoothit,
                                const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title), _niter(niter), _smoothit(smoothit)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  // The regularisation parameter is niter (number of iterations).
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldBayesT<Hist,Hist2D>* RooUnfoldBayesT<Hist,Hist2D>::Clone (const char* newname) const
{
  // Creates a copy of the RooUnfoldBayesT<Hist,Hist2D> object
  RooUnfoldBayesT<Hist,Hist2D>* unfold= new RooUnfoldBayesT<Hist,Hist2D>(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::Init()
{
  this->_nc= this->_ne= 0;
  this->_nbartrue= this->_N0C= 0.0;
  GetSettings();
}

template<class Hist,class Hist2D> void 
RooUnfoldBayesT<Hist,Hist2D>::Reset()
{
  Init();
  RooUnfoldT<Hist,Hist2D>::Reset();
}

template<class Hist,class Hist2D> void 
RooUnfoldBayesT<Hist,Hist2D>::Assign (const RooUnfoldBayesT<Hist,Hist2D>& rhs)
{
  RooUnfoldT<Hist,Hist2D>::Assign (rhs);
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::CopyData (const RooUnfoldBayesT<Hist,Hist2D>& rhs)
{
  this->_niter=    rhs._niter;
  this->_smoothit= rhs._smoothit;
}

template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::Unfold()
{
  this->setup();
  if (this->verbose() >= 2) {
    Print();
    printMatrix(this->_Nji,"RooUnfoldBayesT<Hist,Hist2D> response matrix (Nji)");
  }
  if (this->verbose() >= 1) cout << "Now unfolding..." << endl;
  unfold();
  if (this->verbose() >= 2) Print();
  this->_rec.ResizeTo(this->_nc);
  this->_rec = this->_nbarCi;
  this->_rec.ResizeTo(this->_nt);  // drop fakes in final bin
  this->_unfolded= true;
  this->_haveCov=  false;
}

template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::GetCov()
{
  getCovariance();
  this->_cov.ResizeTo (this->_nt, this->_nt);  // drop fakes in final bin
  this->_haveCov= true;
}

template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::GetSettings()
{
    this->_minparm=1;
    this->_maxparm=15;
    this->_stepsizeparm=1;
    this->_defaultparm=4;
}

//-------------------------------------------------------------------------
template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::setup()
{
  this->_nc = this->_nt;
  this->_ne = this->_nm;

  this->_nEstj.ResizeTo(this->_ne);
  this->_nEstj= this->Vmeasured();

  this->_nCi.ResizeTo(this->_nt);
  this->_nCi= this->_res->Vtruth();

  this->_Nji.ResizeTo(this->_ne,this->_nt);
  h2m (this->_res->Hresponse(), this->_Nji, this->_overflow);   // don't normalise, which is what this->_res->Mresponse() would give us

  if (this->_res->FakeEntries()) {
    TVectorD fakes= this->_res->Vfakes();
    double nfakes= fakes.Sum();
    if (this->verbose()>=0) cout << "Add truth bin for " << nfakes << " fakes" << endl;
    this->_nc++;
    this->_nCi.ResizeTo(this->_nc);
    this->_nCi[this->_nc-1]= nfakes;
    this->_Nji.ResizeTo(this->_ne,this->_nc);
    for (int i= 0; i<this->_nm; i++) this->_Nji(i,this->_nc-1)= fakes[i];
  }

  this->_nbarCi.ResizeTo(this->_nc);
  this->_efficiencyCi.ResizeTo(this->_nc);
  this->_Mij.ResizeTo(this->_nc,this->_ne);
  this->_P0C.ResizeTo(this->_nc);
  this->_UjInv.ResizeTo(this->_ne);
#ifndef OLDERRS
  if (this->_dosys!=2) this->_dnCidnEj.ResizeTo(this->_nc,this->_ne);
#endif
  if (this->_dosys)    this->_dnCidPjk.ResizeTo(this->_nc,this->_ne*this->_nc);

  // Initial distribution
  this->_N0C= this->_nCi.Sum();
  if (this->_N0C!=0.0) {
    this->_P0C= this->_nCi;
    this->_P0C *= 1.0/this->_N0C;
  }
}

//-------------------------------------------------------------------------
template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::unfold()
{
  // Calculate the unfolding matrix.
  // _niter = number of iterations to perform (3 by default).
  // _smoothit = smooth the matrix in between iterations (default false).

  TMatrixD PEjCi(this->_ne,this->_nc), PEjCiEff(this->_ne,this->_nc);
  for (int i = 0 ; i < this->_nc ; i++) {
    if (this->_nCi[i] <= 0.0) { this->_efficiencyCi[i] = 0.0; continue; }
    double eff = 0.0;
    for (int j = 0 ; j < this->_ne ; j++) {
      double response = this->_Nji(j,i) / this->_nCi[i];
      PEjCi(j,i) = PEjCiEff(j,i) = response;  // efficiency of detecting the cause Ci in Effect Ej
      eff += response;
    }
    this->_efficiencyCi[i] = eff;
    double effinv = eff > 0.0 ? 1.0/eff : 0.0;   // reset PEjCiEff if eff=0
    for (int j = 0 ; j < this->_ne ; j++) PEjCiEff(j,i) *= effinv;
  }

  TVectorD PbarCi(this->_nc);

  for (int kiter = 0 ; kiter < this->_niter; kiter++) {

    if (this->verbose()>=1) cout << "Iteration : " << kiter << endl;

    // update prior from previous iteration
    if (kiter>0) {
      this->_P0C = PbarCi;
      this->_N0C = this->_nbartrue;
    }

    for (int j = 0 ; j < this->_ne ; j++) {
      double Uj = 0.0;
      for (int i = 0 ; i < this->_nc ; i++)
        Uj += PEjCi(j,i) * this->_P0C[i];
      this->_UjInv[j] = Uj > 0.0 ? 1.0/Uj : 0.0;
    }

    // Unfolding matrix M
    this->_nbartrue = 0.0;
    for (int i = 0 ; i < this->_nc ; i++) {
      double nbarC = 0.0;
      for (int j = 0 ; j < this->_ne ; j++) {
        double Mij = this->_UjInv[j] * PEjCiEff(j,i) * this->_P0C[i];
        this->_Mij(i,j) = Mij;
        nbarC += Mij * this->_nEstj[j];
      }
      this->_nbarCi[i] = nbarC;
      this->_nbartrue += nbarC;  // best estimate of true number of events
    }

    // new estimate of true distribution
    PbarCi= this->_nbarCi;
    PbarCi *= 1.0/this->_nbartrue;

#ifndef OLDERRS
    if (this->_dosys!=2) {
      if (kiter <= 0) {
        this->_dnCidnEj= this->_Mij;
      } else {
#ifndef OLDMULT
        TVectorD en(this->_nc), nr(this->_nc);
        for (int i = 0 ; i < this->_nc ; i++) {
          if (this->_P0C[i]<=0.0) continue;
          double ni= 1.0/(this->_N0C*this->_P0C[i]);
          en[i]= -ni*this->_efficiencyCi[i];
          nr[i]=  ni*this->_nbarCi[i];
        }
        TMatrixD M1= this->_dnCidnEj;
        M1.NormByColumn(nr,"M");
        TMatrixD M2 (TMatrixD::kTransposed, this->_Mij);
        M2.NormByColumn(this->_nEstj,"M");
        M2.NormByRow(en,"M");
        TMatrixD M3 (M2, TMatrixD::kMult, this->_dnCidnEj);
        this->_dnCidnEj.Mult (this->_Mij, M3);
        this->_dnCidnEj += this->_Mij;
        this->_dnCidnEj += M1;
#else /* OLDMULT */
        TVectorD ksum(this->_ne);
        for (int j = 0 ; j < this->_ne ; j++) {
          for (int k = 0 ; k < this->_ne ; k++) {
            double sum = 0.0;
            for (int l = 0 ; l < this->_nc ; l++) {
              if (this->_P0C[l]>0.0) sum += this->_efficiencyCi[l]*this->_Mij(l,k)*this->_dnCidnEj(l,j)/this->_P0C[l];
            }
            ksum[k]= sum;
          }
          for (int i = 0 ; i < this->_nc ; i++) {
            double dsum = this->_P0C[i]>0 ? this->_dnCidnEj(i,j)*this->_nbarCi[i]/this->_P0C[i] : 0.0;
            for (int k = 0 ; k < this->_ne ; k++) {
              dsum -= this->_Mij(i,k)*this->_nEstj[k]*ksum[k];
            }
            // update dnCidnEj. Note that we can do this in-place due to the ordering of the accesses.
            this->_dnCidnEj(i,j) = this->_Mij(i,j) + dsum/this->_N0C;
          }
        }
#endif
      }
    }
#endif

    if (this->_dosys) {
#ifndef OLDERRS2
      if (kiter > 0) {
        TVectorD mbyu(this->_ne);
        for (int j = 0 ; j < this->_ne ; j++) {
          mbyu[j]= this->_UjInv[j]*this->_nEstj[j]/this->_N0C;
        }
        TMatrixD A= this->_Mij;
        A.NormByRow (mbyu, "M");
        TMatrixD B(A, TMatrixD::kMult, PEjCi);
        TMatrixD dnCidPjkUpd (B, TMatrixD::kMult, this->_dnCidPjk);
        int nec= this->_ne*this->_nc;
        for (int i = 0 ; i < this->_nc ; i++) {
          if (this->_P0C[i]<=0.0) continue;  // skip loop: dnCidPjkUpd(i,jk) will also be 0 because this->_Mij(i,j) will be 0
          double r= PbarCi[i]/this->_P0C[i];
          for (int jk= 0; jk<nec; jk++)
            this->_dnCidPjk(i,jk)= r*this->_dnCidPjk(i,jk) - dnCidPjkUpd(i,jk);
        }
      }
#else  /* OLDERRS2 */
      if (kiter == this->_niter-1)   // used to only calculate this->_dnCidPjk for the final iteration
#endif
      for (int j = 0 ; j < this->_ne ; j++) {
        if (this->_UjInv[j]==0.0) continue;
        double mbyu= this->_UjInv[j]*this->_nEstj[j];
        int j0= j*this->_nc;
        for (int i = 0 ; i < this->_nc ; i++) {
          double b= -mbyu * this->_Mij(i,j);
          for (int k = 0 ; k < this->_nc ; k++) this->_dnCidPjk(i,j0+k) += b*this->_P0C[k];
          if (this->_efficiencyCi[i]!=0.0)
            this->_dnCidPjk(i,j0+i) += (this->_P0C[i]*mbyu - this->_nbarCi[i]) / this->_efficiencyCi[i];
        }
      }
    }

    // no need to smooth the last iteraction
    if (this->_smoothit && kiter < (this->_niter-1)) smooth(PbarCi);

    // Chi2 based on Poisson errors
    double chi2 = getChi2(PbarCi, this->_P0C, this->_nbartrue);
    if (this->verbose()>=1) cout << "Chi^2 of change " << chi2 << endl;

    // and repeat
  }
}

//-------------------------------------------------------------------------
template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::getCovariance()
{
  if (this->_dosys!=2) {
    if (this->verbose()>=1) cout << "Calculating covariances due to number of measured events" << endl;

    // Create the covariance matrix of result from that of the measured distribution
    this->_cov.ResizeTo (this->_nc, this->_nc);
#ifdef OLDERRS
    const TMatrixD& Dprop= this->_Mij;
#else
    const TMatrixD& Dprop= this->_dnCidnEj;
#endif
    if (this->_haveCovMes) {
      ABAT (Dprop, this->GetMeasuredCov(), this->_cov);
    } else {
      TVectorD v= this->Emeasured();
      v.Sqr();
      ABAT (Dprop, v, this->_cov);
    }
  }

  if (this->_dosys) {
    if (this->verbose()>=1) cout << "Calculating covariance due to unfolding matrix..." << endl;

    const TMatrixD& Eres= this->_res->Eresponse();
    TVectorD Vjk(this->_ne*this->_nc);           // vec(Var(j,k))
    for (int j = 0 ; j < this->_ne ; j++) {
      int j0= j*this->_nc;
      for (int i = 0 ; i < this->_nc ; i++) {
        double e= Eres(j,i);
        Vjk[j0+i]= e*e;
      }
    }

    if (this->_dosys!=2) {
      TMatrixD covres(this->_nc,this->_nc);
      ABAT (this->_dnCidPjk, Vjk, covres);
      this->_cov += covres;
    } else {
      this->_cov.ResizeTo (this->_nc, this->_nc);
      ABAT (this->_dnCidPjk, Vjk, this->_cov);
    }
  }
}

//-------------------------------------------------------------------------
template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::smooth(TVectorD& PbarCi) const
{
  // Smooth unfolding distribution. PbarCi is the array of proababilities
  // to be smoothed PbarCi; nevts is the numbers of events
  // (needed to calculate suitable errors for the smearing).
  // PbarCi is returned with the smoothed distribution.

  if (this->_res->GetDimensionTruth() != 1) {
    cerr << "Smoothing only implemented for 1-D distributions" << endl;
    return;
  }
  if (this->verbose()>=1) cout << "Smoothing." << endl;
  TH1::SmoothArray (this->_nc, PbarCi.GetMatrixArray(), 1);
  return;
}

//-------------------------------------------------------------------------
template<class Hist,class Hist2D> double 
RooUnfoldBayesT<Hist,Hist2D>::getChi2(const TVectorD& prob1,
                                 const TVectorD& prob2,
                                 double nevents) const
{
  // calculate the chi^2. prob1 and prob2 are the probabilities
  // and nevents is the number of events used to calculate the probabilities
  double chi2= 0.0;
  int n= prob1.GetNrows();
  if (this->verbose()>=2) cout << "chi2 " << n << " " << nevents << endl;
  for (int i = 0 ; i < n ; i++) {
    double psum  = (prob1[i] + prob2[i])*nevents;
    double pdiff = (prob1[i] - prob2[i])*nevents;
    if (psum > 1.0) {
      chi2 = chi2 + (pdiff*pdiff)/psum;
    } else {
      chi2 = chi2 + (pdiff*pdiff);
    }
  }
  return(chi2);
}

//-------------------------------------------------------------------------
template<class Hist,class Hist2D> void
RooUnfoldBayesT<Hist,Hist2D>::Print(Option_t* option) const
{
  RooUnfoldT<Hist,Hist2D>::Print (option);
  if (this->_nc<=0 || this->_ne<=0) return;

  // Print out some useful info of progress so far

  cout << "-------------------------------------------" << endl;
  cout << "Unfolding Algorithm" << endl;
  cout << "Generated (Training):" << endl;
  cout << "  Total Number of bins   : " << this->_nc << endl;
  cout << "  Total Number of events : " << this->_nCi.Sum() << endl;

  cout << "Measured (Training):" << endl;
  cout << "  Total Number of bins   : " << this->_ne << endl;

  cout << "Input (for unfolding):" << endl;
  cout << "  Total Number of events : " << this->_nEstj.Sum() << endl;

  cout << "Output (unfolded):" << endl;
  cout << "  Total Number of events : " << this->_nbarCi.Sum() <<endl;

  cout << "-------------------------------------------\n" << endl;

  if (((this->_nEstj.Sum())!=0) || ((this->_nCi.Sum())!=0)) {
    int iend = min(this->_nCi.GetNrows(),this->_nEstj.GetNrows());
    cout << "    \tTrain \tTest\tUnfolded"<< endl;
    cout << "Bin \tTruth \tInput\tOutput"<< endl;
    for (int i=0; i < iend ; i++) {
      if ((this->_nCi[i] == 0) && (this->_nEstj[i] == 0) &&
          (this->_nEstj[i] == 0) && (this->_nbarCi[i]==0)) continue;
      cout << i << "\t" << this->_nCi[i]                                      \
           << "\t " << this->_nEstj[i] << "\t " << this->_nbarCi[i] << endl;
    }

    // if the number of bins is different
    if (this->_nCi.GetNrows() > this->_nEstj.GetNrows() ) {
      for (int i=iend; i < this->_nCi.GetNrows() ; i++) {
        cout << i << "\t " << this->_nCi[i] << endl;
      }
    }

    cout << "--------------------------------------------------------" << endl;
    cout << " \t" << (this->_nCi.Sum())
         << "\t " << (this->_nEstj.Sum()) << "\t " << (this->_nbarCi.Sum()) << endl;
    cout << "--------------------------------------------------------\n" << endl;
  }
}

template<class Hist,class Hist2D>
RooUnfoldBayesT<Hist,Hist2D>::RooUnfoldBayesT()
  : RooUnfoldT<Hist,Hist2D>()
{
  // Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldBayesT<Hist,Hist2D>::RooUnfoldBayesT (const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldBayesT<Hist,Hist2D>::RooUnfoldBayesT (const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldBayesT<Hist,Hist2D>& RooUnfoldBayesT<Hist,Hist2D>::operator= (const RooUnfoldBayesT<Hist,Hist2D>& rhs)
{
  // Assignment operator for copying RooUnfoldBayes settings.
  Assign(rhs);
  return *this;
}


template<class Hist,class Hist2D>
void RooUnfoldBayesT<Hist,Hist2D>::SetIterations (int niter)
{
  // Set regularisation parameter (number of iterations)
  this->_niter= niter;
}

template<class Hist,class Hist2D>
void RooUnfoldBayesT<Hist,Hist2D>::SetSmoothing (Bool_t smoothit)
{
  // Enable smoothing
  this->_smoothit= smoothit;
}

template<class Hist,class Hist2D>
int RooUnfoldBayesT<Hist,Hist2D>::GetIterations() const
{
  // Return regularisation parameter (number of iterations)
  return this->_niter;
}

template<class Hist,class Hist2D>
int RooUnfoldBayesT<Hist,Hist2D>::GetSmoothing()  const
{
  // Return smoothing setting
  return this->_smoothit;
}

template<class Hist,class Hist2D>
const TMatrixD& RooUnfoldBayesT<Hist,Hist2D>::UnfoldingMatrix() const
{
  // Access unfolding matrix (Mij)
  return this->_Mij;
}

template<class Hist,class Hist2D>
void  RooUnfoldBayesT<Hist,Hist2D>::SetRegParm (double parm)
{
  // Set regularisation parameter (number of iterations)
  SetIterations(int(parm+0.5));
}

template<class Hist,class Hist2D>
double RooUnfoldBayesT<Hist,Hist2D>::GetRegParm() const
{
  // Return regularisation parameter (number of iterations)
  return GetIterations();
}

template class RooUnfoldBayesT<TH1,TH2>;
ClassImp (RooUnfoldBayes);

template class RooUnfoldBayesT<RooAbsReal,RooAbsReal>;
ClassImp (RooFitUnfoldBayes);
