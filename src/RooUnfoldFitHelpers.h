#ifndef ROOUNFOLDHELPERS_ROOABSREAL_HH
#define ROOUNFOLDHELPERS_ROOABSREAL_HH

#include "RooUnfoldHelpers.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"

namespace RooUnfolding {
  template<> struct Variable<RooAbsReal> {
    RooRealVar* _var;
    Variable(int nBins,double min,double max,const char* name);
    Variable(RooRealVar* var);
  };
}

#endif
