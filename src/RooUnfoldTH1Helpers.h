#ifndef ROOUNFOLDHELPERS_TH1_HH
#define ROOUNFOLDHELPERS_TH1_HH

#include "RooUnfoldHelpers.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace RooUnfolding {
  template<> struct Variable<TH1> {
    int _nBins;
    double _min;
    double _max;
    Variable(int nBins,double min,double max,const char* name) : _nBins(nBins),_min(min),_max(max){};
  };
  template<> struct Variable<TH2> {
    int _nBins;
    double _min;
    double _max;
    Variable(int nBins,double min,double max,const char* name) : _nBins(nBins),_min(min),_max(max){};
  };
  template<> struct Variable<TH3> {
    int _nBins;
    double _min;
    double _max;
    Variable(int nBins,double min,double max,const char* name) : _nBins(nBins),_min(min),_max(max){};
  };    
}

#endif
