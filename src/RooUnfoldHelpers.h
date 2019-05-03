#ifndef ROOUNFOLDHELPERS_HH
#define ROOUNFOLDHELPERS_HH

#include <TVectorD.h>

namespace RooUnfolding {
  enum Dimension { X, Y, Z };

  template<class Hist> Hist* createHist(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname);  
  template<class Hist2D> Hist2D* createHist(const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, int nbinsy, double ymin, double ymax, const char* yname);
  template<class Hist> Hist* createHist(const TVectorD& vec, const char* name, const char* title, int nbinsx, double xmin, double xmax, const char* xname, bool overflow=false);
  
}

#endif
