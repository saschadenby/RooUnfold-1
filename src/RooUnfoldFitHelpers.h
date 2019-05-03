#ifndef ROOUNFOLDFITHELPERS_HH
#define ROOUNFOLDFITHELPERS_HH

class RooAbsReal;

#include "RooUnfoldHelpers.h"

namespace RooUnfolding {
  void reset(RooAbsReal* r);
  int findBin(const RooAbsReal* h, double x, RooUnfolding::Dimension d);
  double min(const RooAbsReal* hist, RooUnfolding::Dimension d);
  double max(const RooAbsReal* hist, RooUnfolding::Dimension d);
  int sumW2N(const RooAbsReal* hist);
  void add(RooAbsReal* hista, RooAbsReal* histb);
  void projectY(RooAbsReal* _res, RooAbsReal* _tru, bool overflow);
  void projectX(RooAbsReal* _res, RooAbsReal* _mes, bool overflow);
  void subtractProjectX(RooAbsReal* _res, RooAbsReal* _mes, RooAbsReal* _fak, bool overflow);
  int fill(RooAbsReal* hist, double x, double w);
  int fill(RooAbsReal* hist, double x, double y, double w);
  int fill(RooAbsReal* hist, double x, double y, double z, double w);
  RooAbsReal* copy(const RooAbsReal* r, bool reset, const char* name = 0, const char* title = 0);
  int entries(const RooAbsReal* hist);
  int dim(const RooAbsReal* hist);
  int nBins(const RooAbsReal* hist);
  int nBins(const RooAbsReal* hist, RooUnfolding::Dimension d);
  double binCenter(const RooAbsReal*h, int i, RooUnfolding::Dimension d);
  double binWidth(const RooAbsReal*h, int i, RooUnfolding::Dimension d);
  double binHighEdge(const RooAbsReal*h, int i, RooUnfolding::Dimension d);
  double binLowEdge(const RooAbsReal*h, int i, RooUnfolding::Dimension d);
  void binXYZ(const RooAbsReal* tru, int i, int& jx, int& jy, int& jz);
  double binError(const RooAbsReal* h, int i, bool overflow);
  double binContent (const RooAbsReal* h, int i, bool overflow);
  void setBinContent (RooAbsReal* h, int i, double val, bool overflow);
  void setBinContent (RooAbsReal* h, int i, int j, double val, bool overflow);
  RooAbsReal* h2h1d(RooAbsReal* h, int nb);
  RooAbsReal* copyHistogram(const RooAbsReal* h, bool includeOverflow);
  const char* varname(const RooAbsReal* h, Dimension d);
}  

#endif
