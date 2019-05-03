#ifndef ROOUNFOLDTH1HELPERS_HH
#define ROOUNFOLDTH1HELPERS_HH

class TH1;
class TH2;

#include "RooUnfoldHelpers.h"

namespace RooUnfolding {
  void reset(TH1* h);
  int findBin(const TH1* h, double x, RooUnfolding::Dimension d);
  double min(const TH1* hist, RooUnfolding::Dimension d);
  double max(const TH1* hist, RooUnfolding::Dimension d);
  int sumW2N(const TH1* hist);
  int entries(const TH1* hist);
  int dim(const TH1* hist);
  int nBins(const TH1* hist);
  int nBins(const TH1* hist, RooUnfolding::Dimension d);
  int bin(const TH1* h, int i, bool overflow);
  int bin(const TH1* h, int i, int j, bool overflow);
  double binCenter(const TH1*h, int i, RooUnfolding::Dimension d);
  double binWidth(const TH1*h, int i, RooUnfolding::Dimension d);
  double binHighEdge(const TH1*h, int i, RooUnfolding::Dimension d);
  double binLowEdge(const TH1*h, int i, RooUnfolding::Dimension d);
  void add(TH1* hista, const TH1* histb);
  void projectY(TH2* _res, TH1* _tru, bool overflow);
  void projectX(TH2* _res, TH1* _mes, bool overflow);
  void subtractProjectX(TH2* _res, TH1* _mes, TH1* _fak, bool overflow);
  int fill(TH1* hist, double x, double w);
  int fill(TH1* hist, double x, double y, double w);
  int fill(TH1* hist, double x, double y, double z, double w);
  int fill(TH2* hist, double x, double y, double w);
  TH1* copy(const TH1* orighist, bool reset, const char* name = 0, const char* title = 0);
  TH2* copy(const TH2* orighist, bool reset, const char* name = 0, const char* title = 0);
  void binXYZ(const TH1* tru, int i, int& jx, int& jy, int& jz);
  double binError(const TH1* h, int i, bool overflow);
  double binContent (const TH1* h, int i, bool overflow);
  void setBinContent (TH1* h, int i, double val, bool overflow);
  void setBinContent (TH1* h, int i, int j, double val, bool overflow);
  TH1* h2h1d(TH1* h, int nb);  
  TH1* h2h1d(TH2* h, int nb);
  TH2* copyHistogram(const TH2* h, bool includeOverflow);
  const char* varname(const TH1* h, Dimension d);
}
#endif
