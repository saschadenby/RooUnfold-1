#ifndef ROOUNFOLDHELPERS_ROOABSREAL_HH
#define ROOUNFOLDHELPERS_ROOABSREAL_HH

#include "RooUnfoldHelpers.h"

#ifndef NOROOFIT

#include "RooRealVar.h"
#include "RooAbsData.h"

class RooHistFunc;

namespace RooUnfolding {
  class RooFitHist { 
  public:
    RooFitHist();
    RooFitHist(const RooFitHist* h);
    RooFitHist(RooAbsReal* f, const std::vector<RooRealVar*>& obs);
    RooFitHist(RooAbsReal* f, RooRealVar* obs);
    RooFitHist(RooAbsReal* f, RooRealVar* obs1, RooRealVar* obs2);
    RooFitHist(RooAbsReal* f, const std::vector<RooRealVar*>& obs, const std::vector<RooRealVar*>& nps);
    RooFitHist(RooAbsReal* f, RooRealVar* obs, const std::vector<RooRealVar*>& nps);
    RooFitHist(RooAbsReal* f, RooRealVar* obs1, RooRealVar* obs2, const std::vector<RooRealVar*>& nps);
    RooFitHist(RooDataHist* f, const std::vector<RooRealVar*>& obs);
    RooFitHist(RooDataHist* f, RooRealVar* obs);
    RooFitHist(RooDataHist* f, RooRealVar* obs1, RooRealVar* obs2);
    
    virtual const char* name() const;
    virtual const char* title() const;
    virtual RooRealVar* obs(size_t) const;
    virtual size_t dim() const;

    virtual double error() const;
    virtual double value() const;
    virtual bool weighted() const;
    virtual int bin() const;
    virtual const std::vector<RooRealVar*>& nps() const;

    virtual void printHist() const;
    virtual RooAbsReal* func() const;
    virtual bool checkValidity() const;
    virtual void replace(const RooAbsCollection& newlist);
  protected:
    RooAbsReal* _func;
    std::vector<RooRealVar*> _obs;
    std::vector<RooRealVar*> _gamma;
    ClassDef(RooFitHist,1)
  };

  template<> struct Variable<RooFitHist> {
    RooRealVar* _var;
    Variable(int nBins,double min,double max,const char* name);
    Variable(RooRealVar* var);
  };
  
  const RooArgSet* getObservables(const RooHistFunc* f);
  void setGammaUncertainties(RooWorkspace* ws);
  RooHistFunc* makeHistFunc(RooDataHist* dhist, const std::vector<RooRealVar*>& obs);  
  RooRealVar* findLeafServer(RooAbsArg* rr, const char* name);
  void importToWorkspace(RooWorkspace* ws, RooAbsReal* object);
  void importToWorkspace(RooWorkspace* ws, RooAbsData* object);
  RooArgSet allVars(RooWorkspace* ws, const char* pattern);
  std::vector<RooAbsArg*> matchingObjects(const RooAbsCollection* c, const char* pattern);
  void printClients(const RooAbsArg* obj);
  void printServers(const RooAbsArg* obj);
}

#endif

#endif
