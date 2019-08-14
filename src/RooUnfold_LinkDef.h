#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class RooUnfold-;
#pragma link C++ class RooUnfoldBayes+;
#pragma link C++ class RooUnfoldSvd-;
#pragma link C++ class RooUnfoldSvd::SVDUnfold-;
#pragma link C++ class RooUnfoldBinByBin+;
#pragma link C++ class RooUnfoldIds+;
#ifndef NOTUNFOLD
#pragma link C++ class RooUnfoldTUnfold+;
#endif
#pragma link C++ class TUnfold+;
#pragma link C++ class TUnfoldSys+;
#pragma link C++ class TUnfoldBinning+;
#pragma link C++ class TUnfoldDensity+;
#pragma link C++ class TUnfoldBinningXML+;
#pragma link C++ class RooUnfoldResponseT<TH1,TH2>-;
#pragma link C++ class RooUnfoldResponse+;
#ifndef NOROOFIT
#pragma link C++ class RooUnfolding::RooFitHist+;
#pragma link C++ class RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>-;
#pragma link C++ class RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>-;
#pragma link C++ class RooUnfolding::RooFitWrapper<RooAbsReal>+;
#pragma link C++ class RooUnfolding::RooFitWrapper<RooAbsPdf>+;
#pragma link C++ class RooUnfoldFunc+;
#pragma link C++ class RooUnfoldPdf+;
#pragma link C++ class RooUnfoldSpec+;
#pragma link C++ class RooUnfoldInvertT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>+;
#pragma link C++ class RooUnfoldBayesT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>+;
#pragma link C++ class RooUnfoldSvdT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>-;
#pragma link C++ class RooUnfoldSvdT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::SVDUnfold-;
#pragma link C++ class RooUnfoldBinByBinT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>+;
#pragma link C++ class RooUnfoldIdsT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>+;
#pragma link C++ class RooUnfoldTUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>+;
#pragma link C++ class RooFitUnfoldResponse+;
#endif
#pragma link C++ class RooUnfoldErrors+;
#pragma link C++ class RooUnfoldParms+;
#pragma link C++ class RooUnfoldInvert+;
#ifdef HAVE_DAGOSTINI
#pragma link C++ class RooUnfoldDagostini+;
#endif

#endif
