#!/usr/bin/env python
# ==============================================================================
#  File and Version Information:
#       $Id$
#
#  Description:
#       Simple example usage of the RooUnfold package using toy MC.
#
#  Author: Tim Adye <T.J.Adye@rl.ac.uk>
#
# ==============================================================================

# Smearing parameters
bias = -2.5
sigma = 0.3

# Binning
truthBins = 40
recoBins = 50

# Kinematic range
xmin = -10.0
xmax = 10.0

# Number of events
nevents = 100000

# In-/Exclude overflow
overflow = False





def smear(xt):
  from ROOT import gRandom

  xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  #  efficiency
  x= gRandom.Rndm();

  if x>xeff: return None;

  xsmear= gRandom.Gaus(bias,sigma);     #  bias and smear

  return xt+xsmear;

def prepare():

  import ROOT

  # Define the response matrix as a 2D histogram.
  response= ROOT.RooUnfoldResponse (recoBins, xmin, xmax, truthBins, xmin, xmax);
  
  # The histograms that are used to validate the unfolding result.
  sig_theory= ROOT.TH1D ("sig_theory", "Test Truth",    truthBins, xmin, xmax);
  sig_theory_alt= ROOT.TH1D ("sig_theory_alt", "Test Truth",    truthBins, xmin, xmax);

  # The histograms that will be unfolded.
  sigbkg_reco= ROOT.TH1D ("sigbkg_reco", "Test Measured", recoBins, xmin, xmax);
  sigbkg_reco_alt= ROOT.TH1D ("sigbkg_reco_alt", "Test Measured", recoBins, xmin, xmax);
  
  # The histogram containing the background that is needed for the unfolding.
  bkg_reco= ROOT.TH1D ("bkg", "Test Bkg",    recoBins, xmin, xmax);


  #  Train with a Breit-Wigner, mean 0.3 and width 2.5.
  for i in range(nevents):
    x_truth = ROOT.gRandom.BreitWigner (0.3, 2.5);
    x_reco = smear (x_truth);
    if x_reco!=None:
      response.Fill (x_reco, x_truth);
    else:
      response.Miss (x_truth);

    x_truth2 = ROOT.gRandom.BreitWigner (0.3, 2.5);
    x_reco2 = smear(x_truth2)
    sig_theory.Fill(x_truth2)

    if x_reco2!=None: sigbkg_reco.Fill(x_reco2)
      
  #  Test with a Gaussian, mean 0 and width 2.
  for i in range(nevents):
    x_truth = ROOT.gRandom.Gaus (0.0, 2.0)
    x_reco = smear(x_truth);
    sig_theory_alt.Fill(x_truth);
    if x_reco!=None: sigbkg_reco_alt.Fill(x_reco);

  for i in range(nevents):
    x_bkg = ROOT.gRandom.Uniform (-10,10.)
    bkg_reco.Fill(x_bkg)
    x_bkg2 = ROOT.gRandom.Uniform (-10,10.)    
    sigbkg_reco.Fill(x_bkg2)
    sigbkg_reco_alt.Fill(x_bkg2)

  histograms = {"sig_theory":sig_theory,
                "sig_theory_alt":sig_theory_alt,
                "sigbkg_reco":sigbkg_reco,
                "sigbkg_reco_alt":sigbkg_reco_alt,
                "sig_response":response.Hresponse(),
                "sig_reco":response.Hmeasured(),
                "sig_truth":response.Htruth(),
                "bkg_reco":bkg_reco}

  return histograms

    
def makePlots(ws):
  import ROOT

  
  # Get the RooUnfoldPdf(RooAbsPdf) object that contains all the 
  # unfolding results.
  unfoldpdf = ws.pdf("unfold")

  # Get the RooUnfold object.
  unfolding = unfoldpdf.unfolding()

  # Get all the input and output distributions as RooAbsReal objects.
  sig_response = ws.function("response")
  sig_reco = ws.function("sig_reco")
  sig_theory = ws.function("sig_theory")
  # sigbkg_reco = ws.function("sigbkg_reco")
  # sig_truth = ws.function("sig_truth")
  # sig_reco = ws.function("unfold_data_minus_bkg")
  # bkg_reco = ws.function("bkg")

  # Get the kinematic variable.
  obs_truth = ws.var("obs_truth")
  
  # Create a RooPlot object.
  plot_truth = obs_truth.frame()
  nexp = unfoldpdf.expectedEvents(0)

  allVars = ROOT.RooArgList(ws.allVars())

  prefit_all = ROOT.RooFitResult.prefitResult(allVars)


  sig_theory.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sig_theory_graph"))
  unfoldpdf.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("unfold_graph"),ROOT.RooFit.MarkerColor(ROOT.kBlack),ROOT.RooFit.MarkerSize(1),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.DrawOption("P"),ROOT.RooFit.VisualizeError(prefit_all),ROOT.RooFit.Normalization(nexp,ROOT.RooAbsReal.NumEvent),ROOT.RooFit.NormRange("full"))



  #unfoldpdf.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("unfold_graph"))

  canvas_truth = ROOT.TCanvas("unfolded","unfolded")
  plot_truth.Draw()
  leg_truth = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
  leg_truth.SetLineWidth(0)
  leg_truth.SetFillStyle(0)
  leg_truth.AddEntry( plot_truth.findObject("sig_theory_graph"), "Theory Prediction", "l" )
  leg_truth.AddEntry( plot_truth.findObject("unfold_graph"), "Unfolded Signal", "pe" )
  # leg_truth.AddEntry( plot_truth.findObject("unfold_graph_sys"), "Unfolded Signal, Sys. Uncertainty", "f" )  
  leg_truth.Draw()
  canvas_truth.SaveAs("unfolded.pdf")


  # Get the RooUnfold object.

  # Get all the histograms in the form of a RooAbsReal

  # Define the histogram variable as a RooRealVar 
  
  # Create a RooPlot object.



  # obs_truth = ws.var(obsname+"_truth")
  # plot_truth = obs_truth.frame(ROOT.RooFit.Title(obs_truth.GetTitle() +", "+label))

  # plot_truth.SetYTitle("Number of Events")
  # plot_truth.SetMinimum(0)
  # copyBinning(hist_truth,plot_truth)  


  # unfoldpdf = ws.pdf("unfold")
  # unfolding = unfoldpdf.unfolding()

  # sig_truth = ws.function("sig_theory")
  # sig_reco = ws.function("sig_reco")
  # data_minus_bkg_reco = ws.function("unfold_data_minus_bkg")
  # bkg_reco = ws.function("bkg_reco")
  # asm_reco = ws.function("asm_reco")

  # allVars = ROOT.RooArgList(ws.allVars())
  # sysVars = ROOT.RooArgList()

  # obs_var = ws.var("obs_truth")
  
  # for v in allVars:

  #   if "alpha_" in v.GetName():
  #     sysVars.add(v)

  # prefit_all = ROOT.RooFitResult.prefitResult(allVars)
  # prefit_sys = ROOT.RooFitResult.prefitResult(sysVars)


  # nexp = unfoldpdf.expectedEvents(0)
  # unfoldpdf.plotOn(plot_truth,ROOT.RooFit.Invisible(),ROOT.RooFit.DrawOption("P"))
  # unfoldpdf.plotOn(plot_truth,ROOT.RooFit.LineColor(0),ROOT.RooFit.FillColor(ROOT.kGray),ROOT.RooFit.FillStyle(1001),ROOT.RooFit.Name("unfold_graph_sys"),ROOT.RooFit.VisualizeError(prefit_sys),
  #                  ROOT.RooFit.Normalization(nexp,ROOT.RooAbsReal.NumEvent),ROOT.RooFit.NormRange("full"))
  # unfoldpdf.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("unfold_graph"),ROOT.RooFit.MarkerColor(ROOT.kBlack),ROOT.RooFit.MarkerSize(1),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.DrawOption("P"),ROOT.RooFit.VisualizeError(prefit_all),
  #                  ROOT.RooFit.Normalization(nexp,ROOT.RooAbsReal.NumEvent),ROOT.RooFit.NormRange("full"))
  # sig_truth.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.FillColor(ROOT.kRed),ROOT.RooFit.FillStyle(3345),ROOT.RooFit.Name("sig_truth_graph"),
  #                  ROOT.RooFit.Normalization(nexp,ROOT.RooAbsReal.NumEvent),ROOT.RooFit.NormRange("full"))
  # canvas_truth = ROOT.TCanvas("unfolded","unfolded")
  # plot_truth.Draw()
  # leg_truth = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
  # leg_truth.SetLineWidth(0)
  # leg_truth.SetFillStyle(0)
  # leg_truth.AddEntry( plot_truth.findObject("sig_truth_graph"), "Theory Prediction", "l" )
  # leg_truth.AddEntry( plot_truth.findObject("unfold_graph"), "Unfolded Signal", "pe" )
  # leg_truth.AddEntry( plot_truth.findObject("unfold_graph_sys"), "Unfolded Signal, Sys. Uncertainty", "f" )  
  # leg_truth.Draw()
  # canvas_truth.SaveAs(pjoin(outdir,"unfolded.pdf"))
  # canvas_truth.SaveAs(pjoin(outdir,"unfolded.root"))
  # canvas_truth.SaveAs(pjoin(outdir,"unfolded.C"))
  # canvas_truth.SaveAs(pjoin(outdir,"unfolded.png"))
  
  # obs_reco = ws.var(obsname+"_reco")
  # plot_reco = obs_reco.frame(ROOT.RooFit.Title(obs_reco.GetTitle() +", "+label))
  # plot_reco.SetYTitle("Number of Events")
  # plot_reco.SetMinimum(0)
  # copyBinning(hist_reco,plot_reco)
  # data_minus_bkg_reco.plotOn(plot_reco,ROOT.RooFit.Name("data_minus_bkg_reco_graph"),
  #                            ROOT.RooFit.MarkerColor(1),ROOT.RooFit.MarkerSize(1),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.DrawOption("P"),ROOT.RooFit.VisualizeError(prefit_all))
  # sig_reco.plotOn           (plot_reco,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sig_reco_graph"),ROOT.RooFit.Normalization(hist_reco.Integral(),ROOT.RooAbsReal.NumEvent),ROOT.RooFit.NormRange("full"))
  # canvas_reco = ROOT.TCanvas("reco","reco")
  # plot_reco.Draw()
  # leg_reco = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
  # leg_reco.SetLineWidth(0)
  # leg_reco.SetFillStyle(0)
  # leg_reco.AddEntry( plot_reco.findObject("data_minus_bkg_reco_graph"), "Background Subtracted Data", "pe" )
  # leg_reco.AddEntry( plot_reco.findObject("sig_reco_graph"), "Signal Reco Monte Carlo", "l" )
  # leg_reco.Draw()
  # canvas_reco.SaveAs(pjoin(outdir,"reco.pdf"))
  # canvas_reco.SaveAs(pjoin(outdir,"reco.C"))    
  # canvas_reco.SaveAs(pjoin(outdir,"reco.root"))  
  # canvas_reco.SaveAs(pjoin(outdir,"reco.png"))



def algorithm(method):
  import ROOT
  alg = None
  if method == "bayes":
    alg= ROOT.RooUnfolding.kBayes;
  elif method == "bbb":
    alg= ROOT.RooUnfolding.kBinByBin;
  elif method == "inv":
    alg= ROOT.RooUnfolding.kInvert;
  elif method == "svd":
    alg= ROOT.RooUnfolding.kSVD;
  elif method == "root":
    alg= ROOT.RooUnfolding.kTUnfold;
  elif method == "ids":
    alg= ROOT.RooUnfolding.kIDS;
  return alg

def main(args):
  from os.path import exists  
  import ROOT

  histograms = prepare()

  spec = ROOT.RooUnfoldSpec("unfold","unfold",histograms["sig_truth"],"obs_truth",histograms["sig_reco"],"obs_reco",histograms["sig_response"],histograms["bkg_reco"],histograms["sigbkg_reco"],overflow,0.0005,False)
    
  # This line instantiates one of the subclasses specific to the unfolding algorithm.
  if args.regparm:
    pdf = spec.makePdf(algorithm(args.method),args.regparm)
  else:
    pdf = spec.makePdf(algorithm(args.method))
    
  theory = pdf.unfolding().response().makeHistFuncTruth(histograms["sig_theory"])

  # required to avoid python garbage collector messing up the RooDataHists added to gDirectory
  ROOT.gDirectory.Clear()
  
  test_truth = spec.makeHistogram(histograms["sig_theory"])
  
  # Here the first unfolding is performed.
  pdf.unfolding().PrintTable(ROOT.cout, test_truth)

  #unfolded_hist = pdf.unfolding().Hreco()
  
  ws = ROOT.RooWorkspace("workspace","workspace")
  
  getattr(ws,"import")(pdf)
  
  getattr(ws,"import")(theory)
  
  # print(type(histograms["sig_theory"]))
  # print(type(theory))

  # pdf.unfolding().Hmeasured().printHist()
  # pdf.unfolding().Hreco().printHist()


  #getattr(ws,"import")(unfolded_hist)

  #ws.Print("t")

  makePlots(ws)

  pdf.Delete()
  theory.Delete()


  

  
if __name__=="__main__":
  from argparse import ArgumentParser
  parser = ArgumentParser(description="RooUnfold testing script")
  parser.add_argument("method",default="bbb",type=str)
  parser.add_argument("--regparm",type=float)
  main(parser.parse_args())
