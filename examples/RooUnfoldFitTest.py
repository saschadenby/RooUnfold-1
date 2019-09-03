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
bias = -1.5
sigma = 0.3

# Binning
truthBins = 40
recoBins = 50

# Kinematic range
xmin = -10.0
xmax = 10.0

# Number of events
nevents = 100000




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

  # Get all the input and output distributions as RooAbsReal objects.
  sig_theory = ws.function("truth_hist")
  data_minus_bkg_reco = ws.function("unfold_data_minus_bkg")

  # Get the kinematic variable.
  obs_truth = ws.var("obs_truth")
  
  # Create a RooPlot object.
  plot_truth = obs_truth.frame()
  nexp = unfoldpdf.expectedEvents(0)

  # Get all the variables of the histograms.
  allVars = ROOT.RooArgList(ws.allVars())

  # Get the errors.
  prefit_all = ROOT.RooFitResult.prefitResult(allVars)

  # Plot the theory prediction.
  sig_theory.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sig_theory_graph"))

  # Plot the unfolded result.
  unfoldpdf.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.Name("unfold_graph"),ROOT.RooFit.MarkerColor(ROOT.kBlack),ROOT.RooFit.MarkerSize(1),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.DrawOption("P"),ROOT.RooFit.VisualizeError(prefit_all),ROOT.RooFit.Normalization(nexp,ROOT.RooAbsReal.NumEvent),ROOT.RooFit.NormRange("full"))

  # Save the plot to a pdf file.
  canvas_truth = ROOT.TCanvas("unfolded","unfolded")
  plot_truth.Draw()
  leg_truth = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
  leg_truth.SetLineWidth(0)
  leg_truth.SetFillStyle(0)
  leg_truth.AddEntry( plot_truth.findObject("sig_theory_graph"), "Theory Prediction", "l" )
  leg_truth.AddEntry( plot_truth.findObject("unfold_graph"), "Unfolded Data", "pe" )
  leg_truth.AddEntry( plot_truth.findObject("data_minus_bkg_reco_graph"), "Reconstructed Data", "pe" )
  leg_truth.Draw()
  canvas_truth.SaveAs("unfolded.pdf")



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

  # Prepare the histograms and response matrix.
  histograms = prepare()

  # Initiate the unfolding setup.
  spec = ROOT.RooUnfoldSpec("unfold","unfold",histograms["sig_truth"],"obs_truth",histograms["sig_reco"],"obs_reco",histograms["sig_response"],histograms["bkg_reco"],histograms["sigbkg_reco"],args.overflow,0.0005,False)
    
  # Create the object that will perform the unfolding. Pass a
  # regularization parameter if given in the command line.
  if args.regparm:
    pdf = spec.makePdf(algorithm(args.method),args.regparm)
  else:
    pdf = spec.makePdf(algorithm(args.method))
  
  # required to avoid python garbage collector messing up the RooDataHists added to gDirectory
  ROOT.gDirectory.Clear()
  
  # Create a RooFitHist object as input for the unfolding
  # results printing.
  test_truth = spec.makeHistogram(histograms["sig_theory"])
  
  # Do the unfolding and print the unfolding results.
  pdf.unfolding().PrintTable(ROOT.cout, test_truth)

  if not args.plot:
    return

  # Create a workspace.
  ws = ROOT.RooWorkspace("workspace","workspace")
  
  # Save the unfolding to the workspace.
  getattr(ws,"import")(pdf)
  
  # Print the saved workspace.
  ws.Print("t")

  # Make a plot the unfolding result and theory prediction.
  makePlots(ws)


  

  
if __name__=="__main__":
  from argparse import ArgumentParser
  parser = ArgumentParser(description="RooUnfold testing script")
  parser.add_argument("method",default="bbb",type=str)
  parser.add_argument("--plot",default=False,type=bool)
  parser.add_argument("--regparm",type=float)
  parser.add_argument("--overflow",default=True,type=bool)
  main(parser.parse_args())
