#!/usr/bin/env python
#1;95;0c ==============================================================================
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
truthXMin = -10.0
truthXMax = 10.0

recoXMin = -10.0
recoXMax = 10.0

# Number of events
nevents = 100000


def smear(xt):
  from ROOT import gRandom

  xeff= 0.3 + (1.0-0.3)/(truthXMax - truthXMin)*(xt-truthXMin);  #  efficiency
  x= gRandom.Rndm();

  if x>xeff: return None;

  xsmear= gRandom.Gaus(bias,sigma);     #  bias and smear

  return xt+xsmear;

def plot_hists(theory, reco, data, unfolded):
  
  print('Plotting unfolding results...')

  import ROOT

  canvas_truth = ROOT.TCanvas("input","input")
  ROOT.gStyle.SetOptStat(0)
  theory.SetLineColor(ROOT.kRed)
  theory.Draw()
  reco.SetLineStyle(7)
  reco.Draw("SAME")
  unfolded.SetLineStyle(7)
  unfolded.SetLineColor(ROOT.kGreen)
  unfolded.Draw("SAME")
  data.SetMarkerStyle(20)
  data.SetMarkerSize(1)
  data.Draw("SAMEP0")
  leg_truth = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
  leg_truth.SetLineWidth(0)
  leg_truth.SetFillStyle(0)
  leg_truth.AddEntry( theory, "Truth MC", "l" )
  leg_truth.AddEntry( reco, "Reco MC", "l" )
  leg_truth.AddEntry( data, "Measured Data", "p" )
  leg_truth.AddEntry( unfolded, "Unfolded Data", "l" )
  leg_truth.Draw()
  canvas_truth.SaveAs("result.pdf")


def prepare_bimodal():

  import ROOT

  # Define the response matrix as a 2D histogram.
  response= ROOT.RooUnfoldResponse (recoBins, recoXMin, recoXMax, truthBins, truthXMin, truthXMax);  

  truth = ROOT.TH1D ("Truth", "Truth", truthBins, truthXMin, truthXMax);
  reco = ROOT.TH1D ("Reco", "Reco", recoBins, recoXMin, recoXMax);
  data = ROOT.TH1D ("Data", "Data", recoBins, recoXMin, recoXMax);

  for i in range(nevents):

    peak_loc = 0.7
    if (i < nevents/2):
      peak_loc = 0.3

    x_truth = ROOT.gRandom.Gaus(truthXMin + (truthXMax-truthXMin)*peak_loc, (truthXMax-truthXMin)*0.1)
    x_reco = smear(x_truth)

    if x_reco!=None:
      response.Fill (x_reco, x_truth);
    else:
      response.Miss (x_truth);

    x_truth = ROOT.gRandom.Gaus(truthXMin + (truthXMax-truthXMin)*peak_loc, (truthXMax-truthXMin)*0.1)
    x_reco = smear(x_truth)
 
    if x_reco!=None: reco.Fill(x_reco)
      
    truth.Fill(x_truth)
     
   
  for i in range(reco.GetNbinsX() + 2):
    x_data = ROOT.gRandom.Poisson(reco.GetBinContent(i))
    data.SetBinContent(i, x_data)


  histograms = {"sig_truth_train":response.Htruth(),
                "sig_reco_train":response.Hmeasured(),
                "sig_truth_test":truth,
                "sigbkg_reco_test":reco,
                "data":data,
                "response":response.Hresponse(),
                "bkg_reco":0
              } 

  return histograms

def prepare_falling():

  import ROOT

  # Define the response matrix as a 2D histogram.
  response= ROOT.RooUnfoldResponse (recoBins, recoXMin, recoXMax, truthBins, truthXMin, truthXMax);  

  truth = ROOT.TH1D ("Truth", "Truth", truthBins, truthXMin, truthXMax);
  reco = ROOT.TH1D ("Reco", "Reco", recoBins, recoXMin, recoXMax);
  data = ROOT.TH1D ("Data", "Data", recoBins, recoXMin, recoXMax);

  for i in range(nevents):


    #x_truth = - ROOT.TMath.Log( 1 - ROOT.gRandom.Uniform( 1 - ROOT.TMath.Exp(-1), xmax - xmax * ROOT.TMath.Exp(-10) )) - 10
    x_truth = ROOT.gRandom.Exp((truthXMax - truthXMin)/4) - (truthXMax - truthXMin)/2 + 2
    x_reco = smear(x_truth)

    if x_reco!=None:
      response.Fill (x_reco, x_truth);
    else:
      response.Miss (x_truth);

    x_truth = ROOT.gRandom.Exp((truthXMax - truthXMin)/4) - (truthXMax - truthXMin)/2 + 2
    x_reco = smear(x_truth)

    truth.Fill(x_truth)

    if x_reco!=None: reco.Fill(x_reco)
    
  for i in range(reco.GetNbinsX() + 2):
    x_data = ROOT.gRandom.Poisson(reco.GetBinContent(i))
    data.SetBinContent(i, x_data)


  histograms = {"sig_truth_train":response.Htruth(),
                "sig_reco_train":response.Hmeasured(),
                "sig_truth_test":truth,
                "sigbkg_reco_test":reco,
                "data":data,
                "response":response.Hresponse(),
                "bkg_reco":0
              } 


  return histograms
      

def prepare():

  import ROOT

  # Define the response matrix as a 2D histogram.
  response= ROOT.RooUnfoldResponse (recoBins, recoXMin, recoXMax, truthBins, truthXMin, truthXMax);  

  # The histograms that are used to validate the unfolding result.
  sig_theory= ROOT.TH1D ("sig_theory", "Test Truth",    truthBins, truthXMin, truthXMax);
  sig_theory_alt= ROOT.TH1D ("sig_theory_alt", "Test Truth",    truthBins, truthXMin, truthXMax);

  # The histograms that will be unfolded.
  sigbkg_reco= ROOT.TH1D ("sigbkg_reco", "Test Measured", recoBins, recoXMin, recoXMax);
  sigbkg_reco_alt= ROOT.TH1D ("sigbkg_reco_alt", "Test Measured", recoBins, recoXMin, recoXMax);
  
  # The histogram containing the background that is needed for the unfolding.
  bkg_reco= ROOT.TH1D ("bkg", "Test Bkg",    recoBins, recoXMin, recoXMax);


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
  elif method == "gp":
    alg= ROOT.RooUnfolding.kGP;
  return alg

def main(args):
  from os.path import exists  
  import ROOT

  # Prepare the histograms and response matrix.
  #histograms = prepare()
  #histograms = prepare_bimodal()
  histograms = prepare_falling()

  # Initiate the unfolding setup.
  spec = ROOT.RooUnfoldSpec("unfold","unfold",histograms["sig_truth_train"],"obs_truth",histograms["sig_reco_train"],"obs_reco",histograms["response"],histograms["bkg_reco"],histograms["data"],args.overflow,0.0005,False)
    
  # Create the object that will perform the unfolding. Pass a
  # regularization parameter if given in the command line.
  if args.regparm:
    func = spec.makeFunc(algorithm(args.method),args.regparm)
  else:
    func = spec.makeFunc(algorithm(args.method))

  ROOT.gDirectory.Clear()
  
  # Create a RooFitHist object as input for the unfolding
  # results printing.
  truth_hist = spec.makeHistogram(histograms["sig_truth_test"])
  
  # Do the unfolding.
  unfold = func.unfolding()
  
  # Print the unfolding results and compare to a truth histogram.
  unfold.PrintTable(ROOT.cout, truth_hist)

  
  if not args.plot:
    return

  unf_hist = unfold.TH1reco()
  
  plot_hists(histograms["sig_truth_test"], histograms["sigbkg_reco_test"], histograms["sigbkg_reco_test"], unf_hist)

  
if __name__=="__main__":
  from argparse import ArgumentParser
  parser = ArgumentParser(description="RooUnfold testing script")
  parser.add_argument("method",default="bbb",type=str)
  parser.add_argument("--plot",action='store_true')
  parser.add_argument("--regparm",type=float)
  parser.add_argument("--overflow",action='store_true',default=False)
  main(parser.parse_args())
