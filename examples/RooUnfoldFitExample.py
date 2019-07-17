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

infname = "hist_tmp.root"

def prepare(smear,write):

  import ROOT

  response= ROOT.RooUnfoldResponse (40, -10.0, 10.0);

  nevents = 100000

  sig_theory= ROOT.TH1D ("sig_theory", "Test Truth",    40, -10.0, 10.0);
  sig_theory_alt= ROOT.TH1D ("sig_theory_alt", "Test Truth",    40, -10.0, 10.0);
  sigbkg_reco= ROOT.TH1D ("sigbkg_reco", "Test Measured", 40, -10.0, 10.0);
  sigbkg_reco_alt= ROOT.TH1D ("sigbkg_reco_alt", "Test Measured", 40, -10.0, 10.0);
  bkg_reco= ROOT.TH1D ("bkg", "Test Bkg",    40, -10.0, 10.0);


  #  Train with a Breit-Wigner, mean 0.3 and width 2.5.
  for i in range(nevents):
    xt= ROOT.gRandom.BreitWigner (0.3, 2.5);
    x= smear (xt);
    if x!=None:
      response.Fill (x, xt);
    else:
      response.Miss (xt);

    xt2= ROOT.gRandom.BreitWigner (0.3, 2.5);
    x2 = smear(xt2)
    sig_theory.Fill(xt2)

    if x2!=None: sigbkg_reco.Fill(x2)
      
  #  Test with a Gaussian, mean 0 and width 2.
  for i in range(nevents):
    xt= ROOT.gRandom.Gaus (0.0, 2.0)
    x= smear (xt);
    sig_theory_alt.Fill(xt);
    if x!=None: sigbkg_reco_alt.Fill(x);

  for i in range(nevents):
    x= ROOT.gRandom.Uniform (-10,10.)
    bkg_reco.Fill(x)
    xd= ROOT.gRandom.Uniform (-10,10.)    
    sigbkg_reco.Fill(xd)
    sigbkg_reco_alt.Fill(xd)

  histograms = {"sig_theory":sig_theory,
                "sig_theory_alt":sig_theory_alt,
                "sigbkg_reco":sigbkg_reco,
                "sigbkg_reco_alt":sigbkg_reco_alt,
                "sig_response":response.Hresponse(),
                "sig_reco":response.Hmeasured(),
                "sig_truth":response.Htruth(),
                "bkg_reco":bkg_reco}

  if write:
    infile = ROOT.TFile.Open(infname,"RECREATE")
    clones = []
    for k in histograms.keys():
      clone = histograms[k].Clone()
      clone.SetName(k)
      clone.SetDirectory(infile)
      clones.append(clone)
    infile.Write()
    infile.Close()

  for h in histograms.values():
    h.SetDirectory(0)

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
  return alg

def makeAuxWorkspace():
  # Create the measurement
  import ROOT
  meas = ROOT.RooStats.HistFactory.Measurement("unfolding", "unfolding")
 
  meas.SetPOI( "mu" )
  meas.AddConstantParam("Lumi")
 
  meas.SetLumi( 1.0 )
  meas.SetLumiRelErr( 0.10 )
  meas.SetExportOnly( False )

  ##############################################################################
  # Create the reco signal region
  srt = ROOT.RooStats.HistFactory.Channel( "SRtruth" )
  srt.SetStatErrorConfig( 0.005, "Poisson" )
  # Create the signal sample
  signal_truth = ROOT.RooStats.HistFactory.Sample( "sig_truth", "sig_truth", infname )
  signal_truth.GetStatError().Activate(True);  
  signal_truth.AddNormFactor( "mu", 1, 0, 3 )
  srt.AddSample( signal_truth )
  meas.AddChannel( srt )
  ##############################################################################    

  ##############################################################################
  # Create the response mapping
  response = ROOT.RooStats.HistFactory.Channel( "response" )
  sample = ROOT.RooStats.HistFactory.Sample( "sig_response", "sig_response", infname )
  sample.GetStatError().Activate(False);  
  response.AddSample( sample )  
  meas.AddChannel( response )
  
  truthRegion = ROOT.RooStats.HistFactory.Channel( "truth" )
  sig_theory = ROOT.RooStats.HistFactory.Sample( "sig_theory", "sig_theory", infname )
  sig_theory.GetStatError().Activate(False);  
  truthRegion.AddSample( sig_theory )  
  sig_theory_alt = ROOT.RooStats.HistFactory.Sample( "sig_theory_alt", "sig_theory_alt", infname )
  sig_theory_alt.GetStatError().Activate(False);  
  truthRegion.AddSample( sig_theory_alt )  
  meas.AddChannel( truthRegion )
  ##############################################################################  
 
  # Collect the histograms from their files,
  # print some output, 
  meas.CollectHistograms()
 
  # Now, do the measurement
  ws = ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast( meas );
  return ws

def makeWorkspace():
  # Create the measurement
  import ROOT
  meas = ROOT.RooStats.HistFactory.Measurement("unfolding", "unfolding")
 
  meas.SetPOI( "mu" )
  meas.AddConstantParam("Lumi")
 
  meas.SetLumi( 1.0 )
  meas.SetLumiRelErr( 0.10 )
  meas.SetExportOnly( False )

  ##############################################################################
  # Create the reco signal region
  sr = ROOT.RooStats.HistFactory.Channel( "SRreco" )
  sr.SetStatErrorConfig( 0.005, "Poisson" )
  # Create the signal sample
  signal = ROOT.RooStats.HistFactory.Sample( "sig_reco", "sig_reco", infname )
  signal.GetStatError().Activate(True);  
  signal.AddNormFactor( "mu", 1, 0, 3 )
  sr.AddSample( signal )
  # Background 
  background = ROOT.RooStats.HistFactory.Sample( "bkg_reco", "bkg_reco", infname )
  sr.AddSample( background )
  # Data 
  sr.SetData( "sigbkg_reco_alt", infname )
  meas.AddChannel( sr )
  ##############################################################################

  # Collect the histograms from their files,
  # print some output, 
  meas.CollectHistograms()
 
  # Now, do the measurement
  ws = ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast( meas );
  return ws

def loadws(fname,wsname):
  import ROOT
  f = ROOT.TFile.Open(fname,"READ")
  ws = f.Get(wsname).Clone()
  f.Close()
  return ws

def wsimport(ws,obj):
  import ROOT
  getattr(ws,"import")(obj,
                       ROOT.RooFit.RenameVariable("obs_x_response","obs_reco"),
                       ROOT.RooFit.RenameVariable("obs_y_response","obs_truth"),
                       ROOT.RooFit.RenameVariable("obs_x_SRtruth","obs_truth"),
                       ROOT.RooFit.RenameVariable("obs_x_truth","obs_truth"),
                       ROOT.RooFit.RenameVariable("obs_x_SRreco","obs_reco"))
  

def makeFinalWorkspace(method,mainws,auxws,regparm):
  import ROOT
  helpws = ROOT.RooWorkspace("helper","helper")
  
  allvars = ROOT.RooArgSet()
  allvars.add(auxws.obj("obs_x_SRtruth"))
  sig_truth_dataset = auxws.pdf("model_SRtruth").generate(allvars,100,True)
  sig_truth_dataset.SetName("sig_truth_dataset")

  sig_response = auxws.function("sig_response_response_nominal")
  sig_response.SetName("sig_response")
  sig_theory = auxws.function("sig_theory_alt_truth_nominal")
  sig_theory.SetName("sig_theory")
  sigbkg_reco_dataset = mainws.data("obsData")
  sigbkg_reco_dataset.SetName("sigbkg_reco_dataset")
  sig_truth = auxws.function("L_x_sig_truth_SRtruth_overallSyst_x_StatUncert")
  sig_truth.SetName("sig_truth")
  sig_reco = mainws.function("L_x_sig_reco_SRreco_overallSyst_x_StatUncert")
  sig_reco.SetName("sig_reco")
  bkg_reco = mainws.function("L_x_bkg_reco_SRreco_overallSyst_x_Exp")
  bkg_reco.SetName("bkg_reco")
  mc = auxws.obj("ModelConfig")

  wsimport(helpws,sig_theory)
  wsimport(helpws,sig_response)
  wsimport(helpws,sig_truth)
  wsimport(helpws,sig_reco)
  wsimport(helpws,bkg_reco)
  wsimport(helpws,sigbkg_reco_dataset)
  wsimport(helpws,sig_truth_dataset)

  sig_response = helpws.function("sig_response")
  sig_theory = helpws.function("sig_theory")
  sigbkg_reco_dataset = helpws.data("sigbkg_reco_dataset").binnedClone()
  sigbkg_reco_dataset.SetName("sigbkg_reco_dataset")
  sig_truth_dataset = helpws.data("sig_truth_dataset")
  sig_truth = helpws.function("sig_truth")
  sig_reco = helpws.function("sig_reco")
  bkg_reco = helpws.function("bkg_reco")
  
  unfold_response = ROOT.RooFitUnfoldResponse("unfold_response","unfold_response",sig_response,sig_truth,sig_reco,0,helpws.var("obs_truth"),helpws.var("obs_reco"))

  sigbkg_reco = unfold_response.makeHistFunc(sigbkg_reco_dataset)
  sigbkg_reco.SetName("sigbkg_reco_dataset")
  sig_reco_dataset = unfold_response.makeHistSum(sigbkg_reco,bkg_reco,1.,-1.)
  sig_reco_dataset.func().SetName("sig_reco_dataset")
  
  if method == ROOT.RooUnfolding.kBayes:
    unfold= ROOT.RooFitUnfoldBayes     (unfold_response, sig_reco_dataset, regparm);     #  OR
  elif method == ROOT.RooUnfolding.kBinByBin:
    unfold= ROOT.RooFitUnfoldBinByBin     (unfold_response, sig_reco_dataset);     #  OR    
  elif method  == ROOT.RooUnfolding.kInvert:
    unfold= ROOT.RooFitUnfoldInvert     (unfold_response, sig_reco_dataset);     #  OR    
  elif method == ROOT.RooUnfolding.kSVD:
    unfold= ROOT.RooFitUnfoldSvd     (unfold_response, sig_reco_dataset, regparm);     #  OR
  elif method == ROOT.RooUnfolding.kTUnfold:
    unfold= ROOT.RooFitUnfoldTUnfold     (unfold_response, sig_reco_dataset);     #  OR    
  elif method == ROOT.RooUnfolding.kIDs:
    unfold= ROOT.RooFitUnfoldIds     (unfold_response, sig_reco_dataset, regparm);     #  OR      
  
  sig_theory_hist = unfold_response.makeHist(sig_theory)

  unfold.PrintTable(ROOT.cout,sig_theory_hist)

  unfoldpdf = ROOT.RooUnfoldPdf("unfolding","unfolding",unfold)
  components = ROOT.RooArgList()
  components.add(unfoldpdf)
  for constr in ROOT.RooUnfolding.matchingObjects(mainws.allPdfs(),"gamma_stat_.*_constraint"):
    components.add(constr)
  for constr in ROOT.RooUnfolding.matchingObjects(auxws.allPdfs(),"gamma_stat_.*_constraint"):
    components.add(constr)
  pdf = ROOT.RooProdPdf("mainPdf","mainPdf",components)
  
  ws = ROOT.RooWorkspace("workspace","workspace")
  ROOT.RooUnfolding.importToWorkspace(ws,sig_theory)
  ROOT.RooUnfolding.importToWorkspace(ws,sig_reco)
  ROOT.RooUnfolding.importToWorkspace(ws,pdf)
  ROOT.RooUnfolding.importToWorkspace(ws,sig_truth_dataset)
  ROOT.RooUnfolding.importToWorkspace(ws,sigbkg_reco_dataset)
  ws.var("obs_truth").setConstant(True)

  return ws

def makePlots_HistFactory(ws):
  import ROOT

  unfoldpdf = ws.pdf("unfolding")
  unfolding = unfoldpdf.unfolding()

  sig_response = ws.function("sig_response")
  sig_theory = ws.function("sig_theory")
  sigbkg_reco_dataset = ws.data("sigbkg_reco_dataset")
  sig_truth_dataset = ws.data("sig_truth_dataset")
  sig_truth = ws.function("sig_truth")
  sig_reco = ws.function("sig_reco")
  bkg_reco = ws.function("bkg_reco")

  obs_truth = ws.var("obs_truth")
  plot_truth = obs_truth.frame()
  sig_theory.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sig_theory"))
  unfoldpdf.plotOn(plot_truth,ROOT.RooFit.Name("unfoldedpdf"))
  canvas_truth = ROOT.TCanvas("unfolded","unfolded")
  plot_truth.Draw()
  leg_truth = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
  leg_truth.AddEntry( plot_truth.findObject("sig_theory"), "Theory Prediction", "l" )
  leg_truth.AddEntry( plot_truth.findObject("unfoldedpdf"), "Unfolded Signal", "l" )
  leg_truth.Draw()
  canvas_truth.SaveAs("unfolded.pdf")
  canvas_truth.SaveAs("unfolded.png")
  
  sig_reco_hist = unfolding.response().Hmeasured()
  sigbkg_reco_hist = unfolding.response().makeHistSum(sig_reco_hist.func(),bkg_reco,1.,1.)
  obs_reco = sig_reco_hist.obs(0)
  plot_reco = obs_reco.frame()
  sigbkg_reco_dataset.plotOn(plot_reco,ROOT.RooFit.Name("data"))
  sigbkg_reco_hist.func().plotOn(plot_reco,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sigbkg_reco_hist"))
  canvas_reco = ROOT.TCanvas("reco","reco")
  plot_reco.Draw()
  leg_reco = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
  leg_reco.AddEntry( plot_reco.findObject("data"), "Data", "l" )
  leg_reco.AddEntry( plot_reco.findObject("sigbkg_reco_hist"), "Reco Monte Carlo", "l" )
  leg_reco.Draw()
  canvas_reco.SaveAs("reco.pdf")
  canvas_reco.SaveAs("reco.png")

def makePlots_RooUnfoldSpec(ws):
  import ROOT

  unfoldpdf = ws.pdf("unfold")
  unfolding = unfoldpdf.unfolding()

  sig_response = ws.function("response")
  sig_theory = ws.function("sig_theory")
  sigbkg_reco = ws.function("sigbkg_reco")
  sig_truth = ws.function("sig_truth")
  sig_reco = ws.function("unfold_data_minus_bkg")
  bkg_reco = ws.function("bkg")

  obs_truth = ws.var("obs_truth")
  plot_truth = obs_truth.frame()
  sig_theory.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sig_theory_graph"))
  unfoldpdf.plotOn(plot_truth,ROOT.RooFit.Name("unfoldedpdf_graph"))
  canvas_truth = ROOT.TCanvas("unfolded","unfolded")
  plot_truth.Draw()
  leg_truth = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
  leg_truth.AddEntry( plot_truth.findObject("sig_theory_graph"), "Theory Prediction", "l" )
  leg_truth.AddEntry( plot_truth.findObject("unfoldedpdf_graph"), "Unfolded Signal", "l" )
  leg_truth.Draw()
  canvas_truth.SaveAs("unfolded.pdf")
  canvas_truth.SaveAs("unfolded.png")
  
  obs_reco = ws.var("obs_reco")
  plot_reco = obs_reco.frame()
  sigbkg_reco.plotOn(plot_reco,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("sigbkg_reco_graph"))
  bkg_reco.plotOn(plot_reco,ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.Name("bkg_reco_graph"))
  canvas_reco = ROOT.TCanvas("reco","reco")
  plot_reco.Draw()
  leg_reco = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
  leg_reco.AddEntry( plot_reco.findObject("bkg_reco_graph"), "Background", "l" )
  leg_reco.AddEntry( plot_reco.findObject("sigbkg_reco_graph"), "Signal + Background", "l" )
  leg_reco.Draw()
  canvas_reco.SaveAs("reco.pdf")
  canvas_reco.SaveAs("reco.png")

def runFit(ws):
  mu = ws.var("mu")
  
  for var in ROOT.RooUnfolding.matchingObjects(ws.allVars(),"gamma_stat_.*"):
    var.setConstant(True)
  
  mu.setVal(1.02)
  unfolding = ws.pdf("unfolding")
  data = ws.data("asimovData")
  
  nps = ROOT.RooUnfolding.allVars(ws,"gamma_.*")
  unfolding.fitTo(data,ROOT.RooFit.GlobalObservables(nps))

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
  
  def smear(xt):
    from ROOT import gRandom
    xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  #  efficiency
    x= gRandom.Rndm();
    if x>xeff: return None;
    xsmear= gRandom.Gaus(args.bias,args.smear);     #  bias and smear
    return xt+xsmear;

  import ROOT
  if args.mode == "HistFactory":
    prepare(smear,True)
  
    if exists("ws-aux.root") and exists("ws.root"):
      mainws = loadws("ws.root","combined")
      auxws = loadws("ws-aux.root","combined")
    else:
      mainws = makeWorkspace()
      auxws = makeAuxWorkspace()
      mainws.writeToFile("ws.root")
      auxws.writeToFile("ws-aux.root")
  
    ROOT.RooUnfoldResponse
    
    ROOT.RooUnfolding.setGammaUncertainties(mainws)
    ROOT.RooUnfolding.setGammaUncertainties(auxws)
  
    import sys
    ws = makeFinalWorkspace(algorithm(args.method),mainws,auxws,4)
  else:
    histograms = prepare(smear,True)

    spec = ROOT.RooUnfoldSpec("unfold","unfold",histograms["sig_theory"],"obs_truth",histograms["sig_reco"],"obs_reco",histograms["sig_response"],histograms["bkg_reco"],histograms["sigbkg_reco"],False,0.0005)
    
    # This line instantiates one of the subclasses specific to the unfolding algorithm.
    pdf = spec.makePdf(algorithm(args.method),args.regparm)

    theory = pdf.unfolding().response().makeHistFuncTruth(histograms["sig_theory"])

    # required to avoid python garbage collector messing up the RooDataHists added to gDirectory
    ROOT.gDirectory.Clear()

    test_truth = spec.makeHistogram(histograms["sig_theory"])

    pdf.unfolding().PrintTable(ROOT.cout, test_truth)

    ws = ROOT.RooWorkspace("workspace","workspace")

    getattr(ws,"import")(pdf)

    getattr(ws,"import")(theory)

    ws.Print("t")

    # Explicitly free the memory that was allocated. PyROOT does not do this
    # for you.
    pdf.Delete()
    theory.Delete()

  ws.writeToFile("unfolding.root")

  if args.mode == "HistFactory":
    makePlots_HistFactory(ws)
  else:
    makePlots_RooUnfoldSpec(ws)


  
if __name__=="__main__":
  from argparse import ArgumentParser
  parser = ArgumentParser(description="RooUnfold testing script")
  parser.add_argument("method",default="bbb",type=str)
  parser.add_argument("--regparm",default=1,type=float)
  parser.add_argument("--mode",choices=["RooUnfoldSpec","HistFactory"],default="RooUnfoldSpec",type=str)
  parser.add_argument("--bias",default=-2.5,type=float)
  parser.add_argument("--smear",default=.2,type=float)
  main(parser.parse_args())
