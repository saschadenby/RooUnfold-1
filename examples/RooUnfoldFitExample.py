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

def prepare(smear):
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

  infile = ROOT.TFile.Open(infname,"RECREATE")
  clones = []
  for k in histograms.keys():
    clone = histograms[k].Clone()
    clone.SetName(k)
    clone.SetDirectory(infile)
    clones.append(clone)
  infile.Write()
  infile.Close()


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

def trySetName(obj,name):
  if obj:
    obj.SetName(name)
    obj.SetTitle(name)

def renameObservables(ws,observables):
  import ROOT
  for obs in observables.keys():
    name = observables[obs]
    trySetName(ws.var(obs),name)
    for ds in ws.allEmbeddedData():
      trySetName(ds.get(0).find(obs),name)
    for ds in ws.allData():
      trySetName(ds.get(0).find(obs),name)
    for hf in ROOT.RooUnfolding.matchingObjects(ws.allFunctions(),".*"):
      if hf.InheritsFrom(ROOT.RooHistFunc.Class()):
        trySetName(hf.dataHist().get(0).find(obs),name)
        trySetName(ROOT.RooUnfolding.getObservables(hf).find(obs),name)
      if hf.InheritsFrom(ROOT.ParamHistFunc.Class()):
        trySetName(hf.get(0).find(obs),name)

def makeFinalWorkspace(method,mainws,auxws):
  import ROOT
  observables = {"obs_x_response":"obs_reco","obs_y_response":"obs_truth","obs_x_SRtruth":"obs_truth","obs_x_truth":"obs_truth","obs_x_SRreco":"obs_reco"}
  
  helpws = ROOT.RooWorkspace("helper","helper")
  
  allvars = ROOT.RooArgSet()
  allvars.add(auxws.obj("obs_x_SRtruth"))
  sig_truth_dataset = auxws.pdf("model_SRtruth").generate(allvars,100,True)
  sig_truth_dataset.SetName("sig_truth_dataset")
  
  renameObservables(auxws,observables)
  renameObservables(mainws,observables)
  
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
  
  ROOT.RooUnfolding.importToWorkspace(helpws,sig_theory)
  ROOT.RooUnfolding.importToWorkspace(helpws,sig_response)
  ROOT.RooUnfolding.importToWorkspace(helpws,sig_truth)
  ROOT.RooUnfolding.importToWorkspace(helpws,sig_reco)
  ROOT.RooUnfolding.importToWorkspace(helpws,bkg_reco)
  ROOT.RooUnfolding.importToWorkspace(helpws,sigbkg_reco_dataset)
  ROOT.RooUnfolding.importToWorkspace(helpws,sig_truth_dataset)
  
  sig_response = helpws.function("sig_response")
  sig_theory = helpws.function("sig_theory")
  sigbkg_reco_dataset = helpws.data("sigbkg_reco_dataset").binnedClone()
  sigbkg_reco_dataset.SetName("sigbkg_reco_dataset")
  sig_truth_dataset = helpws.data("sig_truth_dataset")
  sig_truth = helpws.function("sig_truth")
  sig_reco = helpws.function("sig_reco")
  bkg_reco = helpws.function("bkg_reco")
  
  renameObservables(helpws,observables)
  
  unfold_response = ROOT.RooFitUnfoldResponse("unfold_response","unfold_response",sig_response,sig_truth,sig_reco,0,helpws.var("obs_truth"),helpws.var("obs_reco"))

  sigbkg_reco = unfold_response.makeHistFunc(sigbkg_reco_dataset)
  sigbkg_reco.SetName("sigbkg_reco_dataset")
  sig_reco_dataset = unfold_response.makeHistSum(sigbkg_reco,bkg_reco,1.,-1.)
  sig_reco_dataset.func().SetName("sig_reco_dataset")
  
  if method == "bayes":
    unfold= ROOT.RooFitUnfoldBayes     (unfold_response, sig_reco_dataset, 4);     #  OR
  elif method == "bbb":
    unfold= ROOT.RooFitUnfoldBinByBin     (unfold_response, sig_reco_dataset);     #  OR    
  elif method == "inv":
    unfold= ROOT.RooFitUnfoldInvert     (unfold_response, sig_reco_dataset);     #  OR    
  elif method == "svd":
    unfold= ROOT.RooFitUnfoldSvd     (unfold_response, sig_reco_dataset, 20);     #  OR
  elif method == "root":
    unfold= ROOT.RooFitUnfoldTUnfold     (unfold_response, sig_reco_dataset);     #  OR    
  elif method == "ids":
    unfold= ROOT.RooFitUnfoldIds     (unfold_response, sig_reco_dataset, 3);     #  OR      
  
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

def makePlots(ws):
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

def runFit(ws):
  mu = ws.var("mu")
  
  for var in ROOT.RooUnfolding.matchingObjects(ws.allVars(),"gamma_stat_.*"):
    var.setConstant(True)
  
  mu.setVal(1.02)
  unfolding = ws.pdf("unfolding")
  data = ws.data("asimovData")
  
  nps = ROOT.RooUnfolding.allVars(ws,"gamma_.*")
  unfolding.fitTo(data,ROOT.RooFit.GlobalObservables(nps))


def main(args):
  from os.path import exists
  
  def smear(xt):
    from ROOT import gRandom
    xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  #  efficiency
    x= gRandom.Rndm();
    if x>xeff: return None;
    xsmear= gRandom.Gaus(args.bias,args.smear);     #  bias and smear
    return xt+xsmear;

  prepare(smear)
  
  if exists("ws-aux.root") and exists("ws.root"):
    mainws = loadws("ws.root","combined")
    auxws = loadws("ws-aux.root","combined")
  else:
    mainws = makeWorkspace()
    auxws = makeAuxWorkspace()
    mainws.writeToFile("ws.root")
    auxws.writeToFile("ws-aux.root")
  
  import ROOT
  ROOT.RooUnfoldResponse
  
  ROOT.RooUnfolding.setGammaUncertainties(mainws)
  ROOT.RooUnfolding.setGammaUncertainties(auxws)
  
  import sys
  ws = makeFinalWorkspace(args.method,mainws,auxws)

  ws.writeToFile("unfolding.root")
  
  makePlots(ws)

  #runFit(ws)
  
if __name__=="__main__":
  from argparse import ArgumentParser
  parser = ArgumentParser(description="RooUnfold testing script")
  parser.add_argument("method",default="bayes",type=str)
  parser.add_argument("--bias",default=-2.5,type=float)
  parser.add_argument("--smear",default=.2,type=float)
  main(parser.parse_args())
