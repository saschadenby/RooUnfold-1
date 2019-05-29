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

import sys
method = "bayes"
if len(sys.argv) > 1: method = sys.argv[1]

from ROOT import gRandom, TH1, TH1D, TCanvas, cout
import ROOT

# ==============================================================================
#  Gaussian smearing, systematic translation, and variable inefficiency
# ==============================================================================

def smear(xt):
  xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  #  efficiency
  x= gRandom.Rndm();
  if x>xeff: return None;
  xsmear= gRandom.Gaus(-2.5,0.2);     #  bias and smear
  return xt+xsmear;

# ==============================================================================
#  Example Unfolding
# ==============================================================================

infname = "hist_tmp.root"


def prepare():
  response= ROOT.RooUnfoldResponse (40, -10.0, 10.0);
  
  #  Train with a Breit-Wigner, mean 0.3 and width 2.5.
  for i in range(100000):
    xt= gRandom.BreitWigner (0.3, 2.5);
    x= smear (xt);
    if x!=None:
      response.Fill (x, xt);
    else:
      response.Miss (xt);
      
  hTrue= TH1D ("true", "Test Truth",    40, -10.0, 10.0);
  hMeas= TH1D ("meas", "Test Measured", 40, -10.0, 10.0);
  #  Test with a Gaussian, mean 0 and width 2.
  for i in range(10000):
    xt= gRandom.Gaus (0.0, 2.0)
    x= smear (xt);
    hTrue.Fill(xt);
    if x!=None: hMeas.Fill(x);

  hBkg= TH1D ("bkg", "Test Bkg",    40, -10.0, 10.0);
  for i in range(10000):
    x= gRandom.Uniform (-10,10.)
    hBkg.Fill(x)
    xd= gRandom.Uniform (-10,10.)    
    hMeas.Fill(xd)

  histograms = {"data_truth":hTrue,"measurement":hMeas,"sig_response":response.Hresponse(),"sig":response.Hmeasured(),"sig_truth":response.Htruth(),"bkg":hBkg}

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
  signal_truth = ROOT.RooStats.HistFactory.Sample( "signal", "sig_truth", infname )
  signal_truth.GetStatError().Activate(True);  
  signal_truth.AddNormFactor( "mu", 1, 0, 3 )
  srt.AddSample( signal_truth )
  meas.AddChannel( srt )
  ##############################################################################    

  ##############################################################################
  # Create the response mapping
  response = ROOT.RooStats.HistFactory.Channel( "response" )
  sample = ROOT.RooStats.HistFactory.Sample( "signal", "sig_response", infname )
  sample.GetStatError().Activate(False);  
  response.AddSample( sample )  
  meas.AddChannel( response )
  
  truthRegion = ROOT.RooStats.HistFactory.Channel( "truth" )
  dataTruth = ROOT.RooStats.HistFactory.Sample( "data_truth", "data_truth", infname )
  dataTruth.GetStatError().Activate(False);  
  truthRegion.AddSample( dataTruth )  
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
  signal = ROOT.RooStats.HistFactory.Sample( "signal", "sig", infname )
  signal.GetStatError().Activate(True);  
  signal.AddNormFactor( "mu", 1, 0, 3 )
  sr.AddSample( signal )
  # Background 
  background = ROOT.RooStats.HistFactory.Sample( "background", "bkg", infname )
  sr.AddSample( background )
  # Data 
  sr.SetData( "measurement", infname )
  meas.AddChannel( sr )
  ##############################################################################

  # Collect the histograms from their files,
  # print some output, 
  meas.CollectHistograms()
 
  # Now, do the measurement
  ws = ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast( meas );
  return ws

def loadws(fname,wsname):
  f = ROOT.TFile.Open(fname,"READ")
  ws = f.Get(wsname).Clone()
  f.Close()
  return ws

def trySetName(obj,name):
  if obj:
    obj.SetName(name)
    obj.SetTitle(name)

def renameObservables(ws,observables):
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

from os.path import exists

prepare()

if exists("ws-aux.root") and exists("ws.root"):
  mainws = loadws("ws.root","combined")
  auxws = loadws("ws-aux.root","combined")
else:
  mainws = makeWorkspace()
  auxws = makeAuxWorkspace()
  mainws.writeToFile("ws.root")
  auxws.writeToFile("ws-aux.root")


observables = {"obs_x_response":"obs_reco","obs_y_response":"obs_truth","obs_x_SRtruth":"obs_truth","obs_x_truth":"obs_truth","obs_x_SRreco":"obs_reco"}

ROOT.RooUnfoldResponse

ROOT.RooUnfolding.setGammaUncertainties(mainws)
ROOT.RooUnfolding.setGammaUncertainties(auxws)

helpws = ROOT.RooWorkspace("helper","helper")

allvars = ROOT.RooArgSet()
allvars.add(auxws.obj("obs_x_SRtruth"))
asimov = auxws.pdf("model_SRtruth").generate(allvars,100,True)
asimov.SetName("asimovData")

renameObservables(auxws,observables)
renameObservables(mainws,observables)

response = auxws.function("signal_response_nominal")
data_truth = auxws.function("data_truth_truth_nominal")
data = mainws.data("obsData")
truth = auxws.function("L_x_signal_SRtruth_overallSyst_x_StatUncert")
auxws.obj("obs_truth").setVal(0)
auxws.obj("obs_truth").setVal(10)

reco = mainws.function("L_x_signal_SRreco_overallSyst_x_StatUncert")
bkg = mainws.function("L_x_background_SRreco_overallSyst_x_Exp")
mc = auxws.obj("ModelConfig")

ROOT.RooUnfolding.importToWorkspace(helpws,data_truth)
ROOT.RooUnfolding.importToWorkspace(helpws,response)
ROOT.RooUnfolding.importToWorkspace(helpws,truth)
ROOT.RooUnfolding.importToWorkspace(helpws,reco)
ROOT.RooUnfolding.importToWorkspace(helpws,bkg)
ROOT.RooUnfolding.importToWorkspace(helpws,data)
ROOT.RooUnfolding.importToWorkspace(helpws,asimov)

response = helpws.function("signal_response_nominal")
data_truth = helpws.function("data_truth_truth_nominal")
truth = helpws.function("L_x_signal_SRtruth_overallSyst_x_StatUncert")
reco = helpws.function("L_x_signal_SRreco_overallSyst_x_StatUncert")
bkg = helpws.function("L_x_background_SRreco_overallSyst_x_Exp")
data = helpws.data("obsData")

renameObservables(helpws,observables)

unfolding = ROOT.RooFitUnfoldResponse("unfolding","unfolding",response,truth,reco,0,helpws.var("obs_truth"),helpws.var("obs_reco"));

hist = unfolding.makeHistSum(unfolding.makeHistFunc(data.binnedClone()),bkg,1.,-1.)

if method == "bayes":
  unfold= ROOT.RooFitUnfoldBayes     (unfolding, hist, 4);     #  OR
elif method == "bbb":
  unfold= ROOT.RooFitUnfoldBinByBin     (unfolding, hist);     #  OR    
elif method == "svd":
  unfold= ROOT.RooFitUnfoldSvd     (unfolding, hist, 20);     #  OR
elif method == "root":
  unfold= ROOT.RooFitUnfoldTUnfold     (unfolding, hist);     #  OR    
elif method == "ids":
  unfold= ROOT.RooFitUnfoldIds     (unfolding, hist, 3);     #  OR      

truthhist = unfolding.makeHist(data_truth)
unfold.PrintTable(cout,truthhist)

unfoldpdf = ROOT.RooUnfoldPdf("unfolding","unfolding",unfold)
components = ROOT.RooArgList()
components.add(unfoldpdf)
for constr in ROOT.RooUnfolding.matchingObjects(mainws.allPdfs(),"gamma_stat_.*_constraint"):
  components.add(constr)
for constr in ROOT.RooUnfolding.matchingObjects(auxws.allPdfs(),"gamma_stat_.*_constraint"):
  components.add(constr)
pdf = ROOT.RooProdPdf("mainPdf","mainPdf",components)

ws = ROOT.RooWorkspace("workspace","workspace")
ROOT.RooUnfolding.importToWorkspace(ws,pdf)
ROOT.RooUnfolding.importToWorkspace(ws,helpws.data("asimovData"))
#ws.var("obs_reco").setConstant(True)
ws.var("obs_truth").setConstant(True)
#ws.Print("t")
ws.writeToFile("unfolding.root")

wspdf = ws.pdf("unfolding")

obs_truth = ws.var("obs_truth")
plot_truth = obs_truth.frame()
data_truth.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kRed))
wspdf.plotOn(plot_truth)
canvas_truth = ROOT.TCanvas("unfolded","unfolded")
plot_truth.Draw()
canvas_truth.SaveAs("unfolded.pdf")

hmeas = unfoldpdf.unfolding().response().Hmeasured()
obs_reco = hmeas.obs(0)
plot_reco = obs_reco.frame()
data.plotOn(plot_reco)
hmeas.func().plotOn(plot_reco,ROOT.RooFit.LineColor(ROOT.kRed))
canvas_reco = ROOT.TCanvas("reco","reco")
plot_reco.Draw()
canvas_reco.SaveAs("reco.pdf")

#
#
#exit(0)
#
#mu = ws.var("mu")
#
#for var in ROOT.RooUnfolding.matchingObjects(ws.allVars(),"gamma_stat_.*"):
#  var.setConstant(True)
#
#mu.setVal(1.02)
#unfolding = ws.pdf("unfolding")
#data = ws.data("asimovData")
#
#nps = ROOT.RooUnfolding.allVars(ws,"gamma_.*")
#
#
#unfolding.fitTo(data,ROOT.RooFit.GlobalObservables(nps))
