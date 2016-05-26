# Example showing confidence intervals with four techniques.
'''
IntervalExamples

Author Kyle Cranmer
date   Sep. 2010

An example that shows confidence intervals with four techniques.
The model is a Normal Gaussian G(x|mu,sigma) with 100 samples of x.
The answer is known analytically, self is a good example to validate
the ROOT.RooStats tools.

expected interval is [-0.162917, 0.229075]
plc  interval is     [-0.162917, 0.229075]
fc   interval is     [-0.17    , 0.23]        # stepsize is 0.01
bc   interval is     [-0.162918, 0.229076]
mcmc interval is     [-0.166999, 0.230224]


'''

#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/ProofConfig.h"
#include "RooStats/ToyMCSampler.h"

#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddition.h"
#include "RooDataHist.h"
#include "RooPoisson.h"
#include "RooPlot.h"

#include "TCanvas.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMath.h"
#include"Math/DistFunc.h"
#include "TH1F.h"
#include "TMarker.h"
#include "TStopwatch.h"

#include <iostream>

# use self order for safety on library loading
import ROOT 
using namespace ROOT.RooStats 


def IntervalExamples(self):

  # ROOT.Time self macro
  ROOT.TStopwatch t
  t.Start()


  # set ROOT.RooFit random seed for reproducible results
  ROOT.RooRandom.randomGenerator().SetSeed(3001)

  # make a simple model via the workspace factory
  wspace = ROOT.RooWorkspace()
  wspace.factory("Gaussian.normal(x[-10,10],mu[-1,1],sigma[1])")
  wspace.defineSet("poi", "mu")
  wspace.defineSet("obs", "x")

  # specify components of model for statistical tools
  modelConfig = ROOT.RooStats.ModelConfig(("Example G(x|mu,1)")
  modelConfig.SetWorkspace(*wspace)
  modelConfig.SetPdf( *wspace.pdf("normal") )
  modelConfig.SetParametersOfInterest( *wspace.set("poi") )
  modelConfig.SetObservables( *wspace.set("obs") )

  # create a toy dataset
  data = wspace.pdf("normal").generate(*wspace.set("obs"),100)
  data.Print()

  # for convenience later on
  x = wspace.var("x")
  mu = wspace.var("mu")

  # set confidence level
  confidenceLevel = 0.95

  # example use profile likelihood calculator
  ROOT.RooStats.ProfileLikelihoodCalculator( plc(*data, *modelConfig)
  plc.SetConfidenceLevel( confidenceLevel)
  plInt = plc.GetInterval()

  # example use of Feldman-Cousins
  ROOT.RooStats.FeldmanCousins( fc(*data, *modelConfig)
  fc.SetConfidenceLevel( confidenceLevel)
  fc.SetNBins(100); # number of points to test per parameter
  fc.UseAdaptiveSampling(True); # make it go faster

  # Here, consider only ensembles with 100 events
  # ROOT.The PDF could be extended and self could be removed
  fc.FluctuateNumDataEntries(False)

  # Proof
  #  ProofConfig pc(*wspace, 4, "workers=4", kFALSE);    # proof-lite
  #ProofConfig pc(w, 8, "localhost");    # proof cluster at "localhost"
  #  toymcsampler = (ToyMCSampler*) fc.GetTestStatSampler()
  #  toymcsampler.SetProofConfig(&pc);     # enable proof

  interval = (PointSetInterval*) fc.GetInterval()


  # example use of BayesianCalculator
  # now we also need to specify a prior in the ROOT.RooStats.ModelConfig(
  wspace.factory("Uniform.prior(mu)")
  modelConfig.SetPriorPdf(*wspace.pdf("prior"))

  # example usage of BayesianCalculator
  BayesianCalculator bc(*data, *modelConfig)
  bc.SetConfidenceLevel( confidenceLevel)
  bcInt = bc.GetInterval()

  # example use of MCMCInterval
  ROOT.RooStats.MCMCCalculator( mc(*data, *modelConfig)
  mc.SetConfidenceLevel( confidenceLevel)
  # special options
  mc.SetNumBins(200);        # bins used internally for representing posterior
  mc.SetNumBurnInSteps(500); # first N steps to be ignored as burn-in
  mc.SetNumIters(100000);    # how long to run chain
  mc.SetLeftSideTailFraction(0.5); # for central interval
  mcInt = mc.GetInterval()

  # for self example we know the expected intervals
  expectedLL = data.mean(*x)
    + ROOT.Math.normal_quantile(  (1-confidenceLevel)/2,1)
    / sqrt(data.numEntries())
  expectedUL = data.mean(*x)
    + ROOT.Math.normal_quantile_c((1-confidenceLevel)/2,1)
    / sqrt(data.numEntries()) 

  # Use the intervals
  std.cout << "expected interval is [" <<
    expectedLL << ", " <<
    expectedUL << "]" << endl

  cout << "plc interval is [" <<
    plInt.LowerLimit(*mu) << ", " <<
    plInt.UpperLimit(*mu) << "]" << endl

  std.cout << "fc interval is ["<<
    interval.LowerLimit(*mu) << " , "  <<
    interval.UpperLimit(*mu) << "]" << endl

  cout << "bc interval is [" <<
    bcInt.LowerLimit() << ", " <<
    bcInt.UpperLimit() << "]" << endl

  cout << "mc interval is [" <<
    mcInt.LowerLimit(*mu) << ", " <<
    mcInt.UpperLimit(*mu) << "]" << endl

  mu.setVal(0)
  cout << "is mu=0 in the interval? " <<
    plInt.IsInInterval(ROOT.RooArgSet(*mu)) << endl


  # make a reasonable style
  gStyle.SetCanvasColor(0)
  gStyle.SetCanvasBorderMode(0)
  gStyle.SetPadBorderMode(0)
  gStyle.SetPadColor(0)
  gStyle.SetCanvasColor(0)
  gStyle.SetTitleFillColor(0)
  gStyle.SetFillColor(0)
  gStyle.SetFrameFillColor(0)
  gStyle.SetStatColor(0)


  # some plots
  canvas = ROOT.TCanvas("canvas")
  canvas.Divide(2,2)

  # plot the data
  canvas.cd(1)
  frame = x.frame()
  data.plotOn(frame)
  data.statOn(frame)
  frame.Draw()

  # plot the profile likeihood
  canvas.cd(2)
  ROOT.RooStats.LikelihoodIntervalPlot( plot(plInt)
  plot.Draw()

  # plot the MCMC interval
  canvas.cd(3)
  mcPlot = ROOT.RooStats.MCMCIntervalPlot((*mcInt)
  mcPlot.SetLineColor(ROOT.kGreen)
  mcPlot.SetLineWidth(2)
  mcPlot.Draw()

  canvas.cd(4)
  ROOT.RooPlot bcPlot = bc.GetPosteriorPlot()
  bcPlot.Draw()

  canvas.Update()

  t.Stop()
  t.Print()


