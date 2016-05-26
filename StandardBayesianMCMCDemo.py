# Standard demo of the Bayesian MCMC calculator

'''

Author: Kyle Cranmer
date: Dec. 2010
updated: July 2011 for 1-sided upper limit and SequentialProposalFunction

This is a standard demo that can be used with any ROOT file
prepared in the standard way.  You specify:
 - name for input ROOT file
 - name of workspace inside ROOT file that holds model and data
 - name of ROOT.RooStats.ModelConfig( that specifies details for calculator tools
 - name of dataset

With default parameters the macro will attempt to run the
standard hist2workspace example and read the ROOT file
that it produces.

The actual heart of the demo is only about 10 lines long.

The ROOT.RooStats.MCMCCalculator( is a Bayesian tool that uses
the Metropolis-Hastings algorithm to efficiently integrate
in many dimensions.  It is not as accurate as the BayesianCalculator
for simple problems, it scales to much more complicated cases.
'''

#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/SequentialProposal.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ProposalHelper.h"
#include "RooFitResult.h"


import ROOT
using namespace ROOT.RooStats

intervalType = 2;          # type of interval (0 is shortest, central, upper limit)
maxPOI = -999;        # force a different value of POI for doing the scan (default is given value)


void StandardBayesianMCMCDemo( infile = "",
                               workspaceName = "combined",
                               modelConfigName = "ModelConfig",
                               dataName = "obsData")
  ##############################/
  # First part is just to access a user-defined file
  # or create the standard example file if it doesn't exist
  ##############################

   

    filename = ""
   if not strcmp(infile, ""):      filename = "results/example_combined_GaussExample_model.root"
      fileExist = not gSystem.AccessPathName(filename); # note opposite return code
      # if file does not exists generate with histfactory
      if not fileExist:#ifdef _WIN32
         cout << "HistFactory file cannot be generated on Windows - exit" << endl
         return
#endif
         # Normally self would be run on the command line
         cout <<"will run standard hist2workspace example"<<endl
         gROOT.ProcessLine(".not  prepareHistFactory .")
         gROOT.ProcessLine(".not  hist2workspace config/example.xml")
         cout <<"\n\n---------------------"<<endl
         cout <<"Done creating example input"<<endl
         cout <<"---------------------\n\n"<<endl



   else:
      filename = infile

   # ROOT.Try to open the file
   ROOT.TFile *file = ROOT.TFile.Open(filename)

   # if input file was specified byt not found, quit
   if not file :      cout <<"StandardRooStatsDemoMacro: Input file " << filename << " is not found" << endl
      return




  ##############################/
  # ROOT.Tutorial starts here
  ##############################

  # get the workspace out of the file
  w = (ROOT.RooWorkspace*) file.Get(workspaceName)
  if not w:    cout <<"workspace not found" << endl
    return


  # get the modelConfig out of the file
  mc = (ModelConfig*) w.obj(modelConfigName)

  # get the modelConfig out of the file
  data = w.data(dataName)

  # make sure ingredients are found
  if not data or not mc:    w.Print()
    cout << "data or ROOT.RooStats.ModelConfig( was not found" <<endl
    return


  # Want an efficient proposal function
  # default is uniform.

  '''
  # self one is based on the covariance matrix of fit
  fit = mc.GetPdf().fitTo(*data, ROOT.RooFit.Save())
  ROOT.RooStats.ProposalHelper( ph
  ph.SetVariables((ROOT.RooArgSet&)fit.floatParsFinal())
  ph.SetCovMatrix(fit.covarianceMatrix())
  ph.SetUpdateProposalParameters(kTRUE); # auto-create mean vars and add mappings
  ph.SetCacheSize(100)
  pf = ph.GetProposalFunction()
  '''

  # self proposal function seems fairly robust
  SequentialProposal sp(0.1)
  ######################/
  # create and use the ROOT.RooStats.MCMCCalculator(
  # to find and plot the 95% credible interval
  # on the parameter of interest as specified
  # in the model config
  ROOT.RooStats.MCMCCalculator( mcmc(*data,*mc)
  mcmc.SetConfidenceLevel(0.95); # 95% interval
  #  mcmc.SetProposalFunction(*pf)
  mcmc.SetProposalFunction(sp)
  mcmc.SetNumIters(1000000);         # Metropolis-Hastings algorithm iterations
  mcmc.SetNumBurnInSteps(50);       # first N steps to be ignored as burn-in

  # default is the shortest interval.  
  if intervalType == 0)  mcmc.SetIntervalType(MCMCInterval.kShortest); # for shortest interval (not really needed:
  if (intervalType == 1)  mcmc.SetLeftSideTailFraction(0.5); # for central interval
  if (intervalType == 2)  mcmc.SetLeftSideTailFraction(0.); # for upper limit                                    

  firstPOI = (ROOT.RooRealVar*) mc.GetParametersOfInterest().first()
  if (maxPOI != -999) 
     firstPOI.setMax(maxPOI)

  interval = mcmc.GetInterval()

  # make a plot
  #TCanvas* c1 =
  ROOT.TCanvas("IntervalPlot")
  ROOT.RooStats.MCMCIntervalPlot( plot(*interval)
  plot.Draw()

  c2 = ROOT.TCanvas("extraPlots")
   list = mc.GetNuisanceParameters()
  if list.getSize()>1:    n = list.getSize()
    ny = ROOT.TMath.CeilNint( sqrt(n) )
    nx = ROOT.TMath.CeilNint(double(n)/ny)
    c2.Divide( nx,ny)


  # draw a scatter plot of chain results for poi vs each nuisance parameters
  it = mc.GetNuisanceParameters().createIterator()
  nuis = NULL
  int iPad=1; # iPad, that's funny
  while( (nuis = (ROOT.RooRealVar*) it.Next() ))    c2.cd(iPad++)
    plot.DrawChainScatter(*firstPOI,*nuis)


  # print out the iterval on the first Parameter of Interest
  cout << "\n95% interval on " <<firstPOI.GetName()<<" is : ["<<
    interval.LowerLimit(*firstPOI) << ", "<<
    interval.UpperLimit(*firstPOI) <<"] "<<endl


