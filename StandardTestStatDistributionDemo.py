'''
StandardTestStatDistributionDemo.C
author Kyle Cranmer
date: summer solstice, 2011

This simple script plots the sampling distribution of the profile likelihood
ratio test statistic based on the input Model File.  ROOT.To do self one needs to
specify the value of the parameter of interest that will be used for evaluating
the test statistic and the value of the parameters used for generating the toy data.
In self case, uses the upper-limit estimated from the ProfileLikleihoodCalculator,
which assumes the asymptotic chi-square distribution for -2 log profile likleihood ratio.
Thus, script is handy for checking to see if the asymptotic approximations are valid.
To aid, comparison, script overlays a chi-square distribution as well.
The most common parameter of interest is a parameter proportional to the signal rate,
and often that has a lower-limit of 0, breaks the standard chi-square distribution.
Thus the script allows the parameter to be negative so that the overlay chi-square is
the correct asymptotic distribution.
'''

#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TF1.h"
#include "TSystem.h"

#include "RooWorkspace.h"
#include "RooAbsData.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"

import ROOT
using namespace ROOT.RooStats


useProof = False;  # flag to control whether to use Proof
nworkers = 0;   # number of workers (default use all available cores)

####################################/
# ROOT.The actual macro

void StandardTestStatDistributionDemo( infile = "",
                                       workspaceName = "combined",
                                       modelConfigName = "ModelConfig",
                                       dataName = "obsData")

  # the number of toy MC used to generate the distribution
  nToyMC = 1000
  # ROOT.The parameter below is needed for asymptotic distribution to be chi-square,
  # but set to False if your model is not numerically stable if mu<0
  bool allowNegativeMu=True


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
  # Now get the data and workspace
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


  mc.Print()
  ##############################/
  # Now find the upper limit based on the asymptotic results
  ##############################
  firstPOI = (ROOT.RooRealVar*) mc.GetParametersOfInterest().first()
  ROOT.RooStats.ProfileLikelihoodCalculator( plc(*data,*mc)
  interval = plc.GetInterval()
  plcUpperLimit = interval.UpperLimit(*firstPOI)
  delete interval
  cout << "\n\n--------------------------------------"<<endl
  cout <<"Will generate sampling distribution at " << firstPOI.GetName() << " = " << plcUpperLimit <<endl
  nPOI = mc.GetParametersOfInterest().getSize()
  if nPOI>1:    cout <<"not sure what to do with other parameters of interest, here are their values"<<endl
    mc.GetParametersOfInterest().Print("v")


  ######################/
  # create thte test stat sampler
  ProfileLikelihoodTestStat ts(*mc.GetPdf())

  # to avoid effects from boundary and simplify asymptotic comparison, min=-max
  if allowNegativeMu:
    firstPOI.setMin(-1*firstPOI.getMax())

  # temporary ROOT.RooArgSet
  ROOT.RooArgSet poi
  poi.add(*mc.GetParametersOfInterest())

  # create and configure the ROOT.ToyMCSampler
  ROOT.ToyMCSampler sampler(ts,nToyMC)
  sampler.SetPdf(*mc.GetPdf())
  sampler.SetObservables(*mc.GetObservables())
  sampler.SetGlobalObservables(*mc.GetGlobalObservables())
  if not mc.GetPdf().canBeExtended() and (data.numEntries()==1):    cout << "tell it to use 1 event"<<endl
    sampler.SetNEventsPerToy(1)

  firstPOI.setVal(plcUpperLimit); # set POI value for generation
  sampler.SetParametersForTestStat(*mc.GetParametersOfInterest()); # set POI value for evaluation

  if (useProof) { 
     ProofConfig pc(*w, nworkers, "",False)
     sampler.SetProofConfig(&pc); # enable proof

     
  firstPOI.setVal(plcUpperLimit)
  ROOT.RooArgSet allParameters
  allParameters.add(*mc.GetParametersOfInterest())
  allParameters.add(*mc.GetNuisanceParameters())
  allParameters.Print("v")

  sampDist = sampler.GetSamplingDistribution(allParameters)
  SamplingDistPlot plot
  plot.AddSamplingDistribution(sampDist)
  plot.GetTH1F(sampDist).GetYaxis().SetTitle(Form("f(-log #lambda(#mu=%.2f) | #mu=%.2f)",plcUpperLimit,plcUpperLimit))
  plot.SetAxisTitle(Form("-log #lambda(#mu=%.2f)",plcUpperLimit))

  c1 = ROOT.TCanvas("c1")
  c1.SetLogy()
  plot.Draw()
  min = plot.GetTH1F(sampDist).GetXaxis().GetXmin()
  max = plot.GetTH1F(sampDist).GetXaxis().GetXmax()

  f = ROOT.TF1("f",Form("2*ROOT.Math.chisquared_pdf(2*x,%d,0)",nPOI),min,max)
  f.Draw("same")
  c1.SaveAs("standard_test_stat_distribution.pdf")


