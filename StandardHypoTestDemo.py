# Standard tutorial macro for hypothesis test (for computing the discovery significance) using all
# ROOT.RooStats hypotheiss tests calculators and test statistics
#
#
#Author:  L. Moneta
#
# Usage:
#
# root>.L StandardHypoTestDemo.C
# root> StandardHypoTestDemo("fileName", "workspace name", "S+B modelconfig name", "B model name", "data set name", type, statistic type, #                             number of toys)
#
#
# type = 0 Freq calculator
# type = 1 Hybrid calculator
# type = 2 Asymptotic calculator
# type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)
#
# testStatType = 0 LEP
#              = 1 ROOT.Tevatron
#              = 2 Profile Likelihood
#              = 3 Profile Likelihood one sided (i.e. = 0 if mu_hat < 0)
#


#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSystem.h"
#include "TROOT.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

import ROOT
using namespace ROOT.RooStats


noSystematics = False;              # force all systematics to be off (i.e. set all nuisance parameters as constat
nToysRatio = 4;                   # ratio Ntoys Null/ntoys ALT
poiValue = -1;                    # change poi snapshot value for S+B model (needed for expected p0 values)
int  printLevel=0
generateBinned = False;             # for binned generation

void StandardHypoTestDemo( infile = "",
                           workspaceName = "combined",
                           modelSBName = "ModelConfig",
                           modelBName = "",
                           dataName = "obsData",
                          calcType = 0, # 0 freq 1 hybrid, asymptotic
                          testStatType = 3,   # 0 LEP, ROOT.TeV, LHC, LHC - one sided
                          ntoys = 5000,
                          useNC = False,
                           char nuisPriorName = 0)

'''

  Other Parameter to pass in tutorial
  apart from standard for filename, ws, and data

  type = 0 Freq calculator
  type = 1 Hybrid calculator
  type = 2 Asymptotic calculator

  testStatType = 0 LEP
  = 1 ROOT.Tevatron
  = 2 Profile Likelihood
  = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)

  ntoys:         number of toys to use

  useNumberCounting:  set to ROOT.True when using number counting events

  nuisPriorName:   name of prior for the nnuisance. ROOT.This is often expressed as constraint term in the global model
  It is needed only when using the HybridCalculator (type=1)
  If not given by default the prior pdf from ROOT.RooStats.ModelConfig( is used.

  extra options are available as global paramwters of the macro. ROOT.They major ones are:

  generateBinned       generate binned data sets for toys (default is False) - be careful not to activate with
  a too large (>=3) number of observables
  nToyRatio            ratio of S+B/B toys (default is 2)
  printLevel

'''

   # disable - can cause some problems
   #ToyMCSampler.SetAlwaysUseMultiGen(True)

   SimpleLikelihoodRatioTestStat.SetAlwaysReuseNLL(True)
   ProfileLikelihoodTestStat.SetAlwaysReuseNLL(True)
   RatioOfProfiledLikelihoodsTestStat.SetAlwaysReuseNLL(True)

   #RooRandom.randomGenerator().SetSeed(0)

   # to change minimizers
   # ROOT.Math.MinimizerOptions.SetDefaultStrategy(0)
   # ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
   # ROOT.Math.MinimizerOptions.SetDefaultTolerance(1)

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

  w.Print()

  # get the modelConfig out of the file
  sbModel = (ModelConfig*) w.obj(modelSBName)


  # get the modelConfig out of the file
  data = w.data(dataName)

  # make sure ingredients are found
  if not data or not sbModel:    w.Print()
    cout << "data or ROOT.RooStats.ModelConfig( was not found" <<endl
    return

  # make b model
  bModel = (ModelConfig*) w.obj(modelBName)


   # case of no systematics
   # remove nuisance parameters from model
   if noSystematics:       ROOT.RooArgSet nuisPar = sbModel.GetNuisanceParameters()
      if nuisPar and nuisPar.getSize() > 0:         std.cout << "StandardHypoTestInvDemo" << "  -  Switch off all systematics by setting them constant to their initial values" << std.endl
         ROOT.RooStats.SetAllConstant(*nuisPar)

      if bModel:          ROOT.RooArgSet bnuisPar = bModel.GetNuisanceParameters()
         if bnuisPar:
            ROOT.RooStats.SetAllConstant(*bnuisPar)




  if not bModel :      Info("StandardHypoTestInvDemo", "The background model %s does not exist",modelBName)
      Info("StandardHypoTestInvDemo", "Copy it from ROOT.RooStats.ModelConfig( %s and set POI to zero",modelSBName)
      bModel = (ModelConfig*) sbModel.Clone()
      bModel.SetName(TString(modelSBName)+TString("B_only"))
      ROOT.RooRealVar var = dynamic_cast<RooRealVar*>(bModel.GetParametersOfInterest().first())
      if (not var) return
      oldval = var.getVal()
      var.setVal(0)
      #bModel.SetSnapshot( ROOT.RooArgSet(*var, *w.var("lumi"))  )
      bModel.SetSnapshot( ROOT.RooArgSet(*var)  )
      var.setVal(oldval)


   if not sbModel.GetSnapshot() or poiValue > 0:      Info("StandardHypoTestDemo", "Model %s has no snapshot  - make one using model poi",modelSBName)
      ROOT.RooRealVar var = dynamic_cast<RooRealVar*>(sbModel.GetParametersOfInterest().first())
      if (not var) return
      oldval = var.getVal()
      if poiValue > 0)  var.setVal(poiValue:
      #sbModel.SetSnapshot( ROOT.RooArgSet(*var, *w.var("lumi") ) )
      sbModel.SetSnapshot( ROOT.RooArgSet(*var) )
      if poiValue > 0) var.setVal(oldval:
      #sbModel.SetSnapshot( *sbModel.GetParametersOfInterest() )






   # part 1, testing
   SimpleLikelihoodRatioTestStat slrts = SimpleLikelihoodRatioTestStat(*bModel.GetPdf(), *sbModel.GetPdf())
   # null parameters must includes snapshot of poi plus the nuisance values
   ROOT.RooArgSet nullParams(*bModel.GetSnapshot())
   if bModel.GetNuisanceParameters()) nullParams.add(*bModel.GetNuisanceParameters():

   slrts.SetNullParameters(nullParams)
   ROOT.RooArgSet altParams(*sbModel.GetSnapshot())
   if sbModel.GetNuisanceParameters()) altParams.add(*sbModel.GetNuisanceParameters():
   slrts.SetAltParameters(altParams)


   ProfileLikelihoodTestStat profll = ProfileLikelihoodTestStat(*bModel.GetPdf())


   RatioOfProfiledLikelihoodsTestStat *
      ropl = RatioOfProfiledLikelihoodsTestStat(*bModel.GetPdf(), *sbModel.GetPdf(), sbModel.GetSnapshot())
   ropl.SetSubtractMLE(False)

   if testStatType == 3) profll.SetOneSidedDiscovery(1:
   profll.SetPrintLevel(printLevel)

   # profll.SetReuseNLL(mOptimize)
   # slrts.SetReuseNLL(mOptimize)
   # ropl.SetReuseNLL(mOptimize)

   AsymptoticCalculator.SetPrintLevel(printLevel)

   HypoTestCalculatorGeneric hypoCalc = 0
   # note here Null is B and Alt is S+B
   if calcType == 0) hypoCalc =  FrequentistCalculator(*data, *sbModel, *bModel:
   elif calcType == 1) hypoCalc=  HybridCalculator(*data, *sbModel, *bModel:
   elif calcType == 2) hypoCalc=  AsymptoticCalculator(*data, *sbModel, *bModel:

   if calcType == 0:
       ((FrequentistCalculator*)hypoCalc).SetToys(ntoys, ntoys/nToysRatio)
   if calcType == 1:
       ((HybridCalculator*)hypoCalc).SetToys(ntoys, ntoys/nToysRatio)
   if calcType == 2 :      if testStatType == 3) ((AsymptoticCalculator*) hypoCalc).SetOneSidedDiscovery(True:
      if testStatType != 2 and testStatType != 3:
         Warning("StandardHypoTestDemo", "Only the PL test statistic can be used with AsymptoticCalculator - use by default a two-sided PL")





   # check for nuisance prior pdf in case of nuisance parameters
   if calcType == 1 and (bModel.GetNuisanceParameters() or sbModel.GetNuisanceParameters() ):         ROOT.RooAbsPdf nuisPdf = 0
         if nuisPriorName) nuisPdf = w.pdf(nuisPriorName:
         # use prior defined first in bModel (then in SbModel)
         if not nuisPdf:            Info("StandardHypoTestDemo", "No nuisance pdf given for the HybridCalculator - try to deduce  pdf from the   model")
            if bModel.GetPdf() and bModel.GetObservables() :
               nuisPdf = ROOT.RooStats.MakeNuisancePdf(*bModel, "nuisancePdf_bmodel")
            else:
               nuisPdf = ROOT.RooStats.MakeNuisancePdf(*sbModel, "nuisancePdf_sbmodel")

         if not nuisPdf :            if bModel.GetPriorPdf():               nuisPdf = bModel.GetPriorPdf()
               Info("StandardHypoTestDemo", "No nuisance pdf given - try to use %s that is defined as a prior pdf in the B model",nuisPdf.GetName())

            else:
               Error("StandardHypoTestDemo", "Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified or can be derived")
               return


         assert(nuisPdf)
         Info("StandardHypoTestDemo", "Using as nuisance Pdf ... " )
         nuisPdf.Print()

          ROOT.RooArgSet nuisParams = (bModel.GetNuisanceParameters() ) ? bModel.GetNuisanceParameters() : sbModel.GetNuisanceParameters()
         ROOT.RooArgSet np = nuisPdf.getObservables(*nuisParams)
         if np.getSize() == 0:            Warning("StandardHypoTestDemo", "Prior nuisance does not depend on nuisance parameters. ROOT.They will be smeared in their full range")

         delete np

         ((HybridCalculator*)hypoCalc).ForcePriorNuisanceAlt(*nuisPdf)
         ((HybridCalculator*)hypoCalc).ForcePriorNuisanceNull(*nuisPdf)


   # hypoCalc.ForcePriorNuisanceAlt(*sbModel.GetPriorPdf())
   # hypoCalc.ForcePriorNuisanceNull(*bModel.GetPriorPdf())

   ROOT.ToyMCSampler sampler = (ToyMCSampler *)hypoCalc.GetTestStatSampler()

   if sampler and (calcType == 0 or calcType == 1) :
      # look if pdf is number counting or extended
      if sbModel.GetPdf().canBeExtended() :         if useNC)   Warning("StandardHypoTestDemo", "Pdf is extended: but number counting flag is set: ignore it ":

      else:
         # for not extended pdf
         if not useNC:            nEvents = data.numEntries()
            Info("StandardHypoTestDemo", "Pdf is not extended: number of events to generate taken  from observed data set is %d",nEvents)
            sampler.SetNEventsPerToy(nEvents)

         else:
            Info("StandardHypoTestDemo", "using a number counting pdf")
            sampler.SetNEventsPerToy(1)



      if data.isWeighted() and not generateBinned:         Info("StandardHypoTestDemo", "Data set is weighted, nentries = %d and sum weights = %8.1f but toy generation is unbinned - it would be faster to set generateBinned to ROOT.True\n",data.numEntries(), data.sumEntries())

      if generateBinned)  sampler.SetGenerateBinned(generateBinned:


      # set the test statistic
      if testStatType == 0) sampler.SetTestStatistic(slrts:
      if testStatType == 1) sampler.SetTestStatistic(ropl:
      if testStatType == 2 or testStatType == 3) sampler.SetTestStatistic(profll:



   HypoTestResult htr = hypoCalc.GetHypoTest()
   htr.SetPValueIsRightTail(True)
   htr.SetBackgroundAsAlt(False)
   htr.Print(); # how to get meaningfull CLs at self point?

   delete sampler
   delete slrts
   delete ropl
   delete profll

   if calcType != 2:      HypoTestPlot plot = HypoTestPlot(*htr,100)
      plot.SetLogYaxis(True)
      plot.Draw()

   else:
      std.cout << "Asymptotic results " << std.endl



   # look at expected significances
   # found median of S+B distribution
   if calcType != 2:
      SamplingDistribution altDist = htr.GetAltDistribution()
      HypoTestResult htExp("Expected Result")
      htExp.Append(htr)
      # find quantiles in alt (S+B) distribution
      double p[5]
      double q[5]
      for (i = 0; i < 5; ++i)         sig = -2  + i
         p[i] = ROOT.Math.normal_cdf(sig,1)

      std.vector<double> values = altDist.GetSamplingDistribution()
      ROOT.TMath.Quantiles( values.size(), 5, &values[0], q, p, False)

      for (i = 0; i < 5; ++i)         htExp.SetTestStatisticData( q[i] )
         sig = -2  + i
         std.cout << " Expected p -value and significance at " << sig << " sigma = "
                   << htExp.NullPValue() << " significance " << htExp.Significance() << " sigma " << std.endl



   else:
      # case of asymptotic calculator
      for (i = 0; i < 5; ++i)         sig = -2  + i
         # sigma is inverted here
         pval = AsymptoticCalculator.GetExpectedPValues( htr.NullPValue(), htr.AlternatePValue(), -sig, False)
         std.cout << " Expected p -value and significance at " << sig << " sigma = "
                   << pval << " significance " << ROOT.Math.normal_quantile_c(pval,1) << " sigma " << std.endl






