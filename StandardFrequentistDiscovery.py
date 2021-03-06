# StandardFrequentistDiscovery

'''
 StandardFrequentistDiscovery

 Author: Sven Kreiss, Cranmer
 date: May 2012

 ROOT.This is a standard demo that can be used with any ROOT file
 prepared in the standard way.  You specify:
 - name for input ROOT file
 - name of workspace inside ROOT file that holds model and data
 - name of ROOT.RooStats.ModelConfig( that specifies details for calculator tools
 - name of dataset

 With default parameters the macro will attempt to run the
 standard hist2workspace example and read the ROOT file
 that it produces.
 '''

#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStopwatch.h"

#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRandom.h"
#include "RooRealSumPdf.h"
#include "RooNumIntConfig.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCImportanceSampler.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/FrequentistCalculator.h"
#include "TSystem.h"

#include <vector>

import ROOT
using namespace ROOT.RooStats




double StandardFrequentistDiscovery(
    infile = "",
    workspaceName = "channel1",
    modelConfigNameSB = "ModelConfig",
    dataName = "obsData",
   toys = 1000,
   poiValueForBackground = 0.0,
   poiValueForSignal = 1.0
)
   # ROOT.The workspace contains the model for s+b. ROOT.The b model is "autogenerated"
   # by copying s+b and setting the one parameter of interest to zero.
   # ROOT.To keep the script simple, parameters of interest or different
   # functional forms of the b model are not supported.

   # for now, there is only one parameter of interest, these are
   # its values:

   ##############################/
   # First part is just to access a user-defined file
   # or create the standard example file if it doesn't exist
   ##############################
    filename = ""
   if not strcmp(infile, ""):      filename = "results/example_channel1_GammaExample_model.root"
      fileExist = not gSystem.AccessPathName(filename); # note opposite return code
      # if file does not exists generate with histfactory
      if not fileExist:#ifdef _WIN32
         cout << "HistFactory file cannot be generated on Windows - exit" << endl
         return -1
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
      return -1



   ##############################/
   # ROOT.Tutorial starts here
   ##############################

   ROOT.TStopwatch *mn_t = ROOT.TStopwatch
   mn_t.Start()

   # get the workspace out of the file
   w = (ROOT.RooWorkspace*) file.Get(workspaceName)
   if not w:      cout << "workspace not found" << endl
      return -1.0


   # get the modelConfig out of the file
   mc = (ModelConfig*) w.obj(modelConfigNameSB)

   # get the data out of the file
   data = w.data(dataName)

   # make sure ingredients are found
   if not data or not mc:      w.Print()
      cout << "data or ROOT.RooStats.ModelConfig( was not found" << endl
      return -1.0



   firstPOI = (ROOT.RooRealVar*) mc.GetParametersOfInterest().first()
   firstPOI.setVal(poiValueForSignal)
   mc.SetSnapshot(*mc.GetParametersOfInterest())
   # create null model
   ROOT.RooStats.ModelConfig( *mcNull = mc.Clone("ModelConfigNull")
   firstPOI.setVal(poiValueForBackground)
   mcNull.SetSnapshot(*(ROOT.RooArgSet*)mcNull.GetParametersOfInterest().snapshot())



   # ----------------------------------------------------
   # Configure a ProfileLikelihoodTestStat and a SimpleLikelihoodRatioTestStat
   # to use simultaneously with ROOT.ToyMCSampler
   plts =  ProfileLikelihoodTestStat(*mc.GetPdf())
   plts.SetOneSidedDiscovery(True)
   plts.SetVarName( "q_{0}/2" )

   # ----------------------------------------------------
   # configure the ROOT.ToyMCImportanceSampler with two test statistics
   ROOT.ToyMCSampler toymcs(*plts, 50)



   # Since self tool needs to throw toy MC the PDF needs to be
   # extended or the tool needs to know how many entries in a dataset
   # per pseudo experiment.
   # In the 'number counting form' where the entries in the dataset
   # are counts, not values of discriminating variables, the
   # datasets typically only have one entry and the PDF is not
   # extended.
   if not mc.GetPdf().canBeExtended():      if data.numEntries() == 1:         toymcs.SetNEventsPerToy(1)
      } else cout << "Not sure what to do about self model" << endl


   # We can use PROOF to speed things along in parallel
   # ProofConfig pc(*w, 2, "user@yourfavoriteproofcluster", False)
   ProofConfig pc(*w, 2, "", False)
   #toymcs.SetProofConfig(&pc);    # enable proof


   # instantiate the calculator
   FrequentistCalculator freqCalc(*data, *mc, *mcNull, &toymcs)
   freqCalc.SetToys( toys,toys ); # null toys, toys

   # Run the calculator and print result
   freqCalcResult = freqCalc.GetHypoTest()
   freqCalcResult.GetNullDistribution().SetTitle( "b only" )
   freqCalcResult.GetAltDistribution().SetTitle( "s+b" )
   freqCalcResult.Print()
   pvalue = freqCalcResult.NullPValue()

   # stop timing
   mn_t.Stop()
   cout << "total CPU time: " << mn_t.CpuTime() << endl
   cout << "total real time: " << mn_t.RealTime() << endl

   # plot
   c1 = ROOT.TCanvas()
   HypoTestPlot *plot = HypoTestPlot(*freqCalcResult, 100, -0.49, 9.51 )
   plot.SetLogYaxis(True)

   # add chi2 to plot
   nPOI = 1
   f = ROOT.TF1("f", ROOT.TString.Format("1*ROOT.Math.chisquared_pdf(2*x,%d,0)",nPOI), 0,20)
   f.SetLineColor( kBlack )
   f.SetLineStyle( 7 )
   plot.AddTF1( f, ROOT.TString.Format("#chi^{2}(2x,%d)",nPOI) )

   plot.Draw()
   c1.SaveAs("standard_discovery_output.pdf")


   return pvalue



