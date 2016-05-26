# Standard demo of the Feldman-Cousins calculator
'''
StandardFeldmanCousinsDemo

Author: Kyle Cranmer
date: Dec. 2010

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

The ROOT.RooStats.FeldmanCousins( tools is a classical frequentist calculation
based on the Neyman Construction.  ROOT.The test statistic can be
generalized for nuisance parameters by using the profile likeihood ratio.
But unlike the ROOT.RooStats.ProfileLikelihoodCalculator(, tool explicitly
builds the sampling distribution of the test statistic via toy Monte Carlo.
'''

#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TSystem.h"

#include "RooWorkspace.h"
#include "RooAbsData.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"


import ROOT
using namespace ROOT.RooStats

void StandardFeldmanCousinsDemo( infile = "",
                                 workspaceName = "combined",
                                 modelConfigName = "ModelConfig",
                                 dataName = "obsData")
  ##############################/
  # First part is just to access a user-defined file
  # or create the standard example file if it doesn't exist
  ##############################
   filename = ""
  if not strcmp(infile, ""):    filename = "results/example_combined_GaussExample_model.root"
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
  if not file :    cout <<"StandardRooStatsDemoMacro: Input file " << filename << " is not found" << endl
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


  ######################/
  # create and use the ROOT.RooStats.FeldmanCousins( tool
  # to find and plot the 95% confidence interval
  # on the parameter of interest as specified
  # in the model config
  ROOT.RooStats.FeldmanCousins( fc(*data,*mc)
  fc.SetConfidenceLevel(0.95); # 95% interval
  #fc.AdditionalNToysFactor(0.1); # to speed up the result
  fc.UseAdaptiveSampling(True); # speed it up a bit
  fc.SetNBins(10); # set how many points per parameter of interest to scan
  fc.CreateConfBelt(True); # save the information in the belt for plotting

  # Since self tool needs to throw toy MC the PDF needs to be
  # extended or the tool needs to know how many entries in a dataset
  # per pseudo experiment.
  # In the 'number counting form' where the entries in the dataset
  # are counts, not values of discriminating variables, the
  # datasets typically only have one entry and the PDF is not
  # extended.
  if not mc.GetPdf().canBeExtended():    if data.numEntries()==1:
      fc.FluctuateNumDataEntries(False)
    else:
      cout <<"Not sure what to do about self model" <<endl


  # We can use PROOF to speed things along in parallel
  #  ProofConfig pc(*w, 1, "workers=4", kFALSE)
  #  toymcsampler = (ToyMCSampler*) fc.GetTestStatSampler()
  #  toymcsampler.SetProofConfig(&pc); # enable proof


  # Now get the interval
  interval = fc.GetInterval()
  belt = fc.GetConfidenceBelt()

  # print out the iterval on the first Parameter of Interest
  firstPOI = (ROOT.RooRealVar*) mc.GetParametersOfInterest().first()
  cout << "\n95% interval on " <<firstPOI.GetName()<<" is : ["<<
    interval.LowerLimit(*firstPOI) << ", "<<
    interval.UpperLimit(*firstPOI) <<"] "<<endl

  #######################
  # No nice plots yet, plot the belt by hand

  # Ask the calculator which points were scanned
  parameterScan = (ROOT.RooDataSet*) fc.GetPointsToScan()
  ROOT.RooArgSet* tmpPoint

  # make a histogram of parameter vs. threshold
  histOfThresholds = ROOT.TH1F("histOfThresholds", "",
                                    parameterScan.numEntries(),
                                    firstPOI.getMin(),
                                    firstPOI.getMax())

  # loop through the points that were tested and ask confidence belt
  # what the upper/lower thresholds were.
  # For ROOT.RooStats.FeldmanCousins(, lower cut off is always 0
  for(Int_t i=0; i<parameterScan.numEntries(); ++i)    tmpPoint = (ROOT.RooArgSet*) parameterScan.get(i).clone("temp")
    arMax = belt.GetAcceptanceRegionMax(*tmpPoint)
    arMin = belt.GetAcceptanceRegionMax(*tmpPoint)
    poiVal = tmpPoint.getRealValue(firstPOI.GetName()) 
    histOfThresholds.Fill(poiVal,arMax)

  histOfThresholds.SetMinimum(0)
  histOfThresholds.Draw()


