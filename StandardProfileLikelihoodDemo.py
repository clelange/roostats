# Standard demo of the Profile Likelihood calculcator
'''
StandardProfileLikelihoodDemo

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

The ROOT.RooStats.ProfileLikelihoodCalculator( is based on Wilks's theorem
and the asymptotic properties of the profile likeihood ratio
(eg. that it is chi-square distributed for the ROOT.True value).
'''

#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

import ROOT
using namespace ROOT.RooStats

void StandardProfileLikelihoodDemo( infile = "",
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


  ######################/
  # create and use the ROOT.RooStats.ProfileLikelihoodCalculator(
  # to find and plot the 95% confidence interval
  # on the parameter of interest as specified
  # in the model config
  ROOT.RooStats.ProfileLikelihoodCalculator( pl(*data,*mc)
  pl.SetConfidenceLevel(0.95); # 95% interval
  interval = pl.GetInterval()

  # print out the iterval on the first Parameter of Interest
  firstPOI = (ROOT.RooRealVar*) mc.GetParametersOfInterest().first()
  cout << "\n95% interval on " <<firstPOI.GetName()<<" is : ["<<
    interval.LowerLimit(*firstPOI) << ", "<<
    interval.UpperLimit(*firstPOI) <<"] "<<endl

  # make a plot

  cout << "making a plot of the profile likelihood function ....(if it is taking a lot of time use less points or the ROOT.TF1 drawing option)\n"
  ROOT.RooStats.LikelihoodIntervalPlot( plot(interval)
  plot.SetNPoints(50);  # do not use too many points, could become very slow for some models
  plot.Draw("");  # use option ROOT.TF1 if too slow (plot.Draw("tf1")



