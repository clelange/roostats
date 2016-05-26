# Standard demo of the ProfileInspector class
'''
StandardProfileInspectorDemo

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

The ProfileInspector plots the conditional maximum likelihood estimate
of each nuisance parameter in the model vs. the parameter of interest.
(aka. profiled value of nuisance parameter vs. parameter of interest)
(aka. best fit nuisance parameter with p.o.i fixed vs. parameter of interest)

'''

#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileInspector.h"

import ROOT
using namespace ROOT.RooStats

void StandardProfileInspectorDemo( infile = "",
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


  #######################
  # now use the profile inspector
  ProfileInspector p
  list = p.GetListOfProfilePlots(*data,mc)

  # now make plots
  c1 = ROOT.TCanvas("c1", "ProfileInspectorDemo",800,200)
  if list.GetSize()>4:    n = list.GetSize()
    nx = (int)sqrt(n) 
    ny = ROOT.TMath.CeilNint(n/nx)
    nx = ROOT.TMath.CeilNint( sqrt(n) )
    c1.Divide(ny,nx)
  } else:
    c1.Divide(list.GetSize())
  for(int i=0; i<list.GetSize(); ++i)    c1.cd(i+1)
    list.At(i).Draw("al")


  cout << endl

