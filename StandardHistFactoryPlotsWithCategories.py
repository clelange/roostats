'''
 StandardHistFactoryPlotsWithCategories

 Author: Kyle Cranmer
 date: Spring. 2011

 ROOT.This is a standard demo that can be used with any ROOT file
 prepared in the standard way.  You specify:
 - name for input ROOT file
 - name of workspace inside ROOT file that holds model and data
 - name of ROOT.RooStats.ModelConfig( that specifies details for calculator tools
 - name of dataset

 With default parameters the macro will attempt to run the
 standard hist2workspace example and read the ROOT file
 that it produces.

 ROOT.The macro will scan through all the categories in a simPdf find the corresponding
 observable.  For each cateogry, will loop through each of the nuisance parameters
 and plot
 - the data
 - the nominal model (blue)
 - the +Nsigma (red)
 - the -Nsigma (green)

 You can specify how many sigma to vary by changing nSigmaToVary.
 You can also change the signal rate by changing muVal.

 ROOT.The script produces a lot plots, can merge them by doing:
 gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged.pdf `ls *pdf`
 '''

#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileInspector.h"

import ROOT
using namespace ROOT.RooStats
using namespace std

void StandardHistFactoryPlotsWithCategories( infile = "",
                                             workspaceName = "combined",
                                             modelConfigName = "ModelConfig",
                                             dataName = "obsData")

   double nSigmaToVary=5.
   double muVal=0
   bool doFit=False

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
   if not w:      cout <<"workspace not found" << endl
      return


   # get the modelConfig out of the file
   mc = (ModelConfig*) w.obj(modelConfigName)

   # get the modelConfig out of the file
   data = w.data(dataName)

   # make sure ingredients are found
   if not data or not mc:      w.Print()
      cout << "data or ROOT.RooStats.ModelConfig( was not found" <<endl
      return


   #######################
   # now use the profile inspector

   obs = (ROOT.RooRealVar*)mc.GetObservables().first()
   list = ROOT.TList()


   ROOT.RooRealVar firstPOI = dynamic_cast<RooRealVar*>(mc.GetParametersOfInterest().first())

   firstPOI.setVal(muVal)
   #  firstPOI.setConstant()
   if doFit:      mc.GetPdf().fitTo(*data)


   ####################
   ####################
   ####################

   mc.GetNuisanceParameters().Print("v")
   nPlotsMax = 1000
   cout <<" check expectedData by category"<<endl
   ROOT.RooDataSet* simData=NULL
   simPdf = NULL
   if strcmp(mc.GetPdf().ClassName(), "RooSimultaneous")==0:      cout <<"Is a simultaneous PDF"<<endl
      simPdf = (ROOT.RooSimultaneous *)(mc.GetPdf())
   } else:
      cout <<"Is not a simultaneous PDF"<<endl




   if doFit:      channelCat = (ROOT.RooCategory*) (&simPdf.indexCat())
      iter = channelCat.typeIterator() 
      tt = NULL
      tt=(ROOT.RooCatType*) iter.Next()
      pdftmp = ((ROOT.RooSimultaneous*)mc.GetPdf()).getPdf(tt.GetName()) 
      obstmp = pdftmp.getObservables(*mc.GetObservables()) 
      obs = ((ROOT.RooRealVar*)obstmp.first())
      frame = obs.frame()
      cout <<Form("%s==%s.%s",channelCat.GetName(),channelCat.GetName(),tt.GetName())<<endl
      cout << tt.GetName() << " " << channelCat.getLabel() <<endl
      data.plotOn(frame,MarkerSize(1),Cut(Form("%s==%s.%s",channelCat.GetName(),channelCat.GetName(),tt.GetName())), ROOT.RooFit.DataError(ROOT.RooAbsData.None))

      normCount = data.sumEntries(Form("%s==%s.%s",channelCat.GetName(),channelCat.GetName(),tt.GetName())) 

      pdftmp.plotOn(frame, ROOT.RooFit.LineWidth(2.), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent)) 
      frame.Draw()
      cout <<"events = " << mc.GetPdf().expectedEvents(*data.get()) <<endl
      return




   int nPlots=0
   if not simPdf:
      it = mc.GetNuisanceParameters().createIterator()
      var = NULL
      while( (var = (ROOT.RooRealVar*) it.Next()) != NULL)         frame = obs.frame()
         frame.SetYTitle(var.GetName())
         data.plotOn(frame,MarkerSize(1))
         var.setVal(0)
         mc.GetPdf().plotOn(frame, ROOT.RooFit.LineWidth(1.))
         var.setVal(1)
         mc.GetPdf().plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineWidth(1))
         var.setVal(-1)
         mc.GetPdf().plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineWidth(1))
         list.Add(frame)
         var.setVal(0)



   } else:
      channelCat = (ROOT.RooCategory*) (&simPdf.indexCat())
      #    iter = simPdf.indexCat().typeIterator() 
      iter = channelCat.typeIterator() 
      tt = NULL
      while(nPlots<nPlotsMax and (tt=(ROOT.RooCatType*) iter.Next()))
         cout << "on type " << tt.GetName() << " " << endl
         # Get pdf associated with state from simpdf
         pdftmp = simPdf.getPdf(tt.GetName()) 

         # Generate observables defined by the pdf associated with self state
         obstmp = pdftmp.getObservables(*mc.GetObservables()) 
         #      obstmp.Print()


         obs = ((ROOT.RooRealVar*)obstmp.first())

         it = mc.GetNuisanceParameters().createIterator()
         var = NULL
         while(nPlots<nPlotsMax and (var = (ROOT.RooRealVar*) it.Next()))            c2 = ROOT.TCanvas("c2")
            frame = obs.frame()
            frame.SetName(Form("frame%d",nPlots))
            frame.SetYTitle(var.GetName())

            cout <<Form("%s==%s.%s",channelCat.GetName(),channelCat.GetName(),tt.GetName())<<endl
            cout << tt.GetName() << " " << channelCat.getLabel() <<endl
            data.plotOn(frame,MarkerSize(1),Cut(Form("%s==%s.%s",channelCat.GetName(),channelCat.GetName(),tt.GetName())), ROOT.RooFit.DataError(ROOT.RooAbsData.None))

            normCount = data.sumEntries(Form("%s==%s.%s",channelCat.GetName(),channelCat.GetName(),tt.GetName())) 

            if strcmp(var.GetName(), "Lumi")==0:               cout <<"working on lumi"<<endl
               var.setVal(w.var("nominalLumi").getVal())
               var.Print()
            } else:
               var.setVal(0)

            # w.allVars().Print("v")
            # mc.GetNuisanceParameters().Print("v")
            # pdftmp.plotOn(frame, ROOT.RooFit.LineWidth(2.))
            # mc.GetPdf().plotOn(frame, ROOT.RooFit.LineWidth(2.),Slice(*channelCat,tt.GetName()),ProjWData(*data))
            #pdftmp.plotOn(frame, ROOT.RooFit.LineWidth(2.),Slice(*channelCat,tt.GetName()),ProjWData(*data))
            normCount = pdftmp.expectedEvents(*obs)
            pdftmp.plotOn(frame, ROOT.RooFit.LineWidth(2.), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent)) 

            if strcmp(var.GetName(), "Lumi")==0:               cout <<"working on lumi"<<endl
               var.setVal(w.var("nominalLumi").getVal()+0.05)
               var.Print()
            } else:
               var.setVal(nSigmaToVary)

            # pdftmp.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineWidth(2))
            # mc.GetPdf().plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineWidth(2.),Slice(*channelCat,tt.GetName()),ProjWData(*data))
            #pdftmp.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineWidth(2.),Slice(*channelCat,tt.GetName()),ProjWData(*data))
            normCount = pdftmp.expectedEvents(*obs)
            pdftmp.plotOn(frame, ROOT.RooFit.LineWidth(2.), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent)) 

            if strcmp(var.GetName(), "Lumi")==0:               cout <<"working on lumi"<<endl
               var.setVal(w.var("nominalLumi").getVal()-0.05)
               var.Print()
            } else:
               var.setVal(-nSigmaToVary)

            # pdftmp.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineWidth(2))
            # mc.GetPdf().plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineWidth(2),Slice(*channelCat,tt.GetName()),ProjWData(*data))
            #pdftmp.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineWidth(2),Slice(*channelCat,tt.GetName()),ProjWData(*data))
            normCount = pdftmp.expectedEvents(*obs)
            pdftmp.plotOn(frame, ROOT.RooFit.LineWidth(2.), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent)) 



            # set them back to normal
            if strcmp(var.GetName(), "Lumi")==0:               cout <<"working on lumi"<<endl
               var.setVal(w.var("nominalLumi").getVal())
               var.Print()
            } else:
               var.setVal(0)


            list.Add(frame)

            # quit making plots
            ++nPlots

            frame.Draw()
            c2.SaveAs(Form("%s_%s_%s.pdf",tt.GetName(),obs.GetName(),var.GetName()))
            delete c2






   ####################
   ####################
   ####################

   # now make plots
   c1 = ROOT.TCanvas("c1", "ProfileInspectorDemo",800,200)
   if list.GetSize()>4:      n = list.GetSize()
      nx = (int)sqrt(n) 
      ny = ROOT.TMath.CeilNint(n/nx)
      nx = ROOT.TMath.CeilNint( sqrt(n) )
      c1.Divide(ny,nx)
   } else:
      c1.Divide(list.GetSize())
   for(int i=0; i<list.GetSize(); ++i)      c1.cd(i+1)
      list.At(i).Draw()







