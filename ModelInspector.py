###################
# ROOT.RooStats Model Inspector
# Author Kyle Cranmer <cranmer@cern.ch>
# Version 1, 2011
#   - based on tutorial macro by Bertrand Bellenot, Antcheva
# Version 2, 2011
#   - fixes from Bertrand Bellenot <Bertrand.Bellenot@cern.ch> for scrolling window for many parameters
#
#
# Usage:
# ROOT.The usage is the same as the StandardXxxDemo.C macros.
# ROOT.The macro expects a root file containing a workspace with a ROOT.RooStats.ModelConfig( and a dataset
# $ root
# .L ModelInspector.C+
# ModelInspector(fileName, workspaceName, modelConfigName, dataSetName)
#
# Drag the sliders to adjust the parameters of the model.
# the min and max range of the sliders are used to define the upper & lower variation
# the pointer position of the slider is the central blue curve.
#
# Click the FIT button to
#
# ROOT.To Do:
#  - check boxes to specify which nuisance parameters used in making variation
#  - a button to make the profile inspector plots
#  - a check button to use MINOS errors
#  - have fit button show the covariance matrix from the fit
#  - a button to make the log likelihood plots
#  - a dialog to open the desired file
#  - ability to see teh signal and background contributions?
#

#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGLabel.h"
#include "TGTripleSlider.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "TFile.h"
#include "RooArgSet.h"
#include "TList.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "TGButton.h"
#include <map>
#include "RooFitResult.h"
#include "TROOT.h"
#include <TApplication.h>
#include "RooSimultaneous.h"
#include "RooCategory.h"

enum ETestCommandIdentifiers   HId1,
   HId2,
   HId3,
   HCId1,
   HCId2,

   HSId1


import ROOT
using namespace ROOT.RooStats
class ModelInspectorGUI : public ROOT.TGMainFrame
private:
   ROOT.TRootEmbeddedCanvas *fCanvas
   ROOT.TGLayoutHints       *fLcan
   ROOT.TF1                 *fFitFcn
  ROOT.RooPlot* fPlot
  ROOT.RooWorkspace* fWS
  ROOT.TFile* fFile
  ROOT.RooStats.ModelConfig(* fMC
  ROOT.RooAbsData* fData
  ROOT.RooFitResult* fFitRes

  ROOT.TList fSliderList
  ROOT.TList fFrameList
  vector<RooPlot*> fPlotList
  map<TGTripleHSlider*, fSliderMap
  map<TGTripleHSlider*, fLabelMap

  ROOT.TGButton* fFitButton
  ROOT.TGTextButton *fExitButton

   # BB: a ROOT.TGCanvas and a vertical frame are needed for using scrollbars
   ROOT.TGCanvas *fCan
   ROOT.TGVerticalFrame *fVFrame

   ROOT.TGHorizontalFrame   *fHframe0, *fHframe1, *fHframe2
   ROOT.TGLayoutHints       *fBly, *fBfly1, *fBfly2, *fBfly3
   ROOT.TGTripleHSlider     *fHslider1
  #   ROOT.TGTextEntry         *fTeh1, *fTeh2, *fTeh3
   ROOT.TGTextBuffer        *fTbh1, *fTbh2, *fTbh3
   ROOT.TGCheckButton       *fCheck1, *fCheck2

public:
  ModelInspectorGUI(ROOT.RooWorkspace*, ROOT.RooStats.ModelConfig(*, ROOT.RooAbsData*)
   virtual ~ModelInspectorGUI()

   void CloseWindow()
   void DoText( char *text)
   void DoSlider()
   void DoSlider( char*)
   void DoFit()
   void DoExit()
   void HandleButtons()

   ClassDef(ModelInspectorGUI, 0)


#______________________________________________________________________________
ModelInspectorGUI.ModelInspectorGUI(ROOT.RooWorkspace* w, mc, data)
  : ROOT.TGMainFrame(gClient.GetRoot(), 100, 100)

  ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)
  fWS = w
  fMC = mc
  fData = data
  simPdf = NULL
  numCats = 1
  if strcmp(fMC.GetPdf().ClassName(), "RooSimultaneous")==0:    cout <<"Is a simultaneous PDF"<< endl
    simPdf = (ROOT.RooSimultaneous*)(fMC.GetPdf())
    channelCat = (ROOT.RooCategory*) (&simPdf.indexCat())
    cout <<" with " << channelCat.numTypes() << " categories"<<endl
    numCats = channelCat.numTypes()
  } else:
    cout <<"Is not a simultaneous PDF"<<endl

  fFitRes=0


  #char buf[32]
   SetCleanup(kDeepCleanup)

   # Create an embedded canvas and add to the main frame, in x and y
   # and with 30 pixel margins all around
   fCanvas = ROOT.TRootEmbeddedCanvas("Canvas", self, 600, 400)
   fLcan = ROOT.TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 10)
   AddFrame(fCanvas, fLcan)
   fPlotList.resize(numCats)
   if numCats>1:     fCanvas.GetCanvas().Divide(numCats)
     for(int i=0; i<numCats; ++i)       #   fCanvas.GetCanvas().SetFillColor(33)
       #   fCanvas.GetCanvas().SetFrameFillColor(41)
       fCanvas.GetCanvas().cd(i+1).SetBorderMode(0)
       fCanvas.GetCanvas().cd(i+1).SetGrid()
       #   fCanvas.GetCanvas().SetLogy()



   fHframe0 = ROOT.TGHorizontalFrame(self, 0, 0, 0)

   fCheck1 = ROOT.TGCheckButton(fHframe0, "&Constrained", HCId1)
   fCheck2 = ROOT.TGCheckButton(fHframe0, "&Relative", HCId2)
   fCheck1.SetState(kButtonUp)
   fCheck2.SetState(kButtonUp)
   fCheck1.SetToolTipText("Pointer position constrained to slider sides")
   fCheck2.SetToolTipText("Pointer position relative to slider position")

   fHframe0.Resize(200, 50)


   fHframe2 = ROOT.TGHorizontalFrame(self, 0, 0, 0)

   fFitButton = ROOT.TGTextButton(fHframe2, "Fit")
   fFitButton.Connect("Clicked()", "ModelInspectorGUI",self, "DoFit()")
   fExitButton = ROOT.TGTextButton(fHframe2, "Exit ")
   fExitButton.Connect("Clicked()", "ModelInspectorGUI", self, "DoExit()")


   fCheck1.Connect("Clicked()", "ModelInspectorGUI", self,
                    "HandleButtons()")
   fCheck2.Connect("Clicked()", "ModelInspectorGUI", self,
                    "HandleButtons()")

   fHframe2.Resize(100, 25)

   #--- layout for buttons: top align, expand horizontally
   fBly = ROOT.TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5)

   #--- layout for the frame: place at bottom, aligned
   fBfly1 = ROOT.TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 5, 5)
   fBfly2 = ROOT.TGLayoutHints(kLHintsTop | kLHintsLeft,    5, 5, 5, 5)
   fBfly3 = ROOT.TGLayoutHints(kLHintsTop | kLHintsRight,   5, 5, 5, 5)

   #   fHframe0.AddFrame(fCheck1, fBfly2)
   #   fHframe0.AddFrame(fCheck2, fBfly2)
   fHframe2.AddFrame(fFitButton, fBfly2)
   fHframe2.AddFrame(fExitButton, fBfly3)

   AddFrame(fHframe0, fBly)
   AddFrame(fHframe2, fBly)


   # Loop over POI & NP, slider
   # need maps of NP.slider? or just slider.NP
   ROOT.RooArgSet parameters
   parameters.add(*fMC.GetParametersOfInterest())
   parameters.add(*fMC.GetNuisanceParameters())
   it = parameters.createIterator()
   ROOT.RooRealVar* param=NULL

   # BB: ROOT.This is the part needed in order to have scrollbars
   fCan = ROOT.TGCanvas(self, 100, 100, kFixedSize)
   AddFrame(fCan, ROOT.TGLayoutHints(kLHintsExpandY | kLHintsExpandX))
   fVFrame = ROOT.TGVerticalFrame(fCan.GetViewPort(), 10, 10)
   fCan.SetContainer(fVFrame)
   # And that't it!
   # Obviously, parent of other subframes is now fVFrame instead of "self"...

   while( (param=(ROOT.RooRealVar*)it.Next()) )     cout <<"Adding Slider for "<<param.GetName()<<endl
     hframek = ROOT.TGHorizontalFrame(fVFrame, 0, 0, 0)

     hlabel = ROOT.TGLabel(hframek,Form("%s = %.3f +%.3f",param.GetName(),param.getVal(),param.getError()))
     hsliderk = ROOT.TGTripleHSlider(hframek, 190, kDoubleScaleBoth, HSId1,
                                                     kHorizontalFrame,
                                                     GetDefaultFrameBackground(),
                                                     kFALSE, kFALSE, kFALSE, kFALSE)
     hsliderk.Connect("PointerPositionChanged()", "ModelInspectorGUI",
                       self, "DoSlider()")
     hsliderk.Connect("PositionChanged()", "ModelInspectorGUI",
                       self, "DoSlider()")
     hsliderk.SetRange(param.getMin(),param.getMax())

     hframek.Resize(200, 25)
     fSliderList.Add(hsliderk)
     fFrameList.Add(hframek)

     hsliderk.SetPosition(param.getVal()-param.getError(),param.getVal()+param.getError())
     hsliderk.SetPointerPosition(param.getVal())

     hframek.AddFrame(hlabel, fBly)
     hframek.AddFrame(hsliderk, fBly)
     fVFrame.AddFrame(hframek, fBly)
     fSliderMap[hsliderk]=param.GetName()
     fLabelMap[hsliderk]=hlabel


   # Set main frame name, sub windows (buttons), layout
   # algorithm via Resize() and map main frame
   SetWindowName("RooFit/RooStats Model Inspector")
   MapSubwindows()
   Resize(GetDefaultSize())
   MapWindow()

   DoSlider()


#______________________________________________________________________________
ModelInspectorGUI.~ModelInspectorGUI()
   # Clean up

   Cleanup()


#______________________________________________________________________________
def CloseWindow(self):
   # Called when window is closed via the window manager.

   delete self


#______________________________________________________________________________
def DoText(self, * '''text'''):
   # Handle text entry widgets.

   ROOT.TGTextEntry *te = (TGTextEntry *) gTQSender
   id = te.WidgetId()

   switch (id)      case HId1:
         fHslider1.SetPosition(atof(fTbh1.GetString()),
                                fHslider1.GetMaxPosition())
         break
      case HId2:
         fHslider1.SetPointerPosition(atof(fTbh2.GetString()))
         break
      case HId3:
         fHslider1.SetPosition(fHslider1.GetMinPosition(),
                                atof(fTbh1.GetString()))
         break
      default:
         break

   DoSlider()


#______________________________________________________________________________
def DoFit(self):
  fFitRes = fMC.GetPdf().fitTo(*fData, ROOT.RooFit.Save())
  map<TGTripleHSlider*, it
  it = fSliderMap.begin()
  for(; it!=fSliderMap.end(); ++it)    ROOT.RooRealVar* param=fWS.var(it.second)
    param = (ROOT.RooRealVar*) fFitRes.floatParsFinal().find(it.second)
    it.first.SetPosition(param.getVal()-param.getError(),param.getVal()+param.getError())
    it.first.SetPointerPosition(param.getVal())

  DoSlider()



#______________________________________________________________________________
def DoSlider(self, text):
  cout << "." << text <<endl


#______________________________________________________________________________
def DoSlider(self):
   # Handle slider widgets.

   #char buf[32]

  simPdf = NULL
  numCats = 0
  if strcmp(fMC.GetPdf().ClassName(), "RooSimultaneous")==0:    simPdf = (ROOT.RooSimultaneous*)(fMC.GetPdf())
    channelCat = (ROOT.RooCategory*) (&simPdf.indexCat())
    numCats = channelCat.numTypes()
  } else:



  ######################/
  if not simPdf:    ######################/
    # if not SimPdf
    ######################/

    # pre loop
    map<TGTripleHSlider*, it
    delete fPlot
    fPlot = ((ROOT.RooRealVar*)fMC.GetObservables().first()).frame()
    fData.plotOn(fPlot)
    double normCount

    # high loop
    it = fSliderMap.begin()
    for(; it!=fSliderMap.end(); ++it)       name = it.second
      fWS.var(name).setVal(it.first.GetMaxPosition())
      param = fWS.var(name)
      fLabelMap[it.first].SetText(Form("%s = %.3f [%.3f,%.3f]",param.GetName(),it.first.GetPointerPosition(),it.first.GetMinPosition(),it.first.GetMaxPosition()))

    normCount = fMC.GetPdf().expectedEvents(*fMC.GetObservables())
    fMC.GetPdf().plotOn(fPlot, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent))

    # low loop
    it = fSliderMap.begin()
    for(; it!=fSliderMap.end(); ++it)       name = it.second
      fWS.var(name).setVal(it.first.GetMinPosition())

    normCount = fMC.GetPdf().expectedEvents(*fMC.GetObservables())
    fMC.GetPdf().plotOn(fPlot, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent))

    # central oop
    it = fSliderMap.begin()
    for(; it!=fSliderMap.end(); ++it)       name = it.second
      fWS.var(name).setVal(it.first.GetPointerPosition())

    normCount = fMC.GetPdf().expectedEvents(*fMC.GetObservables())
    fMC.GetPdf().plotOn(fPlot, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent))
    fPlot.Draw()

    fCanvas.GetCanvas().Modified()
    fCanvas.GetCanvas().Update()
    ######################################
  } else:
    ######################################
    # else (is simpdf)
    ######################################
    channelCat = (ROOT.RooCategory*) (&simPdf.indexCat())
    #    iter = simPdf.indexCat().typeIterator() 
    iter = channelCat.typeIterator() 
    tt = NULL
    frameIndex = 0
    while((tt=(ROOT.RooCatType*) iter.Next()))      ++frameIndex
      fCanvas.GetCanvas().cd(frameIndex)

      # pre loop
      pdftmp = simPdf.getPdf(tt.GetName()) 
      obstmp = pdftmp.getObservables(*fMC.GetObservables()) 
      obs = ((ROOT.RooRealVar*)obstmp.first())

      fPlot = fPlotList.at(frameIndex-1)
      if(fPlot) delete fPlot
      fPlot = obs.frame()
      fPlotList.at(frameIndex-1) = fPlot

      msglevel = ROOT.RooMsgService.instance().globalKillBelow()
      ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
      fData.plotOn(fPlot,MarkerSize(1),Cut(Form("%s==%s.%s",channelCat.GetName(),channelCat.GetName(),tt.GetName())), ROOT.RooFit.DataError(ROOT.RooAbsData.None))
      ROOT.RooMsgService.instance().setGlobalKillBelow(msglevel)


      map<TGTripleHSlider*, it
      double normCount

      # high loop
      it = fSliderMap.begin()
      for(; it!=fSliderMap.end(); ++it)          name = it.second
         fWS.var(name).setVal(it.first.GetMaxPosition())
         param = fWS.var(name)
         fLabelMap[it.first].SetText(Form("%s = %.3f [%.3f,%.3f]",param.GetName(),it.first.GetPointerPosition(),it.first.GetMinPosition(),it.first.GetMaxPosition()))

      normCount = pdftmp.expectedEvents(*obs)
      pdftmp.plotOn(fPlot, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineWidth(2.), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent)) 


      # low loop
      it = fSliderMap.begin()
      for(; it!=fSliderMap.end(); ++it)          name = it.second
         fWS.var(name).setVal(it.first.GetMinPosition())
         param = fWS.var(name)
         fLabelMap[it.first].SetText(Form("%s = %.3f [%.3f,%.3f]",param.GetName(),it.first.GetPointerPosition(),it.first.GetMinPosition(),it.first.GetMaxPosition()))

      normCount = pdftmp.expectedEvents(*obs)
      pdftmp.plotOn(fPlot, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineWidth(2.), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent)) 

      # central loop
      it = fSliderMap.begin()
      for(; it!=fSliderMap.end(); ++it)          name = it.second
         fWS.var(name).setVal(it.first.GetPointerPosition())
         param = fWS.var(name)
         fLabelMap[it.first].SetText(Form("%s = %.3f [%.3f,%.3f]",param.GetName(),it.first.GetPointerPosition(),it.first.GetMinPosition(),it.first.GetMaxPosition()))

      normCount = pdftmp.expectedEvents(*obs)
      if not fFitRes:
         pdftmp.plotOn(fPlot, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(2.), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent)) 
      else:
         pdftmp.plotOn(fPlot, ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.VisualizeError(*fFitRes,*fMC.GetNuisanceParameters()), ROOT.RooFit.FillColor(ROOT.kYellow)) 
         pdftmp.plotOn(fPlot, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineWidth(2.), ROOT.RooFit.Normalization(normCount, ROOT.RooAbsReal.NumEvent)) 
         msglevel = ROOT.RooMsgService.instance().globalKillBelow()
         ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
         fData.plotOn(fPlot,MarkerSize(1),Cut(Form("%s==%s.%s",channelCat.GetName(),channelCat.GetName(),tt.GetName())), ROOT.RooFit.DataError(ROOT.RooAbsData.None))
         ROOT.RooMsgService.instance().setGlobalKillBelow(msglevel)

      fPlot.Draw()

    fCanvas.GetCanvas().Modified()
    fCanvas.GetCanvas().Update()
    #####################/
    # end if simPdf:





#______________________________________________________________________________
def HandleButtons(self):
   # Handle different buttons.

   ROOT.TGButton *btn = (TGButton *) gTQSender
   id = btn.WidgetId()

   switch (id)      case HCId1:
         fHslider1.SetConstrained(fCheck1.GetState())
         break
      case HCId2:
         fHslider1.SetRelative(fCheck2.GetState())
         break
      default:
         break


def DoExit(self):
   printf("Exit application...")
   gApplication.Terminate(0)



void ModelInspector( infile = "",
                     workspaceName = "combined",
                     modelConfigName = "ModelConfig",
                     dataName = "obsData")
#ifdef __CINT__
  cout <<"You must use ACLIC for self.  Use ModelInspector.C+"<<endl
  return
#endif

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


  ModelInspectorGUI(w,mc,data)


