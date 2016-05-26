''' -*- mode: c++ -*- '''
# Standard tutorial macro for performing an inverted  hypothesis test for computing an interval
#
# ROOT.This macro will perform a scan of the p-values for computing the interval or limit
#
#Author:  L. Moneta
#
# Usage:
#
# root>.L StandardHypoTestInvDemo.C
# root> StandardHypoTestInvDemo("fileName", "workspace name", "S+B modelconfig name", "B model name", "data set name", type, statistic type, CLS,
#                                number of points, xmin, xmax, of toys, number counting)
#
#
# type = 0 Freq calculator
# type = 1 Hybrid calculator
# type = 2 Asymptotic calculator
# type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)
#
# testStatType = 0 LEP
#              = 1 ROOT.Tevatron
#              = 2 Profile Likelihood two sided
#              = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)
#              = 4 Profile Likelihood signed ( pll = -pll if mu < mu_hat)
#              = 5 Max Likelihood Estimate as test statistic
#              = 6 Number of observed event as test statistic
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
#include "TROOT.h"
#include "TSystem.h"

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
#include "RooStats/NumEventsTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

import ROOT
using namespace ROOT.RooStats
using namespace std

plotHypoTestResult = ROOT.True;          # plot test statistic result at each point
writeResult = ROOT.True;                 # write HypoTestInverterResult in a file
TString resultFileName;                  # file with results (by default is built automatically using the workspace input file name)
optimize = ROOT.True;                    # optmize evaluation of test statistic
useVectorStore = ROOT.True;              # convert data to use roofit data store
generateBinned = False;             # generate binned data sets
noSystematics = False;              # force all systematics to be off (i.e. set all nuisance parameters as constat
                                         # to their nominal values)
nToysRatio = 2;                   # ratio Ntoys S+b/ntoysB
maxPOI = -1;                      # max value used of POI (in case of auto scan)
useProof = False;                   # use Proof Lite when using toys (for freq or hybrid)
nworkers = 0;                        # number of worker for ProofLite (default use all available cores)
enableDetailedOutput = False;       # enable detailed output with all fit information for each toys (output will be written in result file)
rebuild = False;                    # re-do extra toys for computing expected limits and rebuild test stat
                                         # distributions (N.B self requires much more CPU (factor is equivalent to nToyToRebuild)
nToyToRebuild = 100;                 # number of toys used to rebuild
int rebuildParamValues=0;                # = 0   do a profile of all the parameters on the B (alt snapshot) before performing a rebuild operation (default)
                                         # = 1   use initial workspace parameters with B snapshot values
                                         # = 2   use all initial workspace parameters with B
                                         # Otherwise the rebuild will be performed using
initialFit = -1;                     # do a first  fit to the model (-1 : default, skip fit, do always fit)
randomSeed = -1;                     # random seed (if = -1: use default value, if = 0 always random )
                                         # NOTE: Proof uses automatically a random seed

nAsimovBins = 0;                     # number of bins in observables used for Asimov data sets (0 is the default and it is given by workspace, is 100)

reuseAltToys = False;                # reuse same toys for alternate hypothesis (if set one gets more stable bands)
confidenceLevel = 0.95;            # confidence level value



massValue = "";              # extra string to tag output file of result
minimizerType = "";                  # minimizer type (default is what is in ROOT.Math.MinimizerOptions.DefaultMinimizerType()
printLevel = 0;                    # print level for debugging PL test statistics and calculators

useNLLOffset = False;               # use NLL offset when fitting (self increase stability of fits) 


# internal class to run the inverter and more

namespace ROOT.RooStats
   class HypoTestInvTool
   public:
      HypoTestInvTool()
      ~HypoTestInvTool(){

      HypoTestInverterResult *
      RunInverter(ROOT.RooWorkspace * w,
                   char * modelSBName, * modelBName,
                   char * dataName,
                  int type, testStatType,
                  bool useCLs,
                  int npoints, poimin, poimax, ntoys,
                  useNumberCounting = False,
                   char nuisPriorName = 0)



      void
      AnalyzeResult( HypoTestInverterResult * r,
                     int calculatorType,
                     int testStatType,
                     bool useCLs,
                     int npoints,
                      char fileNameBase = 0 )

      void SetParameter( char * name, * value)
      void SetParameter( char * name, value)
      void SetParameter( char * name, value)
      void SetParameter( char * name, value)

   private:

      bool mPlotHypoTestResult
      bool mWriteResult
      bool mOptimize
      bool mUseVectorStore
      bool mGenerateBinned
      bool mUseProof
      bool mRebuild
      bool mReuseAltToys
      bool mEnableDetOutput
      int     mNWorkers
      int     mNToyToRebuild
      int     mRebuildParamValues
      int     mPrintLevel
      int     mInitialFit
      int     mRandomSeed
      double  mNToysRatio
      double  mMaxPoi
      int mAsimovBins
      std.string mMassValue
      std.string mMinimizerType;                  # minimizer type (default is what is in ROOT.Math.MinimizerOptions.DefaultMinimizerType()
      ROOT.TString     mResultFileName


} # end namespace ROOT.RooStats

RooStats.HypoTestInvTool.HypoTestInvTool() : mPlotHypoTestResult(True),
                                               mWriteResult(False),
                                               mOptimize(True),
                                               mUseVectorStore(True),
                                               mGenerateBinned(False),
                                               mUseProof(False),
                                               mEnableDetOutput(False),
                                               mRebuild(False),
                                               mReuseAltToys(False),
                                               mNWorkers(4),
                                               mNToyToRebuild(100),
                                               mRebuildParamValues(0),
                                               mPrintLevel(0),
                                               mInitialFit(-1),
                                               mRandomSeed(-1),
                                               mNToysRatio(2),
                                               mMaxPoi(-1),
                                               mAsimovBins(0),
                                               mMassValue(""),
                                               mMinimizerType(""),
                                               mResultFileName()



void
RooStats.HypoTestInvTool.SetParameter( char * name, value)   #
   # set boolean parameters
   #

   std.string s_name(name)

   if (s_name.find("PlotHypoTestResult") != std.string.npos) mPlotHypoTestResult = value
   if (s_name.find("WriteResult") != std.string.npos) mWriteResult = value
   if (s_name.find("Optimize") != std.string.npos) mOptimize = value
   if (s_name.find("UseVectorStore") != std.string.npos) mUseVectorStore = value
   if (s_name.find("GenerateBinned") != std.string.npos) mGenerateBinned = value
   if (s_name.find("UseProof") != std.string.npos) mUseProof = value
   if (s_name.find("EnableDetailedOutput") != std.string.npos) mEnableDetOutput = value
   if (s_name.find("Rebuild") != std.string.npos) mRebuild = value
   if (s_name.find("ReuseAltToys") != std.string.npos) mReuseAltToys = value

   return




void
RooStats.HypoTestInvTool.SetParameter( char * name, value)   #
   # set integer parameters
   #

   std.string s_name(name)

   if (s_name.find("NWorkers") != std.string.npos) mNWorkers = value
   if (s_name.find("NToyToRebuild") != std.string.npos) mNToyToRebuild = value
   if (s_name.find("RebuildParamValues") != std.string.npos) mRebuildParamValues = value
   if (s_name.find("PrintLevel") != std.string.npos) mPrintLevel = value
   if (s_name.find("InitialFit") != std.string.npos) mInitialFit = value
   if (s_name.find("RandomSeed") != std.string.npos) mRandomSeed = value
   if (s_name.find("AsimovBins") != std.string.npos) mAsimovBins = value

   return




void
RooStats.HypoTestInvTool.SetParameter( char * name, value)   #
   # set double precision parameters
   #

   std.string s_name(name)

   if (s_name.find("NToysRatio") != std.string.npos) mNToysRatio = value
   if (s_name.find("MaxPOI") != std.string.npos) mMaxPoi = value

   return




void
RooStats.HypoTestInvTool.SetParameter( char * name, * value)   #
   # set string parameters
   #

   std.string s_name(name)

   if s_name.find("MassValue") != std.string.npos) mMassValue.assign(value:
   if s_name.find("MinimizerType") != std.string.npos) mMinimizerType.assign(value:
   if (s_name.find("ResultFileName") != std.string.npos) mResultFileName = value

   return




void
StandardHypoTestInvDemo( char infile = 0,
                         char wsName = "combined",
                         char modelSBName = "ModelConfig",
                         char modelBName = "",
                         char dataName = "obsData",
                        calculatorType = 0,
                        testStatType = 0,
                        useCLs = ROOT.True ,
                        npoints = 6,
                        poimin = 0,
                        poimax = 5,
                        int ntoys=1000,
                        useNumberCounting = False,
                         char nuisPriorName = 0)'''

  Other Parameter to pass in tutorial
  apart from standard for filename, ws, and data

  type = 0 Freq calculator
  type = 1 Hybrid calculator
  type = 2 Asymptotic calculator
  type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)

  testStatType = 0 LEP
  = 1 ROOT.Tevatron
  = 2 Profile Likelihood
  = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)
  = 4 Profiel Likelihood signed ( pll = -pll if mu < mu_hat)
  = 5 Max Likelihood Estimate as test statistic
  = 6 Number of observed event as test statistic

  useCLs          scan for CLs (otherwise for CLs+b)

  npoints:        number of points to scan , autoscan npoints = -1

  poimin, min/max value to scan in case of fixed scans
  (if min >  max, to find automatically)

  ntoys:         number of toys to use

  useNumberCounting:  set to ROOT.True when using number counting events

  nuisPriorName:   name of prior for the nnuisance. ROOT.This is often expressed as constraint term in the global model
  It is needed only when using the HybridCalculator (type=1)
  If not given by default the prior pdf from ROOT.RooStats.ModelConfig( is used.

  extra options are available as global paramwters of the macro. ROOT.They major ones are:

  plotHypoTestResult   plot result of tests at each point (TS distributions) (defauly is ROOT.True)
  useProof             use Proof   (default is ROOT.True)
  writeResult          write result of scan (default is ROOT.True)
  rebuild              rebuild scan for expected limits (require extra toys) (default is False)
  generateBinned       generate binned data sets for toys (default is False) - be careful not to activate with
  a too large (>=3) number of observables
  nToyRatio            ratio of S+B/B toys (default is 2)


'''



   ROOT.TString filename(infile)
   if filename.IsNull():      filename = "results/example_combined_GaussExample_model.root"
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




   HypoTestInvTool calc

   # set parameters
   calc.SetParameter("PlotHypoTestResult", plotHypoTestResult)
   calc.SetParameter("WriteResult", writeResult)
   calc.SetParameter("Optimize", optimize)
   calc.SetParameter("UseVectorStore", useVectorStore)
   calc.SetParameter("GenerateBinned", generateBinned)
   calc.SetParameter("NToysRatio", nToysRatio)
   calc.SetParameter("MaxPOI", maxPOI)
   calc.SetParameter("UseProof", useProof)
   calc.SetParameter("EnableDetailedOutput", enableDetailedOutput)
   calc.SetParameter("NWorkers", nworkers)
   calc.SetParameter("Rebuild", rebuild)
   calc.SetParameter("ReuseAltToys", reuseAltToys)
   calc.SetParameter("NToyToRebuild", nToyToRebuild)
   calc.SetParameter("RebuildParamValues", rebuildParamValues)
   calc.SetParameter("MassValue", massValue.c_str())
   calc.SetParameter("MinimizerType", minimizerType.c_str())
   calc.SetParameter("PrintLevel", printLevel)
   calc.SetParameter("InitialFit",initialFit)
   calc.SetParameter("ResultFileName",resultFileName)
   calc.SetParameter("RandomSeed",randomSeed)
   calc.SetParameter("AsimovBins",nAsimovBins)

   # enable offset for all roostats
   if (useNLLOffset) ROOT.RooStats.UseNLLOffset(True); 

   ROOT.RooWorkspace w = dynamic_cast<RooWorkspace*>( file.Get(wsName) )
   HypoTestInverterResult r = 0
   std.cout << w << "\t" << filename << std.endl
   if w != NULL:      r = calc.RunInverter(w, modelSBName, modelBName,
                           dataName, calculatorType, testStatType, useCLs,
                           npoints, poimin, poimax,
                           ntoys, useNumberCounting, nuisPriorName )
      if not r:         std.cerr << "Error running the HypoTestInverter - Exit " << std.endl
         return


   else:
      # case workspace is not present look for the inverter result
      std.cout << "Reading an HypoTestInverterResult with name " << wsName << " from file " << filename << std.endl
      r = dynamic_cast<HypoTestInverterResult*>( file.Get(wsName) ); #
      if not r:         std.cerr << "File " << filename << " does not contain a workspace or an HypoTestInverterResult - Exit "
                   << std.endl
         file.ls()
         return



   calc.AnalyzeResult( r, calculatorType, testStatType, useCLs, npoints, infile )

   return




void
RooStats.HypoTestInvTool.AnalyzeResult( HypoTestInverterResult * r,
                                          int calculatorType,
                                          int testStatType,
                                          bool useCLs,
                                          int npoints,
                                           char * fileNameBase )
   # analyze result produced by the inverter, save it in a file


   lowerLimit = 0
   llError = 0
#if defined ROOT_SVN_VERSION and  ROOT_SVN_VERSION >= 44126
   if r.IsTwoSided():      lowerLimit = r.LowerLimit()
      llError = r.LowerLimitEstimatedError()

#else:
   lowerLimit = r.LowerLimit()
   llError = r.LowerLimitEstimatedError()
#endif

   upperLimit = r.UpperLimit()
   ulError = r.UpperLimitEstimatedError()

   #std.cout << "DEBUG : [ " << lowerLimit << " , " << upperLimit << "  ] " << std.endl

   if lowerLimit < upperLimit*(1.- 1.E-4) and lowerLimit != 0:
      std.cout << "The computed lower limit is: " << lowerLimit << " +/- " << llError << std.endl
   std.cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std.endl


   # compute expected limit
   std.cout << "Expected upper limits, the B (alternate) model : " << std.endl
   std.cout << " expected limit (median) " << r.GetExpectedUpperLimit(0) << std.endl
   std.cout << " expected limit (-1 sig) " << r.GetExpectedUpperLimit(-1) << std.endl
   std.cout << " expected limit (+1 sig) " << r.GetExpectedUpperLimit(1) << std.endl
   std.cout << " expected limit (-2 sig) " << r.GetExpectedUpperLimit(-2) << std.endl
   std.cout << " expected limit (+2 sig) " << r.GetExpectedUpperLimit(2) << std.endl


   # detailed output
   if mEnableDetOutput:      mWriteResult=True
      Info("StandardHypoTestInvDemo", "detailed output will be written in output result file")



   # write result in a file
   if r != NULL and mWriteResult:
      # write to a file the results
       char calcType = (calculatorType == 0) ? "Freq" : (calculatorType == 1) ? "Hybr" : "Asym"
       char limitType = (useCLs) ? "CLs" : "Cls+b"
       char scanType = (npoints < 0) ? "auto" : "grid"
      if mResultFileName.IsNull():         mResultFileName = ROOT.TString.Format("%s_%s_%s_ts%d_",calcType,limitType,scanType,testStatType)
         #strip the / from the filename
         if mMassValue.size()>0:            mResultFileName += mMassValue.c_str()
            mResultFileName += "_"


         name = fileNameBase
         name.Replace(0, name.Last('/')+1, "")
         mResultFileName += name


      # get (if existing) rebuilt UL distribution
      uldistFile = "RULDist.root"
      ROOT.TObject ulDist = 0
      existULDist = not gSystem.AccessPathName(uldistFile)
      if existULDist:         ROOT.TFile fileULDist = ROOT.TFile.Open(uldistFile)
         if fileULDist) ulDist= fileULDist.Get("RULDist":



      ROOT.TFile fileOut = ROOT.TFile(mResultFileName, "RECREATE")
      r.Write()
      if ulDist) ulDist.Write(:
      Info("StandardHypoTestInvDemo", "HypoTestInverterResult has been written in the file %s",mResultFileName.Data())

      fileOut.Close()



   # plot the result ( p values vs scan points)
   typeName = ""
   if calculatorType == 0 :
      typeName = "Frequentist"
   if calculatorType == 1 :
      typeName = "Hybrid"
   elif calculatorType == 2 or calculatorType == 3:      typeName = "Asymptotic"
      mPlotHypoTestResult = False


    char resultName = r.GetName()
   plotTitle = ROOT.TString.Format("%s CL Scan for workspace %s",typeName.c_str(),resultName)
   HypoTestInverterPlot *plot = HypoTestInverterPlot("HTI_Result_Plot",plotTitle,r)

   # plot in a canvas with style
   c1Name = ROOT.TString.Format("%s_Scan",typeName.c_str())
   ROOT.TCanvas c1 = ROOT.TCanvas(c1Name)
   c1.SetLogy(False)

   plot.Draw("CLb 2CL");  # plot all and Clb

   # if useCLs:
   #    plot.Draw("CLb 2CL");  # plot all and Clb
   # else:
   #    plot.Draw("");  # plot all and Clb

    nEntries = r.ArraySize()

   # plot test statistics distributions for the two hypothesis
   if mPlotHypoTestResult:      ROOT.TCanvas c2 = ROOT.TCanvas()
      if nEntries > 1:         ny = ROOT.TMath.CeilNint(TMath.Sqrt(nEntries))
         nx = ROOT.TMath.CeilNint(double(nEntries)/ny)
         c2.Divide( nx,ny)

      for (int i=0; i<nEntries; i++)         if nEntries > 1) c2.cd(i+1:
         SamplingDistPlot pl = plot.MakeTestStatPlot(i)
         pl.SetLogYaxis(True)
         pl.Draw()








# internal routine to run the inverter
HypoTestInverterResult *
RooStats.HypoTestInvTool.RunInverter(ROOT.RooWorkspace * w,
                                        char * modelSBName, * modelBName,
                                        char * dataName, type, testStatType,
                                       bool useCLs, npoints, poimin, poimax,
                                       int ntoys,
                                       bool useNumberCounting,
                                        char * nuisPriorName )
   std.cout << "Running HypoTestInverter on the workspace " << w.GetName() << std.endl

   w.Print()


   ROOT.RooAbsData data = w.data(dataName)
   if not data:      Error("StandardHypoTestDemo", "Not existing data %s",dataName)
      return 0

   else:
      std.cout << "Using data set " << dataName << std.endl

   if mUseVectorStore:      ROOT.RooAbsData.setDefaultStorageType(ROOT.RooAbsData.Vector)
      data.convertToVectorStore() 



   # get models from WS
   # get the modelConfig out of the file
   bModel = (ModelConfig*) w.obj(modelBName)
   sbModel = (ModelConfig*) w.obj(modelSBName)

   if not sbModel:      Error("StandardHypoTestDemo", "Not existing ROOT.RooStats.ModelConfig( %s",modelSBName)
      return 0

   # check the model
   if not sbModel.GetPdf():      Error("StandardHypoTestDemo", "Model %s has no pdf ",modelSBName)
      return 0

   if not sbModel.GetParametersOfInterest():      Error("StandardHypoTestDemo", "Model %s has no poi ",modelSBName)
      return 0

   if not sbModel.GetObservables():      Error("StandardHypoTestInvDemo", "Model %s has no observables ",modelSBName)
      return 0

   if not sbModel.GetSnapshot() :      Info("StandardHypoTestInvDemo", "Model %s has no snapshot  - make one using model poi",modelSBName)
      sbModel.SetSnapshot( *sbModel.GetParametersOfInterest() )


   # case of no systematics
   # remove nuisance parameters from model
   if noSystematics:       ROOT.RooArgSet nuisPar = sbModel.GetNuisanceParameters()
      if nuisPar and nuisPar.getSize() > 0:         std.cout << "StandardHypoTestInvDemo" << "  -  Switch off all systematics by setting them constant to their initial values" << std.endl
         ROOT.RooStats.SetAllConstant(*nuisPar)

      if bModel:          ROOT.RooArgSet bnuisPar = bModel.GetNuisanceParameters()
         if bnuisPar:
            ROOT.RooStats.SetAllConstant(*bnuisPar)



   if not bModel or bModel == sbModel:      Info("StandardHypoTestInvDemo", "The background model %s does not exist",modelBName)
      Info("StandardHypoTestInvDemo", "Copy it from ROOT.RooStats.ModelConfig( %s and set POI to zero",modelSBName)
      bModel = (ModelConfig*) sbModel.Clone()
      bModel.SetName(TString(modelSBName)+TString("_with_poi_0"))
      ROOT.RooRealVar var = dynamic_cast<RooRealVar*>(bModel.GetParametersOfInterest().first())
      if (not var) return 0
      oldval = var.getVal()
      var.setVal(0)
      bModel.SetSnapshot( ROOT.RooArgSet(*var)  )
      var.setVal(oldval)

   else:
      if not bModel.GetSnapshot() :         Info("StandardHypoTestInvDemo", "Model %s has no snapshot  - make one using model poi and 0 values ",modelBName)
         ROOT.RooRealVar var = dynamic_cast<RooRealVar*>(bModel.GetParametersOfInterest().first())
         if var:            oldval = var.getVal()
            var.setVal(0)
            bModel.SetSnapshot( ROOT.RooArgSet(*var)  )
            var.setVal(oldval)

         else:
            Error("StandardHypoTestInvDemo", "Model %s has no valid poi",modelBName)
            return 0




   # check model  has global observables when there are nuisance pdf
   # for the hybrid case the globobs are not needed
   if type != 1 :      hasNuisParam = (sbModel.GetNuisanceParameters() and sbModel.GetNuisanceParameters().getSize() > 0)
      hasGlobalObs = (sbModel.GetGlobalObservables() and sbModel.GetGlobalObservables().getSize() > 0)
      if hasNuisParam and not hasGlobalObs :         # try to see if model has nuisance parameters first
         ROOT.RooAbsPdf constrPdf = ROOT.RooStats.MakeNuisancePdf(*sbModel, "nuisanceConstraintPdf_sbmodel")
         if constrPdf:            Warning("StandardHypoTestInvDemo", "Model %s has nuisance parameters but no global observables associated",sbModel.GetName())
            Warning("StandardHypoTestInvDemo", "\tThe effect of the nuisance parameters will not be treated correctly ")




   # save all initial parameters of the model including the global observables
   ROOT.RooArgSet initialParameters
   ROOT.RooArgSet allParams = sbModel.GetPdf().getParameters(*data)
   allParams.snapshot(initialParameters)
   delete allParams

   # run first a data fit

    ROOT.RooArgSet poiSet = sbModel.GetParametersOfInterest()
   ROOT.RooRealVar *poi = (ROOT.RooRealVar*)poiSet.first()

   std.cout << "StandardHypoTestInvDemo : POI initial value:   " << poi.GetName() << " = " << poi.getVal()   << std.endl

   # fit the data first (need to use constraint )
   ROOT.TStopwatch tw

   doFit = initialFit
   if (testStatType == 0 and initialFit == -1) doFit = False;  # case of LEP test statistic
   if (type == 3  and initialFit == -1) doFit = False;         # case of Asymptoticcalculator with nominal Asimov
   poihat = 0

   if minimizerType.size()==0) minimizerType = ROOT.Math.MinimizerOptions.DefaultMinimizerType(:
   else:
      ROOT.Math.MinimizerOptions.SetDefaultMinimizer(minimizerType.c_str())

   Info("StandardHypoTestInvDemo", "Using %s as minimizer for computing the test statistic",
        ROOT.Math.MinimizerOptions.DefaultMinimizerType().c_str() )

   if doFit:
      # do the fit : By doing a fit the POI snapshot (for S+B)  is set to the fit value
      # and the nuisance parameters nominal values will be set to the fit value.
      # ROOT.This is relevant when using LEP test statistics

      Info( "StandardHypoTestInvDemo", " Doing a first fit to the observed data ")
      ROOT.RooArgSet constrainParams
      if sbModel.GetNuisanceParameters() ) constrainParams.add(*sbModel.GetNuisanceParameters():
      ROOT.RooStats.RemoveConstantParameters(&constrainParams)
      tw.Start()
      ROOT.RooFitResult fitres = sbModel.GetPdf().fitTo(*data,InitialHesse(False), Hesse(False),
                                                       Minimizer(minimizerType.c_str(), "Migrad"), Strategy(0), PrintLevel(mPrintLevel), Constrain(constrainParams), ROOT.RooFit.Save(True), Offset(ROOT.RooStats.IsNLLOffset()) )
      if fitres.status() != 0:         Warning("StandardHypoTestInvDemo", "Fit to the model failed - try with strategy 1 and perform first an Hesse computation")
         fitres = sbModel.GetPdf().fitTo(*data,InitialHesse(True), Hesse(False),Minimizer(minimizerType.c_str(), "Migrad"), Strategy(1), PrintLevel(mPrintLevel+1), Constrain(constrainParams),
                                           Save(True), Offset(ROOT.RooStats.IsNLLOffset()) )

      if fitres.status() != 0:
         Warning("StandardHypoTestInvDemo", " Fit still failed - continue anyway.....")


      poihat  = poi.getVal()
      std.cout << "StandardHypoTestInvDemo - Best Fit value : " << poi.GetName() << " = "
                << poihat << " +/- " << poi.getError() << std.endl
      std.cout << "Time for fitting : "; tw.Print()

      #save best fit value in the poi snapshot
      sbModel.SetSnapshot(*sbModel.GetParametersOfInterest())
      std.cout << "StandardHypoTestInvo: snapshot of S+B Model " << sbModel.GetName()
                << " is set to the best fit value" << std.endl



   # print a message in case of LEP test statistics because it affects result by doing or not doing a fit
   if testStatType == 0:      if not doFit:
         Info("StandardHypoTestInvDemo", "Using LEP test statistic - an initial fit is not done and the ROOT.TS will use the nuisances at the model value")
      else:
         Info("StandardHypoTestInvDemo", "Using LEP test statistic - an initial fit has been done and the ROOT.TS will use the nuisances at the best fit value")



   # build test statistics and hypotest calculators for running the inverter

   SimpleLikelihoodRatioTestStat slrts(*sbModel.GetPdf(),*bModel.GetPdf())

   # null parameters must includes snapshot of poi plus the nuisance values
   ROOT.RooArgSet nullParams(*sbModel.GetSnapshot())
   if sbModel.GetNuisanceParameters()) nullParams.add(*sbModel.GetNuisanceParameters():
   if sbModel.GetSnapshot()) slrts.SetNullParameters(nullParams:
   ROOT.RooArgSet altParams(*bModel.GetSnapshot())
   if bModel.GetNuisanceParameters()) altParams.add(*bModel.GetNuisanceParameters():
   if bModel.GetSnapshot()) slrts.SetAltParameters(altParams:
   if mEnableDetOutput) slrts.EnableDetailedOutput(:

   # ratio of profile likelihood - need to pass snapshot for the alt
   RatioOfProfiledLikelihoodsTestStat
      ropl(*sbModel.GetPdf(), *bModel.GetPdf(), bModel.GetSnapshot())
   ropl.SetSubtractMLE(False)
   if testStatType == 11) ropl.SetSubtractMLE(True:
   ropl.SetPrintLevel(mPrintLevel)
   ropl.SetMinimizer(minimizerType.c_str())
   if mEnableDetOutput) ropl.EnableDetailedOutput(:

   ProfileLikelihoodTestStat profll(*sbModel.GetPdf())
   if testStatType == 3) profll.SetOneSided(True:
   if testStatType == 4) profll.SetSigned(True:
   profll.SetMinimizer(minimizerType.c_str())
   profll.SetPrintLevel(mPrintLevel)
   if mEnableDetOutput) profll.EnableDetailedOutput(:

   profll.SetReuseNLL(mOptimize)
   slrts.SetReuseNLL(mOptimize)
   ropl.SetReuseNLL(mOptimize)

   if mOptimize:      profll.SetStrategy(0)
      ropl.SetStrategy(0)
      ROOT.Math.MinimizerOptions.SetDefaultStrategy(0)


   if (mMaxPoi > 0) poi.setMax(mMaxPoi);  # increase limit

   MaxLikelihoodEstimateTestStat maxll(*sbModel.GetPdf(),*poi)
   NumEventsTestStat nevtts

   AsymptoticCalculator.SetPrintLevel(mPrintLevel)

   # create the HypoTest calculator class
   HypoTestCalculatorGeneric hc = 0
   if type == 0) hc = FrequentistCalculator(*data, *bModel, *sbModel:
   elif type == 1) hc = HybridCalculator(*data, *bModel, *sbModel:
   # elif type == 2 ) hc = AsymptoticCalculator(*data, *bModel, *sbModel, False, mAsimovBins:
   # elif (type == 3 ) hc = AsymptoticCalculator(*data, *bModel, *sbModel, ROOT.True, mAsimovBins);  # for using Asimov data generated with nominal values
   elif type == 2 ) hc = AsymptoticCalculator(*data, *bModel, *sbModel, :
   elif (type == 3 ) hc = AsymptoticCalculator(*data, *bModel, *sbModel, ROOT.True );  # for using Asimov data generated with nominal values
   else:
      Error("StandardHypoTestInvDemo", "Invalid - type = %d supported values are only :\n\t\t\t 0 (Frequentist) , 1 (Hybrid) , 2 (Asymptotic) ",type)
      return 0


   # set the test statistic
   ROOT.TestStatistic testStat = 0
   if (testStatType == 0) testStat = &slrts
   if (testStatType == 1 or testStatType == 11) testStat = &ropl
   if (testStatType == 2 or testStatType == 3 or testStatType == 4) testStat = &profll
   if (testStatType == 5) testStat = &maxll
   if (testStatType == 6) testStat = &nevtts

   if testStat == 0:      Error("StandardHypoTestInvDemo", "Invalid - test type = %d supported values are only :\n\t\t\t 0 (SLR) , 1 (Tevatron) , 2 (PLR), 3 (PLR1), 4(MLE)",testStatType)
      return 0



   ROOT.ToyMCSampler *toymcs = (ToyMCSampler*)hc.GetTestStatSampler()
   if toymcs and (type == 0 or type == 1) :      # look if pdf is number counting or extended
      if sbModel.GetPdf().canBeExtended() :         if useNumberCounting)   Warning("StandardHypoTestInvDemo", "Pdf is extended: but number counting flag is set: ignore it ":

      else:
         # for not extended pdf
         if not useNumberCounting  :            nEvents = data.numEntries()
            Info("StandardHypoTestInvDemo", "Pdf is not extended: number of events to generate taken  from observed data set is %d",nEvents)
            toymcs.SetNEventsPerToy(nEvents)

         else:
            Info("StandardHypoTestInvDemo", "using a number counting pdf")
            toymcs.SetNEventsPerToy(1)



      toymcs.SetTestStatistic(testStat)

      if data.isWeighted() and not mGenerateBinned:         Info("StandardHypoTestInvDemo", "Data set is weighted, nentries = %d and sum weights = %8.1f but toy generation is unbinned - it would be faster to set mGenerateBinned to ROOT.True\n",data.numEntries(), data.sumEntries())

      toymcs.SetGenerateBinned(mGenerateBinned)

      toymcs.SetUseMultiGen(mOptimize)

      if mGenerateBinned and  sbModel.GetObservables().getSize() > 2:         Warning("StandardHypoTestInvDemo", "generate binned is activated but the number of ovservable is %d. ROOT.Too much memory could be needed for allocating all the bins",sbModel.GetObservables().getSize() )


      # set the random seed if needed
      if mRandomSeed >= 0) ROOT.RooRandom.randomGenerator().SetSeed(mRandomSeed:



   # specify if need to re-use same toys
   if reuseAltToys:      hc.UseSameAltToys()


   if type == 1:      HybridCalculator *hhc = dynamic_cast<HybridCalculator*> (hc)
      assert(hhc)

      hhc.SetToys(ntoys,ntoys/mNToysRatio); # can use less ntoys for b hypothesis

      # remove global observables from ROOT.RooStats.ModelConfig( (self is probably not needed anymore in 5.32)
      bModel.SetGlobalObservables(ROOT.RooArgSet() )
      sbModel.SetGlobalObservables(ROOT.RooArgSet() )


      # check for nuisance prior pdf in case of nuisance parameters
      if bModel.GetNuisanceParameters() or sbModel.GetNuisanceParameters() :
         # fix for using multigen (does not work in self case)
         toymcs.SetUseMultiGen(False)
         ROOT.ToyMCSampler.SetAlwaysUseMultiGen(False)

         ROOT.RooAbsPdf nuisPdf = 0
         if nuisPriorName) nuisPdf = w.pdf(nuisPriorName:
         # use prior defined first in bModel (then in SbModel)
         if not nuisPdf:            Info("StandardHypoTestInvDemo", "No nuisance pdf given for the HybridCalculator - try to deduce  pdf from the model")
            if bModel.GetPdf() and bModel.GetObservables() :
               nuisPdf = ROOT.RooStats.MakeNuisancePdf(*bModel, "nuisancePdf_bmodel")
            else:
               nuisPdf = ROOT.RooStats.MakeNuisancePdf(*sbModel, "nuisancePdf_sbmodel")

         if not nuisPdf :            if bModel.GetPriorPdf():               nuisPdf = bModel.GetPriorPdf()
               Info("StandardHypoTestInvDemo", "No nuisance pdf given - try to use %s that is defined as a prior pdf in the B model",nuisPdf.GetName())

            else:
               Error("StandardHypoTestInvDemo", "Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified or can be derived")
               return 0


         assert(nuisPdf)
         Info("StandardHypoTestInvDemo", "Using as nuisance Pdf ... " )
         nuisPdf.Print()

          ROOT.RooArgSet nuisParams = (bModel.GetNuisanceParameters() ) ? bModel.GetNuisanceParameters() : sbModel.GetNuisanceParameters()
         ROOT.RooArgSet np = nuisPdf.getObservables(*nuisParams)
         if np.getSize() == 0:            Warning("StandardHypoTestInvDemo", "Prior nuisance does not depend on nuisance parameters. ROOT.They will be smeared in their full range")

         delete np

         hhc.ForcePriorNuisanceAlt(*nuisPdf)
         hhc.ForcePriorNuisanceNull(*nuisPdf)




   elif type == 2 or type == 3:      if testStatType == 3) ((AsymptoticCalculator*) hc).SetOneSided(True:
      if testStatType != 2 and testStatType != 3:
         Warning("StandardHypoTestInvDemo", "Only the PL test statistic can be used with AsymptoticCalculator - use by default a two-sided PL")

   elif type == 0 or type == 1:      ((FrequentistCalculator*) hc).SetToys(ntoys,ntoys/mNToysRatio)
      # store also the fit information for each poi point used by calculator based on toys
      if mEnableDetOutput) ((FrequentistCalculator*) hc).StoreFitInfo(True:


   # Get the result
   ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)



   HypoTestInverter calc(*hc)
   calc.SetConfidenceLevel(confidenceLevel)


   calc.UseCLs(useCLs)
   calc.SetVerbose(True)

   # can speed up using proof-lite
   if mUseProof:      ProofConfig pc(*w, mNWorkers, "", kFALSE)
      toymcs.SetProofConfig(&pc);    # enable proof



   if npoints > 0:      if poimin > poimax:         # if no min/max given scan between MLE and +4 sigma
         poimin = int(poihat)
         poimax = int(poihat +  4 * poi.getError())

      std.cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std.endl
      calc.SetFixedScan(npoints,poimin,poimax)

   else:
      #poi.setMax(10*int( (poihat+ 10 *poi.getError() )/10 ) )
      std.cout << "Doing an  automatic scan  in interval : " << poi.getMin() << " , " << poi.getMax() << std.endl


   tw.Start()
   HypoTestInverterResult r = calc.GetInterval()
   std.cout << "Time to perform limit scan \n"
   tw.Print()

   if mRebuild:
      std.cout << "\n***************************************************************\n"
      std.cout << "Rebuild the upper limit distribution by re-generating set of pseudo-experiment and re-compute for each of them a upper limit\n\n"


      allParams = sbModel.GetPdf().getParameters(*data)

      # define on which value of nuisance parameters to do the rebuild
      # default is best fit value for bmodel snapshot



      if mRebuildParamValues != 0:         # set all parameters to their initial workspace values
         *allParams = initialParameters

      if mRebuildParamValues == 0 or mRebuildParamValues == 1 :          ROOT.RooArgSet constrainParams
          if sbModel.GetNuisanceParameters() ) constrainParams.add(*sbModel.GetNuisanceParameters():
          ROOT.RooStats.RemoveConstantParameters(&constrainParams)

           ROOT.RooArgSet poiModel = sbModel.GetParametersOfInterest()
          bModel.LoadSnapshot()

          # do a profile using the B model snapshot
          if mRebuildParamValues == 0 :
             ROOT.RooStats.SetAllConstant(*poiModel,True)

             sbModel.GetPdf().fitTo(*data,InitialHesse(False), Hesse(False),
                                      Minimizer(minimizerType.c_str(), "Migrad"), Strategy(0), PrintLevel(mPrintLevel), Constrain(constrainParams), Offset(ROOT.RooStats.IsNLLOffset()) )


             std.cout << "rebuild using fitted parameter value for B-model snapshot" << std.endl
             constrainParams.Print("v")

             ROOT.RooStats.SetAllConstant(*poiModel,False)


      std.cout << "StandardHypoTestInvDemo: Initial parameters used for rebuilding: "
      ROOT.RooStats.PrintListContent(*allParams, std.cout)
      delete allParams

      calc.SetCloseProof(1)
      tw.Start()
      SamplingDistribution limDist = calc.GetUpperLimitDistribution(True,mNToyToRebuild)
      std.cout << "Time to rebuild distributions " << std.endl
      tw.Print()

      if limDist:         std.cout << "Expected limits after rebuild distribution " << std.endl
         std.cout << "expected upper limit  (median of limit distribution) " << limDist.InverseCDF(0.5) << std.endl
         std.cout << "expected -1 sig limit (0.16% quantile of limit dist) " << limDist.InverseCDF(ROOT.Math.normal_cdf(-1)) << std.endl
         std.cout << "expected +1 sig limit (0.84% quantile of limit dist) " << limDist.InverseCDF(ROOT.Math.normal_cdf(1)) << std.endl
         std.cout << "expected -2 sig limit (.025% quantile of limit dist) " << limDist.InverseCDF(ROOT.Math.normal_cdf(-2)) << std.endl
         std.cout << "expected +2 sig limit (.975% quantile of limit dist) " << limDist.InverseCDF(ROOT.Math.normal_cdf(2)) << std.endl

         # Plot the upper limit distribution
         SamplingDistPlot limPlot( (mNToyToRebuild < 200) ? 50 : 100)
         limPlot.AddSamplingDistribution(limDist)
         limPlot.GetTH1F().SetStats(True); # display statistics
         limPlot.SetLineColor(ROOT.kBlue)
         ROOT.TCanvas("limPlot", "Upper Limit Distribution")
         limPlot.Draw()

         #/ save result in a file
         limDist.SetName("RULDist")
         ROOT.TFile fileOut = ROOT.TFile("RULDist.root", "RECREATE")
         limDist.Write()
         fileOut.Close()


         #update r to a updated result object containing the rebuilt expected p-values distributions
         # (it will not recompute the expected limit)
         if (r) delete r;  # need to delete previous object since GetInterval will return a cloned copy
         r = calc.GetInterval()


      else:
         std.cout << "ERROR : failed to re-build distributions " << std.endl


   return r




def ReadResult(self, * fileName, * resultName="", useCLs=True):   # read a previous stored result from a file given the result name

   StandardHypoTestInvDemo(fileName, resultName, "", "", "",0,0,useCLs)



#ifdef USE_AS_MAIN
def main(self):    StandardHypoTestInvDemo()

#endif




