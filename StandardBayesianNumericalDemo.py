# Standard demo of the numerical Bayesian calculator

'''


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

The BayesianCalculator is based on Bayes's theorem
and performs the integration using ROOT's numeric integration utilities
'''

#include "TFile.h"
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"

#include "RooUniform.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/RooStatsUtils.h"
#include "RooPlot.h"
#include "TSystem.h"

import ROOT
using namespace ROOT.RooStats


integrationType = "";  # integration ROOT.Type (default is adaptive (numerical integration)
                               # possible values are "TOYMC" (toy MC integration, when nuisances have a constraints pdf)
                               #  "VEGAS" , "MISER", or "PLAIN"  (these are all possible MC integration)
nToys = 10000;             # number of toys used for the MC integrations - for Vegas should be probably set to an higher value
scanPosterior = False;    # flag to compute interval by scanning posterior (it is more robust but maybe less precise)
nScanPoints = 20;          # number of points for scanning the posterior (scanPosterior = False it is used only for plotting). Use by default a low value to speed-up tutorial
intervalType = 1;          # type of interval (0 is shortest, central, upper limit)
maxPOI = -999;        # force a different value of POI for doing the scan (default is given value)
nSigmaNuisance = -1;   # force integration of nuisance parameters to be withing nSigma of their error (do first a model fit to find nuisance error)

void StandardBayesianNumericalDemo( infile = "",
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
  # create and use the BayesianCalculator
  # to find and plot the 95% credible interval
  # on the parameter of interest as specified
  # in the model config

  # before we do that, must specify our prior
  # it belongs in the model config, it may not have
  # been specified
  ROOT.RooUniform prior("prior", "",*mc.GetParametersOfInterest())
  getattr(w, 'import')(prior)
  mc.SetPriorPdf(*w.pdf("prior"))

  # do without systematics
  #mc.SetNuisanceParameters(ROOT.RooArgSet() )
  if nSigmaNuisance > 0:     ROOT.RooAbsPdf pdf = mc.GetPdf()
     assert(pdf); 
     ROOT.RooFitResult res = pdf.fitTo(*data, ROOT.RooFit.Save(True), Minimizer(ROOT.Math.MinimizerOptions.DefaultMinimizerType().c_str()), Hesse(True),
                                     PrintLevel(ROOT.Math.MinimizerOptions.DefaultPrintLevel()-1) )

     res.Print()
     ROOT.RooArgList nuisPar(*mc.GetNuisanceParameters())
     for (i = 0; i < nuisPar.getSize(); ++i)        ROOT.RooRealVar v = dynamic_cast<RooRealVar*> (&nuisPar[i] )
        assert( v)
        v.setMin( ROOT.TMath.Max( v.getMin(), v.getVal() - nSigmaNuisance * v.getError() ) ); 
        v.setMax( ROOT.TMath.Min( v.getMax(), v.getVal() + nSigmaNuisance * v.getError() ) )
        std.cout << "setting interval for nuisance  " << v.GetName() << " : [ " << v.getMin() << " , " << v.getMax() << " ]" << std.endl




  BayesianCalculator bayesianCalc(*data,*mc)
  bayesianCalc.SetConfidenceLevel(0.95); # 95% interval

  # default of the calculator is central interval.  here use shortest , or upper limit depending on input
  # doing a shortest interval might require a longer time since it requires a scan of the posterior function
  if (intervalType == 0)  bayesianCalc.SetShortestInterval(); # for shortest interval
  if (intervalType == 1)  bayesianCalc.SetLeftSideTailFraction(0.5); # for central interval
  if (intervalType == 2)  bayesianCalc.SetLeftSideTailFraction(0.); # for upper limit

  if not integrationType.IsNull() :     bayesianCalc.SetIntegrationType(integrationType); # set integrationType
     bayesianCalc.SetNumIters(nToys); # set number of ietrations (i.e. number of toys for MC integrations)


  # in case of toyMC make a nnuisance pdf
  if integrationType.Contains("TOYMC") :    ROOT.RooAbsPdf nuisPdf = ROOT.RooStats.MakeNuisancePdf(*mc, "nuisance_pdf")
    cout << "using ROOT.TOYMC integration: make nuisance pdf from the model " << std.endl
    nuisPdf.Print()
    bayesianCalc.ForceNuisancePdf(*nuisPdf)
    scanPosterior = ROOT.True; # for ROOT.ToyMC the posterior is scanned anyway so used given points


  # compute interval by scanning the posterior function
  if scanPosterior:
     bayesianCalc.SetScanOfPosterior(nScanPoints)

  poi = (ROOT.RooRealVar*) mc.GetParametersOfInterest().first()
  if maxPOI != -999 and  maxPOI > poi.getMin():
    poi.setMax(maxPOI)


  interval = bayesianCalc.GetInterval()

  # print out the iterval on the first Parameter of Interest
  cout << "\n95% interval on " << poi.GetName()<<" is : ["<<
    interval.LowerLimit() << ", "<<
    interval.UpperLimit() <<"] "<<endl


  # make a plot
  # since plotting may take a long time (it requires evaluating
  # the posterior in many points) self command will speed up
  # by reducing the number of points to plot - do 50

  cout << "\nDrawing plot of posterior function....." << endl

  bayesianCalc.SetScanOfPosterior(nScanPoints)

  ROOT.RooPlot plot = bayesianCalc.GetPosteriorPlot()
  plot.Draw()


