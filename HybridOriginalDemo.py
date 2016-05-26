# Example on how to use the HybridCalculatorOriginal class
#
# Author: Gregory Schott
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooConstVar.h"

#ifndef __CINT__  # problem including self file with CINT
#include "RooGlobalFunc.h"
#endif

#include "RooStats/HybridCalculatorOriginal.h"
#include "RooStats/HybridResult.h"
#include "RooStats/HybridPlot.h"

def HybridOriginalDemo(self, ntoys = 1000):
  #***********************************************************************#
  # ROOT.This macro show an example on how to use ROOT.RooStats/HybridCalculatorOriginal #
  #***********************************************************************#
  #
  # With self example, should CL_sb = 0.130 CL_b = 0.946
  # (if data had -2lnQ = -3.0742). You can compare to the expected plot:
  # http:#www-ekp.physik.uni-karlsruhe.de/~schott/roostats/hybridplot_example.png

  import ROOT
  using namespace ROOT.RooStats

  #/ set ROOT.RooFit random seed
  ROOT.RooRandom.randomGenerator().SetSeed(3007)

  #/ build the models for background and signal+background
  x = ROOT.RooRealVar("x", "",-3,3)
  ROOT.RooArgList observables(x); # variables to be generated

  # gaussian signal
#  ROOT.RooRealVar sig_mean("sig_mean", "",0)
#  ROOT.RooRealVar sig_sigma("sig_sigma", "",0.8)
#  sig = ROOT.RooGaussian_pdf("sig_pdf", "",x,sig_mean,sig_sigma)
  sig = ROOT.RooGaussian_pdf("sig_pdf", "",x, ROOT.RooFit.RooConst(0.0), ROOT.RooFit.RooConst(0.8))
  ROOT.RooRealVar sig_yield("sig_yield", "",20,0,300)

  # flat background (extended PDF)
#  ROOT.RooRealVar bkg_slope("bkg_slope", "",0)
#  ROOT.RooPolynomial bkg_pdf("bkg_pdf", "",x,bkg_slope)
  ROOT.RooPolynomial bkg_pdf("bkg_pdf", "", x, ROOT.RooFit.RooConst(0))
  ROOT.RooRealVar bkg_yield("bkg_yield", "",40,0,300)
  ROOT.RooExtendPdf bkg_ext_pdf("bkg_ext_pdf", "",bkg_pdf,bkg_yield)

#  bkg_yield.setConstant(kTRUE)
  sig_yield.setConstant(kTRUE)

  # total sig+bkg (extended PDF)
  ROOT.RooAddPdf tot_pdf("tot_pdf", "", ROOT.RooArgList(sig_pdf,bkg_pdf), ROOT.RooArgList(sig_yield,bkg_yield))

  #/ build the prior PDF on the parameters to be integrated
  # gaussian contraint on the background yield ( N_B = 40 +/- 10  ie. 25% )
  ROOT.RooGaussian bkg_yield_prior("bkg_yield_prior", "",bkg_yield, ROOT.RooFit.RooConst(bkg_yield.getVal()), ROOT.RooFit.RooConst(10.))

  ROOT.RooArgSet nuisance_parameters(bkg_yield); # variables to be integrated

  #/ generate a data sample
  data = tot_pdf.generate(observables, ROOT.RooFit.Extended())

  #***********************************************************************#

  #/ run HybridCalculator on those inputs

  # use interface from HypoTest calculator by default

  HybridCalculatorOriginal myHybridCalc(*data, tot_pdf , bkg_ext_pdf ,
                                   &nuisance_parameters, &bkg_yield_prior)

  # here I use the default test statistics: 2*lnQ (optional)
  myHybridCalc.SetTestStatistic(1)
  #myHybridCalc.SetTestStatistic(3); # profile likelihood ratio

  myHybridCalc.SetNumberOfToys(ntoys)
  myHybridCalc.UseNuisance(True)

  # for speed up generation (do binned data)
  myHybridCalc.SetGenerateBinned(False)

  # calculate by running ntoys for the S+B and B hypothesis and retrieve the result
  myHybridResult = myHybridCalc.GetHypoTest()

  if not  myHybridResult:     std.cerr << "\nError returned from Hypothesis test" << std.endl
     return


  #/ run 1000 toys without gaussian prior on the background yield
  #myHybridResult = myHybridCalc.Calculate(*data,1000,False)

  #/ save the toy-MC results to file, way splitting into sub-batch jobs is possible
  #TFile fileOut("some_hybridresult.root", "RECREATE")
  #fileOut.cd()
  #myHybridResult.Write()
  #fileOut.Close()

  #/ read the results from a file
  #TFile fileIn("some_hybridresult.root")
  #myOtherHybridResult = (HybridResult*) fileIn.Get("myHybridCalc")

  #/ example on how to merge with toy-MC results obtained in another job
  #mergedHybridResult = HybridResult("mergedHybridResult", "self object holds merged results")
  #mergedHybridResult.Add(myHybridResult)
  #mergedHybridResult.Add(myOtherHybridResult)
  #/ or
  #myHybridResult.Add(myOtherHybridResult)

  #/ nice plot of the results
  myHybridPlot = myHybridResult.GetPlot("myHybridPlot", "Plot of results with HybridCalculatorOriginal",100)
  myHybridPlot.Draw()

  #/ recover and display the results
  clsb_data = myHybridResult.CLsplusb()
  clb_data = myHybridResult.CLb()
  cls_data = myHybridResult.CLs()
  data_significance = myHybridResult.Significance()
  min2lnQ_data = myHybridResult.GetTestStat_data()

  #/ compute the mean expected significance from toys
  mean_sb_toys_test_stat = myHybridPlot.GetSBmean()
  myHybridResult.SetDataTestStatistics(mean_sb_toys_test_stat)
  toys_significance = myHybridResult.Significance()

  std.cout << "Completed HybridCalculatorOriginal example:\n"
  std.cout << " - -2lnQ = " << min2lnQ_data << endl
  std.cout << " - CL_sb = " << clsb_data << std.endl
  std.cout << " - CL_b  = " << clb_data << std.endl
  std.cout << " - CL_s  = " << cls_data << std.endl
  std.cout << " - significance data = " << data_significance << std.endl
  std.cout << " - mean significance toys = " << toys_significance << std.endl



