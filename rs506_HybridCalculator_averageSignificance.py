#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/HybridResult.h"
#include "RooStats/HybridPlot.h"

#include "TFile.h"
#include "TStopwatch.h"
#include "TCanvas.h"

using namespace ROOT.RooFit
import ROOT

def rs506_HybridCalculator_averageSignificance(self, fname="WS_GaussOverFlat_withSystematics.root", ntoys=2500, outputplot="hc_averagesign_counting_syst.ps"):
  ROOT.RooRandom.randomGenerator().SetSeed(100)
  ROOT.TStopwatch t
  t.Start()

  ROOT.TFile* file =new ROOT.TFile(fname)
  my_WS = (ROOT.RooWorkspace*) file.Get("myWS")
  if (not my_WS) return

  #Import the objects needed
  ROOT.RooAbsPdf* model=my_WS.pdf("model")
  ROOT.RooAbsPdf* priorNuisance=my_WS.pdf("priorNuisance")

  # ROOT.RooArgSet* paramInterestSet=my_WS.set("paramInterest")
  #RooRealVar* paramInterest=(ROOT.RooRealVar*) paramInterestSet.first()
  ROOT.RooAbsPdf* modelBkg=my_WS.pdf("modelBkg")
   ROOT.RooArgSet* nuisanceParam=my_WS.set("parameters")
  ROOT.RooArgList observable(*(my_WS.set("observables") ) )

  HybridCalculator * hc=new HybridCalculator(*model,*modelBkg,observable)
  hc.SetNumberOfToys(ntoys)
  hc.SetTestStatistic(1)
  bool usepriors=False
  if priorNuisance!=0:    hc.UseNuisance(kTRUE)
    hc.SetNuisancePdf(*priorNuisance)
    usepriors=True
    nuisanceParam.Print()
    hc.SetNuisanceParameters(*nuisanceParam)
  }else:
    hc.UseNuisance(kFALSE)


  ROOT.RooRandom.randomGenerator().SetSeed(0)
  HybridResult* hcresult=hc.Calculate(ntoys,usepriors)
  hcplot = hcresult.GetPlot("HybridPlot", "Toy MC Q ",100)
  mean_sb_toys_test_stat = hcplot.GetSBmean()
  hcresult.SetDataTestStatistics(mean_sb_toys_test_stat)

  mean_significance = hcresult.Significance()

  cout<<"significance:" <<mean_significance<<endl
  hcplot = hcresult.GetPlot("HybridPlot", "Toy MC -2ln Q ",100)

  ROOT.TCanvas*c1=new ROOT.TCanvas()
  c1.cd()
  hcplot.Draw(outputplot)
  c1.Print(outputplot)
  file.Close()
  t.Stop()
  t.Print()

