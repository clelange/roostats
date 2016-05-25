# /
#
# 'Hypothesis ROOT.Test Inversion' ROOT.RooStats tutorial macro #801
# author: Gregory Schott
# date Sep 2009
#
# ROOT.This tutorial shows an example of using the HypoTestInverterOriginal class
#
# /


import ROOT


def rs801_HypoTestInverterOriginal():
    # prepare the model
    lumi = ROOT.RooRealVar("lumi", "luminosity", 1)
    r = ROOT.RooRealVar("r", "cross-section ratio", 3.74, 0, 50)
    ns = ROOT.RooFormulaVar("ns", "1*r*lumi", ROOT.RooArgList(lumi, r))
    nb = ROOT.RooRealVar("nb", "background yield", 1)
    x = ROOT.RooRealVar("x", "dummy observable", 0, 1)
    p0 = ROOT.RooConstVar(ROOT.RooFit.RooConst(0))
    flatPdf = ROOT.RooPolynomial("flatPdf", "flat PDF", x, ROOT.RooArgList(p0))
    totPdf = ROOT.RooAddPdf(
        "totPdf", "S+B model", ROOT.RooArgList(flatPdf, flatPdf), ROOT.RooArgList(ns, nb))
    bkgPdf = ROOT.RooExtendPdf("bkgPdf", "B-only model", flatPdf, nb)
    data = totPdf.generate(ROOT.RooArgSet(x), 1)

    # prepare the calculator
    myhc = ROOT.RooStats.HybridCalculatorOriginal(data, totPdf, bkgPdf, 0, 0)
    myhc.SetTestStatistic(2)
    myhc.SetNumberOfToys(1000)
    myhc.UseNuisance(False)

    # run the hypothesis-test invertion
    myInverter = ROOT.RooStats.HypoTestInverterOriginal(myhc, r)
    myInverter.SetTestSize(0.10)
    myInverter.UseCLs(True)
    # myInverter.RunFixedScan(5,1,6)
    # scan for a 95% UL
    myInverter.RunAutoScan(3., 5, myInverter.Size() / 2, 0.005)
    # run an alternative autoscan algorithm
    # myInverter.RunAutoScan(1,6,myInverter.Size()/2,0.005,1)
    # myInverter.RunOnePoint(3.9)

    results = myInverter.GetInterval()

    myInverterPlot = ROOT.RooStats.HypoTestInverterPlot(
        "myInverterPlot", "", results)
    c = ROOT.TCanvas("rs801_HypoTestInverterOriginal",
                     "rs801_HypoTestInverterOriginal", 800, 600)
    gr1 = myInverterPlot.MakePlot()
    gr1.Draw("ALP")

    ulError = results.UpperLimitEstimatedError()

    upperLimit = results.UpperLimit()
    print "The computed upper limit is: ", upperLimit
    print "an estimated error on self upper limit is: ", ulError
    # expected result: 4.10

    c.SaveAs("rs801_HypoTestInverterOriginal.png")


if __name__ == "__main__":
    rs801_HypoTestInverterOriginal()
