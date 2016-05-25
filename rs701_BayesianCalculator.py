# /
#
# 'Bayesian Calculator' ROOT.RooStats tutorial macro #701
# author: Gregory Schott
# date Sep 2009
#
# ROOT.This tutorial shows an example of using the BayesianCalculator class
#
# /


import ROOT


def rs701_BayesianCalculator(useBkg=True, confLevel=0.90):

    w = ROOT.RooWorkspace("w", True)
    w.factory("SUM::pdf(s[0.001,15]*Uniform(x[0,1]),b[1,0,2]*Uniform(x))")
    w.factory("Gaussian::prior_b(b,1,1)")
    w.factory("PROD::model(pdf,prior_b)")
    model = w.pdf("model")  # pdf*priorNuisance
    nuisanceParameters = ROOT.RooArgSet(w.var("b"))

    POI = w.var("s")
    priorPOI = w.factory("Uniform::priorPOI(s)")
    priorPOI2 = w.factory("GenericPdf::priorPOI2('1/sqrt(@0)',s)")

    w.factory("n[3]")  # observed number of events
    # create a data set with n observed events
    data = ROOT.RooDataSet("data", "", ROOT.RooArgSet(
        (w.var("x")), (w.var("n"))), "n")
    data.add(ROOT.RooArgSet((w.var("x"))), w.var("n").getVal())

    # to suppress messgaes when pdf goes to zero
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    nuisPar = ROOT.RooArgSet()
    if (useBkg):
        nuisPar = nuisanceParameters
    # if not useBkg) ((ROOT.RooRealVar *)w.var("b")).setVal(0:

    size = 1. - confLevel
    print "\nBayesian Result using a Flat prior "
    bcalc = ROOT.RooStats.BayesianCalculator(data, model, ROOT.RooArgSet(POI), priorPOI, nuisPar)
    bcalc.SetTestSize(size)
    interval = bcalc.GetInterval()
    cl = bcalc.ConfidenceLevel()
    print cl, "% CL central interval: [ ", interval.LowerLimit(), " - ", interval.UpperLimit(), " ] or ", cl + (1. - cl) / 2, "% CL limits\n"
    plot = ROOT.RooPlot = bcalc.GetPosteriorPlot()
    c1 = ROOT.TCanvas("c1", "Bayesian Calculator Result")
    c1.Divide(1, 2)
    c1.cd(1)
    plot.Draw()
    c1.Update()

    print "\nBayesian Result using a 1/sqrt(s) prior  "
    bcalc2 = ROOT.RooStats.BayesianCalculator(data, model, ROOT.RooArgSet(POI), priorPOI2, nuisPar)
    bcalc2.SetTestSize(size)
    interval2 = bcalc2.GetInterval()
    cl = bcalc2.ConfidenceLevel()
    print cl, "% CL central interval: [ ", interval2.LowerLimit(), " - ", interval2.UpperLimit(), " ] or ", cl + (1. - cl) / 2, "% CL limits\n"

    plot2 = bcalc2.GetPosteriorPlot()
    c1.cd(2)
    plot2.Draw()
    ROOT.gPad.SetLogy()
    c1.Update()

    c1.SaveAs("rs701_BayesianCalculator.png")

    # observe one event while expecting one background event . the 95% CL upper limit on s is 4.10
    # observe one event while expecting zero background event . the 95% CL
    # upper limit on s is 4.74


if __name__ == "__main__":
    rs701_BayesianCalculator()
