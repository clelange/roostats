# /
#
# 'Bernstein Correction' ROOT.RooStats tutorial macro
# author: Kyle Cranmer
# date March. 2009
#
# ROOT.This tutorial shows usage of a the BernsteinCorrection utility in ROOT.RooStats.
# ROOT.The idea is that one has a distribution coming either from data or Monte Carlo
# (called "reality" in the macro) and a nominal model that is not sufficiently
# flexible to take into account the real distribution.  One wants to take into
# account the systematic associated with self imperfect modeling by augmenting
# the nominal model with some correction term (in self case a polynomial).
# ROOT.The BernsteinCorrection utility will import into your workspace a corrected model
# given by nominal(x) * poly_N(x), poly_N is an n-th order polynomial in
# the Bernstein basis.  ROOT.The degree N of the polynomial is chosen by specifying the tolerance
# one has in adding an extra term to the polynomial.
# ROOT.The Bernstein basis is nice because it only has positive-definite terms
# and works well with PDFs.
# Finally, macro makes a plot of:
#  - the data (drawn from 'reality'),
#  - the best fit of the nominal model (blue)
#  - and the best fit corrected model.
# /


import ROOT


# ____________________________________
def rs_bernsteinCorrection():
    # set range of observable
    lowRange = -1
    highRange = 5

    # make a ROOT.RooRealVar for the observable
    x = ROOT.RooRealVar("x", "x", lowRange, highRange)

    # ROOT.True model
    narrow = ROOT.RooGaussian(
        "narrow", "", x, ROOT.RooFit.RooConst(0.), ROOT.RooFit.RooConst(.8))
    wide = ROOT.RooGaussian(
        "wide", "", x, ROOT.RooFit.RooConst(0.), ROOT.RooFit.RooConst(2.))
    reality = ROOT.RooAddPdf("reality", "", ROOT.RooArgList(
        narrow, wide), ROOT.RooArgList(ROOT.RooFit.RooConst(0.8)))

    data = reality.generate(ROOT.RooArgSet(x), 1000)

    # nominal model
    sigma = ROOT.RooRealVar("sigma", "", 1., 0, 10)
    nominal = ROOT.RooGaussian(
        "nominal", "", x, ROOT.RooFit.RooConst(0.), sigma)

    wks = ROOT.RooWorkspace("myWorksspace")

    getattr(wks, 'import')(data, ROOT.RooFit.Rename("data"))
    getattr(wks, 'import')(nominal)

    # use Minuit2
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")

    # ROOT.The tolerance sets the probability to add an unnecessary term.
    # lower tolerance will add fewer terms, higher tolerance
    # will add more terms and provide a more flexible function.
    tolerance = 0.05
    bernsteinCorrection = ROOT.RooStats.BernsteinCorrection(tolerance)
    degree = bernsteinCorrection.ImportCorrectedPdf(
        wks, "nominal", "x", "data")

    if degree < 0:
        ROOT.RooStats.Error("rs_bernsteinCorrection",
                            "Bernstein correction failed not  ")
        return

    print " Correction based on Bernstein Poly of degree ", degree

    frame = x.frame()
    data.plotOn(frame)
    # plot the best fit nominal model in blue
    minimType = ROOT.Math.MinimizerOptions.DefaultMinimizerType()
    nominal.fitTo(data, ROOT.RooFit.PrintLevel(0),
                  ROOT.RooFit.Minimizer(minimType))
    nominal.plotOn(frame)

    # plot the best fit corrected model in red
    corrected = wks.pdf("corrected")
    if (not corrected):
        return

    # fit corrected model
    corrected.fitTo(data, ROOT.RooFit.PrintLevel(0),
                    ROOT.RooFit.Minimizer(minimType))
    corrected.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))

    # plot the correction term (* norm constant) in dashed green
    # should make norm constant just be 1, depend on binning of data
    poly = wks.pdf("poly")
    if (poly):
        poly.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGreen),
                    ROOT.RooFit.LineStyle(ROOT.kDashed))

    # self is a switch to check the sampling distribution
    # of -2 log LR for two comparisons:
    # the first is for n-1 vs. n degree polynomial corrections
    # the second is for n vs. n+1 degree polynomial corrections
    # Here we choose n to be the one chosen by the tolerance
    # critereon above, eg. n = "degree" in the code.
    # Setting self to True is takes about 10 min.
    checkSamplingDist = True
    numToyMC = 20  # increse self value for sensible results

    c1 = ROOT.TCanvas()
    if checkSamplingDist:
        c1.Divide(1, 2)
        c1.cd(1)

    frame.Draw()
    ROOT.gPad.Update()

    if checkSamplingDist:    # check sampling dist
        ROOT.Math.MinimizerOptions.SetDefaultPrintLevel(-1)
        samplingDist = ROOT.TH1F("samplingDist", "", 20, 0, 10)
        samplingDistExtra = ROOT.TH1F("samplingDistExtra", "", 20, 0, 10)
        bernsteinCorrection.CreateQSamplingDist(
            wks, "nominal", "x", "data", samplingDist, samplingDistExtra, degree, numToyMC)

        c1.cd(2)
        samplingDistExtra.SetLineColor(ROOT.kRed)
        samplingDistExtra.Draw()
        samplingDist.Draw("same")
    c1.SaveAs("rs_bernsteinCorrection.png")


if __name__ == "__main__":
    rs_bernsteinCorrection()
