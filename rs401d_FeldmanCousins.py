# /
#
# 'Neutrino Oscillation Example from Feldman & Cousins'
# author: Kyle Cranmer
# date March 2009
#
# ROOT.This tutorial shows a more complex example using the ROOT.RooStats.FeldmanCousins( utility
# to create a confidence interval for a toy neutrino oscillation experiment.
# ROOT.The example attempts to faithfully reproduce the toy example described in Feldman & Cousins'
# original paper, Phys.Rev.D57:3873-3889,1998.
#
# to run it:
# .x tutorials/roostats/rs401d_FeldmanCousins.C+
# /


import ROOT


def rs401d_FeldmanCousins(doFeldmanCousins=False, doMCMC=True):

    # to time the macro
    t = ROOT.TStopwatch()
    t.Start()

    '''
    ROOT.Taken from Feldman & Cousins paper, Phys.Rev.D57:3873-3889,1998.
    e-Print: physics/9711021 (see page 13.)

    Quantum mechanics dictates that the probability of such a transformation is given by the formula
    P (\nu \mu \to \nu e ) = sin^2 (2\theta) sin^2 (1.27 \Delta m^2 L /E )
    where P is the probability for a \nu \mu to transform into a \nu e , is the distance in km between
    the creation of the neutrino from meson decay and its interaction in the detector, is the
    neutrino energy in GeV, and \Delta m^2 = |m^2- m^2 | in (eV/c^2 )^2 .

    ROOT.To demonstrate how self works in practice, how it compares to alternative approaches
    that have been used, consider a toy model of a typical neutrino oscillation experiment.
    ROOT.The toy model is defined by the following parameters: Mesons are assumed to decay to
    neutrinos uniformly in a region 600 m to 1000 m from the detector. ROOT.The expected background
    from conventional \nu e interactions and misidentified \nu \mu interactions is assumed to be 100
    events in each of 5 energy bins which span the region from 10 to 60 GeV. We assume that
    the \nu \mu flux is such that if P (\nu \mu \to \nu e ) = 0.01 averaged over any bin, that bin would
    have an expected additional contribution of 100 events due to \nu \mu \to \nu e oscillations.
   '''

    # Make signal model model
    E = ROOT.RooRealVar("E", "", 15, 10, 60, "GeV")
    # need these units in formula
    L = ROOT.RooRealVar("L", "", .800, .600, 1.0, "km")
    deltaMSq = ROOT.RooRealVar(
        "deltaMSq", "#Delta m^{2}", 40, 1, 300, "eV/c^{2}")
    sinSq2theta = ROOT.RooRealVar(
        "sinSq2theta", "sin^{2}(2#theta)", .006, .0, .02)
    # RooRealVar deltaMSq("deltaMSq", "#Delta m^{2}",40,20,70, "eV/c^{2}")
    #  ROOT.RooRealVar sinSq2theta("sinSq2theta", "sin^{2}(2#theta)", .006,.001,.01)
    # PDF for oscillation only describes deltaMSq dependence, goes into sigNorm
    # 1) ROOT.The code for self PDF was created by issuing these commands
    #    root [0] ROOT.RooClassFactory x
    # root [1] x.makePdf("NuMuToNuE_Oscillation", "L,E,deltaMSq", "",
    # "pow(sin(1.27*deltaMSq*L/E),2)")
    x = ROOT.RooClassFactory()
    x.makePdf("NuMuToNuE_Oscillation", "L,E,deltaMSq", "", "pow(sin(1.27*deltaMSq*L/E),2)")
    # This is the way to handle user defined pdf (generated with RooClassFactory).
    # Compile once and for all in ROOT and then add it as a library
    ROOT.gROOT.ProcessLineSync(".x NuMuToNuE_Oscillation.cxx+")
    # ROOT.gSystem.Load("NuMuToNuE_Oscillation_cxx.so")

    PnmuTone = ROOT.NuMuToNuE_Oscillation(
        "PnmuTone", "P(#nu_{#mu} #rightarrow #nu_{e}", L, E, deltaMSq)

    # only E is observable, create the signal model by integrating out L
    sigModel = PnmuTone.createProjection(L)

    # create   \int dE' dL' P(E',L' | \Delta m^2).
    # Given ROOT.RooFit will renormalize the PDF in the range of the observables,
    # the average probability to oscillate in the experiment's acceptance
    # needs to be incorporated into the extended term in the likelihood.
    # Do self by creating a ROOT.RooAbsReal representing the integral and divide by
    # the area in the E-L plane.
    # ROOT.The integral should be over "primed" observables, we need
    # an independent copy of PnmuTone not to interfere with the original.

    # Independent copy for Integral
    EPrime = ROOT.RooRealVar("EPrime", "", 15, 10, 60, "GeV")
    # need these units in formula
    LPrime = ROOT.RooRealVar("LPrime", "", .800, .600, 1.0, "km")
    PnmuTonePrime = ROOT.NuMuToNuE_Oscillation("PnmuTonePrime", "P(#nu_{#mu} #rightarrow #nu_{e}",
                                                        LPrime, EPrime, deltaMSq)
    intProbToOscInExp = PnmuTonePrime.createIntegral(
        ROOT.RooArgSet(EPrime, LPrime))

    # Getting the flux is a bit tricky.  It is more celear to include a cross section term that is not
    # explicitly refered to in the text, eg.
    # # events bin = flux * cross-section for nu_e interaction in E bin * average prob nu_mu osc. to nu_e in bin
    # maxEventsInBin = flux * cross-section for nu_e interaction in E bin
    # maxEventsInBin * 1% chance bin =  100 events / bin
    # maxEventsInBin = 10,000.
    # for 5 bins, maxEventsTot = 50,000
    maxEventsTot = ROOT.RooConstVar(
        "maxEventsTot", "maximum number of sinal events", 50000)
    inverseArea = ROOT.RooConstVar("inverseArea", "1/(#Delta E #Delta L)",
                                   1. / (EPrime.getMax() - EPrime.getMin()) / (LPrime.getMax() - LPrime.getMin()))

    # sigNorm = maxEventsTot * (\int dE dL prob to oscillate in experiment /
    # Area) * sin^2(2\theta)
    sigNorm = ROOT.RooProduct("sigNorm", "", ROOT.RooArgSet(
        maxEventsTot, intProbToOscInExp, inverseArea, sinSq2theta))
    # bkg = 5 bins * 100 events / bin
    bkgNorm = ROOT.RooConstVar("bkgNorm", "normalization for background", 500)

    # flat background (0th order polynomial, no arguments for coefficients)
    bkgEShape = ROOT.RooPolynomial("bkgEShape", "flat bkg shape", E)

    # total model
    model = ROOT.RooAddPdf("model", "", ROOT.RooArgList(sigModel, bkgEShape),
                           ROOT.RooArgList(sigNorm, bkgNorm))

    # for debugging, model tree
    #  model.printCompactTree()
    #  model.graphVizTree("model.dot")

    # turn off some messages
    ROOT.RooMsgService.instance().setStreamStatus(0, ROOT.kFALSE)
    ROOT.RooMsgService.instance().setStreamStatus(1, ROOT.kFALSE)
    ROOT.RooMsgService.instance().setStreamStatus(2, ROOT.kFALSE)

    #######################
    # n events in data to data, sum of sig+bkg
    nEventsData = bkgNorm.getVal() + sigNorm.getVal()
    print "generate toy data nEvents = ", nEventsData
    # adjust random seed to get a toy dataset similar to one in paper.
    # Found by trial and error (3 trials, not very "fine tuned")
    ROOT.RooRandom.randomGenerator().SetSeed(3)
    # create a toy dataset
    data = model.generate(ROOT.RooArgSet(E), nEventsData)

    # /
    # make some plots
    dataCanvas = ROOT.TCanvas("dataCanvas")
    dataCanvas.Divide(2, 2)

    # plot the PDF
    dataCanvas.cd(1)
    hh = PnmuTone.createHistogram("hh", E, ROOT.RooFit.Binning(
        40), ROOT.RooFit.YVar(L, ROOT.RooFit.Binning(40)), ROOT.RooFit.Scaling(ROOT.kFALSE))
    hh.SetLineColor(ROOT.kBlue)
    hh.SetTitle("True Signal Model")
    hh.Draw("surf")

    # plot the data with the best fit
    dataCanvas.cd(2)
    Eframe = E.frame()
    data.plotOn(Eframe)
    model.fitTo(data, ROOT.RooFit.Extended())
    model.plotOn(Eframe)
    ras_sigModel = ROOT.RooArgSet(sigModel)
    ras_bkgEShape = ROOT.RooArgSet(bkgEShape)
    model.plotOn(Eframe, ROOT.RooFit.Components(
        ras_sigModel), ROOT.RooFit.LineColor(ROOT.kRed))
    model.plotOn(Eframe, ROOT.RooFit.Components(
        ras_bkgEShape), ROOT.RooFit.LineColor(ROOT.kGreen))
    model.plotOn(Eframe)
    Eframe.SetTitle("toy data with best fit model (and sig+bkg components)")
    Eframe.Draw()

    # plot the likelihood function
    dataCanvas.cd(3)
    nll = ROOT.RooNLLVar("nll", "nll", model, data, ROOT.RooFit.Extended())
    pll = ROOT.RooProfileLL(
        "pll", "", nll, ROOT.RooArgSet(deltaMSq, sinSq2theta))
    # hhh = nll.createHistogram("hhh",sinSq2theta, ROOT.RooFit.Binning(40),
    # ROOT.RooFit.YVar(deltaMSq, ROOT.RooFit.Binning(40)))
    hhh = pll.createHistogram("hhh", sinSq2theta, ROOT.RooFit.Binning(
        40), ROOT.RooFit.YVar(deltaMSq, ROOT.RooFit.Binning(40)), ROOT.RooFit.Scaling(ROOT.kFALSE))
    hhh.SetLineColor(ROOT.kBlue)
    hhh.SetTitle("Likelihood Function")
    hhh.Draw("surf")

    dataCanvas.Update()

    #############################
    # show use of Feldman-Cousins utility in ROOT.RooStats
    # set the distribution creator, encodes the test statistic
    parameters = ROOT.RooArgSet(deltaMSq, sinSq2theta)
    w = ROOT.RooWorkspace()

    modelConfig = ROOT.RooStats.ModelConfig()
    modelConfig.SetWorkspace(w)
    modelConfig.SetPdf(model)
    modelConfig.SetParametersOfInterest(parameters)

    fc = ROOT.RooStats.FeldmanCousins(data, modelConfig)
    fc.SetTestSize(.1)  # set size of test
    fc.UseAdaptiveSampling(True)
    fc.SetNBins(10)  # number of points to test per parameter

    # use the Feldman-Cousins tool
    interval = 0
    if doFeldmanCousins:
        interval = fc.GetInterval()

    # /
    # / show use of ProfileLikeihoodCalculator utility in ROOT.RooStats
    plc = ROOT.RooStats.ProfileLikelihoodCalculator(data, modelConfig)
    plc.SetTestSize(.1)

    plcInterval = plc.GetInterval()

    # /
    # / show use of ROOT.RooStats.MCMCCalculator( utility in ROOT.RooStats
    mcInt = 0

    if doMCMC:      # turn some messages back on
        ROOT.RooMsgService.instance().setStreamStatus(0, ROOT.kTRUE)
        ROOT.RooMsgService.instance().setStreamStatus(1, ROOT.kTRUE)

        mcmcWatch = ROOT.TStopwatch()
        mcmcWatch.Start()

        axisList = ROOT.RooArgList(deltaMSq, sinSq2theta)
        mc = ROOT.RooStats.MCMCCalculator(data, modelConfig)
        mc.SetNumIters(5000)
        mc.SetNumBurnInSteps(100)
        mc.SetUseKeys(True)
        mc.SetTestSize(.1)
        # set which is x and y axis in posterior histogram
        mc.SetAxes(axisList)
        # mc.SetNumBins(50)
        mcInt = mc.GetInterval()

        mcmcWatch.Stop()
        mcmcWatch.Print()

    ######################
    # make plot of resulting interval

    dataCanvas.cd(4)

    # first plot a small dot for every point tested
    if doFeldmanCousins:
        parameterScan = fc.GetPointsToScan()
        hist = parameterScan.createHistogram("sinSq2theta:deltaMSq", 30, 30)
        #  hist.Draw()
        forContour = hist.Clone()

        # now loop through the points and put a marker if it's in the interval
        tmpPoint = ROOT.RooArgSet()
        # loop over points to test
        for i in range(parameterScan.numEntries()):
            # get a parameter point from the list of points to test.
            tmpPoint = parameterScan.get(i).clone("temp")

            if interval:
                if interval.IsInInterval(tmpPoint):
                    forContour.SetBinContent(hist.FindBin(tmpPoint.getRealValue("sinSq2theta"),
                                                          tmpPoint.getRealValue("deltaMSq")), 1)
                else:
                    forContour.SetBinContent(hist.FindBin(tmpPoint.getRealValue("sinSq2theta"),
                                                          tmpPoint.getRealValue("deltaMSq")), 0)

        if interval:
            level = 0.5
            forContour.SetContour(1, level)
            forContour.SetLineWidth(2)
            forContour.SetLineColor(ROOT.kRed)
            forContour.Draw("cont2,same")

    mcPlot = 0
    if mcInt:
        print "MCMC actual confidence level: ", mcInt.GetActualConfidenceLevel()
        mcPlot = ROOT.RooStats.MCMCIntervalPlot(mcInt)
        mcPlot.SetLineColor(ROOT.kMagenta)
        mcPlot.Draw()

    dataCanvas.Update()

    plotInt = ROOT.RooStats.LikelihoodIntervalPlot(plcInterval)
    plotInt.SetTitle("90% Confidence Intervals")
    if mcInt:
        plotInt.Draw("same")
    else:
        plotInt.Draw()
    dataCanvas.Update()

    # / print timing info
    t.Stop()
    t.Print()


if __name__ == "__main__":
    rs401d_FeldmanCousins()
