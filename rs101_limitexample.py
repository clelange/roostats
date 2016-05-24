# /
#
# 'Limit Example' ROOT.RooStats tutorial macro #101
# author: Kyle Cranmer
# date June. 2009
#
# ROOT.This tutorial shows an example of creating a simple
# model for a number counting experiment with uncertainty
# on both the background rate and signal efficeincy. We then
# use a Confidence Interval Calculator to set a limit on the signal.
#
#
# /


import ROOT


def rs101_limitexample():
    # /
    # An example of setting a limit in a number counting experiment with uncertainty on background and signal
    # /

    # to time the macro
    t = ROOT.TStopwatch()
    t.Start()

    # /
    # ROOT.The Model building stage
    # /
    wspace = ROOT.RooWorkspace()
    # counting model
    wspace.factory(
        "Poisson::countingModel(obs[150,0,300], sum(s[50,0,120]*ratioSigEff[1.,0,3.],b[100]*ratioBkgEff[1.,0.,3.]))")
    #  wspace.factory("Gaussian.sigConstraint(ratioSigEff,1,0.05)"); # 5% signal efficiency uncertainty
    # wspace.factory("Gaussian.bkgConstraint(ratioBkgEff,1,0.1)"); # 10%
    # background efficiency uncertainty
    # 5% signal efficiency uncertainty
    wspace.factory("Gaussian::sigConstraint(gSigEff[1,0,3],ratioSigEff,0.05)")
    # 10% background efficiency uncertainty
    wspace.factory("Gaussian::bkgConstraint(gSigBkg[1,0,3],ratioBkgEff,0.2)")
    # product of terms
    wspace.factory(
        "PROD::modelWithConstraints(countingModel,sigConstraint,bkgConstraint)")
    wspace.Print()

    modelWithConstraints = wspace.pdf("modelWithConstraints")  # get the model
    obs = wspace.var("obs")  # get the observable
    s = wspace.var("s")  # get the signal we care about
    # get the background and set it to a constant.  Uncertainty included in
    # ratioBkgEff
    b = wspace.var("b")
    b.setConstant()

    # get uncertaint parameter to constrain
    ratioSigEff = wspace.var("ratioSigEff")
    # get uncertaint parameter to constrain
    ratioBkgEff = wspace.var("ratioBkgEff")
    # need to constrain these in the fit (should change default behavior)
    constrainedParams = ROOT.RooArgSet(ratioSigEff, ratioBkgEff)

    # global observables for signal efficiency
    gSigEff = wspace.var("gSigEff")
    # global obs for background efficiency
    gSigBkg = wspace.var("gSigBkg")
    gSigEff.setConstant()
    gSigBkg.setConstant()

    # Create an example dataset with 160 observed events
    obs.setVal(160.)
    data = ROOT.RooDataSet("exampleData", "exampleData", ROOT.RooArgSet(obs))
    data.add(ROOT.RooArgSet(obs))

    all = ROOT.RooArgSet(s, ratioBkgEff, ratioSigEff)

    # not necessary
    modelWithConstraints.fitTo(data, ROOT.RooFit.Constrain(
        ROOT.RooArgSet(ratioSigEff, ratioBkgEff)))

    # Now let's make some confidence intervals for s, parameter of interest
    paramOfInterest = ROOT.RooArgSet(s)

    modelConfig = ROOT.RooStats.ModelConfig(wspace)
    modelConfig.SetPdf(modelWithConstraints)
    modelConfig.SetParametersOfInterest(paramOfInterest)
    modelConfig.SetNuisanceParameters(constrainedParams)
    modelConfig.SetObservables(ROOT.RooArgSet(obs))
    modelConfig.SetGlobalObservables(ROOT.RooArgSet(gSigEff, gSigBkg))
    modelConfig.SetName("ModelConfig")
    getattr(wspace, 'import')(modelConfig)
    getattr(wspace, 'import')(data)
    wspace.SetName("w")
    wspace.writeToFile("rs101_ws.root")

    # First, let's use a Calculator based on the Profile Likelihood Ratio
    # ProfileLikelihoodCalculator plc(*data, *modelWithConstraints,
    # paramOfInterest)
    plc = ROOT.RooStats.ProfileLikelihoodCalculator(data, modelConfig)
    plc.SetTestSize(.05)
    lrint = plc.GetInterval()  # that was easy.

    # Let's make a plot
    dataCanvas = ROOT.TCanvas("dataCanvas")
    dataCanvas.Divide(2, 1)

    dataCanvas.cd(1)
    plotInt = ROOT.RooStats.LikelihoodIntervalPlot(lrint)
    plotInt.SetTitle("Profile Likelihood Ratio and Posterior for S")
    plotInt.Draw()

    # Second, a Calculator based on the Feldman Cousins technique
    fc = ROOT.RooStats.FeldmanCousins(data, modelConfig)
    fc.UseAdaptiveSampling(True)
    # number counting analysis: dataset always has 1 entry with N events
    # observed
    fc.FluctuateNumDataEntries(False)
    fc.SetNBins(100)  # number of points to test per parameter
    fc.SetTestSize(.05)
    #  fc.SaveBeltToFile(True); # optional
    fcint = 0
    fcint = fc.GetInterval()  # that was easy.

    fit = modelWithConstraints.fitTo(data, ROOT.RooFit.Save(True))

    # ROOT.Third, a Calculator based on Markov Chain monte carlo
    # Before configuring the calculator, let's make a ProposalFunction
    # that will achieve a high acceptance rate
    ph = ROOT.RooStats.ProposalHelper()
    ph.SetVariables(ROOT.RooArgSet(fit.floatParsFinal()))
    ph.SetCovMatrix(fit.covarianceMatrix())
    ph.SetUpdateProposalParameters(True)
    ph.SetCacheSize(100)
    pdfProp = ph.GetProposalFunction()  # that was easy

    mc = ROOT.RooStats.MCMCCalculator(data, modelConfig)
    mc.SetNumIters(20000)  # steps to propose in the chain
    mc.SetTestSize(.05)  # 95% CL
    mc.SetNumBurnInSteps(40)  # ignore first N steps in chain as "burn in"
    mc.SetProposalFunction(pdfProp)
    mc.SetLeftSideTailFraction(0.5)  # find a "central" interval
    mcInt = mc.GetInterval()  # that was easy

    # Get Lower and Upper limits from Profile Calculator
    print "Profile lower limit s = ", lrint.LowerLimit(s)
    print "Profile upper limit s = ", lrint.UpperLimit(s)

    # Get Lower and Upper limits from FeldmanCousins with profile construction
    if fcint:
        fcul = fcint.UpperLimit(s)
        fcll = fcint.LowerLimit(s)
        print "FC lower limit s = ", fcll
        print "FC upper limit s = ", fcul
        fcllLine = ROOT.TLine(fcll, 0, fcll, 1)
        fculLine = ROOT.TLine(fcul, 0, fcul, 1)
        fcllLine.SetLineColor(ROOT.kRed)
        fculLine.SetLineColor(ROOT.kRed)
        fcllLine.Draw("same")
        fculLine.Draw("same")
        dataCanvas.Update()

    # Plot MCMC interval and print some statistics
    mcPlot = ROOT.RooStats.MCMCIntervalPlot(mcInt)
    mcPlot.SetLineColor(ROOT.kMagenta)
    mcPlot.SetLineWidth(2)
    mcPlot.Draw("same")

    mcul = mcInt.UpperLimit(s)
    mcll = mcInt.LowerLimit(s)
    print "MCMC lower limit s = ", mcll
    print "MCMC upper limit s = ", mcul
    print "MCMC Actual confidence level: ", mcInt.GetActualConfidenceLevel()

    # 3-d plot of the parameter points
    dataCanvas.cd(2)
    # also plot the points in the markov chain
    chainData = mcInt.GetChainAsDataSet()

    assert(chainData)
    print "plotting the chain data - nentries = ", chainData.numEntries()
    chain = ROOT.RooStats.GetAsTTree(
        "chainTreeData", "chainTreeData", chainData)
    assert(chain)
    chain.SetMarkerStyle(6)
    chain.SetMarkerColor(ROOT.kRed)

    chain.Draw("s:ratioSigEff:ratioBkgEff", "nll_MarkovChain_local_",
               "box")  # 3-d box proporional to posterior

    # the points used in the profile construction
    parScanData = fc.GetPointsToScan()
    assert(parScanData)
    print "plotting the scanned points used in the frequentist construction - npoints = ", parScanData.numEntries()
    # getting the tree and drawing it -crashes (very strange....)
    # parameterScan =  ROOT.RooStats.GetAsTTree("parScanTreeData", "parScanTreeData",*parScanData)
    # assert(parameterScan)
    # parameterScan.Draw("s:ratioSigEff:ratioBkgEff", "", "candle goff")
    gr = ROOT.TGraph2D(parScanData.numEntries())
    for ievt in range(parScanData.numEntries()):
        evt = parScanData.get(ievt)
        x = evt.getRealValue("ratioBkgEff")
        y = evt.getRealValue("ratioSigEff")
        z = evt.getRealValue("s")
        gr.SetPoint(ievt, x, y, z)
        # std.print ievt << "  " << x << "  " << y << "  " << z << std.endl

    gr.SetMarkerStyle(24)
    gr.Draw("P SAME")

    dataCanvas.SaveAs("rs101_limitexample.png")

    # print timing info
    t.Stop()
    t.Print()


if __name__ == "__main__":
    rs101_limitexample()
