# /
#
# 'Number Counting Example' ROOT.RooStats tutorial macro #100
# author: Kyle Cranmer
# date Nov. 2008
#
# ROOT.This tutorial shows an example of a combination of
# two searches using number counting with background uncertainty.
#
# ROOT.The macro uses a ROOT.RooStats "factory" to construct a PDF
# that represents the two number counting analyses with background
# uncertainties.  ROOT.The uncertainties are taken into account by
# considering a sideband measurement of a size that corresponds to the
# background uncertainty.  ROOT.The problem has been studied in these references:
#   http:#arxiv.org/abs/physics/0511028
#   http:#arxiv.org/abs/physics/0702156
#   http:#cdsweb.cern.ch/record/1099969?ln=en
#
# After using the factory to make the model, use a ROOT.RooStats
# ROOT.RooStats.ProfileLikelihoodCalculator( for a Hypothesis test and a confidence interval.
# ROOT.The calculator takes into account systematics by eliminating nuisance parameters
# with the profile likelihood.  ROOT.This is equivalent to the method of MINOS.
#
# /


import ROOT
from array import array


######################
# main driver to choose one
def rs_numberCountingCombination(flag=1):
    if flag == 1:
        rs_numberCountingCombination_expected()
    if flag == 2:
        rs_numberCountingCombination_observed()
    if flag == 3:
        rs_numberCountingCombination_observedWithTau()


# /
def rs_numberCountingCombination_expected():

    # /
    # An example of a number counting combination with two channels.
    # We consider both hypothesis testing and the equivalent confidence interval.
    # /

    # /
    # ROOT.The Model building stage
    # /

    # Step 1, arrays with signal & bkg expectations and background
    # uncertainties
    s = array('d', [20., 10.])           # expected signal
    b = array('d', [100., 100.])         # expected background
    db = array('d', [.0100, .0100])     # fractional background uncertainty

    # Step 2, a ROOT.RooStats factory to build a PDF for a
    # number counting combination and add it to the workspace.
    # We need to give the signal expectation to relate the masterSignal
    # to the signal contribution in the individual channels.
    # The model neglects correlations in background uncertainty,
    # but they could be added without much change to the example.
    f = ROOT.RooStats.NumberCountingPdfFactory()
    wspace = ROOT.RooWorkspace()
    f.AddModel(s, 2, wspace, "TopLevelPdf", "masterSignal")

    # Step 3, a RooStats factory to add datasets to the workspace.
    # Step 3a.
    # Add the expected data to the workspace
    f.AddExpData(s, b, db, 2, wspace, "ExpectedNumberCountingData")

    # see below for a printout of the workspace
    wspace.Print()  # uncomment to see structure of workspace

    # /
    # ROOT.The Hypothesis testing stage:
    # /
    # Step 4, the null hypothesis for the calculator
    # Here you need to know the name of the variables corresponding to
    # hypothesis.
    mu = wspace.var("masterSignal")
    poi = ROOT.RooArgSet(mu)
    nullParams = ROOT.RooArgSet("nullParams")
    _ = nullParams.addClone(mu)
    # here we explicitly set the value of the parameters for the null
    nullParams.setRealValue("masterSignal", 0)

    # Step 5, a calculator for doing the hypothesis test.
    # because self is a
    plc = ROOT.RooStats.ProfileLikelihoodCalculator(wspace.data(
        "ExpectedNumberCountingData"), wspace.pdf("TopLevelPdf"), poi, 0.05, nullParams)

    # Step 6, the Calculator to get a HypoTestResult
    htr = plc.GetHypoTest()
    assert(htr != 0)
    print "-------------------------------------------------"
    print "The p-value for the null is ", htr.NullPValue()
    print "Corresponding to a signifcance of ", htr.Significance()
    print "-------------------------------------------------\n\n"

    ''' expected case should return:
     -------------------------------------------------
     ROOT.The p-value for the null is 0.015294
     Corresponding to a signifcance of 2.16239
     -------------------------------------------------
  '''

    #####################
    # Confidence Interval Stage

    # Step 8, we re-use the ROOT.RooStats.ProfileLikelihoodCalculator( to return a confidence interval.
    # We need to specify what are our parameters of interest
    paramsOfInterest = nullParams  # they are the same as before in self case
    plc.SetParameters(paramsOfInterest)
    lrint = plc.GetInterval()  # that was easy.
    lrint.SetConfidenceLevel(0.95)

    # Step 9, a plot of the likelihood ratio and the interval obtained
    # paramsOfInterest.setRealValue("masterSignal",1.)
    # find limits
    lower = lrint.LowerLimit(mu)
    upper = lrint.UpperLimit(mu)

    c = ROOT.TCanvas("rs_numberCountingCombination_expected",
                     "rs_numberCountingCombination_expected", 800, 600)

    lrPlot = ROOT.RooStats.LikelihoodIntervalPlot(lrint)
    lrPlot.SetMaximum(3.)
    lrPlot.Draw()

    c.SaveAs("rs_numberCountingCombination_expected.png")

    # Step 10a. Get upper and lower limits
    print "lower limit on signal = ", lower
    print "upper limit on signal = ", upper

    # Step 10b, if masterSignal=0 is in the interval.
    # Note, is equivalent to the question of a 2-sigma hypothesis test:
    # "is the parameter point masterSignal=0 inside the 95% confidence interval?"
    # Since the signficance of the Hypothesis test was > 2-sigma it should not be:
    # eg. we exclude masterSignal=0 at 95% confidence.
    paramsOfInterest.setRealValue("masterSignal", 0.)
    print "-------------------------------------------------"
    print "Consider self parameter point:"
    paramsOfInterest.first().Print()
    if lrint.IsInInterval(paramsOfInterest):
        print "It IS in the interval."
    else:
        print "It is NOT in the interval."
    print "-------------------------------------------------\n\n"

    # Step 10c, also ask about the parameter point masterSignal=2, is inside
    # the interval.
    paramsOfInterest.setRealValue("masterSignal", 2.)
    print "-------------------------------------------------"
    print "Consider self parameter point:"
    paramsOfInterest.first().Print()
    if lrint.IsInInterval(paramsOfInterest):
        print "It IS in the interval."
    else:
        print "It is NOT in the interval."
    print "-------------------------------------------------\n\n"

    return

    '''
    # Here's an example of what is in the workspace
    #  wspace.Print()
    ROOT.RooWorkspace(NumberCountingWS) Number Counting WS contents

    variables
    ---------
    (x_0,masterSignal,expected_s_0,b_0,y_0,tau_0,x_1,expected_s_1,b_1,y_1,tau_1)

    p.d.f.s
    -------
    ROOT.RooProdPdf.joint[ pdfs=(sigRegion_0,sideband_0,sigRegion_1,sideband_1) ] = 2.20148e-08
    ROOT.RooPoisson.sigRegion_0[ x=x_0 mean=splusb_0 ] = 0.036393
    ROOT.RooPoisson.sideband_0[ x=y_0 mean=bTau_0 ] = 0.00398939
    ROOT.RooPoisson.sigRegion_1[ x=x_1 mean=splusb_1 ] = 0.0380088
    ROOT.RooPoisson.sideband_1[ x=y_1 mean=bTau_1 ] = 0.00398939

    functions
    --------
    ROOT.RooAddition.splusb_0[ set1=(s_0,b_0) set2=() ] = 120
    ROOT.RooProduct.s_0[ compRSet=(masterSignal,expected_s_0) compCSet=() ] = 20
    ROOT.RooProduct.bTau_0[ compRSet=(b_0,tau_0) compCSet=() ] = 10000
    ROOT.RooAddition.splusb_1[ set1=(s_1,b_1) set2=() ] = 110
    ROOT.RooProduct.s_1[ compRSet=(masterSignal,expected_s_1) compCSet=() ] = 10
    ROOT.RooProduct.bTau_1[ compRSet=(b_1,tau_1) compCSet=() ] = 10000

    datasets
    --------
    ROOT.RooDataSet.ExpectedNumberCountingData(x_0,y_0,x_1,y_1)

    embedded precalculated expensive components
    -------------------------------------------
    '''


def rs_numberCountingCombination_observed():

    # /
    # ROOT.The same example with observed data in a main
    # measurement and an background-only auxiliary
    # measurement with a factor tau more background
    # than in the main measurement.

    # /
    # ROOT.The Model building stage
    # /

    # Step 1, arrays with signal & bkg expectations and background uncertainties
    # We still need the expectation to relate signal in different channels
    # with the master signal
    s = array('d', [20., 10.])           # expected signal

    # Step 2, a ROOT.RooStats factory to build a PDF for a
    # number counting combination and add it to the workspace.
    # We need to give the signal expectation to relate the masterSignal
    # to the signal contribution in the individual channels.
    # ROOT.The model neglects correlations in background uncertainty,
    # but they could be added without much change to the example.
    f = ROOT.RooStats.NumberCountingPdfFactory()
    wspace = ROOT.RooWorkspace()
    f.AddModel(s, 2, wspace, "TopLevelPdf", "masterSignal")

    # Step 3, a ROOT.RooStats factory to add datasets to the workspace.
    # Add the observed data to the workspace
    mainMeas = array('d', [123., 117.])      # observed main measurement
    bkgMeas = array('d', [111.23, 98.76])    # observed background
    # observed fractional background uncertainty
    dbMeas = array('d', [.011, .0095])
    f.AddData(mainMeas, bkgMeas, dbMeas, 2,
              wspace, "ObservedNumberCountingData")

    # see below for a printout of the workspace
    #  wspace.Print();  #uncomment to see structure of workspace

    # /
    # ROOT.The Hypothesis testing stage:
    # /
    # Step 4, the null hypothesis for the calculator
    # Here you need to know the name of the variables corresponding to
    # hypothesis.
    mu = wspace.var("masterSignal")
    poi = ROOT.RooArgSet(mu)
    nullParams = ROOT.RooArgSet("nullParams")
    _ = nullParams.addClone(mu)
    # here we explicitly set the value of the parameters for the null
    nullParams.setRealValue("masterSignal", 0)

    # Step 5, a calculator for doing the hypothesis test.
    # because self is a
    plc = ROOT.RooStats.ProfileLikelihoodCalculator(wspace.data("ObservedNumberCountingData"),
                                                    wspace.pdf("TopLevelPdf"), poi, 0.05, nullParams)

    wspace.var("tau_0").Print()
    wspace.var("tau_1").Print()

    # Step 7, the Calculator to get a HypoTestResult
    htr = plc.GetHypoTest()
    print "-------------------------------------------------"
    print "The p-value for the null is ", htr.NullPValue()
    print "Corresponding to a signifcance of ", htr.Significance()
    print "-------------------------------------------------\n\n"

    ''' observed case should return:
     -------------------------------------------------
     ROOT.The p-value for the null is 0.0351669
     Corresponding to a signifcance of 1.80975
     -------------------------------------------------
  '''

    #####################
    # Confidence Interval Stage

    # Step 8, we re-use the ROOT.RooStats.ProfileLikelihoodCalculator( to return a confidence interval.
    # We need to specify what are our parameters of interest
    paramsOfInterest = nullParams  # they are the same as before in self case
    plc.SetParameters(paramsOfInterest)
    lrint = plc.GetInterval()  # that was easy.
    lrint.SetConfidenceLevel(0.95)

    # Step 9c. Get upper and lower limits
    print "lower limit on signal = ",   lrint.LowerLimit(mu)
    print "upper limit on signal = ",   lrint.UpperLimit(mu)


def rs_numberCountingCombination_observedWithTau():

    # /
    # ROOT.The same example with observed data in a main
    # measurement and an background-only auxiliary
    # measurement with a factor tau more background
    # than in the main measurement.

    # /
    # ROOT.The Model building stage
    # /

    # Step 1, arrays with signal & bkg expectations and background uncertainties
    # We still need the expectation to relate signal in different channels
    # with the master signal
    s = array('d', [20., 10.])           # expected signal

    # Step 2, a ROOT.RooStats factory to build a PDF for a
    # number counting combination and add it to the workspace.
    # We need to give the signal expectation to relate the masterSignal
    # to the signal contribution in the individual channels.
    # ROOT.The model neglects correlations in background uncertainty,
    # but they could be added without much change to the example.
    f = ROOT.RooStats.NumberCountingPdfFactory()
    wspace = ROOT.RooWorkspace()
    f.AddModel(s, 2, wspace, "TopLevelPdf", "masterSignal")

    # Step 3, a ROOT.RooStats factory to add datasets to the workspace.
    # Add the observed data to the workspace in the on-off problem.
    mainMeas = array('d', [123., 117.])      # observed main measurement
    sideband = array('d', [11123., 9876.])    # observed sideband
    # ratio of bkg in sideband to bkg in main measurement, experimental design.
    tau = array('d', [100., 100.])
    f.AddDataWithSideband(mainMeas, sideband, tau, 2,
                          wspace, "ObservedNumberCountingDataWithSideband")

    # see below for a printout of the workspace
    #  wspace.Print();  #uncomment to see structure of workspace

    # /
    # ROOT.The Hypothesis testing stage:
    # /
    # Step 4, the null hypothesis for the calculator
    # Here you need to know the name of the variables corresponding to
    # hypothesis.
    mu = wspace.var("masterSignal")
    poi = ROOT.RooArgSet(mu)
    nullParams = ROOT.RooArgSet("nullParams")
    _ = nullParams.addClone(mu)
    # here we explicitly set the value of the parameters for the null
    nullParams.setRealValue("masterSignal", 0)

    # Step 5, a calculator for doing the hypothesis test.
    # because self is a
    plc = ROOT.RooStats.ProfileLikelihoodCalculator(wspace.data("ObservedNumberCountingDataWithSideband"),
                                                    wspace.pdf("TopLevelPdf"), poi, 0.05, nullParams)

    # Step 7, the Calculator to get a HypoTestResult
    htr = plc.GetHypoTest()
    print "-------------------------------------------------"
    print "The p-value for the null is ", htr.NullPValue()
    print "Corresponding to a signifcance of ", htr.Significance()
    print "-------------------------------------------------\n\n"

    ''' observed case should return:
     -------------------------------------------------
     ROOT.The p-value for the null is 0.0352035
     Corresponding to a signifcance of 1.80928
     -------------------------------------------------
  '''

    #####################
    # Confidence Interval Stage

    # Step 8, we re-use the ROOT.RooStats.ProfileLikelihoodCalculator( to return a confidence interval.
    # We need to specify what are our parameters of interest
    paramsOfInterest = nullParams  # they are the same as before in self case
    plc.SetParameters(paramsOfInterest)
    lrint = plc.GetInterval()  # that was easy.
    lrint.SetConfidenceLevel(0.95)

    # Step 9c. Get upper and lower limits
    print "lower limit on signal = ",   lrint.LowerLimit(mu)
    print "upper limit on signal = ",   lrint.UpperLimit(mu)


if __name__ == "__main__":
    rs_numberCountingCombination()
