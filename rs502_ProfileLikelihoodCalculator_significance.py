#####################################
#
# ROOT.RooStats tutorial macro #502
# 2009/08 - Nils Ruthmann, Schott
#
# Show how to run the ROOT.RooStats classes to perform specific tasks. ROOT.The
# ROOT file containing a workspace holding the models, and other
# objects needed to run can be prepared with any of the rs500*.C
# tutorial macros.
#
# Compute with ROOT.RooStats.ProfileLikelihoodCalculator( a 68% CL interval and the
# signal significance for the given data.
#
#####################################


import ROOT


def rs502_ProfileLikelihoodCalculator_significance(fileName="WS_GaussOverFlat.root"):
    # Open the ROOT file and import from the workspace the objects needed for
    # self tutorial
    file = ROOT.TFile(fileName)
    myWS = file.Get("myWS")
    modelTmp = myWS.pdf("model")
    data = myWS.data("data")
    priorNuisance = myWS.pdf("priorNuisance")
    POI = myWS.set("POI")
    parameterOfInterest = (POI.first())
    assert(parameterOfInterest)

    # If there are nuisance parameters, their prior distribution to the full
    # model
    model = modelTmp
    if (priorNuisance != 0):
        model = ROOT.RooProdPdf("constrainedModel", "Model with nuisance parameters", ROOT.RooArgList(modelTmp, priorNuisance))
    # myWS.var("B").setConstant()

    # Set up the ROOT.RooStats.ProfileLikelihoodCalculator(
    plc = ROOT.RooStats.ProfileLikelihoodCalculator(data, model, POI)
    # ROOT.The 68% CL ROOT.RooStats.ProfileLikelihoodCalculator( interval
    # correspond to test size of 0.32
    plc.SetTestSize(0.32)

    # Compute the confidence interval: a fit is needed first in order to
    # locate the minimum of the -log(likelihood) and ease the upper limit
    # computation

    model.fitTo(data, ROOT.RooFit.SumW2Error(ROOT.kFALSE))
    # Pointer to the confidence interval
    interval = plc.GetInterval()

    MLE = parameterOfInterest.getVal()
    lowerLimit = interval.LowerLimit(parameterOfInterest)
    parameterOfInterest.setVal(MLE)
    upperLimit = interval.UpperLimit(parameterOfInterest)
    parameterOfInterest.setVal(MLE)

    # Make a plot of the profile-likelihood and confidence interval
    c1 = ROOT.TCanvas("rs502_ProfileLikelihoodCalculator_significance",
                      "rs502_ProfileLikelihoodCalculator_significance", 800, 600)
    plot = ROOT.RooStats.LikelihoodIntervalPlot(interval)
    plot.Draw()

    # Get the significance using the GetHypoTest function: (a plot is not
    # possible)

    # create a copy of the POI parameters to set the values to zero
    nullparams = ROOT.RooArgSet()
    nullClone_ = nullparams.addClone(parameterOfInterest)
    nullparams.first().setVal(0)
    plc.SetNullParameters(nullparams)

    testresult = plc.GetHypoTest()
    significance = testresult.Significance()

    #   # Another way is to use directly
    #   # the profile-log-likelihood function
    #   parameterOfInterest.setVal(MLE)
    #   profile = interval.GetLikelihoodRatio()
    #   profile.getVal()
    #   #Go to the Bkg Hypothesis
    #   myarg = (ROOT.RooRealVar *) profile.getVariables().find(parameterOfInterest.GetName()); # <-- cloned!
    #   myarg.setVal(0)
    #   delta_nll = profile.getVal()
    #   significance = sqrt(fabs(2*delta_nll))
    #   if (delta_nll<0) significance *= -1


    print "68% CL interval: [ ", lowerLimit, " ; ", upperLimit, " ]\n"
    print "significance estimation: ", significance

    c1.SaveAs("rs502_ProfileLikelihoodCalculator_significance.png")

    # closing the file will delete the workspace
    file.Close()


if __name__ == "__main__":
    rs502_ProfileLikelihoodCalculator_significance()
