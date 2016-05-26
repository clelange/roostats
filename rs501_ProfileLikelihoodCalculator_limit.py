#####################################
#
# ROOT.RooStats tutorial macro #501
# 2009/08 - Nils Ruthmann, Schott
#
# Show how to run the ROOT.RooStats classes to perform specific tasks. ROOT.The
# ROOT file containing a workspace holding the models, and other
# objects needed to run can be prepared with any of the rs500*.C
# tutorial macros.
#
# Compute with ROOT.RooStats.ProfileLikelihoodCalculator( a 95% CL upper limit on
# the parameter of interest for the given data.
#
#####################################


import ROOT


def rs501_ProfileLikelihoodCalculator_limit(fileName="WS_GaussOverFlat.root"):
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
        print modelTmp, type(modelTmp)
        model = ROOT.RooProdPdf(
            "constrainedModel", "Model with nuisance parameters", ROOT.RooArgList(modelTmp, priorNuisance))

    # Set up the ROOT.RooStats.ProfileLikelihoodCalculator(
    plc = ROOT.RooStats.ProfileLikelihoodCalculator(data, model, POI)
    # ROOT.RooStats.ProfileLikelihoodCalculator( usually make intervals: the
    # 95% CL one-sided upper-limit is the same as the two-sided upper-limit of
    # a 90% CL interval
    plc.SetTestSize(0.10)

    # Pointer to the confidence interval
    model.fitTo(data, ROOT.RooFit.SumW2Error(ROOT.kFALSE))  # <-- problem
    interval = plc.GetInterval()

    # Compute the upper limit: a fit is needed first in order to locate the
    # minimum of the -log(likelihood) and ease the upper limit computation
    model.fitTo(data, ROOT.RooFit.SumW2Error(ROOT.kFALSE))  # <-- problem
    upperLimit = interval.UpperLimit(parameterOfInterest)  # <-- to simplify

    file.Close()

    # Make a plot of the profile-likelihood and confidence interval
    c1 = ROOT.TCanvas("rs501_ProfileLikelihoodCalculator_limit",
                      "rs501_ProfileLikelihoodCalculator_limit", 800, 600)
    plot = ROOT.RooStats.LikelihoodIntervalPlot(interval)
    plot.Draw()

    print "One sided upper limit at 95% CL: ", upperLimit
    c1.SaveAs("rs501_ProfileLikelihoodCalculator_limit.png")


if __name__ == "__main__":
    rs501_ProfileLikelihoodCalculator_limit()
