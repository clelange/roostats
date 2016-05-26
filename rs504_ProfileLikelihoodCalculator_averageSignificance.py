#####################################
#
# ROOT.RooStats tutorial macro #504
# 2009/08 - Nils Ruthmann, Schott
#
# /


import ROOT


def rs504_ProfileLikelihoodCalculator_averageSignificance(fname="WS_GaussOverFlat_withSystematics.root", ntoys=100, outputplot="rs504_ProfileLikelihoodCalculator_averageSignificance.png"):

    t = ROOT.TStopwatch()
    t.Start()

    file = ROOT.TFile(fname)
    my_WS = file.Get("myWS")
    # Import the objects needed
    model_naked = my_WS.pdf("model")
    priorNuisance = my_WS.pdf("priorNuisance")
    paramInterestSet = my_WS.set("POI")
    paramInterest = paramInterestSet.first()
    observable = my_WS.set("observables")
    nuisanceParam = my_WS.set("parameters")

    # If there are nuisance parameters present, their prior to the model
    model = model_naked
    if priorNuisance != 0:
        model = ROOT.RooProdPdf(
            "constrainedModel", "Model with nuisance parameters", model_naked, priorNuisance)
        # From now work with the product of both

    # Save the default values of the parameters:
    parameters = model.getVariables()
    default_parameters = ROOT.RooArgSet("default_parameters")
    it = parameters.createIterator()
    currentparam = it.Next()
    while(type(currentparam) != ROOT.TObject):
        print "param:", currentparam
        default_parameters.addClone(currentparam, False)
        currentparam = it.Next()

    nll_nuisance = 0
    if priorNuisance != 0:
        nll_nuisance = ROOT.RooFormulaVar(
            "nllSyst", "-TMath.Log(@0)", ROOT.RooArgList(priorNuisance))
    else:
        nll_nuisance = ROOT.RooFormulaVar(
            "nllSyst", "0", ROOT.RooArgList(priorNuisance))

    ROOT.RooRandom.randomGenerator().SetSeed(110)

    #--------------------------------------------------------------------
    # ROOMCSTUDY

    # For simplicity use ROOT.RooMCStudy.
    mcs = 0
    if nuisanceParam:
        mcs = ROOT.RooMCStudy(model, observable, ROOT.RooFit.Extended(ROOT.kTRUE),
                              ROOT.RooFit.FitOptions(ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.PrintEvalErrors(-1), ROOT.RooFit.Minos(ROOT.kFALSE)), ROOT.RooFit.Constrain(nuisanceParam))
    else:
        mcs = ROOT.RooMCStudy(model, observable, ROOT.RooFit.Extended(ROOT.kTRUE),
                              ROOT.RooFit.FitOptions(ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minos(ROOT.kFALSE), ROOT.RooFit.PrintEvalErrors(-1)))

    # Adding a module which allows to compute the significance in each toy
    # experiment
    sigModule = ROOT.RooDLLSignificanceMCSModule(paramInterest, 0)
    # If there are nuisance parameters present, should be generated according to their pdf for every toy experiment.
    # self is done using a MCSModule
    mcs.addModule(sigModule)
    mcs.generateAndFit(ntoys)

    signstring = "significance_nullhypo_"
    mcssign_histo = ROOT.RooAbsData.createHistogram(
        mcs.fitParDataSet(), signstring + paramInterest.GetName())

    c2 = ROOT.TCanvas()
    c2.Divide(2, 2)
    c2.cd(1)
    mcssign_histo.Draw()
    c2.cd(2)
    mcs.plotPull(paramInterest).Draw()
    c2.cd(3)
    if my_WS.var("B"):
        mcs.plotParam(my_WS.var("B")).Draw()
        c2.cd(4)
        mcs.plotError(my_WS.var("B")).Draw()

    c2.Print(outputplot)
    print "The average significance after ", ntoys, " toys is: ", mcssign_histo.GetMean()

    # file.Close()
    t.Stop()
    t.Print()


if __name__ == "__main__":
    rs504_ProfileLikelihoodCalculator_averageSignificance()
