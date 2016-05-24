# /
#
# rs102_hypotestwithshapes for ROOT.RooStats project
# Author: Kyle Cranmer <cranmer@cern.ch>
#
# Modified from version of February 29, 2008
#
# ROOT.This tutorial macro shows a typical search for a particle
# by studying an invariant mass distribution.
# ROOT.The macro creates a simple signal model and two background models,
# which are added to a ROOT.RooWorkspace.
# ROOT.The macro creates a toy dataset, then uses a ROOT.RooStats
# ProfileLikleihoodCalculator to do a hypothesis test of the
# background-only and signal+background hypotheses.
# In self example, uncertainties are not taken into account, but
# normalization uncertainties are.
#
# /

# use self order for safety on library loading
import ROOT

# ____________________________________


def rs102_hypotestwithshapes():
    # ROOT.The main macro.

    # Create a workspace to manage the project.
    wspace = ROOT.RooWorkspace("myWS")

    # add the signal and background models to the workspace
    AddModel(wspace)

    # add some toy data to the workspace
    AddData(wspace)

    # inspect the workspace if you wish
    #  wspace.Print()

    # do the hypothesis test
    DoHypothesisTest(wspace)

    # make some plots
    MakePlots(wspace)


# ____________________________________
def AddModel(wks):
    # Make models for signal (Higgs) and background (Z+jets and QCD)
    # In real life, part requires an intellegent modeling
    # of signal and background -- self is only an example.

    # set range of observable
    lowRange = 60
    highRange = 200

    # make a ROOT.RooRealVar for the observable
    invMass = ROOT.RooRealVar("invMass", "M_{inv}", lowRange, highRange, "GeV")

    # /
    # make a simple signal model.
    mH = ROOT.RooRealVar("mH", "Higgs Mass", 130, 90, 160)
    sigma1 = ROOT.RooRealVar("sigma1", "Width of Gaussian", 12., 2, 100)
    sigModel = ROOT.RooGaussian(
        "sigModel", "Signal Model", invMass, mH, sigma1)
    # we will test self specific mass point for the signal
    mH.setConstant()
    # and we assume we know the mass resolution
    sigma1.setConstant()

    # /
    # make zjj model.  Just like signal model
    mZ = ROOT.RooRealVar("mZ", "Z Mass", 91.2, 0, 100)
    sigma1_z = ROOT.RooRealVar("sigma1_z", "Width of Gaussian", 10., 6, 100)
    zjjModel = ROOT.RooGaussian(
        "zjjModel", "Z+jets Model", invMass, mZ, sigma1_z)
    # we know Z mass
    mZ.setConstant()
    # assume we know resolution too
    sigma1_z.setConstant()

    #######################
    # make QCD model
    a0 = ROOT.RooRealVar("a0", "a0", 0.26, -1, 1)
    a1 = ROOT.RooRealVar("a1", "a1", -0.17596, -1, 1)
    a2 = ROOT.RooRealVar("a2", "a2", 0.018437, -1, 1)
    # a3 = ROOT.RooRealVar("a3", "a3", 0.02, -1, 1)
    qcdModel = ROOT.RooChebychev(
        "qcdModel", "A  Polynomial for QCD", invMass, ROOT.RooArgList(a0, a1, a2))

    # let's assume self shape is known, the normalization is not
    a0.setConstant()
    a1.setConstant()
    a2.setConstant()

    #######################
    # combined model

    # Setting the fraction of Zjj to be 40% for initial guess.
    fzjj = ROOT.RooRealVar(
        "fzjj", "fraction of zjj background events", .4, 0., 1)

    # Set the expected fraction of signal to 20%.
    fsigExpected = ROOT.RooRealVar(
        "fsigExpected", "expected fraction of signal events", .2, 0., 1)
    fsigExpected.setConstant()  # use mu as main parameter, fix self.

    # Introduce mu: the signal strength in units of the expectation.
    # eg. mu = 1 is the SM, mu = 0 is no signal, mu=2 is 2x the SM
    mu = ROOT.RooRealVar(
        "mu", "signal strength in units of SM expectation", 1, 0., 2)

    # Introduce ratio of signal efficiency to nominal signal efficiency.
    # ROOT.This is useful if you want to do limits on cross section.
    ratioSigEff = ROOT.RooRealVar(
        "ratioSigEff", "ratio of signal efficiency to nominal signal efficiency", 1., 0., 2)
    ratioSigEff.setConstant(ROOT.kTRUE)

    # finally the signal fraction is the product of the terms above.
    fsig = ROOT.RooProduct("fsig", "fraction of signal events",
                           ROOT.RooArgList(mu, ratioSigEff, fsigExpected))

    # full model
    model = ROOT.RooAddPdf("model", "sig+zjj+qcd background shapes",
                           ROOT.RooArgList(sigModel, zjjModel, qcdModel), ROOT.RooArgList(fsig, fzjj))

    # interesting for debugging and visualizing the model
    #  model.printCompactTree("", "fullModel.txt")
    #  model.graphVizTree("fullModel.dot")

    getattr(wks, 'import')(model)


# ____________________________________
def AddData(wks):  # Add a toy dataset

    nEvents = 150
    model = wks.pdf("model")
    invMass = wks.var("invMass")

    data = model.generate(ROOT.RooArgSet(invMass), nEvents)

    getattr(wks, 'import')(data, ROOT.RooFit.Rename("data"))


# ____________________________________
def DoHypothesisTest(wks):

    # Use a ROOT.RooStats ProfileLikleihoodCalculator to do the hypothesis
    # test.
    model = ROOT.RooStats.ModelConfig()
    model.SetWorkspace(wks)
    model.SetPdf("model")

    # plc.SetData("data")

    plc = ROOT.RooStats.ProfileLikelihoodCalculator()
    plc.SetData((wks.data("data")))

    # here we explicitly set the value of the parameters for the null.
    # We want no signal contribution, eg. mu = 0
    mu = wks.var("mu")
    #   nullParams = ROOT.RooArgSet("nullParams")
    #   nullParams.addClone(mu)
    poi = ROOT.RooArgSet(mu)
    nullParams = ROOT.RooArgSet(poi.snapshot())
    nullParams.setRealValue("mu", 0)

    # plc.SetNullParameters(nullParams)
    plc.SetModel(model)
    # NOTE: using snapshot will import nullparams
    # in the WS and merge with existing "mu"
    # model.SetSnapshot(nullParams)

    # use instead setNuisanceParameters
    plc.SetNullParameters(nullParams)

    # We get a HypoTestResult out of the calculator, we can query it.
    htr = plc.GetHypoTest()
    print "-------------------------------------------------"
    print "The p-value for the null is ", htr.NullPValue()
    print "Corresponding to a signifcance of ", htr.Significance()
    print "-------------------------------------------------\n\n"


# ____________________________________
def MakePlots(wks):
    # Make plots of the data and the best fit model in two cases:
    # first the signal+background case
    # second the background-only case.

    # get some things out of workspace
    model = wks.pdf("model")
    sigModel = wks.pdf("sigModel")
    zjjModel = wks.pdf("zjjModel")
    qcdModel = wks.pdf("qcdModel")

    mu = wks.var("mu")
    invMass = wks.var("invMass")
    data = wks.data("data")

    #############################
    # Make plots for the Alternate hypothesis, eg. let mu float

    mu.setConstant(ROOT.kFALSE)

    model.fitTo(data, ROOT.RooFit.Save(ROOT.kTRUE),
                ROOT.RooFit.Minos(ROOT.kFALSE), ROOT.RooFit.Hesse(ROOT.kFALSE), ROOT.RooFit.PrintLevel(-1))

    # plot sig candidates, model, individual componenets
    cdata = ROOT.TCanvas()
    frame = invMass.frame()
    data.plotOn(frame)
    model.plotOn(frame)
    ras_sigModel = ROOT.RooArgSet(sigModel)
    ras_zjjModel = ROOT.RooArgSet(zjjModel)
    ras_qcdModel = ROOT.RooArgSet(qcdModel)
    model.plotOn(frame, ROOT.RooFit.Components(ras_sigModel),
                 ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))
    model.plotOn(frame, ROOT.RooFit.Components(ras_zjjModel),
                 ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlack))
    model.plotOn(frame, ROOT.RooFit.Components(ras_qcdModel),
                 ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen))

    frame.SetTitle("An example fit to the signal + background model")
    frame.Draw()
    cdata.SaveAs("rs102_hypotestwithshapes_alternateFit.png")

    #############################
    # Do Fit to the Null hypothesis.  Eg. fix mu=0

    mu.setVal(0)  # set signal fraction to 0
    mu.setConstant(ROOT.kTRUE)  # set constant

    model.fitTo(data, ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Minos(
        ROOT.kFALSE), ROOT.RooFit.Hesse(ROOT.kFALSE), ROOT.RooFit.PrintLevel(-1))

    # plot signal candidates with background model and components
    cbkgonly = ROOT.TCanvas()
    xframe2 = invMass.frame()
    data.plotOn(xframe2, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    model.plotOn(xframe2)
    model.plotOn(xframe2, ROOT.RooFit.Components(ras_zjjModel),
                 ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlack))
    model.plotOn(xframe2, ROOT.RooFit.Components(ras_qcdModel),
                 ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen))

    xframe2.SetTitle("An example fit to the background-only model")
    xframe2.Draw()

    cbkgonly.SaveAs("rs102_hypotestwithshapes_nullFit.png")


if __name__ == "__main__":
    rs102_hypotestwithshapes()
