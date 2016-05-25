# /
#
# SPlot tutorial
# author: Kyle Cranmer
# date Dec. 2008
#
# ROOT.This tutorial shows an example of using SPlot to unfold two distributions.
# ROOT.The physics context for the example is that we want to know
# the isolation distribution for real electrons from Z events
# and fake electrons from QCD.  Isolation is our 'control' variable
# ROOT.To unfold them, need a model for an uncorrelated variable that
# discriminates between Z and QCD.  ROOT.To do self, use the invariant
# mass of two electrons.  We model the Z with a Gaussian and the QCD
# with a falling exponential.
#
# Note, we don't have real data in self tutorial, need to generate
# toy data.  ROOT.To do that we need a model for the isolation variable for
# both Z and QCD.  ROOT.This is only used to generate the toy data, would
# not be needed if we had real data.
# /


import ROOT


def rs301_splot():

    # Create a workspace to manage the project.
    wspace = ROOT.RooWorkspace("myWS")

    # add the signal and background models to the workspace.
    # Inside self function you will find a discription our model.
    AddModel(wspace)

    # add some toy data to the workspace
    AddData(wspace)

    # inspect the workspace if you wish
    #  wspace.Print()

    # do sPlot.
    # This wil make a dataset with sWeights added for every event.
    DoSPlot(wspace)

    # Make some plots showing the discriminating variable and
    # the control variable after unfolding.
    MakePlots(wspace)


# ____________________________________
def AddModel(ws):
    # Make models for signal (Higgs) and background (Z+jets and QCD)
    # In real life, part requires an intellegent modeling
    # of signal and background -- self is only an example.

    # set range of observable
    lowRange = 00
    highRange = 200

    # make a ROOT.RooRealVar for the observables
    invMass = ROOT.RooRealVar("invMass", "M_{inv}", lowRange, highRange, "GeV")
    isolation = ROOT.RooRealVar("isolation", "isolation", 0., 20., "GeV")

    # /
    # make 2-d model for Z including the invariant mass
    # distribution  and an isolation distribution which we want to
    # unfold from QCD.
    print "make z model"
    # mass model for Z
    mZ = ROOT.RooRealVar("mZ", "Z Mass", 91.2, lowRange, highRange)
    sigmaZ = ROOT.RooRealVar("sigmaZ", "Width of Gaussian", 2, 0, 10, "GeV")
    mZModel = ROOT.RooGaussian("mZModel", "Z+jets Model", invMass, mZ, sigmaZ)
    # we know Z mass
    mZ.setConstant()
    # we leave the width of the Z free during the fit in self example.

    # isolation model for Z.  Only used to generate toy MC.
    # the exponential is of the form exp(c*x).  If we want
    # the isolation to decay an e-fold every R GeV, use
    # c = -1/R.
    zIsolDecayConst = ROOT.RooConstVar("zIsolDecayConst",
                                       "z isolation decay  constant", -1)
    zIsolationModel = ROOT.RooExponential("zIsolationModel", "z isolation model",
                                          isolation, zIsolDecayConst)

    # make the combined Z model
    zModel = ROOT.RooProdPdf("zModel", "4-d model for Z",
                             ROOT.RooArgList(mZModel, zIsolationModel))

    #######################
    # make QCD model

    print "make qcd model"
    # mass model for QCD.
    # the exponential is of the form exp(c*x).  If we want
    # the mass to decay an e-fold every R GeV, use
    # c = -1/R.
    # We can leave self parameter free during the fit.
    qcdMassDecayConst = ROOT.RooRealVar("qcdMassDecayConst",
                                        "Decay  for QCD mass spectrum",
                                        -0.01, -100, 100, "1/GeV")
    qcdMassModel = ROOT.RooExponential("qcdMassModel", "qcd Mass Model",
                                       invMass, qcdMassDecayConst)

    # isolation model for QCD.  Only used to generate toy MC
    # the exponential is of the form exp(c*x).  If we want
    # the isolation to decay an e-fold every R GeV, use
    # c = -1/R.
    qcdIsolDecayConst = ROOT.RooConstVar("qcdIsolDecayConst",
                                         "Et resolution constant", -.1)
    qcdIsolationModel = ROOT.RooExponential("qcdIsolationModel", "QCD isolation model",
                                            isolation, qcdIsolDecayConst)

    # make the 2-d model
    qcdModel = ROOT.RooProdPdf("qcdModel", "2-d model for QCD",
                               ROOT.RooArgList(qcdMassModel, qcdIsolationModel))

    #######################
    # combined model

    # ROOT.These variables represent the number of Z or QCD events
    # ROOT.They will be fitted.
    zYield = ROOT.RooRealVar("zYield", "fitted yield for Z", 50, 0., 1000)
    qcdYield = ROOT.RooRealVar(
        "qcdYield", "fitted yield for QCD", 100, 0., 1000)

    # now make the combined model
    print "make full model"
    model = ROOT.RooAddPdf("model", "z+qcd background models",
                           ROOT.RooArgList(zModel, qcdModel),
                           ROOT.RooArgList(zYield, qcdYield))

    # interesting for debugging and visualizing the model
    model.graphVizTree("fullModel.dot")

    print "import model"

    getattr(ws, 'import')(model)


# ____________________________________
def AddData(ws):  # Add a toy dataset

    # how many events do we want?
    nEvents = 1000

    # get what we need out of the workspace to make toy data
    model = ws.pdf("model")
    invMass = ws.var("invMass")
    isolation = ws.var("isolation")

    # make the toy data
    print "make data set and import to workspace"
    data = model.generate(ROOT.RooArgSet(invMass, isolation), nEvents)

    # import data into workspace
    getattr(ws, 'import')(data, ROOT.RooFit.Rename("data"))


# ____________________________________
def DoSPlot(ws):

    print "Calculate sWeights"

    # get what we need out of the workspace to do the fit
    model = ws.pdf("model")
    zYield = ws.var("zYield")
    qcdYield = ws.var("qcdYield")
    data = ws.data("data")

    # fit the model to the data.
    model.fitTo(data, ROOT.RooFit.Extended())

    # ROOT.The sPlot technique requires that we fix the parameters
    # of the model that are not yields after doing the fit.
    sigmaZ = ws.var("sigmaZ")
    qcdMassDecayConst = ws.var("qcdMassDecayConst")
    sigmaZ.setConstant()
    qcdMassDecayConst.setConstant()

    ROOT.RooMsgService.instance().setSilentMode(True)

    # Now we use the SPlot class to add SWeights to our data set
    # based on our model and our yield variables
    sData = ROOT.RooStats.SPlot("sData", "An SPlot",
                                data, model, ROOT.RooArgList(zYield, qcdYield))

    # Check that our weights have the desired properties

    print "Check SWeights:"

    print "Yield of Z is ", zYield.getVal(), ".  From sWeights it is ", sData.GetYieldFromSWeight("zYield")

    print "Yield of QCD is ", qcdYield.getVal(), ".  From sWeights it is ", sData.GetYieldFromSWeight("qcdYield")

    for i in range(10):
        print "z Weight   ", sData.GetSWeight(i, "zYield"), "   qcd Weight   ", sData.GetSWeight(i, "qcdYield"), "  ROOT.Total Weight   ", sData.GetSumOfEventSWeight(i)

    print ""

    # import self dataset with sWeights
    print "import dataset with sWeights"
    getattr(ws, 'import')(data, ROOT.RooFit.Rename("dataWithSWeights"))


def MakePlots(ws):
    # Here we make plots of the discriminating variable (invMass) after the fit
    # and of the control variable (isolation) after unfolding with sPlot.
    print "make plots"

    # make our canvas
    cdata = ROOT.TCanvas("sPlot", "sPlot demo", 400, 600)
    cdata.Divide(1, 3)

    # get what we need out of the workspace
    model = ws.pdf("model")
    zModel = ws.pdf("zModel")
    qcdModel = ws.pdf("qcdModel")

    isolation = ws.var("isolation")
    invMass = ws.var("invMass")

    # note, get the dataset with sWeights
    data = ws.data("dataWithSWeights")

    # self shouldn't be necessary, to fix something with workspace
    # do self to set parameters back to their fitted values.
    model.fitTo(data, ROOT.RooFit.Extended())

    # plot invMass for data with full model and individual componenets overlayed
    #  cdata = ROOT.TCanvas()
    cdata.cd(1)
    frame = invMass.frame()
    data.plotOn(frame)
    model.plotOn(frame)
    ras_zModel = ROOT.RooArgSet(zModel)
    ras_qcdModel = ROOT.RooArgSet(qcdModel)
    model.plotOn(frame, ROOT.RooFit.Components(
        ras_zModel), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))
    model.plotOn(frame, ROOT.RooFit.Components(
        ras_qcdModel), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen))

    frame.SetTitle("Fit of model to discriminating variable")
    frame.Draw()

    # Now use the sWeights to show isolation distribution for Z and QCD.
    # ROOT.The SPlot class can make self easier, here we demonstrait in more
    # detail how the sWeights are used.  ROOT.The SPlot class should make self
    # very easy and needs some more development.

    # Plot isolation for Z component.
    # Do self by plotting all events weighted by the sWeight for the Z component.
    # ROOT.The SPlot class adds a variable that has the name of the corresponding
    # yield + "_sw".
    cdata.cd(2)

    # create weightfed data set
    dataw_z = ROOT.RooDataSet(data.GetName(), data.GetTitle(),
                           data, data.get(), "", "zYield_sw")

    frame2 = isolation.frame()
    dataw_z.plotOn(frame2, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

    frame2.SetTitle("isolation distribution for Z")
    frame2.Draw()

    # Plot isolation for QCD component.
    # Eg. plot all events weighted by the sWeight for the QCD component.
    # ROOT.The SPlot class adds a variable that has the name of the corresponding
    # yield + "_sw".
    cdata.cd(3)
    dataw_qcd = ROOT.RooDataSet(data.GetName(), data.GetTitle(), data, data.get(), "", "qcdYield_sw")
    frame3 = isolation.frame()
    dataw_qcd.plotOn(frame3, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

    frame3.SetTitle("isolation distribution for QCD")
    frame3.Draw()

    cdata.SaveAs("rs301_splot.png")


if __name__ == "__main__":
    rs301_splot()
