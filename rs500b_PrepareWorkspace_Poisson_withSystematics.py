#####################################
#
# ROOT.RooStats tutorial macro #500b
# 2009/08 - Nils Ruthmann, Schott
#
# Prepare a workspace (stored in a ROOT file) containing a models,
# data and other objects needed to run statistical classes in
# ROOT.RooStats.
#
# In self macro a PDF model is built for a counting analysis.  A
# certain number of events are observed (self can be enforced or left
# free) while a number of background events is expected.  It is also
# assumed there is a systematic uncertainty on the number of expected
# background events.  ROOT.The parameter of interest is the signal yield
# and we assume for it a flat prior.  All needed objects are stored
# in a ROOT file (within a ROOT.RooWorkspace container); self ROOT file
# can then be fed as input to various statistical methods.
#
# root -q -x -l 'rs500b_PrepareWorkspace_Poisson_withSystematics.C()'
#
#####################################


import ROOT


# prepare the workspace
# type = 0 : binned data with fixed number of total events
# type = 1 : binned data with N with POisson fluctuations
# type = 2 : binned data without any bin-by bin fluctuations (Asimov data)

def rs500b_PrepareWorkspace_Poisson_withSystematics(fileName="WS_Poisson_withSystematics.root", type=1):

    # use a ROOT.RooWorkspace to store the pdf models, informations, of
    # parameters,...
    myWS = ROOT.RooWorkspace("myWS")

    # Observable
    myWS.factory("x[0,0,1]")

    # Pdf in observable,
    myWS.factory("Uniform::sigPdf(x)")
    myWS.factory("Uniform::bkgPdf(x)")
    myWS.factory("SUM::model(S[100,0,1500]*sigPdf,B[1000,0,3000]*bkgPdf)")

    # Background only pdf
    myWS.factory("ExtendPdf::modelBkg(bkgPdf,B)")

    # Priors
    myWS.factory("Gaussian::priorNuisance(B,1000,200)")
    myWS.factory("Uniform::priorPOI(S)")

    # Definition of observables and parameters of interest
    myWS.defineSet("observables", "x")
    myWS.defineSet("POI", "S")
    myWS.defineSet("parameters", "B")

    # Generate data
    data = 0
    # binned data with fixed number of events
    if (type == 0):
        data = myWS.pdf("model").generateBinned(myWS.set("observables"),
                                                myWS.var("S").getVal(), ROOT.RooFit.Name("data"))
    # binned data with Poisson fluctuations
    if (type == 1):
        data = myWS.pdf("model").generateBinned(myWS.set("observables"),
                                                ROOT.RooFit.Extended(), ROOT.RooFit.Name("data"))
    # Asimov data: binned data without any fluctuations (average case)
    if (type == 2):
        data = myWS.pdf("model").generateBinned(myWS.set("observables"),
                                                ROOT.RooFit.Name("data"), ROOT.RooFit.ExpectedData())
    getattr(myWS, 'import')(data)

    myWS.writeToFile(fileName)
    print "\nRooFit model initialized and stored in ", fileName

    # control plot of the generated data
    c1 = ROOT.TCanvas("rs500b_PrepareWorkspace_Poisson_withSystematics",
                      "rs500b_PrepareWorkspace_Poisson_withSystematics", 800, 600)
    plot = myWS.var("x").frame()
    data.plotOn(plot)
    plot.Draw()
    c1.SaveAs("rs500b_PrepareWorkspace_Poisson_withSystematics.png")


if __name__ == "__main__":
    rs500b_PrepareWorkspace_Poisson_withSystematics()
