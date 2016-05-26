#####################################
#
# ROOT.RooStats tutorial macro #500a
# 2009/08 - Nils Ruthmann, Schott
#
# Prepare a workspace (stored in a ROOT file) containing a models,
# data and other objects needed to run statistical classes in
# ROOT.RooStats.
#
# In self macro a PDF model is built for a counting analysis.  A
# certain number of events are observed (self can be enforced or left
# free) while a number of background events is expected.  In self
# macro, systematic uncertainty is considered (see
# rs500b_PrepareWorkspace_Poisson_withSystematics.C) ROOT.The parameter of
# interest is the signal yield and we assume for it a flat prior.
# All needed objects are stored in a ROOT file (within a ROOT.RooWorkspace
# container); self ROOT file can then be fed as input to various
# statistical methods.
#
# root -q -x -l 'rs500a_PrepareWorkspace_Poisson.C()'
#
#####################################


import ROOT


def rs500a_PrepareWorkspace_Poisson(fileName="WS_Poisson.root", type=1):

    # use a ROOT.RooWorkspace to store the pdf models, informations, of
    # parameters,...
    myWS = ROOT.RooWorkspace("myWS")

    # Observable
    myWS.factory("x[0,0,1]")

    # Pdf in observable,
    myWS.factory("Uniform::sigPdf(x)")
    myWS.factory("Uniform::bkgPdf(x)")
    myWS.factory("SUM::model(S[0,0,60]*sigPdf,B[10]*bkgPdf)")
    # Background only pdf
    myWS.factory("ExtendPdf::modelBkg(bkgPdf,B)")

    # Prior
    myWS.factory("Uniform::priorPOI(S)")

    # Definition of observables and parameters of interest
    myWS.defineSet("observables", "x")
    myWS.defineSet("POI", "S")

    # Generate data
    data = 0
    # binned data with fixed number of events
    if (type == 0):
        data = myWS.pdf("model").generateBinned(
            myWS.set("observables"), myWS.var("S").getVal(), ROOT.RooFit.Name("data"))
    # binned data with Poisson fluctuations
    if (type == 1):
        data = myWS.pdf("model").generateBinned(
            myWS.set("observables"), ROOT.RooFit.Extended(), ROOT.RooFit.Name("data"))
    # Asimov data: binned data without any fluctuations (average case)
    if (type == 2):
        data = myWS.pdf("model").generateBinned(
            myWS.set("observables"), ROOT.RooFit.Name("data"), ROOT.RooFit.ExpectedData())
    getattr(myWS, 'import')(data)

    # control plot of the generated data
    c1 = ROOT.TCanvas("rs500a_PrepareWorkspace_Poisson",
                      "rs500a_PrepareWorkspace_Poisson", 800, 600)
    plot = myWS.var("x").frame()
    data.plotOn(plot)
    plot.Draw()


    myWS.writeToFile(fileName)
    print "\nRooFit model initialized and stored in ", fileName
    c1.SaveAs("rs500a_PrepareWorkspace_Poisson.png")

if __name__ == "__main__":
    rs500a_PrepareWorkspace_Poisson()
