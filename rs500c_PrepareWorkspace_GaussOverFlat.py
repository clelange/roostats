#####################################
#
# ROOT.RooStats tutorial macro #500c
# 2009/08 - Nils Ruthmann, Schott
#
# Prepare a workspace (stored in a ROOT file) containing a models,
# data and other objects needed to run statistical classes in
# ROOT.RooStats.
#
# In self macro a PDF model is built assuming signal has a Gaussian
# PDF and the background a flat PDF.  ROOT.The parameter of interest is
# the signal yield and we assume for it a flat prior.  In self macro,
# no systematic uncertainty is considered (see
# rs500d_PrepareWorkspace_GaussOverFlat_withSystematics.C).  All needed
# objects are stored in a ROOT file (within a ROOT.RooWorkspace
# container); self ROOT file can then be fed as input to various
# statistical methods.
#
# root -q -x -l 'rs500c_PrepareWorkspace_GaussOverFlat.C()'
#
#####################################


import ROOT


# prepare the workspace
# type = 0 : unbinned data with total number of events fluctuating using Poisson statistics
# type = 1 : binned data with total number of events fluctuating using Poisson statistics
# type = 2 : binned data without any bin-bin fuctuation (Asimov data)

def rs500c_PrepareWorkspace_GaussOverFlat(fileName="WS_GaussOverFlat.root", type=1):
    # use a ROOT.RooWorkspace to store the pdf models, informations, of
    # parameters,...
    myWS = ROOT.RooWorkspace("myWS")

    # Observable
    myWS.factory("mass[0,500]")

    # Pdf in observable,
    myWS.factory("Gaussian::sigPdf(mass,200,50)")
    myWS.factory("Uniform::bkgPdf(mass)")
    myWS.factory("SUM::model(S[20,0,60]*sigPdf,B[10]*bkgPdf)")
    # Background only pdf
    myWS.factory("ExtendPdf::modelBkg(bkgPdf,B)")
    # Priors
    myWS.factory("Uniform::priorPOI(S)")

    # Definition of observables and parameters of interest
    myWS.defineSet("observables", "mass")
    myWS.defineSet("POI", "S")

    # Generate data
    data = 0
    # unbinned data with Poisson fluctuations for total number of events
    if (type == 0):
        data = myWS.pdf("model").generate(myWS.set("observables"),
                                          ROOT.RooFit.Extended(), ROOT.RooFit.Name("data"))
    # binned data with Poisson fluctuations for total number of events
    if (type == 1):
        data = myWS.pdf("model").generateBinned(
            myWS.set("observables"), ROOT.RooFit.Extended(), ROOT.RooFit.Name("data"))
    # binned without any fluctuations (average case)
    if (type == 2):
        data = myWS.pdf("model").generateBinned(
            myWS.set("observables"), ROOT.RooFit.Name("data"), ROOT.RooFit.ExpectedData())

    getattr(myWS, 'import')(data)

    myWS.writeToFile(fileName)
    print "\nRooFit model initialized and stored in ", fileName

    # control plot of the generated data
    c1 = ROOT.TCanvas("rs500c_PrepareWorkspace_GaussOverFlat",
                      "rs500c_PrepareWorkspace_GaussOverFlat", 800, 600)
    plot = myWS.var("mass").frame()
    data.plotOn(plot)
    plot.Draw()
    c1.SaveAs("rs500c_PrepareWorkspace_GaussOverFlat.png")


if __name__ == "__main__":
    rs500c_PrepareWorkspace_GaussOverFlat()
