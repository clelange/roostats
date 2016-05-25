# /
#
# 'High Level Factory Example' ROOT.RooStats tutorial macro #602
# author: Danilo Piparo
# date August. 2009
#
# ROOT.This tutorial shows an example of creating a combined
# model using the High Level model Factory.
#
#
# /


import ROOT


def rs602_HLFactoryCombinationexample():

    # create a card
    card_name = "HLFactoryCombinationexample.rs"
    ofile = open(card_name, "w")
    ofile.write("// ROOT.The simplest card for combination\n\n")
    ofile.write("gauss1 = Gaussian(x[0,100],mean1[50,0,100],4);\n")
    ofile.write("flat1 = Polynomial(x,0);\n")
    ofile.write(
        "sb_model1 = SUM(nsig1[120,0,300]*gauss1 , nbkg1[100,0,1000]*flat1);\n")
    ofile.write("gauss2 = Gaussian(x,mean2[80,0,100],5);\n")
    ofile.write("flat2 = Polynomial(x,0);\n")
    ofile.write(
        "sb_model2 = SUM(nsig2[90,0,400]*gauss2 , nbkg2[80,0,1000]*flat2);\n")

    ofile.close()

    hlf = ROOT.RooStats.HLFactory(
        "HLFactoryCombinationexample", card_name, False)

    hlf.AddChannel("model1", "sb_model1", "flat1")
    hlf.AddChannel("model2", "sb_model2", "flat2")
    pdf = hlf.GetTotSigBkgPdf()
    thecat = hlf.GetTotCategory()
    x = hlf.GetWs().arg("x")

    data = pdf.generate(ROOT.RooArgSet(x, thecat), ROOT.RooFit.Extended())

    # --- Perform extended ML fit of composite PDF to toy data ---
    pdf.fitTo(data)

    # --- Plot toy data and composite PDF overlaid ---
    c = ROOT.TCanvas("rs602_HLFactoryCombinationexample",
                     "rs602_HLFactoryCombinationexample", 800, 600)
    xframe = x.frame()

    data.plotOn(xframe)

    # thecat.setIndex(0)
    pdf.plotOn(xframe, ROOT.RooFit.Slice(thecat, "0"),
               ROOT.RooFit.ProjWData(ROOT.RooArgSet(thecat), data))

    # thecat.setIndex(1)
    pdf.plotOn(xframe, ROOT.RooFit.Slice(thecat, "1"),
               ROOT.RooFit.ProjWData(ROOT.RooArgSet(thecat), data))

    ROOT.gROOT.SetStyle("Plain")
    xframe.Draw()

    c.SaveAs("rs602_HLFactoryCombinationexample.png")


if __name__ == "__main__":
    rs602_HLFactoryCombinationexample()
