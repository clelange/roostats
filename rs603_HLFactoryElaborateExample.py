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


def rs603_HLFactoryElaborateExample():

    # --- Prepare the 2 needed datacards for self example ---

    card_name = "rs603_card_WsMaker.rs"
    ofile = open(card_name, "w")
    ofile.write("// ROOT.The simplest card for combination\n\n")
    ofile.write("gauss1 = Gaussian(x[0,100],mean1[50,0,100],4);\n")
    ofile.write("flat1 = Polynomial(x,0);\n")
    ofile.write(
        "sb_model1 = SUM(nsig1[120,0,300]*gauss1 , nbkg1[100,0,1000]*flat1);\n\n")
    ofile.write("echo In the middle! ;\n\n")
    ofile.write("gauss2 = Gaussian(x,mean2[80,0,100],5);\n")
    ofile.write("flat2 = Polynomial(x,0);\n")
    ofile.write(
        "sb_model2 = SUM(nsig2[90,0,400]*gauss2 , nbkg2[80,0,1000]*flat2);\n\n")
    ofile.write("echo At the end! ;\n")
    ofile.close()

    card_name2 = "rs603_card.rs"
    ofile2 = open(card_name2, "w")
    ofile2.write("// The simplest card for combination\n\n")
    ofile2.write("gauss1 = Gaussian(x[0,100],mean1[50,0,100],4);\n")
    ofile2.write("flat1 = Polynomial(x,0);\n")
    ofile2.write(
        "sb_model1 = SUM(nsig1[120,0,300]*gauss1 , nbkg1[100,0,1000]*flat1);\n\n")
    ofile2.write("echo In the middle! ;\n\n")
    ofile2.write("gauss2 = Gaussian(x,mean2[80,0,100],5);\n")
    ofile2.write("flat2 = Polynomial(x,0);\n")
    ofile2.write(
        "sb_model2 = SUM(nsig2[90,0,400]*gauss2 , nbkg2[80,0,1000]*flat2);\n\n")
    ofile2.write("# include rs603_included_card.rs;\n\n")
    ofile2.write("echo At the end! ;\n")
    ofile2.close()

    card_name3 = "rs603_included_card.rs"
    ofile3 = open(card_name3, "w")
    ofile3.write("echo Now reading the included file! ;\n\n")
    ofile3.write(
        "echo Including datasets in a Workspace in a ROOT.Root file...;\n")
    ofile3.write("data1 = import(rs603_infile.root,\n")
    ofile3.write("               rs603_ws,\n")
    ofile3.write("               data1);\n\n")
    ofile3.write("data2 = import(rs603_infile.root,\n")
    ofile3.write("               rs603_ws,\n")
    ofile3.write("               data2);\n")
    ofile3.close()

    # --- Produce the two separate datasets into a WorkSpace ---

    hlf = ROOT.RooStats.HLFactory(
        "HLFactoryComplexExample", "rs603_card_WsMaker.rs", False)

    x = hlf.GetWs().arg("x")
    pdf1 = hlf.GetWs().pdf("sb_model1")
    pdf2 = hlf.GetWs().pdf("sb_model2")

    w = ROOT.RooWorkspace("rs603_ws")

    data1 = pdf1.generate(ROOT.RooArgSet(x), ROOT.RooFit.Extended())
    data1.SetName("data1")
    getattr(w, 'import')(data1)

    data2 = pdf2.generate(ROOT.RooArgSet(x), ROOT.RooFit.Extended())
    data2.SetName("data2")
    getattr(w, 'import')(data2)

    # --- Write the WorkSpace into a rootfile ---

    outfile = ROOT.TFile.Open("rs603_infile.root", "RECREATE")
    w.Write()
    outfile.Close()

    print "-------------------------------------------------------------------"
    print " ROOT.Rootfile and Workspace prepared"
    print "-------------------------------------------------------------------"

    hlf_2 = ROOT.RooStats.HLFactory(
        "HLFactoryElaborateExample", "rs603_card.rs", False)

    x = hlf_2.GetWs().var("x")
    pdf1 = hlf_2.GetWs().pdf("sb_model1")
    pdf2 = hlf_2.GetWs().pdf("sb_model2")

    hlf_2.AddChannel("model1", "sb_model1", "flat1", "data1")
    hlf_2.AddChannel("model2", "sb_model2", "flat2", "data2")

    data = hlf_2.GetTotDataSet()
    pdf = hlf_2.GetTotSigBkgPdf()
    thecat = hlf_2.GetTotCategory()

    # --- Perform extended ML fit of composite PDF to toy data ---
    pdf.fitTo(data)

    # --- Plot toy data and composite PDF overlaid ---
    c = ROOT.TCanvas("rs603_HLFactoryElaborateExample",
                     "rs603_HLFactoryElaborateExample", 800, 600)
    xframe = x.frame()

    data.plotOn(xframe)
    # thecat.setIndex(0)
    pdf.plotOn(xframe, ROOT.RooFit.Slice(thecat, "0"),
               ROOT.RooFit.ProjWData(ROOT.RooArgSet(thecat), data))

    thecat.setIndex(1)
    pdf.plotOn(xframe, ROOT.RooFit.Slice(thecat, "1"),
               ROOT.RooFit.ProjWData(ROOT.RooArgSet(thecat), data))

    ROOT.gROOT.SetStyle("Plain")

    xframe.Draw()

    c.SaveAs("rs603_HLFactoryElaborateExample.png")


if __name__ == "__main__":
    rs603_HLFactoryElaborateExample()
