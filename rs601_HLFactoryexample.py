# /
#
# 'High Level Factory Example' ROOT.RooStats tutorial macro #601
# author: Danilo Piparo
# date August. 2009
#
# ROOT.This tutorial shows an example of creating a simple
# model using the High Level model Factory.
#
#
# /


import ROOT


def rs601_HLFactoryexample():
    # --- Build the datacard and dump to file---

    card_name = "HLFactoryexample.rs"
    ofile = open(card_name, 'w')
    ofile.write("// The simplest card\n\n")
    ofile.write(
        "gauss = Gaussian(mes[5.20,5.30],mean[5.28,5.2,5.3],width[0.0027,0.001,1]);\n")
    ofile.write("argus = ArgusBG(mes,5.291,argpar[-20,-100,-1]);\n")
    ofile.write(
        "summe = SUM(nsig[200,0,10000]*gauss,nbkg[800,0,10000]*argus);\n\n")

    ofile.close()

    hlf = ROOT.RooStats.HLFactory("HLFactoryexample", card_name, False)

    # --- ROOT.Take elements out of the internal workspace ---

    w = hlf.GetWs()

    mes = w.arg("mes")
    summe = w.pdf("summe")
    argus = w.pdf("argus")
#    mean = dynamic_cast<RooRealVar*>(w.arg("mean"))
#    argpar = dynamic_cast<RooRealVar*>(w.arg("argpar"))

    # --- Generate a toyMC sample from composite PDF ---
    data = summe.generate(ROOT.RooArgSet(mes), 2000)

    # --- Perform extended ML fit of composite PDF to toy data ---
    summe.fitTo(data)

    # --- Plot toy data and composite PDF overlaid ---
    c = ROOT.TCanvas("rs601_HLFactoryexample",
                     "rs601_HLFactoryexample", 800, 400)
    mesframe = mes.frame()
    data.plotOn(mesframe)
    summe.plotOn(mesframe)
    ras_argus = ROOT.RooArgSet(argus)
    summe.plotOn(mesframe, ROOT.RooFit.Components(
        ras_argus), ROOT.RooFit.LineStyle(ROOT.kDashed))

    ROOT.gROOT.SetStyle("Plain")
    mesframe.Draw()

    c.SaveAs("rs601_HLFactoryexample.png")


if __name__ == "__main__":
    rs601_HLFactoryexample()
