import ROOT


def rs505_HybridCalculator_significance(fname="WS_GaussOverFlat_withSystematics.root", ntoys=5000, outputplot="rs505_HybridCalculator_significance.png"):
    ROOT.RooRandom.randomGenerator().SetSeed(100)
    t = ROOT.TStopwatch()
    t.Start()
    file = ROOT.TFile(fname)
    my_WS = file.Get("myWS")
    if (not my_WS):
        return
    # Import the objects needed
    model = my_WS.pdf("model")
    data = my_WS.data("data")
    priorNuisance = my_WS.pdf("priorNuisance")
    # ROOT.RooArgSet* paramInterestSet=my_WS.set("POI")
    # paramInterest = (ROOT.RooRealVar*) paramInterestSet.first()
    modelBkg = my_WS.pdf("modelBkg")
    # ROOT.RooArgSet* observable=my_WS.set("observables")
    nuisanceParam = my_WS.set("parameters")

    print type(model)
    hc = ROOT.RooStats.HybridCalculator(data, model, modelBkg)
    hc.SetNumberOfToys(ntoys)
    usepriors = False

    if priorNuisance != 0:
        hc.UseNuisance(ROOT.kTRUE)
        hc.SetNuisancePdf(priorNuisance)
        usepriors = True
        print "The following nuisance parameters are taken into account:"
        nuisanceParam.Print()
        hc.SetNuisanceParameters(nuisanceParam)
    else:
        hc.UseNuisance(ROOT.kFALSE)

    hcresult = hc.GetHypoTest()
    clsb_data = hcresult.CLsplusb()
    clb_data = hcresult.CLb()
    cls_data = hcresult.CLs()
    data_significance = hcresult.Significance()

    print "CL_b:", clb_data
    print "CL_s:", cls_data
    print "CL_sb:", clsb_data
    print "significance:", data_significance

    hcPlot = hcresult.GetPlot("hcPlot", "p Values Plot", 100)
    c1 = ROOT.TCanvas()
    c1.cd()
    hcPlot.Draw()
    c1.Print(outputplot)
    file.Close()
    t.Stop()
    t.Print()


if __name__ == "__main__":
    rs505_HybridCalculator_significance()
