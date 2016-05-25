# /
#
# Produces an interval on the mean signal in a number counting
# experiment with known background using the Feldman-Cousins technique.
# date Jan. 2009
# updated June 2010
#
# Using the ROOT.RooStats ROOT.RooStats.FeldmanCousins( tool with 200 bins
# it takes 1 min and the interval is [0.2625, 10.6125]
# with a step size of 0.075.
# ROOT.The interval in Feldman & Cousins's original paper is [.29, 10.81]
#  Phys.Rev.D57:3873-3889,1998.
# /


import ROOT


def rs401c_FeldmanCousins():
    # to time the macro... about 30 s
    t = ROOT.TStopwatch()
    t.Start()

    # make a simple model
    x = ROOT.RooRealVar("x", "", 1, 0, 50)
    mu = ROOT.RooRealVar("mu", "", 2.5, 0, 15)  # with a limit on mu>=0
    b = ROOT.RooConstVar("b", "", 3.)
    mean = ROOT.RooAddition("mean", "", ROOT.RooArgList(mu, b))
    pois = ROOT.RooPoisson("pois", "", x, mean)
    parameters = ROOT.RooArgSet(mu)

    # create a toy dataset
    data = pois.generate(ROOT.RooArgSet(x), 1)
    data.Print("v")

    dataCanvas = ROOT.TCanvas("dataCanvas")
    frame = x.frame()
    data.plotOn(frame)
    frame.Draw()
    dataCanvas.Update()

    w = ROOT.RooWorkspace()
    modelConfig = ROOT.RooStats.ModelConfig("poissonProblem", w)
    modelConfig.SetPdf(pois)
    modelConfig.SetParametersOfInterest(parameters)
    modelConfig.SetObservables(ROOT.RooArgSet(x))
    w.Print()

    # show use of Feldman-Cousins
    fc = ROOT.RooStats.FeldmanCousins(data, modelConfig)
    fc.SetTestSize(.05)  # set size of test
    fc.UseAdaptiveSampling(True)
    # number counting analysis: dataset always has 1 entry with N events
    # observed
    fc.FluctuateNumDataEntries(False)
    fc.SetNBins(100)  # number of points to test per parameter

    # use the Feldman-Cousins tool
    interval = fc.GetInterval()

    # make a canvas for plots
    intervalCanvas = ROOT.TCanvas("intervalCanvas")

    print "is self point in the interval? ", interval.IsInInterval(parameters)

    print "interval is [", interval.LowerLimit(mu), ", ", interval.UpperLimit(mu), "]"

    # using 200 bins it takes 1 min and the answer is
    # interval is [0.2625, 10.6125] with a step size of .075
    # ROOT.The interval in Feldman & Cousins's original paper is [.29, 10.81]
    #  Phys.Rev.D57:3873-3889,1998.

    # No dedicated plotting class yet, do it by hand:

    parameterScan = fc.GetPointsToScan()
    hist = parameterScan.createHistogram("mu", 30)
    hist.Draw()

    tmpPoint = ROOT.RooArgSet()
    mark = []
    # loop over points to test
    for i in range(parameterScan.numEntries()):
        # print "on parameter point ", i, " out of ", parameterScan.numEntries()
        # get a parameter point from the list of points to test.
        tmpPoint = parameterScan.get(i).clone("temp")

        mark.append(ROOT.TMarker(tmpPoint.getRealValue("mu"), 1, 25))
        if interval.IsInInterval(tmpPoint):
            mark[i].SetMarkerColor(ROOT.kBlue)
        else:
            mark[i].SetMarkerColor(ROOT.kRed)

        mark[i].Draw("s")
        # delete tmpPoint
        #    delete mark

    intervalCanvas.SaveAs("rs401c_FeldmanCousins.png")
    t.Stop()
    t.Print()


if __name__ == "__main__":
    rs401c_FeldmanCousins()
