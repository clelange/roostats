
import ROOT
def ROOT.TestNonCentral(self):
  ROOT.RooWorkspace w("w")
  # k <2, use sum
  w.factory("NonCentralChiSquare.nc(x[0,50],k[1.99,0,5],lambda[5])")
  # kk > 2 can use bessel
  w.factory("NonCentralChiSquare.ncc(x,kk[2.01,0,5],lambda)")
  # kk > 2, sum
  w.factory("NonCentralChiSquare.nccc(x,kk,lambda)")
  ((ROOT.RooNonCentralChiSquare*)w.pdf("nccc")).SetForceSum(True)

  # a normal "central" chi-square for comparision when lambda.0
  w.factory("ChiSquarePdf.cs(x,k)")

  #w.var("kk").setVal(4.); # test a large kk

  ncdata = w.pdf("nc").generate(*w.var("x"),100)
  csdata = w.pdf("cs").generate(*w.var("x"),100)
  plot = w.var("x").frame()
  ncdata.plotOn(plot,MarkerColor(ROOT.kRed))
  csdata.plotOn(plot,MarkerColor(ROOT.kBlue))
  w.pdf("nc").plotOn(plot, ROOT.RooFit.LineColor(ROOT.kRed))
  w.pdf("ncc").plotOn(plot, ROOT.RooFit.LineColor(ROOT.kGreen))
  w.pdf("nccc").plotOn(plot, ROOT.RooFit.LineColor(ROOT.kYellow), ROOT.RooFit.LineStyle(ROOT.kDashed))
  w.pdf("cs").plotOn(plot, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(ROOT.kDotted))
  plot.Draw()


