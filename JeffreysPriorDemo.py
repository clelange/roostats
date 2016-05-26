# tutorial demonstrating and validates the ROOT.RooJeffreysPrior class

'''
JeffreysPriorDemo.C

author Kyle Cranmer
date   Dec. 2010

This tutorial demonstraites and validates the ROOT.RooJeffreysPrior class

Jeffreys's prior is an 'objective prior' based on formal rules.
It is calculated from the Fisher information matrix.

Read more:
http:#en.wikipedia.org/wiki/Jeffreys_prior

The analytic form is not known for most PDFs, it is for
simple cases like the Poisson mean, mean, sigma.

This class uses numerical tricks to calculate the Fisher Information Matrix
efficiently.  In particular, takes advantage of a property of the
'Asimov data' as described in
Asymptotic formulae for likelihood-based tests of physics
Glen Cowan, Cranmer, Gross, Vitells
http:#arxiv.org/abs/arXiv:1007.1727

This Demo has four parts:
TestJeffreysPriorDemo -- validates Poisson mean case 1/sqrt(mu)
TestJeffreysGaussMean -- validates Gaussian mean case
TestJeffreysGaussSigma -- validates Gaussian sigma case 1/sigma
TestJeffreysGaussMeanAndSigma -- demonstraites 2-d example

'''

#include "RooJeffreysPrior.h"

#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "TMatrixDSym.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooNumIntConfig.h"
#include "TH1F.h"

import ROOT

def JeffreysPriorDemo(self):  ROOT.RooWorkspace w("w")
  w.factory("Uniform.u(x[0,1])")
  w.factory("mu[100,1,200]")
  w.factory("ExtendPdf.p(u,mu)")

  #  w.factory("Poisson.pois(n[0,inf],mu)")

  asimov = w.pdf("p").generateBinned(*w.var("x"),ExpectedData())
  #  asimov2 = w.pdf("pois").generateBinned(*w.var("n"),ExpectedData())

  res = w.pdf("p").fitTo(*asimov, ROOT.RooFit.Save(),SumW2Error(kTRUE))

  asimov.Print()
  res.Print()
  cov = res.covarianceMatrix()
  cout << "variance = " << (cov.Determinant()) << endl
  cout << "stdev = " << sqrt(cov.Determinant()) << endl
  cov.Invert()
  cout << "jeffreys = " << sqrt(cov.Determinant()) << endl

  w.defineSet("poi", "mu")
  w.defineSet("obs", "x")
  #  w.defineSet("obs2", "n")

  ROOT.RooJeffreysPrior pi("jeffreys", "jeffreys",*w.pdf("p"),*w.set("poi"),*w.set("obs"))
  #  pi.specialIntegratorConfig(kTRUE).method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")  
  #  pi.specialIntegratorConfig(kTRUE).getConfigSection("RooIntegrator1D").setRealValue("maxSteps",10)

  #  JeffreysPrior pi2("jeffreys2", "jeffreys",*w.pdf("pois"),*w.set("poi"),*w.set("obs2"))

  #  return
  test = ROOT.RooGenericPdf("test", "test", "1./sqrt(mu)",*w.set("poi"))

  c1 = ROOT.TCanvas
  plot = w.var("mu").frame()
  #  pi.plotOn(plot, Normalization(1, ROOT.RooAbsReal.Raw),Precision(.1))
  pi.plotOn(plot)
  #  pi2.plotOn(plot, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDotted))
  test.plotOn(plot, ROOT.RooFit.LineColor(ROOT.kRed))
  plot.Draw()




#_________________________________________________
def ROOT.TestJeffreysGaussMean(self):  ROOT.RooWorkspace w("w")
  w.factory("Gaussian.g(x[0,-20,20],mu[0,-5,5],sigma[1,0,10])")
  w.factory("n[10,.1,200]")
  w.factory("ExtendPdf.p(g,n)")
  w.var("sigma").setConstant()
  w.var("n").setConstant()

  asimov = w.pdf("p").generateBinned(*w.var("x"),ExpectedData())

  res = w.pdf("p").fitTo(*asimov, ROOT.RooFit.Save(),SumW2Error(kTRUE))

  asimov.Print()
  res.Print()
  cov = res.covarianceMatrix()
  cout << "variance = " << (cov.Determinant()) << endl
  cout << "stdev = " << sqrt(cov.Determinant()) << endl
  cov.Invert()
  cout << "jeffreys = " << sqrt(cov.Determinant()) << endl

  #  w.defineSet("poi", "mu,sigma")
  w.defineSet("poi", "mu")
  w.defineSet("obs", "x")

  ROOT.RooJeffreysPrior pi("jeffreys", "jeffreys",*w.pdf("p"),*w.set("poi"),*w.set("obs"))
  #  pi.specialIntegratorConfig(kTRUE).method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")  
  #  pi.specialIntegratorConfig(kTRUE).getConfigSection("RooIntegrator1D").setRealValue("maxSteps",3)

   temp = w.set("poi")
  pi.getParameters(*temp).Print()

  #  return
  test = ROOT.RooGenericPdf("test", "test", "1",*w.set("poi"))

  c1 = ROOT.TCanvas
  plot = w.var("mu").frame()
  pi.plotOn(plot)
  test.plotOn(plot, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDotted))
  plot.Draw()




#_________________________________________________
def ROOT.TestJeffreysGaussSigma(self):  # self one is VERY sensitive
  # if the Gaussian is narrow ~ range(x)/nbins(x) then the peak isn't resolved
  #   and you get really bizzare shapes
  # if the Gaussian is too wide range(x) ~ sigma then PDF gets renormalized
  #   and the PDF falls off too fast at high sigma
  ROOT.RooWorkspace w("w")
  w.factory("Gaussian.g(x[0,-20,20],mu[0,-5,5],sigma[1,1,5])")
  w.factory("n[100,.1,2000]")
  w.factory("ExtendPdf.p(g,n)")
  #  w.var("sigma").setConstant()
  w.var("mu").setConstant()
  w.var("n").setConstant()
  w.var("x").setBins(301)

  asimov = w.pdf("p").generateBinned(*w.var("x"),ExpectedData())

  res = w.pdf("p").fitTo(*asimov, ROOT.RooFit.Save(),SumW2Error(kTRUE))

  asimov.Print()
  res.Print()
  cov = res.covarianceMatrix()
  cout << "variance = " << (cov.Determinant()) << endl
  cout << "stdev = " << sqrt(cov.Determinant()) << endl
  cov.Invert()
  cout << "jeffreys = " << sqrt(cov.Determinant()) << endl


  #  w.defineSet("poi", "mu,sigma")
  #w.defineSet("poi", "mu,sigma,n")
  w.defineSet("poi", "sigma")
  w.defineSet("obs", "x")

  ROOT.RooJeffreysPrior pi("jeffreys", "jeffreys",*w.pdf("p"),*w.set("poi"),*w.set("obs"))
  #  pi.specialIntegratorConfig(kTRUE).method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")  
  pi.specialIntegratorConfig(kTRUE).getConfigSection("RooIntegrator1D").setRealValue("maxSteps",3)

   temp = w.set("poi")
  pi.getParameters(*temp).Print()
  #  return

  #  return
  test = ROOT.RooGenericPdf("test", "test", "sqrt(2.)/sigma",*w.set("poi"))

  c1 = ROOT.TCanvas
  plot = w.var("sigma").frame()
  pi.plotOn(plot)
  test.plotOn(plot, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDotted))
  plot.Draw()





#_________________________________________________
def ROOT.TestJeffreysGaussMeanAndSigma(self):  # self one is VERY sensitive
  # if the Gaussian is narrow ~ range(x)/nbins(x) then the peak isn't resolved
  #   and you get really bizzare shapes
  # if the Gaussian is too wide range(x) ~ sigma then PDF gets renormalized
  #   and the PDF falls off too fast at high sigma
  ROOT.RooWorkspace w("w")
  w.factory("Gaussian.g(x[0,-20,20],mu[0,-5,5],sigma[1,1,5])")
  w.factory("n[100,.1,2000]")
  w.factory("ExtendPdf.p(g,n)")
  #  w.var("sigma").setConstant()
  #  w.var("mu").setConstant()
  w.var("n").setConstant()
  w.var("x").setBins(301)

  asimov = w.pdf("p").generateBinned(*w.var("x"),ExpectedData())

  res = w.pdf("p").fitTo(*asimov, ROOT.RooFit.Save(),SumW2Error(kTRUE))

  asimov.Print()
  res.Print()
  cov = res.covarianceMatrix()
  cout << "variance = " << (cov.Determinant()) << endl
  cout << "stdev = " << sqrt(cov.Determinant()) << endl
  cov.Invert()
  cout << "jeffreys = " << sqrt(cov.Determinant()) << endl


  w.defineSet("poi", "mu,sigma")
  #w.defineSet("poi", "mu,sigma,n")
  #  w.defineSet("poi", "sigma")
  w.defineSet("obs", "x")

  ROOT.RooJeffreysPrior pi("jeffreys", "jeffreys",*w.pdf("p"),*w.set("poi"),*w.set("obs"))
  #  pi.specialIntegratorConfig(kTRUE).method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D")  
  pi.specialIntegratorConfig(kTRUE).getConfigSection("RooIntegrator1D").setRealValue("maxSteps",3)

   temp = w.set("poi")
  pi.getParameters(*temp).Print()
  #  return

  c1 = ROOT.TCanvas
  Jeff2d = pi.createHistogram("2dJeffreys",*w.var("mu"), ROOT.RooFit.Binning(10), ROOT.RooFit.YVar(*w.var("sigma"), ROOT.RooFit.Binning(10)))
  Jeff2d.Draw("surf")




