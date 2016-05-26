# comparison of MCMC and PLC in a multi-variate gaussian problem

#include "RooGlobalFunc.h"
#include <stdlib.h>
#include "TMatrixDSym.h"
#include "RooMultiVarGaussian.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "RooAbsReal.h"
#include "RooFitResult.h"
#include "TStopwatch.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MetropolisHastings.h"
#include "RooStats/MarkovChain.h"
#include "RooStats/ConfInterval.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ProposalFunction.h"
#include "RooStats/PdfProposal.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/LikelihoodInterval.h"

using namespace std
import ROOT
using namespace ROOT.RooStats


def MultivariateGaussianTest(self, dim = 4, nPOI = 2):
  '''
    Authors: Kevin Belasco and Kyle Cranmer.

    ROOT.This tutorial produces an N-dimensional multivariate Gaussian
    with a non-trivial covariance matrix.  By default N=4 (called "dim").

    A subset of these are considered parameters of interest.
    ROOT.This problem is tractable analytically.

    We use self mainly as a test of Markov Chain Monte Carlo
    and we compare the result to the profile likelihood ratio.

    We use the proposal helper to create a customized
    proposal function for self problem.

    For N=4 and 2 parameters of interest it takes about 10-20 seconds
    and the acceptance rate is 37%
   '''

  # let's time self challenging example
  ROOT.TStopwatch t
  t.Start()

   ROOT.RooArgList xVec
   ROOT.RooArgList muVec
   ROOT.RooArgSet poi

   # make the observable and means
   Int_t i,j
   ROOT.RooRealVar* x
   ROOT.RooRealVar* mu_x
   for (i = 0; i < dim; i++)      name = Form("x%d", i)
      x = ROOT.RooRealVar(name, name, 0, -3,3)
      xVec.add(*x)

      mu_name = Form("mu_x%d",i)
      mu_x = ROOT.RooRealVar(mu_name, mu_name, 0, -2,2)
      muVec.add(*mu_x)


   # put them into the list of parameters of interest
   for (i = 0; i < nPOI; i++)      poi.add(*muVec.at(i))


   # make a covariance matrix that is all 1's
   ROOT.TMatrixDSym cov(dim)
   for (i = 0; i < dim; i++)      for (j = 0; j < dim; j++)         if (i == j) cov(i,j) = 3.
         else        cov(i,j) = 1.0



   # now make the multivariate Gaussian
   ROOT.RooMultiVarGaussian mvg("mvg", "mvg", xVec, muVec, cov)

   #####################/
   # make a toy dataset
   data = mvg.generate(ROOT.RooArgSet(x)Vec, 100)

   #######################/
   # now create the model config for self problem
   w = ROOT.RooWorkspace("MVG")
   ROOT.RooStats.ModelConfig( modelConfig(w)
   modelConfig.SetPdf(mvg)
   modelConfig.SetParametersOfInterest(poi)

   #####################/
   # Setup calculators

   # MCMC
   # we want to setup an efficient proposal function
   # using the covariance matrix from a fit to the data
   fit = mvg.fitTo(*data, ROOT.RooFit.Save(True))
   ROOT.RooStats.ProposalHelper( ph
   ph.SetVariables((ROOT.RooArgSet&)fit.floatParsFinal())
   ph.SetCovMatrix(fit.covarianceMatrix())
   ph.SetUpdateProposalParameters(True)
   ph.SetCacheSize(100)
   pdfProp = ph.GetProposalFunction()

   # now create the calculator
   ROOT.RooStats.MCMCCalculator( mc(*data, modelConfig)
   mc.SetConfidenceLevel(0.95)
   mc.SetNumBurnInSteps(100)
   mc.SetNumIters(10000)
   mc.SetNumBins(50)
   mc.SetProposalFunction(*pdfProp)

   mcInt = mc.GetInterval()
   poiList = mcInt.GetAxes()

   # now setup the profile likelihood calculator
   ROOT.RooStats.ProfileLikelihoodCalculator( plc(*data, modelConfig)
   plc.SetConfidenceLevel(0.95)
   plInt = (LikelihoodInterval*)plc.GetInterval()


   #####################/
   # make some plots
   ROOT.RooStats.MCMCIntervalPlot( mcPlot(*mcInt)

   c1 = ROOT.TCanvas()
   mcPlot.SetLineColor(ROOT.kGreen)
   mcPlot.SetLineWidth(2)
   mcPlot.Draw()

   ROOT.RooStats.LikelihoodIntervalPlot( plPlot(plInt)
   plPlot.Draw("same")

   if poiList.getSize() == 1:      p = (ROOT.RooRealVar*)poiList.at(0)
      ll = mcInt.LowerLimit(*p)
      ul = mcInt.UpperLimit(*p)
      cout << "MCMC interval: [" << ll << ", " << ul << "]" << endl


   if poiList.getSize() == 2:      p0 = (ROOT.RooRealVar*)poiList.at(0)
      p1 = (ROOT.RooRealVar*)poiList.at(1)
      scatter = ROOT.TCanvas()
      ll = mcInt.LowerLimit(*p0)
      ul = mcInt.UpperLimit(*p0)
      cout << "MCMC interval on p0: [" << ll << ", " << ul << "]" << endl
      ll = mcInt.LowerLimit(*p0)
      ul = mcInt.UpperLimit(*p0)
      cout << "MCMC interval on p1: [" << ll << ", " << ul << "]" << endl

      #MCMC interval on p0: [-0.2, 0.6]
      #MCMC interval on p1: [-0.2, 0.6]

      mcPlot.DrawChainScatter(*p0, *p1)
      scatter.Update()


  t.Print()



