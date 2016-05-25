# This example is a generalization of the on/off problem.

'''
FourBin Instructional Tutorial:
 authors:
 Kyle Cranmer <cranmer@cern.ch>
 ROOT.Tanja Rommerskirchen <tanja.rommerskirchen@cern.ch>

 date: June 1, 2010

This example is a generalization of the on/off problem.
It's a common setup for SUSY searches.  Imagine that one has two
variables "x" and "y" (eg. missing ET and SumET), figure.
The signal region has high values of both of these variables (top right).
One can see low values of "x" or "y" acting as side-bands.  If we
just used "y" as a sideband, would have the on/off problem.
 - In the signal region we observe non events and expect s+b events.
 - In the region with low values of "y" (bottom right)
   we observe noff events and expect tau*b events.
Note the significance of tau.  In the background only case:
   tau ~ <expectation off> / <expectation on>
If tau is known, model is sufficient, often tau is not known exactly.
So one can use low values of "x" as an additional constraint for tau.
Note that self technique critically depends on the notion that the
joint distribution for "x" and "y" can be factorized.
Generally, regions have many events, it the ratio can be
measured very precisely there.  So we extend the model to describe the
left two boxes... denoted with "bar".
  - In the upper left we observe nonbar events and expect bbar events
  - In the bottom left we observe noffbar events and expect tau bbar events
Note again we have:
   tau ~ <expecation off bar> / <expectation on bar>
One can further expand the model to account for the systematic associated
to assuming the distribution of "x" and "y" factorizes (eg. that
tau is the same for off/on and offbar/onbar). ROOT.This can be done in several
ways, here we introduce an additional parameter rho, so that
one set of models will use tau and the other tau*rho. ROOT.The choice is arbitary,
but it has consequences on the numerical stability of the algorithms.
The "bar" measurements typically have more events (& smaller relative errors).
If we choose <expectation noffbar> = tau * rho * <expectation noonbar>, the
product tau*rho will be known very precisely (~1/sqrt(bbar)) and the contour
in those parameters will be narrow and have a non-trivial tau~1/rho shape.
However, we choose to put rho on the non/noff measurements (where the
product will have an error ~1/sqrt(b)), contours will be more ameanable
to numerical techniques.  ROOT.Thus, we choose to define
   tau := <expecation off bar> / (<expectation on bar>)
   rho := <expecation off> / (<expectation on> * tau)

^ y
|
|---------------------------+
|               |           |
|     nonbar    |    non    |
|      bbar     |    s+b    |
|               |           |
|---------------+-----------|
|               |           |
|    noffbar    |    noff   |
|    tau bbar   | tau b rho |
|               |           |
+----------------------------. x


Left in self way, problem is under-constrained.  However, may
have some auxiliary measurement (usually based on Monte Carlo) to
constrain rho.  Let us call self auxiliary measurement that gives
the nominal value of rho "rhonom".  ROOT.Thus, is a 'constraint' term in
the full model: P(rhonom | rho).  In self case, consider a Gaussian
constraint with standard deviation sigma.

In the example, initial values of the parameters are
  - s    = 40
  - b    = 100
  - tau  = 5
  - bbar = 1000
  - rho  = 1
  (sigma rho = 20%)
and in the toy dataset:
   - non = 139
   - noff = 528
   - nonbar = 993
   - noffbar = 4906
   - rhonom = 1.27824

Note, covariance matrix of the parameters has large off-diagonal terms.
Clearly s, are anti-correlated.  Similary, noffbar >> nonbar, would
expect bbar, to be anti-correlated.

This can be seen below.
            GLOBAL      b    bbar   rho      s     tau
        b  0.96820   1.000  0.191 -0.942 -0.762 -0.209
     bbar  0.91191   0.191  1.000  0.000 -0.146 -0.912
      rho  0.96348  -0.942  0.000  1.000  0.718 -0.000
        s  0.76250  -0.762 -0.146  0.718  1.000  0.160
      tau  0.92084  -0.209 -0.912 -0.000  0.160  1.000

Similarly, tau*rho appears as a product, expect rho,tau
to be anti-correlated. When the error on rho is significantly
larger than 1/sqrt(bbar), is essentially known and the
correlation is minimal (tau mainly cares about bbar, rho about b,s).
In the alternate parametrizaton (bbar* tau * rho) the correlation coefficient
for rho, is large (and negative).


The code below uses best-practices for ROOT.RooFit & ROOT.RooStats as of
June 2010.  It proceeds as follows:
 - create a workspace to hold the model
 - use workspace factory to quickly create the terms of the model
 - use workspace factory to define total model (a prod pdf)
 - create a ROOT.RooStats ROOT.RooStats.ModelConfig( to specify observables, of interest
 - add to the ROOT.RooStats.ModelConfig( a prior on the parameters for Bayesian techniques
   - note, pdf it is factorized for parameters of interest & nuisance params
 - visualize the model
 - write the workspace to a file
 - use several of ROOT.RooStats IntervalCalculators & compare results
'''


import ROOT


def FourBinInstructional(doBayesian=False, doFeldmanCousins=False, doMCMC=False):
    # let's time self challenging example
    t = ROOT.TStopwatch()
    t.Start()

    # set ROOT.RooFit random seed for reproducible results
    _rnd = ROOT.RooRandom.randomGenerator().SetSeed(4357)

    # make model
    wspace = ROOT.RooWorkspace("wspace")
    wspace.factory("Poisson::on(non[0,1000], sum::splusb(s[40,0,100],b[100,0,300]))")
    wspace.factory("Poisson::off(noff[0,5000], prod::taub(b,tau[5,3,7],rho[1,0,2]))")
    wspace.factory("Poisson::onbar(nonbar[0,10000], bbar[1000,500,2000])")
    wspace.factory("Poisson::offbar(noffbar[0,1000000], prod::lambdaoffbar(bbar, tau))")
    wspace.factory("Gaussian::mcCons(rhonom[1.,0,2], rho, sigma[.2])")
    wspace.factory("PROD::model(on,off,onbar,offbar,mcCons)")
    wspace.defineSet("obs","non,noff,nonbar,noffbar,rhonom")

    wspace.factory("Uniform::prior_poi({s})")
    wspace.factory("Uniform::prior_nuis({b,bbar,tau, rho})")
    wspace.factory("PROD::prior(prior_poi,prior_nuis)")

    # /
    # Control some interesting variations
    # define parameers of interest
    # for 1-d plots
    wspace.defineSet("poi", "s")
    wspace.defineSet("nuis", "b,tau,rho,bbar")
    # for 2-d plots to inspect correlations:
    #  wspace.defineSet("poi", "s,rho")

    # test simpler cases where parameters are known.
    #  wspace.var("tau").setConstant()
    #  wspace.var("rho").setConstant()
    #  wspace.var("b").setConstant()
    #  wspace.var("bbar").setConstant()

    # inspect workspace
    #  wspace.Print()

    ##############################
    # Generate toy data
    # generate toy data assuming current value of the parameters
    # import into workspace.
    # add Verbose() to see how it's being generated
    data = wspace.pdf("model").generate(wspace.set("obs"), 1)
    #  data.Print("v")
    getattr(wspace, 'import')(data)

    # /
    # Now the statistical tests
    # model config
    modelConfig = ROOT.RooStats.ModelConfig(("FourBins"))
    modelConfig.SetWorkspace(wspace)
    modelConfig.SetPdf(wspace.pdf("model"))
    modelConfig.SetPriorPdf(wspace.pdf("prior"))
    modelConfig.SetParametersOfInterest(wspace.set("poi"))
    modelConfig.SetNuisanceParameters(wspace.set("nuis"))
    getattr(wspace, 'import')(modelConfig)
    wspace.writeToFile("FourBin.root")

    #########################
    # If you want to see the covariance matrix uncomment
    #  wspace.pdf("model").fitTo(data)

    # use ProfileLikelihood
    plc = ROOT.RooStats.ProfileLikelihoodCalculator(data, modelConfig)
    plc.SetConfidenceLevel(0.95)
    plInt = plc.GetInterval()
    msglevel = ROOT.RooMsgService.instance().globalKillBelow()
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    plInt.LowerLimit(wspace.var("s"))  # get ugly print out of the way. Fix.
    ROOT.RooMsgService.instance().setGlobalKillBelow(msglevel)

    # use FeldmaCousins (takes ~20 min)
    fc = ROOT.RooStats.FeldmanCousins(data, modelConfig)
    fc.SetConfidenceLevel(0.95)
    # number counting: dataset always has 1 entry with N events observed
    fc.FluctuateNumDataEntries(False)
    fc.UseAdaptiveSampling(True)
    fc.SetNBins(40)
    fcInt = 0
    if(doFeldmanCousins):  # takes 7 minutes
        fcInt = fc.GetInterval()  # fix cast

    # use BayesianCalculator (only 1-d parameter of interest, for self problem)
    bc = ROOT.RooStats.BayesianCalculator(data, modelConfig)
    bc.SetConfidenceLevel(0.95)
    bInt = 0
    if doBayesian and wspace.set("poi").getSize() == 1:
        bInt = bc.GetInterval()
    else:
        print "Bayesian Calc. only supports on parameter of interest"

    # use ROOT.RooStats.MCMCCalculator(  (takes about 1 min)
    # Want an efficient proposal function, derive it from covariance
    # matrix of fit
    fit = wspace.pdf("model").fitTo(data, ROOT.RooFit.Save())
    ph = ROOT.RooStats.ProposalHelper()
    ph.SetVariables(fit.floatParsFinal())
    ph.SetCovMatrix(fit.covarianceMatrix())
    # auto-create mean vars and add mappings
    ph.SetUpdateProposalParameters(ROOT.kTRUE)
    ph.SetCacheSize(100)
    pf = ph.GetProposalFunction()

    mc = ROOT.RooStats.MCMCCalculator(data, modelConfig)
    mc.SetConfidenceLevel(0.95)
    mc.SetProposalFunction(pf)
    mc.SetNumBurnInSteps(500)  # first N steps to be ignored as burn-in
    mc.SetNumIters(50000)
    mc.SetLeftSideTailFraction(0.5)  # make a central interval
    mcInt = 0
    if doMCMC:
        mcInt = mc.GetInterval()

    ###################
    # Make some  plots
    c1 = ROOT.gROOT.Get("c1")
    if not c1:
        c1 = ROOT.TCanvas("c1")

    if doBayesian and doMCMC:
        c1.Divide(3)
        c1.cd(1)

    elif doBayesian or doMCMC:
        c1.Divide(2)
        c1.cd(1)

    lrplot = ROOT.RooStats.LikelihoodIntervalPlot(plInt)
    lrplot.Draw()

    if doBayesian and wspace.set("poi").getSize() == 1:
        c1.cd(2)
        # the plot takes a long time and print lots of error
        # using a scan it is better
        bc.SetScanOfPosterior(20)
        plot = bc.GetPosteriorPlot()
        plot.Draw()

    if doMCMC:
        if doBayesian and wspace.set("poi").getSize() == 1:
            c1.cd(3)
        else:
            c1.cd(2)
        mcPlot = ROOT.RooStats.MCMCIntervalPlot(mcInt)
        mcPlot.Draw()

    ##################
    # querry intervals
    print "Profile Likelihood interval s = [", plInt.LowerLimit(wspace.var("s")), ", ", plInt.UpperLimit(wspace.var("s")), "]"
    # Profile Likelihood interval s = [12.1902, 88.6871]

    if doBayesian and wspace.set("poi").getSize() == 1:
        print "Bayesian interval s = [", bInt.LowerLimit(), ", ",  bInt.UpperLimit(), "]"

    if doFeldmanCousins:
        print "Feldman Cousins interval s = [", fcInt.LowerLimit(wspace.var("s")), ", ", fcInt.UpperLimit(wspace.var("s")), "]"
    # Feldman Cousins interval s = [18.75 +/- 2.45, 83.75 +/- 2.45]

    if doMCMC:
        print "MCMC interval s = [", mcInt.LowerLimit(wspace.var("s")), ", ", mcInt.UpperLimit(wspace.var("s")), "]"
        # MCMC interval s = [15.7628, 84.7266]

    t.Print()
    c1.SaveAs("FourBinInstructional.png")


if __name__ == "__main__":
    FourBinInstructional(False, False, False)
