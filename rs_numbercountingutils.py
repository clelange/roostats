# /
#
# 'Number Counting Utils' ROOT.RooStats tutorial
# author: Kyle Cranmer
# date June. 2009
#
# ROOT.This tutorial shows an example of the ROOT.RooStats standalone
# utilities that calculate the p-value or Z value (eg. significance in
# 1-sided Gaussian standard deviations) for a number counting experiment.
# ROOT.This is a hypothesis test between background only and signal-plus-background.
# ROOT.The background estimate has uncertainty derived from an auxiliary or sideband
# measurement.
#
# Documentation for these utilities can be found here:
# http:#root.cern.ch/root/html/RooStats__NumberCountingUtils.html
#
#
# ROOT.This problem is often called a proto-type problem for high energy physics.
# In some references it is referred to as the on/off problem.
#
# ROOT.The problem is treated in a fully frequentist fashion by
# interpreting the relative background uncertainty as
# being due to an auxiliary or sideband observation
# that is also Poisson distributed with only background.
# Finally, considers the test as a ratio of Poisson means
# where an interval is well known based on the conditioning on the total
# number of events and the binomial distribution.
# For more on self, see
#  http:#arxiv.org/abs/0905.3831
#  http:#arxiv.org/abs/physics/physics/0702156
#  http:#arxiv.org/abs/physics/0511028
#
# /


import ROOT


def rs_numbercountingutils():

    # From the root prompt, can see the full list of functions by using
    # tab-completion

    # root [0] ROOT.RooStats.NumberCountingUtils.  <tab>
    # BinomialExpZ
    # BinomialWithTauExpZ
    # BinomialObsZ
    # BinomialWithTauObsZ
    # BinomialExpP
    # BinomialWithTauExpP
    # BinomialObsP
    # BinomialWithTauObsP

    # For each of the utilities you can inspect the arguments by tab completion

    # root [1] NumberCountingUtils.BinomialExpZ( <tab>
    # Double_t BinomialExpZ(Double_t sExp, bExp, fractionalBUncertainty)

    # /
    # Here we see common usages where the experimenter
    # has a relative background uncertainty, without
    # explicit reference to the auxiliary or sideband
    # measurement

    # /
    # Expected p-values and significance with background uncertainty
    ##########################
    sExpected = 50
    bExpected = 100
    relativeBkgUncert = 0.1

    pExp = ROOT.RooStats.NumberCountingUtils.BinomialExpP(
        sExpected, bExpected, relativeBkgUncert)
    zExp = ROOT.RooStats.NumberCountingUtils.BinomialExpZ(
        sExpected, bExpected, relativeBkgUncert)
    print "expected p-value =", pExp, "  Z value (Gaussian sigma) = ", zExp

    # /
    # Expected p-values and significance with background uncertainty
    ##########################
    observed = 150
    pObs = ROOT.RooStats.NumberCountingUtils.BinomialObsP(
        observed, bExpected, relativeBkgUncert)
    zObs = ROOT.RooStats.NumberCountingUtils.BinomialObsZ(
        observed, bExpected, relativeBkgUncert)
    print "observed p-value =", pObs, "  Z value (Gaussian sigma) = ", zObs

    # /
    # Here we see usages where the experimenter has knowledge
    # about the properties of the auxiliary or sideband
    # measurement.  In particular, ratio tau of background
    # in the auxiliary measurement to the main measurement.
    # Large values of tau mean small background uncertainty
    # because the sideband is very constraining.

    # Usage:
    # root [0] ROOT.RooStats.NumberCountingUtils.BinomialWithTauExpP(
    # Double_t BinomialWithTauExpP(Double_t sExp, bExp, tau)

    # /
    # Expected p-values and significance with background uncertainty
    ##########################
    tau = 1

    pExpWithTau = ROOT.RooStats.NumberCountingUtils.BinomialWithTauExpP(
        sExpected, bExpected, tau)
    zExpWithTau = ROOT.RooStats.NumberCountingUtils.BinomialWithTauExpZ(
        sExpected, bExpected, tau)
    print "expected p-value =", pExpWithTau, "  Z value (Gaussian sigma) = ", zExpWithTau

    # /
    # Expected p-values and significance with background uncertainty
    ##########################
    pObsWithTau = ROOT.RooStats.NumberCountingUtils.BinomialWithTauObsP(
        observed, bExpected, tau)
    zObsWithTau = ROOT.RooStats.NumberCountingUtils.BinomialWithTauObsZ(
        observed, bExpected, tau)
    print "observed p-value =", pObsWithTau, "  Z value (Gaussian sigma) = ", zObsWithTau


if __name__ == "__main__":
    rs_numbercountingutils()
