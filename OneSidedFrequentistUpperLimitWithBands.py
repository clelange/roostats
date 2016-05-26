'''
OneSidedFrequentistUpperLimitWithBands

Author: Kyle Cranmer,
Contributions from Haichen Wang and Daniel Whiteson
date: Dec. 2010 - Feb. 2011
v1. Jan 28

This is a standard demo that can be used with any ROOT file
prepared in the standard way.  You specify:
 - name for input ROOT file
 - name of workspace inside ROOT file that holds model and data
 - name of ROOT.RooStats.ModelConfig( that specifies details for calculator tools
 - name of dataset

With default parameters the macro will attempt to run the
standard hist2workspace example and read the ROOT file
that it produces.

The first ~100 lines define a test statistic, the main macro starts.
You may want to control:
  double confidenceLevel=0.95
  nPointsToScan = 30
  nToyMC = 200

This uses a modified version of the profile likelihood ratio as
a test statistic for upper limits (eg. stat = 0 if muhat>mu).

Based on the observed data, defines a set of parameter points
to be tested based on the value of the parameter of interest
and the conditional MLE (eg. profiled) values of the nuisance parameters.

At each parameter point, pseudo-experiments are generated using self
fixed reference model and then the test statistic is evaluated.
Note, nuisance parameters are floating in the fits.  For each point,
the threshold that defines the 95% acceptance region is found.  ROOT.This
forms a "Confidence Belt".

After constructing the confidence belt, can find the confidence
interval for any particular dataset by finding the intersection
of the observed test statistic and the confidence belt.  First
this is done on the observed data to get an observed 1-sided upper limt.

Finally, expected limit and bands (from background-only) are
formed by generating background-only data and finding the upper limit.
This is done by hand for now, later be part of the ROOT.RooStats tools.

On a technical note, technique is NOT the Feldman-Cousins technique,
because that is a 2-sided interval BY DEFINITION.  However, the
Feldman-Cousins technique self is a Neyman-Construction.  For technical
reasons the easiest way to implement self right now is to use the
FeldmanCousins tool and then change the test statistic that it is using.

Building the confidence belt can be computationally expensive.  Once it is built,
one could save it to a file and use it in a separate step.

We can use PROOF to speed things along in parallel, however,
the test statistic has to be installed on the workers
so either turn off PROOF or include the modified test statistic
in your $ROOTSYS/roofit/roostats/inc directory,
add the additional line to the LinkDef.h file,
and recompile root.

Note, you have a boundary on the parameter of interest (eg. cross-section)
the threshold on the one-sided test statistic starts off very small because we
are only including downward fluctuations.  You can see the threshold in these printouts:

[#0] PROGRESS:Generation -- generated toys: 500 / 999
NeymanConstruction: Prog: 12/50 MC = 39 self stat = 0
 SigXsecOverSM=0.69 alpha_syst1=0.136515 alpha_syst3=0.425415 beta_syst2=1.08496 [-1e+30, 0.011215]  interval = 1

this tells you the values of the parameters being used to generate the pseudo-experiments
and the threshold in self case is 0.011215.  One would expect for 95% that the threshold
would be ~1.35 once the cross-section is far enough away from 0 that it is essentially
unaffected by the boundary.  As one reaches the last points in the scan, the
theshold starts to get artificially high.  ROOT.This is because the range of the parameter in
the fit is the same as the range in the scan.  In the future, should be independently
controled, they are not now.  As a result the ~50% of pseudo-experiments that have an
upward fluctuation end up muhat = muMax.  Because of self, upper range of the
parameter should be well above the expected upper limit... but not too high or one will
need a very large value of nPointsToScan to resolve the relevant region.  ROOT.This can be
improved, self is the first version of self script.

Important note: when the model includes external constraint terms, a Gaussian
constraint to a nuisance parameter centered around some nominal value there is
a subtlety.  ROOT.The asymptotic results are all based on the assumption that all the
measurements fluctuate... including the nominal values from auxiliary measurements.
If these do not fluctuate, corresponds to an "conditional ensemble".  ROOT.The
result is that the distribution of the test statistic can become very non-chi^2.
This results in thresholds that become very large. ROOT.This can be seen in the following
thought experiment.  Say the model is
   Pois(N | s + b) * G(b0|b,sigma)
where G(b0|b,sigma) is the external constraint and b0 is 100.  If N is also 100
then the profiled value of b given s is going to be some trade off betwen 100-s and b0.
If sigma is \sqrt(N), the profiled value of b is probably 100 - s/2   So for
s=60 we are going to have a profiled value of b~70.  Now when we generate pseudo-experiments
for s=60, b=70 we will have N~130 and the average shat will be 30, 60.  In practice,
this is only an issue for values of s that are very excluded.  For values of s near the 95%
limit self should not be a big effect.  ROOT.This can be avoided if the nominal values of the constraints also fluctuate, that requires that those parameters are ROOT.RooRealVars in the model.
This version does not deal with self issue, it will be addressed in a future version.
'''

#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "RooWorkspace.h"
#include "RooSimultaneous.h"
#include "RooAbsData.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"

#include "RooStats/RooStatsUtils.h"
#include "RooStats/ProfileLikelihoodTestStat.h"

import ROOT
using namespace ROOT.RooStats

useProof = False;  # flag to control whether to use Proof
nworkers = 0;   # number of workers (default use all available cores)

####################################/
# ROOT.The actual macro

void OneSidedFrequentistUpperLimitWithBands( infile = "",
                                             workspaceName = "combined",
                                             modelConfigName = "ModelConfig",
                                             dataName = "obsData")

#ifdef __CINT__
  cout << "DO NOT RUN WITH CINT: we are using a custom test statistic "
  cout << "which requires that self tutorial must be compiled "
  cout << "with ACLIC" << endl
  return
#endif

  double confidenceLevel=0.95
  nPointsToScan = 20
  nToyMC = 200

  ##############################/
  # First part is just to access a user-defined file
  # or create the standard example file if it doesn't exist
  ##############################
    filename = ""
   if not strcmp(infile, ""):      filename = "results/example_combined_GaussExample_model.root"
      fileExist = not gSystem.AccessPathName(filename); # note opposite return code
      # if file does not exists generate with histfactory
      if not fileExist:#ifdef _WIN32
         cout << "HistFactory file cannot be generated on Windows - exit" << endl
         return
#endif
         # Normally self would be run on the command line
         cout <<"will run standard hist2workspace example"<<endl
         gROOT.ProcessLine(".not  prepareHistFactory .")
         gROOT.ProcessLine(".not  hist2workspace config/example.xml")
         cout <<"\n\n---------------------"<<endl
         cout <<"Done creating example input"<<endl
         cout <<"---------------------\n\n"<<endl



   else:
      filename = infile

   # ROOT.Try to open the file
   ROOT.TFile *file = ROOT.TFile.Open(filename)

   # if input file was specified byt not found, quit
   if not file :      cout <<"StandardRooStatsDemoMacro: Input file " << filename << " is not found" << endl
      return



  ##############################/
  # Now get the data and workspace
  ##############################

  # get the workspace out of the file
  w = (ROOT.RooWorkspace*) file.Get(workspaceName)
  if not w:    cout <<"workspace not found" << endl
    return


  # get the modelConfig out of the file
  mc = (ModelConfig*) w.obj(modelConfigName)

  # get the modelConfig out of the file
  data = w.data(dataName)

  # make sure ingredients are found
  if not data or not mc:    w.Print()
    cout << "data or ROOT.RooStats.ModelConfig( was not found" <<endl
    return


  ##############################/
  # Now get the POI for convenience
  # you may want to adjust the range of your POI
  ##############################
  firstPOI = (ROOT.RooRealVar*) mc.GetParametersOfInterest().first()
  #  firstPOI.setMin(0)
  #  firstPOI.setMax(10)

  ######################/
  # create and use the ROOT.RooStats.FeldmanCousins( tool
  # to find and plot the 95% confidence interval
  # on the parameter of interest as specified
  # in the model config
  # REMEMBER, will change the test statistic
  # so self is NOT a Feldman-Cousins interval
  ROOT.RooStats.FeldmanCousins( fc(*data,*mc)
  fc.SetConfidenceLevel(confidenceLevel)
  #  fc.AdditionalNToysFactor(0.25); # degrade/improve sampling that defines confidence belt
  #  fc.UseAdaptiveSampling(True); # speed it up a bit, don't use for expectd limits
  fc.SetNBins(nPointsToScan); # set how many points per parameter of interest to scan
  fc.CreateConfBelt(True); # save the information in the belt for plotting

  ######################/
  # Feldman-Cousins is a unified limit by definition
  # but the tool takes care of a few things for us like which values
  # of the nuisance parameters should be used to generate toys.
  # so let's just change the test statistic and realize self is
  # no longer "Feldman-Cousins" but is a fully frequentist Neyman-Construction.
  #  ProfileLikelihoodTestStatModified onesided(*mc.GetPdf())
  #  fc.GetTestStatSampler().SetTestStatistic(&onesided)
  # ((ToyMCSampler*) fc.GetTestStatSampler()).SetGenerateBinned(True)
  toymcsampler = (ToyMCSampler*) fc.GetTestStatSampler()
  testStat = dynamic_cast<ProfileLikelihoodTestStat*>(toymcsampler.GetTestStatistic())
  testStat.SetOneSided(True)

  # Since self tool needs to throw toy MC the PDF needs to be
  # extended or the tool needs to know how many entries in a dataset
  # per pseudo experiment.
  # In the 'number counting form' where the entries in the dataset
  # are counts, not values of discriminating variables, the
  # datasets typically only have one entry and the PDF is not
  # extended.
  if not mc.GetPdf().canBeExtended():    if data.numEntries()==1:
      fc.FluctuateNumDataEntries(False)
    else:
      cout <<"Not sure what to do about self model" <<endl


  # We can use PROOF to speed things along in parallel
  # However, test statistic has to be installed on the workers
  # so either turn off PROOF or include the modified test statistic
  # in your $ROOTSYS/roofit/roostats/inc directory,
  # add the additional line to the LinkDef.h file,
  # and recompile root.
  if (useProof) { 
     ProofConfig pc(*w, nworkers, "", False)
     toymcsampler.SetProofConfig(&pc); # enable proof


  if mc.GetGlobalObservables():     cout << "will use global observables for unconditional ensemble"<<endl
     mc.GetGlobalObservables().Print()
     toymcsampler.SetGlobalObservables(*mc.GetGlobalObservables())



  # Now get the interval
  interval = fc.GetInterval()
  belt = fc.GetConfidenceBelt()

  # print out the iterval on the first Parameter of Interest
  cout << "\n95% interval on " <<firstPOI.GetName()<<" is : ["<<
    interval.LowerLimit(*firstPOI) << ", "<<
    interval.UpperLimit(*firstPOI) <<"] "<<endl

  # get observed UL and value of test statistic evaluated there
  ROOT.RooArgSet tmpPOI(*firstPOI)
  observedUL = interval.UpperLimit(*firstPOI)
  firstPOI.setVal(observedUL)
  obsTSatObsUL = fc.GetTestStatSampler().EvaluateTestStatistic(*data,tmpPOI)


  # Ask the calculator which points were scanned
  parameterScan = (ROOT.RooDataSet*) fc.GetPointsToScan()
  ROOT.RooArgSet* tmpPoint

  # make a histogram of parameter vs. threshold
  histOfThresholds = ROOT.TH1F("histOfThresholds", "",
                                    parameterScan.numEntries(),
                                    firstPOI.getMin(),
                                    firstPOI.getMax())
  histOfThresholds.GetXaxis().SetTitle(firstPOI.GetName())
  histOfThresholds.GetYaxis().SetTitle("Threshold")

  # loop through the points that were tested and ask confidence belt
  # what the upper/lower thresholds were.
  # For ROOT.RooStats.FeldmanCousins(, lower cut off is always 0
  for(Int_t i=0; i<parameterScan.numEntries(); ++i)    tmpPoint = (ROOT.RooArgSet*) parameterScan.get(i).clone("temp")
    #cout <<"get threshold"<<endl
    arMax = belt.GetAcceptanceRegionMax(*tmpPoint)
    poiVal = tmpPoint.getRealValue(firstPOI.GetName()) 
    histOfThresholds.Fill(poiVal,arMax)

  c1 = ROOT.TCanvas()
  c1.Divide(2)
  c1.cd(1)
  histOfThresholds.SetMinimum(0)
  histOfThresholds.Draw()
  c1.cd(2)

  ##############################/
  # Now we generate the expected bands and power-constriant
  ##############################

  # First: find parameter point for mu=0, conditional MLEs for nuisance parameters
  nll = mc.GetPdf().createNLL(*data)
  profile = nll.createProfile(*mc.GetParametersOfInterest())
  firstPOI.setVal(0.)
  profile.getVal(); # self will do fit and set nuisance parameters to profiled values
  poiAndNuisance = ROOT.RooArgSet()
  if mc.GetNuisanceParameters():
    poiAndNuisance.add(*mc.GetNuisanceParameters())
  poiAndNuisance.add(*mc.GetParametersOfInterest())
  w.saveSnapshot("paramsToGenerateData",*poiAndNuisance)
  paramsToGenerateData = (ROOT.RooArgSet*) poiAndNuisance.snapshot()
  cout << "\nWill use these parameter points to generate pseudo data for bkg only" << endl
  paramsToGenerateData.Print("v")


  ROOT.RooArgSet unconditionalObs
  unconditionalObs.add(*mc.GetObservables())
  unconditionalObs.add(*mc.GetGlobalObservables()); # comment self out for the original conditional ensemble

  double CLb=0
  double CLbinclusive=0

  # Now we generate background only and find distribution of upper limits
  histOfUL = ROOT.TH1F("histOfUL", "",100,0,firstPOI.getMax())
  histOfUL.GetXaxis().SetTitle("Upper Limit (background only)")
  histOfUL.GetYaxis().SetTitle("Entries")
  for(int imc=0; imc<nToyMC; ++imc)
    # set parameters back to values for generating pseudo data
    #    cout << "\n get current nuis, vals, again" << endl
    w.loadSnapshot("paramsToGenerateData")
    #    poiAndNuisance.Print("v")

    toyData = 0
    # now generate a toy dataset
    if not mc.GetPdf().canBeExtended():      if data.numEntries()==1:
         toyData = mc.GetPdf().generate(*mc.GetObservables(),1)
      else:
         cout <<"Not sure what to do about self model" <<endl
    } else:
       #      cout << "generating extended dataset"<<endl
       toyData = mc.GetPdf().generate(*mc.GetObservables(), ROOT.RooFit.Extended())


    # generate global observables
    # need to be careful for simpdf
    #    globalData = mc.GetPdf().generate(*mc.GetGlobalObservables(),1)

    simPdf = dynamic_cast<RooSimultaneous*>(mc.GetPdf())
    if not simPdf:      ROOT.RooDataSet *one = mc.GetPdf().generate(*mc.GetGlobalObservables(), 1)
       ROOT.RooArgSet *values = one.get()
      ROOT.RooArgSet *allVars = mc.GetPdf().getVariables()
      *allVars = *values
      delete allVars
      delete values
      delete one
    } else:

      #try fix for sim pdf
      iter = simPdf.indexCat().typeIterator() 
      tt = NULL
      while((tt=(ROOT.RooCatType*) iter.Next()))
         # Get pdf associated with state from simpdf
         pdftmp = simPdf.getPdf(tt.GetName()) 

         # Generate only global variables defined by the pdf associated with self state
         globtmp = pdftmp.getObservables(*mc.GetGlobalObservables()) 
         tmp = pdftmp.generate(*globtmp,1) 

         # ROOT.Transfer values to output placeholder
         *globtmp = *tmp.get(0) 

         # Cleanup
         delete globtmp 
         delete tmp 



    #    globalData.Print("v")
    #    unconditionalObs = *globalData.get()
    #    mc.GetGlobalObservables().Print("v")
    #    delete globalData
    #    cout << "data = " << endl
    #    toyData.get().Print("v")

    # get test stat at observed UL in observed data
    firstPOI.setVal(observedUL)
    toyTSatObsUL = fc.GetTestStatSampler().EvaluateTestStatistic(*toyData,tmpPOI)
    #    toyData.get().Print("v")
    #    cout <<"obsTSatObsUL " <<obsTSatObsUL << "toyTS " << toyTSatObsUL << endl
    if(obsTSatObsUL < toyTSatObsUL) # not sure about <= part yet
      CLb+= (1.)/nToyMC
    if(obsTSatObsUL <= toyTSatObsUL) # not sure about <= part yet
      CLbinclusive+= (1.)/nToyMC


    # loop over points in belt to find upper limit for self toy data
    thisUL = 0
    for(Int_t i=0; i<parameterScan.numEntries(); ++i)      tmpPoint = (ROOT.RooArgSet*) parameterScan.get(i).clone("temp")
      arMax = belt.GetAcceptanceRegionMax(*tmpPoint)
      firstPOI.setVal( tmpPoint.getRealValue(firstPOI.GetName()) )
      #   thisTS = profile.getVal()
      thisTS = fc.GetTestStatSampler().EvaluateTestStatistic(*toyData,tmpPOI)

      #   cout << "poi = " << firstPOI.getVal()
      # << " max is " << arMax << " profile = " << thisTS << endl
      #      cout << "thisTS = " << thisTS<<endl
      if thisTS<=arMax:         thisUL = firstPOI.getVal()
      } else:
         break





    '''
    # loop over points in belt to find upper limit for self toy data
    thisUL = 0
    for(Int_t i=0; i<histOfThresholds.GetNbinsX(); ++i)      tmpPoint = (ROOT.RooArgSet*) parameterScan.get(i).clone("temp")
      cout <<"----------------  "<<i<<endl
      tmpPoint.Print("v")
      cout << "from hist " << histOfThresholds.GetBinCenter(i+1) <<endl
      arMax = histOfThresholds.GetBinContent(i+1)
      # cout << " threhold Hist = aMax " << arMax<<endl
      # arMax2 = belt.GetAcceptanceRegionMax(*tmpPoint)
      # cout << "from arMax2 = "<< arMax2 << endl; # not the same due to ROOT.TH1F not ROOT.TH1D
      # cout << "scan - hist" << arMax2-arMax << endl
      firstPOI.setVal( histOfThresholds.GetBinCenter(i+1))
      #   thisTS = profile.getVal()
      thisTS = fc.GetTestStatSampler().EvaluateTestStatistic(*toyData,tmpPOI)

      #   cout << "poi = " << firstPOI.getVal()
      # << " max is " << arMax << " profile = " << thisTS << endl
      #      cout << "thisTS = " << thisTS<<endl

      # NOTE: need to add a small epsilon term for single precision vs. double precision
      if thisTS<=arMax + 1e-7:         thisUL = firstPOI.getVal()
      } else:
         break


    '''

    histOfUL.Fill(thisUL)

    # for few events, is often the same, UL is often the same
    #    cout << "thisUL = " << thisUL<<endl

    delete toyData

  histOfUL.Draw()
  c1.SaveAs("one-sided_upper_limit_output.pdf")

  # if you want to see a plot of the sampling distribution for a particular scan point:
  '''
  SamplingDistPlot sampPlot
  indexInScan = 0
  tmpPoint = (ROOT.RooArgSet*) parameterScan.get(indexInScan).clone("temp")
  firstPOI.setVal( tmpPoint.getRealValue(firstPOI.GetName()) )
  toymcsampler.SetParametersForTestStat(tmpPOI)
  samp = toymcsampler.GetSamplingDistribution(*tmpPoint)
  sampPlot.AddSamplingDistribution(samp)
  sampPlot.Draw()
   '''

  # Now find bands and power constraint
  bins = histOfUL.GetIntegral()
  cumulative = (TH1F*) histOfUL.Clone("cumulative")
  cumulative.SetContent(bins)
  double band2sigDown, band1sigDown, bandMedian, band1sigUp,band2sigUp
  for(int i=1; i<=cumulative.GetNbinsX(); ++i)    if bins[i]<RooStats.SignificanceToPValue(2):
      band2sigDown=cumulative.GetBinCenter(i)
    if bins[i]<RooStats.SignificanceToPValue(1):
      band1sigDown=cumulative.GetBinCenter(i)
    if bins[i]<0.5:
      bandMedian=cumulative.GetBinCenter(i)
    if bins[i]<RooStats.SignificanceToPValue(-1):
      band1sigUp=cumulative.GetBinCenter(i)
    if bins[i]<RooStats.SignificanceToPValue(-2):
      band2sigUp=cumulative.GetBinCenter(i)

  cout << "-2 sigma  band " << band2sigDown << endl
  cout << "-1 sigma  band " << band1sigDown << " [Power Constriant)]" << endl
  cout << "median of band " << bandMedian << endl
  cout << "+1 sigma  band " << band1sigUp << endl
  cout << "+2 sigma  band " << band2sigUp << endl

  # print out the iterval on the first Parameter of Interest
  cout << "\nobserved 95% upper-limit "<< interval.UpperLimit(*firstPOI) <<endl
  cout << "CLb strict [P(toy>obs|0)] for observed 95% upper-limit "<< CLb <<endl
  cout << "CLb inclusive [P(toy>=obs|0)] for observed 95% upper-limit "<< CLbinclusive <<endl

  delete profile
  delete nll


