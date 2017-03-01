################################################################
## In principle all you have to setup is defined in this file ##
################################################################
from configManager import configMgr
from ROOT import kBlack,kWhite,kGray,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,kGreen,kSpring,kYellow,kOrange
from configWriter import fitConfig,Measurement,Channel,Sample
from systematic import Systematic
from math import sqrt
import time
from logger import Logger

import os

######################
# Results data
######################

results = {}
## key=region name => [exp, expUnc, obs]

# RJigsaw Moriond 2017 results on v115
results["SRJigsawSRS1a"] = [274, 35, 300] 
results["SRJigsawSRS1b"] = [190, 22, 221] 
results["SRJigsawSRS2a"] = [73, 10, 85] 
results["SRJigsawSRS2b"] = [56, 7, 67] 
results["SRJigsawSRS3a"] = [40.4, 3.5, 52] 
results["SRJigsawSRS3b"] = [26.7, 2.0, 32] 

results["SRJigsawSRG1a"] = [46.65, 4.88, 38] 
results["SRJigsawSRG1b"] = [15.0, 1.8, 13] 
results["SRJigsawSRG2a"] = [19.9, 2.2, 29] 
results["SRJigsawSRG2b"] = [6.06, 1.0, 10] 
results["SRJigsawSRG3a"] = [4.2, 0.8, 8] 
results["SRJigsawSRG3b"] = [1.3, 0.4, 4] 

results["SRJigsawSRC1"] = [36.8, 3.66, 31] 
results["SRJigsawSRC2"] = [30.25, 2.96, 25] 
results["SRJigsawSRC3"] = [20.38, 2.54, 12] 
results["SRJigsawSRC4"] = [25.95, 3.65, 21] 
results["SRJigsawSRC5"] = [7.51, 1.56, 8] 


##########################

log = Logger('SimpleUL')
try:
    pickedSRs
except NameError:
    log.fatal("No region specified!")    

if len(pickedSRs) == 0:
    log.fatal("No region specified!")    

for SR in pickedSRs:

    if SR not in results:
        log.warning("SR %s not found in results dict!")
        time.sleep(3)

    ##########################

    # Set observed and expected number of events in counting experiment
    ndata     =  float(results[SR][2]) # Number of events observed in data
    if ndata == 0.0:
        ndata = 0.001

    nbkg      =  float(results[SR][0]) # Number of predicted bkg events
    nbkgErr   =  float(results[SR][1]) # (Absolute) Statistical error on bkg estimate

    lumiError = 0.032 	# Relative luminosity uncertainty

    ucb = Systematic("ucb", configMgr.weights, 1 + nbkgErr/nbkg, 1 - nbkgErr/nbkg, "user","userOverallSys")

    ##########################

    # Setting the parameters of the hypothesis test
    #configMgr.nTOYs=5000
    configMgr.calculatorType=2 # 2=asymptotic calculator, 0=frequentist calculator
    configMgr.testStatType=3   # 3=one-sided profile likelihood test statistic (LHC default)
    configMgr.nPoints=20       # number of values scanned of signal-strength for upper-limit determination of signal strength.

    ##########################

    # Give the analysis a name
    configMgr.analysisName = "SimpleUL_%s" % SR
    configMgr.outputFileName = "results/%s_Output.root" % configMgr.analysisName

    # Define cuts
    configMgr.cutsDict["UserRegion"] = "1."

    # Define weights
    configMgr.weights = "1."

    # Define samples
    bkgSample = Sample("Bkg",kGreen-9)
    bkgSample.setStatConfig(False)
    bkgSample.buildHisto([nbkg],"UserRegion","cuts")
    #bkgSample.buildStatErrors([nbkgErr],"UserRegion","cuts")
    #bkgSample.addSystematic(corb)
    bkgSample.addSystematic(ucb)

    dataSample = Sample("Data",kBlack)
    dataSample.setData()
    dataSample.buildHisto([ndata],"UserRegion","cuts")

    # Define top-level
    ana = configMgr.addFitConfig("SPlusB")
    ana.addSamples([bkgSample,dataSample])
    #ana.setSignalSample(sigSample)

    # Define measurement
    meas = ana.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=lumiError)
    meas.addPOI("mu_SIG")
    meas.addParamSetting("Lumi",True,1)

    # Add the channel
    chan = ana.addChannel("cuts",["UserRegion"],1,0.,1.)
    chan.addDiscoverySamples(["SIG"], [1.], [0.], [1000.], [kMagenta])
    ana.setSignalChannels([chan])

    # These lines are needed for the user analysis to run
    # Make sure file is re-made when executing HistFactory
    if configMgr.executeHistFactory:
        if os.path.isfile("data/%s.root" % configMgr.analysisName):
            os.remove("data/%s.root" % configMgr.analysisName) 
