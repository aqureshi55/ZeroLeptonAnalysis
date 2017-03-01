#!/usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libSusyFitter.so")

from optparse import OptionParser
from ROOT import TFile, RooWorkspace, TObject, TString, RooAbsReal, RooRealVar, RooFitResult, RooDataSet, RooAddition, RooArgSet, RooFormulaVar, RooAbsData, RooRandom
from ROOT import *
from ROOT import TMath, TMap, RooExpandedFitResult
#from ROOT import Util, TMath, TMap, RooExpandedFitResult
import sys, os, string, shutil, pickle, subprocess, copy, time

ROOT.gROOT.SetBatch(True)

doBlind=False
lumi=36.1

################################################################################

import os
import sys

parser = OptionParser()
parser.add_option("--asymptotic", default=False, action="store_true", help="asymptotic")
parser.add_option("-o", "--output-dir", default="results/",
                  help="output dir under which files will be stored", metavar="DIR")
(options, args) = parser.parse_args()
options.output_dir += "/" #to be sure

muRanges = {}
# muRanges["RJR-S1a"] = 160
# muRanges["RJR-S1b"] = 120
# muRanges["RJR-S2a"] = 60
# muRanges["RJR-S2b"] = 50
# muRanges["RJR-S3a"] = 40
# muRanges["RJR-S3b"] = 30
# muRanges["RJR-S4"] = 30
# muRanges["RJR-G1a"] = 40
# muRanges["RJR-G1b"] = 25
# muRanges["RJR-G2a"] = 30
# muRanges["RJR-G2b"] = 30
# muRanges["RJR-G3a"] = 20
# muRanges["RJR-G3b"] = 20
# muRanges["RJR-G4"] = 20
# muRanges["RJR-C1"] = 25
# muRanges["RJR-C2"] = 50
# muRanges["RJR-C3"] = 60
# muRanges["RJR-C4"] = 25
# muRanges["RJR-C5"] = 20
muRanges["SRJigsawSRS1a"] = 160
muRanges["SRJigsawSRS1b"] = 120
muRanges["SRJigsawSRS2a"] = 60
muRanges["SRJigsawSRS2b"] = 50
muRanges["SRJigsawSRS3a"] = 40
muRanges["SRJigsawSRS3b"] = 30
muRanges["SRJigsawSRS4"] = 30
muRanges["SRJigsawSRG1a"] = 40
muRanges["SRJigsawSRG1b"] = 25
muRanges["SRJigsawSRG2a"] = 30
muRanges["SRJigsawSRG2b"] = 30
muRanges["SRJigsawSRG3a"] = 20
muRanges["SRJigsawSRG3b"] = 20
muRanges["SRJigsawSRG4"] = 20
muRanges["SRJigsawSRC1"] = 25
muRanges["SRJigsawSRC2"] = 50
muRanges["SRJigsawSRC3"] = 60
muRanges["SRJigsawSRC4"] = 25
muRanges["SRJigsawSRC5"] = 20

regions = sorted(muRanges.keys())

#loop over analysis and compute UL,p0,...
for anaShortname in regions:
    cmd_HF = "HistFitter.py -F disc -t -w -f -r {0} {1}/ToolKit/simpleUL/SimpleUL.py &> {0}.log".format(anaShortname, os.getenv("ZEROLEPTONFITTER")) 
    #print cmd_HF
    #subprocess.call(cmd_HF, shell=True)

    fileName = "%s/SimpleUL_%s/SPlusB_combined_NormalMeasurement_model.root" % (options.output_dir, anaShortname)

    outName="UL_%s.tex" % (anaShortname)
    outNameTMP=outName+".tmp"

    nPoints = 40
    nToys = 1000
    muRange = muRanges[anaShortname]
    
    option = "-N %d -R %d" % (nPoints, muRange) 

    if options.asymptotic:
        option += " -a"#asymptotic instead of toys

    cmd_UL = "python $HISTFITTER/scripts/UpperLimitTable.py {0} -c combined -n {1:d} -p mu_SIG -w {2} -l {3:.3f} -o {4} >& UL_{5}.log".format(option, nToys, fileName, lumi, outNameTMP, anaShortname)   
    cmd_pval = "cd %s/.. && HistFitter.py -z -r %s $ZEROLEPTONFITTER/ToolKit/simpleUL/SimpleUL.py && cd -" % (options.output_dir, anaShortname)
    
    cmd_rename ="cat "+outNameTMP+"| sed -e 's/combined/"+anaShortname+"/g' > "+outName

    cmd = "{0} && {1} && {2} &".format(cmd_HF, cmd_UL, cmd_rename)
    subprocess.call(cmd, shell=True)
        

######################################
#merge
######################################

debut="""\\begin{table}
\\begin{center}
\\setlength{\\tabcolsep}{0.0pc}
\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lccccc}
\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}
{\\bf Signal channel}                        & $\\langle\\epsilon{\\rm \\sigma}\\rangle_{\\rm obs}^{95}$[fb]  &  $S_{\\rm obs}^{95}$  & $S_{\\rm exp}^{95}$ & $CL_{B}$ & $p(s=0)$  \\\\
\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}
%%
"""

cmd = """cat  UL_*.tex | grep RJR | grep \" &\" | sort -n"""
res = os.popen(cmd).readlines()
milieu=""
for aLine in res:
    milieu+=aLine

fin="""
%
\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}
\\end{tabular*}
\\end{center}
\\caption[Breakdown of upper limits.]{
Left to right: 95\\% CL upper limits on the visible cross section
($\\langle\\epsilon\\sigma\\rangle_{\\rm obs}^{95}$) and on the number of
signal events ($S_{\\rm obs}^{95}$ ).  The third column
($S_{\\rm exp}^{95}$) shows the 95\\% CL upper limit on the number of
signal events, given the expected number (and $\\pm 1\\sigma$
excursions on the expectation) of background events.
The last two columns
indicate the $CL_B$ value, i.e. the confidence level observed for
the background-only hypothesis, and the discovery $p$-value ($p(s = 0)$). 
\\label{table.results.exclxsec.pval.upperlimit}}
\\end{table}"""

f=open("ULp0.tex","w")
f.write(debut+milieu+fin)
f.close()
