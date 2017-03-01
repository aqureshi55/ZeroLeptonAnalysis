#!/usr/bin/env python

from optparse import OptionParser
import sys, os, string, shutil, pickle, subprocess, copy, time

import os
import sys
import shutil

from ChannelsDict import *

from ZLFitterConfig import *
zlFitterConfig = ZLFitterConfig()

doBlind = zlFitterConfig.blindSR
lumi = zlFitterConfig.luminosity

parser = OptionParser()

parser.add_option("-a","--asymptotic", default=False, action="store_true", help="asymptotic")
parser.add_option("-m", "--merge", default=False, action="store_true", help="Merge results")
parser.add_option("-r", "--region", default="SRJigsawSRC1", help="Region")
parser.add_option("-o", "--output-dir", default="results/",
                  help="output dir under which files can be found", metavar="DIR")
(options, args) = parser.parse_args()

options.output_dir = os.path.abspath(options.output_dir)

# loop over analysis and compute UL,p0,...
for anaName in finalChannelsDict.keys():
#anaName = options.region

# names
    fileName = os.path.join(options.output_dir, "ZL_"+anaName+"_Discovery/Fit__Discovery_combined_NormalMeasurement_model.root")
    outName="UL_%s.tex" % (anaName)
    outNameTMP=outName+".tmp"

# options    
    nPoints = 30
    muRange = 100
    nToys = 1000

    if not("SRS3b" in anaName):
        continue
    
    if ("SRC" in anaName):
        muRange = 40
    if ("SRS1a" in anaName):
        muRange = 200
    elif ("SRS1b" in anaName):
        muRange = 100
    elif ("SRS" in anaName):
        muRange = 50
    else:
        muRange = 60
    

    option = "-N %d -R %d" % (nPoints, muRange)  
    if options.asymptotic:
        option += " -a"
    else:
        option += " -n %d"%(nToys) 

    cmd = "python $HISTFITTER/scripts/UpperLimitTable.py %s -c combined -p mu_SIG -w %s -l %f -o %s" % (option, fileName, lumi, outNameTMP)   

    cmd+="  >&"+anaName+"_UL.log &"

    print cmd
        
    if not options.merge:
        subprocess.call(cmd, shell=True)
    else:
        if os.path.exists(outNameTMP): 

            cmd="cat "+outNameTMP+"| sed -e 's/combined/"+anaName+"/g' > "+outName
            subprocess.call(cmd, shell=True)
        else:
            print "WARNING: file %s is missing! Waiting for a few seconds." % outNameTMP
            time.sleep(1)


######################################
#merge
######################################
if options.merge:
        
    debut="""\\begin{table}
    \\centering
    \\setlength{\\tabcolsep}{0.0pc}
    \\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lccccc}
    \\noalign{\\smallskip}\\hline\\noalign{\\smallskip}
    {\\bf Signal channel}                        & $\\langle\\epsilon{\\rm \\sigma}\\rangle_{\\rm obs}^{95}$[fb]  &  $S_{\\rm obs}^{95}$  & $S_{\\rm exp}^{95}$ & $CL_{B}$ & $p(s=0)$  \\\\
    \\noalign{\\smallskip}\\hline\\noalign{\\smallskip}
    %%
    """

    cmd = """cat  UL_SR*.tex | grep SR | grep \" &\" | sort -n"""
    res = os.popen(cmd).readlines()
    milieu=""
    for aLine in res:
        milieu+=aLine

    fin="""
    %
    \\noalign{\\smallskip}\\hline\\noalign{\\smallskip}
    \\end{tabular*}
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
    \\end{table}\n\n"""

    f=open("ULp0.tex","w")
    f.write(debut+milieu+fin)
    f.close()
