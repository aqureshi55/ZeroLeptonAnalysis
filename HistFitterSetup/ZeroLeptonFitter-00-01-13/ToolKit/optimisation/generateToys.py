#!/usr/bin/env python

import argparse
import datetime
import itertools
import json
import math
import os
import socket
import sys
import textwrap
import time
from optparse import OptionParser

if os.getenv("ZEROLEPTONFITTER") is None:
    print("Cannot run without ZeroLeptonFitter setup!")
    sys.exit()

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import copy
from pprint import pprint

from ChannelConfig import createChannelConfigFromString
from zerolepton.utils import chunks
from zerolepton.utils import cartesian_product
from zerolepton.contours.utils import buildCutString

from zerolepton.colors import colors

from ZLFitterConfig import *
zlFitterConfig = ZLFitterConfig()

doBlind = zlFitterConfig.blindSR
lumi = zlFitterConfig.luminosity

parser = OptionParser()

parser.add_option("-o", "--output-dir", default="results/",
                  help="output dir under which files can be found", metavar="DIR")
(options, args) = parser.parse_args()
options.output_dir = os.path.join(os.path.abspath(options.output_dir),"results/")

from zerolepton.grids.config import GridConfig

SignalRegions = [

"SRJigsawSRG1a",
"SRJigsawSRG1b",
"SRJigsawSRG2a",
"SRJigsawSRG2b",
"SRJigsawSRG3a",
"SRJigsawSRG3b",
"SRJigsawSRG4",

"SRJigsawSRS1a",
"SRJigsawSRS1b",
"SRJigsawSRS2a",
"SRJigsawSRS2b",
"SRJigsawSRS3a",
"SRJigsawSRS3b",
"SRJigsawSRS4",

"SRJigsawSRC1",
"SRJigsawSRC2",
"SRJigsawSRC3",
"SRJigsawSRC4",
"SRJigsawSRC5",

]

# options                                                                                                                                                     
nPoints = 30
muRange = 100
nToys = 1000


commands = []
for SignalRegion in SignalRegions:
    fileName = os.path.join(options.output_dir, "ZL_"+SignalRegion+"_Discovery/Fit__Discovery_combined_NormalMeasurement_model.root")
    outName="UL_%s.tex" % (SignalRegion)
    outNameTMP=outName+".tmp"
    if ("SRC" in SignalRegion):
        muRange = 40
    if ("SRS1a" in SignalRegion):
        muRange = 200
    elif ("SRS1b" in SignalRegion):
        muRange = 100
    elif ("SRS" in SignalRegion):
        muRange = 50
    else:
        muRange = 60

    option = "-N %d -R %d" % (nPoints, muRange)
    option += " -n %d"%(nToys)

    cmd = "python $HISTFITTER/scripts/UpperLimitTable.py %s -c combined -p mu_SIG -w %s -l %f -o %s" % (option, fileName, lumi, outNameTMP)

    commands.append( (SignalRegion, SignalRegion, cmd) )

data = {}
timestamp = time.time()
data['commands'] = commands
data['timestamp'] = timestamp

suffix = datetime.datetime.fromtimestamp(timestamp).strftime("%Y%m%d-%H%M%S")

# store 2 files. we do this to prevent accidental deletion.
defaultDir = os.path.join(os.getenv('ZEROLEPTONFITTER'), 'toys')
if not os.path.exists(defaultDir): os.makedirs(defaultDir)
filenames = ["{0}/toys-{1}.json".format(defaultDir, suffix)]
if os.getenv('PWD') != defaultDir and os.getenv('PWD') != os.getenv('ZEROLEPTONFITTER'):
    filenames.append("{0}/toys-{1}.json".format(os.getenv('PWD'), suffix))

for filename in filenames:
	with open(filename, 'w') as outfile:
		json.dump(data, outfile)

	print(colors.OKGREEN + "Wrote {1} toy commands to {0}".format(filename, len(commands)) + colors.ENDC)
