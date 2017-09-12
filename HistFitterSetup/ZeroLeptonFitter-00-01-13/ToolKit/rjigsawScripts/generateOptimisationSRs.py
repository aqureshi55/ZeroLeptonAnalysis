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

if os.getenv("ZEROLEPTONFITTER") is None:
    print("Cannot run without ZeroLeptonFitter setup!")
    sys.exit()

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import copy
from pprint import pprint

# from ChannelConfig import createChannelConfigFromString
#from ContourUtils import buildCutString
#from Utils import chunks
#from Utils import cartesian_product
from zerolepton.utils import chunks
from zerolepton.utils import cartesian_product
from zerolepton.contours.utils import buildCutString

from zerolepton.colors import colors

# simple script to generate a list of optimisation commands that would get run
# we take a set of passed cuts and generate a command from it

# example command:
# ./generateOptimisationSRs.py --grid=SM_GG_direct --point 1350_250 --nJets 2 --MET 200 --jetpt1 200 300 --jetpt2 200 300 --jetpt3 60 100 150 --jetpt4 60 100 150 --dPhi123 0.4 0.8 1.2 --aplanarity 0.02 --METovermeffNj 0.2 0.25 0.3 --meff 800 1000 1200 1400

# NOTE: new variables must always be added to zerolepton.contours.utils' buildCutString() method!

def loadGridPoints(grid):
    import pickle
    with open(os.path.join(os.getenv("ZEROLEPTONFITTER"), "ToolKit", "signalPointPickle.pkl"), "rb") as f:
        pointDict = pickle.load(f)

    if not grid in pointDict:
        print("Grid name {0} not in points dictionary".format(grid))
        sys.exit()

    return pointDict[grid].values()

parser = argparse.ArgumentParser()
parser.add_argument("--grid", default=None, type=str)
parser.add_argument("--point", default=[], type=str, action='append')
parser.add_argument("--pointsPerCommand", default=10, type=int)
parser.add_argument("--entire-grid", default=False, action="store_true", help="use entire grid")
parser.add_argument("--outputSuffix", default="", type=str, help="output filename suffix. A timestamp is used by default.")
parser.add_argument('--mode', choices=['discovery', 'exclusion', 'exclusionUL', 'all'], default='exclusion', help="Run discovery, exclusion or exclusion ULs to optimise")

args = parser.parse_args()
modeMap = {"discovery": "-z", "exclusion" : "-p", "exclusionUL" : "-p -l", "all" : "-z -p -l"}

grid = args.grid
points = args.point
mode = args.mode
print mode
outputSuffix = args.outputSuffix
pointsPerCommand = args.pointsPerCommand

# if grid is None:
#     print("Can't run without a grid!")
#     sys.exit()

if not mode in modeMap:
    print("Don't know HistFitter argument for optimisation mode {0}".format(mode))
    sys.exit()

from zerolepton.grids.config import GridConfig

# build a dictionary out of the cuts
cuts = vars(args)
# remove the useless args
cuts.pop('grid', None)
cuts.pop('point', None)
cuts.pop('mode', None)
cuts.pop('outputSuffix', None)
cuts.pop('entire_grid', None)
cuts.pop('pointsPerCommand', None)

#print cuts
#for var in cuts:
#    print var, len(cuts[var])

commands = []
cutStrings = set()
i = 0
nChunks = 1


SignalRegions = [

"SRJigsawSRG1a",
"SRJigsawSRG1b",
"SRJigsawSRG2a",
"SRJigsawSRG2b",
"SRJigsawSRG3a",
"SRJigsawSRG3b",
"SRJigsawSRG4", #MOR

"SRJigsawSRS1a",
"SRJigsawSRS1b",
"SRJigsawSRS2a",
"SRJigsawSRS2b",
"SRJigsawSRS3a",
"SRJigsawSRS3b",
"SRJigsawSRS4", #MOR

"SRJigsawSRC1",
"SRJigsawSRC2",
"SRJigsawSRC3",
"SRJigsawSRC4",
"SRJigsawSRC5",

]


for SignalRegion in SignalRegions:

    # for subsetPoints in chunks(points, pointsPerCommand):
    i = 1

    myCmds = []

    if mode == "exclusion":
        cmd = "../../../HistFitter-00-00-52/scripts/HistFitter.py -t -w -d -f  -F bkg  -V  -r {2} {3}/analysis/ZeroLepton_Run2_RJigsaw.py".format( 0,0, SignalRegion, os.getenv('ZEROLEPTONFITTER')   )
    elif mode == "discovery":
        cmd = "../../../HistFitter-00-00-52/scripts/HistFitter.py -t -w -f  -z -F disc -r {2} {3}/analysis/ZeroLepton_Run2_RJigsaw.py".format( 0,0, SignalRegion, os.getenv('ZEROLEPTONFITTER')   )

    commands.append( (SignalRegion, SignalRegion, cmd) )



timestamp = time.time()

data = {}
data['argv'] = " ".join(sys.argv)
data['mode'] = mode
data['timestamp'] = timestamp
data['grid'] = grid
data['points'] = points
data['commands'] = commands
data['cutStrings'] = list(cutStrings)

suffix = datetime.datetime.fromtimestamp(timestamp).strftime("%Y%m%d-%H%M%S")
if outputSuffix is not None and outputSuffix != "":
    suffix = outputSuffix

# store 2 files. we do this to prevent accidental deletion.
defaultDir = os.path.join(os.getenv('ZEROLEPTONFITTER'), 'optimisation')
if not os.path.exists(defaultDir): os.makedirs(defaultDir)
filenames = ["{0}/optimisation-{1}-{2}.json".format(defaultDir, grid, suffix)]
if os.getenv('PWD') != defaultDir and os.getenv('PWD') != os.getenv('ZEROLEPTONFITTER'):
    filenames.append("{0}/optimisation-{1}-{2}.json".format(os.getenv('PWD'), grid, suffix))

for filename in filenames:
	with open(filename, 'w') as outfile:
		json.dump(data, outfile)

	print(colors.OKGREEN + "Wrote {1} optimisation commands to {0}".format(filename, i) + colors.ENDC)




