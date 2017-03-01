#!/bin/env python

import os,sys,commands;
import ROOT;
from optparse import OptionParser;

# set batch mode
ROOT.gROOT.SetBatch(True);

RegionName = "SRAll"

parser = OptionParser();
parser.add_option("-i","--inputDir" ,action="store",dest="inputDir" ,default="/n/atlasfs/atlasdata/atlasdata1/crogan/SAMPLES_ATLAS/NTUPLES/0LEP/v115_sys/");
parser.add_option("-o","--outputname",action="store",dest="outputname",default="makeSignalPointPickle.py");

inputSamples=["GG_direct","SS_direct","GG_onestepCC","SS_onestepCC","SM_GG_N2","GG_onestepN2C","SS_onestepN2C"];

(options, args) = parser.parse_args();


outtxtfile = open(options.outputname,"w");

header="\
#!/usr/bin/env python\n\
\n\
import pickle\n\
\n\
outfile = open('signalPointPickle.pkl','wb')\n\
pointdict = {};\n\
\n\
";

outtxtfile.write(header);

for sample in inputSamples :

  outtxtfile.write("#"+sample+"\n");
  outtxtfile.write("pointdict[\'"+sample+"\'] = {\n");
  
  infile=ROOT.TFile(options.inputDir+"/"+sample+".root");
  infile.ls(sample+"*_SRAll");
  keyList = infile.GetListOfKeys();
  keyNames=[];
  for i in range(keyList.GetEntries()):
    keyNames.append((keyList.At(i).GetName()));
  keyNames.sort();
  
  #print "keyNames =",keyNames;
  tempCounter = 1
  for keyName in keyNames:
   #print keyName;
   if keyName.startswith(sample) and keyName.endswith(RegionName):
     tmpstr  = keyName.replace(sample,"");
     tmpstr  = tmpstr.replace(RegionName,"");
     tmpstr  = tmpstr.strip("_");
     print "mass string : "+tmpstr;
     massesStr = tmpstr.split("_");
     #masses    = [ int(massStr) for massStr in massesStr ] ;
  
     # get run number
     print "treename = "+keyName;
     tree = infile.Get(keyName);
     tree.GetEntry(0);
     runNumber=tree.RunNumber;
     runNumber = tempCounter
     if len(massesStr)==2 : 
       outtxtfile.write("  "+(str)(runNumber)+": ("+massesStr[0]+","+massesStr[1]+"),\n");
     elif len(massesStr)==3 :
       outtxtfile.write("  "+(str)(runNumber)+": ("+massesStr[0]+","+massesStr[1]+","+massesStr[2]+"),\n");
       pass;
     
     pass; # end of matched tree selection
     tempCounter = tempCounter + 1
  outtxtfile.write("}\n"); 
  outtxtfile.write("\n"); 
  pass; # end of loop over signal samples


footer="pickle.dump(pointdict,outfile)";
outtxtfile.write(footer);

outtxtfile.close();
  

