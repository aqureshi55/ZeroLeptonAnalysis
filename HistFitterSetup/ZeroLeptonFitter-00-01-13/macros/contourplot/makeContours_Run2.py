#!/usr/bin/env python
# usage : ./makeCountours_Run2.py --all --grid <gridname>

import pprint
import time
import ROOT
import socket
import glob
import sys
import os

from array import array
from ROOT import *
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

if not os.getenv("HISTFITTER"): 
    print "$HISTFITTER is not defined! Exiting now."
    sys.exit()
    pass;

if not os.getenv("ZEROLEPTONFITTER"): 
    print "$ZEROLEPTONFITTER is not defined! Exiting now."
    sys.exit()
    pass;
    
from ROOT import TGraph
ROOT.gSystem.Load("libSusyFitter.so");

ROOT.gROOT.LoadMacro("$ZEROLEPTONFITTER/macros/contourplot/ContourUtils.C");
ROOT.gROOT.LoadMacro("$ZEROLEPTONFITTER/macros/contourplot/contourmacros/GetSRName.C");
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Import Modules                                                             # 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
import sys, os, string, shutil,pickle,subprocess
from ChannelsDict import *
from ZLFitterConfig import *


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Some global variables
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# CHECK:
#INPUTDIR="%s/results" % (os.getenv("ZEROLEPTONFITTER"))
#INPUTDIR="%s/optimisation/optimisation-GG_onestepCC-20170218-224212/results" % (os.getenv("ZEROLEPTONFITTER"))
INPUTDIR="%s/optimisation/optimisation-SS_onestepCC-20170219-002955/results" % (os.getenv("ZEROLEPTONFITTER"))
#INPUTDIR="%s/optimisation/optimisation-GG_direct-20170214-190428/results" % (os.getenv("ZEROLEPTONFITTER"))

#INPUTDIR="%s/optimisation/optimisation-SS_direct-20170214-124632/results" % (os.getenv("ZEROLEPTONFITTER"))

# OUTPUTDIR is where the combination of these using hadd will go, as well as 
#           all the list files for the contours, the histograms and the plots
OUTPUTDIR = "Outputs/"

# Cross sections to use. (Up, down is the theory uncertainty)
# Our plotting always uses Nominal for exp+obs+yellow band, up and down only for two extra obs curves
#allXS=["Nominal"]
allXS=["Nominal", "Up", "Down"]

###########################################################################
# useful functions
###########################################################################

def parseCmdLine(args):
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-m", dest="doMerge", help="merge output root files", action='store_true', default=False)
    parser.add_option("-c", dest="makeContours", help="create contours", action='store_true', default=False)
    parser.add_option("-o", dest="doOring", help="Analysis Oring", action='store_true', default=False)
    parser.add_option("-p", dest="makePlots", help="create plots", action='store_true', default=False)
    parser.add_option("--inputDir", dest="inputDir", help="input directory", default=INPUTDIR)
    parser.add_option("--outputDir", dest="outputDir", help="output directory", default=OUTPUTDIR)
    parser.add_option("--all", dest="doAll", help="do all steps", action='store_true', default=False)
    parser.add_option("--grid", dest="grid", help="grid name SS_direct, GG_direct or GG_onestepCC or SS_onestepCC, SM_GG_N2", default="GG_direct")
    parser.add_option("--suffix", dest="suffix", help="suffix to append after grid name in output files (default empty)", default="")
    parser.add_option("--match", dest="match", help="name to match input files against", default="")
    parser.add_option("--filter", dest="filter", help="name to filter input files against", default="")
    parser.add_option("--ul", dest="makeUL", help="do upper limits", action='store_true', default=False)
    parser.add_option("--merge-ul-only", help="only merge UL files", action='store_true', default=False)
    parser.add_option("--discovery", dest="discovery", help="use discovery lines", action='store_true', default=False)

    (config, args) = parser.parse_args(args)

    # make output directory
    try:
      os.mkdir(config.outputDir)
    except:
      pass
    if config.outputDir[-1] != '/':
      config.outputDir+='/'
      pass;
    print "Output directory : %s" % (config.outputDir);

    if config.filter != "" and config.match != "" and (config.filter.find(config.match) != -1 or config.match.find(config.filter) != -1):
        print "--match and --filter overlap and negate each other -> zero input files by default! Exiting now."
        sys.exit()
        pass;
    
    if config.makeUL and config.discovery:
        print "--ul and --discovery cannot be used simultaneously!"
        sys.exit()
        pass;

    if not config.doAll and not config.doMerge and not config.makeContours and not config.doOring and not (config.makePlots or config.makeATLASplots):
        print "No step to execute specified!"
        pass;

    print "Grid name: ", config.grid, ", Output name: config.outputName"    
        
    config.outputName = config.grid + config.suffix

    return config
# end of def parseCmdLine()

def wait(sec):
    os.system('setterm -cursor off')
    while sec > 0:
        sys.stdout.write(str(sec) + '...     \r')
        sys.stdout.flush()
        sec -= 1
        try:
            time.sleep(1)
        except KeyboardInterrupt:
            os.system('setterm -cursor on')
            print
            sys.exit()
            pass;
        pass;
    os.system('setterm -cursor on')
# end of def wait()

def getFileList(inputdir):
    cachefile = os.getenv("ZEROLEPTONFITTER")+"/macros/contourplot/results.cache"

    readDir = False

    print "Read filename from cache file : "+cachefile;
    # Add file names in cachefile(results.cache) to instance "filenames"
    if os.path.exists(cachefile) and os.path.isfile(cachefile):
        f = open(cachefile)
        filenames = [l.strip() for l in f.readlines() if l.strip()]
        f.close()

        # if filenames exist in cachefile(results.cache), skip read directory
        if len(filenames) == 0: readDir = True
        else                  : return filenames
    else:
        readDir = True
        pass;

    print "Read filename existing in "+inputdir;
    # Read directory in inputdir
    if readDir:
        dirnames = os.listdir(inputdir)
        filenames = []
        i=1
        print "Found %d directories, reading them all..." % len(dirnames)
        for d in dirnames:
            sys.stdout.write('%d / %d \r' % (i, len(dirnames)))
            sys.stdout.flush()
            
            # get file name in inputdir/*/
            fnames = os.listdir(os.path.join(inputdir, d))
            for f in fnames:
                filenames.append(os.path.join(inputdir, d, f))
                pass;
            i+=1
            pass;

        # write file names to cachefile(results.cache)
        f = open(cachefile, "w")
        for n in filenames:  f.write("%s\n" % n)
        f.close()
    else:
        print "INFO: read %d lines from results.cache (file modified: %s)" % (len(filenames), time.ctime(os.path.getmtime(cachefile)))
        print "Waiting 3 seconds in case this is not correct and you want to delete the file..."
        wait(3)
        pass;

    return filenames

# end of def getFileList()

def MergeFiles(config):
    # get file names in inputdir/*/ or in cachefile
    filenames = getFileList(config.inputDir) 
   
    grid_name=config.grid # SS_direct, GG_direct, GG_onestepCC, SS_onestepCC, SM_GG_N2, GG_onestep_mlsp60
    name_match=config.match   # special input file match-selectoin by argument (--match)
    name_filter=config.filter # special input file filtre-selectoin by argument (--filter)

    # used : (GG_direct, SS_direct, GG_onestepCC, SS_onestepCC, SM_GG_N2)
    if grid_name in ["GG_direct" , "SS_direct" , "GG_onestepCC"  , "SS_onestepCC" , "SM_GG_N2", "GG_onestep_mlsp60"] :
        # set anaList (signal region list)
        config.anaList = finalChannelsDict.keys()
    else:
        print "correct grid_name is not defined. (grid_name=%s)" % (grid_name)
        config.anaList = finalChannelsDict.keys()
        pass;

    #merge histfitter output root files in each SRs
    for ana in config.anaList:
        print "  MergeFiles : start merge of %s" % (ana);
        thisAnaFilenames = []

        # First, filter filenames to suitable files
        for filename in filenames:
            # skip files other than root file
            if not filename.endswith(".root"):
                continue;
                pass;
            # select filename containing ana (SR name)
            if filename.find(ana+"_") == -1: 
                continue;
                pass;
            # select filename containing grid_name (GG_direct/SS_direct/GG_onestepCC/SS_onestepCC/SM_GG_N2)
            if filename.find(grid_name+"_") == -1:
              if not ( grid_name=="GG_onestep_mlsp60" and "GG_onestepCC" in filename ) : 
                continue;
                pass;
              pass;
            # select filename containing match
            if name_match != "" and filename.find(name_match) == -1:
                continue;
                pass;
            # select filename not containing filter
            if name_filter != "" and filename.find(name_filter) != -1: 
                continue;
                pass;

            # append only suitable file
            thisAnaFilenames.append(filename)

            pass;
      
        # initialize thisFilenames[xs] for cross-section
        # loop over cross-section on (nominal or up or down)
        thisFilenames = {}
        for xs in allXS:
            thisFilenames[xs] = []
            pass;

        thisULFilenames = []
        # Now, split filenames into UL files and each xs of hypotest files
        for filename in thisAnaFilenames:
            # get upperlimit files
            if filename.find("upperlimit") != -1:
                thisULFilenames.append(filename)
            # get hypotest files
            else:
                for xs in allXS:
                    # check cross-section nominal/up/down
                    if filename.find(xs) != -1:
                        thisFilenames[xs].append(filename)
                        break;
                        pass;
                    pass;
                pass;

        # merge hypotest files
        for xs in allXS:
            # skip if only merge upperlimit option is declared
            if config.merge_ul_only:
                print "Not merging %s, %s, %s, --merge-only-ul used" % (config.grid, ana, xs)
            # skip if the number of hypotest files is 0
            elif len(thisFilenames[xs]) == 0:
                print "No files for %s, %s, %s -> skipping" % (config.grid, ana, xs)
            # proceed to use hypotest files
            else:
                maxfilenum=200
                # default
                if len(thisFilenames[xs])<maxfilenum:
                    inputFiles = " ".join(thisFilenames[xs])
                    # outputName=config.grid+config.suffix
                    outputFilename = config.outputDir+config.outputName+"_"+ana+"_fixSigXSec"+xs+".root"
                    # hadd
                    cmd="hadd -f "+outputFilename+" "+inputFiles
                    print "  MergeFiles :",cmd;
                    subprocess.call(cmd, shell=True)
                # just separate merge process to 2step when files are too much
                # (make former/later merged files and after that merge them to one file) 
                else:
                    inputFiles = " ".join(thisFilenames[xs][:maxfilenum])
                    outputFilename1 = config.outputDir+config.outputName+"1_"+ana+"_fixSigXSec"+xs+".root"
                    cmd="hadd -f "+outputFilename1+" "+inputFiles 
                    print "  MergeFiles :",cmd;
                    subprocess.call(cmd, shell=True)
                    inputFiles = " ".join(thisFilenames[xs][maxfilenum:])
                    outputFilename2 = config.outputDir+config.outputName+"2_"+ana+"_fixSigXSec"+xs+".root"
                    cmd="hadd -f "+outputFilename2+" "+inputFiles 
                    print "  MergeFiles :",cmd;
                    subprocess.call(cmd, shell=True)
                    outputFilename = config.outputDir+config.outputName+"_"+ana+"_fixSigXSec"+xs+".root"            
                    cmd="hadd -f "+outputFilename+" "+outputFilename1+" " +outputFilename2
                    print "  MergeFiles :",cmd;
                    subprocess.call(cmd, shell=True)
                    # remove temporary former/later merged files
                    os.remove(outputFilename1)
                    os.remove(outputFilename2)
                    pass;
                pass;
            pass;

        # merge upperlimit files
        if config.makeUL:
            if len(thisULFilenames) == 0:
                print "No UL files for %s, %s -> skipping" % (config.grid, ana)
            else:
                maxfilenum2=200
                # default
                if len(thisULFilenames)<maxfilenum2:
                    inputFiles = " ".join(thisULFilenames)
                    # outputName=config.grid+config.suffix
                    outputFilename = config.outputDir+config.outputName+"_"+ana+"_upperlimit.root"
                    # hadd
                    cmd="hadd -f "+outputFilename+" "+inputFiles 
                    print "  MergeFiles :",cmd;
                    subprocess.call(cmd, shell=True)
                # just separate merge process to 2step when files are too much
                # (make former/later merged files and after that merge them to one file) 
                else:
                    inputFiles = " ".join(thisULFilenames[:maxfilenum2])
                    outputFilename1 = config.outputDir+config.outputName+"1_"+ana+"_upperlimit.root"
                    cmd="hadd -f "+outputFilename1+" "+inputFiles 
                    print "  MergeFiles :",cmd;
                    subprocess.call(cmd, shell=True)
                    inputFiles = " ".join(thisULFilenames[maxfilenum2:])
                    outputFilename2 = config.outputDir+config.outputName+"2_"+ana+"_upperlimit.root"
                    cmd="hadd -f "+outputFilename2+" "+inputFiles 
                    print "  MergeFiles :",cmd;
                    subprocess.call(cmd, shell=True)
                    outputFilename = config.outputDir+config.outputName+"_"+ana+"_upperlimit.root"            
                    cmd="hadd -f "+outputFilename+" "+outputFilename1+" " +outputFilename2
                    print "  MergeFiles :",cmd;
                    subprocess.call(cmd, shell=True)
                    # remove temporary former/later merged files
                    os.remove(outputFilename1)
                    os.remove(outputFilename2)
                    pass;
                pass;
            pass;
        pass;

    print " MergeFiles : END!";

# end of def MergeFiles()


def MakeContours(config):
    automaticRejection = False # automatic rejection of bad points in HistFitter

    # loop over config.anaList ( SRs )
    for ana in config.anaList:        
        print "MakeContours:", ana

        # loop over cross-section nominal or up or down
        for xs in allXS:
            # skip up/down when --discovery
            if config.discovery and xs != "Nominal":
                print "--discovery only uses nominal values, ignoring xsec %s" % xs
                continue;
                pass;
            
            # hypotest root file base name
            basename = config.outputDir+config.outputName+"_"+ana+"_fixSigXSec"+xs
            # upperlimit root file base name
            basenameUL = config.outputDir+config.outputName+"_"+ana+"_upperlimit"  
            
            # default cut
            cutStr = "1"; # accept everything
            
            # Format of the hypotests, normally follows e.g. hypo_<grid>_2000_1000 
            # default : used in GG/SS_direct ..
            format     = "hypo_"+config.grid+"_%f_%f";
            interpretation = "m0:m12";            
            
            if config.grid.find("SS_onestepCC")!=-1:
                format     = "hypo_SS_onestepCC_%f_%f_%f";
                interpretation = "msquark:mchargino:mlsp";
                pass;
            if config.grid.find("GG_onestepCC")!=-1:
                format     = "hypo_GG_onestepCC_%f_%f_%f";
                interpretation = "mgluino:mchargino:mlsp";
                pass;
            if config.grid.find("SM_GG_N2")!=-1:
                format     = "hypo_SM_GG_N2_%f_%f_%f";
                interpretation = "mgluino:mlsp2:mlsp";
                pass;
            if config.grid.find("GG_onestep_mlsp60")!=-1:
                format     = "hypo_GG_onestepCC_%f_%f_%f";
                interpretation = "mgluino:mchargino:mlsp";
                pass;
            if config.discovery:
                format = format.replace("hypo", "hypo_discovery")
                pass;
            
            print "INFO: format set to %s, %s" % (format, interpretation)

            listSuffix = "__1_harvest_list"

            # modification of mass filter
            # onestep grids want mlsp != 60
            print config.grid;
            if (config.grid.find("GG_onestepCC")!=-1) or  (config.grid.find("SS_onestepCC")!=-1):
                cutStr = "mlsp!=60"
                listSuffix = "__mlspNE60_harvest_list"
                print "removed mlsp!=60 tentatively" 
                pass;
            if (config.grid.find("GG_onestep_mlsp60")!=-1):
                cutStr = "mlsp==60"
                listSuffix = "__mlspEE60_harvest_list"
                print "select only mlsp==60 tentatively" 
                pass;

            if config.discovery:                
                listSuffix = "_discovery_1_harvest_list"
                pass;

            listSuffix+=".json"
                
            inputfile = basename+".root"
            print "MakeContours: inputfile name", inputfile
            if os.path.isfile(inputfile):
                fitResultFormat = format;
                print "CollectAndWriteHypoTestResults( %s, %s, %s, %s, %i, %s )" % ( inputfile, fitResultFormat, interpretation, cutStr, int(automaticRejection), config.outputDir ) ;
                CollectAndWriteHypoTestResults( inputfile, fitResultFormat,  interpretation, str(cutStr), int(automaticRejection), config.outputDir ) ;
                """
                cmd = "root -l -b -q \"$ZEROLEPTONFITTER/macros/contourplot/makelistfiles.C(\\\""+inputfile+"\\\",\\\""+fitResultFormat+"\\\",\\\""+interpretation+"\\\",\\\""+str(cutStr)+"\\\","+str(int((automaticRejection)))+",\\\""+config.outputDir+"\\\")\"";
                print cmd;
                ret=subprocess.check_call(cmd)
                ret=subprocess.call(cmd, shell=True)
                print "process result :",ret;
                """
                

            # get extra information from upper limits computation
            # and merge the files
            if config.makeUL and xs == "Nominal":
                inputfile = basenameUL+".root"
                if os.path.isfile(inputfile):
                    # by defition, ULs are not discovery -> don't care about passing -d
                    CollectAndWriteHypoTestResults( inputfile, format, interpretation, cutStr, int(automaticRejection), config.outputDir ) ;
                    
                    # merge the output files for hypotest and upperlimit in 1 output file
                    mergeFileList(config, basename+listSuffix, basenameUL+listSuffix)
                    
                    pass;
                pass;
            

            if not os.path.exists(basename+listSuffix):
                print "INFO: file %s does not exist, skipping call to makecontourhists.C" % (basename+listSuffix)
                continue;
            
            # convert : json file --> plain text file
            cmd = "GenerateTreeDescriptionFromJSON.py -f %s" % (basename+listSuffix)
            print cmd;
            subprocess.call(cmd, shell=True)
            
            # mv summary_harvest_tree_description.* from config.outputDir to here
            cmd="cp -v "+config.outputDir+"summary_harvest_tree_description.py ./summary_harvest_tree_description.py"
            subprocess.call(cmd, shell=True)
            cmd="cp -v "+config.outputDir+"summary_harvest_tree_description.h  ./summary_harvest_tree_description.h"
            subprocess.call(cmd, shell=True)
           
            # convert : plain text file --> root file
            listSuffix = "".join(listSuffix.split(".")[:-1]) # remove ".json"
            cmd = "root -l -b -q \"$ZEROLEPTONFITTER/macros/contourplot/makecontourhists.C(\\\""+basename+listSuffix+"\\\",\\\""+config.grid+"\\\")\""
            print cmd
            subprocess.call(cmd, shell=True)
                        
            cmd = "mv -v *_list.root "+config.outputDir
            print cmd
            subprocess.call(cmd, shell=True)

            pass;
        pass;
    pass;
# end of def MakeContoures()


# Merge regions into combination based on best region (simply pick lowest p-value)
def Oring(config):
    # print which SRs you have to put in GetSRName.C
    for indx,ana in enumerate(config.anaList):
        print indx, ana
        pass;

    import sys, os, string, shutil, pickle, subprocess    
    import ROOT

    from summary_harvest_tree_description import treedescription
    dummy,description = treedescription()
    allpar = description.split(':')
    print allpar;

    ROOT.gROOT.Reset()
    ROOT.gROOT.SetBatch(True)

    # For all xsecs, merge on best selectpar (normally expected CLs)
    for xsecStr in allXS:
        myMap = {}

        # skip when xsec != Nominal in discovery mode
        if config.discovery and xsecStr != "Nominal":
            continue;

        # if no UL's available, merge on CLsexp. Otherwise, use expectedUpperLimit
        # you should use CLsexp even if it's for observed limit because it can be bias?
        selectpar = "CLsexp"
        #selectpar = "CLs" 

        if config.makeUL and xsecStr == "Nominal": # ULs only available for Nominal
            selectpar = "expectedUpperLimit"
            pass;

        # select best SR based on p0exp    
        if config.discovery:
            selectpar = "p0exp"
            pass;

        print "Using selectpar = %s" % selectpar

        # default
        par1_s = "m0"
        par2_s = "m12"
        par3_s = ""
        listSuffix = "__1_harvest_list"

        # modifification depending on signal grid
        if config.grid.find("GG_onestepCC")!=-1:
            print "Ordering GG_onestepCC"
            par1_s = "mgluino"
            par2_s = "mchargino"
            par3_s = "mlsp"
            listSuffix = "__mlspNE60_harvest_list"
            pass;
        if config.grid.find("SS_onestepCC")!=-1:
            print "Ordering SS_onestepCC"
            par1_s = "msquark"
            par2_s = "mchargino"
            par3_s = "mlsp"
            listSuffix = "__mlspNE60_harvest_list"
            pass;
        if config.grid.find("SM_GG_N2")!=-1 :
            print "Ordering SM_GG_N2"
            par1_s = "mgluino"
            par2_s = "mlsp2"
            par3_s = "mlsp"
            pass;
        if config.grid.find("GG_onestep_mlsp60")!=-1:
            print "Ordering GG_onestep_mlsp60"
            par1_s = "mgluino"
            par2_s = "mchargino"
            par3_s = "mlsp"
            listSuffix = "__mlspEE60_harvest_list"
            pass;


        if config.discovery:
            listSuffix = "_discovery_1_harvest_list"
            pass;
            
        infoFilename = config.outputDir+config.outputName+"_combined_fixSigXSec"+xsecStr+listSuffix+"_infoFile"
        file_info = open(infoFilename,"w")
     
        # loop over ana (SRs)
        print " Oring :",config.anaList;
        for indx,ana in enumerate(config.anaList):

            filename = config.outputDir+config.outputName+"_"+ana+"_fixSigXSec"+xsecStr+listSuffix
            print filename
            infoline="%i : %s\n"%(indx+1,ana)
            file_info.write(infoline)
            if not os.path.exists(filename):
                print "file does not exist -> skip"
                continue;
            
            infile = open(filename,'r')
            for line in infile.readlines():
                vals = line.strip().split(' ')
                if len(allpar) != len(vals): 
                    print 'PRB!!!!!!!!!!!!!!!!!!!!'
                    print "summary file says %d components; file has %d per line" % (len(allpar),len(vals))
                    continue;
                
                vals[allpar.index("fID/C")]=(indx+1)
 
                pval = float( vals[allpar.index(selectpar+"/F")])
                par1 = float( vals[allpar.index(par1_s+"/F")])
                par2 = float( vals[allpar.index(par2_s+"/F")])
                key = "%d_%d" % (par1, par2)
                
                if config.grid.find("GG_onestepCC")!=-1: 
                    par3 = float( vals[allpar.index(par3_s+"/F")])
                    key += "_%d" % par3
                    pass;
                if config.grid.find("SS_onestepCC")!=-1: 
                    par3 = float( vals[allpar.index(par3_s+"/F")])
                    key += "_%d" % par3
                    pass;
                if config.grid.find("SM_GG_N2")!=-1: 
                    par3 = float( vals[allpar.index(par3_s+"/F")])
                    key += "_%d" % par3
                    pass;
                if config.grid.find("GG_onestep_mlsp60")!=-1: 
                    par3 = float( vals[allpar.index(par3_s+"/F")])
                    key += "_%d" % par3
                    pass;
           
                #print "DEBUG: %s selectpar=%.2e for %s" % (ana, pval, key)
                # ignore negative pvalue
                if pval < 0:
                    print "INFO: %s removing negative selectpar (%s = %.2e) for %s" % (ana, selectpar, pval, key)
                    continue;
                # ignore expectedUpperLimit < 0.00001 
                if selectpar == "expectedUpperLimit" and pval < 0.00001:
                    print "INFO: %s removing %s < 0.00001 for %s" % (ana, selectpar, key)
                    continue;
                # 110 is 20 times our default step -> this is almost certainly a bug
                if selectpar == "expectedUpperLimit" and pval == 110.0:
                    print "INFO: %s removing expUL==110.0 for %s" % (ana, key)
                    continue;
                # if -1sig, -2sig == 0 and +1sig, 2sig == 100 -> almost certainly a bug too
                if selectpar == "expectedUpperLimit" and float(vals[allpar.index("expectedUpperLimitMinus1Sig")]) == 0.0 and float(vals[allpar.index("expectedUpperLimitMinus2Sig")]) == 0.0 and float(vals[allpar.index("expectedUpperLimitPlus1Sig")]) == 100.0 and float(vals[allpar.index("expectedUpperLimitPlus2Sig")]) == 100.0:
                    print "INFO: %s removing point %s with expULMinus1Sig == expULMinus2Sig == 0 and expULPlus1Sig == expULPlus2Sig == 100" % (ana, key)
                    continue;
                # ignore observed upperLimit=0 when merging on UL
                if selectpar == "expectedUpperLimit" and float(vals[allpar.index("upperLimit")]) == 0.0:
                    print "INFO: %s removing obsUL=0.0 for %s" % (ana, key)
                    continue;

                # throw away points with CLsexp > 0.99 and UL < 1.0 and CLs=-1 and UL<1 when merging on UL                  
                CLsExp = float( vals[allpar.index("CLsexp/F")])
                if selectpar == "expectedUpperLimit" and pval < 1.0 and (CLsExp>0.99 or CLsExp<0) and float( vals[allpar.index("upperLimit")])<1:                   
                    if CLsExp>0.99: print "INFO: %s replacing CLsexp with 0.01 since UL < 1.0  and CLsexp=1 for %s" % (ana, key)
                    elif CLsExp<0: print "INFO: %s replacing CLsexp with 0.01 since UL < 1.0  and CLsexp=-1 for %s" % (ana, key)
                    vals[allpar.index("CLsexp/F")] = str(0.01)
                    vals[allpar.index("CLs/F")] = str(0.01)
                    vals[allpar.index("clsu1s/F")] = str(0.01)
                    vals[allpar.index("clsd1s/F")] = str(0.01)
                    vals[allpar.index("p1/F")] = str(0.01)

                    pass;
                 
                newline=""
                #print vals;
                for val in vals:
                    newline += (str)(val)+" "
                    pass;
                newline = newline.rstrip(" ")
                newline += "\n"

                # float strings -> so make them float, then an int to throw away .0000 and then to bool
                failedcov =  bool(int(float(vals[allpar.index("failedcov/F")])))  # Mediocre cov matrix quality
                covqual = int(float(vals[allpar.index("covqual/F")]))             # covqual
                failedfit = bool(int(float(vals[allpar.index("failedfit/F")])))   # Fit failure
                failedp0 = bool(int(float(vals[allpar.index("failedp0/F")])))     # Base p0 ~ 0.5 (this can reject good fits)!
                fitstatus = bool(int(float(vals[allpar.index("fitstatus/F")])))   # Fit status from Minuit
                nofit = bool(int(float(vals[allpar.index("nofit/F")])))           # Whether there's a fit present
                # ignore some checker
                """
                # ignore failed fit
                if failedfit:
                    print "INFO: %s removing failedfit=true for %s" % (ana, key)
                    continue
                # ignore bad mediocre cov matrix quality
                if failedcov:
                    print "INFO: %s removing failedcov=true for %s" % (ana, key)
                # ignore if covqual<3 & covqual!=-1
                if covqual < 3 and covqual != -1:
                    print "INFO: %s removing check if (covqual<3 and covqual!=-1) for %s (found covqual=%d)" % (ana, key, covqual)
                    continue
                """
    
                key = (par1,par2)
                if config.grid.find("GG_onestepCC")!=-1: 
                    key = (par1,par2,par3)
                    pass;
                if config.grid.find("SS_onestepCC")!=-1: 
                    key = (par1,par2,par3)
                    pass;
                if config.grid.find("SM_GG_N2")!=-1: 
                    key = (par1,par2,par3)
                    pass;
                if config.grid.find("GG_onestep_mlsp60")!=-1: 
                    key = (par1,par2,par3)
                    pass;
                
                if key not in myMap.keys():
                    myMap[key] = [pval,newline]
                else:
                    if pval < myMap[key][0] and pval>=0:
                        print "DEBUG: %s found new best value - selectpar=%.2e < previous=%e for %s" % (ana, pval, myMap[key][0], key)
                        myMap[key][0] = pval
                        myMap[key][1] = newline
                        pass;
                    pass;
                pass;
            infile.close()

        #print myMap
        combined_filename = config.outputDir+config.outputName+"_combined_fixSigXSec"+xsecStr+listSuffix
        combinedfile = open(combined_filename,"w")
        for key,info in myMap.items():
            print key,": ",info[0],newline;
            combinedfile.write(info[1])
            pass;
        combinedfile.close()
        file_info.close()
        cmd="root -b -q \"$ZEROLEPTONFITTER/macros/contourplot/makecontourhists.C(\\\""+combined_filename+"\\\",\\\""+config.grid+"\\\")\""
        print cmd
        subprocess.call(cmd, shell=True)
        cmd="mv *_list.root "+config.outputDir
        subprocess.call(cmd, shell=True)
    pass;
# end of def Oring()



###########################################################################
#Main
###########################################################################

def main():
    config = parseCmdLine(sys.argv[1:])

    print "## Start makeContours_Run2.py ##"
    # default anaList
    anaList = [];
    anaList = finalChannelsDict.keys();
            
    # anaList=[] in doAll
    config.anaList=anaList
    print config.anaList
    
    # merge result rootfiles of HistFitter to Outputs/*fixSigXSec*.root in each SRs
    if config.doMerge or config.doAll:
        # merge file and set config.anaList ( in GG/SS_direct or GG_onestepCC,  anaList=[SR2jl/m/t,...] )
        print "## Start MergeFiles ##"
        MergeFiles(config)
        pass;
        
    if config.makeContours or config.doAll:
        print "## Start MakeContours ##"
        MakeContours(config)
        pass;

    if config.doOring or config.doAll:
        print "## Start Oring ##"
        Oring(config)
        pass;
  
if __name__ == "__main__":
    main()
