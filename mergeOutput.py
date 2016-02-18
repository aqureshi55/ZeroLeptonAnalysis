#!/usr/bin/env python

########### Initialization ######################################
##
##

import ROOT
import logging
import shutil
import os, sys
import itertools
import array

import discoverInput

logging.basicConfig(level=logging.INFO)
from optparse import OptionParser

import atexit
@atexit.register
def quite_exit():
	ROOT.gSystem.Exit(0)


logging.info("loading packages...")
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")

##
##
########### Configuration ######################################
##
##

lumi = 1000  ## in pb-1
if len(sys.argv)==1:
	search_directories = ["/afs/cern.ch/user/l/leejr/work/public/fromGrid/",]
	if 'bnl' in os.getenv('HOSTNAME'):
		search_directories = ["/pnfs/usatlas.bnl.gov/users/russsmith/photonTruthStudies_zG_test/"]
else:
	search_directories = [str(sys.argv[i]) for i in xrange(1,len(sys.argv))]

print search_directories

tmpOutputDirectory = "tmpOutput"
outputDirectory = "output"
#treePrefix = ""

selection = "1"

def main():

	##
	##
	########### Gather input ######################################
	##
	##

	logging.info("creating new sample handler")
	sh_all = ROOT.SH.SampleHandler()

	discoverInput.discover(sh_all, search_directories)
	print len(sh_all)

	logging.info("adding my tags defined in discoverInput.py")
	discoverInput.addTags(sh_all)

	ROOT.SH.readSusyMetaDir(sh_all,"$ROOTCOREBIN/data/SUSYTools")

	## Split up samplehandler into per-BG SH's based on tag metadata

	sh_data = sh_all.find("data")
	sh_signal = sh_all.find("signal")
	sh_bg = {}

	sh_bg["qcd"] = sh_all.find("qcd")
	sh_bg["top"] = sh_all.find("top")
	sh_bg["wjets"] = sh_all.find("wjets")
	sh_bg["zjets"] = sh_all.find("zjets")
	sh_bg["gamma"] = sh_all.find("gamma")

	######This will be done per samplehandler ############################
	##
	##

	#Creation of output directory names
	outputFileNames = {
		sh_all: "All",
		sh_data: "Data",
		sh_signal: "Signal",
		sh_bg["qcd"]: "QCD",
		sh_bg["top"]: "Top",
		sh_bg["wjets"]: "WJets",
		sh_bg["zjets"]: "ZJets",
		sh_bg["gamma"]: "Gamma",
		}

	treesToProcess = set()

	for mysamplehandler in [
								# sh_all,
								# sh_data,
								# sh_bg["qcd"],
								# sh_bg["top"],
								# sh_bg["wjets"],
#								sh_bg["zjets"]
								sh_bg["gamma"]
							]:

		# print mysamplehandler

		filesToEventuallyHadd = []

		for sample in mysamplehandler:

			sample_name = sample.getMetaString("sample_name")
			# print sample_name
			dsid = sample_name.split(".")[2]

			if not treesToProcess :
				treesToProcess = getListOfTreeNames(sample)
			print treesToProcess

			attachCounters(sample)

			mydir = tmpOutputDirectory

			try:
				os.stat(mydir)
			except:
				os.mkdir(mydir)

			outputSampleFileName = "%s/%s.root"%(tmpOutputDirectory, dsid)
			filesToEventuallyHadd.append(outputSampleFileName)

			outputSampleFile = ROOT.TFile(outputSampleFileName,"RECREATE")
			print outputSampleFile
			print treesToProcess

			for itree in treesToProcess:
				mysamplehandler.setMetaString("nc_tree", itree)
				print itree
				mytree =  sample.makeTChain().Clone()
#				mytree.Print()
				outputTree = ROOT.addBranch( mytree, getNormFactor(sample) , selection)
				outputTree.Write()

			print "Saved tree %s with %s events . . ." % ( outputTree.GetName(), outputTree.GetEntries() )

			outputSampleFile.Close()


		mydir = outputDirectory
		try:
			os.stat(mydir)
		except:
			os.mkdir(mydir)
		print 'hadd -O -f %s/%s.root %s'%(outputDirectory, outputFileNames[mysamplehandler], " ".join(filesToEventuallyHadd) )
		os.system('hadd -O -f %s/%s.root %s'%
			(outputDirectory, outputFileNames[mysamplehandler], " ".join(filesToEventuallyHadd) )
			)



#To scale the histograms in the files after the event loop is done...
def getNormFactor(sample):

	tempxs = sample.getMetaDouble("nc_xs") * sample.getMetaDouble("kfactor") * sample.getMetaDouble("filter_efficiency")

	print "Norm weight for %s is %f/(%f or %f)"%(sample.getMetaString("short_name"), tempxs, sample.getMetaDouble("nc_nevt"), sample.getMetaDouble("nc_sumw"))
	m_eventscaling = tempxs
	if sample.getMetaDouble("nc_nevt"):
		m_eventscaling /= sample.getMetaDouble("nc_nevt") if "jetjet" in sample.getMetaString("short_name") else sample.getMetaDouble("nc_sumw")
	else:
		m_eventscaling = 0.
	return m_eventscaling


addBranchCode = """
TTree * addBranch(TTree* tree, float normalization, TString selection="1"){

                std::cout << "copying tree" << std::endl;
                std::cout << tree->GetEntries() << std::endl;

		TTree * newtree = tree->CopyTree(selection);
		int nevents = newtree->GetEntries();

                std::cout << "got entries" << std::endl;

		double normweight = normalization;
                std::cout << "normweight : " << normweight << std::endl;
		TBranch * bnormweight = newtree->Branch("normweight",&normweight,"normweight/D");

                std::cout << "looping entries" << std::endl;

		for (Long64_t i=0;i<nevents;i++) {
			newtree->GetEntry(i);
			bnormweight->Fill();
			if(i%10000==0) cout<< i << " of " << nevents << endl;
		}
                std::cout << "returning new tree" << std::endl;

		return newtree;
}
"""

ROOT.gInterpreter.Declare(addBranchCode)

# This python function is replaced by the c++ function above for speed
# def addLeaves(tree,normalization,selection="1"):
# 	return 0
# 	events = tree.GetEntries()
# 	leaves = "normweight/D"
# 	leafValues = array.array("d", [0.])
# 	# newtree = tree.CloneTree(0)
# 	newtree = tree.CopyTree(selection)
# 	newBranch = newtree.Branch( "normweight" , leafValues, leaves )
# 	for i in xrange(events):
# 		tree.GetEntry(i)
# 		leafValues[0] = ROOT.Double(normalization)
# 		newtree.Fill()
# 		if i % 10000 == 0:
# 			print "%s of %s: %s" % (i,events,leafValues)
# 	return newtree


def attachCounters(sample):

		m_nevt = 0
		m_sumw = 0

		#Go to the grid and get the metadata output
		sh_metadata = ROOT.SH.SampleHandler()
		discoverInput.discover(sh_metadata, search_directories, sample.getMetaString("sample_name").replace("trees","metadata") )
		assert len(sh_metadata) == 1
		metadata_sample = sh_metadata[0]
		for myfile in [ROOT.TFile(ifilepath) for ifilepath in metadata_sample.makeFileList() ]:
			# print myfile
			try:
				m_nevt += myfile.Get("Counter_JobBookeeping_JobBookeeping").GetBinContent(1)
				m_sumw += myfile.Get("Counter_JobBookeeping_JobBookeeping").GetBinContent(2)
			except:
				pass

		sample.setMetaDouble("nc_nevt",m_nevt)
		sample.setMetaDouble("nc_sumw",m_sumw)



def getListOfTreeNames(sample):
	f = ROOT.TFile(sample.fileName(0) )
	print f
	listOfTrees = set([key.GetName() for key in f.GetListOfKeys() if ('TTree' in f.Get(key.GetName()).ClassName())])
	print listOfTrees
	return listOfTrees




if __name__ == "__main__":
	main()

