import ROOT
import logging
import shutil
import os


def discover(sh, search_directories, pattern="*"):

	# scan for datasets in the given directories
	for directory in search_directories:
	    print directory
	    ROOT.SH.ScanDir().samplePattern(pattern).scan(sh, directory)


	logging.info("%d different datasets found scanning all directories", len(sh))

	return sh



def addTags(sh_all):
	for sample in sh_all:
		short_name = sample.getMetaString("sample_name").split(".")[5]
		dsid = sample.getMetaString("sample_name").split(".")[4]
		sample.setMetaString( "short_name" , sample.getMetaString("sample_name") )

		if ('truth' in sample.getMetaString("sample_name")) :
			sample.addTag("truth")
		else :
			sample.addTag("reco")

		if ('truth' in sample.getMetaString("sample_name")) :
			if(('LO') in sample.getMetaString("sample_name") or
			   "1Gam" in short_name or
			   ('lo') in sample.getMetaString("sample_name") or
			   ('gamma') in sample.getMetaString("sample_name")
			   ) :
				sample.addTag("lo")
			else : sample.addTag('nlo')

		print short_name

		if "physics_" in short_name:
			sample.addTag("data")

		if "GG_direct" in short_name:
			sample.addTag("signal")
			sample.addTag("gg_direct")
		if "SS_direct" in short_name:
			sample.addTag("signal")
			sample.addTag("ss_direct")

		if "ttbar" in short_name:
			sample.addTag("top")
			sample.addTag("ttbar")
		if "Wt_inclusive" in short_name:
			sample.addTag("top")
			sample.addTag("singletop")
		if "PowHPEvG_singletop" in short_name:
			sample.addTag("top")
			sample.addTag("singletop")

		if "jetjet" in short_name:
			sample.addTag("qcd")

		if "1Gam" in short_name:
			sample.addTag("gamma")

		if "PowHP8EvG_W" in short_name:
			sample.addTag("wjets")
		if "PowHP8EvG_Z" in short_name:
			sample.addTag("zjets")
		if "Znunu" in short_name:
			sample.addTag("zjets")
		if "Sherpa_Wqq" in short_name:
			sample.addTag("zjets")
		if "_Zmumu_" in short_name:
			sample.addTag("zjets")
		if "ggH" in short_name:
			sample.addTag("signal")
