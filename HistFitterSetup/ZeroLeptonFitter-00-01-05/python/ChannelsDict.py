from ChannelConfig import * 

###########################################################
# definition of the control regions
###########################################################

cleaningCut="( abs(jetEta[0])>2.4 || jet1Chf/jetFracSamplingMax[0]>0.1)"

regionDict={}
regionDict["SR"] = Region("SR", "SRAll", [cleaningCut], [])

regionDict["CRW"] = Region("CRW", "CRWT", ["nBJet==0"], ["bTagWeight"])
regionDict["CRT"] = Region("CRT", "CRWT", ["nBJet>0"], ["bTagWeight"])
regionDict["CRWT"] = Region("CRWT", "CRWT", ["nBJet>=0"], ["bTagWeight"])


regionDict["VRWf"] = Region("VRWf", "CRWT", ["nBJet==0"], ["bTagWeight"]) #ATT: systWeights[0] is a proxy for the lepton weight
regionDict["VRTf"] = Region("VRTf", "CRWT", ["nBJet>0"], ["bTagWeight"])

regionDict["VRWMf"] = Region("VRWMf", "VRWT", ["nBJet==0"], ["bTagWeight"])
regionDict["VRTMf"] = Region("VRTMf", "VRWT", ["nBJet>0"], ["bTagWeight"])
regionDict["VRWM"] = Region("VRWM", "VRWT", ["nBJet==0"], ["bTagWeight"])
regionDict["VRTM"] = Region("VRTM", "VRWT", ["nBJet>0"], ["bTagWeight"])

regionDict["CRZ"] = Region("CRZ", "CRZ", [], [])

regionDict["VRWTplus"] = Region("VRWTplus", "CRWT", ["lep1sign>0"], ["bTagWeight"])
regionDict["VRWTminus"] = Region("VRWTminus", "CRWT", ["lep1sign<0"], ["bTagWeight"])
regionDict["VRWTfplus"] = Region("VRWTfplus", "CRWT", ["lep1sign>0"], ["bTagWeight"])
regionDict["VRWTfminus"] = Region("VRWTfminus", "CRWT", ["lep1sign<0"], ["bTagWeight"])

regionDict["CRY"] = Region("CRY", "CRY", ["(phSignal==1 && phPt>130.)"], [])#extra weights should be applied only to gamma+jets
regionDict["VRYf"] = Region("VRY", "CRY", ["(phSignal==1 && phPt>130.)"], [])#extra weights should be applied only to gamma+jets

regionDict["CRQ"] = Region("CRQ", "SRAll", [cleaningCut])#ATT: qcd weight
regionDict["VRQ1"] = Region("VRQ1", "SRAll", [cleaningCut])#ATT: qcd weight
regionDict["VRQ2"] = Region("VRQ2", "SRAll", [cleaningCut])#ATT: qcd weight
regionDict["VRQ3"] = Region("VRQ3", "SRAll", [cleaningCut])#ATT: qcd weight
regionDict["VRQ4"] = Region("VRQ4", "SRAll", [cleaningCut])#ATT: qcd weight

regionDict["VRZ"] = Region("VRZ", "CRZ", [], [])
regionDict["VRZf"] = Region("VRZf", "CRZ", [], [])

regionDict["VRT2L"] = Region("VRT2L", "CRZ_VR1b", ["(mll>116000 &&  lep1Pt<200000 && lep2Pt<100000)"], [])#ATT: qcd weight


##for data-driven BG estimation##
regionDict["VRW"] = Region("VRW", "CRWT", ["nBJet==0"], ["bTagWeight"])
regionDict["CRWL"] = Region("VRWL", "CRWT", ["nBJet==0"], ["bTagWeight"])
regionDict["CRWVL"] = Region("VRWVL", "CRWT", ["nBJet==0"], ["bTagWeight"])
regionDict["VRT"] = Region("VRT", "CRWT", ["nBJet>0"], ["bTagWeight"])
regionDict["CRTL"] = Region("VRTL", "CRWT", ["nBJet>0"], ["bTagWeight"])
regionDict["CRTVL"] = Region("VRTVL", "CRWT", ["nBJet>0"], ["bTagWeight" ])
regionDict["CRYL"] = Region("CRYL", "CRY", ["(phSignal==1 && phPt>130."], [])#extra weights should be applied only to gamma+jets
regionDict["CRZL"] = Region("CRZL", "CRZ", [], [])
regionDict["CRZVL"] = Region("CRZVL", "CRZ", [], [])
regionDict["SRZVL"] = Region("SRZVL", "SRAll", ["( abs(jetEta[0])>2.4 || jet1Chf/jetFracSamplingMax[0]>0.1)"], [])



###########################################################
# definition of the final analysis channels
###########################################################

# FAR SRs
finalChannelsDict = {}

#SS_direct_800_400 2jl
# SR2jbase-MeSig15-Meff1200-sljetpt200-dphi0.8-ap0.00
anaSRFAR = ChannelConfig(name="SR2jl", regionDict=regionDict, fullname="SR2jbase-MeSig15-Meff1200-sljetpt200-dphi0.8-ap0.00")
anaSRFAR.nJets = 2
anaSRFAR.dPhi = 0.8
anaSRFAR.MET_over_meffNj = 0.0
anaSRFAR.METsig = 15
anaSRFAR.jetpt1 = 200.
anaSRFAR.jetpt2 = 200.
anaSRFAR.meffIncl = 1200
finalChannelsDict[anaSRFAR.name] = anaSRFAR

#GG_direct_750_650, GG_onestep_825_785_745 2jm
# SR2jbase-MeSig20-Meff1800-sljetpt60-dphi0.4-ap0.00
anaSRFAR = ChannelConfig(name="SR2jm", regionDict=regionDict, fullname="SR2jbase-MeSig20-Meff1800-sljetpt50-dphi0.4-ap0.00")
anaSRFAR.nJets = 2
anaSRFAR.dPhi = 0.4
anaSRFAR.MET_over_meffNj = 0.0
anaSRFAR.METsig = 20
anaSRFAR.jetpt1 = 200.
anaSRFAR.jetpt2 = 50.
anaSRFAR.meffIncl = 1800
finalChannelsDict[anaSRFAR.name] = anaSRFAR

#SS_direct_1200_0 2jt
# SR2jbase-MeSig20-Meff2000-sljetpt200-dphi0.8-ap0.00
anaSRFAR = ChannelConfig(name="SR2jt", regionDict=regionDict, fullname="SR2jbase-MeSig20-Meff2000-sljetpt200-dphi0.8-ap0.00")
anaSRFAR.nJets = 2
anaSRFAR.dPhi = 0.8
anaSRFAR.MET_over_meffNj = 0.0
anaSRFAR.METsig = 20
anaSRFAR.jetpt1 = 200.
anaSRFAR.jetpt2 = 200.
anaSRFAR.meffIncl = 2000
finalChannelsDict[anaSRFAR.name] = anaSRFAR

# GG_direct_1400_0 4jt
anaSRFAR = ChannelConfig(name="SR4jt", regionDict=regionDict, fullname="SR4jbase-MetoMeff0.2-Meff2200-sljetpt100-34jetpt100-dphi0.4-ap0.04")
anaSRFAR.nJets = 4
anaSRFAR.dPhi = 0.4
anaSRFAR.dPhiR = 0.2
anaSRFAR.MET_over_meffNj = 0.2
anaSRFAR.MET = 200
anaSRFAR.METsig = 0
anaSRFAR.Ap = 0.04
anaSRFAR.jetpt1 = 200.
anaSRFAR.jetpt2 = 100.
anaSRFAR.jetpt3 = 100.
anaSRFAR.jetpt4 = 100.
anaSRFAR.meffIncl = 2200
finalChannelsDict[anaSRFAR.name] = anaSRFAR

#GG_onestepCC_1265_945_625 at 4/fb 5j
# SR5jbase-MetoMeff0.25-Meff1600-sljetpt100-34jetpt100-dphi0.4-ap0.02
anaSRFAR = ChannelConfig(name="SR5j", regionDict=regionDict, fullname="SR5jbase-MetoMeff0.25-Meff1600-sljetpt100-34jetpt100-dphi0.4-ap0.04")
anaSRFAR.nJets = 5
anaSRFAR.dPhi = 0.4
anaSRFAR.dPhiR = 0.2
anaSRFAR.MET_over_meffNj = 0.25
anaSRFAR.METsig = 0
anaSRFAR.Ap = 0.04
anaSRFAR.jetpt1 = 200.
anaSRFAR.jetpt2 = 100.
anaSRFAR.jetpt3 = 100.
anaSRFAR.jetpt4 = 100.
anaSRFAR.jetpt5 = 50.
anaSRFAR.meffIncl = 1600
finalChannelsDict[anaSRFAR.name] = anaSRFAR

#GG_onestepCC_1265_945_625 at 4/fb 6jm
# SR6jbase-MetoMeff0.25-Meff1600-sljetpt100-34jetpt100-dphi0.4-ap0.02
anaSRFAR = ChannelConfig(name="SR6jm", regionDict=regionDict, fullname="SR6jbase-MetoMeff0.25-Meff1600-sljetpt100-34jetpt100-dphi0.4-ap0.04")
anaSRFAR.nJets = 6
anaSRFAR.dPhi = 0.4
anaSRFAR.dPhiR = 0.2
anaSRFAR.MET_over_meffNj = 0.25
anaSRFAR.METsig = 0
anaSRFAR.Ap = 0.04
anaSRFAR.jetpt1 = 200.
anaSRFAR.jetpt2 = 100.
anaSRFAR.jetpt3 = 100.
anaSRFAR.jetpt4 = 100.
anaSRFAR.jetpt5 = 50.
anaSRFAR.jetpt6 = 50.
anaSRFAR.meffIncl = 1600
finalChannelsDict[anaSRFAR.name] = anaSRFAR

#GG_onestepCC_1545_785_25 at 4/fb
# SR6jbase-MetoMeff0.2-Meff1800-sljetpt100-34jetpt100-dphi0.4-ap0.02
anaSRFAR = ChannelConfig(name="SR6jt", regionDict=regionDict, fullname="SR6jbase-MetoMeff0.2-Meff1800-sljetpt100-34jetpt100-dphi0.4-ap0.04")
anaSRFAR.nJets = 6
anaSRFAR.dPhi = 0.4
anaSRFAR.dPhiR = 0.2
anaSRFAR.MET_over_meffNj = 0.2
anaSRFAR.METsig = 0
anaSRFAR.Ap = 0.04
anaSRFAR.jetpt1 = 200.
anaSRFAR.jetpt2 = 100.
anaSRFAR.jetpt3 = 100.
anaSRFAR.jetpt4 = 100.
anaSRFAR.jetpt5 = 50.
anaSRFAR.jetpt6 = 50.
anaSRFAR.meffIncl = 2000
finalChannelsDict[anaSRFAR.name] = anaSRFAR


###########################################################
# definition of the channels for plotting purpose
###########################################################

channelsForPlottingDict = {}

#SRs used for plotting
anaSR2jrun1 = ChannelConfig(name="SR2jrun1", regionDict=regionDict)
anaSR2jrun1.nJets = 2
anaSR2jrun1.jetpt1 = 130
anaSR2jrun1.jetpt2 = 60
anaSR2jrun1.MET = 160
anaSR2jrun1.dPhi = 0.4
#anaSR2jrun1.MET_over_meff  =  0.1
channelsForPlottingDict[anaSR2jrun1.name] = anaSR2jrun1

anaSR3jrun1 = ChannelConfig(name="SR3jrun1", regionDict=regionDict)
anaSR3jrun1.nJets = 3
anaSR3jrun1.jetpt1 = 130
anaSR3jrun1.jetpt2 = 60
anaSR3jrun1.jetpt3 = 60
anaSR3jrun1.MET = 160
anaSR3jrun1.dPhi = 0.4
#anaSR3jrun1.MET_over_meff  =  0.1
channelsForPlottingDict[anaSR3jrun1.name] = anaSR3jrun1

anaSR4jrun1 = ChannelConfig(name="SR4jrun1", regionDict=regionDict)
anaSR4jrun1.nJets = 4
anaSR4jrun1.jetpt1 = 130
anaSR4jrun1.jetpt2 = 60
anaSR4jrun1.jetpt3 = 60
anaSR4jrun1.jetpt4 = 60
anaSR4jrun1.MET = 160
anaSR4jrun1.dPhi = 0.4
#anaSR4jrun1.MET_over_meff  =  0.1
channelsForPlottingDict[anaSR4jrun1.name] = anaSR4jrun1

anaSR2jEPS = ChannelConfig(name="SR2jEPS", regionDict=regionDict)
anaSR2jEPS.nJets = 2
anaSR2jEPS.jetpt1 = 100
anaSR2jEPS.jetpt2 = 60
anaSR2jEPS.MET = 100
#anaSR2jEPS.dPhi = 0.4
channelsForPlottingDict[anaSR2jEPS.name] = anaSR2jEPS


###########################################################
# special channel to compute kappa
# very loose selection to compute kappa correction
# CRT and CRW should be add as constraining regions
# CRZ and VRZ should be add as validation regions
###########################################################
#finalChannelsDict={}
anaVL=ChannelConfig(name="VLForKappa",regionDict=regionDict)
anaVL.nJets=2
anaVL.dPhi=0
anaVL.MET=200    
anaVL.MET_upper=300  
anaVL.jetpt1 = 200.
#finalChannelsDict[anaVL.name]=anaVL

###########################################################
# SR for testing the fit
# Upper cut on meff are applied
###########################################################

testChannelsDict={}

anaSR = ChannelConfig(name="SR2jTest", regionDict=regionDict)
anaSR.nJets = 2
anaSR.dPhi = 0.4
anaSR.METsig = 10
anaSR.jetpt1 = 200.
anaSR.meffIncl = 800
anaSR.meffInclUpperCut = 1000
testChannelsDict[anaSR.name] = anaSR


anaSR = ChannelConfig(name="SR3jTest", regionDict=regionDict)
anaSR.nJets = 3
anaSR.dPhi = 0.4
anaSR.METsig = 10
anaSR.jetpt1 = 200.
anaSR.meffIncl = 800
anaSR.meffInclUpperCut = 1000
testChannelsDict[anaSR.name] = anaSR


anaSR = ChannelConfig(name="SR4jTest", regionDict=regionDict)
anaSR.nJets = 4
anaSR.dPhi = 0.4
anaSR.dPhiR = 0.2
anaSR.METsig = 10
#anaSR.MET_over_meffNj = 0.15
anaSR.jetpt1 = 200.
anaSR.meffIncl = 800
anaSR.meffInclUpperCut = 1000
testChannelsDict[anaSR.name] = anaSR


anaSR = ChannelConfig(name="SR5jTest", regionDict=regionDict)
anaSR.nJets = 5
anaSR.dPhi = 0.4
anaSR.dPhiR = 0.2
anaSR.METsig = 10
#anaSR.MET_over_meffNj = 0.2
anaSR.jetpt1 = 200.
anaSR.meffIncl = 800
anaSR.meffInclUpperCut = 1000
testChannelsDict[anaSR.name] = anaSR



anaSR = ChannelConfig(name="SR6jTest", regionDict=regionDict)
anaSR.nJets = 6
anaSR.dPhi = 0.4
anaSR.dPhiR = 0.2
anaSR.METsig = 10
#anaSR.MET_over_meffNj = 0.2
anaSR.jetpt1 = 200.
anaSR.meffIncl = 800
anaSR.meffInclUpperCut = 1000
testChannelsDict[anaSR.name] = anaSR



#This is temporary and will be removed with the actual SR to find SUSY
#to be removed
#finalChannelsDict=testChannelsDict
#finalChannelsDict.update(testChannelsDict)





###########################################################
# all channels
###########################################################

allChannelsDict = finalChannelsDict.copy()
allChannelsDict.update(channelsForPlottingDict)




