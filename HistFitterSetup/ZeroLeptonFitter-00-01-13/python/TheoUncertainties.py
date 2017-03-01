#######################################
# The uncertainties in this file are the
# uncertainties on the transfer factor
#######################################

def getError(channelName, regionName, sysDict):
    error = 0
    #print "toto"
    #print sysDict
    if (channelName,regionName) in sysDict.keys():
        error = sysDict[(channelName,regionName)]
    elif (channelName,"default") in sysDict.keys():
        error = sysDict[(channelName,"default")]
    elif ("default","default") in sysDict.keys():
        error = sysDict[("default","default")]
    return error


#######################################
#Diboson
#######################################
dibosonFlatSysDict = {}
dibosonFlatSysDict[("default","default")] = 0.2

#######################################
#Z+jets
#######################################
# Sherpa vs MadGraph
zTheoSysGeneratorDict = {}
zTheoSysGeneratorDict[("default","default")] = 0.11
zTheoSysGeneratorDict[("SRJigsawSRS1a","default")] = 0.0835866
zTheoSysGeneratorDict[("SRJigsawSRS1b","default")] = 0.0763378
zTheoSysGeneratorDict[("SRJigsawSRS2a","default")] = 0.0689389
zTheoSysGeneratorDict[("SRJigsawSRS2b","default")] = 0.0756263
zTheoSysGeneratorDict[("SRJigsawSRS3a","default")] = 0.0664213
zTheoSysGeneratorDict[("SRJigsawSRS3b","default")] = 0.0645755
zTheoSysGeneratorDict[("SRJigsawSRS4","default")] = 0.0597674
zTheoSysGeneratorDict[("SRJigsawSRG1a","default")] = 0.0813217
zTheoSysGeneratorDict[("SRJigsawSRG1b","default")] = 0.0837482
zTheoSysGeneratorDict[("SRJigsawSRG2a","default")] = 0.0672614
zTheoSysGeneratorDict[("SRJigsawSRG2b","default")] = 0.0998166
zTheoSysGeneratorDict[("SRJigsawSRG3a","default")] = 0.0828165
zTheoSysGeneratorDict[("SRJigsawSRG3b","default")] = 0.208771
zTheoSysGeneratorDict[("SRJigsawSRG4","default")] = 0.141553
zTheoSysGeneratorDict[("SRJigsawSRC1","default")] = 0.0801382
zTheoSysGeneratorDict[("SRJigsawSRC2","default")] = 0.0832055
zTheoSysGeneratorDict[("SRJigsawSRC3","default")] = 0.11795
zTheoSysGeneratorDict[("SRJigsawSRC4","default")] = 0.111056
zTheoSysGeneratorDict[("SRJigsawSRC5","default")] = 0.164243

#######################################
#W+jets
#######################################
# SHerpa vs MadGraph
wTheoSysGeneratorDict = {}
wTheoSysGeneratorDict[("default","default")] = 0.10
wTheoSysGeneratorDict[("SRJigsawSRS1a","default")] = 0.110786
wTheoSysGeneratorDict[("SRJigsawSRS1b","default")] = 0.115062
wTheoSysGeneratorDict[("SRJigsawSRS2a","default")] = 0.131817
wTheoSysGeneratorDict[("SRJigsawSRS2b","default")] = 0.13979
wTheoSysGeneratorDict[("SRJigsawSRS3a","default")] = 0.120151
wTheoSysGeneratorDict[("SRJigsawSRS3b","default")] = 0.124532
wTheoSysGeneratorDict[("SRJigsawSRS4","default")] = 0.139939
wTheoSysGeneratorDict[("SRJigsawSRG1a","default")] = 0.111579
wTheoSysGeneratorDict[("SRJigsawSRG1b","default")] = 0.0968154
wTheoSysGeneratorDict[("SRJigsawSRG2a","default")] = 0.0737987
wTheoSysGeneratorDict[("SRJigsawSRG2b","default")] = 0.0954486
wTheoSysGeneratorDict[("SRJigsawSRG3a","default")] = 0.116864
wTheoSysGeneratorDict[("SRJigsawSRG3b","default")] = 0.358518
wTheoSysGeneratorDict[("SRJigsawSRG4","default")] = 0.268703
wTheoSysGeneratorDict[("SRJigsawSRC1","default")] = 0.168484
wTheoSysGeneratorDict[("SRJigsawSRC2","default")] = 0.166257
wTheoSysGeneratorDict[("SRJigsawSRC3","default")] = 0.100981
wTheoSysGeneratorDict[("SRJigsawSRC4","default")] = 0.118663
wTheoSysGeneratorDict[("SRJigsawSRC5","default")] = 0.139578


#######################################
#Top
#######################################
# PowhegPythia vs aMcAtNloHerwigppE
topTheoSysGeneratorDict = {}
topTheoSysGeneratorDict[("default","default")] = 0.10
topTheoSysGeneratorDict[("SR2jl","SR")]= -0.0270272272382
topTheoSysGeneratorDict[("SR2jm","SR")]= 0.349636497148
topTheoSysGeneratorDict[("SR2jt","SR")]= 0.461538513093
topTheoSysGeneratorDict[("SR4jt","SR")]= 1.0
topTheoSysGeneratorDict[("SR5j","SR")]= -0.762753262497
topTheoSysGeneratorDict[("SR6jm","SR")]= -0.517948585602
topTheoSysGeneratorDict[("SR6jt","SR")]= -1.55128217279


#Additionnal radiation
topTheoSysRadDict= {}
topTheoSysRadDict[("default","default")] = (0.1,0.1)
topTheoSysRadDict[("SR2jl","SR")]=( 0.0142020364542 , -0.201620192291 )
topTheoSysRadDict[("SR2jm","SR")]=( 0.186026896402 , -0.245470948202 )
topTheoSysRadDict[("SR2jt","SR")]=( 0.125000113877 , 0.0569949678427 )
topTheoSysRadDict[("SR4jt","SR")]=( 0.346846851333 , -0.237373681417 )
topTheoSysRadDict[("SR5j","SR")]=( 0.188578327358 , -0.188673048593 )
topTheoSysRadDict[("SR6jm","SR")]=( 0.106187637855 , -0.453252653384 )
topTheoSysRadDict[("SR6jt","SR")]=( -0.479054176483 , -0.353741257356 )


#PowhegPythia vs PowhegHerwig
topTheoSysPowhegHerwigDict= {}
topTheoSysPowhegHerwigDict[("default","default")] = 0.1
topTheoSysPowhegHerwigDict[("SR2jl","SR")]= 0.00132307230229
topTheoSysPowhegHerwigDict[("SR2jm","SR")]= -0.266853972654
topTheoSysPowhegHerwigDict[("SR2jt","SR")]= -0.105263077168
topTheoSysPowhegHerwigDict[("SR4jt","SR")]= 0.664750955402
topTheoSysPowhegHerwigDict[("SR5j","SR")]= 0.028225738458
topTheoSysPowhegHerwigDict[("SR6jm","SR")]= 0.244110848105
topTheoSysPowhegHerwigDict[("SR6jt","SR")]= 0.416849812194


# Tune P2012 vs A14
topTheoSysA14Dict= {}
topTheoSysA14Dict[("default","default")] = 0.1
topTheoSysA14Dict[("SR2jl","SR")]= -0.0729149044195
topTheoSysA14Dict[("SR2jm","SR")]= -0.126971311478
topTheoSysA14Dict[("SR2jt","SR")]= -0.10526311211
topTheoSysA14Dict[("SR4jt","SR")]= 0.151515168029
topTheoSysA14Dict[("SR5j","SR")]= 0.250880763878
topTheoSysA14Dict[("SR6jm","SR")]= 0.209849710805
topTheoSysA14Dict[("SR6jt","SR")]= 0.0245097704868


