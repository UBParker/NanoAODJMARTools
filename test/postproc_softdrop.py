#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODJMARTools.postprocessing.modules.jme.softdrop import *

redir = 'root://cmsxrootd.fnal.gov/'

files=[#redir +
    #"/store/user/asparker/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZandJetSkimNANAODreclusterDY1JetsToLLM-50TuneCP513TeV-madgraphMLM-pythia8/180207_192924/0000/94XNanoV0-DYtoLL-nanoTrees_96.root",
    "/uscms_data/d3/aparker/nanoAod/NanoAODJMAR_V0/CMSSW_9_4_4/src/PhysicsTools/NanoAODJMAR/test/94X_JMARNANOSkim3/ZandJetSkimofNANAODreclusterDY1JetsToLLM-50-trees-90per.root"
   ]

import random
random.seed(12345)

p1=PostProcessor(".",files,'','',[sdb0()],provenance=False, noOut=True, histFileName='sdb0.root', histDirName='sdb0', postfix='plots')

p1.run()
