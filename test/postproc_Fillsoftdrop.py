#!/usr/bin/env python
import os, sys
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODJMARTools.postprocessing.modules.jme.FillSDhistosForUnfolding import *

redir = 'root://cmsxrootd.fnal.gov/'

files=[#redir +
    #"/store/user/asparker/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZandJetSkimNANAODreclusterDY1JetsToLLM-50TuneCP513TeV-madgraphMLM-pythia8/180207_192924/0000/94XNanoV0-DYtoLL-nanoTrees_96.root",
   redir + "/store/user/asparker/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZplusJetSelection94XNANAODreclusterDY1JetsToLLM-50TuneCP513TeV-madgraphMLM-pythia8/180212_155223/0000/94XNanoV0-DYtoLL-nanoTrees_7.root",
   redir + "/store/user/asparker/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZplusJetSelection94XNANAODreclusterDY1JetsToLLM-50TuneCP513TeV-madgraphMLM-pythia8/180212_155223/0000/94XNanoV0-DYtoLL-nanoTrees_8.root",
   redir + "/store/user/asparker/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZplusJetSelection94XNANAODreclusterDY1JetsToLLM-50TuneCP513TeV-madgraphMLM-pythia8/180212_155223/0000/94XNanoV0-DYtoLL-nanoTrees_9.root",
   redir + "/store/user/asparker/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZplusJetSelection94XNANAODreclusterDY1JetsToLLM-50TuneCP513TeV-madgraphMLM-pythia8/180212_155223/0000/94XNanoV0-DYtoLL-nanoTrees_10.root",
   redir + "/store/user/asparker/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZplusJetSelection94XNANAODreclusterDY1JetsToLLM-50TuneCP513TeV-madgraphMLM-pythia8/180212_155223/0000/94XNanoV0-DYtoLL-nanoTrees_11.root",
   redir + "/store/user/asparker/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/ZplusJetSelection94XNANAODreclusterDY1JetsToLLM-50TuneCP513TeV-madgraphMLM-pythia8/180212_155223/0000/94XNanoV0-DYtoLL-nanoTrees_12.root",




   ]

import random
random.seed(12345)

p1=PostProcessor(".",files,'','',[sdb0()],provenance=False, noOut=True, histFileName='sdb0.root', histDirName='sdb0', postfix='plots')

p1.run()
