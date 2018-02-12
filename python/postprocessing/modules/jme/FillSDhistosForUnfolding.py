import ROOT
import math, os
import numpy as np
import array as array
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import matchFastJetObjectCollection


class softDropProducer(Module):
    def __init__(self, beta=0.0, zcut=0.1, bname="sdb0", R=0.8, ptmin=200.):
        self.writeHistFile=True
        self.verbose =  True
        
        self.beta = beta
        self.zcut = zcut
        self.R = R
        self.ptmin = ptmin

        self.jetBranchName = "FatJet"
        self.genJetBranchName = "GenJetAK8"
        self.rhoBranchName = "fixedGridRhoFastjetAll"
        self.pfCandsBranchName = "PFCandsAK8"
        self.genCandsBranchName = "GenPartAK8"
        self.bname = bname
        
        print "Load C++ Recluster worker module"                
        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")
                
        

    def beginJob(self, histFile, histDirName):
        self.sd = ROOT.SoftDropWrapper(self.beta,self.zcut, self.R, self.ptmin)
        Module.beginJob(self, histFile, histDirName)
        self.addObject(ROOT.TH1F('h_ak8sdm_'+self.bname,   'h_ak8sdm_'+self.bname, 25, 0, 250) )
        self.binsGen = array.array('d', [0., 1., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 100.]) 
        self.nGen = len(self.binsGen) - 1
        self.binsDet = array.array('d', [0., 0.5, 1., 3., 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30., 32.5, 35., 37.5, 40., 42.5, 45., 47.5, 50., 75., 100.])
        self.nDet = len(self.binsDet) - 1
        self.nDetSD =  len(self.binsDet) - 1 #18
        self.nGenSD =  len(self.binsGen) - 1 #9
        
        self.resp0 =  ROOT.TH2D('h_response_'+self.bname,     'h_response_'+self.bname ,   self.nDetSD, self.binsDet, self.nGenSD, self.binsGen) 
        self.reco0 = ROOT.TH1D('h_reco_'+self.bname,         'h_reco_'+self.bname,       self.nDetSD, self.binsDet)
        self.gen0 = ROOT.TH1D('h_gen_'+self.bname,          'h_gen_'+self.bname,        self.nGenSD, self.binsGen)
        self.fake0 = ROOT.TH1D('h_fake_'+self.bname,         'h_fake_'+self.bname,       self.nDetSD, self.binsDet)
        self.miss0 =  ROOT.TH1D('h_miss_'+self.bname,         'h_miss_'+self.bname,       self.nGenSD, self.binsGen)
        
        
        self.addObject( self.resp0 )
        self.addObject( self.reco0 )
        self.addObject(  self.gen0 )
        self.addObject( self.fake0 )
        self.addObject( self.miss0 )
        

    def endJob(self):
        Module.endJob(self)
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        #jets = Collection(event, self.jetBranchName )
        #genJets = Collection(event, self.genJetBranchName )
        pfCands = Collection(event, self.pfCandsBranchName )
        genCands = Collection(event, self.genCandsBranchName )

        pfCandsVec = ROOT.vector("TLorentzVector")()
        for p in pfCands :
            pfCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
        sdjets = self.sd.result( pfCandsVec )
        #sdjets.sort(key=lambda x:x.p4().Perp(),reverse=True)

        genCandsVec = ROOT.vector("TLorentzVector")()
        for p in genCands :
            genCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
        gensdjets = self.sd.result( genCandsVec )
        #gensdjets.sort(key=lambda x:x.p4().Perp(),reverse=True)
        

        if len(pfCands) == 0 :
            return False
        # Read the flags written by event selector 
        # 
        if event.miss < 1 and event.reco < 1 and event.fake < 1 :
            return False
        print 'Event : ', event.event
        
               
        if not event.goodreco and not event.goodgen :
            return False
        if event.goodreco and event.goodgen :
            recoToGen = matchFastJetObjectCollection( sdjets, gensdjets, dRmax=0.05)
            for reco,gen in recoToGen.iteritems():
                if reco == None or reco.perp() < 200. or reco.m() < 1.:
                    continue
                if event.reco != 0 : 
                    self.reco0.Fill(reco.m())
                    if self.verbose : print "Filling reco histo with SD jet of mass {:3.0f} GeV and Pt of {:3.0f} GeV".format(reco.m(), reco.perp())
                if gen != None and gen.perp() > 200. * 0.8 and  gen.m() > 1. :
                    if event.response != 0 and event.gen != 0 : 
                        #print "event.response is {}".format(event.response) 
                        self.resp0.Fill(reco.m(), gen.m() )
                        self.gen0.Fill(gen.m()) 
                        if self.verbose : print "Filling response and gen histo with gen SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV".format(gen.m(), gen.perp())
                        #return True
                if event.fake != 0 and event.gen == 0 : 
                    self.fake0.Fill(reco.m())
                    if self.verbose : print "Filling fake histo with SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV ".format(reco.m(), reco.perp())
                    #return True
            for igen,gen in enumerate(gensdjets):
                if gen != None and gen not in recoToGen.values() :
                    if event.miss and gen.perp() > 200. and  gen.m() > 1. : 
                        self.gen0.Fill(gen.m())
                        self.resp0.Fill(-1., gen.m() )
                        self.miss0.Fill(gen.m())        
                        if self.verbose : print "Filling miss/gen/response histo with gen SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV".format(gen.m(), gen.perp())        
                        #return True
        elif (event.goodreco and not event.goodgen ) :
            #Fake
            if len(sdjets) < 1 : return False
            if event.fake == 0 : return False
            for reco in sdjets:

                if reco.perp() < 200. or  reco.m() < 1.: continue
                self.fake0.Fill(reco.m())
                self.reco0.Fill(reco.m())
                if self.verbose : print "Filling reco fake histo with SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV".format(reco.m(), reco.perp())
                #self.resp0.Fill( reco.p4().M(), -1. ) 
                #return True
        elif  (not event.goodreco and event.goodgen ) :
            #Miss
            if event.miss == 0 : return False
            if len(gensdjets) < 1 : return False
            for gen in gensdjets :
                if gen.perp() < 200.  or  gen.m() < 1. :continue 
                self.miss0.Fill(gen.m())
                self.gen0.Fill(gen.m())
                self.resp0.Fill( -1.0, gen.m() )
                if self.verbose : print "Filling miss/gen/response histo with gen SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV".format(gen.m(), gen.perp())
                #return True
        
        typeofill = ''
        if event.miss : typeofill = 'miss'
        if event.fake :typeofill = 'fake'
        if event.reco and not event.fake : typeofill = 'reco'
        if event.reco and event.gen : typeofill = 'gen'
        ## print "NANOAOD PF jets:"
        ## for i,jet in enumerate(jets):
        ##     print ' %5d: %6.2f %6.2f %6.2f %6.2f %6.2f' % (i, jet.pt, jet.eta, jet.phi, jet.mass, jet.msoftdrop)


        if typeofill == 'fake' or  typeofill == 'reco':
            print '-------'
            print "filling Soft Drop {} histo with reco jet below".format(typeofill)      
            print '-------'            
            #print "On-the-fly PF jets:"
            for i,sdjet in enumerate(sdjets):
                print ' %5d: %6.2f %6.2f %6.2f %6.2f' % (i, sdjet.perp(), sdjet.eta(), sdjet.phi(), sdjet.m())
                #print 'Subjets:'
                #subs = sdjet.pieces()
                #for j,sub in enumerate(subs):
                #    print '      : %6.2f %6.2f %6.2f %6.2f' % (sub.perp(), sub.eta(), sub.phi(), sub.m())

        if typeofill =='miss' or  typeofill == 'gen':
            print '-------'
            print "filling Soft Drop {} histo with gen jet below".format(typeofill)
            ## print "NANOAOD Gen jets:"
            ## for i,jet in enumerate(genJets):
            ##     print ' %5d: %6.2f %6.2f %6.2f %6.2f' % (i, jet.pt, jet.eta, jet.phi, jet.mass)
            print '-------'            
            #print "On-the-fly Gen jets:"
            for i,sdjet in enumerate(sdjets):
                print ' %5d: %6.2f %6.2f %6.2f %6.2f' % (i, sdjet.perp(), sdjet.eta(), sdjet.phi(), sdjet.m())
                #print 'Subjets:'
            #subs = sdjet.pieces()
            #for j,sub in enumerate(subs):
            #    print '      : %6.2f %6.2f %6.2f %6.2f' % (sub.perp(), sub.eta(), sub.phi(), sub.m())
            print '-------'
         
            
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

sdb0 = lambda : softDropProducer(beta=0.0, zcut=0.1, bname="sdb0")
sdb1 = lambda : softDropProducer(beta=1.0, zcut=0.1, bname="sdb1")
