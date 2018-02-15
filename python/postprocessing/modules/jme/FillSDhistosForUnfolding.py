import ROOT
import math, os
import numpy as np
import array as array
import fastjet
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.toolsFJ import matchFastJetObjectCollection


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
        #self.addObject(ROOT.TH1F('h_ak8sdm_'+self.bname,   'h_ak8sdm_'+self.bname, 25, 0, 250) )
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
        print "Create histograms"

    def endJob(self):
        Module.endJob(self)
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    
    def analyze(self, event):
        #print "analyze"
        """process event, return True (go to next module) or False (fail, go to next event)"""
        #jets = Collection(event, self.jetBranchName )
        #genJets = Collection(event, self.genJetBranchName )
        pfCands = Collection(event, self.pfCandsBranchName )
        genCands = Collection(event, self.genCandsBranchName )


        if len(pfCands) == 0 and len(genCands) ==0 :
            return False
        # Read the flags written by event selector 
        # 
        if event.miss < 1 and event.reco < 1  :
            return False


        print 'Event : {} goodgen : {} goodreco : {}  gen: {} miss: {} reco : {} fake : {}  response : {}  '.format( event.event , event.goodgen , event.goodreco, event.gen, event.miss, event.reco, event.fake, event.response )
        
               
        if not event.goodreco and not event.goodgen :
            return False
        if event.goodreco and event.goodgen :
            
            goodrecoP4 = ROOT.TLorentzVector( event.goodrecojet0_pt , event.goodrecojet0_eta , event.goodrecojet0_phi , event.goodrecojet0_m )
            
            pfCandsVec = ROOT.vector("TLorentzVector")()
            for p in pfCands :
                t = ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                drt = goodrecoP4.DeltaR(t) 
                if drt < 0.8 :
                   pfCandsVec.push_back( t )
                
            sdjets = self.sd.result( pfCandsVec )
            gsdjets = [ x for x in sdjets if x.perp() > 200.  and abs(x.eta()) < 2.5 ]
            gsdjets.sort(key=lambda x:x.perp(),reverse=True)

            goodgenP4 = ROOT.TLorentzVector( event.goodgenjet0_pt , event.goodgenjet0_eta , event.goodgenjet0_phi , event.goodgenjet0_m )
            genCandsVec = ROOT.vector("TLorentzVector")()
            for p in genCands :
                t = ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                drt = goodgenP4.DeltaR(t)
                if drt < 0.8 : genCandsVec.push_back( t )
            gensdjets = self.sd.result( genCandsVec )
            ggensdjets= [ x for x in gensdjets if x.perp() > 200.*0.8   ]
            ggensdjets.sort(key=lambda x:x.perp(),reverse=True)

            if len(ggensdjets) < 1 and len(gsdjets) < 1 :
                return False
          
            recoToGen = matchFastJetObjectCollection( gsdjets, ggensdjets, dRmax=0.05)
            print recoToGen
            for reco,gen in recoToGen.iteritems():
          
                if reco == None :
                    continue
                if event.reco > 0 : 
                    self.reco0.Fill(reco.m())
                    if self.verbose : print "Filling reco histo with SD jet of mass {:3.0f} GeV and Pt of {:3.0f} GeV".format(reco.m(), reco.perp())
                if gen != None :          
                    if event.response > 0  : 
                        self.resp0.Fill(reco.m(), gen.m() )
                        self.gen0.Fill(gen.m()) 
                        if self.verbose : print "Filling response and gen histo with gen SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV".format(gen.m(), gen.perp())
                  
                if event.fake > 0  :
                  
                    self.fake0.Fill(reco.m())
                    if self.verbose : print "Filling fake histo with SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV ".format(reco.m(), reco.perp())
                  
            for igen,gen in enumerate(gensdjets):
                if gen != None and gen not in recoToGen.values() :
                    if event.miss > 0 :
                        self.gen0.Fill(gen.m())
                        self.resp0.Fill(-1., gen.m() )
                        self.miss0.Fill(gen.m())        
                        if self.verbose : print "Filling miss/gen/response histo with gen SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV".format(gen.m(), gen.perp())        
                  
           
        elif (event.goodreco and not event.goodgen ) :
            #Fake
            goodrecoP4 = ROOT.TLorentzVector( event.goodrecojet0_pt , event.goodrecojet0_eta , event.goodrecojet0_phi , event.goodrecojet0_m )
            
            pfCandsVec = ROOT.vector("TLorentzVector")()
            for p in pfCands :
                t = ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                drt = goodrecoP4.DeltaR(t) 
                if drt < 0.8 :
                   pfCandsVec.push_back( t )
                
            sdjets = self.sd.result( pfCandsVec )
            gsdjets = [ x for x in sdjets if x.perp() > 200.  and abs(x.eta()) < 2.5 ]
            gsdjets.sort(key=lambda x:x.perp(),reverse=True)
            
            if len(gsdjets) < 1 : return False
            if event.fake < 0  : return False
            for reco in gsdjets:
                #if reco.m() < .01 : continue #return False
                self.fake0.Fill(reco.m())
                self.reco0.Fill(reco.m())
                if self.verbose : print "Filling reco fake histo with SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV".format(reco.m(), reco.perp())
                #self.resp0.Fill( reco.p4().M(), -1. ) 
                return True
        elif  (not event.goodreco and event.goodgen ) :
            #Miss
            goodgenP4 = ROOT.TLorentzVector( event.goodgenjet0_pt , event.goodgenjet0_eta , event.goodgenjet0_phi , event.goodgenjet0_m )
            genCandsVec = ROOT.vector("TLorentzVector")()
            for p in genCands :
                t = ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E())
                drt = goodgenP4.DeltaR(t)
                if drt < 0.8 : genCandsVec.push_back( t )
            gensdjets = self.sd.result( genCandsVec )
            ggensdjets= [ x for x in gensdjets if x.perp() > 200.*0.8   ]
            ggensdjets.sort(key=lambda x:x.perp(),reverse=True)
            
            if event.miss <  0 : return False
            if len(ggensdjets) < 1 : return False
            for gen in ggensdjets :
                #if gen.m() < .01 : continue
                self.miss0.Fill(gen.m())
                self.gen0.Fill(gen.m())
                self.resp0.Fill( -1.0, gen.m() )
                if self.verbose : print "Filling miss/gen/response histo with gen SD jet of mass {:3.0f} GeV  and Pt of {:3.0f} GeV".format(gen.m(), gen.perp())
                return True
        
        typeofill = ''
        if event.miss> 0 : typeofill = 'miss'
        if event.fake >0:typeofill = 'fake'
        if event.reco>0 and not event.fake<0 : typeofill = 'reco'
        if event.reco>0 and event.gen>0 : typeofill = 'gen'
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
