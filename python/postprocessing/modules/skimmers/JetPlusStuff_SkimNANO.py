
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
import numpy
import random
import array
import copy


def myclosest(obj,collection,presel=lambda x,y: True):
    ret = None; drMin = 999
    for x in collection:
        if not presel(obj,x): continue
        dr = obj.DeltaR(x)
        if dr < drMin: 
            ret = x; drMin = dr
    return (ret,drMin)
    
algNames = [
    "sdB0", "sdB1", "sdBn1"
]

class ZPlusJet_SkimNANO(Module):
    def __init__(self ):
        self.writeHistFile = True
        self.verbose = False
        self.R = 0.8
        self.ptmin = 0.
        
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)

        ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")
        
        self.sdBn1Z0p05 = ROOT.SoftDropWrapper(-1. ,0.05, self.R, self.ptmin)
        self.sdBn1 = ROOT.SoftDropWrapper(-1. ,0.1, self.R, self.ptmin)
        self.sdBn1Z0p15 = ROOT.SoftDropWrapper(-1. ,0.15, self.R, self.ptmin)


        self.sdB0Z0p05 = ROOT.SoftDropWrapper(0. ,0.05, self.R, self.ptmin)
        self.sdB0 = ROOT.SoftDropWrapper(0. ,0.1, self.R, self.ptmin)
        self.sdB0Z0p15 = ROOT.SoftDropWrapper(0. ,0.15, self.R, self.ptmin)

        self.sdB1Z0p05 = ROOT.SoftDropWrapper(1. ,0.05, self.R, self.ptmin)
        self.sdB1 = ROOT.SoftDropWrapper(1. ,0.1, self.R, self.ptmin)
        self.sdB1Z0p15 = ROOT.SoftDropWrapper(1. ,0.15, self.R, self.ptmin)

        self.algsToRun = [ self.sdB0, self.sdB1, self.sdBn1 ]
        self.algNames = copy.copy(algNames)
        
        ### Kinematics Cuts ###
        ### This is for ANY analysis where the signal contains :
       
        ###  >=1  high Pt (> 200 GeV) AK8 Jet (Reco or Gen level)
        ###  >=1 Muon (Pt > 30 GeV) or >=1 Electron (Pt > 35 )     

        ### Considering either
        ###Z - > mu+ mu-  + Jet


        self.minMu0pt = 30.
        # Trigger: HLT_IsoTkMu22 || HLT_IsoMu22

        self.minEl0pt = 35.
        #Trigger: HLT_Ele32_eta2p1_WPTight_Gsf

        self.minJetPt = 200.

        self.maxObjEta = 2.5



        ### Control plots of observed particles ###


        ### h_ histograms contain muons and electrons
        self.addObject( ROOT.TH1D('h_lep0pt',          'h_lep0pt',        40, 0, 200 ) )
        self.addObject( ROOT.TH1D('h_lep0eta',         'h_lep0eta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_lep0phi',         'h_lep0phi',      100, -5, 5 ) )


        ### he_* histograms contain Electrons
        self.addObject( ROOT.TH1D('he_lep0pt',          'he_lep0pt',        40, 0, 200 ) )
        self.addObject( ROOT.TH1D('he_lep0eta',         'he_lep0eta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('he_lep0phi',         'he_lep0phi',      100, -5, 5 ) )

        ### hm_* histograms contain Muons
        self.addObject( ROOT.TH1D('hm_lep0pt',          'hm_lep0pt',        40, 0, 200 ) )
        self.addObject( ROOT.TH1D('hm_lep0eta',         'hm_lep0eta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('hm_lep0phi',         'hm_lep0phi',      100, -5, 5 ) )

        self.addObject( ROOT.TH1D('h_genjetpt',          'h_genjetpt',   100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_genjeteta',         'h_genjeteta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_genjetphi',         'h_genjetphi',      100, -5, 5 ) )
        self.addObject( ROOT.TH1D('h_genjetmass',        'h_genjetmass',      300, 0, 300 ) )
        #self.addObject( ROOT.TH1D('h_genjettau1',        'h_genjettau1',      100, 0, 1 ) )
        #self.addObject( ROOT.TH1D('h_genjettau2',        'h_genjettau2',      100, 0, 1 ) )
        #self.addObject( ROOT.TH1D('h_genjettau3',        'h_genjettau3',      100, 0, 1 ) )
        #self.addObject( ROOT.TH1D('h_genjettau4',        'h_genjettau4',      100, 0, 1 ) )
        #self.addObject( ROOT.TH1D('h_genjetn2b1',        'h_genjetn2b1',      100, 0, 1 ) )
        #self.addObject( ROOT.TH1D('h_genjetn3b1',        'h_genjetn3b1',      100, 0, 1 ) )


        self.addObject( ROOT.TH1D('h_recojetpt',          'h_recojetpt',  3000/5, 0, 3000 ) )
        self.addObject( ROOT.TH1D('h_recojeteta',         'h_recojeteta',      48, -3, 3 ) )
        self.addObject( ROOT.TH1D('h_recojetphi',         'h_recojetphi',      100, -5, 5 ) )

        self.addObject( ROOT.TH1D('h_recojetmass_u',        'h_recojetmass_u',      1000/5, 0, 1000 ) )

        self.addObject( ROOT.TH1D('h_recojettau1',        'h_recojettau1',      100, 0, 1 ) )
        self.addObject( ROOT.TH1D('h_recojettau2',        'h_recojettau2',      100, 0, 1 ) )
        self.addObject( ROOT.TH1D('h_recojettau3',        'h_recojettau3',      100, 0, 1 ) )
        self.addObject( ROOT.TH1D('h_recojettau4',        'h_recojettau4',      100, 0, 1 ) )
        self.addObject( ROOT.TH1D('h_recojetn2b1',        'h_recojetn2b1',      100, 0, 1 ) )
        self.addObject( ROOT.TH1D('h_recojetn3b1',        'h_recojetn3b1',      100, 0, 1 ) )

        self.addObject( ROOT.TH1D('h_drGenReco',    'h_drGenReco',   40, 0, 0.8) )
        self.addObject( ROOT.TH1D('h_drGenGroomed', 'h_drGenGroomed',40, 0, 0.8) )
                          
        self.sdrecoMass = {}
        self.sdgenMass= {}

        for ialg,alg in enumerate(self.algsToRun):
            h = ROOT.TH1D('h_recojetmass'+ self.algNames[ialg] ,        'h_recojetmass'+ self.algNames[ialg] ,      1000, 0, 1000 )
            self.sdrecoMass[self.algNames[ialg] ] = h
            self.addObject( h )

            j = ROOT.TH1D('h_genjetmass'+ self.algNames[ialg] ,        'h_genjetmass'+ self.algNames[ialg] ,      1000, 0, 1000 )
            self.sdgenMass[self.algNames[ialg] ] = j
            self.addObject( j )
  
    def endJob(self):
        Module.endJob(self)
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Zlep_isMu",  "I")
        self.out.branch("fake",  "I")
        self.out.branch("miss",  "I")
        self.out.branch("gen",  "I")
        self.out.branch("breco",  "I")
        self.out.branch("response",  "I")
        self.out.branch("goodgen",  "I")
        self.out.branch("goodreco",  "I") 

        self.out.branch("goodrecojet0_pt",  "F")
        self.out.branch("goodrecojet0_eta",  "F")
        self.out.branch("goodrecojet0_phi",  "F")
        self.out.branch("goodrecojet0_mass",  "F")

        self.out.branch("goodgenjet0_pt",  "F")
        self.out.branch("goodgenjet0_eta",  "F")
        self.out.branch("goodgenjet0_phi",  "F")
        self.out.branch("goodgenjet0_mass",  "F")


        for ialg,alg in enumerate(self.algsToRun):
            self.out.branch("goodreco" + self.algNames[ialg] + "jet0_pt",  "F")
            self.out.branch("goodreco" + self.algNames[ialg] + "jet0_eta",  "F")
            self.out.branch("goodreco" + self.algNames[ialg] + "jet0_phi",  "F")
            self.out.branch("goodreco" + self.algNames[ialg] + "jet0_mass",  "F")

            self.out.branch("goodgen" + self.algNames[ialg] + "jet0_pt",  "F")
            self.out.branch("goodgen" + self.algNames[ialg] + "jet0_eta",  "F")
            self.out.branch("goodgen" + self.algNames[ialg] + "jet0_phi",  "F")
            self.out.branch("goodgen" + self.algNames[ialg] + "jet0_mass",  "F")

        pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def getSubjets(self, p4, subjets, dRmax=0.8):
        ret = []
        for subjet in subjets :
            if p4.DeltaR(subjet.p4()) < dRmax and len(ret) < 2 :
                ret.append(subjet.p4())
        return ret

    def printP4( self, c ):
        if hasattr( c, "p4"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.p4().Perp(), c.p4().Eta(), c.p4().Phi(), c.p4().M() )
        elif hasattr( c, "Perp"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.Perp(), c.Eta(), c.Phi(), c.M() )
        else:
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.perp(), c.eta(), c.phi(), c.m() )
        return s
    def printCollection(self,coll):
        for ic,c in enumerate(coll):
            s = self.printP4( c )
            print ' %3d : %s' % ( ic, s )
            
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        weight = 1.0

        
        isMC = event.run == 1
        if self.verbose:
            print '------------------------ ', event.event
        genjets = []
        #print "len(genjets) = {}".format(len(genjets))
        if isMC:
            goodgen= False
            ###### Get gen Z candidate ######
            genleptons = Collection(event, "GenDressedLepton")
            genjet0 = None
            genjet0Groomed = []
            for il,genl in enumerate(genleptons) :
              
                if len(genleptons) < 1 : 
                    continue
                    #return False
                if abs(genl.pdgId) != 13 and  abs(genl.pdgId) != 11 :
                    continue
                    #return False
                if self.verbose :
                    print '----'
                    print 'Gen leptons:'
                    self.printCollection( genleptons )
           
                goodgen = True
                if self.verbose:
                    print '-----'
                    print 'Gen Z:'
                    print self.printP4( Zboson )
                if goodgen : break
            if goodgen :
                ###### Get list of gen jets #######
                allgenjets = list(Collection(event, "GenJetAK8"))
                allgenparts = list(Collection(event, "GenPartAK8"))

                genCandsVec = ROOT.vector("TLorentzVector")()
                for p in allgenparts :
                    genCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )

                

                # List of gen jets:
                genjets = [ x for x in allgenjets if x.p4().Perp() > self.minJetPt * 0.8  ]
                genjetsGroomed = {}
                
                
                # Create the map of ungroomed to groomed gen jets
                if len(genjets) < 1 :
                    goodgen = False
                else :
                    for igenjet,genjet in enumerate(genjets):
                        genjetsGroomed[igenjet] = {}
                        # Cluster only the particles near the appropriate jet to save time
                        constituents = ROOT.vector("TLorentzVector")()
                        for x in genCandsVec:
                            if genjet.p4().DeltaR( x ) < 0.8:
                                constituents.push_back(x)
                        for ialg,alg in enumerate(self.algsToRun):
                            groomedjetsFJ = alg.result( constituents )
                            groomedjets = [ ROOT.TLorentzVector( x.px(), x.py(), x.pz(), x.e() ) for x in groomedjetsFJ]                            
                            if len(groomedjets) > 0 : 
                                genjetsGroomed[igenjet][ialg] = groomedjets[0]
                            else:
                                print 'Grooming failed. Inputs are:'
                                self.printCollection( constituents )
                                genjetsGroomed[igenjet][ialg] = None
                    if self.verbose:
                        print '-----'
                        print ' genjets:'
                        for igen, gen in enumerate( genjets ):
                            print " %12d: %25s" % ( igen, self.printP4(gen))
                            for ialg,gengroomed in genjetsGroomed[igen].iteritems() :
                                if gengroomed :
                                    print " ---- alg %3d: %25s" % ( ialg, self.printP4(gengroomed) )
                    
                                
                            

        goodreco = False    

        ###### Get reco Z candidate #######
        # List of reco muons
        allmuons = Collection(event, "Muon")
        # Select reco muons:
        muons = [ x for x in allmuons if (x.mediumId  and x.p4().Perp() > self.minMu0pt and abs(x.p4().Eta()) < self.maxObjEta)]
 
        # List of reco electrons
        allelectrons = Collection(event, "Electron")
        # Select reco muons:
        electrons = [ x for x in allelectrons if (   x.p4().Perp() > self.minEl0pt and abs(x.p4().Eta()) < self.maxObjEta  )]


        if len(muons) < 1 and len(electrons) < 1:
            goodreco = False
        else : 
            goodreco = True             
       
        lep0 = None

        Zismu = 0   
        if len(muons) >= 1 and len(electrons) < 1:           
            lep0 = muons[0].p4()
            Zismu =1
        elif len(muons) < 1 and len(electrons) >= 1:           
            lep0 = electrons[0].p4()
        elif  len(muons) >= 1 and len(electrons) >= 1:
            if  muons[0].p4().Perp() >  electrons[0].p4().Perp() :
                Zismu =1
                lep0 = muons[0].p4()
    
            else :

                lep0 = electrons[0].p4()

        if goodreco :
            self.h_lep0pt.Fill(lep0.Perp())
            self.h_lep0eta.Fill(lep0.Eta())
            self.h_lep0phi.Fill(lep0.Phi())

            if Zismu == 0 :
                self.he_lep0pt.Fill(lep0.Perp())
                self.he_lep0eta.Fill(lep0.Eta())
                self.he_lep0phi.Fill(lep0.Phi())
            elif Zismu == 1 :
                self.hm_lep0pt.Fill(lep0.Perp())
                self.hm_lep0eta.Fill(lep0.Eta())
                self.hm_lep0phi.Fill(lep0.Phi())

        #if self.verbose:
        #    print '-----'
        #    print ' recoZ:', self.printP4( Zcand )
        
        ###### Get list of reco jets #######
        # List of reco jets:
        allrecojets = list(Collection(event, "FatJet"))
        allrecoparts = list(Collection(event, "PFCandsAK8"))

        recoCandsVec = ROOT.vector("TLorentzVector")()
        for p in allrecoparts :
            recoCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
                    

        if goodreco :
            recojets = [ x for x in allrecojets if x.p4().Perp() > self.minJetPt and  abs(x.eta) < self.maxObjEta  and x.p4().DeltaR(lep0) >0.8  ]
            recojets.sort(key=lambda x:x.p4().Perp(),reverse=True)
            if self.verbose:
                print '----'
                print ' recojets:'
                self.printCollection( recojets )
            recojetsGroomed = {}
            # Create the map of ungroomed to groomed reco jets
            if len(recojets) < 1 : goodreco = False
            else :
                for irecojet,recojet in enumerate(recojets):
                    recojetsGroomed[irecojet] = {}
                    # Cluster only the particles near the appropriate jet to save time
                    constituents = ROOT.vector("TLorentzVector")()
                    for x in recoCandsVec:
                        if recojet.p4().DeltaR( x ) < 0.8:
                            constituents.push_back(x)
                    for ialg,alg in enumerate(self.algsToRun):
                        groomedjetsFJ = alg.result( constituents )
                        groomedjets = [ ROOT.TLorentzVector( x.px(), x.py(), x.pz(), x.e() ) for x in groomedjetsFJ]
                        if len(groomedjets) > 0 : 
                            recojetsGroomed[irecojet][ialg] = groomedjets[0]
                        else:
                            recojetsGroomed[irecojet][ialg] = None
                            print 'RECO Grooming failed. Constituents:'
                            self.printCollection( constituents )
                if self.verbose:
                    print '-----'
                    print ' recojets:'
                    for ireco, reco in enumerate( recojets ):
                        print " %12d: %25s" % ( ireco, self.printP4(reco))
                        for ialg,recogroomed in recojetsGroomed[ireco].iteritems() :
                            if recogroomed :
                                print " ---- alg %3d: %25s" % ( ialg, self.printP4(recogroomed) )
                                
        else :
            recojets = []
            recojetsGroomed = []
   
       
        if isMC == False :
            goodgen = False
             
        # Dictionary to hold reco--> gen matching
        recoToGen =  None
        if isMC and goodgen and goodreco :
            recoToGen = matchObjectCollection( recojets, genjets, dRmax=0.05 )
    
                                
        if not goodgen and not goodreco :
            self.out.fillBranch("goodreco",  -2)
            self.out.fillBranch("goodgen",  -2)
            return False
        
        if not goodreco or len(recojets) < 1 or recojets[0] == None:
            self.out.fillBranch("goodreco",  0)
            self.out.fillBranch("goodrecojet0_pt",  -1.0 )
            self.out.fillBranch("goodrecojet0_eta", -1.0 )
            self.out.fillBranch("goodrecojet0_phi",  -1.0 )
            self.out.fillBranch("goodrecojet0_mass", -1.0 )
            for ialg,alg in enumerate(self.algsToRun):
                self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_pt",  -1.0 )
                self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_eta", -1.0 )
                self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_phi", -1.0 )
                self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_mass",   -1.0 )
                
        else:
            ireco = 0
            reco = recojets[ireco]
            self.out.fillBranch("goodreco",  1)
            self.out.fillBranch("goodrecojet0_pt",  reco.p4().Perp() )
            self.out.fillBranch("goodrecojet0_eta",  reco.p4().Eta() )
            self.out.fillBranch("goodrecojet0_phi",  reco.p4().Phi() )
            self.out.fillBranch("goodrecojet0_mass",  reco.p4().M() )
            self.h_recojetpt.Fill( reco.p4().Perp() )
            self.h_recojeteta.Fill( reco.p4().Eta() )
            self.h_recojetphi.Fill( reco.p4().Phi() )
            self.h_recojetmass_u.Fill( reco.p4().M()  )
            
            #########################################
            # TO DO: Need to apply jet corrections!!!
            #########################################

            
            for ialg,alg in enumerate(self.algsToRun):
                recoGroomed = recojetsGroomed[ireco][ialg] if ireco in recojetsGroomed and ialg in recojetsGroomed[ireco] else None
                if recoGroomed != None : 
                    self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_pt",   recoGroomed.Perp() )
                    self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_eta",  recoGroomed.Eta() )
                    self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_phi",  recoGroomed.Phi() )
                    self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_mass",    recoGroomed.M() )
                    self.sdrecoMass[ self.algNames[ialg] ].Fill(recoGroomed.M() )
                else : 
                    self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_pt",  -1.0 )
                    self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_eta", -1.0 )
                    self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_phi", -1.0 )
                    self.out.fillBranch("goodreco" + self.algNames[ialg] + "jet0_mass",   -1.0 )

        if not goodgen or len(genjets) < 1 or genjets[0] == None :
            self.out.fillBranch("goodgen",  0)
            self.out.fillBranch("goodgenjet0_pt",  -1.0 )
            self.out.fillBranch("goodgenjet0_eta", -1.0 )
            self.out.fillBranch("goodgenjet0_phi",  -1.0 )
            self.out.fillBranch("goodgenjet0_mass", -1.0 )
            for ialg,alg in enumerate(self.algsToRun):
                self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_pt",  -1.0 )
                self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_eta", -1.0 )
                self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_phi", -1.0 )
                self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_mass",   -1.0 )
        else :
            igen = 0
            gen = genjets[igen]
            self.out.fillBranch("goodgen",  1)            
            #print "filling Ungroomed histos---------"
            self.out.fillBranch("goodgenjet0_pt",  gen.p4().Perp() )
            self.out.fillBranch("goodgenjet0_eta",  gen.p4().Eta() )
            self.out.fillBranch("goodgenjet0_phi",  gen.p4().Phi() )
            self.out.fillBranch("goodgenjet0_mass",  gen.p4().M() )
            for ialg,alg in enumerate(self.algsToRun):
                genGroomed = genjetsGroomed[igen][ialg] if igen in genjetsGroomed and ialg in genjetsGroomed[igen] else None
                if genGroomed != None : 
                    self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_pt",   genGroomed.Perp() )
                    self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_eta",  genGroomed.Eta() )
                    self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_phi",  genGroomed.Phi() )
                    self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_mass",    genGroomed.M() )
                    self.sdgenMass[ self.algNames[ialg] ].Fill(genGroomed.M() )

                else : 
                    self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_pt",  -1.0 )
                    self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_eta", -1.0 )
                    self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_phi", -1.0 )
                    self.out.fillBranch("goodgen" + self.algNames[ialg] + "jet0_mass",   -1.0 )
            
            self.h_genjetpt.Fill( gen.p4().Perp() )                      
            self.h_genjeteta.Fill( gen.p4().Eta() )
            self.h_genjetphi.Fill( gen.p4().Phi() )
            self.h_genjetmass.Fill( gen.p4().M()  )
            #self.h_genjettau1.Fill( gen.tau1  )
            #self.h_genjettau2.Fill( gen.tau2  )
            #self.h_genjettau3.Fill( gen.tau3 )
            #self.h_genjettau4.Fill( gen.tau4  )
            #self.h_genjetn2b1.Fill( gen.n2b1 )
            #self.h_genjetn3b1.Fill( gen.n3b1 ) 


        return True
# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

zplusjetsxs = lambda : ZPlusJet_SkimNANO() 
