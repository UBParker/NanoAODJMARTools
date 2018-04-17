# NanoAODJMARTools
Tools for using the NANOAOD postprocessing framework for JMAR. 


## Obtaining JMAR NANOAOD

You should first follow the directions outlined in the [NanoAODJMAR](https://github.com/cms-jet/NanoAODJMAR) repository. This repository assumes the NANOAOD files follow that structure. 

## JMAR NANOAOD Analysis : With CMSSW

First, set up a new fastjet and fastjet-contrib. You need at least `fastjet-3.3.0`. This is in 
`/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/fastjet/3.3.0/` and `/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/fastjet-contrib/1.033/`. 

Now make a `CMSSW` working area and get this code:
```
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
git clone https://github.com/UBParker/NanoAODJMARTools.git PhysicsTools/NanoAODJMARTools
```


You can then use the following XML files to go into `$CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/`

```
cp $CMSSW_BASE/src/PhysicsTools/NanoAODJMARTools/xmlfiles/* $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/
scram setup fastjet
scram setup fastjet-contrib
```

Add fastjet to your python path:
```
setenv PYTHONPATH /uscms_data/d2/rappocc/fastjet/bare/install_330/lib/python2.7/site-packages:$PYTHONPATH

(or export in bash)
```


Compile and run:
```
scram b -j 10
cd PhysicsTools/NanoAODJMARTools/test
python postproc_softdrop.py
```

If you just want the library for the `SoftDrop` python interface, do this in python (after scram above):

```python 
import fastjet
import ROOT

ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")
beta=0.0
zcut=0.1
R=0.8
ptmin=200.
sd = ROOT.SoftDropWrapper(beta,zcut, R, ptmin)

#### IN YOUR EVENT LOOP:
        pfCandsVec = ROOT.vector("TLorentzVector")()
        for p in pfCands :
            pfCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
        sdjets = self.sd.result( pfCandsVec )
```



## JMAR NANOAOD Analysis : Without CMSSW 

Coming soon. 


## Technical details

This assumes you have [fastjet 3.3.0](http://fastjet.fr/repo/doxygen-3.3.0/), which implements the python front-end to fastjet. The `fastjet-contrib` packages do not yet have a python implementation, so this is implemented [here](https://github.com/cms-jet/NanoAODJMARTools/blob/master/src/Recluster.cc) and [here](https://github.com/cms-jet/NanoAODJMARTools/blob/master/interface/Recluster.h).



