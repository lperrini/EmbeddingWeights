# EmbeddingWeights
```
cmsrel CMSSW_7_6_3/src/
cd CMSSW_7_6_3/src
cmsenv
git clone https://github.com/lperrini/EmbeddingWeights.git
mv EmbeddingWeights/* $CMSSW_BASE
scram b -j 9
```

run the test code MCEmbeddingWeights/MCEmbeddingWeights/bin/testMCEmbeddingWeights.cc by simply calling:
```
testMCEmbeddingWeights
```

the test program creates 1000 Z -> mumu toy MC events and fills histograms of the MCEmbeddingWeights without (evtWeight_woTauSpinner) and with (evtWeight_wTauSpinner) taking tau polarization effects into account
If you want to plug-in the class in your own analysis code, do (you find the entire example in the test programm):
```
#include "MCEmbeddingWeights/MCEmbeddingWeights/interface/MCEmbeddingWeights.h"
int decayMode1 = -999;
int decayMode2 = -999;
if(isem){ decayMode1=MCEmbeddingWeights::kTauToElecDecay; decayMode2=MCEmbeddingWeights::kTauToMuDecay;}
if(iset){ decayMode1=MCEmbeddingWeights::kTauToElecDecay; decayMode2=MCEmbeddingWeights::kTauToHadDecay;}
if(ismt){ decayMode1=MCEmbeddingWeights::kTauToMuDecay;   decayMode2=MCEmbeddingWeights::kTauToHadDecay;}
if(istt){ decayMode1=MCEmbeddingWeights::kTauToHadDecay;  decayMode2=MCEmbeddingWeights::kTauToHadDecay;}
MCEmbeddingWeights mcEmbeddingWeights(decayMode1, pt1, eta1, decayMode2, pt2, eta2, sqrtS,verbosity);
//compute weights
double evtWeight_woTauSpinner = mcEmbeddingWeights(muon1P4, q_1, muon2P4, q_2, false);
double evtWeight_wTauSpinner  = mcEmbeddingWeights(muon1P4, q_1, muon2P4, q_2, true);
```






