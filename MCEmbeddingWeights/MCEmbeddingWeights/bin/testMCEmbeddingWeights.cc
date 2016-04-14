
/**
   \class MCEmbeddingWeights MCEmbeddingWeights.cc "MCEmbeddingWeights/MCEmbeddingWeights/bin/MCEmbeddingWeights.cc"
   \brief Basic example of the use of the MCEmbeddingWeights class
*/

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "MCEmbeddingWeights/MCEmbeddingWeights/interface/MCEmbeddingWeights.h"

#include <Math/VectorUtil.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TBenchmark.h>

const double muonMass          = 0.105658369;
const double nomMassZ          = 91.1876;
const double breitWignerWidthZ = 2.4952;

namespace
{
  double square(double x)
  {
    return x*x;
  }
}

int main(int argc, char* argv[]) 
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);

  // load framework libraries
  gSystem->Load("libFWCoreFWLite");
  FWLiteEnabler::enable();

  // CV: mutau channel with muonPt > 19 GeV and tauPt > 20 GeV
  int decayMode1 = MCEmbeddingWeights::kTauToMuDecay;
  double ptMin1 = 19.;
  double absEtaMax1 = 2.1;
  int decayMode2 = MCEmbeddingWeights::kTauToHadDecay;
  double ptMin2 = 20.;
  double absEtaMax2 = 2.3;

  // CV: set center-of-mass energy to 13 TeV (LHC run 2)
  double sqrtS = 13.e+3; 

  int verbosity = 0;

  MCEmbeddingWeights mcEmbeddingWeights(decayMode1, ptMin1, absEtaMax1, decayMode2, ptMin2, absEtaMax2, sqrtS, verbosity);

  TRandom3 rnd;

  TH1* histogram_evtWeight_woTauSpinner = new TH1D("evtWeight_woTauSpinner", "evtWeight_woTauSpinner", 101, -0.005, +1.005);
  TH1* histogram_evtWeight_wTauSpinner  = new TH1D("evtWeight_wTauSpinner",  "evtWeight_wTauSpinner",  101, -0.005, +1.005);

  TBenchmark* clock = new TBenchmark();
  clock->Start("<testMCEmbeddingWeights>");

  const int numToys = 1000;
  for ( int iToy = 0; iToy < numToys; ++iToy ) {    
    double zPt_lab  = TMath::Abs(rnd.Gaus(0., 50.));
    double zEta_lab = rnd.Uniform(-3., +3.);
    double zPhi_lab = rnd.Uniform(-TMath::Pi(), +TMath::Pi());
    double zMass    = rnd.BreitWigner(nomMassZ, breitWignerWidthZ);
    double zP_lab   = zPt_lab*TMath::CosH(zEta_lab); // CV: taken from https://en.wikipedia.org/wiki/Pseudorapidity
    double zEn_lab  = TMath::Sqrt(square(zP_lab) + square(zMass));
    double zPx_lab  = zPt_lab*TMath::Cos(zPhi_lab);
    double zPy_lab  = zPt_lab*TMath::Sin(zPhi_lab);
    double zPz_lab  = zPt_lab*TMath::SinH(zEta_lab); // CV: taken from https://en.wikipedia.org/wiki/Pseudorapidity
    reco::Candidate::LorentzVector zP4_lab(zPx_lab, zPy_lab, zPz_lab, zEn_lab);

    double muon1En_rf = 0.5*zMass;
    double muon1Px_rf, muon1Py_rf, muon1Pz_rf;
    rnd.Sphere(muon1Px_rf, muon1Py_rf, muon1Pz_rf, muon1En_rf);
    reco::Candidate::LorentzVector muon1P4_rf(muon1Px_rf, muon1Py_rf, muon1Pz_rf, muon1En_rf);
    double muon2En_rf = 0.5*zMass;
    double muon2Px_rf = -muon1Px_rf;
    double muon2Py_rf = -muon1Py_rf;
    double muon2Pz_rf = -muon1Pz_rf;
    reco::Candidate::LorentzVector muon2P4_rf(muon2Px_rf, muon2Py_rf, muon2Pz_rf, muon2En_rf);

    ROOT::Math::Boost boost_to_rf(zP4_lab.BoostToCM());
    ROOT::Math::Boost boost_to_lab(boost_to_rf.Inverse());

    reco::Particle::LorentzVector muon1P4_lab = boost_to_lab(muon1P4_rf);
    reco::Particle::LorentzVector muon2P4_lab = boost_to_lab(muon2P4_rf);

    if ( !((muon1P4_lab.pt() > 19 && TMath::Abs(muon2P4_lab.eta()) < 2.1 && muon1P4_lab.pt() > 20 && TMath::Abs(muon2P4_lab.eta()) < 2.3) ||
	   (muon1P4_lab.pt() > 20 && TMath::Abs(muon2P4_lab.eta()) < 2.3 && muon1P4_lab.pt() > 19 && TMath::Abs(muon2P4_lab.eta()) < 2.1)) ) {
      continue;
    }
	 
    double evtWeight_woTauSpinner = mcEmbeddingWeights(muon1P4_lab, +1, muon2P4_lab, -1, false);
    histogram_evtWeight_woTauSpinner->Fill(evtWeight_woTauSpinner);
    double evtWeight_wTauSpinner = mcEmbeddingWeights(muon1P4_lab, +1, muon2P4_lab, -1, true);
    //double evtWeight_wTauSpinner = -1.;
    histogram_evtWeight_wTauSpinner->Fill(evtWeight_wTauSpinner);

    std::cout << "toy #" << iToy << ":" << std::endl;
    //std::cout << "rf:" << std::endl;
    //reco::Particle::LorentzVector zP4_rf = boost_to_rf(zP4_lab);
    //std::cout << " Z: pT = " << zP4_rf.pt() << ", eta = " << zP4_rf.eta() << ", phi = " << zP4_rf.phi() << ", mass = " << zP4_rf.mass() << std::endl;
    //std::cout << " mu+: pT = " << muon1P4_rf.pt() << ", eta = " << muon1P4_rf.eta() << ", phi = " << muon1P4_rf.phi() << ", mass = " << muon1P4_rf.mass() << std::endl;
    //std::cout << " mu-: pT = " << muon2P4_rf.pt() << ", eta = " << muon2P4_rf.eta() << ", phi = " << muon2P4_rf.phi() << ", mass = " << muon2P4_rf.mass() << std::endl;
    //reco::Particle::LorentzVector dimuonP4_rf = muon1P4_rf + muon2P4_rf;
    //std::cout << " dimuon: pT = " << dimuonP4_rf.pt() << ", eta = " << dimuonP4_rf.eta() << ", phi = " << dimuonP4_rf.phi() << ", mass = " << dimuonP4_rf.mass() << std::endl;
    //std::cout << "lab:" << std::endl;
    std::cout << " Z: pT = " << zP4_lab.pt() << ", eta = " << zP4_lab.eta() << ", phi = " << zP4_lab.phi() << ", mass = " << zP4_lab.mass() << std::endl;
    std::cout << " mu+: pT = " << muon1P4_lab.pt() << ", eta = " << muon1P4_lab.eta() << ", phi = " << muon1P4_lab.phi() << ", mass = " << muon1P4_lab.mass() << std::endl;
    std::cout << " mu-: pT = " << muon2P4_lab.pt() << ", eta = " << muon2P4_lab.eta() << ", phi = " << muon2P4_lab.phi() << ", mass = " << muon2P4_lab.mass() << std::endl;
    reco::Particle::LorentzVector dimuonP4_lab = muon1P4_lab + muon2P4_lab;
    std::cout << " dimuon: pT = " << dimuonP4_lab.pt() << ", eta = " << dimuonP4_lab.eta() << ", phi = " << dimuonP4_lab.phi() << ", mass = " << dimuonP4_lab.mass() << std::endl;
    std::cout << "--> weight: woTauSpinner = " << evtWeight_woTauSpinner << ", wTauSpinner = " << evtWeight_wTauSpinner << std::endl;    
  }

  clock->Stop("<testMCEmbeddingWeights>");
  clock->Show("<testMCEmbeddingWeights>");
  delete clock;
  
  TFile* outputFile = new TFile("testMCEmbeddingWeights.root", "RECREATE");
  histogram_evtWeight_woTauSpinner->Write();
  histogram_evtWeight_wTauSpinner->Write();
  delete outputFile;

  return 0;
}
