#include "MCEmbeddingWeights/MCEmbeddingWeights/interface/MCEmbeddingWeights.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "GeneratorInterface/TauolaInterface/interface/TauolaFactory.h"
#include "MCEmbeddingWeights/MCEmbeddingWeights/interface/read_particles_from_HepMC.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include "Tauola/Tauola.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "TauSpinner/SimpleParticle.h"
#include "TauSpinner/Tauola_wrapper.h"

#include "HepPDT/TableBuilder.hh"
#include "CLHEP/Random/JamesRandom.h"
#include "LHAPDF/LHAPDF.h"

#include <Math/VectorUtil.h>
#include <TMath.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <queue>

const double tauMass           = 1.77690;
const double muonMass          = 0.105658369;
const double electronMass      = 0.00051099893; // GeV
const double protonMass        = 0.938272;
const double nomMassZ          = 91.1876;
const double breitWignerWidthZ = 2.4952;

bool MCEmbeddingWeights::tauola_isInitialized_ = false;

namespace
{
  HepPDT::ParticleDataTable* loadParticleDataTable(const std::string& inputFileName)
  {  
    std::ifstream inputFile(inputFileName);
    if( !inputFile ) { 
      std::cerr << "Failed to open inputFile = " << inputFileName << "!!" << std::endl;
      assert(0);
    }
    HepPDT::ParticleDataTable* pdgTable = new HepPDT::ParticleDataTable();
    HepPDT::TableBuilder* tableBuilder = new HepPDT::TableBuilder(*pdgTable);
    bool status = addParticleTable(inputFile, *tableBuilder, true);
    if( !status ) { 
      std::cout << "Failed to read ParticleDataTable from file = " << inputFileName << "!!" << std::endl;
      assert(0);
    }
    delete tableBuilder;
    return pdgTable;
  }
}

MCEmbeddingWeights::MCEmbeddingWeights(int decayMode1, double ptMin1, double absEtaMax1, int decayMode2, double ptMin2, double absEtaMax2, double sqrtS, int verbosity)
  : decayMode1_(decayMode1),
    ptMin1_(ptMin1),
    absEtaMax1_(absEtaMax1),
    decayMode2_(decayMode2),
    ptMin2_(ptMin2),
    absEtaMax2_(absEtaMax2),
    motherParticleID_(23),
    beamEnergy_(0.5*sqrtS),
    generatorMode_("Tauola"),
    maxNumberOfAttempts_(1000),
    targetParticle1Mass_(tauMass),
    targetParticle1AbsPdgID_(15),
    targetParticle2Mass_(tauMass),
    targetParticle2AbsPdgID_(15),    
    tauola_(0),
    pdgTable_(0),
    rnd_(0),
    verbosity_(verbosity)
{
  // set TAUOLA mdtau parameter, defined at https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTauolaInterface
  int mdtau = -1;
  if      ( decayMode1 == kTauToElecDecay && decayMode2 == kTauToElecDecay ) mdtau = 121;
  else if ( decayMode1 == kTauToElecDecay && decayMode2 == kTauToMuDecay   ) mdtau = 123;
  else if ( decayMode1 == kTauToElecDecay && decayMode2 == kTauToHadDecay  ) mdtau = 115;
  else if ( decayMode1 == kTauToMuDecay   && decayMode2 == kTauToElecDecay ) mdtau = 123;
  else if ( decayMode1 == kTauToMuDecay   && decayMode2 == kTauToMuDecay   ) mdtau = 122;
  else if ( decayMode1 == kTauToMuDecay   && decayMode2 == kTauToHadDecay  ) mdtau = 116;
  else if ( decayMode1 == kTauToHadDecay  && decayMode2 == kTauToElecDecay ) mdtau = 115;
  else if ( decayMode1 == kTauToHadDecay  && decayMode2 == kTauToMuDecay   ) mdtau = 116;
  else if ( decayMode1 == kTauToHadDecay  && decayMode2 == kTauToHadDecay  ) mdtau = 132;

  edm::FileInPath inputFileName_pdgTable("SimGeneral/HepPDTESSource/data/pythiaparticle.tbl"); // CV: as defined in SimGeneral/HepPDTESSource/python/pythiapdt_cfi.py
  //pdgTable_ = new HepPDT::ParticleDataTable(inputFileName_pdgTable.fullPath());
  pdgTable_ = loadParticleDataTable(inputFileName_pdgTable.fullPath());
  if ( verbosity_ >= 2 ) {
    pdgTable_->writeParticleData(std::cout);
  }
  rnd_ = new CLHEP::HepJamesRandom(12345);

  edm::ParameterSet cfgInputCards;
  cfgInputCards.addParameter<int>("mdtau", mdtau);
  cfgInputCards.addParameter<int>("pjak1", 0);
  cfgInputCards.addParameter<int>("pjak2", 0);
  edm::ParameterSet cfgTauola;
  cfgTauola.addParameter<bool>("UseTauolaPolarization", false);
  cfgTauola.addParameter<edm::ParameterSet>("InputCards", cfgInputCards);
  tauola_ = TauolaFactory::get()->create("Tauolapp114", cfgTauola);
  tauola_->setRandomEngine(rnd_);
  tauola_->init(pdgTable_);
  tauola_isInitialized_ = true;

  //Tauolapp::Tauola::initialize();
  LHAPDF::initPDFSetByName("MSTW2008nnlo90cl.LHgrid");

  TauSpinner::initialize_spinner(true, 0, 0, 0, 1.3e+4);

  if ( decayMode1_ == decayMode2_ ) { // CV: consider both combinations of "leading" and "subleading" particles if both taus decay to particles of same type
    VisPtAndEtaCut visPtAndEtaCut1_1(decayMode1_, 1, ptMin1_, absEtaMax1_);
    VisPtAndEtaCut visPtAndEtaCut2_1(decayMode2_, 0, ptMin2_, absEtaMax2_);
    VisPtAndEtaCutCombination visPtAndEtaCut12;
    visPtAndEtaCut12.cuts_.push_back(visPtAndEtaCut1_1);
    visPtAndEtaCut12.cuts_.push_back(visPtAndEtaCut2_1);
    visPtAndEtaCuts_.push_back(visPtAndEtaCut12);
    VisPtAndEtaCut visPtAndEtaCut1_2(decayMode1_, 0, ptMin1_, absEtaMax1_);
    VisPtAndEtaCut visPtAndEtaCut2_2(decayMode2_, 1, ptMin2_, absEtaMax2_);
    VisPtAndEtaCutCombination visPtAndEtaCut21;
    visPtAndEtaCut21.cuts_.push_back(visPtAndEtaCut1_2);
    visPtAndEtaCut21.cuts_.push_back(visPtAndEtaCut2_2);
    visPtAndEtaCuts_.push_back(visPtAndEtaCut21);
  } else {
    VisPtAndEtaCut visPtAndEtaCut1(decayMode1_, 0, ptMin1_, absEtaMax1_);
    VisPtAndEtaCut visPtAndEtaCut2(decayMode2_, 0, ptMin2_, absEtaMax2_);
    VisPtAndEtaCutCombination visPtAndEtaCut;
    visPtAndEtaCut.cuts_.push_back(visPtAndEtaCut1);
    visPtAndEtaCut.cuts_.push_back(visPtAndEtaCut2);
    visPtAndEtaCuts_.push_back(visPtAndEtaCut);
  }
}

MCEmbeddingWeights::~MCEmbeddingWeights()
{
  delete tauola_;
  delete pdgTable_;
}

namespace
{
  double square(double x)
  {
    return x*x;
  }
}

void MCEmbeddingWeights::transformMuMu2LepLep(reco::Particle* muon1, reco::Particle* muon2)
{
//--- transform a muon pair into an electron/tau pair,
//    taking into account the difference between muon and electron/tau mass

  reco::Particle::LorentzVector muon1P4_lab = muon1->p4();
  reco::Particle::LorentzVector muon2P4_lab = muon2->p4();
  reco::Particle::LorentzVector zP4_lab = muon1P4_lab + muon2P4_lab;

  ROOT::Math::Boost boost_to_rf(zP4_lab.BoostToCM());
  ROOT::Math::Boost boost_to_lab(boost_to_rf.Inverse());

  reco::Particle::LorentzVector zP4_rf = boost_to_rf(zP4_lab);

  reco::Particle::LorentzVector muon1P4_rf = boost_to_rf(muon1P4_lab);
  reco::Particle::LorentzVector muon2P4_rf = boost_to_rf(muon2P4_lab);

  double muon1P_rf2 = square(muon1P4_rf.px()) + square(muon1P4_rf.py()) + square(muon1P4_rf.pz());
  double lep1Mass2 = square(targetParticle1Mass_);
  double lep1En_rf = 0.5*zP4_rf.E();
  double lep1P_rf2  = square(lep1En_rf) - lep1Mass2;
  if ( lep1P_rf2 < 0. ) lep1P_rf2 = 0.;
  float scaleFactor1 = sqrt(lep1P_rf2/muon1P_rf2);
  reco::Particle::LorentzVector lep1P4_rf = reco::Particle::LorentzVector(
    scaleFactor1*muon1P4_rf.px(), scaleFactor1*muon1P4_rf.py(), scaleFactor1*muon1P4_rf.pz(), lep1En_rf);

  double muon2P_rf2 = square(muon2P4_rf.px()) + square(muon2P4_rf.py()) + square(muon2P4_rf.pz());
  double lep2Mass2 = square(targetParticle2Mass_);
  double lep2En_rf = 0.5*zP4_rf.E();
  double lep2P_rf2  = square(lep2En_rf) - lep2Mass2;
  if ( lep2P_rf2 < 0. ) lep2P_rf2 = 0.;
  float scaleFactor2 = sqrt(lep2P_rf2/muon2P_rf2);
  reco::Particle::LorentzVector lep2P4_rf = reco::Particle::LorentzVector(
    scaleFactor2*muon2P4_rf.px(), scaleFactor2*muon2P4_rf.py(), scaleFactor2*muon2P4_rf.pz(), lep2En_rf);

  reco::Particle::LorentzVector lep1P4_lab = boost_to_lab(lep1P4_rf);
  reco::Particle::LorentzVector lep2P4_lab = boost_to_lab(lep2P4_rf);

  // perform additional checks:
  // the following tests guarantee a deviation of less than 0.1% 
  // for the following values of the original muons and the embedded electrons/taus in terms of:
  //  - invariant mass
  //  - transverse momentum
  if ( !(fabs(zP4_lab.mass() - (lep1P4_lab + lep2P4_lab).mass()) < (1.e-3*zP4_lab.mass()) &&
	 fabs(zP4_lab.pt()   - (lep1P4_lab + lep2P4_lab).pt())   < (1.e-3*zP4_lab.pt())) ) {
    edm::LogError ("MCEmbeddingWeights::transformMuMu2LepLep") 
      << "The kinematics of muons and embedded taus differ by more than 0.1%:" << std::endl 
      << " mass(muon1 + muon2) = " << zP4_lab.mass() << ", mass(lep1 + lep2) = " << (lep1P4_lab + lep2P4_lab).mass() << std::endl
      << " Pt(muon1 + muon2) = " << zP4_lab.pt() << ", Pt(lep1 + lep2) = " << (lep1P4_lab + lep2P4_lab).pt() << " --> please CHECK !!" << std::endl;
  }
  
  muon1->setP4(lep1P4_lab);
  muon2->setP4(lep2P4_lab);

  muon1->setPdgId(targetParticle1AbsPdgID_*muon1->pdgId()/abs(muon1->pdgId())); 
  muon2->setPdgId(targetParticle2AbsPdgID_*muon2->pdgId()/abs(muon2->pdgId()));

  muon1->setStatus(1);
  muon2->setStatus(1);
}

namespace
{
  struct pT_and_etaType
  {
    pT_and_etaType(double pT, double absEta)
      : pT_(pT),
	absEta_(absEta)
    {}
    ~pT_and_etaType() {}
    double pT_;
    double absEta_;
  };

  bool isHigherPt(const pT_and_etaType& pT_and_eta1, const pT_and_etaType& pT_and_eta2)
  {
    return (pT_and_eta1.pT_ > pT_and_eta2.pT_);
  }
}

bool MCEmbeddingWeights::testEvent(HepMC::GenEvent* genEvt)
{
  std::vector<pT_and_etaType> muon_pT_and_eta;
  std::vector<pT_and_etaType> electron_pT_and_eta;
  std::vector<pT_and_etaType> hadTau_pT_and_eta;

  int genParticleIdx = 0;
  for ( HepMC::GenEvent::particle_iterator genParticle = genEvt->particles_begin();
	genParticle != genEvt->particles_end(); ++genParticle ) {
    if ( verbosity_ >= 1 ) {
      std::cout << "genParticle #" << (genParticleIdx + 1) << " (pdgId = " << (*genParticle)->pdg_id() << "):" 
		<< " pT = " << (*genParticle)->momentum().perp() << ", eta = " << (*genParticle)->momentum().eta() << ", phi = " << (*genParticle)->momentum().phi() 
		<< std::endl;
    }
    if ( abs((*genParticle)->pdg_id()) == 15 && (*genParticle)->end_vertex() ) {
      reco::Candidate::LorentzVector visP4;
      std::queue<const HepMC::GenParticle*> decayProducts;
      decayProducts.push(*genParticle);
      int type = kTauToHadDecay;
      int decayProductIdx = 0;
      while ( !decayProducts.empty() && decayProductIdx < 100 ) { // CV: protection against entering infinite loop in case of corrupt particle relations
	const HepMC::GenParticle* decayProduct = decayProducts.front();
	if ( verbosity_ >= 1 ) {
	  std::cout << "decayProduct #" << (decayProductIdx + 1) << " (pdgId = " << decayProduct->pdg_id() << "):" 
		    << " pT = " << decayProduct->momentum().perp() << ", eta = " << decayProduct->momentum().eta() << ", phi = " << decayProduct->momentum().phi() 
		    << std::endl;
	}
	decayProducts.pop();
	if ( !decayProduct->end_vertex() ) { // stable decay product
	  int absPdgId = abs(decayProduct->pdg_id());
	  if ( !(absPdgId == 12 || absPdgId == 14 || absPdgId == 16) ) visP4 += (reco::Candidate::LorentzVector)decayProduct->momentum();
          if ( absPdgId == 11 ) type = kTauToElecDecay;
	  if ( absPdgId == 13 ) type = kTauToMuDecay;
	} else { // decay product decays further...
	  HepMC::GenVertex* decayVtx = decayProduct->end_vertex();
	  for ( HepMC::GenVertex::particles_out_const_iterator daughter = decayVtx->particles_out_const_begin();
		daughter != decayVtx->particles_out_const_end(); ++daughter ) {
	    decayProducts.push(*daughter);
	  }
	}
	++decayProductIdx;
      }

      double visPt = visP4.pt();
      if ( TMath::IsNaN(visPt) ) {
	edm::LogWarning ("MCEmbeddingWeights::testEvent")
	  << " Transverse momentum of visible tau decay products is NaN !!" << std::endl;
	return false;
      }
      double visAbsEta = fabs(visP4.eta());
      if      ( type == kTauToMuDecay   ) muon_pT_and_eta.push_back(pT_and_etaType(visPt, visAbsEta));
      else if ( type == kTauToElecDecay ) electron_pT_and_eta.push_back(pT_and_etaType(visPt, visAbsEta));
      else if ( type == kTauToHadDecay  ) hadTau_pT_and_eta.push_back(pT_and_etaType(visPt, visAbsEta));
      if ( verbosity_ >= 1 ) {
        std::string type_string = "";
        if      ( type == kTauToMuDecay   ) type_string = "mu";
        else if ( type == kTauToElecDecay ) type_string = "elec";
        else if ( type == kTauToHadDecay  ) type_string = "had";
        std::cout << "visLeg #" << (genParticleIdx + 1) << " (type = " << type_string << "):" 
		  << " Pt = " << visP4.pt() << ", eta = " << visP4.eta() << ", phi = " << visP4.phi() 
		  << " (X = " << (visP4.energy()/(*genParticle)->momentum().e()) << ")" << std::endl;
      }
      ++genParticleIdx;
    }
  }

  std::sort(electron_pT_and_eta.begin(), electron_pT_and_eta.end(), isHigherPt);
  std::sort(muon_pT_and_eta.begin(), muon_pT_and_eta.end(), isHigherPt);
  std::sort(hadTau_pT_and_eta.begin(), hadTau_pT_and_eta.end(), isHigherPt);
    
  // check if visible decay products pass Pt cuts
  //
  // NOTE: return value = True if leg1 > lowerThreshold[i] && leg2 > lowerThreshold[i] for **any** path i
  //      (e.g. (leg1Pt > 10 && leg2Pt > 20) || (leg1Pt > 20 && leg2Pt > 10), consistent with logic used by HLT)
  //
  for ( std::vector<VisPtAndEtaCutCombination>::const_iterator visPtAndEtaCut = visPtAndEtaCuts_.begin();
	visPtAndEtaCut != visPtAndEtaCuts_.end(); ++visPtAndEtaCut ) {
    if ( verbosity_ >= 2 ) {
      visPtAndEtaCut->print(std::cout);
    }

    bool passesVisPtAndEtaCut = true;
   
    for ( std::vector<VisPtAndEtaCut>::const_iterator cut = visPtAndEtaCut->cuts_.begin();
	  cut != visPtAndEtaCut->cuts_.end(); ++cut ) {
      std::vector<pT_and_etaType>* collection = 0;
      switch ( cut->type_ ) {
      case kTauToElecDecay:
	collection = &electron_pT_and_eta; 
	break;
      case kTauToMuDecay: 
	collection = &muon_pT_and_eta; 
	break;
      case kTauToHadDecay: 
	collection = &hadTau_pT_and_eta; 
	break;
      }
      assert(collection);
      
      // j-th tau decay product fails pT cut on visible tau decay products
      if ( !(cut->index_ < collection->size() && (*collection)[cut->index_].pT_ > cut->minPt_ && (*collection)[cut->index_].absEta_ < cut->maxAbsEta_) ) {
	passesVisPtAndEtaCut = false;
	break;
      }
    }

    // all visible tau decay products satisfy pT cuts for i-th path 
    if ( verbosity_ >= 2 ) {
      std::cout << "passes vis. pT and eta cuts = " << passesVisPtAndEtaCut << std::endl;
    }
    if ( passesVisPtAndEtaCut ) return true;
  }

  // visible pT and eta cuts failed for all paths
  return false;
}

//-------------------------------------------------------------------------------
// CV: auxiliary functions copied from TauAnalysis/MCEmbeddingTools/src/embeddingAuxFunctions.cc,
//     define here to make MCEmbeddingWeights independent of TauAnalysis/MCEmbeddingTools package
namespace{
  
void repairBarcodes(HepMC::GenEvent* genEvt)
{
  // AB: Note that we cannot do the barcode re-assignment "inline" without first
  //     creating a copy of the vertex and particle collections, because otherwise
  //     changing a barcode might invalidate the iterator which is used for
  //     iterating over the collection.
  const std::vector<HepMC::GenVertex*> vertices(genEvt->vertices_begin(), genEvt->vertices_end());
  const std::vector<HepMC::GenParticle*> particles(genEvt->particles_begin(), genEvt->particles_end());

  int next_genVtx_barcode = 1;
  for ( std::vector<HepMC::GenVertex*>::const_iterator iter = vertices.begin(); 
	iter != vertices.end(); ++iter, ++next_genVtx_barcode ) {
    while ( !(*iter)->suggest_barcode(-next_genVtx_barcode) ) {
      ++next_genVtx_barcode;
    }
  }

  int next_genParticle_barcode = 1;
  for ( std::vector<HepMC::GenParticle*>::const_iterator iter = particles.begin(); 
	iter != particles.end(); ++iter, ++next_genParticle_barcode ) {
    while ( !(*iter)->suggest_barcode(next_genParticle_barcode) ) {
      ++next_genParticle_barcode;
    }
  }
}

}
//-------------------------------------------------------------------------------

double MCEmbeddingWeights::operator()(const reco::Candidate::LorentzVector& muon1P4, int muon1Charge, const reco::Candidate::LorentzVector& muon2P4, int muon2Charge, bool addTauSpinnerWeight)
{
  int muon1PdgId = -13*muon1Charge/abs(muon1Charge);
  reco::Particle::Point muon1Vtx(0.,0.,0.);
  int muon2PdgId = -13*muon2Charge/abs(muon2Charge);
  reco::Particle::Point muon2Vtx(0.,0.,0.);
  reco::Particle embedLepton1(muon1Charge, muon1P4, muon1Vtx, muon1PdgId, 0, true);
  reco::Particle embedLepton2(muon2Charge, muon2P4, muon2Vtx, muon2PdgId, 0, true);
  transformMuMu2LepLep(&embedLepton1, &embedLepton2);
  std::vector<reco::Particle> embedParticles;	
  embedParticles.push_back(embedLepton1);
  embedParticles.push_back(embedLepton2);

  reco::Particle::LorentzVector genMotherP4;
  double ppCollisionPosX = 0.;
  double ppCollisionPosY = 0.;
  double ppCollisionPosZ = 0.;
  int idx = 0;
  for ( std::vector<reco::Particle>::const_iterator embedParticle = embedParticles.begin();
	embedParticle != embedParticles.end(); ++embedParticle ) {
    if ( verbosity_ >= 2 ) {
      std::cout << "embedParticle #" << idx << ": Pt = " << embedParticle->pt() << "," 
		<< " eta = " << embedParticle->eta() << ", phi = " << embedParticle->phi() << ", mass = " << embedParticle->mass() << std::endl;
      std::cout << "(production vertex: x = " << embedParticle->vertex().x() << ", y = " << embedParticle->vertex().y() << ", z = " << embedParticle->vertex().z() << ")" << std::endl;
    }
    genMotherP4 += embedParticle->p4();
    const reco::Particle::Point& embedParticleVertex = embedParticle->vertex();
    ppCollisionPosX += embedParticleVertex.x();
    ppCollisionPosY += embedParticleVertex.y();
    ppCollisionPosZ += embedParticleVertex.z();
    ++idx;
  }

  int numEmbedParticles = embedParticles.size();
  if ( numEmbedParticles > 0 ) {
    ppCollisionPosX /= numEmbedParticles;
    ppCollisionPosY /= numEmbedParticles;
    ppCollisionPosZ /= numEmbedParticles;
  }

  HepMC::GenVertex* ppCollisionVtx = new HepMC::GenVertex(HepMC::FourVector(ppCollisionPosX*10., ppCollisionPosY*10., ppCollisionPosZ*10., 0.)); // convert from cm to mm
  double protonEn = beamEnergy_;
  double protonPz = sqrt(square(protonEn) - square(protonMass));
  ppCollisionVtx->add_particle_in(new HepMC::GenParticle(HepMC::FourVector(0., 0., +protonPz, protonEn), 2212, 3));
  ppCollisionVtx->add_particle_in(new HepMC::GenParticle(HepMC::FourVector(0., 0., -protonPz, protonEn), 2212, 3));

  HepMC::GenVertex* genMotherDecayVtx = new HepMC::GenVertex(HepMC::FourVector(ppCollisionPosX*10., ppCollisionPosY*10., ppCollisionPosZ*10., 0.)); // Z-bosons decay immediately
  HepMC::GenParticle* genMother = new HepMC::GenParticle((HepMC::FourVector)genMotherP4, motherParticleID_, (generatorMode_ == "Pythia" ? 3 : 2), HepMC::Flow(), HepMC::Polarization(0,0));
  ppCollisionVtx->add_particle_out(genMother);

  genMotherDecayVtx->add_particle_in(genMother);
  for ( std::vector<reco::Particle>::const_iterator embedParticle = embedParticles.begin();
	embedParticle != embedParticles.end(); ++embedParticle ) {
    genMotherDecayVtx->add_particle_out(new HepMC::GenParticle((HepMC::FourVector)embedParticle->p4(), embedParticle->pdgId(), 1, HepMC::Flow(), HepMC::Polarization(0,0)));
  }

  HepMC::GenEvent* genEvt_output = new HepMC::GenEvent();
  genEvt_output->add_vertex(ppCollisionVtx);
  genEvt_output->add_vertex(genMotherDecayVtx);

  repairBarcodes(genEvt_output);

//--- compute probability for visible tau decay products to pass pT and eta cuts
  numEvents_tried_   = 0;
  sumWeights_tried_  = 0.;
  numEvents_passed_  = 0;
  sumWeights_passed_ = 0.;      
  HepMC::IO_HEPEVT conv;
  for ( int iTrial = 0; iTrial < maxNumberOfAttempts_; ++iTrial ) {
    HepMC::GenEvent* tempEvt = new HepMC::GenEvent(*genEvt_output);
    if ( generatorMode_ == "Pythia" )	{ // Pythia
      throw cms::Exception("Configuration") 
	<< "Pythia is currently not supported !!\n";
    } else if ( generatorMode_ == "Tauola" ) { // TAUOLA
      conv.write_event(tempEvt);      
      tempEvt = tauola_->decay(tempEvt);      
      if ( verbosity_ >= 2 ) {
	tempEvt->write(std::cout);
      }
    }

    double evtWeight = 1.;
    if ( addTauSpinnerWeight ) {
      TauSpinner::SimpleParticle mother, tau1, tau2;
      std::vector<TauSpinner::SimpleParticle> tau1_daughters, tau2_daughters;
      readParticlesFromHepMC(tempEvt, mother, tau1, tau2, tau1_daughters, tau2_daughters);
      evtWeight = TauSpinner::calculateWeightFromParticlesH(mother, tau1, tau2, tau1_daughters, tau2_daughters);
    }
    
    ++numEvents_tried_;
    sumWeights_tried_ += evtWeight;

    bool passesVisPtCuts = testEvent(tempEvt);
    if ( passesVisPtCuts ) {
      ++numEvents_passed_;
      sumWeights_passed_ += evtWeight;
    } 

    delete tempEvt;
  }

  delete genEvt_output;

  double retVal = ( numEvents_tried_ > 0 ) ? (sumWeights_passed_/sumWeights_tried_) : -1.;
  if ( verbosity_ >= 1 ) {
    std::cout << "--> retVal = " << retVal << std::endl;
  }
  return retVal;
}
