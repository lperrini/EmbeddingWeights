#ifndef MCEmbeddingWeights_MCEmbeddingWeights_MCEmbeddingWeights_h
#define MCEmbeddingWeights_MCEmbeddingWeights_MCEmbeddingWeights_h

/** \class MCEmbeddingWeights
 *
 * Compute probability of Z -> tautau event corresponding to given Z -> mumu event
 * to pass pT and eta cuts applied on visible tau decay products.
 *
 * \author Christian Veelken, Tallinn
 *
 */

#include "GeneratorInterface/TauolaInterface/interface/TauolaInterfaceBase.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "HepPDT/ParticleDataTable.hh"
#include "CLHEP/Random/RandomEngine.h"
#include "HepMC/IO_HEPEVT.h"

#include <vector>

class MCEmbeddingWeights
{
 public:
  enum { kUndefined, kTauToElecDecay, kTauToMuDecay, kTauToHadDecay };
  explicit MCEmbeddingWeights(int, double, double, int, double, double, double, int = 0);
  virtual ~MCEmbeddingWeights();

  double operator()(const reco::Candidate::LorentzVector&, int, const reco::Candidate::LorentzVector&, int, bool = true);

 protected:
  void transformMuMu2LepLep(reco::Particle*, reco::Particle*);
  bool testEvent(HepMC::GenEvent*);

  int decayMode1_;
  double ptMin1_;
  double absEtaMax1_;
  int decayMode2_;
  double ptMin2_;
  double absEtaMax2_;

  int motherParticleID_;
  double beamEnergy_;
  std::string generatorMode_;
  int maxNumberOfAttempts_;

  double targetParticle1Mass_;
  int targetParticle1AbsPdgID_;
  double targetParticle2Mass_;
  int targetParticle2AbsPdgID_;

  gen::TauolaInterfaceBase* tauola_;
  // keep track if TAUOLA/PHOTOS interface has already been initialized,
  // to avoid multiple initializations of TAUOLA interface, which makes TAUOLA crash.
  static bool tauola_isInitialized_;
  HepPDT::ParticleDataTable* pdgTable_;
  CLHEP::HepRandomEngine* rnd_;

  struct VisPtAndEtaCut 
  { 
    VisPtAndEtaCut(int type, unsigned index, double minPt, double maxAbsEta)
      : type_(type),
	index_(index),
	minPt_(minPt),
	maxAbsEta_(maxAbsEta)
    {}
    ~VisPtAndEtaCut() {}
    int type_; 
    unsigned index_; 
    double minPt_; 
    double maxAbsEta_; 
    void print(std::ostream& stream) const
    {
      std::string type_string = "undefined";
      if      ( type_ == kTauToElecDecay ) type_string = "elec";
      else if ( type_ == kTauToMuDecay   ) type_string = "mu";
      else if ( type_ == kTauToHadDecay  ) type_string = "had";
      stream << type_string << " #" << index_ << ": minPt = " << minPt_ << ", maxAbsEta = " << maxAbsEta_ << std::endl;
    }
  };
  struct VisPtAndEtaCutCombination
  {
    std::vector<VisPtAndEtaCut> cuts_;
    void print(std::ostream& stream) const
    {
      for ( std::vector<VisPtAndEtaCut>::const_iterator cut = cuts_.begin();
            cut != cuts_.end(); ++cut ) {
        cut->print(stream);
      }
    }
  };
  std::vector<VisPtAndEtaCutCombination> visPtAndEtaCuts_;

  unsigned numEvents_tried_;
  double sumWeights_tried_;
  unsigned numEvents_passed_;
  double sumWeights_passed_;

  int verbosity_;
};

#endif
