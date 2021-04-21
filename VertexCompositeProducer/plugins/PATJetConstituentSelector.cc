/* \class PATJetConstituentSelector
 */

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

class PATJetConstituentSelector : public edm::stream::EDProducer<> {
public:

  using JetsOutput = std::vector<pat::Jet>;
  using ConstituentsOutput = std::vector<pat::PackedCandidate>;
  
  PATJetConstituentSelector(edm::ParameterSet const& params) :
    srcToken_{consumes<edm::View<pat::Jet>>(params.getParameter<edm::InputTag>("src"))},
    selector_{params.getParameter<std::string>("cut")}
  {
//    produces<JetsOutput>();
    produces<ConstituentsOutput>();
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions)
  {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("InputTag used for retrieving jets in event.");
    desc.add<std::string>("cut")->setComment("Cut used by which to select jets.  For example:\n"
                                             "  \"pt > 100.0 && abs(rapidity()) < 2.4\".");

    // addDefault must be used here instead of add unless this function is specialized
    // for different sets of template parameter types. Each specialization would need
    // a different module label. Otherwise the generated cfi filenames will conflict
    // for the different plugins.
    descriptions.addDefault(desc);
  }

  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override
  {
    auto jets = std::make_unique<JetsOutput>();
    auto candsOut = std::make_unique<ConstituentsOutput>();

    edm::Handle<edm::View<pat::Jet>> h_jets;
    iEvent.getByToken(srcToken_, h_jets);
    
    // Now set the Ptrs with the orphan handles.
    for (auto const& jet : *h_jets) {
      // Check the selection
      if (selector_(jet)) {
	// Add the jets that pass to the output collection
	jets->push_back(jet);

        std::vector<const pat::PackedCandidate*> daughters;
        for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
          auto const* cand = jet.daughter(i);
          auto packed_cand = dynamic_cast<const pat::PackedCandidate*>(cand);
          if (packed_cand) {
            auto packed_cand = dynamic_cast<const pat::PackedCandidate*>(cand);
            daughters.push_back(packed_cand);
          }
        }
      }
    }

//    iEvent.put(std::move(jets));
    iEvent.put(std::move(candsOut));
  }

private:
  edm::EDGetTokenT<edm::View<pat::Jet>> const srcToken_;
  StringCutObjectSelector<pat::Jet> const selector_;
};

DEFINE_FWK_MODULE(PATJetConstituentSelector);