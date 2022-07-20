// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TMath.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

//#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
//#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
//#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


//
// class decleration
//

#define MAXTRK 2000

using namespace std;

class TrackTreeProducer : public edm::EDAnalyzer {
public:
  explicit TrackTreeProducer(const edm::ParameterSet&);
  ~TrackTreeProducer();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&) ;
  virtual void fillGEN(const edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  virtual void initTree();

  // ----------member data ---------------------------
    
    edm::Service<TFileService> fs;

    TTree* TrackNtuple;
    
    bool   saveTree_;

    //options
    bool doRecoNtuple_;
    bool doGenNtuple_;  

    //cuts 
    double multMax_;
    double multMin_;    
    double cutPtMin_;
    double cutEtaMax_;
    double cutPtErrOverPt_;
    double cutDzOverDzSigma_;
    double cutDxyOverDxySigma_;

    //tree branches
    //event info
    float bestvx;
    float bestvy;
    float bestvz;
    
    //Track info
    int nTrk;
    float pt[MAXTRK];
    float eta[MAXTRK];
    float phi[MAXTRK];
    float trkPtErrOverPt[MAXTRK];
    float trkDzOverDzSigma[MAXTRK];
    float trkDxyOverDxySigma[MAXTRK];

    // gen info    
    int nTrk_gen;
    float pt_gen[MAXTRK];
    float eta_gen[MAXTRK];
    float phi_gen[MAXTRK];

    //tokens
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;

    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

TrackTreeProducer::TrackTreeProducer(const edm::ParameterSet& iConfig)
{
    //options
    doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
    doGenNtuple_ = iConfig.getUntrackedParameter<bool>("doGenNtuple");
    saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", 1000);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", 0);
    cutPtMin_ = iConfig.getUntrackedParameter<double>("cutPtMin", 0.3);
    cutEtaMax_ = iConfig.getUntrackedParameter<double>("cutEtaMax", 2.4);
    cutPtErrOverPt_ = iConfig.getUntrackedParameter<double>("cutPtErrOverPt", 0.1);
    cutDzOverDzSigma_ = iConfig.getUntrackedParameter<double>("cutDzOverDzSigma", 5.0);
    cutDxyOverDxySigma_ = iConfig.getUntrackedParameter<double>("cutDxyOverDxySigma", 5.0);

    //input tokens
    tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
    tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

}


TrackTreeProducer::~TrackTreeProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
TrackTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup&
iSetup)
{
    using std::vector;
    using namespace edm;
    using namespace reco;

    if(doGenNtuple_) fillGEN(iEvent,iSetup);
    if(doRecoNtuple_) fillRECO(iEvent,iSetup);

    if(saveTree_) TrackNtuple->Fill();
}

void
TrackTreeProducer::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //get collections
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);
    
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);
    
    //best vertex
    bestvz=-999.9; bestvx=-999.9; bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    //const reco::Vertex & vtx = (*vertices)[0];
    //bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    //bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();

     unsigned int daughter = 0;
     int greatestvtx = 0;
     for (unsigned int i = 0 ; i< vertices->size(); ++i){
      daughter = (*vertices)[i].tracksSize();
      if( daughter > (*vertices)[greatestvtx].tracksSize()) greatestvtx = i;
     }
     if(vertices->size()>0){
      bestvx = (*vertices)[greatestvtx].position().x();
      bestvy = (*vertices)[greatestvtx].position().y();
      bestvz = (*vertices)[greatestvtx].position().z();
      bestvxError = (*vertices)[greatestvtx].xError();
      bestvyError = (*vertices)[greatestvtx].yError();
      bestvzError = (*vertices)[greatestvtx].zError();
     }
    
    nTrk = 0;
    if(multMax_>1 && multMin_<100)
    {
      for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];
        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(fabs(trk.dzError()*trk.dzError()+bestvzError*bestvzError));
        double dxyerror = sqrt(fabs(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError));
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>cutPtErrOverPt_) continue;
        if(fabs(dzvtx/dzerror) > cutDzOverDzSigma_) continue;
        if(fabs(dxyvtx/dxyerror) > cutDxyOverDxySigma_) continue;
        if(fabs(trk.eta())>=cutEtaMax_) continue;
        if(trk.pt()<=cutPtMin_) continue;

        pt[nTrk] = trk.pt();
        eta[nTrk] = trk.eta();
        phi[nTrk] = trk.phi();
        trkPtErrOverPt[nTrk] = fabs(trk.ptError())/trk.pt();
        trkDzOverDzSigma[nTrk] = fabs(dzvtx/dzerror);
        trkDxyOverDxySigma[nTrk] = fabs(dxyvtx/dxyerror);

        nTrk++;
      }
    }
    if(nTrk==0) return;

}

void
TrackTreeProducer::fillGEN(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(tok_genParticle_,genpars);

    nTrk_gen = 0;
    for(unsigned it=0; it<genpars->size(); ++it){

        const reco::GenParticle & trk = (*genpars)[it];

        //int id = trk.pdgId();

        if(trk.charge()==0) continue;
        if(trk.status()!=1) continue;
        if(fabs(trk.eta())>=cutEtaMax_) continue;
        if(trk.pt()<=cutPtMin_) continue;
        pt_gen[nTrk_gen] = trk.pt();
        eta_gen[nTrk_gen] = trk.eta();
        phi_gen[nTrk_gen] = trk.phi();
        nTrk_gen+=1;

    }
}

// ------------ method called once each job just before starting event
//loop  ------------
void
TrackTreeProducer::beginJob()
{
    TH1D::SetDefaultSumw2();
    
    if(!doRecoNtuple_ && !doGenNtuple_)
    {
        cout<<"No output for either RECO or GEN!! Fix config!!"<<endl; return;
    }

    if(saveTree_) initTree();
}

void 
TrackTreeProducer::initTree()
{ 
    TrackNtuple = fs->make< TTree>("TrackNtuple","TrackNtuple");
    
    if(doRecoNtuple_) 
    { 
  
    // Event info
    TrackNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
    TrackNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
    TrackNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");

    // Track info
    TrackNtuple->Branch("nTrk", &nTrk, "nTrk/I"); 
    TrackNtuple->Branch("pT",&pt,"pT[nTrk]/F");
    TrackNtuple->Branch("eta",&eta,"eta[nTrk]/F");
    TrackNtuple->Branch("phi",&phi,"phi[nTrk]/F");
    TrackNtuple->Branch("trkPtErrOverPt",&trkPtErrOverPt,"trkPtErrOverPt[nTrk]/F");
    TrackNtuple->Branch("trkDzOverDzSigma",&trkDzOverDzSigma,"trkDzOverDzSigma[nTrk]/F");
    TrackNtuple->Branch("trkDxyOverDxySigma",&trkDxyOverDxySigma,"trkDxyOverDxySigma[nTrk]/F");

    } // doRecoNtuple_

    if(doGenNtuple_)
    {
        TrackNtuple->Branch("nTrk_gen",&nTrk_gen,"nTrk_gen/I");
        TrackNtuple->Branch("pt_gen",&pt_gen,"pt_gen[nTrk_gen]/F");
        TrackNtuple->Branch("eta_gen",&eta_gen,"eta_gen[nTrk_gen]/F");
        TrackNtuple->Branch("phi_gen",&phi_gen,"phi_gen[nTrk_gen]/F");

    }
}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
TrackTreeProducer::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackTreeProducer);
