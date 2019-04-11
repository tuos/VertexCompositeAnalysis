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
#include <TVector3.h>
#include <TMath.h>

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
// constants, enums and typedefs
//

#define PI 3.1416
#define MAXCAN 50000
#define MAXDAU 3
#define MAXGDAU 2
#define MAXTRG 1024
#define MAXSEL 100

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
typedef ROOT::Math::SVector<double, 6> SVector6;

//
// class decleration
//

class PATCompositeTreeProducer : public edm::EDAnalyzer {
public:
  explicit PATCompositeTreeProducer(const edm::ParameterSet&);
  ~PATCompositeTreeProducer();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(const edm::Event&, const edm::EventSetup&);
  virtual void fillGEN(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void initHistogram();
  virtual void initTree();
  reco::GenParticleRef findMother(const reco::GenParticleRef&);

  // ----------member data ---------------------------

  edm::Service<TFileService> fs;

  TTree* PATCompositeNtuple;
  TH2F*  hMassVsMVA[6][10];
  TH2F*  hpTVsMVA[6][10];
  TH2F*  hetaVsMVA[6][10];
  TH2F*  hyVsMVA[6][10];
  TH2F*  hVtxProbVsMVA[6][10];
  TH2F*  h3DCosPointingAngleVsMVA[6][10];
  TH2F*  h3DPointingAngleVsMVA[6][10];
  TH2F*  h2DCosPointingAngleVsMVA[6][10];
  TH2F*  h2DPointingAngleVsMVA[6][10];
  TH2F*  h3DDecayLengthSignificanceVsMVA[6][10];
  TH2F*  h3DDecayLengthVsMVA[6][10];
  TH2F*  h2DDecayLengthSignificanceVsMVA[6][10];
  TH2F*  h2DDecayLengthVsMVA[6][10];
  TH2F*  h3DDCAVsMVA[6][10];
  TH2F*  h2DDCAVsMVA[6][10];
  TH2F*  hzDCASignificanceDaugtherVsMVA[MAXDAU][6][10];
  TH2F*  hxyDCASignificanceDaugtherVsMVA[MAXDAU][6][10];
  TH2F*  hNHitDVsMVA[MAXDAU][6][10];
  TH2F*  hpTDVsMVA[MAXDAU][6][10];
  TH2F*  hpTerrDVsMVA[MAXDAU][6][10];
  TH2F*  hEtaDVsMVA[MAXDAU][6][10];
  TH2F*  hdedxHarmonic2DVsMVA[MAXDAU][6][10];
  TH2F*  hdedxHarmonic2DVsP[MAXDAU][6][10];

  bool   saveTree_;
  bool   saveHistogram_;
  bool   saveAllHistogram_;
  double massHistPeak_;
  double massHistWidth_;
  int    massHistBins_;

  //options
  bool doRecoNtuple_;
  bool doGenNtuple_;   
  bool doGenMatching_;
  bool doGenMatchingTOF_;
  bool hasSwap_;
  bool decayInGen_;
  bool twoLayerDecay_;
  bool threeProngDecay_;
  bool doMuon_;
  bool doMuonFull_;
  int PID_;
  const std::vector<int> PID_dau_;
  const ushort NDAU_;
  const ushort NGDAU_ = MAXGDAU;

  //cut variables
  double multMax_;
  double multMin_;
  double deltaR_; //deltaR for Gen matching

  std::vector<double> pTBins_;
  std::vector<double> yBins_;

  //tree branches
  //event info
  uint  runNb;
  uint  eventNb;
  uint  lsNb;
  short trigPrescale[MAXTRG];
  short centrality;
  int   Ntrkoffline;
  int   Npixel;
  short nPV;
  ushort candSize;
  bool  trigHLT[MAXTRG];
  bool  evtSel[MAXSEL];
  float HFsumETPlus;
  float HFsumETMinus;
  float ZDCPlus;
  float ZDCMinus;
  float bestvx;
  float bestvy;
  float bestvz;
  float ephfpAngle[3];
  float ephfmAngle[3];
  float ephfpQ[3];
  float ephfmQ[3];
  float ephfpSumW;
  float ephfmSumW;

  //Composite candidate info
  float mva[MAXCAN];
  float pt[MAXCAN];
  float eta[MAXCAN];
  float phi[MAXCAN];
  float flavor[MAXCAN];
  float y[MAXCAN];
  float mass[MAXCAN];
  float VtxProb[MAXCAN];
  float dlos[MAXCAN];
  float dl[MAXCAN];
  float dlerror[MAXCAN];
  float agl[MAXCAN];
  float vtxChi2[MAXCAN];
  float ndf[MAXCAN];
  float agl_abs[MAXCAN];
  float agl2D[MAXCAN];
  float agl2D_abs[MAXCAN];
  float dlos2D[MAXCAN];
  float dl2D[MAXCAN];
  bool isSwap[MAXCAN];
  bool matchGEN[MAXCAN];
  int idmom_reco[MAXCAN];
    
  //dau candidate info
  float grand_mass[MAXCAN];
  float grand_VtxProb[MAXCAN];
  float grand_dlos[MAXCAN];
  float grand_dl[MAXCAN];
  float grand_dlerror[MAXCAN];
  float grand_agl[MAXCAN];
  float grand_vtxChi2[MAXCAN];
  float grand_ndf[MAXCAN];
  float grand_agl_abs[MAXCAN];
  float grand_agl2D[MAXCAN];
  float grand_agl2D_abs[MAXCAN];
  float grand_dlos2D[MAXCAN];

  //dau info
  float dzos[MAXDAU][MAXCAN];
  float dxyos[MAXDAU][MAXCAN];
  float nhit[MAXDAU][MAXCAN];
  bool  trkquality[MAXDAU][MAXCAN];
  float ptDau[MAXDAU][MAXCAN];
  float ptErr[MAXDAU][MAXCAN];
  float pDau[MAXDAU][MAXCAN];
  float etaDau[MAXDAU][MAXCAN];
  float phiDau[MAXDAU][MAXCAN];
  short chargeDau[MAXDAU][MAXCAN];
  int   pid[MAXDAU][MAXCAN];
  float tof[MAXDAU][MAXCAN];
  float H2dedx[MAXDAU][MAXCAN];
  float T4dedx[MAXDAU][MAXCAN];
  float trkChi[MAXDAU][MAXCAN];
   
  //grand-dau info
  float grand_dzos[MAXGDAU][MAXCAN];
  float grand_dxyos[MAXGDAU][MAXCAN];
  float grand_nhit[MAXGDAU][MAXCAN];
  bool  grand_trkquality[MAXGDAU][MAXCAN];
  float grand_pt[MAXGDAU][MAXCAN];
  float grand_ptErr[MAXGDAU][MAXCAN];
  float grand_p[MAXGDAU][MAXCAN];
  float grand_eta[MAXGDAU][MAXCAN];
  short grand_charge[MAXGDAU][MAXCAN];
  float grand_H2dedx[MAXGDAU][MAXCAN];
  float grand_T4dedx[MAXGDAU][MAXCAN];
  float grand_trkChi[MAXGDAU][MAXCAN];
    
  //dau muon info
  bool  onestmuon[MAXDAU][MAXCAN];
  bool  pfmuon[MAXDAU][MAXCAN];
  bool  glbmuon[MAXDAU][MAXCAN];
  bool  trkmuon[MAXDAU][MAXCAN];
  bool  tightmuon[MAXDAU][MAXCAN];
  bool  softmuon[MAXDAU][MAXCAN];
  bool  hybridmuon[MAXDAU][MAXCAN];
  bool  hpmuon[MAXDAU][MAXCAN];
  bool  trgmuon[MAXDAU][MAXTRG][MAXCAN];
  short nmatchedst[MAXDAU][MAXCAN];
  short ntrackerlayer[MAXDAU][MAXCAN];
  short npixellayer[MAXDAU][MAXCAN];
  short npixelhit[MAXDAU][MAXCAN];
  short nmuonhit[MAXDAU][MAXCAN];
  float glbtrkchi[MAXDAU][MAXCAN];
  float muonbestdxy[MAXDAU][MAXCAN];
  float muonbestdz[MAXDAU][MAXCAN];
  float muondxy[MAXDAU][MAXCAN];
  float muondz[MAXDAU][MAXCAN];
  short nmatchedch[MAXDAU][MAXCAN];
  float matchedenergy[MAXDAU][MAXCAN];
  float dx_seg[MAXDAU][MAXCAN];
  float dy_seg[MAXDAU][MAXCAN];
  float dxSig_seg[MAXDAU][MAXCAN];
  float dySig_seg[MAXDAU][MAXCAN];
  float ddxdz_seg[MAXDAU][MAXCAN];
  float ddydz_seg[MAXDAU][MAXCAN];
  float ddxdzSig_seg[MAXDAU][MAXCAN];
  float ddydzSig_seg[MAXDAU][MAXCAN];

  // gen info
  float weight_gen;
  ushort candSize_gen;
  float pt_gen[MAXCAN];
  float eta_gen[MAXCAN];
  short status_gen[MAXCAN];
  int idmom[MAXCAN];
  float y_gen[MAXCAN];
  int iddau[MAXDAU][MAXCAN];
  short chargedau[MAXDAU][MAXCAN];
  float ptdau[MAXDAU][MAXCAN];
  float etadau[MAXDAU][MAXCAN];
  float phidau[MAXDAU][MAXCAN];
  float massdau[MAXDAU][MAXCAN];
    
  bool useAnyMVA_;
  bool isSkimMVA_;
  bool isCentrality_;
  bool isEventPlane_;
  bool useDeDxData_;

  //tokens
  edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> patCompositeCandidateCollection_Token_;
  edm::EDGetTokenT<MVACollection> MVAValues_Token_;

  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
  edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
  edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
  edm::EDGetTokenT<GenEventInfoProduct> tok_genInfo_;

  edm::EDGetTokenT<int> tok_centBinLabel_;
  edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

  edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;

  //trigger
  const std::vector<std::string> triggerNames_;
  const std::vector<std::string> filterNames_;
  edm::EDGetTokenT<edm::TriggerResults> tok_triggerResults_;
  const ushort NTRG_;

  //event selection
  const std::vector<std::string> eventFilters_;
  edm::EDGetTokenT<edm::TriggerResults> tok_filterResults_;
  const ushort NSEL_;

  //prescale provider
  HLTPrescaleProvider hltPrescaleProvider_;
};

//
// static data member definitions
//

//
// constructors and destructor
//

PATCompositeTreeProducer::PATCompositeTreeProducer(const edm::ParameterSet& iConfig) :
  PID_dau_(iConfig.getUntrackedParameter<std::vector<int> >("PID_dau")),
  NDAU_(PID_dau_.size()>MAXDAU ? MAXDAU : PID_dau_.size()),
  patCompositeCandidateCollection_Token_(consumes<pat::CompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"))),
  triggerNames_(iConfig.getUntrackedParameter<std::vector<std::string> >("triggerPathNames")),
  filterNames_(iConfig.getUntrackedParameter<std::vector<std::string> >("triggerFilterNames")),
  NTRG_(triggerNames_.size()>MAXTRG ? MAXTRG : triggerNames_.size()),
  eventFilters_(iConfig.getUntrackedParameter<std::vector<std::string> >("eventFilterNames")),
  NSEL_(eventFilters_.size()>MAXSEL ? MAXSEL : eventFilters_.size()),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
  //options
  doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
  doGenNtuple_ = iConfig.getUntrackedParameter<bool>("doGenNtuple");
  twoLayerDecay_ = iConfig.getUntrackedParameter<bool>("twoLayerDecay");
  threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
  doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
  doGenMatchingTOF_ = iConfig.getUntrackedParameter<bool>("doGenMatchingTOF");
  hasSwap_ = iConfig.getUntrackedParameter<bool>("hasSwap");
  decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
  doMuon_ = iConfig.getUntrackedParameter<bool>("doMuon");
  doMuonFull_ = iConfig.getUntrackedParameter<bool>("doMuonFull");
  PID_ = iConfig.getUntrackedParameter<int>("PID");

  saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
  saveHistogram_ = iConfig.getUntrackedParameter<bool>("saveHistogram");
  saveAllHistogram_ = iConfig.getUntrackedParameter<bool>("saveAllHistogram");
  massHistPeak_ = iConfig.getUntrackedParameter<double>("massHistPeak");
  massHistWidth_ = iConfig.getUntrackedParameter<double>("massHistWidth");
  massHistBins_ = iConfig.getUntrackedParameter<int>("massHistBins");

  //cut variables
  multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
  multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
  deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.03);

  pTBins_ = iConfig.getUntrackedParameter<std::vector<double> >("pTBins");
  yBins_  = iConfig.getUntrackedParameter<std::vector<double> >("yBins");

  //input tokens
  tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
  MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
  tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));
  tok_genInfo_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

  useDeDxData_ = (iConfig.exists("useDeDxData") ? iConfig.getParameter<bool>("useDeDxData") : false);
  if(useDeDxData_)
  {
    Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
    Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));
  }

  isCentrality_ = (iConfig.exists("isCentrality") ? iConfig.getParameter<bool>("isCentrality") : false);
  if(isCentrality_)
  {
    tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
    tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  }

  isEventPlane_ = (iConfig.exists("isEventPlane") ? iConfig.getParameter<bool>("isEventPlane") : false);
  if(isEventPlane_) tok_eventplaneSrc_ = consumes<reco::EvtPlaneCollection>(iConfig.getParameter<edm::InputTag>("eventplaneSrc"));

  useAnyMVA_ = (iConfig.exists("useAnyMVA") ? iConfig.getParameter<bool>("useAnyMVA") : false);
  if(useAnyMVA_) MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
  isSkimMVA_ = iConfig.getUntrackedParameter<bool>("isSkimMVA");
  
  tok_triggerResults_ = consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultCollection"));
  tok_filterResults_ = consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("FilterResultCollection"));
}


PATCompositeTreeProducer::~PATCompositeTreeProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
PATCompositeTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(doGenNtuple_) fillGEN(iEvent,iSetup);
  if(doRecoNtuple_) fillRECO(iEvent,iSetup);
  if(saveTree_) PATCompositeNtuple->Fill();
}


void
PATCompositeTreeProducer::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Check inputs
  if((!threeProngDecay_ && NDAU_!=2) || (threeProngDecay_ && NDAU_!=3))
  {
    edm::LogError("PATCompositeAnalyzer") << "Want threeProngDecay but PID daughter vector size is: " << NDAU_ << " !" << std::endl;
    return;
  }

  //get collections
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);
  if(!vertices.isValid()) { edm::LogError("PATCompositeAnalyzer") << "Primary vertices  collection not found!" << std::endl; return; }

  edm::Handle<pat::CompositeCandidateCollection> v0candidates;
  iEvent.getByToken(patCompositeCandidateCollection_Token_, v0candidates);
  if(!v0candidates.isValid()) { edm::LogError("PATCompositeAnalyzer") << "V0 candidate collection not found!" << std::endl; return; }

  edm::Handle<MVACollection> mvavalues;
  if(useAnyMVA_)
  {
    iEvent.getByToken(MVAValues_Token_, mvavalues);
    if(!mvavalues.isValid()) { edm::LogError("PATCompositeAnalyzer") << "MVA collection not found!" << std::endl; return; }
    assert( (*mvavalues).size() == v0candidates->size() );
  }

  edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1 , dEdxHandle2;
  if(useDeDxData_)
  {
    iEvent.getByToken(Dedx_Token1_, dEdxHandle1);
    iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
  }

  runNb = iEvent.id().run();
  eventNb = iEvent.id().event();
  lsNb = iEvent.luminosityBlock();

  //Trigger Information
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(tok_triggerResults_, triggerResults);
  if(triggerNames_.size()>0)
  {
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
    for (ushort iTr=0; iTr<triggerNames_.size(); iTr++) {
      trigHLT[iTr] = false;
      const auto& triggerIndex = triggerNames.triggerIndex(triggerNames_.at(iTr));
      if(triggerIndex>=triggerNames.size()) continue;
      bool isTriggerFired = false;
      if(triggerResults->accept(triggerIndex)) isTriggerFired = true;
      int prescaleValue = -1;
      if (hltPrescaleProvider_.hltConfigProvider().inited() && hltPrescaleProvider_.prescaleSet(iEvent,iSetup)>=0)
      {
        const auto& presInfo = hltPrescaleProvider_.prescaleValuesInDetail(iEvent, iSetup, triggerNames_.at(iTr));
        const auto& hltPres = presInfo.second;
        const short& l1Pres = ((presInfo.first.size()==1) ? presInfo.first.at(0).second : ((presInfo.first.size()>1) ? 1 : -1));
        prescaleValue = hltPres*l1Pres;
      }
      trigPrescale[iTr] = prescaleValue;
      if(isTriggerFired) trigHLT[iTr] = true;
    }
  }

  //Event selection information
  edm::Handle<edm::TriggerResults> filterResults;
  iEvent.getByToken(tok_filterResults_, filterResults);
  if(eventFilters_.size()>0)
  {
    const edm::TriggerNames& filterNames = iEvent.triggerNames(*filterResults);
    for(ushort iFr=0; iFr<eventFilters_.size(); ++iFr)
    {
      evtSel[iFr] = false;
      const auto& index = filterNames.triggerIndex(eventFilters_.at(iFr));
      if(index < filterNames.size()) evtSel[iFr] = filterResults->accept(index);
    }
  }

  centrality = -1;
  if(isCentrality_)
  {
    edm::Handle<reco::Centrality> cent;
    iEvent.getByToken(tok_centSrc_, cent);
    HFsumETPlus = (cent.isValid() ? cent->EtHFtowerSumPlus() : -1.);
    HFsumETMinus = (cent.isValid() ? cent->EtHFtowerSumMinus() : -1.);
    Npixel = (cent.isValid() ? cent->multiplicityPixel() : -1);
    ZDCPlus = (cent.isValid() ? cent->zdcSumPlus() : -1.);
    ZDCMinus = (cent.isValid() ? cent->zdcSumMinus() : -1.);
    Ntrkoffline = (cent.isValid() ? cent->Ntracks() : -1);
      
    edm::Handle<int> cbin;
    iEvent.getByToken(tok_centBinLabel_, cbin);
    if(!cbin.isValid()) { edm::LogError("PATCompositeAnalyzer") << "Centrality bin not valid!" << std::endl; return; }
    centrality = *cbin;
  }

  if(isEventPlane_)
  {
    edm::Handle<reco::EvtPlaneCollection> eventplanes;
    iEvent.getByToken(tok_eventplaneSrc_, eventplanes);
    if(!eventplanes.isValid()) { edm::LogError("PATCompositeAnalyzer") << "Event plane collection not valid!" << std::endl; return; }

    const reco::EvtPlane& ephfp1 = (*eventplanes)[0];
    const reco::EvtPlane& ephfm1 = (*eventplanes)[1];
    const reco::EvtPlane& ephfp2 = (*eventplanes)[6];
    const reco::EvtPlane& ephfm2 = (*eventplanes)[7];
    const reco::EvtPlane& ephfp3 = (*eventplanes)[13];
    const reco::EvtPlane& ephfm3 = (*eventplanes)[14];

    ephfpAngle[0] = ephfp1.angle(2);
    ephfpAngle[1] = ephfp2.angle(2);
    ephfpAngle[2] = ephfp3.angle(2);

    ephfmAngle[0] = ephfm1.angle(2);
    ephfmAngle[1] = ephfm2.angle(2);
    ephfmAngle[2] = ephfm3.angle(2);

    ephfpQ[0] = ephfp1.q(2);
    ephfpQ[1] = ephfp2.q(2);
    ephfpQ[2] = ephfp3.q(2);

    ephfmQ[0] = ephfm1.q(2);
    ephfmQ[1] = ephfm2.q(2);
    ephfmQ[2] = ephfm3.q(2);

    ephfpSumW = ephfp2.sumw();
    ephfmSumW = ephfm2.sumw();
  }

  nPV = vertices->size();
  //best vertex
  const reco::Vertex& vtx = (*vertices)[0];
  bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
  const math::XYZPoint bestvtx(bestvx, bestvy, bestvz);
  const double& bestvzError = vtx.zError(), bestvxError = vtx.xError(), bestvyError = vtx.yError();

  //Gen info for matching
  std::vector< std::vector< std::vector<double> > > pVectDau;
  std::vector<int> pVectIDmom;

  if(doGenMatching_)
  {
    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(tok_genParticle_, genpars);
    if(!genpars.isValid()) { edm::LogError("PATCompositeAnalyzer") << "Gen matching cannot be done without Gen collection!" << std::endl; return; }

    for(uint idx=0; idx<genpars->size(); ++idx)
    {
      const auto& trk = reco::GenParticleRef(genpars, idx);

      if (trk.isNull()) continue; //check gen particle ref
      if(fabs(trk->pdgId())!=PID_) continue; //check is target
      if(trk->status()!=2) continue; //check status

      const ushort& nDau = trk->numberOfDaughters();
      if (nDau!=NDAU_) continue;

      bool skipGen = false;
      std::vector<int> PIDvec = PID_dau_;
      for(ushort iDau=0; iDau<nDau; iDau++)
      {
        const auto& Dd = *(trk->daughter(iDau));
        bool found = false;
        for(ushort iPID=0; iPID<PIDvec.size(); iPID++)
        {
          if(fabs(Dd.pdgId())==PIDvec[iPID])
          {
            PIDvec.erase(PIDvec.begin()+iPID);
            found=true;
            break;
          }
        }
        if(!found) { skipGen = (PIDvec.size()>0); break; }
      }
      if(skipGen) continue;

      const auto& mom = findMother(trk);
      const int& idmom_tmp = (mom->pdgId()!=trk->pdgId() ? mom->pdgId() : -77);

      std::vector< std::vector<double> > pVect;
      for(ushort iDau=0; iDau<nDau; iDau++)
      {
        const auto& Dd = *(trk->daughter(iDau));
        std::vector<double> Dvector;
        Dvector.push_back(Dd.pt());
        Dvector.push_back(Dd.eta());
        Dvector.push_back(Dd.phi());
        Dvector.push_back(Dd.charge());
        Dvector.push_back(Dd.mass());
        pVect.push_back(Dvector);
      }
      pVectDau.push_back(pVect);
      pVectIDmom.push_back(idmom_tmp);
    }
  }

  //RECO Candidate info
  candSize = v0candidates->size();
  for(ushort it=0; it<candSize; ++it)
  { 
    const auto& trk = (*v0candidates)[it];

    pt[it] = trk.pt();
    eta[it] = trk.eta();
    phi[it] = trk.phi();
    mass[it] = trk.mass();
    y[it] = trk.rapidity();
    flavor[it] = trk.pdgId()/fabs(trk.pdgId());

    mva[it] = 0.0;
    if(useAnyMVA_) mva[it] = (*mvavalues)[it];

    //vtxChi2
    vtxChi2[it] = trk.userFloat("vertexChi2");
    ndf[it] = trk.userFloat("vertexNdof");
    VtxProb[it] = TMath::Prob(vtxChi2[it],ndf[it]);

    const double& secvz = trk.vz(), secvx = trk.vx(), secvy = trk.vy();
    const double& px = trk.px(), py = trk.py(), pz = trk.pz();

    //PAngle
    const TVector3 ptosvec(secvx-bestvx, secvy-bestvy, secvz-bestvz);
    const TVector3 secvec(px, py, pz);
    const TVector3 ptosvec2D(secvx-bestvx, secvy-bestvy, 0);
    const TVector3 secvec2D(px,py,0);

    agl[it] = std::cos(secvec.Angle(ptosvec));
    agl_abs[it] = secvec.Angle(ptosvec);
    agl2D[it] = std::cos(secvec2D.Angle(ptosvec2D));
    agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
        
    //Decay length 3D
    const SMatrixSym3D& trkCovMat = *trk.userData<reco::Vertex::CovarianceMatrix>("vertexCovariance");
    const SMatrixSym3D& totalCov = vtx.covariance() + trkCovMat;
    const SVector3 distanceVector(secvx-bestvx, secvy-bestvy, secvz-bestvz);

    dl[it] = ROOT::Math::Mag(distanceVector);
    dlerror[it] = std::sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl[it];
    dlos[it] = dl[it]/dlerror[it];

    //Decay length 2D
    const SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1), vtx.covariance(1,1), 0, 0, 0);
    const SVector6 v2(trkCovMat(0,0), trkCovMat(0,1), trkCovMat(1,1), 0, 0, 0);
    const SMatrixSym3D& totalCov2D = SMatrixSym3D(v1) + SMatrixSym3D(v2);
    const SVector3 distanceVector2D(secvx-bestvx, secvy-bestvy, 0);

    dl2D[it] = ROOT::Math::Mag(distanceVector2D);
    const double& dl2Derror = std::sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D[it];
    dlos2D[it] = dl2D[it]/dl2Derror;

    const ushort& nDau = trk.numberOfDaughters();
    if(nDau!=NDAU_) { edm::LogError("PATCompositeAnalyzer") << "Expected " << NDAU_ << " daughters but V0 candidate has " << nDau << " daughters!" << std::endl; return; }

    //Gen match
    if(doGenMatching_)
    {
      matchGEN[it] = false;
      isSwap[it] = false;
      idmom_reco[it] = -77;
      for(ushort iGen=0; iGen<pVectDau.size(); iGen++)
      {
        bool foundGen = false;
        std::vector<double> dmassGEN;
        auto pVect = pVectDau[iGen];
        for(ushort iDau=0; iDau<nDau; iDau++)
        {
          const auto& Dd = *(trk.daughter(iDau));
          bool found = false;
          for(ushort iDG=0; iDG<pVect.size(); iDG++)
          {
            const auto& Dvector = pVect[iDG];
            if(Dd.charge()!=Dvector.at(3)) continue;
            if(fabs((Dd.pt()-Dvector.at(0))/Dd.pt()) > 0.5) continue; //check deltaPt matching
            const double& deltaR = reco::deltaR(Dd.eta(), Dd.phi(), Dvector.at(1), Dvector.at(2));
            if(deltaR <= deltaR_)
            {
              dmassGEN.push_back(Dvector.at(4));
              pVect.erase(pVect.begin()+iDG);
              found=true;
              break;
            }
          }
          if(found) { foundGen = true; break; }
        }
        if(foundGen)
        {
          matchGEN[it] = true; //matched gen
          //check swap
          for(ushort iDau=0; iDau<nDau; iDau++)
          {
            if(fabs(dmassGEN[iDau] - trk.daughter(iDau)->mass())>0.01) { isSwap[it] = true; break; }
          }
          //check prompt & record mom id
          idmom_reco[it] = pVectIDmom[iGen];
          break;
        }
      }
    }

    for(ushort iDau=0; iDau<nDau; iDau++)
    {
      const auto& dau = *(trk.daughter(iDau));

      ptDau[iDau][it] = dau.pt();
      pDau[iDau][it] = dau.p();
      etaDau[iDau][it] = dau.eta();
      phiDau[iDau][it] = dau.phi();
      chargeDau[iDau][it] = dau.charge();

      pid[iDau][it] = -99999;
      if(doGenMatchingTOF_)
      {
        const auto& ppDau = dynamic_cast<const pat::PATObject<reco::Candidate>*>(trk.daughter(iDau));
        const auto& trk = (ppDau ? ppDau->genParticleRef() : reco::GenParticleRef());
        if(trk.isNonnull()) { pid[iDau][it] = trk->pdgId(); }
      }

      //trk info
      if(!twoLayerDecay_)
      {
        const auto& dtrk = dau.get<reco::TrackRef>();

        //trk quality
        trkquality[iDau][it] = (dtrk.isNonnull() ? dtrk->quality(reco::TrackBase::highPurity) : false);

        //trk dEdx
        H2dedx[iDau][it] = -999.9;
        if(dtrk.isNonnull() && dEdxHandle1.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle1.product();
          H2dedx[iDau][it] = dEdxTrack[dtrk].dEdx();
        }
        T4dedx[iDau][it] = -999.9;
        if(dtrk.isNonnull() && dEdxHandle2.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle2.product();
          T4dedx[iDau][it] = dEdxTrack[dtrk].dEdx();
        }

        //track Chi2
        trkChi[iDau][it] = (dtrk.isNonnull() ? dtrk->normalizedChi2() : 99.);

        //track pT error
        ptErr[iDau][it] = (dtrk.isNonnull() ? dtrk->ptError() : -1.);

        //trkNHits
        nhit[iDau][it] = (dtrk.isNonnull() ? dtrk->numberOfValidHits() : -1);

        //DCA
        dzos[iDau][it] = 99.;
        dxyos[iDau][it] = 99.;
        if (dtrk.isNonnull())
        {
          const double& dzbest = dtrk->dz(bestvtx);
          const double& dxybest = dtrk->dxy(bestvtx);
          const double& dzerror = std::sqrt(dtrk->dzError()*dtrk->dzError() + bestvzError*bestvzError);
          const double& dxyerror = std::sqrt(dtrk->d0Error()*dtrk->d0Error() + bestvxError*bestvyError);
          dzos[iDau][it] = dzbest/dzerror;
          dxyos[iDau][it] = dxybest/dxyerror;
        }
      }
 
      if(doMuon_)
      {
        const auto& muon = (dau.isMuon() ? *dynamic_cast<const pat::Muon*>(trk.daughter(iDau)) : pat::Muon());

        // Tight ID Muon POG Run 2
        glbmuon[iDau][it] = (dau.isMuon() ? muon.isGlobalMuon() : false);
        pfmuon[iDau][it]  = (dau.isMuon() ? muon.isPFMuon() : false);
        glbtrkchi[iDau][it] = (muon.globalTrack().isNonnull() ? muon.globalTrack()->normalizedChi2() : 99.);
        nmuonhit[iDau][it] = (muon.globalTrack().isNonnull() ? muon.globalTrack()->hitPattern().numberOfValidMuonHits() : -1);
        nmatchedst[iDau][it] = (dau.isMuon() ? muon.numberOfMatchedStations() : -1);
        npixelhit[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().numberOfValidPixelHits() : -1);
        ntrackerlayer[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() : -1);
        muonbestdxy[iDau][it] = (muon.muonBestTrack().isNonnull() ? muon.muonBestTrack()->dxy(bestvtx) : 99.);
        muonbestdz[iDau][it] = (muon.muonBestTrack().isNonnull() ? muon.muonBestTrack()->dz(bestvtx) : 99.);
        tightmuon[iDau][it] = (
                               glbmuon[iDau][it] &&
                               pfmuon[iDau][it] &&
                               (glbtrkchi[iDau][it] < 10.) &&
                               (nmuonhit[iDau][it] > 0) &&
                               (nmatchedst[iDau][it] > 1) &&
                               (npixelhit[iDau][it] > 0) &&
                               (ntrackerlayer[iDau][it] > 5) &&
                               (fabs(muonbestdxy[iDau][it]) < 0.2) &&
                               (fabs(muonbestdz[iDau][it]) < 0.5)
                               );

        // Soft ID Muon POG Run 2
        onestmuon[iDau][it] = (dau.isMuon() ? muon::isGoodMuon(muon, muon::SelectionType::TMOneStationTight) : false);
        npixellayer[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->hitPattern().pixelLayersWithMeasurement() : -1);
        hpmuon[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->quality(reco::TrackBase::highPurity) : false);
        muondxy[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->dxy(bestvtx) : 99.);
        muondz[iDau][it] = (muon.innerTrack().isNonnull() ? muon.innerTrack()->dz(bestvtx) : 99.);
        softmuon[iDau][it] = (
                              onestmuon[iDau][it] &&
                              (ntrackerlayer[iDau][it] > 5) &&
                              (npixellayer[iDau][it] > 0) &&
                              hpmuon[iDau][it] &&
                              (fabs(muondxy[iDau][it]) < 0.3) &&
                              (fabs(muondz[iDau][it]) < 20.)
                              );

        // Hybrid Soft ID HIN PAG Run 2 PbPb
        trkmuon[iDau][it] = (dau.isMuon() ? muon.isTrackerMuon() : false);
        hybridmuon[iDau][it] = (
                                glbmuon[iDau][it] &&
                                (ntrackerlayer[iDau][it] > 5) &&
                                (npixellayer[iDau][it] > 0) &&
                                (fabs(muondxy[iDau][it]) < 0.3) &&
                                (fabs(muondz[iDau][it]) < 20.)
                                );

        // Muon Trigger Matching
        if(dau.isMuon())
        {
          for(ushort iTr=0; iTr<filterNames_.size(); iTr++)
          {
            trgmuon[iDau][iTr][it] = false;
            const auto& muHLTMatchesFilter = muon.triggerObjectMatchesByFilter(filterNames_.at(iTr));
            if(muHLTMatchesFilter.size()>0) trgmuon[iDau][iTr][it] = true;
          }
        }

        if(doMuonFull_)
        {
          nmatchedch[iDau][it] = (dau.isMuon() ? muon.numberOfMatches() : -1);
          matchedenergy[iDau][it] = (dau.isMuon() ? muon.calEnergy().hadMax : -99.);

          dx_seg[iDau][it] = 999.9;
          dy_seg[iDau][it] = 999.9;
          dxSig_seg[iDau][it] = 999.9;
          dySig_seg[iDau][it] = 999.9;
          ddxdz_seg[iDau][it] = 999.9;
          ddydz_seg[iDau][it] = 999.9;
          ddxdzSig_seg[iDau][it] = 999.9;
          ddydzSig_seg[iDau][it] = 999.9;
          const std::vector<reco::MuonChamberMatch>& muchmatches = muon.matches();
          for(ushort ich=0; ich<muchmatches.size(); ich++)
          {
            const double& x_exp = muchmatches[ich].x;
            const double& y_exp = muchmatches[ich].y;
            const double& xerr_exp = muchmatches[ich].xErr;
            const double& yerr_exp = muchmatches[ich].yErr;
            const double& dxdz_exp = muchmatches[ich].dXdZ;
            const double& dydz_exp = muchmatches[ich].dYdZ;
            const double& dxdzerr_exp = muchmatches[ich].dXdZErr;
            const double& dydzerr_exp = muchmatches[ich].dYdZErr;

            const std::vector<reco::MuonSegmentMatch>& musegmatches = muchmatches[ich].segmentMatches;
            for(ushort jseg=0; jseg<musegmatches.size(); jseg++)
            {
              const double& x_seg = musegmatches[jseg].x;
              const double& y_seg = musegmatches[jseg].y;
              const double& xerr_seg = musegmatches[jseg].xErr;
              const double& yerr_seg = musegmatches[jseg].yErr;
              const double& dxdz_seg = musegmatches[jseg].dXdZ;
              const double& dydz_seg = musegmatches[jseg].dYdZ;
              const double& dxdzerr_seg = musegmatches[jseg].dXdZErr;
              const double& dydzerr_seg = musegmatches[jseg].dYdZErr;

              const double& dseg = std::sqrt((x_seg-x_exp)*(x_seg-x_exp) + (y_seg-y_exp)*(y_seg-y_exp));
              const double& dxerr_seg = std::sqrt(xerr_seg*xerr_seg + xerr_exp*xerr_exp);
              const double& dyerr_seg = std::sqrt(yerr_seg*yerr_seg + yerr_exp*yerr_exp);
              const double& ddxdzerr_seg = std::sqrt(dxdzerr_seg*dxdzerr_seg + dxdzerr_exp*dxdzerr_exp);
              const double& ddydzerr_seg = std::sqrt(dydzerr_seg*dydzerr_seg + dydzerr_exp*dydzerr_exp);

              if(dseg < std::sqrt(dx_seg[iDau][it]*dx_seg[iDau][it] + dy_seg[iDau][it]*dy_seg[iDau][it]))
              {
                dx_seg[iDau][it] = x_seg - x_exp;
                dy_seg[iDau][it] = y_seg - y_exp;
                dxSig_seg[iDau][it] = dx_seg[iDau][it] / dxerr_seg;
                dySig_seg[iDau][it] = dy_seg[iDau][it] / dyerr_seg;
                ddxdz_seg[iDau][it] = dxdz_seg - dxdz_exp;
                ddydz_seg[iDau][it] = dydz_seg - dydz_exp;
                ddxdzSig_seg[iDau][it] = ddxdz_seg[iDau][it] / ddxdzerr_seg;
                ddydzSig_seg[iDau][it] = ddydz_seg[iDau][it] / ddydzerr_seg;
              }
            }
          }
        }
      }
    }
 
    if(twoLayerDecay_)
    {
      const auto& d = *(trk.daughter(0));
      grand_mass[it] = d.mass();
      for(ushort iGDau=0; iGDau<NGDAU_; iGDau++)
      {
        if(!d.daughter(iGDau)) continue;
        const auto& gd = *(d.daughter(iGDau));
        const auto& gdau = gd.get<reco::TrackRef>();

        //trk quality
        grand_trkquality[iGDau][it] = (gdau.isNonnull() ? gdau->quality(reco::TrackBase::highPurity) : false);

        //trk dEdx
        grand_H2dedx[iGDau][it] = -999.9;
        if(gdau.isNonnull() && dEdxHandle1.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle1.product();
          grand_H2dedx[iGDau][it] = dEdxTrack[gdau].dEdx();
        }   
        grand_T4dedx[iGDau][it] = -999.9;
        if(gdau.isNonnull() && dEdxHandle2.isValid())
        {
          const edm::ValueMap<reco::DeDxData>& dEdxTrack = *dEdxHandle2.product();
          grand_T4dedx[iGDau][it] = dEdxTrack[gdau].dEdx();
        }

        //track pt
        grand_pt[iGDau][it] = gd.pt();

        //track momentum
        grand_p[iGDau][it] = gd.p();

        //track eta
        grand_eta[iGDau][it] = gd.eta();

        //track charge
        grand_charge[iGDau][it] = gd.charge();

        //track Chi2
        grand_trkChi[iGDau][it] = (gdau.isNonnull() ? gdau->normalizedChi2() : 99.);

        //track pT error
        grand_ptErr[iGDau][it] = (gdau.isNonnull() ? gdau->ptError() : -1.);

        //trkNHits
        grand_nhit[iGDau][it] = (gdau.isNonnull() ? gdau->numberOfValidHits() : -1);

        //DCA
        grand_dzos[iGDau][it] = 99.;
        grand_dxyos[iGDau][it] = 99.;
        if(gdau.isNonnull())
        {
          const double& gdzbest = gdau->dz(bestvtx);
          const double& gdxybest = gdau->dxy(bestvtx);
          const double& gdzerror = std::sqrt(gdau->dzError()*gdau->dzError() + bestvzError*bestvzError);
          const double& gdxyerror = std::sqrt(gdau->d0Error()*gdau->d0Error() + bestvxError*bestvyError);
          grand_dzos[iGDau][it] = gdzbest/gdzerror;
          grand_dxyos[iGDau][it] = gdxybest/gdxyerror;
        }
      }
   
      //vtxChi2
      grand_vtxChi2[it] = d.vertexChi2();
      grand_ndf[it] = d.vertexNdof();
      grand_VtxProb[it] = TMath::Prob(grand_vtxChi2[it], grand_ndf[it]);

      //PAngle
      const double& secvz = d.vz(), secvx = d.vx(), secvy = d.vy();
      const TVector3 ptosvec(secvx-bestvx, secvy-bestvy, secvz-bestvz);
      const TVector3 secvec(d.px(), d.py(), d.pz());            
      const TVector3 ptosvec2D(secvx-bestvx, secvy-bestvy, 0);
      const TVector3 secvec2D(d.px(), d.py(), 0);

      grand_agl[it] = std::cos(secvec.Angle(ptosvec));
      grand_agl_abs[it] = secvec.Angle(ptosvec);
      grand_agl2D[it] = std::cos(secvec2D.Angle(ptosvec2D));
      grand_agl2D_abs[it] = secvec2D.Angle(ptosvec2D);

      //Decay length 3D
      const SMatrixSym3D& totalCov = vtx.covariance() + d.vertexCovariance();
      const SVector3 distanceVector(secvx-bestvx, secvy-bestvy, secvz-bestvz);

      grand_dl[it] = ROOT::Math::Mag(distanceVector);
      grand_dlerror[it] = std::sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/grand_dl[it];
      grand_dlos[it] = grand_dl[it]/grand_dlerror[it];

      //Decay length 2D
      const SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1), vtx.covariance(1,1), 0, 0, 0);
      const SVector6 v2(d.vertexCovariance(0,0), d.vertexCovariance(0,1), d.vertexCovariance(1,1), 0, 0, 0);
      const SMatrixSym3D totalCov2D = SMatrixSym3D(v1) + SMatrixSym3D(v2);
      const SVector3 distanceVector2D(secvx-bestvx, secvy-bestvy, 0);

      const double& gdl2D = ROOT::Math::Mag(distanceVector2D);
      const double& gdl2Derror = std::sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/gdl2D;
      grand_dlos2D[it] = gdl2D/gdl2Derror;
    }

    if(saveHistogram_)
    {
      for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
      {
        for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
        {
          if(pt[it]<pTBins_[ipt+1] && pt[it]>pTBins_[ipt] && y[it]<yBins_[iy+1] && y[it]>yBins_[iy])
          {
            hMassVsMVA[iy][ipt]->Fill(mva[it], mass[it]);

            if(saveAllHistogram_)
            {
              hpTVsMVA[iy][ipt]->Fill(mva[it], pt[it]);
              hetaVsMVA[iy][ipt]->Fill(mva[it], eta[it]);
              hyVsMVA[iy][ipt]->Fill(mva[it], y[it]);
              hVtxProbVsMVA[iy][ipt]->Fill(mva[it], VtxProb[it]);
              h3DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl[it]);
              h3DPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl_abs[it]);
              h2DCosPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl2D[it]);
              h2DPointingAngleVsMVA[iy][ipt]->Fill(mva[it], agl2D_abs[it]);
              h3DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it], dlos[it]);
              h3DDecayLengthVsMVA[iy][ipt]->Fill(mva[it], dl[it]);
              h2DDecayLengthSignificanceVsMVA[iy][ipt]->Fill(mva[it], dlos2D[it]);
              h2DDecayLengthVsMVA[iy][ipt]->Fill(mva[it], dl2D[it]);
              for (ushort iDau=0; iDau<NDAU_; iDau++)
              {
                hzDCASignificanceDaugtherVsMVA[iDau][iy][ipt]->Fill(mva[it], dzos[iDau][it]);
                hxyDCASignificanceDaugtherVsMVA[iDau][iy][ipt]->Fill(mva[it], dxyos[iDau][it]);
                hNHitDVsMVA[iDau][iy][ipt]->Fill(mva[it], nhit[iDau][it]);
                hpTDVsMVA[iDau][iy][ipt]->Fill(mva[it], ptDau[iDau][it]);
                hpTerrDVsMVA[iDau][iy][ipt]->Fill(mva[it], ptErr[iDau][it]/ptDau[iDau][it]);
                hEtaDVsMVA[iDau][iy][ipt]->Fill(mva[it], etaDau[iDau][it]);
                hdedxHarmonic2DVsMVA[iDau][iy][ipt]->Fill(mva[it], H2dedx[iDau][it]);
                hdedxHarmonic2DVsP[iDau][iy][ipt]->Fill(pDau[iDau][it], H2dedx[iDau][it]);
              }
            }
          }
        }
      }
    }
  }
}


void
PATCompositeTreeProducer::fillGEN(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<GenEventInfoProduct> geninfo;
  iEvent.getByToken(tok_genInfo_, geninfo);
  weight_gen = (geninfo.isValid() ? geninfo->weight() : -1.0);

  edm::Handle<reco::GenParticleCollection> genpars;
  iEvent.getByToken(tok_genParticle_, genpars);
  if(!genpars.isValid()) { edm::LogError("PATCompositeAnalyzer") << "Gen matching cannot be done without Gen collection!" << std::endl; return; }

  candSize_gen = 0;
  for(uint idx=0; idx<genpars->size(); ++idx)
  {
    const auto& trk = reco::GenParticleRef(genpars, idx);

    if (trk.isNull()) continue; //check gen particle ref
    if(fabs(trk->pdgId())!=PID_) continue; //check is target

    pt_gen[idx] = trk->pt();
    eta_gen[idx] = trk->eta();
    status_gen[idx] = trk->status();
    y_gen[idx] = trk->rapidity();
    
    const auto& mom = findMother(trk);
    idmom[idx] = (mom->pdgId()!=trk->pdgId() ? mom->pdgId() : -77);

    if(decayInGen_)
    {
      for(ushort iDau=0; iDau<NDAU_; iDau++)
      {
        const auto& Dd = trk->daughter(iDau);
        iddau[iDau][candSize_gen] = (Dd ? Dd->pdgId() : 99999);
        chargedau[iDau][candSize_gen] = (Dd ? Dd->charge() : 9);
        ptdau[iDau][candSize_gen] = (Dd ? Dd->pt() : -1.);
        etadau[iDau][candSize_gen] = (Dd ? Dd->eta() : 9.);
        phidau[iDau][candSize_gen] = (Dd ? Dd->phi() : 9.);
        massdau[iDau][candSize_gen] = (Dd ? Dd->mass() : -1.);
      }
    }

    candSize_gen++;
  }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
PATCompositeTreeProducer::beginJob()
{
  TH1D::SetDefaultSumw2();

  if(!doRecoNtuple_ && !doGenNtuple_) { edm::LogError("PATCompositeAnalyzer") << "No output for either RECO or GEN!! Fix config!!" << std::endl; return; }
  if(twoLayerDecay_ && doMuon_) { edm::LogError("PATCompositeAnalyzer") << "Muons cannot be coming from two layer decay!! Fix config!!" << std::endl; return; }

  if(saveHistogram_) initHistogram();
  if(saveTree_) initTree();
}


void
PATCompositeTreeProducer::initHistogram()
{
  for(unsigned int ipt=0;ipt<pTBins_.size()-1;ipt++)
  {
    for(unsigned int iy=0;iy<yBins_.size()-1;iy++)
    {
      hMassVsMVA[iy][ipt] = fs->make<TH2F>(Form("hMassVsMVA_y%d_pt%d",iy,ipt),";mva;mass(GeV)",100,-1.,1.,massHistBins_,massHistPeak_-massHistWidth_,massHistPeak_+massHistWidth_);
      if(saveAllHistogram_)
      {
        hpTVsMVA[iy][ipt] = fs->make<TH2F>(Form("hpTVsMVA_y%d_pt%d",iy,ipt),";mva;pT;",100,-1,1,100,0,10);
        hetaVsMVA[iy][ipt] = fs->make<TH2F>(Form("hetaVsMVA_y%d_pt%d",iy,ipt),";mva;eta;",100,-1.,1.,40,-4,4);
        hyVsMVA[iy][ipt] = fs->make<TH2F>(Form("hyVsMVA_y%d_pt%d",iy,ipt),";mva;y;",100,-1.,1.,40,-4,4);
        hVtxProbVsMVA[iy][ipt] = fs->make<TH2F>(Form("hVtxProbVsMVA_y%d_pt%d",iy,ipt),";mva;VtxProb;",100,-1.,1.,100,0,1);
        h3DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DCosPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;3DCosPointingAngle;",100,-1.,1.,100,-1,1);
        h3DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;3DPointingAngle;",100,-1.,1.,50,-3.14,3.14);
        h2DCosPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DCosPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;2DCosPointingAngle;",100,-1.,1.,100,-1,1);
        h2DPointingAngleVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DPointingAngleVsMVA_y%d_pt%d",iy,ipt),";mva;2DPointingAngle;",100,-1.,1.,50,-3.14,3.14);
        h3DDecayLengthSignificanceVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDecayLengthSignificanceVsMVA_y%d_pt%d",iy,ipt),";mva;3DDecayLengthSignificance;",100,-1.,1.,300,0,30);
        h2DDecayLengthSignificanceVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDecayLengthSignificanceVsMVA_y%d_pt%d",iy,ipt),";mva;2DDecayLengthSignificance;",100,-1.,1.,300,0,30);
        h3DDecayLengthVsMVA[iy][ipt] = fs->make<TH2F>(Form("h3DDecayLengthVsMVA_y%d_pt%d",iy,ipt),";mva;3DDecayLength;",100,-1.,1.,300,0,30);
        h2DDecayLengthVsMVA[iy][ipt] = fs->make<TH2F>(Form("h2DDecayLengthVsMVA_y%d_pt%d",iy,ipt),";mva;2DDecayLength;",100,-1.,1.,300,0,30);
        for(ushort d=1; d<=NDAU_; d++)
        {
          hzDCASignificanceDaugtherVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hzDCASignificanceDaugther%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;zDCASignificanceDaugther%d;",d),100,-1.,1.,100,-10,10);
          hxyDCASignificanceDaugtherVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hxyDCASignificanceDaugther%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;xyDCASignificanceDaugther%d;",d),100,-1.,1.,100,-10,10);
          hNHitDVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hNHitD%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;NHitD%d;",d),100,-1.,1.,100,0,100);
          hpTDVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hpTD%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;pTD%d;",d),100,-1.,1.,100,0,10);
          hpTerrDVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hpTerrD%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;pTerrD%d;",d),100,-1.,1.,50,0,0.5);
          hEtaDVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hEtaD%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;EtaD%d;",d),100,-1.,1.,40,-4,4);
          if(useDeDxData_)
          {
            hdedxHarmonic2DVsMVA[d-1][iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D%dVsMVA_y%d_pt%d",d,iy,ipt),Form(";mva;dedxHarmonic2D%d;",d),100,-1.,1.,100,0,10);
            hdedxHarmonic2DVsP[d-1][iy][ipt] = fs->make<TH2F>(Form("hdedxHarmonic2D%dVsP_y%d_pt%d",d,iy,ipt),Form(";p (GeV);dedxHarmonic2D%d",d),100,0,10,100,0,10);
          }
        }
      }
    }
  }
}


void 
PATCompositeTreeProducer::initTree()
{ 
  PATCompositeNtuple = fs->make< TTree>("VertexCompositeNtuple","VertexCompositeNtuple");

  if(doRecoNtuple_)
  {
    // Event info
    
    PATCompositeNtuple->Branch("RunNb",&runNb,"RunNb/i");
    PATCompositeNtuple->Branch("LSNb",&lsNb,"LSNb/i");
    PATCompositeNtuple->Branch("EventNb",&eventNb,"EventNb/i");
    PATCompositeNtuple->Branch("nPV",&nPV,"nPV/S");
    PATCompositeNtuple->Branch("bestvtxX",&bestvx,"bestvtxX/F");
    PATCompositeNtuple->Branch("bestvtxY",&bestvy,"bestvtxY/F");
    PATCompositeNtuple->Branch("bestvtxZ",&bestvz,"bestvtxZ/F");
    PATCompositeNtuple->Branch("candSize",&candSize,"candSize/s");
    if(isCentrality_) 
    {
      PATCompositeNtuple->Branch("centrality",&centrality,"centrality/S");
      PATCompositeNtuple->Branch("Npixel",&Npixel,"Npixel/I");
      PATCompositeNtuple->Branch("HFsumETPlus",&HFsumETPlus,"HFsumETPlus/F");
      PATCompositeNtuple->Branch("HFsumETMinus",&HFsumETMinus,"HFsumETMinus/F");
      PATCompositeNtuple->Branch("ZDCPlus",&ZDCPlus,"ZDCPlus/F");
      PATCompositeNtuple->Branch("ZDCMinus",&ZDCMinus,"ZDCMinus/F");
      PATCompositeNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
    }
    if(isEventPlane_) {
      PATCompositeNtuple->Branch("ephfpAngle",&ephfpAngle,"ephfpAngle[3]/F");
      PATCompositeNtuple->Branch("ephfmAngle",&ephfmAngle,"ephfmAngle[3]/F");
      PATCompositeNtuple->Branch("ephfpQ",&ephfpQ,"ephfpQ[3]/F");
      PATCompositeNtuple->Branch("ephfmQ",&ephfmQ,"ephfmQ[3]/F");
      PATCompositeNtuple->Branch("ephfpSumW",&ephfpSumW,"ephfpSumW/F");
      PATCompositeNtuple->Branch("ephfmSumW",&ephfmSumW,"ephfmSumW/F");
    }
    PATCompositeNtuple->Branch("trigPrescale",&trigPrescale,Form("trigPrescale[%d]/S",NTRG_));
    PATCompositeNtuple->Branch("trigHLT",&trigHLT,Form("trigHLT[%d]/O",NTRG_));
    PATCompositeNtuple->Branch("evtSel",&evtSel,Form("evtSel[%d]/O",NSEL_));

    // particle info
    PATCompositeNtuple->Branch("pT",&pt,"pT[candSize]/F");
    PATCompositeNtuple->Branch("eta",&eta,"eta[candSize]/F");
    PATCompositeNtuple->Branch("phi",&phi,"phi[candSize]/F");
    PATCompositeNtuple->Branch("mass",&mass,"mass[candSize]/F");
    PATCompositeNtuple->Branch("y",&y,"y[candSize]/F");
    if(useAnyMVA_) PATCompositeNtuple->Branch("mva",&mva,"mva[candSize]/F");

    if(!isSkimMVA_)
    {
      //Composite candidate info RECO
      PATCompositeNtuple->Branch("flavor",&flavor,"flavor[candSize]/F");
      PATCompositeNtuple->Branch("VtxProb",&VtxProb,"VtxProb[candSize]/F");
      PATCompositeNtuple->Branch("3DCosPointingAngle",&agl,"3DCosPointingAngle[candSize]/F");
      PATCompositeNtuple->Branch("3DPointingAngle",&agl_abs,"3DPointingAngle[candSize]/F");
      PATCompositeNtuple->Branch("2DCosPointingAngle",&agl2D,"2DCosPointingAngle[candSize]/F");
      PATCompositeNtuple->Branch("2DPointingAngle",&agl2D_abs,"2DPointingAngle[candSize]/F");
      PATCompositeNtuple->Branch("3DDecayLengthSignificance",&dlos,"3DDecayLengthSignificance[candSize]/F");
      PATCompositeNtuple->Branch("3DDecayLength",&dl,"3DDecayLength[candSize]/F");
      PATCompositeNtuple->Branch("3DDecayLengthError",&dlerror,"3DDecayLengthError[candSize]/F");
      PATCompositeNtuple->Branch("2DDecayLengthSignificance",&dlos2D,"2DDecayLengthSignificance[candSize]/F");
      PATCompositeNtuple->Branch("2DDecayLength",&dl2D,"2DDecayLength[candSize]/F");

      if(doGenMatching_)
      {
        PATCompositeNtuple->Branch("isSwap",&isSwap,"isSwap[candSize]/O");
        PATCompositeNtuple->Branch("idmom_reco",&idmom_reco,"idmom_reco[candSize]/I");
        PATCompositeNtuple->Branch("matchGEN",&matchGEN,"matchGEN[candSize]/O");
      }
 
      if(doGenMatchingTOF_)
      {
        for(ushort iDau=1; iDau<=NDAU_; iDau++)
        {
          PATCompositeNtuple->Branch(Form("PIDD%d",iDau),&pid[iDau-1],Form("PIDD%d[candSize]/I",iDau));
          PATCompositeNtuple->Branch(Form("TOFD%d",iDau),&tof[iDau-1],Form("TOFD%d[candSize]/F",iDau));
        }
      }

      //daugther & grand daugther info
      if(twoLayerDecay_)
      {
        PATCompositeNtuple->Branch("massDaugther1",&grand_mass,"massDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("VtxProbDaugther1",&grand_VtxProb,"VtxProbDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DCosPointingAngleDaugther1",&grand_agl,"3DCosPointingAngleDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DPointingAngleDaugther1",&grand_agl_abs,"3DPointingAngleDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("2DCosPointingAngleDaugther1",&grand_agl2D,"2DCosPointingAngleDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("2DPointingAngleDaugther1",&grand_agl2D_abs,"2DPointingAngleDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DDecayLengthSignificanceDaugther1",&grand_dlos,"3DDecayLengthSignificanceDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DDecayLengthDaugther1",&grand_dl,"3DDecayLengthDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("3DDecayLengthErrorDaugther1",&grand_dlerror,"3DDecayLengthErrorDaugther1[candSize]/F");
        PATCompositeNtuple->Branch("2DDecayLengthSignificanceDaugther1",&grand_dlos2D,"2DDecayLengthSignificanceDaugther1[candSize]/F");
        for(ushort iGDau=1; iGDau<=NGDAU_; iGDau++)
        {
          PATCompositeNtuple->Branch(Form("zDCASignificanceGrandDaugther%d",iGDau),&grand_dzos[iGDau-1],Form("zDCASignificanceGrandDaugther%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("xyDCASignificanceGrandDaugther%d",iGDau),&grand_dxyos[iGDau-1],Form("xyDCASignificanceGrandDaugther%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("NHitGrandD%d",iGDau),&grand_nhit[iGDau-1],Form("NHitGrandD%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("HighPurityGrandDaugther%d",iGDau),&grand_trkquality[iGDau-1],Form("HighPurityGrandDaugther%d[candSize]/O",iGDau));
          PATCompositeNtuple->Branch(Form("pTGrandD%d",iGDau),&grand_pt[iGDau-1],Form("pTGrandD%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("pTerrGrandD%d",iGDau),&grand_ptErr[iGDau-1],Form("pTerrGrandD%d[candSize]/F",iGDau));
          PATCompositeNtuple->Branch(Form("EtaGrandD%d",iGDau),&grand_eta[iGDau-1],Form("EtaGrandD%d[candSize]/F",iGDau));
          if(useDeDxData_)
          {
            PATCompositeNtuple->Branch(Form("dedxPixelHarmonic2GrandD%d",iGDau),&grand_T4dedx[iGDau-1],Form("dedxPixelHarmonic2GrandD%d[candSize]/F",iGDau));
            PATCompositeNtuple->Branch(Form("dedxHarmonic2GrandD%d",iGDau),&grand_H2dedx[iGDau-1],Form("dedxHarmonic2GrandD%d[candSize]/F",iGDau));
          }
        }
      }
      for(ushort iDau=1; iDau<=NDAU_; iDau++)
      {
        PATCompositeNtuple->Branch(Form("zDCASignificanceDaugther%d",iDau),&dzos[iDau-1],Form("zDCASignificanceDaugther%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("xyDCASignificanceDaugther%d",iDau),&dxyos[iDau-1],Form("xyDCASignificanceDaugther%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("NHitD%d",iDau),&nhit[iDau-1],Form("NHitD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("HighPurityDaugther%d",iDau),&trkquality[iDau-1],Form("HighPurityDaugther%d[candSize]/O",iDau));
        PATCompositeNtuple->Branch(Form("pTD%d",iDau),&ptDau[iDau-1],Form("pTD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("pTerrD%d",iDau),&ptErr[iDau-1],Form("pTerrD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("EtaD%d",iDau),&etaDau[iDau-1],Form("EtaD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("PhiD%d",iDau),&phiDau[iDau-1],Form("PhiD%d[candSize]/F",iDau));
        PATCompositeNtuple->Branch(Form("chargeD%d",iDau),&chargeDau[iDau-1],Form("chargeD%d[candSize]/S",iDau));
        if(useDeDxData_)
        {
          PATCompositeNtuple->Branch(Form("dedxPixelHarmonic2D%d",iDau),&T4dedx[iDau-1],Form("dedxPixelHarmonic2D%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("dedxHarmonic2D%d",iDau),&H2dedx[iDau-1],Form("dedxHarmonic2D%d[candSize]/F",iDau));
        }
      }
 
      if(doMuon_)
      {
        for(ushort iDau=1; iDau<=NDAU_; iDau++)
        {
          if(fabs(PID_dau_[iDau-1])!=13) continue;
          PATCompositeNtuple->Branch(Form("OneStMuon%d",iDau),&onestmuon[iDau-1],Form("OneStMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("PFMuon%d",iDau),&pfmuon[iDau-1],Form("PFMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("GlbMuon%d",iDau),&glbmuon[iDau-1],Form("GlbMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("trkMuon%d",iDau),&trkmuon[iDau-1],Form("trkMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("tightMuon%d",iDau),&tightmuon[iDau-1],Form("tightMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("softMuon%d",iDau),&softmuon[iDau-1],Form("softMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("hybridMuon%d",iDau),&hybridmuon[iDau-1],Form("hybridMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("HPMuon%d",iDau),&hpmuon[iDau-1],Form("hybridMuon%d[candSize]/O",iDau));
          PATCompositeNtuple->Branch(Form("trigMuon%d",iDau),trgmuon[iDau-1],Form("trigMuon%d[%d][candSize]/O",iDau,NTRG_));
          PATCompositeNtuple->Branch(Form("nMatchedStationD%d",iDau),&nmatchedst[iDau-1],Form("nMatchedStationD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("nTrackerLayerD%d",iDau),&ntrackerlayer[iDau-1],Form("nTrackerLayerD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("nPixelLayerD%d",iDau),&npixellayer[iDau-1],Form("nPixelLayerD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("nPixelHitD%d",iDau),&npixelhit[iDau-1],Form("nPixelHitD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("nMuonHitD%d",iDau),&nmuonhit[iDau-1],Form("nMuonHitD%d[candSize]/S",iDau));
          PATCompositeNtuple->Branch(Form("GlbTrkChiD%d",iDau),&glbtrkchi[iDau-1],Form("GlbTrkChiD%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("muondXYD%d",iDau),&muonbestdxy[iDau-1],Form("muondXYD%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("muondZD%d",iDau),&muonbestdz[iDau-1],Form("muondZD%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("dXYD%d",iDau),&muondxy[iDau-1],Form("dXYD%d[candSize]/F",iDau));
          PATCompositeNtuple->Branch(Form("dZD%d",iDau),&muondz[iDau-1],Form("dZD%d[candSize]/F",iDau));
          if(doMuonFull_)
          {
            PATCompositeNtuple->Branch(Form("nMatchedChamberD%d",iDau),&nmatchedch[iDau-1],Form("nMatchedChamberD%d[candSize]/S",iDau));
            PATCompositeNtuple->Branch(Form("EnergyDepositionD%d",iDau),&matchedenergy[iDau-1],Form("EnergyDepositionD%d[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("dx%d_seg",iDau),        &dx_seg[iDau-1], Form("dx%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("dy%d_seg",iDau),        &dy_seg[iDau-1], Form("dy%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("dxSig%d_seg",iDau),     &dxSig_seg[iDau-1], Form("dxSig%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("dySig%d_seg",iDau),     &dySig_seg[iDau-1], Form("dySig%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("ddxdz%d_seg",iDau),     &ddxdz_seg[iDau-1], Form("ddxdz%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("ddydz%d_seg",iDau),     &ddydz_seg[iDau-1], Form("ddydz%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("ddxdzSig%d_seg",iDau),  &ddxdzSig_seg[iDau-1], Form("ddxdzSig%d_seg[candSize]/F",iDau));
            PATCompositeNtuple->Branch(Form("ddydzSig%d_seg",iDau),  &ddydzSig_seg[iDau-1], Form("ddydzSig%d_seg[candSize]/F",iDau));
          }
        }
      }
    }
  } // doRecoNtuple_

  if(doGenNtuple_)
  {
    PATCompositeNtuple->Branch("weight_gen",&weight_gen,"weight_gen/F");
    PATCompositeNtuple->Branch("candSize_gen",&candSize_gen,"candSize_gen/I");
    PATCompositeNtuple->Branch("pT_gen",&pt_gen,"pT_gen[candSize_gen]/F");
    PATCompositeNtuple->Branch("eta_gen",&eta_gen,"eta_gen[candSize_gen]/F");
    PATCompositeNtuple->Branch("y_gen",&y_gen,"y_gen[candSize_gen]/F");
    PATCompositeNtuple->Branch("status_gen",&status_gen,"status_gen[candSize_gen]/I");
    PATCompositeNtuple->Branch("MotherID_gen",&idmom,"MotherID_gen[candSize_gen]/I");

    if(decayInGen_)
    {
      for(ushort iDau=0; iDau<NDAU_; iDau++)
      {
        PATCompositeNtuple->Branch(Form("DauID%d_gen",iDau),&iddau[iDau],Form("DauID%d_gen[candSize_gen]/I",iDau));
        PATCompositeNtuple->Branch(Form("chargeD%d_gen",iDau),&chargedau[iDau],Form("chargeD%d_gen[candSize_gen]/I",iDau));
        PATCompositeNtuple->Branch(Form("pTD%d_gen",iDau),&ptdau[iDau],Form("pTD%d_gen[candSize_gen]/F",iDau));
        PATCompositeNtuple->Branch(Form("EtaD%d_gen",iDau),&etadau[iDau],Form("EtaD%d_gen[candSize_gen]/F",iDau));
        PATCompositeNtuple->Branch(Form("PhiD%d_gen",iDau),&phidau[iDau],Form("PhiD%d_gen[candSize_gen]/F",iDau));
        PATCompositeNtuple->Branch(Form("massDaugther%d_gen",iDau),&massdau[iDau],Form("massDaugther%d_gen[candSize_gen]/F",iDau));
      }
    }
  }
}


//--------------------------------------------------------------------------------------------------
void 
PATCompositeTreeProducer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = true;
  EDConsumerBase::Labels triggerResultsLabel;
  EDConsumerBase::labelsForToken(tok_triggerResults_, triggerResultsLabel);
  hltPrescaleProvider_.init(iRun, iSetup, triggerResultsLabel.process, changed);
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
PATCompositeTreeProducer::endJob()
{    
}

reco::GenParticleRef
PATCompositeTreeProducer::findMother(const reco::GenParticleRef& genParRef)
{
  reco::GenParticleRef genMomRef = genParRef;
  int pdg = genParRef->pdgId(); const int pdg_OLD = pdg;
  while(pdg==pdg_OLD && genMomRef->numberOfMothers()>0) {
    genMomRef = genMomRef->motherRef(0);
    pdg = genMomRef->pdgId();
  }
  return ( (pdg_OLD==pdg) ? genParRef : genMomRef );
}


//define this as a plug-in
DEFINE_FWK_MODULE(PATCompositeTreeProducer);