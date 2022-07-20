import FWCore.ParameterSet.Config as cms

trackana = cms.EDAnalyzer('TrackTreeProducer',
  doRecoNtuple = cms.untracked.bool(True),
  VertexCollection = cms.untracked.InputTag("offlinePrimaryVertices"),
  TrackCollection = cms.untracked.InputTag("generalTracks"),
  cutPtMin_ = cms.untracked.double(0.3),
  cutEtaMax_ = cms.untracked.double(2.4),
  cutPtErrOverPt_ = cms.untracked.double(0.1),
  cutDzOverDzSigma_ = cms.untracked.double(5.0),
  cutDxyOverDxySigma_ = cms.untracked.double(5.0),
  doGenNtuple = cms.untracked.bool(False),
  GenParticleCollection = cms.untracked.InputTag("genParticles"),
  saveTree = cms.untracked.bool(True)
                         )

