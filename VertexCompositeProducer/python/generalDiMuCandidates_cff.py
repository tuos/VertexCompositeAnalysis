import FWCore.ParameterSet.Config as cms

from VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cfi import *

HybridSoftIdReco2018 = "(isGlobalMuon && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0)"
SoftIdReco = "(innerTrack.isNonnull && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0 && innerTrack.quality(\"highPurity\"))"
TightIdReco = "(isGlobalMuon && isPFMuon && globalTrack.normalizedChi2 < 10 && globalTrack.hitPattern.numberOfValidMuonHits > 0 && numberOfMatchedStations > 1 && track.hitPattern.trackerLayersWithMeasurement > 5 && track.hitPattern.numberOfValidPixelHits > 0)"

generalMuMuMassMin0Candidates = generalDiMuCandidates.clone(
    mllCutMin = cms.double(0.0),
    muonSelection = cms.string("(p > 2.5 && abs(eta) < 2.4) && ("+HybridSoftIdReco2018+" || "+SoftIdReco+" || "+TightIdReco+")"),
    candidateSelection = cms.string("mass > 0.0")
)
