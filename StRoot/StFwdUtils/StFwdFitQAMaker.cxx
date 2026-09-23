#include "StFwdUtils/StFwdFitQAMaker.h"
#include "StFwdTrackMaker/Common.h"

#include "TMath.h"
#include "TVector3.h"

#include "StFwdTrackMaker/StFwdTrackMaker.h"
#include "StFwdTrackMaker/include/Tracker/FwdTracker.h"
#include "StBFChain/StBFChain.h"

#include "StEvent/StEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StHelixModel.h"
#include "StEvent/StPrimaryTrack.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StPrimaryVertex.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StTrackDetectorInfo.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StTriggerData.h"
#include "StEvent/StFstHitCollection.h"
#include "StEvent/StFstHit.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StChain/StChainOpt.h"

#include "StEventUtilities/StEventHelper.h"

#include "StRoot/StEpdUtil/StEpdGeom.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrackCollection.h"


#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_g2t_event_Table.h"

#include "StarMagField/StarMagField.h"

#include "St_base/StMessMgr.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"


#include "TROOT.h"
#include "TLorentzVector.h"
#include "StEvent/StFwdTrack.h"
#include "StFcsDbMaker/StFcsDb.h"

//________________________________________________________________________
StFwdFitQAMaker::StFwdFitQAMaker() : StMaker("fwdFitQA"), mOutputFilename("StFwdFitQAMaker.root"){};
//________________________________________________________________________
int StFwdFitQAMaker::Init() {
    LOG_DEBUG << "StFwdFitQAMaker::Init" << endm;

    addHist( new TH1F("McQ", ";McQ; counts", 3, -2, 1) );
    addHist( new TH1F("McNumFst", ";NumFST; counts", 10, 0, 10) );
    addHist( new TH1F("McNumFtt", ";NumFTT; counts", 10, 0, 10) );
    addHist( new TH1F("McNumHits", ";NumHits; counts", 10, 0, 10) );
    addHist( new TH1F("McPt", ";McPt; counts", 100, 0, 5) );
    addHist( new TH1F("McEta", ";McEta; counts", 50, 0, 5) );
    addHist( new TH1F("McPhi", ";McPhi; counts", 32, -2*3.1415926, 2*3.1415926) );
    addHist( new TH1F("McStartVtx", ";StartVtx; counts", 10, 0, 10) );
    addHist( new TH1F("McPID", ";GEANT PID; counts", 20, 0, 20) );
    addHist( new TH1F("McPrimaryPt", ";McPrimaryPt; counts", 100, 0, 5) );
    addHist( new TH1F("McPrimaryEta", ";McPrimaryEta; counts", 50, 0, 5) );
    addHist( new TH1F("McPrimaryPhi", ";McPrimaryPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("McPt3FST", ";McPt; counts", 100, 0, 5) );
    addHist( new TH1F("McEta3FST", ";McEta; counts", 50, 0, 5) );
    addHist( new TH1F("McPhi3FST", ";McPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("McPrimaryPt3FST", ";McPrimaryPt (3 FST); counts", 100, 0, 5) );
    addHist( new TH1F("McPrimaryEta3FST", ";McPrimaryEta (3 FST); counts", 50, 0, 5) );
    addHist( new TH1F("McPrimaryPhi3FST", ";McPrimaryPhi (3 FST); counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("RcQ", ";RcQ; counts", 11, -5.5, 5.5) );
    addHist( new TH1F("RcPt", ";RcPt; counts", 100, 0, 5) );
    addHist( new TH1F("RcEta", ";RcEta; counts", 50, 0, 5) );
    addHist( new TH1F("RcPhi", ";RcPhi; counts", 32, -2*3.1415926, 2*3.1415926) );
    addHist( new TH1F("RcPID", "RC;GEANT PID; counts", 20, 0, 20) );

    addHist( new TH1F("RcMatchedMcPt", ";RcMatchedMcPt; counts", 100, 0, 5) );
    addHist( new TH1F("RcMatchedMcEta", ";RcMatchedMcEta; counts", 50, 0, 5) );
    addHist( new TH1F("RcMatchedMcPhi", ";RcMatchedMcPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("RcMatched3FSTMcPtPrimary", ";RcMatched3FSTMcPt; counts", 100, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcEtaPrimary", ";RcMatched3FSTMcEta; counts", 50, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcPhiPrimary", ";RcMatched3FSTMcPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("RcMatched3FSTMcPtGlobal", ";RcMatched3FSTMcPt; counts", 100, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcEtaGlobal", ";RcMatched3FSTMcEta; counts", 50, 0, 5) );
    addHist( new TH1F("RcMatched3FSTMcPhiGlobal", ";RcMatched3FSTMcPhi; counts", 32, -2*3.1415926, 2*3.1415926) );

    addHist( new TH1F("RcQOverPt", ";RcQ/RcPt; counts", 2000, -1, 1) );

    addHist( new TH1F("RcIdTruth", ";RcIdTruth; counts", 102, -2, 100) );
    addHist( new TH1F( "deltaRelPt", ";(pT_{RC} - pT_{MC})/pT_{MC}; count", 400, -2, 2 ) );
    addHist( new TH1F( "deltaRelInvPt", ";(pT_{RC}^{-1} - pT_{MC}^{-1})/pT_{MC}^{-1}; count", 400, -2, 2 ) );
    addHist( new TH1F( "curveResolutionPrim", ";(C_{RC} - C_{MC})/C_{MC}; count", 400, -2, 2 ) );
    addHist( new TH1F( "curveResolutionGlobal", ";(C_{RC} - C_{MC})/C_{MC}; count", 400, -2, 2 ) );
    addHist( new TH2F( "deltaRelInvPtVsPt", ";pT_{MC};(pT_{RC}^{-1} - pT_{MC}^{-1})/pT_{MC}^{-1}; count", 500, 0, 5, 400, -1, 4 ) );

    addHist( new TH2F("RcQMcQ", ";McQ; RcQ", 5, -2.5, 2.5, 5, -2.5, 2.5) );
    addHist( new TH2F("RcPtMcPt", ";McPt; RcPt", 100, 0.0, 1.0, 100, 0.0, 1.0) );
    addHist( new TH2F("RcPtMcPtPrim", "Primary;McPt; RcPt", 100, 0.0, 1.0, 100, 0.0, 1.0) );
    addHist( new TH2F("RcPtMcPtGlobal", "Global;McPt; RcPt", 100, 0.0, 1.0, 100, 0.0, 1.0) );
    addHist( new TH2F( "curveResolutionVsPtPrim", ";Pt_{MC};(C_{RC} - C_{MC})/C_{MC}; ", 100, 0, 5.0, 400, -2, 2 ) );
    addHist( new TH2F( "curveResolutionVsPtGlobal", ";Pt_{MC};(C_{RC} - C_{MC})/C_{MC}; ", 100, 0, 5.0, 400, -2, 2 ) );

    addHist( new TH1F( "RecoTrackStats", "Track Stats (All Reco); Cat;Counts", 10, 0, 10 )  );
    addHist( new TH1F( "GlobalTrackStats", "Track Stats (Global); Cat;Counts", 10, 0, 10 )  );
    addHist( new TH1F( "BeamlineTrackStats", "Track Stats (Beamline); Cat;Counts", 10, 0, 10 )  );
    addHist( new TH1F( "PrimaryTrackStats", "Track Stats (Primary); Cat;Counts", 10, 0, 10 )  );
    addHist( new TH1F( "SecondaryTrackStats", "Track Stats (Secondary); Cat;Counts", 10, 0, 10 )  );
    addHist( new TH1F( "McTrackStats", "MCTrack Stats; Cat;Counts", 10, 0, 10 )  );
    addHist( new TH1F( "TrackSeedStats", "Track Seed Stats; Cat;Counts", 10, 0, 10 )  );

    addHist( new TH1F( "McStartVertex", "Mc Start Vertex; Start Vertex;Counts", 500, -2, 498 )  );

    addHist( new TH1F( "McMult", "Mc Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McPrimaryMult", "Mc Primary Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McSecondaryMult", "Mc Secondary Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdMult", "Mc Fwd Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdPrimaryMult", "Mc Fwd Primary Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdSecondaryMult", "Mc Fwd Secondary Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdFst1Mult", "Mc Fwd Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdFst2Mult", "Mc Fwd Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdFstGE2Mult", "Mc Fwd Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdFst3Mult", "Mc Fwd Mult; N;Counts", 500, 0, 500 )  );

    addHist( new TH1F( "McFwdPrimaryFst1Mult", "Mc Fwd Primary nFST==1 Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdPrimaryFst2Mult", "Mc Fwd Primary nFST==2 Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdPrimaryFstGE2Mult", "Mc Fwd Primary nFST>=2 Mult; N;Counts", 500, 0, 500 )  );
    addHist( new TH1F( "McFwdPrimaryFst3Mult", "Mc Fwd Primary nFST==3 Mult; N;Counts", 500, 0, 500 )  );

    addHist( new TH1F( "McFwdAccMult", "Mc Fwd Primary nFST==3 Pt>0.1 Mult; N;Counts", 500, 0, 500 )  );

    addHist( new TH1F( "GlobalMult", "Global Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "BeamlineMult", "Beamline Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "PrimaryMult", "Primary Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "SecondaryMult", "Secondary Mult; N;Counts", 50, 0, 50 )  );

    addHist( new TH1F( "GlobalAccMult", "Global Acc Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "BeamlineAccMult", "Beamline Acc Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "PrimaryAccMult", "Primary Acc Mult; N;Counts", 50, 0, 50 )  );
    
    addHist( new TH1F( "ConvergedGlobalMult", "Global Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "ConvergedBeamlineMult", "Beamline Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "ConvergedPrimaryMult", "Primary Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "ConvergedSecondaryMult", "Secondary Mult; N;Counts", 50, 0, 50 )  );
    
    addHist( new TH1F( "GoodGlobalMult", "Global Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "GoodBeamlineMult", "Beamline Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "GoodPrimaryMult", "Primary Mult; N;Counts", 50, 0, 50 )  );
    addHist( new TH1F( "GoodSecondaryMult", "Secondary Mult; N;Counts", 50, 0, 50 )  );

    addHist( new TH2F("QidVsPt", "QMatrix; Pt; Qid;", 100, 0, 2.0, 2, -0.5, 1.5 ) );
    addHist( new TH2F("QidVsPtPrim", "QMatrix; Pt; Qid;", 100, 0, 5.0, 2, -0.5, 1.5 ) );
    addHist( new TH2F("QidVsPtGlobal", "QMatrix; Pt; Qid;", 100, 0, 5.0, 2, -0.5, 1.5 ) );
    addHist( new TH2F("QMisidVsPtPrim", "QMatrix; Pt; Qid;", 100, 0, 5.0, 2, -0.5, 1.5 ) );
    addHist( new TH2F("QMisidVsPtGlobal", "QMatrix; Pt; Qid;", 100, 0, 5.0, 2, -0.5, 1.5 ) );

    //Chi2 plots
    makeHistogramSet("Chi2All", 1000, 0, 500, " (All); #Chi^{2}; counts");
    makeHistogramSet("Chi2Good", 1000, 0, 500, " (Good); #Chi^{2}; counts");
    makeHistogramSet("Chi2Bad", 100, 0, 50000, " (Bad); #Chi^{2}; counts");

    // EPD
    addHist( new TH2F("epdId", ";EPD ID (pp); EPD ID (tt)", 200, -0.5, 199.5, 200, -0.5, 199.5), "EPD" );
    addHist( new TH1F("epdEnergy", ";EPD nMIPs; counts", 200, 0, 10), "EPD" );
    addHist( new TH2F("epdEnergyPP", ";EPD ID (pp); EPD nMIPs;",14, -0.5, 13.5, 200, 0, 10), "EPD" );
    addHist( new TH2F("epdEnergyTT", ";EPD ID (tt); EPD nMIPs;",32, -0.5, 31.5, 200, 0, 10), "EPD" );

    addHist( new TH2F("epdTrackdXdY", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100), "EPD" );
    addHist( new TH2F("epdTrackdRdPhi", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()), "EPD" );

    addHist( new TH2F("epdTrackdXdYMixed", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100), "EPD" );
    addHist( new TH2F("epdTrackdRdPhiMixed", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()), "EPD" );

    // EPD vs. Global
    addHist( new TH2F("epdGlobaldXdY", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100), "EPD" );
    addHist( new TH2F("epdGlobaldRdPhi", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()), "EPD" );

    addHist( new TH2F("epdGlobaldXdYMixed", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100), "EPD" );
    addHist( new TH2F("epdGlobaldRdPhiMixed", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()), "EPD" );

    // EPD vs. Beamline
    addHist( new TH2F("epdBeamlinedXdY", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100), "EPD" );
    addHist( new TH2F("epdBeamlinedRdPhi", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()), "EPD" );

    addHist( new TH2F("epdBeamlinedXdYMixed", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100), "EPD" );
    addHist( new TH2F("epdBeamlinedRdPhiMixed", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()), "EPD" );

    // EPD vs. Primary
    addHist( new TH2F("epdPrimarydXdY", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100), "EPD" );
    addHist( new TH2F("epdPrimarydRdPhi", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()), "EPD" );

    addHist( new TH2F("epdPrimarydXdYMixed", ";dX [cm]; dY [cm]", 400, -100, 100, 400, -100, 100), "EPD" );
    addHist( new TH2F("epdPrimarydRdPhiMixed", ";dX [cm]; dY [cm]", 500, 0, 250, 120, -TMath::Pi(), TMath::Pi()), "EPD" );

    // TRACK SEEDS
    makeHistogramSet("SeedFstFtt", 5, 0, 5, 5, 0, 5, " ; nFTT; nFST", "Seed");
    addHist( new TH2F( "SeedFstType", "; type; nFST", 5, 0, 5, 5, 0, 5 ), "Seed");
    addHist( new TH2F( "SeedFttType", "; type; nFTT", 5, 0, 5, 5, 0, 5 ), "Seed");
    addHist( new TH2F( "SeedEpdType", "; type; nEPD", 5, 0, 5, 5, 0, 5 ), "Seed");
    addHist( new TH2F( "SeedQaTruthType", "; type; QA Truth", 5, 0, 5, 101, 0, 101 ), "Seed");


    // EVENT STATS
    addHist( new TH1F( "mNumSeeds", "mNumSeeds", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mAttemptedFits", "mAttemptedFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mGoodFits", "mGoodFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mFailedFits", "mFailedFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mAttemptedGlobalFits", "mAttemptedGlobalFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mFailedGlobalFits", "mFailedGlobalFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mGoodGlobalFits", "mGoodGlobalFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mFailedGlobalRefits", "mFailedGlobalRefits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mGoodGlobalRefits", "mGoodGlobalRefits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mNumFwdVertices", "mNumFwdVertices", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mAttemptedPrimaryFits", "mAttemptedPrimaryFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mGoodPrimaryFits", "mGoodPrimaryFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mFailedPrimaryFits", "mFailedPrimaryFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mFailedPrimaryRefits", "mFailedPrimaryRefits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mGoodPrimaryRefits", "mGoodPrimaryRefits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mAttemptedBeamlineFits", "mAttemptedBeamlineFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mGoodBeamlineFits", "mGoodBeamlineFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mFailedBeamlineFits", "mFailedBeamlineFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mFailedBeamlineRefits", "mFailedBeamlineRefits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mGoodBeamlineRefits", "mGoodBeamlineRefits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mAttemptedSecondaryFits", "mAttemptedSecondaryFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mGoodSecondaryFits", "mGoodSecondaryFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mFailedSecondaryFits", "mFailedSecondaryFits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mFailedSecondaryRefits", "mFailedSecondaryRefits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mGoodSecondaryRefits", "mGoodSecondaryRefits", 100, 0, 100 ), "EventStats" );

    addHist( new TH1F( "numGlobalFoundHits", "numGlobalFoundHits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "numBeamlineFoundHits", "numBeamlineFoundHits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "numPrimaryFoundHits", "numPrimaryFoundHits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "numSecondaryFoundHits", "numSecondaryFoundHits", 100, 0, 100 ), "EventStats" );
    addHist( new TH1F( "mStep1Duration", "mStep1Duration", 500, 0, 500 ), "EventStats" );
    addHist( new TH1F( "mSeedFindingDuration", "mSeedFindingDuration", 500, 0, 500 ), "EventStats" );
    addHist( new TH1F( "mStep2Duration", "mStep2Duration", 500, 0, 500 ), "EventStats" );
    addHist( new TH1F( "mStep3Duration", "mStep3Duration", 500, 0, 500 ), "EventStats" );
    addHist( new TH1F( "mStep4Duration", "mStep4Duration", 500, 0, 500 ), "EventStats" );
    addHist( new TH1F( "mGlobalFitDuration", "mGlobalFitDuration", 500, 0, 500 ), "EventStats" );
    addHist( new TH1F( "mBeamlineFitDuration", "mBeamlineFitDuration", 500, 0, 500 ), "EventStats" );
    addHist( new TH1F( "mPrimaryFitDuration", "mPrimaryFitDuration", 500, 0, 500 ), "EventStats" );
    addHist( new TH1F( "mSecondaryFitDuration", "mSecondaryFitDuration", 500, 0, 500 ), "EventStats" );

    return kStOK;
}
//________________________________________________________________________
int StFwdFitQAMaker::Finish() {

    auto prevDir = gDirectory;

    // output file name
    string name = "StFwdFitQAMaker.root";
    TFile *fOutput = new TFile(mOutputFilename.Data(), "RECREATE");
    fOutput->cd();


    // accumulated analysis
    addHist( (TH1*)(getHist( "RcMatchedMcEta" )->Clone( "EffEta" )) );
    getHist("EffEta")->Divide( (TH1*)getHist( "McEta" ) );

    addHist( (TH1*)(getHist( "RcMatchedMcPt" )->Clone( "EffPt" )) );
    getHist("EffPt")->Divide( (TH1*)getHist( "McPt" ) );

    addHist( (TH1*)(getHist( "RcMatchedMcPhi" )->Clone( "EffPhi" )) );
    getHist("EffPhi")->Divide( (TH1*)getHist( "McPhi" ) );


        // 3 FST hits on MC Track
    addHist( (TH1*)(getHist( "RcMatched3FSTMcEtaPrimary" )->Clone( "Eff3FSTEtaPrimary" )) );
    getHist("Eff3FSTEtaPrimary")->Divide( (TH1*)getHist( "McEta3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPtPrimary" )->Clone( "Eff3FSTPtPrimary" )) );
    getHist("Eff3FSTPtPrimary")->Divide( (TH1*)getHist( "McPt3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPhiPrimary" )->Clone( "Eff3FSTPhiPrimary" )) );
    getHist("Eff3FSTPhiPrimary")->Divide( (TH1*)getHist( "McPhi3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcEtaGlobal" )->Clone( "Eff3FSTEtaGlobal" )) );
    getHist("Eff3FSTEtaGlobal")->Divide( (TH1*)getHist( "McEta3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPtGlobal" )->Clone( "Eff3FSTPtGlobal" )) );
    getHist("Eff3FSTPtGlobal")->Divide( (TH1*)getHist( "McPt3FST" ) );

    addHist( (TH1*)(getHist( "RcMatched3FSTMcPhiGlobal" )->Clone( "Eff3FSTPhiGlobal" )) );
    getHist("Eff3FSTPhiGlobal")->Divide( (TH1*)getHist( "McPhi3FST" ) );

    getHist( "McTrackStats")->GetXaxis()->SetBinLabel( 1, "ALLMC" );
    getHist( "McTrackStats")->GetXaxis()->SetBinLabel( 1, "PRIMMC" );
    getHist( "McTrackStats")->GetXaxis()->SetBinLabel( 2, "MCFWD" );
    getHist( "McTrackStats")->GetXaxis()->SetBinLabel( 3, "MCAcc" );
    
    getHist( "GlobalTrackStats")->GetXaxis()->SetBinLabel( 1, "All" );
    getHist( "GlobalTrackStats")->GetXaxis()->SetBinLabel( 2, "Converge" );    
    getHist( "GlobalTrackStats")->GetXaxis()->SetBinLabel( 3, "Good" );
    getHist( "GlobalTrackStats")->GetXaxis()->SetBinLabel( 4, "q+" );
    getHist( "GlobalTrackStats")->GetXaxis()->SetBinLabel( 5, "q0" );
    getHist( "GlobalTrackStats")->GetXaxis()->SetBinLabel( 6, "q-" );

    getHist( "BeamlineTrackStats")->GetXaxis()->SetBinLabel( 1, "All" );
    getHist( "BeamlineTrackStats")->GetXaxis()->SetBinLabel( 2, "Converge" );    
    getHist( "BeamlineTrackStats")->GetXaxis()->SetBinLabel( 3, "Good" );
    getHist( "BeamlineTrackStats")->GetXaxis()->SetBinLabel( 4, "q+" );
    getHist( "BeamlineTrackStats")->GetXaxis()->SetBinLabel( 5, "q0" );
    getHist( "BeamlineTrackStats")->GetXaxis()->SetBinLabel( 6, "q-" );

    getHist( "PrimaryTrackStats")->GetXaxis()->SetBinLabel( 1, "All" );
    getHist( "PrimaryTrackStats")->GetXaxis()->SetBinLabel( 2, "Converge" );    
    getHist( "PrimaryTrackStats")->GetXaxis()->SetBinLabel( 3, "Good" );
    getHist( "PrimaryTrackStats")->GetXaxis()->SetBinLabel( 4, "q+" );
    getHist( "PrimaryTrackStats")->GetXaxis()->SetBinLabel( 5, "q0" );
    getHist( "PrimaryTrackStats")->GetXaxis()->SetBinLabel( 6, "q-" );

    getHist( "SecondaryTrackStats")->GetXaxis()->SetBinLabel( 1, "All" );
    getHist( "SecondaryTrackStats")->GetXaxis()->SetBinLabel( 2, "Converge" );    
    getHist( "SecondaryTrackStats")->GetXaxis()->SetBinLabel( 3, "Good" );
    getHist( "SecondaryTrackStats")->GetXaxis()->SetBinLabel( 4, "q+" );
    getHist( "SecondaryTrackStats")->GetXaxis()->SetBinLabel( 5, "q0" );
    getHist( "SecondaryTrackStats")->GetXaxis()->SetBinLabel( 6, "q-" );


    for (auto nh : mHists) {
        // check if name containes "Global" or "Beamline" or "Primary" or "Secondary"
        if (nh.first.Contains("Global") ) {
            nh.second->SetLineColor(kRed);
        } else if (nh.first.Contains("Beamline") ) {
            nh.second->SetLineColor(kBlue);
        } else if (nh.first.Contains("Primary") ) {
            nh.second->SetLineColor(kMagenta);
        } else if (nh.first.Contains("Secondary") ) {
            nh.second->SetLineColor(kGreen);
        } else if (nh.first.Contains("Mc") ) {
            nh.second->SetLineColor(kBlack);
        }

        if ( mHistsDirectories.count( nh.first ) > 0 && mHistsDirectories[nh.first] != "" ) {
            fOutput->mkdir( mHistsDirectories[nh.first].Data() );
            fOutput->cd( mHistsDirectories[nh.first].Data() );
        } else {
            fOutput->cd();
            // gDirectory->cd();
        }
        nh.second->SetDirectory(gDirectory);
        nh.second->Write();
    }

    // restore previous directory
    gDirectory = prevDir;

    LOG_INFO << "Writing StFwdFitQAMaker output" << endm;

    return kStOk;
}
//________________________________________________________________________
int StFwdFitQAMaker::Make() {
    LOG_DEBUG << "StFwdFitQAMaker::Make" << endm;
    long long itStart = FwdTrackerUtils::nowNanoSecond();
    if (!mAnalyzeMuDst)
        ProcessFwdTracks();
    else
        ProcessFwdMuTracks();
    
    ProcessData();
    LOG_DEBUG << "Processing Fwd Track Fit QA took: " << (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e6 << " ms" << endm;
    return kStOK;
} // Make
//________________________________________________________________________
void StFwdFitQAMaker::Clear(const Option_t *opts) { LOG_DEBUG << "StFwdFitQAMaker::CLEAR" << endm; }
//________________________________________________________________________
void StFwdFitQAMaker::ProcessData(  ){
    LOG_INFO << "StFwdFitQAMaker::ProcessData" << endm;
    FillEventStats();
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent)
        return;
    

    auto mFcsDb = dynamic_cast<StFcsDb*>(GetDataSet("fcsDb"));
    if (!mFcsDb)
        return;

    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    auto fcsc = stEvent->fcsCollection();
    if (!ftc || !fcsc)
        return;
    
    StEpdGeom epdgeo;

    std::vector<TVector3> mFcsPreHits;
    // LOAD PRESHOWER HITS (EPD)
    for ( int det = 4; det < 6; det ++ ) {

        StSPtrVecFcsHit& hits = stEvent->fcsCollection()->hits(det);
        int nh=stEvent->fcsCollection()->numberOfHits(det);
        LOG_INFO << "Det " << det << " has " << nh << " hits" << endm;
        for ( int i = 0; i < nh; i++ ){
            StFcsHit* hit=hits[i];

            if(det==kFcsPresNorthDetId || det==kFcsPresSouthDetId){ //EPD
                double zepd=375.0;
                int pp,tt,n;
                double x[5],y[5];

                mFcsDb->getEPDfromId(det,hit->id(),pp,tt);
                getHist2( "epdId" )->Fill( pp,tt );
                getHist ( "epdEnergy" )->Fill( hit->energy() );
                getHist2( "epdEnergyPP" )->Fill( pp, hit->energy() );
                getHist2( "epdEnergyTT" )->Fill( tt, hit->energy() );
                if ( hit->energy() < 0.2 ) continue;

                // Get STAR position of EPD tile corners
                epdgeo.GetCorners(100*pp+tt,&n,x,y);
                double x0 = (x[0] + x[1] + x[2] + x[3]) / 4.0;
                double y0 = (y[0] + y[1] + y[2] + y[3]) / 4.0;
                mFcsPreHits.push_back( TVector3( x0, y0, zepd ) );
            } // if det
        } // for i
    } // for det

    size_t nGlobal = 0;
    size_t nBeamline = 0;
    size_t nPrimary = 0;
    size_t nSecondary = 0;

    size_t nBeamlineAcc = 0;

    size_t nConvergedGlobal = 0;
    size_t nConvergedBeamline = 0;
    size_t nConvergedPrimary = 0;
    size_t nConvergedSecondary = 0;

    size_t nGoodGlobal = 0;
    size_t nGoodBeamline = 0;
    size_t nGoodPrimary = 0;
    size_t nGoodSecondary = 0;

    float minChi2 = 1e-5;
    float maxChi2 = 1e11;
    LOG_INFO << "The StFwdFitQAMaker has found number of Fwd Tracks: " << ftc->numberOfTracks() << endm;
    for ( auto fwdTrack : ftc->tracks() ){

        FillTrackSeedHistograms( fwdTrack );

        if (fwdTrack->isGlobalTrack() ){
            nGlobal++;
            getHist( "GlobalTrackStats" )->Fill( 0 );
            getHist( "GlobalChi2All" )->Fill( fwdTrack->chi2() );
            if ( fwdTrack->didFitConverge() ) {
                nConvergedGlobal++;
                getHist( "GlobalTrackStats" )->Fill( 1 );
                if ( fwdTrack->chi2() > minChi2 && fwdTrack->chi2() < maxChi2 ) {
                    getHist( "GlobalChi2Good" )->Fill( fwdTrack->chi2() );
                    nGoodGlobal++;
                    getHist( "GlobalTrackStats" )->Fill( 2 );
                    if ( fwdTrack->charge() == 1 ) 
                        getHist( "GlobalTrackStats" )->Fill( 3 );
                    else if ( fwdTrack->charge() == 0 ) 
                        getHist( "GlobalTrackStats" )->Fill( 4 );
                    else if ( fwdTrack->charge() == -1 )
                        getHist( "GlobalTrackStats" )->Fill( 5 );
                    else {
                        LOG_INFO << "Unknown charge on global: " << fwdTrack->charge() << endm;
                    }
                } else {
                    getHist( "GlobalChi2Bad" )->Fill( fwdTrack->chi2() );
                }
            }
        }
        if (fwdTrack->isBeamLineConstrainedTrack() ){
            nBeamline++;
            getHist( "BeamlineTrackStats" )->Fill( 0 );
            getHist( "BeamlineChi2All" )->Fill( fwdTrack->chi2() );
            if ( fwdTrack->didFitConverge() ) {
                nConvergedBeamline++;
                if ( fwdTrack->momentum().perp() > 0.1 ) nBeamlineAcc++;
                getHist( "BeamlineTrackStats" )->Fill( 1 );
                if ( fwdTrack->chi2() > minChi2 && fwdTrack->chi2() < maxChi2 ) {
                    getHist( "BeamlineChi2Good" )->Fill( fwdTrack->chi2() );
                    nGoodBeamline++;
                    getHist( "BeamlineTrackStats" )->Fill( 2 );
                    if ( fwdTrack->charge() == 1 ) 
                        getHist( "BeamlineTrackStats" )->Fill( 3 );
                    else if ( fwdTrack->charge() == 0 ) 
                        getHist( "BeamlineTrackStats" )->Fill( 4 );
                    else if ( fwdTrack->charge() == -1 )
                        getHist( "BeamlineTrackStats" )->Fill( 5 );
                    else {
                        LOG_INFO << "Unknown charge on beamline: " << fwdTrack->charge() << endm;
                    }
                } else {
                    getHist( "BeamlineChi2Bad" )->Fill( fwdTrack->chi2() );
                }
            }
        }
        if (fwdTrack->isPrimaryTrack() ){
            nPrimary++;
            getHist( "PrimaryTrackStats" )->Fill( 0 );
            getHist( "PrimaryChi2All" )->Fill( fwdTrack->chi2() );
            if ( fwdTrack->didFitConverge() ) {
                nConvergedPrimary++;
                getHist( "PrimaryTrackStats" )->Fill( 1 );
                if ( fwdTrack->chi2() > minChi2 && fwdTrack->chi2() < maxChi2 ) {
                    nGoodPrimary++;
                    getHist( "PrimaryChi2Good" )->Fill( fwdTrack->chi2() );
                    getHist( "PrimaryTrackStats" )->Fill( 2 );
                    if ( fwdTrack->charge() == 1 ) 
                        getHist( "PrimaryTrackStats" )->Fill( 3 );
                    else if ( fwdTrack->charge() == 0 ) 
                        getHist( "PrimaryTrackStats" )->Fill( 4 );
                    else if ( fwdTrack->charge() == -1 )
                        getHist( "PrimaryTrackStats" )->Fill( 5 );
                    else {
                        LOG_INFO << "Unknown charge on primary: " << fwdTrack->charge() << endm;
                    }
                } else {
                    getHist( "PrimaryChi2Bad" )->Fill( fwdTrack->chi2() );
                }
            }
        }
        if (fwdTrack->isFwdVertexConstrainedTrack() ){
            nSecondary++;
            getHist( "SecondaryTrackStats" )->Fill( 0 );
            if ( fwdTrack->didFitConverge() ) {
                nConvergedSecondary++;
                getHist( "SecondaryTrackStats" )->Fill( 1 );
                getHist( "SecondaryChi2All" )->Fill( fwdTrack->chi2() );
                if ( fwdTrack->chi2() > minChi2 && fwdTrack->chi2() < maxChi2 ) {
                    nGoodSecondary++;
                    getHist( "SecondaryChi2Good" )->Fill( fwdTrack->chi2() );
                    getHist( "SecondaryTrackStats" )->Fill( 2 );
                    if ( fwdTrack->charge() == 1 ) 
                        getHist( "SecondaryTrackStats" )->Fill( 3 );
                    else if ( fwdTrack->charge() == 0 ) 
                        getHist( "SecondaryTrackStats" )->Fill( 4 );
                    else if ( fwdTrack->charge() == -1 )
                        getHist( "SecondaryTrackStats" )->Fill( 5 );
                    else {
                        LOG_INFO << "Unknown charge on secondary: " << fwdTrack->charge() << endm;
                    }
                } else {
                    getHist( "SecondaryChi2Bad" )->Fill( fwdTrack->chi2() );
                }
            }
        }

        // get EPD projection
        auto proj = fwdTrack->getProjectionFor(kFcsPresId);
        
        if ( fwdTrack->chi2() < 0.01 ) continue;
        // now loop on EPD hits
        for ( auto epdHit : mFcsPreHits){
            TVector3 track(proj.mXYZ.x(),proj.mXYZ.y(),proj.mXYZ.z());
            double dx = track.X() - epdHit.X();
            double dy = track.Y() - epdHit.Y();
            double dR = TMath::Sqrt(dx*dx + dy*dy);
            double dPhi = track.DeltaPhi(epdHit);

            getHist2( "epdTrackdXdY" )->Fill( dx, dy );
            getHist2( "epdTrackdRdPhi" )->Fill( dR, dPhi );

            if ( fwdTrack->isGlobalTrack() ){
                getHist2( "epdGlobaldXdY" )->Fill( dx, dy );
                getHist2( "epdGlobaldRdPhi" )->Fill( dR, dPhi );
            }
            if ( fwdTrack->isBeamLineConstrainedTrack() ){
                getHist2( "epdBeamlinedXdY" )->Fill( dx, dy );
                getHist2( "epdBeamlinedRdPhi" )->Fill( dR, dPhi );
            }
            if ( fwdTrack->isPrimaryTrack() ){
                getHist2( "epdPrimarydXdY" )->Fill( dx, dy );
                getHist2( "epdPrimarydRdPhi" )->Fill( dR, dPhi );
            }
        }

        // Mixed event
        for ( auto epdHit : mFcsPreHitsLastEvent){
            TVector3 track(proj.mXYZ.x(),proj.mXYZ.y(),proj.mXYZ.z());
            double dx = track.X() - epdHit.X();
            double dy = track.Y() - epdHit.Y();
            double dR = TMath::Sqrt(dx*dx + dy*dy);
            double dPhi = track.DeltaPhi(epdHit);

            getHist2( "epdTrackdXdYMixed" )->Fill( dx, dy );
            getHist2( "epdTrackdRdPhiMixed" )->Fill( dR, dPhi );
        }
        
    }

    mFcsPreHitsLastEvent.clear();
    for ( auto epdHit : mFcsPreHits){
        mFcsPreHitsLastEvent.push_back(epdHit);
    }

    // Reco track stats
    getHist( "RecoTrackStats" )->Fill( 1, nGlobal );
    getHist( "RecoTrackStats" )->Fill( 2, nBeamline );
    getHist( "RecoTrackStats" )->Fill( 3, nPrimary );
    getHist( "RecoTrackStats" )->Fill( 4, nSecondary );
    
    // Multiplicity plots
    getHist( "GlobalMult" )->Fill( nGlobal );
    getHist( "BeamlineMult" )->Fill( nBeamline );
    getHist( "PrimaryMult" )->Fill( nPrimary );
    getHist( "SecondaryMult" )->Fill( nSecondary );

    getHist( "BeamlineAccMult" )->Fill( nBeamlineAcc );

    // Converged Multiplicity plots
    getHist( "ConvergedGlobalMult" )->Fill( nConvergedGlobal );
    getHist( "ConvergedBeamlineMult" )->Fill( nConvergedBeamline );
    getHist( "ConvergedPrimaryMult" )->Fill( nConvergedPrimary );
    getHist( "ConvergedSecondaryMult" )->Fill( nConvergedSecondary );
    
    // Good Multiplicity plots
    getHist( "GoodGlobalMult" )->Fill( nGoodGlobal );
    getHist( "GoodBeamlineMult" )->Fill( nGoodBeamline );
    getHist( "GoodPrimaryMult" )->Fill( nGoodPrimary );
    getHist( "GoodSecondaryMult" )->Fill( nGoodSecondary );

}
//________________________________________________________________________
void StFwdFitQAMaker::ProcessFwdTracks(  ){
    // This is an example of how to process fwd track collection
    LOG_DEBUG << "StFwdFitQAMaker::ProcessFwdTracks" << endm;
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent)
        return;
    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    if (!ftc)
        return;

    // Get geant tracks
    St_g2t_track *g2t_track = (St_g2t_track *)GetDataSet("geant/g2t_track");

    if (!g2t_track)
        return;

    mcTracks.clear();
    LOG_DEBUG << g2t_track->GetNRows() << " mc tracks in geant/g2t_track " << endm;
    for (int irow = 0; irow < g2t_track->GetNRows(); irow++) {
        g2t_track_st *track = (g2t_track_st *)g2t_track->At(irow);

        if (0 == track)
            continue;

        McFwdTrack mct;

        mct.id = track->id;
        mct.p.SetXYZ( track->p[0], track->p[1], track->p[2] );
        mct.q = track->charge;
        mct.g2track = track;
        mct.pid = (int)track->ge_pid;

        mcTracks.push_back( mct );

        getHist( "McPID" )->Fill( mct.pid );
        getHist( "McEta" )->Fill( mct.p.Eta() );
        getHist( "McPt" )->Fill( mct.p.Pt() );
        getHist( "McPhi" )->Fill( mct.p.Phi() );
        getHist( "McQ" )->Fill( mct.p.Pt() );

        if (mct.g2track->start_vertex_p == 1){
            getHist( "McPrimaryEta" )->Fill( mct.p.Eta() );
            getHist( "McPrimaryPt" )->Fill( mct.p.Pt() );
            getHist( "McPrimaryPhi" )->Fill( mct.p.Phi() );
        }
        
    } // get MC tracks

    /***********************/
    St_g2t_fts_hit *g2t_fsi_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_fsi_hit");
    if ( !g2t_fsi_hits ){
        LOG_DEBUG << "No g2t_fts_hits, cannot load FST hits from GEANT" << endm;
        return;
    }
    int nfsi = g2t_fsi_hits->GetNRows();

    for (int i = 0; i < nfsi; i++) {
        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_fsi_hits->At(i);
        if (0 == git) continue; // geant hit
        int track_id = git->track_p - 1;
        if ( track_id >= 0 && track_id < mcTracks.size() ){
            mcTracks[track_id].numFST++;
        }
    } // Get FST hits on MC tracks

    St_g2t_fts_hit *g2t_stg_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_stg_hit");
    if (!g2t_stg_hits){
        LOG_WARN << "geant/g2t_stg_hit is empty" << endm;
        return;
    }
    int nstg = g2t_stg_hits->GetNRows();
    for (int i = 0; i < nstg; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_stg_hits->At(i);
        if (0 == git) continue; // geant hit
        int track_id = git->track_p - 1;
        int volume_id = git->volume_id;

        // only use the hits on the front modules
        if ( volume_id % 2 ==0 )
            continue;
        if ( track_id >= 0 && track_id < mcTracks.size() ){
            mcTracks[track_id].numFTT++;
        }
    }
    
    getHist( "McMult" ) ->Fill( mcTracks.size() );

    auto isFwd    = [](McFwdTrack const& t){ return t.numFST >= 1 || t.numFTT >= 1; };
    auto isFst1   = [](McFwdTrack const& t){ return t.numFST == 1; };
    auto isFst2   = [](McFwdTrack const& t){ return t.numFST == 2; };
    auto isFst3   = [](McFwdTrack const& t){ return t.numFST == 3; };
    auto isSecondary   = [](McFwdTrack const& t){ return t.g2track->start_vertex_p != 1; };
    auto isPrimary     = [](McFwdTrack const& t){ return t.g2track->start_vertex_p == 1; };
    auto isFwdPrimary  = [](McFwdTrack const& t){ return (t.numFST >= 1 || t.numFTT >= 1) && t.g2track->start_vertex_p == 1; };
    auto isFwdSecondary= [](McFwdTrack const& t){ return (t.numFST >= 1 || t.numFTT >= 1) && t.g2track->start_vertex_p != 1; };

    auto isFwdPrimaryFST1  = [](McFwdTrack const& t){ return (t.numFST == 1 ) && t.g2track->start_vertex_p == 1; };
    auto isFwdPrimaryFST2  = [](McFwdTrack const& t){ return (t.numFST == 2 ) && t.g2track->start_vertex_p == 1; };
    auto isFwdPrimaryFSTGE2  = [](McFwdTrack const& t){ return (t.numFST >= 2 ) && t.g2track->start_vertex_p == 1; };
    auto isFwdPrimaryFST3  = [](McFwdTrack const& t){ return (t.numFST == 3 ) && t.g2track->start_vertex_p == 1; };
    auto isFwdAcc  = [](McFwdTrack const& t){ return (t.numFST == 3 ) && t.g2track->start_vertex_p == 1 && t.p.Perp() > 0.1; };
    
    getHist( "McFwdMult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFwd) );
    getHist( "McFwdFst1Mult" ) ->Fill(std::count_if(mcTracks.begin(), mcTracks.end(), isFst1) );
    getHist( "McFwdFst2Mult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFst2) );
    getHist( "McFwdFst3Mult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFst3) );
    getHist( "McPrimaryMult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isPrimary ) );
    getHist( "McSecondaryMult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isSecondary ) );
    getHist( "McFwdPrimaryMult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFwdPrimary ) );
    getHist( "McFwdSecondaryMult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFwdSecondary ) );

    getHist( "McFwdPrimaryFst1Mult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFwdPrimaryFST1 ) );
    getHist( "McFwdPrimaryFst2Mult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFwdPrimaryFST2 ) );
    getHist( "McFwdPrimaryFstGE2Mult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFwdPrimaryFSTGE2 ) );
    getHist( "McFwdPrimaryFst3Mult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFwdPrimaryFST3 ) );
    getHist( "McFwdAccMult" ) ->Fill( std::count_if(mcTracks.begin(), mcTracks.end(), isFwdAcc ) );

    for ( auto mct : mcTracks ){
        LOG_INFO << "Track ID = " << mct.id << " numFST = " << mct.numFST << " numFTT = " << mct.numFTT << endm;
        LOG_INFO 
            << "\t n_fts_hit=" << mct.g2track->n_fts_hit 
            << ", n_stg_hit=" << mct.g2track->n_stg_hit 
            << ", n_pre_hit=" << mct.g2track->hit_stg_p 
            << ", n_wca_hit=" << mct.g2track->n_wca_hit
            << ", n_hca_hit=" << mct.g2track->hit_hca_p 
            << endm;
        LOG_INFO << "\t start_vertex_p=" << mct.g2track->start_vertex_p << endm;

        getHist( "McTrackStats" )->Fill( 1 ); // all
        
        getHist( "McNumFst" )->Fill( mct.numFST );
        getHist( "McNumFtt" )->Fill( mct.numFTT );
        getHist( "McNumHits" )->Fill( mct.numFTT + mct.numFST );

        if ( mct.numFST >= 1 || mct.numFTT >= 1 ){
            getHist( "McTrackStats" )->Fill( 1 );
        }
        if ( mct.numFST >= 3 ){
            getHist( "McEta3FST" )->Fill( mct.p.Eta() );
            getHist( "McPt3FST" )->Fill( mct.p.Pt() );
            getHist( "McPhi3FST" )->Fill( mct.p.Phi() );

            if ( mct.g2track->start_vertex_p == 1 ){
                getHist( "McPrimaryEta3FST" )->Fill( mct.p.Eta() );
                getHist( "McPrimaryPt3FST" )->Fill( mct.p.Pt() );
                getHist( "McPrimaryPhi3FST" )->Fill( mct.p.Phi() );
            }

        }
        if ( mct.numFST >=3 || mct.numFTT >= 4){
            getHist( "McTrackStats" )->Fill( 2 );
        }
    }

    /**************************************************************************
    */
    for ( auto fwdTrack : ftc->tracks() ){
        
        int idx = fwdTrack->idTruth() - 1;
        
        getHist( "RcIdTruth" )->Fill( fwdTrack->idTruth() );
        getHist( "RcEta" )->Fill( fwdTrack->momentum().pseudoRapidity() );
        getHist( "RcPt" )->Fill( fwdTrack->momentum().perp() );
        getHist( "RcPhi" )->Fill( fwdTrack->momentum().phi() );

        LOG_INFO << "idx = " << idx << endm;
        if ( idx < 0 || idx >= mcTracks.size() ) continue;
        McFwdTrack mct = mcTracks[idx];
        

        getHist( "RcMatchedMcEta" )->Fill( mct.p.Eta() );
        getHist( "RcMatchedMcPt" )->Fill( mct.p.Pt() );
        getHist( "RcMatchedMcPhi" )->Fill( mct.p.Phi() );

        getHist( "RcPtMcPt" )->Fill( mct.p.Pt(), fwdTrack->momentum().perp() );

        if ( fwdTrack->vertexIndex() < 100 ){
            getHist( "RcPtMcPtPrim" )->Fill( mct.p.Pt(), fwdTrack->momentum().perp() );
        } else {
            getHist( "RcPtMcPtGlobal" )->Fill( mct.p.Pt(), fwdTrack->momentum().perp() );
        }

        if ( mct.numFST < 1 || mct.numFTT < 1 ) continue;
        // if ( mct.numFTT < 4 ) continue;


        if ( fwdTrack->trackType() == StFwdTrack::kPrimaryVertexConstrained ){
            getHist( "RcMatched3FSTMcEtaPrimary" ) ->Fill( mct.p.Eta() );
            getHist( "RcMatched3FSTMcPtPrimary" )  ->Fill( mct.p.Pt() );
            getHist( "RcMatched3FSTMcPhiPrimary" ) ->Fill( mct.p.Phi() );
        } else {
            getHist( "RcMatched3FSTMcEtaGlobal" ) ->Fill( mct.p.Eta() );
            getHist( "RcMatched3FSTMcPtGlobal" )  ->Fill( mct.p.Pt() );
            getHist( "RcMatched3FSTMcPhiGlobal" ) ->Fill( mct.p.Phi() );
        }
        getHist( "RcPID" )              ->Fill( mct.pid );
        ((TH2*)getHist( "RcQMcQ" ))     ->Fill( mct.q, fwdTrack->charge() );
        double qop = fwdTrack->charge() / fwdTrack->momentum().perp();
        getHist( "RcQOverPt" )           ->Fill( qop );
        getHist( "RcQ" )                ->Fill( fwdTrack->charge() );

        //double dEta  = fwdTrack->momentum().pseudoRapidity() - mct.p.Eta();
        // double dPhi  = fwdTrack->momentum().DeltaPhi( mct.p );
        // double dR    = fwdTrack->momentum().DeltaR( mct.p );
        //double dPt   = fwdTrack->momentum().perp() - mct.p.Pt();
        double dRPt  = (fwdTrack->momentum().perp() - mct.p.Pt()) / mct.p.Pt();
        double iPt   = 1.0 / fwdTrack->momentum().perp();
        double curve   = 1.0 / fwdTrack->momentum().perp();
        double curveMc = 1.0 / mct.p.Pt();
        double iMcPt = 1.0 / mct.p.Pt();
        double dRIPt = (iPt - iMcPt) / iMcPt;
        double dCurve = (curve - curveMc) / curveMc;

        getHist( "deltaRelPt" ) ->Fill( dRPt );
        getHist( "deltaRelInvPt" )->Fill( dRIPt );
        ((TH2*)getHist( "deltaRelInvPtVsPt" ))->Fill( mct.p.Pt(), dRIPt );
        if ( mct.g2track )
            getHist( "McStartVtx" )->Fill( mct.g2track->start_vertex_p );
        
        if ( fwdTrack->trackType() == StFwdTrack::kPrimaryVertexConstrained ){
            getHist( "curveResolutionPrim" )->Fill( dCurve );
            getHist( "curveResolutionVsPtPrim" )->Fill( mct.p.Pt(), dCurve );
            getHist("QidVsPtPrim")->Fill( mct.p.Pt(), mct.q == fwdTrack->charge() ? 1 : 0 );
            getHist("QMisidVsPtPrim")->Fill( mct.p.Pt(), mct.q == fwdTrack->charge() ? 0 : 1 );
        } else {
            getHist( "curveResolutionGlobal" )->Fill( dCurve );
            getHist( "curveResolutionVsPtGlobal" )->Fill( mct.p.Pt(), dCurve );
            getHist("QidVsPtGlobal")->Fill( mct.p.Pt(), mct.q == fwdTrack->charge() ? 1 : 0 );
            getHist("QMisidVsPtGlobal")->Fill( mct.p.Pt(), mct.q == fwdTrack->charge() ? 0 : 1 );
        }
    }
    LOG_INFO << "<<ProcessFwdTracks" << endm;
}

// what do I want to know:
// 2D nFST vs. nFTT for each track type
// 2D: nFST vs. track type
// 2D: nFTT vs. track type
void StFwdFitQAMaker::FillTrackSeedHistograms( StFwdTrack *fwdTrack ){
    if ( fwdTrack->isGlobalTrack() ){
        getHist2( "GlobalSeedFstFtt" )->Fill( fwdTrack->mFTTPoints.size(), fwdTrack->mFSTPoints.size() );
    }
    else if ( fwdTrack->isBeamLineConstrainedTrack() ){
        getHist2( "BeamlineSeedFstFtt" )->Fill( fwdTrack->mFTTPoints.size(), fwdTrack->mFSTPoints.size() );
    }
    else if ( fwdTrack->isPrimaryTrack() ){
        getHist2( "PrimarySeedFstFtt" )->Fill( fwdTrack->mFTTPoints.size(), fwdTrack->mFSTPoints.size() );
    }
    else if ( fwdTrack->isFwdVertexConstrainedTrack() ){
        getHist2( "SecondarySeedFstFtt" )->Fill( fwdTrack->mFTTPoints.size(), fwdTrack->mFSTPoints.size() );
    }
    // now 2d vs type
    getHist2( "SeedFstType" )->Fill( fwdTrack->trackType(), fwdTrack->mFSTPoints.size() );
    getHist2( "SeedFttType" )->Fill( fwdTrack->trackType(), fwdTrack->mFTTPoints.size() );
    //Seed Epd type is filled in EventStats because it is not in StFwdTrack yet

    getHist( "SeedQaTruthType" )->Fill( fwdTrack->trackType(), fwdTrack->qaTruth() );
}

void StFwdFitQAMaker::FillEventStats(){
    // Get the StFwdTrackMaker 
    StFwdTrackMaker *fwdTrackMaker = static_cast<StFwdTrackMaker *>(GetMaker( "fwdTrack" ));
    if ( !fwdTrackMaker ) {
        LOG_WARN << "StFwdTrackMaker not found" << endm;
        return;
    }
    auto stats = fwdTrackMaker->GetEventStats();

    getHist("mNumSeeds") -> Fill( stats.mNumSeeds );
    getHist("mAttemptedFits") -> Fill( stats.mAttemptedFits );
    getHist("mGoodFits") -> Fill( stats.mGoodFits );
    getHist("mFailedFits") -> Fill( stats.mFailedFits );
    getHist("mAttemptedGlobalFits") -> Fill( stats.mAttemptedGlobalFits );
    getHist("mFailedGlobalFits") -> Fill( stats.mFailedGlobalFits );
    getHist("mGoodGlobalFits") -> Fill( stats.mGoodGlobalFits );
    getHist("mFailedGlobalRefits") -> Fill( stats.mFailedGlobalRefits );
    getHist("mGoodGlobalRefits") -> Fill( stats.mGoodGlobalRefits );
    getHist("mNumFwdVertices") -> Fill( stats.mNumFwdVertices );
    getHist("mAttemptedPrimaryFits") -> Fill( stats.mAttemptedPrimaryFits );
    getHist("mGoodPrimaryFits") -> Fill( stats.mGoodPrimaryFits );
    getHist("mFailedPrimaryFits") -> Fill( stats.mFailedPrimaryFits );
    getHist("mFailedPrimaryRefits") -> Fill( stats.mFailedPrimaryRefits );
    getHist("mGoodPrimaryRefits") -> Fill( stats.mGoodPrimaryRefits );
    getHist("mAttemptedBeamlineFits") -> Fill( stats.mAttemptedBeamlineFits );
    getHist("mGoodBeamlineFits") -> Fill( stats.mGoodBeamlineFits );
    getHist("mFailedBeamlineFits") -> Fill( stats.mFailedBeamlineFits );
    getHist("mFailedBeamlineRefits") -> Fill( stats.mFailedBeamlineRefits );
    getHist("mGoodBeamlineRefits") -> Fill( stats.mGoodBeamlineRefits );
    getHist("mAttemptedSecondaryFits") -> Fill( stats.mAttemptedSecondaryFits );
    getHist("mGoodSecondaryFits") -> Fill( stats.mGoodSecondaryFits );
    getHist("mFailedSecondaryFits") -> Fill( stats.mFailedSecondaryFits );
    getHist("mFailedSecondaryRefits") -> Fill( stats.mFailedSecondaryRefits );
    getHist("mGoodSecondaryRefits") -> Fill( stats.mGoodSecondaryRefits );

    for( auto elem : stats.mGlobalNumEpdFoundHits){
        getHist( "SeedEpdType" ) -> Fill( 0.5, elem );
    }
    for( auto elem : stats.mBeamlineNumEpdFoundHits){
        getHist( "SeedEpdType" ) -> Fill( 1.5, elem );
    }
    for( auto elem : stats.mPrimaryNumEpdFoundHits){
        getHist( "SeedEpdType" ) -> Fill( 2.5, elem );
    }
    for( auto elem : stats.mSecondaryNumEpdFoundHits){
        getHist( "SeedEpdType" ) -> Fill( 3.5, elem );
    }

    for( auto elem : stats.numGlobalFoundHits){
        getHist( "numGlobalFoundHits" ) -> Fill( elem );
    }
    for( auto elem : stats.numBeamlineFoundHits){
        getHist( "numBeamlineFoundHits" ) -> Fill( elem );
    }
    for( auto elem : stats.numPrimaryFoundHits){
        getHist( "numPrimaryFoundHits" ) -> Fill( elem );
    }
    for( auto elem : stats.numSecondaryFoundHits){
        getHist( "numSecondaryFoundHits" ) -> Fill( elem );
    }
    for( auto elem : stats.mStep1Duration){
        getHist( "mStep1Duration" ) -> Fill( elem );
    }
    for( auto elem : stats.mSeedFindingDuration){
        getHist( "mSeedFindingDuration" ) -> Fill( elem );
    }
    for( auto elem : stats.mStep2Duration){
        getHist( "mStep2Duration" ) -> Fill( elem );
    }
    for( auto elem : stats.mStep3Duration){
        getHist( "mStep3Duration" ) -> Fill( elem );
    }
    for( auto elem : stats.mStep4Duration){
        getHist( "mStep4Duration" ) -> Fill( elem );
    }
    for( auto elem : stats.mGlobalFitDuration){
        getHist( "mGlobalFitDuration" ) -> Fill( elem );
    }
    for( auto elem : stats.mBeamlineFitDuration){
        getHist( "mBeamlineFitDuration" ) -> Fill( elem );
    }
    for( auto elem : stats.mPrimaryFitDuration){
        getHist( "mPrimaryFitDuration" ) -> Fill( elem );
    }
    for( auto elem : stats.mSecondaryFitDuration){
        getHist( "mSecondaryFitDuration" ) -> Fill( elem );
    }
}

