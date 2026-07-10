#ifndef ST_FWD_ALIGNMENT_MAKER_H
#define ST_FWD_ALIGNMENT_MAKER_H

// StFwdAlignmentMaker -- unbiased ("hit-removed") FST/FTT residuals.
//
// For each good track, and for each FST/FTT hit on it in turn: build a seed
// with that one hit removed, refit, project the refit trajectory to the
// removed hit's plane, and record hit-vs-projection. Removing the hit from
// the fit before measuring its residual is what makes this an UNBIASED
// residual -- StFwdResidualMaker's residuals are BIASED (the hit is part of
// the fit that produced the projection being compared against it), which
// systematically pulls the fit toward the hit and underestimates any real
// misalignment. See proposal_alignment_path.txt for the full writeup.
//
// This is deliberately NOT a replacement for StFwdResidualMaker -- it costs
// roughly (num FST+FTT planes) extra refits per track, so it is not meant to
// run in routine production. It is also deliberately independent of it: this
// maker does not touch StEvent/StMuDst/StPicoDst, does not modify anything
// StFwdTrackMaker produces, and writes its own standalone ROOT TTree. It only
// *reads* the track results StFwdTrackMaker already computed for this event
// (via the small additive GetForwardTrackerBase() accessor added to
// StFwdTrackMaker for this purpose -- see StFwdTrackMaker.h/.cxx) and
// performs additional, independent refits alongside them.
//
// Track types: works for any StFwdTrack::StFwdTrackType via setTrackType().
// BLC/Primary/BLCVtx/FCSConstrained types include a primary-vertex spacepoint
// in their seed; this maker never removes that point (it isn't a detector
// plane being aligned) -- the vertex hit always carries detid==kTpcId
// (confirmed for BLCVtx: FwdTracker.h's doBLCVertexFitting() constructs it via
// mBLCVtxHit.setXYZDetId(..., kTpcId)), which the removal loop's
// isFst()||isFtt() test already excludes regardless of track type, so this
// needed no code change to extend past Global.
//
// The caveat worth keeping in mind for vertex-constrained types (BLCVtx etc.):
// if the vertex itself was built from very few tracks, removing one track's
// FST/FTT hit and refitting *that* track doesn't change the (already-fixed)
// vertex position, but the vertex was itself derived partly from that hit
// through that track's original fit -- a secondary, indirect leakage, diluted
// by however many tracks/events went into the vertex fit. Global tracks have
// no vertex constraint at all, so remain the cleanest case; BLCVtx numbers
// should be read with this in mind, not as fully independent of the hit being
// tested.
//
// Usage (added to the same chain as StFwdTrackMaker, after it):
//   StFwdAlignmentMaker* align = new StFwdAlignmentMaker("fwdAlignment.root");
//   align->setTrackMaker(fwdTrack);      // required
//   chain->AddAfter("fwdTrack", align);
//
// Off by default in the sense that nothing runs unless a driver macro
// explicitly creates and adds this maker -- see the runFwdAlignment flag in
// fwd_afterburner_db.C and script/sim.C.

#include "StChain/StMaker.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

class StFwdTrackMaker;
class ForwardTrackMaker;

class StFwdAlignmentMaker : public StMaker {
public:
    StFwdAlignmentMaker(const char* outFile = "fwdAlignment.root",
                         const char* makerName = "fwdAlignment");
    virtual ~StFwdAlignmentMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    // Required: gives this maker read-only access to the already-fitted
    // tracks/seeds for this event, via StFwdTrackMaker's own tracker
    // instance. Must be called before the chain runs (right after both
    // makers are constructed in the driver macro).
    void setTrackMaker(StFwdTrackMaker* mk) { mFwdTrackMaker = mk; }

    // Track type to study (StFwdTrack::StFwdTrackType: 0=Global 1=BLC
    // 2=Primary 3=FwdVtx 4=BLCVtx 5=FCSConstrained). See class-level note
    // above for the vertex-leakage caveat on non-Global types.
    void setTrackType(UChar_t t) { mTrackType = t; }

    // Skip the leave-one-out test for a given hit if fewer than this many
    // FST+FTT hits would remain afterward -- keeps refits away from GenFit's
    // convergence floor. Default 4 (STAR forward tracks have at most 7
    // FST+FTT points: 3 FST + 4 FTT).
    void setMinRemainingHits(int n) { mMinRemainingHits = n; }

    ClassDef(StFwdAlignmentMaker, 1)

private:
    TString mOutFile;
    TFile*  mFout = nullptr;
    TTree*  mTree = nullptr;

    StFwdTrackMaker* mFwdTrackMaker = nullptr;

    UChar_t mTrackType         = 0;  // Global by default
    int     mMinRemainingHits  = 4;

    // TTree branch buffers -- one row per (track, plane-tested) hit. See
    // proposal_alignment_path.txt section 3(c) for why these columns and not
    // pre-computed deltas.
    Int_t    b_run             = 0;
    Int_t    b_event           = 0;
    Int_t    b_trackId         = 0;
    UChar_t  b_trackType       = 0;
    Int_t    b_detType         = 0;      // 0=FST 1=FTT
    Int_t    b_genfitPlaneIndex = 0;     // raw sensor/plane index (see getPlaneFor()); grouping into
                                          // FST disk 0-2 / FTT plane 0-3 is left to analysis, via hitZ
    Int_t    b_stripDir        = 0;      // FTT only: 0=undetermined/FST, 1=x-strip, 2=y-strip
    Float_t  b_hitX = 0, b_hitY = 0, b_hitZ = 0;    // measured, global frame
    Float_t  b_projX = 0, b_projY = 0, b_projZ = 0; // unbiased refit projection, global frame
    Float_t  b_chi2ndf      = 0;
    Int_t    b_nPointsUsed  = 0;
    Bool_t   b_converged    = false;

    Long64_t mNTracksSeen    = 0;
    Long64_t mNTracksTested  = 0;
    Long64_t mNRowsWritten   = 0;

    void bookTree();

#ifndef __CINT__
    // Runs the leave-one-out loop for one already-fitted track, using the
    // SAME ForwardTrackMaker/TrackFitter instance StFwdTrackMaker built for
    // the nominal fit (via GetForwardTrackerBase()). Fills mTree rows for
    // every hit that clears the min-remaining-hits guard and refits/projects
    // successfully. Declared out-of-line (not inline in the header) because
    // it needs GenfitTrackResult/Seed_t/ForwardTrackMaker, which are CINT-
    // unfriendly types kept out of this header on purpose (same pattern
    // StFwdResidualMaker already uses for its own physics core).
    void processTrack(class GenfitTrackResult& gtr, ForwardTrackMaker* ft,
                       int runId, int eventId, int trackId);
#endif
};

#endif
