#ifndef ST_FWD_RESIDUAL_MAKER_H
#define ST_FWD_RESIDUAL_MAKER_H

#include "StChain/StMaker.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

// Complex headers guarded from CINT
#ifndef __CINT__
#include <vector>
#endif

class StFwdResidualMaker : public StMaker {
public:
    // makerName must be unique if multiple StFwdResidualMaker instances (one
    // per track type) are added to the same chain -- StMaker's own name
    // defaults the same for every instance otherwise, which StChain doesn't
    // tolerate well (name-keyed maker lookups only find the first one).
    StFwdResidualMaker(const char* outFile = "fwdDetResidual.root", const char* makerName = "fwdResidual");
    virtual ~StFwdResidualMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    // Track type used for FST/FTT residuals (StFwdTrack::StFwdTrackType):
    // 0=Global 1=BLC 2=Primary 3=FwdVtx 4=BLCVtx 5=FCSConstrained
    // Global (default) has no vertex constraint -> least biased residuals for MC.
    // For real data with a good TPC vertex (or B=0 straight-track runs), BLC or
    // Primary may be preferable. Plane-usage histograms (all 6 types) are
    // unaffected by this setting.
    void setTrackType(UChar_t t) { mResidualTrackType = t; }

    // Data source for forward tracks:
    //  false (default) = StEvent::fwdTrackCollection() -- freshly refit tracks
    //                     from this event (needed by the afterburner: MuDst-level
    //                     StMuFwdTrackCollection is only ever populated by
    //                     StMuDstMaker's write-mode fillFwdTrack(), which the
    //                     afterburner never runs, so it would otherwise be stale
    //                     production-time data or empty).
    //  true             = StMuDst::muFwdTrackCollection() -- for analyzing an
    //                     already-produced MuDst directly (no live StEvent/refit
    //                     chain), where that collection genuinely was filled by
    //                     a normal write-mode BFC production.
    void setUseMuDst(bool b = true) { mUseMuDst = b; }

    ClassDef(StFwdResidualMaker, 1)

private:
    TString mOutFile;
    TFile*  mFout = nullptr;
    UChar_t mResidualTrackType = 0;   // Global by default
    bool    mUseMuDst = false;        // false = read StEvent (see setUseMuDst)

    // Counters
    int mNGoodTracks = 0;
    int mNFstRes     = 0;
    int mNFttRes     = 0;

    // FST histograms: [disk 0-2][0=vs x,1=vs y,2=vs r,3=vs rphi,4=vs phi]
    // vs phi (raw, not r*phi) is what a transverse (x,y) misalignment shows up
    // in -- see debug/index.html Issue #28 -- so it's kept separate from vs rphi.
    TH1F* hFst[3];
    TH2F* h2Fst[3][5];

    // FTT horizontal (dy): [plane 0-3]
    TH1F* hFttDy[4];
    TH2F* h2FttDy[4][4];

    // FTT vertical (dx): [plane 0-3]
    TH1F* hFttDx[4];
    TH2F* h2FttDx[4][4];

    // Plane-usage / "efficiency" histogram, one per StFwdTrack::StFwdTrackType (0-5):
    // 20 bins — 0:Vtx, 1-3:FST1-3, 4-6:FTT1(x,y,uv), 7-9:FTT2, 10-12:FTT3, 13-15:FTT4,
    // 16:EPD(not impl.), 17:FCSEcal, 18:FCSHcal, 19:AllTrk (always filled)
    TH1F* hPlaneUsage[6];

    // Same 20-bin scheme as hPlaneUsage, but split by quadrant and charge sign,
    // mResidualTrackType only (keeps the count reasonable -- this is meant for a
    // targeted look at the FTT disk2 south-half dead-zone / charge-asymmetry
    // question, not routine per-type bookkeeping like hPlaneUsage above).
    // quadrant: 0=N-top(x<0,y>0) 1=N-bottom(x<0,y<=0) 2=S-top(x>=0,y>0) 3=S-bottom(x>=0,y<=0)
    // charge:   0=positive(q>=0) 1=negative(q<0)
    // Quadrant is from the mean (x,y) of all FST+FTT seed points used by the track.
    TH1F* hPlaneUsageQC[4][2];

    // XY of the hit actually used in the fit, at each disk/orientation, split by
    // charge sign, mResidualTrackType only. FST: disk 0-2, no orientation split
    // (single 2D-strip measurement per disk). FTT: disk 0-3 x orientation
    // (0=x/vertical-strip, 1=y/horizontal-strip, 2=uv/diagonal).
    // 3 FST + 4*3 FTT = 15 plots x 2 charges = 30 total.
    TH2F* h2FstHitXY[3][2];
    TH2F* h2FttHitXY[4][3][2];

    // Same gating (a hit at that disk/orientation was actually used in the fit)
    // and same binning as h2Fst/FttHitXY above, but filled with the track's
    // continuous projected (x,y) at that disk/plane (fstProjs[disk]/fttProjs
    // [plane], already computed for findClosestZ matching) instead of the hit's
    // own position. An FTT strip only measures one coordinate precisely; the
    // hit's stored position along the strip's long axis is a fixed/quantized
    // placeholder, which made h2Fst/FttHitXY show banding along that axis
    // instead of the real, continuous track footprint. This is the fix.
    TH2F* h2FstProjXY[3][2];
    TH2F* h2FttProjXY[4][3][2];

    void bookHistos();

    void fillFst(int disk, float res, float hx, float hy);
    void fillFttDy(int plane, float res, float hx, float hy);
    void fillFttDx(int plane, float res, float hx, float hy);

#ifndef __CINT__
    // Neutral, source-agnostic copies of the pieces of a track this maker needs.
    // Filled from either StMuFwdTrack (MuDst) or StFwdTrack (StEvent) at the top
    // of Make(), so the residual/plane-usage logic below is written -- and
    // verified -- exactly once, for both sources.
    struct FwdProj { float x = 0, y = 0, z = 0; };
    struct FwdSeedPt { float x = 0, y = 0, z = 0; float cov[9] = {0}; };

    Int_t makeFromStEvent();
    Int_t makeFromMuDst();

    void processTrack(UChar_t trackType, char charge, bool hasEcal, bool hasHcal,
                       const std::vector<FwdSeedPt>& fstPts,
                       const std::vector<FwdSeedPt>& fttPts,
                       const std::vector<FwdProj>& fstProjs,
                       const std::vector<FwdProj>& fttProjs);

    void processFstPoints(const std::vector<FwdSeedPt>& fstPts,
                          const std::vector<FwdProj>& fstProjs);
    void processFttPoints(const std::vector<FwdSeedPt>& fttPts,
                          const std::vector<FwdProj>& fttProjs);
    int  findClosestZ(float hz, const std::vector<FwdProj>& projs, int maxIdx);
    void fillPlaneUsage(UChar_t trackType, char charge, bool hasEcal, bool hasHcal,
                        const std::vector<FwdSeedPt>& fstPts,
                        const std::vector<FwdSeedPt>& fttPts,
                        const std::vector<FwdProj>& fstProjs,
                        const std::vector<FwdProj>& fttProjs);
#endif
};

#endif
