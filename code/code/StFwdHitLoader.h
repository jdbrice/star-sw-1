#ifndef ST_FWD_HIT_LOADER_H
#define ST_FWD_HIT_LOADER_H

#ifndef __CINT__
#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StFwdTrackMaker/include/Tracker/FwdDataSource.h"
#include "StEvent/StFstConsts.h"   // kFstrStart[], kFstStripPitchR for rasterizer fix
#else
class FwdHit;
#endif

class StEvent;
class StMuDstMaker;
class St_g2t_fts_hit;
class StFcsDb;


class FstRasterizer {
  public:
    FstRasterizer( double r = 3.0, double phi = 0.0040906154 ) {
        mRasterR = r;
        mRasterPhi = phi;
    }
    void setRPhi( double r, double phi ) {
        mRasterR = r;
        mRasterPhi = phi;
    }
    ~FstRasterizer() {}

    TVector3 raster(TVector3 p0) {
        TVector3 p = p0;
        double r = p.Perp();
        double phi = p.Phi();
        const double minR = 5.0;
        //AAA NOTE: mRasterR=3.0 cm does NOT match kFstStripPitchR=2.875 cm, causing position
        //  biases of 0.06–0.94 cm per strip (growing outward).  Fixing rasterR to 2.875 cm
        //  and sigma to pitch/sqrt12 made things WORSE (hit-search windows too tight with
        //  GEANT path).  Root fix: enable RunFstChain=true and use StFstFastSimMaker STEVENT
        //  path which has correct strip centers (from RSegment[]) and positionError=dr/sqrt12.
        //  Leave rasterR=3.0 cm for now to keep GEANT path stable.
        p.SetPerp(minR + (std::floor((r - minR) / mRasterR) * mRasterR + mRasterR / 2.0));
        p.SetPhi(-TMath::Pi() + (std::floor((phi + TMath::Pi()) / mRasterPhi) * mRasterPhi + mRasterPhi / 2.0));
        return p;
    }
    double mRasterR, mRasterPhi;
};

// Per-wedge phi correction for FST hits loaded from MuDst -- see
// bugreport_StFstHitMaker.txt: the official reconstruction (StFstHitMaker.cxx)
// computes each hit's global position from idealized per-wedge geometry
// constants only, discarding the real per-sensor DB/survey alignment
// correction (StFstDb::getRotations()) before it's ever applied. Since we
// aren't reprocessing raw DAQ data to fix that at the source, this applies an
// empirically-measured per-(disk, wedge) correction at MuDst-load time
// instead -- downstream of production, upstream of track fitting, exactly
// where FstRasterizer (above) already does a similar load-time correction for
// a different issue.
//
// Correction values come from the mean r*dphi (unbiased, hit-removed
// residual) in each of 12 fixed 30-degree phi bins per disk, measured over
// the full 150-file BLCVtx alignment campaign (~117M rows) -- see
// script/computeWedgeCorrection.C. Binned by the hit's own measured phi
// rather than by StFstConsts.h's wedge-index labeling (kFstphiStart etc.):
// each of the 12 phi bins found in data is exactly one physical wedge
// regardless of what numeric wedge index label it carries internally, so
// this sidesteps needing to re-derive that (non-monotonic) mapping.
class FstWedgeAligner {
  public:
    // wedgePhiCorrection[disk][wedgePhiBin], radians -- corrected_phi = phi - correction.
    // Bin 0 covers phi in [-180,-150) deg; bin index increases by 30 deg per
    // step (bin = floor(((phiDeg+180) mod 360) / 30)). Boundary phase
    // confirmed empirically 2026-07-07 by finding the actual step-jump
    // locations in the fine-binned (4.5 deg) FST0 r*dPhi-vs-phi profile from
    // the 150-file baseline: jumps cluster at 0, +-30, +-60, +-90, +-120,
    // +-150, +-180 deg, NOT at the +-15-deg-offset boundaries an earlier
    // version of this code assumed -- that +15 deg (half-wedge) phase was a
    // bug, silently smearing each "wedge" bin across two real physical
    // wedges' edges. Measured 2026-07-07 from the full 150-file BLCVtx
    // alignment campaign (~117M FST rows) via script/computeWedgeCorrection.C
    // -- rerun that script and paste in fresh numbers as more stats come in.
    // Declared here, defined in StFwdHitLoader.cxx -- CINT's rootcint (this
    // header gets a dictionary generated) doesn't understand constexpr for a
    // static array member, so this uses the classic out-of-line-definition
    // pattern instead.
    static const double kWedgePhiCorrection[3][12];

    static int wedgePhiBin(double phi) {
        double phiDeg = phi * 180.0 / TMath::Pi();
        double shifted = phiDeg + 180.0;
        while (shifted < 0)     shifted += 360.0;
        while (shifted >= 360)  shifted -= 360.0;
        return (int)(shifted / 30.0);
    }

    // diskIndex: 0-2. Returns the corrected phi (radians, still in (-pi,pi] after folding).
    static double correctPhi(int diskIndex, double phi) {
        if (diskIndex < 0 || diskIndex > 2) return phi;
        int bin = wedgePhiBin(phi);
        double corrected = phi - kWedgePhiCorrection[diskIndex][bin];
        while (corrected >  TMath::Pi()) corrected -= 2*TMath::Pi();
        while (corrected <= -TMath::Pi()) corrected += 2*TMath::Pi();
        return corrected;
    }
};

class StFwdHitLoader {
 public:
    StFwdHitLoader() : 
    mFttDataSource(StFwdHitLoader::DataSource::STEVENT), 
    mFstDataSource(StFwdHitLoader::DataSource::STEVENT), 
    mEpdDataSource(StFwdHitLoader::DataSource::STEVENT) 
    {}
    ~StFwdHitLoader() {}
    void clear() {
      #if !defined (__CINT__)
        mFwdHitsFtt.clear();
        mFwdHitsFst.clear();
        mFwdHitsEpd.clear();
      #endif

        // clear vectors for visualization OBJ hits
        mSpacepointsFtt.clear();
        mSpacepointsFst.clear();
        mSpacepointsEpd.clear();
    }

  #if !defined(__CINT__) && !defined(__CLING__)
    /********************************/
    // Load hits from the FTT detector (location based on mFttDataSource)
    //  * @param mcTrackMap : map of mc tracks
    //  * @param hitMap : FTT hitmap to populate
    //  * @return number of hits loaded
    //  * @note: This function is called by StFwdTrackMaker::Make()
    int loadFttHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );
    int loadFttPointsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );
    int loadFttPointsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );

    /********************************/
    // Load hits from the FST detector (location based on mFstDataSource)
    //  * @param mcTrackMap : map of mc tracks
    //  * @param hitMap : FST hitmap to populate
    //  * @return number of hits loaded
    //  * @note: This function is called by StFwdTrackMaker::Make()
    int loadFstHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );
    int loadFstHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );
    int loadFstHitsFromMuDst( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );
    int loadFstHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );

    /********************************/
    // Load hits from the EPD detector (location based on mEpdDataSource)
    //  * @param mcTrackMap : map of mc tracks
    //  * @param hitMap : EPD hitmap to populate
    //  * @return number of hits loaded
    //  * @note: This function is called by StFwdTrackMaker::Make()
    int loadEpdHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, StFcsDb *fcsDb = nullptr );
    int loadEpdHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, StFcsDb *fcsDb = nullptr );
  #endif

    // all caps to avoid conflict with types themselves
    enum DataSource { GEANT=0, STEVENT, MUDST, IGNORE };
    /********************************/
    // Specify the data source for FTT, FST and EPD hits
    //  * @param ds : data source (GEANT, StEvent, MuDst)
    //  * @return None
    void setFttDataSource( DataSource ds ) { mFttDataSource = ds; }
    void setFstDataSource( DataSource ds ) { mFstDataSource = ds; }
    void setEpdDataSource( DataSource ds ) { mEpdDataSource = ds; }
    void setDataSource( DataSource ds ) { mFttDataSource = ds; mFstDataSource = ds; mEpdDataSource = ds; }

    void setStEvent( StEvent *stEvent ) { mStEvent = stEvent; }
    void setMuDstMaker( StMuDstMaker *muDstMaker ) { mMuDstMaker = muDstMaker; }
    void setTables( St_g2t_fts_hit *stg_hits, St_g2t_fts_hit *fst_hits, St_g2t_fts_hit *epd_hits ) {
        mGeantFtt = stg_hits;
        mGeantFst = fst_hits;
        mGeantEpd = epd_hits;
    }
    FstRasterizer mFstRasterizer;

    // Off by default -- opt-in per-run via setApplyFstWedgeAlignment(true), same
    // convention as this session's other experimental additions (e.g.
    // StFwdAlignmentMaker's runFwdAlignment). See FstWedgeAligner above and
    // bugreport_StFstHitMaker.txt for why this exists: the official
    // reconstruction never applies the real per-sensor DB alignment, so this
    // corrects for the resulting per-wedge offset at MuDst-load time instead.
    void setApplyFstWedgeAlignment( bool apply ) { mApplyFstWedgeAlignment = apply; }
    bool mApplyFstWedgeAlignment = false;
  protected:
    DataSource mFttDataSource;
    DataSource mFstDataSource;
    DataSource mEpdDataSource;
    float mEpdThreshold = 0.2;

    // Pointers to these are used by StFwdTrackMaker, clear the vectors after each event
    #if !defined (__CINT__)
    vector<FwdHit> mFwdHitsFtt;
    vector<FwdHit> mFwdHitsFst;
    vector<FwdHit> mFwdHitsEpd;

    // this disables logging at compile time
    constexpr static int kLogVerbose = 10;
    constexpr static int kLogInfo = 1;
    constexpr static int kLogSilent = 0;
    constexpr static int kLogLevel = kLogSilent;
    #endif

    vector<TVector3> mSpacepointsFtt;
    vector<TVector3> mSpacepointsFst;
    vector<TVector3> mSpacepointsEpd;

    /// Non-owning pointer: StEvent provided externally via `setStEvent`
    StEvent *mStEvent = nullptr; // pointer to StEvent
    /// Non-owning pointer: StMuDstMaker provided externally via `setMuDstMaker`
    StMuDstMaker *mMuDstMaker = nullptr; // pointer to StMuDstMaker

    /// Non-owning pointers to GEANT hits
    St_g2t_fts_hit *mGeantFtt = nullptr; // pointer to GEANT FTT hits
    St_g2t_fts_hit *mGeantFst = nullptr; // pointer to GEANT FST hits
    St_g2t_fts_hit *mGeantEpd = nullptr; // pointer to GEANT EPD hits

    

    bool verbosity = 0;
};


#endif // ST_FWD_HIT_LOADER_H
