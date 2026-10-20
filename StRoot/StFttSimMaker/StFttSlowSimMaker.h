/***************************************************************************
 * StFttSimHitMaker.h
 ***************************************************************************
 *
 * Description: GEANT-truth -> StFttRawHit injector for testing the REAL
 * StFttClusterMaker/StFttClusterPointMaker/StFttDb reconstruction chain with
 * single-particle-gun (or any GEANT) MC, instead of the GEANT-truth blur
 * used by StFttClusterPointMaker::MakeGeantPoints() / StFttFastSimMaker.
 *
 * Sits where StFttRawHitMaker + StFttHitCalibMaker normally would (both are
 * real-DAQ/calibration specific and are skipped for simulated input): reads
 * St_g2t_fts_hit truth hits, inverts StFttDb::getGloablOffset_ClusterPoint()
 * to get the local (row, strip) address, and injects synthetic StFttRawHits
 * (3-strip charge-sharing cluster) via StFttDb::reverseHardwareMap() so the
 * unmodified, real StFttClusterMaker does the actual clustering.
 *
 * See proposal_ftt_sim_maker.txt (section 6) and status_ftt_sim_maker.txt
 * for the full design/derivation and current status.
 ***************************************************************************/
#ifndef STFTTSIMHITMAKER_H
#define STFTTSIMHITMAKER_H

#include "StMaker.h"

#ifndef __CINT__
#include <map>
#endif

class StEvent;
class StFttCollection;
class StFttDb;

class StFttSimHitMaker : public StMaker {
public:
    StFttSimHitMaker( const char* name = "fttSimHit" );
    ~StFttSimHitMaker();

    int Init();
    int Make();

    void SetDebug( int v = 1 ) { mDebug = v; }

    // Real-data-based charge sharing (status_ftt_sim_maker.txt item 10):
    // per-hit total charge drawn ONCE per truth hit (shared between the
    // intended and diagonal orientation attempts -- they're the same
    // physical foil's two strip layers reading the same charge deposit)
    // from an exponential with mean=mMeanTotalAdc (measured sumAdc
    // mean=RMS=~430-450, exponential-consistent), split over the true hit's
    // own center strip +-1 via a Gaussian(sigma=mChargeShareSigma) response
    // integrated over each strip's bounds (erf). Each candidate strip's
    // SUMMED (across all truth hits landing on it) ADC is kept only if it
    // clears mStripAdcThreshold -- EXCEPT the center (d==0) strip, which is
    // always kept, floored up to mStripAdcThreshold if the draw came in low:
    // a naive threshold on all 3 strips would silently drop ~30% of truth
    // hits entirely (low-Q exponential draws), an unintended efficiency loss
    // this project deliberately keeps at 100% for now (see status page TODO)
    // -- only the neighbor strips' presence/absence should vary with Q and
    // sub-strip position, giving realistic 1/2/3-strip multiplicity. Not
    // meant to reproduce the measured profile beyond +-1 strip (that tail is
    // small and likely pileup/noise-dominated in the busy real sample this
    // was tuned from).
    void setChargeShareSigma( float mm )   { mChargeShareSigma = mm; }
    void setMeanTotalAdc( float adc )      { mMeanTotalAdc = adc; }
    void setStripAdcThreshold( float adc ) { mStripAdcThreshold = adc; }

private:
    static bool computeRowStrip( UChar_t orientation, double local_x, double local_y, int &row, int &centerStrip, double &subStripOffset );
    void addStrip( int plane, int quad, int row, int strip, float adc, UShort_t idTruth, UChar_t expectOrientation );

#ifndef __CINT__
    // A physical channel (plane,quad,row,strip) can receive contributions
    // from more than one truth hit in the same event -- two nearby tracks,
    // or one track's center strip landing on another's 10% neighbor spill
    // -- important for high-occupancy/embedding studies (akio). ADC must be
    // SUMMED per channel, exactly once per real hardware address, not left
    // as separate StFttRawHit records at the same (feb,vmm,ch): a real DAQ
    // channel reports one ADC value per event, and StFttClusterMaker's
    // clustering isn't written to expect (nor reliably tolerate) duplicate
    // entries at the same strip -- see status_ftt_sim_maker.txt item 9.
    struct ChannelAccum {
        int plane = 0, quad = 0, row = 0, strip = 0;
        UChar_t orientation = 4; // kFttUnknownOrientation -- StEnumerations.h not included here
        float adcSum = 0;
        float maxAdc = -1;     // to pick the dominant contributor's truth id
        UShort_t idTruth = 0;
        bool isCenter = false; // true if ANY truth hit this event called this its own d==0 strip -- see setStripAdcThreshold() comment
    };
    static long long packChannelKey( int plane, int quad, int row, int strip );
    void accumulateStrip( std::map<long long, ChannelAccum> &accum, int plane, int quad, int row, int strip, float adc, UShort_t idTruth, UChar_t orientation, bool isCenter );
#endif

    StEvent*         mEvent;         //! pointer to StEvent
    StFttCollection* mFttCollection; //! pointer to StFttCollection
    StFttDb*         mFttDb;         //! pointer to StFttDb

    int   mDebug = 0;
    float mChargeShareSigma  = 2.4;    // mm, fit to the measured nStrips==3 25:48:25 ratio
    float mMeanTotalAdc      = 431.5;  // measured cluster sumAdc mean
    float mStripAdcThreshold = 70.0;   // ~ observed hardware zero-suppression edge

    ClassDef( StFttSimHitMaker, 0 )
};

#endif // STFTTSIMHITMAKER_H
