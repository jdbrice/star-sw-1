/***************************************************************************
 * StFttSlowSimMaker.h
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
#ifndef STFTTSLOWSIMMAKER_H
#define STFTTSLOWSIMMAKER_H

#include "StMaker.h"

#ifndef __CINT__
#include <map>
#endif

class StEvent;
class StFttCollection;
class StFttDb;

class StFttSlowSimMaker : public StMaker {
public:
    StFttSlowSimMaker( const char* name = "fttSlowSim" );
    ~StFttSlowSimMaker();

    int Init();
    int Make();

    void SetDebug( int v = 1 ) { mDebug = v; }

    // Charge model (see ftt_sim_hit_maker/index.html#decharge for the
    // measurement behind it): the total charge of one crossing is
    //
    //     totalAdc = mAdcPerKeV * (GEANT de, summed over the crossing's
    //                              steps) * gainFluctuation()
    //
    // i.e. driven by the energy the particle actually deposited, not drawn
    // from a fixed distribution. In a 3mm gas gap the Landau dE/dx spread
    // dominates and the gas amplification is only a correction -- measured
    // against real p+p clusters, de alone (mGainSpread=0) already matches
    // the spectrum shape better than the exponential this used to draw, and
    // mGainSpread~0.3 is the optimum; larger values get worse again.
    //
    // mAdcPerKeV is anchored on the cluster PEAK-strip ADC, not on sumAdc:
    // the peak strip is the center strip, which this +-1-strip model
    // represents faithfully, whereas 20% of real clusters have 4-7 strips
    // that it cannot produce at all (shallow-incidence tracks crossing many
    // strips). With one parameter the peak-strip mean AND rms both land on
    // the real values; sumAdc then comes out ~15% high, which is that same
    // missing population. See ftt_sim_hit_maker/index.html#decharge.
    // Drawn per crossing in TWO parts, common and per-view -- see
    // setGainSpreadIndependent() -- then split over the true
    // hit's own center strip +-1 via a Gaussian(sigma=mChargeShareSigma) response
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
    void setAdcPerKeV( float adcPerKeV )   { mAdcPerKeV = adcPerKeV; }
    void setGainSpread( float relWidth )   { mGainSpread = relWidth; }  // 0 = no gas amplification spread

    // How much of the single-view spread above is INDEPENDENT between the two
    // strip layers of a foil, the rest being common to both.
    //
    // The two layers read the same avalanche, so the naive expectation is that
    // all of it is common and the measured X/U ratio should be a delta. It is
    // not: in 391k/465k clean 1+1 pairs from 682k Run 23 ZF events, the ratio's
    // lognormal CORE has a relative width of 0.207 (U/X) and 0.213 (V/Y),
    // i.e. 0.21 +- 0.02 including fit-window and background-model systematics.
    // A Gaussian core fit gives 0.207 as well, so the number does not depend on
    // the model. Splitting that equally between the two views (the ratio alone
    // cannot say which view is noisier) gives
    //     sigma_independent = 0.21 / sqrt(2) = 0.15
    // and keeping the single-view total at the tuned mGainSpread then fixes
    //     sigma_common = sqrt(0.30^2 - 0.15^2) = 0.26.
    // So the per-view ADC spectrum mAdcPerKeV/mGainSpread were tuned against is
    // preserved exactly, while X-U goes from perfectly correlated to the
    // measured 0.21. NOTE the RMS of the whole ratio distribution is 0.47, which
    // would imply 0.33 -- that is the ~45% background and tail, not gain.
    // See ftt_sim_hit_maker/charge_sharing.html#xusplit.
    void setGainSpreadIndependent( float relWidth ) { mGainSpreadIndependent = relWidth; }
    void setStripAdcThreshold( float adc ) { mStripAdcThreshold = adc; }

private:
    static bool computeRowStrip( UChar_t orientation, double local_x, double local_y, int &row, int &centerStrip, double &subStripOffset );
    double gainFluctuation( double relWidth );
    // The part of mGainSpread left once the independent part is taken out.
    // The two draws MULTIPLY, so their relative widths combine as
    // (1+t^2) = (1+c^2)(1+i^2), not in quadrature -- solving that for c keeps
    // the single-view total at exactly mGainSpread instead of ~1% over it.
    double gainSpreadCommon() const {
        double t = mGainSpread, i = mGainSpreadIndependent;
        if ( i >= t ) return 0.0;              // caller asked for all of it per-view
        return sqrt( ( t*t - i*i ) / ( 1.0 + i*i ) );
    }
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
    float mAdcPerKeV         = 600.0;  // tuned through the full chain so the cluster PEAK-strip ADC
                                       // matches real data (mean 227, rms 199 vs 226/196)
    float mGainSpread        = 0.3;    // TOTAL relative width of the single-view gain spread; 0 = pure de scaling
    float mGainSpreadIndependent = 0.15; // the part of it independent per strip layer; measured, see the setter
    float mStripAdcThreshold = 70.0;   // ~ observed hardware zero-suppression edge

    ClassDef( StFttSlowSimMaker, 0 )
};

#endif // STFTTSLOWSIMMAKER_H
