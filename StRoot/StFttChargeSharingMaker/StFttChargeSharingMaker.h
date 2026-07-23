/***************************************************************************
 * StFttChargeSharingMaker.h
 ***************************************************************************
 *
 * Description: reads real (MuDst) FTT raw hits and measures the actual
 * per-strip charge-sharing profile within a cluster, to replace
 * StFttSimHitMaker's flat 1:8:1 3-strip model with something realistic.
 *
 * Time cut: StFttDb::getTimeCut() (checked directly against this run's DB)
 * returns mode=kHitCalibratedTime and a constant window [-65,100] for every
 * VMM -- i.e. the cut applies to time = dbcid - anchor(VMM), where the
 * per-VMM anchor is NOT in the DB (only the window is); real reconstruction
 * (StFttHitCalibMaker/HitCalibHelper) learns it online from the first
 * ~200 hits/VMM of the run. This maker reproduces that exact algorithm
 * (anchor = mode of that VMM's dbcid distribution) but over the whole
 * file's statistics instead of an 8-15-event bootstrap, then re-clusters
 * the surviving in-time hits itself (StFttClusterMaker's own algorithm:
 * max-ADC anchor + adjacent/non-increasing-ADC expansion) since MuDst's
 * StMuFttCluster<->StMuFttRawHit link is never populated
 * (StMuFttUtil::rebuildRelationships() is an empty stub) and StMuFttCluster
 * doesn't even carry a row() field.
 *
 * See status_ftt_sim_maker.txt item 10 for the full plan/derivation.
 ***************************************************************************/
#ifndef STFTTCHARGESHARINGMAKER_H
#define STFTTCHARGESHARINGMAKER_H

#include "StMaker.h"

#ifndef __CINT__
#include <map>
#include <vector>
#endif

class StMuDst;
class TFile;
class TH1F;
class TH2F;
class TProfile;

class StFttChargeSharingMaker : public StMaker {
public:
    StFttChargeSharingMaker( const char* name = "fttChargeSharing" );
    ~StFttChargeSharingMaker();

    int Init();
    int Make();
    int Finish();

    void setOutputFile( const char* fn ) { mOutFile = fn; }
    void setWarmupEvents( int n ) { mWarmupEvents = n; }
    void setTimeCutWindow( int lo, int hi ) { mTimeCutLo = lo; mTimeCutHi = hi; }
    void setMinSamplesForReady( int n ) { mMinSamplesForReady = n; }

private:
#ifndef __CINT__
    struct StripHit { int strip; float adc; };

    // per-VMM online time calibration, HitCalibHelper's algorithm (anchor =
    // mode of the dbcid distribution), just over the whole file instead of
    // an early bootstrap -- see header comment above.
    std::map<int, std::map<short,int> > mDbcidHist; // uuid -> {dbcid: count}
    std::map<int, short> mAnchor;                    // uuid -> cached anchor
    bool readyFor( int uuid );
    short anchorFor( int uuid );
    void fillCalib( int uuid, short dbcid );
    static int vmmUuid( int plane, int quad, int feb, int vmm );

    // re-clustering: SAME algorithm as StFttClusterMaker::FindClusters()/
    // SearchClusterEdges(), operating on a plain per-(plane,quad,row,orientation)
    // strip vector instead of StFttRawHit objects. outSums, if given, gets each
    // found cluster's sumAdc appended (for the same-foil combinatoric check).
    void clusterAndFill( std::vector<StripHit>& hits, std::vector<float>* outSums = 0 );
    static long long packGroupKey( int plane, int quad, int row, int orientation );
#endif

    TString mOutFile;
    TFile*  mFout;

    int mWarmupEvents;        // events used only to build up per-VMM anchors, not analyzed
    int mTimeCutLo;           // confirmed from this run's DB window
    int mTimeCutHi;
    int mMinSamplesForReady;  // matches HitCalibHelper::MIN_BCID_SAMPLES

    int mEventCount;
    int mNTotalHits;
    int mNInTimeHits;
    int mNClustersFound;

    TH1F*     mHNStrips;
    TH1F*     mHSumAdc;
    TH1F*     mHSumAdcN1;   // sumAdc for nStrips==1 clusters -- is the low-sumAdc peak isolated single-strip hits?
    TH1F*     mHSumAdcN2p;  // sumAdc for nStrips>=2 clusters
    TH1F*     mHPeakAdc;
    TH1F*     mHRawAdcInTime; // every individual in-time raw hit's own ADC, pre-clustering -- reveals any real per-channel threshold edge
    TProfile* mPProfileAll;
    TProfile* mPProfileByMult[6]; // index = min(nStrips,5); [0],[1] unused

    // Vertical+DiagonalV share the front foil, Horizontal+DiagonalH the back
    // foil (same physical charge deposit, read out by two differently-angled
    // strip layers) -- per (plane,quad,event) total in-time ADC on one layer
    // vs. the other, to check whether they correlate as expected for a
    // shared source. (Sums ALL hits per orientation this event/plane/quad --
    // a first pass; see mH2ClusterXvsU below for the real per-cluster check.)
    TH2F* mH2XvsU; // Vertical vs DiagonalV total ADC, per plane/quad/event
    TH2F* mH2YvsV; // Horizontal vs DiagonalH total ADC, per plane/quad/event

    // Per-CLUSTER version: every (X-cluster,U-cluster) combination within the
    // same plane/quad/event (and same for Y/V) -- true matches should show as
    // a correlated peak on top of the combinatorial background from
    // unrelated pairings, unlike the trivially-correlated event-total sums
    // above. Ratio histograms are the actual candidate cut variable.
    TH2F* mH2ClusterXvsU;
    TH2F* mH2ClusterYvsV;
    TH1F* mHRatioUX; // U/X cluster sumAdc ratio, all combinations
    TH1F* mHRatioVY; // V/Y cluster sumAdc ratio, all combinations

    // Exclusive (no combinatorial ambiguity): only plane/quad/events with
    // EXACTLY one X-cluster and one U-cluster (same for Y/V) -- isolates
    // genuine pairs to check whether a real correlation/ratio peak exists at
    // all before fighting the full-combinatorics background above.
    TH2F* mH2ClusterXvsU_excl;
    TH2F* mH2ClusterYvsV_excl;
    TH1F* mHRatioUX_excl;
    TH1F* mHRatioVY_excl;

    ClassDef( StFttChargeSharingMaker, 0 )
};

#endif // STFTTCHARGESHARINGMAKER_H
