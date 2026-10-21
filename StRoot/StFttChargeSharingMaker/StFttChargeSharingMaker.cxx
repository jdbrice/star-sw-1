#include "StFttChargeSharingMaker.h"

#include <algorithm>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuFttCollection.h"
#include "StMuDSTMaker/COMMON/StMuFttRawHit.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StFttRawHit.h"
#include "StFttDbMaker/StFttDb.h"

ClassImp( StFttChargeSharingMaker )

StFttChargeSharingMaker::StFttChargeSharingMaker( const char* name )
: StMaker( name ),
  mOutFile( "fttChargeSharing.root" ),
  mFout( nullptr ),
  mWarmupEvents( 20 ),
  mTimeCutLo( -65 ),
  mTimeCutHi( 100 ),
  mMinSamplesForReady( 200 ),
  mEventCount( 0 ),
  mNTotalHits( 0 ),
  mNInTimeHits( 0 ),
  mNClustersFound( 0 ), mNPQOSamples(0), mUseDb(false), mFttDb(nullptr), mNFromDbAnchor(0), mNFallback(0),
  mH2ClusterXvsU_cross(nullptr), mH2ClusterYvsV_cross(nullptr),
  mHRatioUX_cross(nullptr), mHRatioVY_cross(nullptr),
  mNPairsXU(0), mNPairsXUCross(0), mNPairsYV(0), mNPairsYVCross(0),
  mNExclXU(0), mNExclXUCross(0), mNExclYV(0), mNExclYVCross(0),
  mCrossTolerance(1.6),
  mHNStrips( nullptr ),
  mHSumAdc( nullptr ),
  mHSumAdcN1( nullptr ),
  mHSumAdcN2p( nullptr ),
  mHPeakAdc( nullptr ),
  mHRawAdcInTime( nullptr ),
  mPProfileAll( nullptr ),
  mH2XvsU( nullptr ),
  mH2YvsV( nullptr ),
  mH2ClusterXvsU( nullptr ),
  mH2ClusterYvsV( nullptr ),
  mHRatioUX( nullptr ),
  mHRatioVY( nullptr ),
  mH2ClusterXvsU_excl( nullptr ),
  mH2ClusterYvsV_excl( nullptr ),
  mHRatioUX_excl( nullptr ),
  mHRatioVY_excl( nullptr )
{
    for ( int io = 0; io < 4; io++ ) mNClusOrient[io] = 0;
    for ( int i = 0; i < 6; i++ ) mPProfileByMult[i] = nullptr;
}

StFttChargeSharingMaker::~StFttChargeSharingMaker()
{ }

Int_t StFttChargeSharingMaker::Init()
{
    mFout = new TFile( mOutFile, "RECREATE" );

    mHNStrips = new TH1F( "hNStrips", "cluster multiplicity;nStrips;clusters", 12, -0.5, 11.5 );
    mHSumAdc  = new TH1F( "hSumAdc",  "cluster total ADC;sumAdc;clusters", 200, 0, 4000 );
    mHSumAdcN1  = new TH1F( "hSumAdcN1",  "cluster total ADC, nStrips==1;sumAdc;clusters", 200, 0, 4000 );
    mHSumAdcN2p = new TH1F( "hSumAdcN2p", "cluster total ADC, nStrips>=2;sumAdc;clusters", 200, 0, 4000 );
    mHPeakAdc = new TH1F( "hPeakAdc", "cluster peak-strip ADC;peakAdc;clusters", 200, 0, 1024 );
    mHRawAdcInTime = new TH1F( "hRawAdcInTime", "in-time raw hit ADC (pre-clustering);adc;hits", 200, 0, 1024 );
    mPProfileAll = new TProfile( "pProfileAll", "mean ADC fraction vs strip offset from peak (all clusters);offset;<ADC fraction>", 9, -4.5, 4.5 );
    const char* multTag[6] = { "", "", "2", "3", "4", "5p" };
    for ( int m = 2; m <= 5; m++ ) {
        mPProfileByMult[m] = new TProfile( Form("pProfile_nStrips%s", multTag[m]),
            Form("mean ADC fraction vs offset, nStrips=%s;offset;<ADC fraction>", multTag[m]),
            9, -4.5, 4.5 );
    }

    mH2XvsU = new TH2F( "h2XvsU", "Vertical(X) vs DiagonalV(U) total ADC, per plane/quad/event;X sumAdc;U sumAdc", 200, 0, 100000, 200, 0, 100000 );
    mH2YvsV = new TH2F( "h2YvsV", "Horizontal(Y) vs DiagonalH(V) total ADC, per plane/quad/event;Y sumAdc;V sumAdc", 200, 0, 100000, 200, 0, 100000 );

    mH2ClusterXvsU = new TH2F( "h2ClusterXvsU", "cluster sumAdc, all X-U combinations same plane/quad/event;X cluster sumAdc;U cluster sumAdc", 150, 0, 3000, 150, 0, 3000 );
    mH2ClusterYvsV = new TH2F( "h2ClusterYvsV", "cluster sumAdc, all Y-V combinations same plane/quad/event;Y cluster sumAdc;V cluster sumAdc", 150, 0, 3000, 150, 0, 3000 );
    mHRatioUX = new TH1F( "hRatioUX", "cluster sumAdc ratio U/X, all combinations;U/X;combinations", 150, 0, 3 );
    mHRatioVY = new TH1F( "hRatioVY", "cluster sumAdc ratio V/Y, all combinations;V/Y;combinations", 150, 0, 3 );

    mH2ClusterXvsU_excl = new TH2F( "h2ClusterXvsU_excl", "cluster sumAdc, exactly 1 X + 1 U this plane/quad/event;X cluster sumAdc;U cluster sumAdc", 150, 0, 3000, 150, 0, 3000 );
    mH2ClusterYvsV_excl = new TH2F( "h2ClusterYvsV_excl", "cluster sumAdc, exactly 1 Y + 1 V this plane/quad/event;Y cluster sumAdc;V cluster sumAdc", 150, 0, 3000, 150, 0, 3000 );
    mHRatioUX_excl = new TH1F( "hRatioUX_excl", "cluster sumAdc ratio U/X, exactly 1 X + 1 U;U/X;pairs", 150, 0, 3 );
    mH2ClusterXvsU_cross = new TH2F( "h2ClusterXvsU_cross", "cluster sumAdc, X-U pairs whose strips CROSS;X cluster sumAdc;U cluster sumAdc", 150, 0, 3000, 150, 0, 3000 );
    mH2ClusterYvsV_cross = new TH2F( "h2ClusterYvsV_cross", "cluster sumAdc, Y-V pairs whose strips CROSS;Y cluster sumAdc;V cluster sumAdc", 150, 0, 3000, 150, 0, 3000 );
    mHRatioUX_cross = new TH1F( "hRatioUX_cross", "cluster sumAdc ratio U/X, crossing pairs;U/X;pairs", 150, 0, 3 );
    mHRatioVY_cross = new TH1F( "hRatioVY_cross", "cluster sumAdc ratio V/Y, crossing pairs;V/Y;pairs", 150, 0, 3 );
    mHRatioVY_excl = new TH1F( "hRatioVY_excl", "cluster sumAdc ratio V/Y, exactly 1 Y + 1 V;V/Y;pairs", 150, 0, 3 );

    return kStOk;
}

// Same convention as StFttHitCalibMaker (fob-based uuid), just computed
// directly from the already-mapped plane/quad/feb/vmm fields every
// StMuFttRawHit already carries -- no need for an StFttDb call at all.
int StFttChargeSharingMaker::vmmUuid( int plane, int quad, int feb, int vmm )
{
    return vmm + 10 * ( feb + 10 * ( quad + 10 * plane ) );
}

bool StFttChargeSharingMaker::readyFor( int uuid )
{
    std::map<int, std::map<short,int> >::iterator it = mDbcidHist.find( uuid );
    if ( it == mDbcidHist.end() ) return false;
    size_t n = 0;
    for ( std::map<short,int>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt ) n += jt->second;
    return n >= (size_t)mMinSamplesForReady;
}

short StFttChargeSharingMaker::anchorFor( int uuid )
{
    std::map<int, short>::iterator cached = mAnchor.find( uuid );
    std::map<int, std::map<short,int> >::iterator it = mDbcidHist.find( uuid );
    if ( it == mDbcidHist.end() || it->second.empty() ) return 0;

    // recompute the mode -- cheap, dbcidHist per VMM only ever holds as many
    // distinct keys as there are distinct dbcid values actually seen (a few
    // hundred at most), not one entry per hit
    short best = it->second.begin()->first;
    int bestCount = it->second.begin()->second;
    for ( std::map<short,int>::iterator jt = it->second.begin(); jt != it->second.end(); ++jt ) {
        if ( jt->second > bestCount ) { bestCount = jt->second; best = jt->first; }
    }
    mAnchor[uuid] = best;
    return best;
}

void StFttChargeSharingMaker::fillCalib( int uuid, short dbcid )
{
    mDbcidHist[uuid][dbcid]++;
}

long long StFttChargeSharingMaker::packGroupKey( int plane, int quad, int row, int orientation )
{
    return ( (long long)plane << 16 ) | ( (long long)quad << 12 ) | ( (long long)row << 4 ) | (long long)orientation;
}

// Re-implements StFttClusterMaker::FindClusters()/SearchClusterEdges()/
// CalculateClusterInfo() exactly (max-ADC anchor, expand while adjacent
// (strip diff<=1) AND ADC non-increasing -- GetThresholdFor() is currently
// a stub returning 0.0 in the real maker too, so any positive ADC passes),
// on a plain per-(plane,quad,row,orientation) strip vector instead of
// StFttRawHit objects (that link isn't available from MuDst -- see header).
void StFttChargeSharingMaker::clusterAndFill( std::vector<StripHit>& hits, int row, std::vector<ClusterInfo>* out )
{
    std::sort( hits.begin(), hits.end(), []( const StripHit& a, const StripHit& b ) { return a.strip < b.strip; } );

    while ( !hits.empty() ) {
        size_t anchor = 0;
        for ( size_t i = 1; i < hits.size(); i++ ) if ( hits[i].adc > hits[anchor].adc ) anchor = i;

        size_t left = anchor, right = anchor;
        while ( right + 1 < hits.size() ) {
            if ( hits[right+1].strip - hits[right].strip > 1 ) break;
            if ( hits[right+1].adc > hits[right].adc ) break;
            right++;
        }
        while ( left > 0 ) {
            if ( hits[left].strip - hits[left-1].strip > 1 ) break;
            if ( hits[left-1].adc > hits[left].adc ) break;
            left--;
        }

        int nStrips = (int)(right - left + 1);
        float sumAdc = 0;
        int peakIdx = (int)left;
        float peakAdc = -1;
        for ( size_t k = left; k <= right; k++ ) {
            sumAdc += hits[k].adc;
            if ( hits[k].adc > peakAdc ) { peakAdc = hits[k].adc; peakIdx = (int)k; }
        }
        int peakStrip = hits[peakIdx].strip;

        mHNStrips->Fill( nStrips );
        mHSumAdc->Fill( sumAdc );
        if ( nStrips == 1 ) mHSumAdcN1->Fill( sumAdc );
        else mHSumAdcN2p->Fill( sumAdc );
        mHPeakAdc->Fill( peakAdc );
        int multIdx = nStrips; if ( multIdx > 5 ) multIdx = 5;
        for ( size_t k = left; k <= right; k++ ) {
            int offset = hits[k].strip - peakStrip;
            float frac = ( sumAdc > 0 ) ? hits[k].adc / sumAdc : 0;
            mPProfileAll->Fill( offset, frac );
            if ( multIdx >= 2 && mPProfileByMult[multIdx] ) mPProfileByMult[multIdx]->Fill( offset, frac );
        }
        mNClustersFound++;
        if ( out ) { ClusterInfo ci; ci.sum = sumAdc; ci.row = row; ci.peakStrip = peakStrip; out->push_back( ci ); }

        hits.erase( hits.begin() + left, hits.begin() + right + 1 );
    }
}

Int_t StFttChargeSharingMaker::Make()
{
    mEventCount++;

    StMuDst* muDst = (StMuDst*)GetInputDS("MuDst");
    if ( !muDst ) {
        StMuDstMaker* muMaker = (StMuDstMaker*)GetMaker("MuDst");
        if ( muMaker ) muDst = muMaker->muDst();
    }
    if ( !muDst ) { LOG_WARN << "StFttChargeSharingMaker::Make() - no MuDst" << endm; return kStOk; }

    StMuFttCollection* ftt = muDst->muFttCollection();
    if ( !ftt ) { LOG_WARN << "StFttChargeSharingMaker::Make() - no muFttCollection" << endm; return kStOk; }

    int nRaw = ftt->numberOfRawHits();
    mNTotalHits += nRaw;

    // pass 1 (always): feed the online per-VMM dbcid calibration, same as
    // production's StFttHitCalibMaker -- every hit, every event, so later
    // events benefit from more statistics than production's own bootstrap
    std::vector<int>   uuidOf( nRaw );
    for ( int j = 0; j < nRaw; j++ ) {
        StMuFttRawHit* mh = ftt->getRawHit( j );
        int uuid = vmmUuid( mh->plane(), mh->quadrant(), mh->feb(), mh->vmm() );
        uuidOf[j] = uuid;
        fillCalib( uuid, (short)mh->dbcid() );
    }

    if ( mUseDb && !mFttDb ) {
        mFttDb = static_cast<StFttDb*>( GetDataSet( "fttDb" ) );
        if ( !mFttDb )
            printf( "StFttChargeSharingMaker: setUseDb(true) but no \"fttDb\" in the chain -- "
                    "add St_db_Maker and StFttDbMaker ahead of this maker. Self-calibrating.\n" );
        else
            printf( "StFttChargeSharingMaker: per-VMM anchor and window from Calibrations/ftt/fttDataWindowsB\n" );
    }

    if ( mEventCount <= mWarmupEvents ) return kStOk; // let anchors build up, don't analyze yet

    // pass 2: select in-time hits (time = dbcid - anchor(VMM) within window),
    // group by (plane,quad,row,orientation)
    std::map<long long, std::vector<StripHit> > groups;

    // Same-foil charge correlation (see header comment on mH2XvsU/mH2YvsV):
    // total in-time ADC per (plane,quad,orientation), collapsed across row,
    // this event -- Vertical+DiagonalV are the same front foil, Horizontal+
    // DiagonalH the same back foil.
    std::map<long long, double> totalByPQO;

    for ( int j = 0; j < nRaw; j++ ) {
        StMuFttRawHit* mh = ftt->getRawHit( j );
        int uuid = uuidOf[j];

        // Prefer the per-VMM anchor and per-VMM window from the calibration;
        // fall back to the self-calibrated anchor and the flat window for any
        // VMM the calibration could not anchor (anchor -1).
        int anchor = -1, lo = mTimeCutLo, hi = mTimeCutHi;
        if ( mFttDb ) {
            // StFttDb keys on the raw hit's own electronics fields, so a scratch
            // StFttRawHit carrying them is all it needs; nothing else is read.
            StFttRawHit probe;
            probe.setRaw( mh->sector(), mh->rdo(), mh->feb(), mh->vmm(),
                          mh->channel(), mh->adc(), mh->bcid(), mh->tb(), mh->dbcid() );
            Short_t dbAnchor = -1;
            if ( mFttDb->getAnchor( &probe, dbAnchor ) && dbAnchor >= 0 ) {
                int mode = 0, dbLo = mTimeCutLo, dbHi = mTimeCutHi;
                mFttDb->getTimeCut( &probe, mode, dbLo, dbHi );
                anchor = dbAnchor; lo = dbLo; hi = dbHi;
                mNFromDbAnchor++;
            }
        }
        if ( anchor < 0 ) {
            if ( !readyFor( uuid ) ) continue;
            anchor = (int)anchorFor( uuid );
            mNFallback++;
        }

        int t = (int)mh->dbcid() - anchor;
        if ( t < lo || t > hi ) continue;
        mNInTimeHits++;
        mHRawAdcInTime->Fill( mh->adc() );

        long long key = packGroupKey( mh->plane(), mh->quadrant(), mh->row(), mh->orientation() );
        StripHit sh; sh.strip = mh->strip(); sh.adc = mh->adc();
        groups[key].push_back( sh );

        long long pqoKey = ( (long long)mh->plane() * 4 + mh->quadrant() ) * 4 + mh->orientation();
        totalByPQO[pqoKey] += mh->adc();
    }

    // Per-cluster sums, collapsed across row but split by (plane,quad,
    // orientation), for the same-foil combinatoric check below -- packGroupKey()
    // is (plane<<16)|(quad<<12)|(row<<4)|orientation, so orientation/quad/plane
    // unpack straight back out of the group key.
    std::map<long long, std::vector<ClusterInfo> > clusterSumsByPQO;
    for ( std::map<long long, std::vector<StripHit> >::iterator it = groups.begin(); it != groups.end(); ++it ) {
        long long key = it->first;
        int orientation = (int)( key & 0xF );
        int quad        = (int)( ( key >> 12 ) & 0xF );
        int plane       = (int)( key >> 16 );
        long long pqoKey = ( (long long)plane * 4 + quad ) * 4 + orientation;
        int row = (int)( ( key >> 4 ) & 0xFF );
        clusterAndFill( it->second, row, &clusterSumsByPQO[pqoKey] );
    }

    for ( int p = 0; p < 4; p++ ) {
        for ( int q = 0; q < 4; q++ ) {
            long long base = ( (long long)p * 4 + q ) * 4;
            std::map<long long,double>::iterator itX = totalByPQO.find( base + kFttVertical );
            std::map<long long,double>::iterator itU = totalByPQO.find( base + kFttDiagonalV );
            std::map<long long,double>::iterator itY = totalByPQO.find( base + kFttHorizontal );
            std::map<long long,double>::iterator itV = totalByPQO.find( base + kFttDiagonalH );
            double x = ( itX != totalByPQO.end() ) ? itX->second : 0.0;
            double u = ( itU != totalByPQO.end() ) ? itU->second : 0.0;
            double y = ( itY != totalByPQO.end() ) ? itY->second : 0.0;
            double v = ( itV != totalByPQO.end() ) ? itV->second : 0.0;
            if ( x > 0 || u > 0 ) mH2XvsU->Fill( x, u );
            if ( y > 0 || v > 0 ) mH2YvsV->Fill( y, v );

            std::vector<ClusterInfo>& xList = clusterSumsByPQO[ base + kFttVertical ];
            std::vector<ClusterInfo>& uList = clusterSumsByPQO[ base + kFttDiagonalV ];
            std::vector<ClusterInfo>& yList = clusterSumsByPQO[ base + kFttHorizontal ];
            std::vector<ClusterInfo>& vList = clusterSumsByPQO[ base + kFttDiagonalH ];
            mNClusOrient[0] += (double)yList.size();   // Horizontal
            mNClusOrient[1] += (double)xList.size();   // Vertical
            mNClusOrient[2] += (double)vList.size();   // DiagonalH
            mNClusOrient[3] += (double)uList.size();   // DiagonalV
            mNPQOSamples++;
            for ( size_t ix = 0; ix < xList.size(); ix++ ) {
                for ( size_t iu = 0; iu < uList.size(); iu++ ) {
                    mH2ClusterXvsU->Fill( xList[ix].sum, uList[iu].sum );
                    if ( xList[ix].sum > 0 ) mHRatioUX->Fill( uList[iu].sum / xList[ix].sum );
                    mNPairsXU++;
                    double cx, cy;
                    if ( mFttDb && mFttDb->stripsCross( xList[ix].row, xList[ix].peakStrip, kFttVertical,
                                                        uList[iu].row, uList[iu].peakStrip, kFttDiagonalV,
                                                        cx, cy, mCrossTolerance ) ) {
                        mNPairsXUCross++;
                        mH2ClusterXvsU_cross->Fill( xList[ix].sum, uList[iu].sum );
                        if ( xList[ix].sum > 0 ) mHRatioUX_cross->Fill( uList[iu].sum / xList[ix].sum );
                    }
                }
            }
            for ( size_t iy = 0; iy < yList.size(); iy++ ) {
                for ( size_t iv = 0; iv < vList.size(); iv++ ) {
                    mH2ClusterYvsV->Fill( yList[iy].sum, vList[iv].sum );
                    if ( yList[iy].sum > 0 ) mHRatioVY->Fill( vList[iv].sum / yList[iy].sum );
                    mNPairsYV++;
                    double cx, cy;
                    if ( mFttDb && mFttDb->stripsCross( yList[iy].row, yList[iy].peakStrip, kFttHorizontal,
                                                        vList[iv].row, vList[iv].peakStrip, kFttDiagonalH,
                                                        cx, cy, mCrossTolerance ) ) {
                        mNPairsYVCross++;
                        mH2ClusterYvsV_cross->Fill( yList[iy].sum, vList[iv].sum );
                        if ( yList[iy].sum > 0 ) mHRatioVY_cross->Fill( vList[iv].sum / yList[iy].sum );
                    }
                }
            }
            double ecx, ecy;
            if ( xList.size() == 1 && uList.size() == 1 ) {
                mH2ClusterXvsU_excl->Fill( xList[0].sum, uList[0].sum );
                if ( xList[0].sum > 0 ) mHRatioUX_excl->Fill( uList[0].sum / xList[0].sum );
                mNExclXU++;
                if ( mFttDb && mFttDb->stripsCross( xList[0].row, xList[0].peakStrip, kFttVertical,
                                                    uList[0].row, uList[0].peakStrip, kFttDiagonalV,
                                                    ecx, ecy, mCrossTolerance ) ) mNExclXUCross++;
            }
            if ( yList.size() == 1 && vList.size() == 1 ) {
                mH2ClusterYvsV_excl->Fill( yList[0].sum, vList[0].sum );
                if ( yList[0].sum > 0 ) mHRatioVY_excl->Fill( vList[0].sum / yList[0].sum );
                mNExclYV++;
                if ( mFttDb && mFttDb->stripsCross( yList[0].row, yList[0].peakStrip, kFttHorizontal,
                                                    vList[0].row, vList[0].peakStrip, kFttDiagonalH,
                                                    ecx, ecy, mCrossTolerance ) ) mNExclYVCross++;
            }
        }
    }

    return kStOk;
}

Int_t StFttChargeSharingMaker::Finish()
{
    printf( "\n=== StFttChargeSharingMaker::Finish() ===\n" );
    printf( "  events processed: %d (warmup=%d, analyzed=%d)\n", mEventCount, mWarmupEvents, mEventCount - mWarmupEvents );
    printf( "  total raw hits seen: %d, in-time (kept for clustering): %d (%.1f%%)\n",
            mNTotalHits, mNInTimeHits, mNTotalHits > 0 ? 100.0 * mNInTimeHits / mNTotalHits : 0.0 );
    printf( "  clusters found: %d\n", mNClustersFound );
    if ( mNPQOSamples > 0 ) {
        if ( mFttDb )
            printf( "  time cut: per-VMM anchor+window from the DB (%d hits on a DB anchor, %d fell back)\n",
                    mNFromDbAnchor, mNFallback );
        else
            printf( "  time cut used: [%d, %d], self-calibrated anchors\n", mTimeCutLo, mTimeCutHi );
        printf( "  clusters per (plane,quad,event), by orientation -- the pairing figure of merit:\n" );
        printf( "      H %.2f   V %.2f   dH %.2f   dV %.2f   (%d samples)\n",
                mNClusOrient[0]/mNPQOSamples, mNClusOrient[1]/mNPQOSamples,
                mNClusOrient[2]/mNPQOSamples, mNClusOrient[3]/mNPQOSamples, mNPQOSamples );
        if ( mNExclXU > 0 || mNExclYV > 0 )
            printf( "  gate EFFICIENCY on the clean 1+1 pairs (must be high, or the gate is wrong): "
                    "X-U %.0f of %.0f (%.1f%%), Y-V %.0f of %.0f (%.1f%%)\n",
                    mNExclXUCross, mNExclXU, mNExclXU > 0 ? 100.0*mNExclXUCross/mNExclXU : 0.0,
                    mNExclYVCross, mNExclYV, mNExclYV > 0 ? 100.0*mNExclYVCross/mNExclYV : 0.0 );
        if ( mNPairsXU > 0 || mNPairsYV > 0 )
            printf( "  geometric crossing gate (tolerance %.1f mm): X-U %.0f of %.0f pairs cross (%.1f%%), "
                    "Y-V %.0f of %.0f (%.1f%%)\n",
                    mCrossTolerance, mNPairsXUCross, mNPairsXU,
                    mNPairsXU > 0 ? 100.0*mNPairsXUCross/mNPairsXU : 0.0,
                    mNPairsYVCross, mNPairsYV,
                    mNPairsYV > 0 ? 100.0*mNPairsYVCross/mNPairsYV : 0.0 );
        printf( "      combinations per (plane,quad,event): X-U %.1f   Y-V %.1f\n",
                (mNClusOrient[1]/mNPQOSamples)*(mNClusOrient[3]/mNPQOSamples),
                (mNClusOrient[0]/mNPQOSamples)*(mNClusOrient[2]/mNPQOSamples) );
    }
    printf( "  distinct VMMs with a computed anchor: %d\n", (int)mAnchor.size() );

    if ( mFout ) {
        mFout->cd();
        mHNStrips->Write();
        mHSumAdc->Write();
        mHSumAdcN1->Write();
        mHSumAdcN2p->Write();
        mHPeakAdc->Write();
        mHRawAdcInTime->Write();
        mH2XvsU->Write();
        mH2YvsV->Write();
        mH2ClusterXvsU->Write();
        mH2ClusterYvsV->Write();
        mHRatioUX->Write();
        mHRatioVY->Write();
        mH2ClusterXvsU_excl->Write();
        mH2ClusterYvsV_excl->Write();
        mHRatioUX_excl->Write();
        mH2ClusterXvsU_cross->Write();
        mH2ClusterYvsV_cross->Write();
        mHRatioUX_cross->Write();
        mHRatioVY_cross->Write();
        mHRatioVY_excl->Write();
        mPProfileAll->Write();
        for ( int m = 2; m <= 5; m++ ) if ( mPProfileByMult[m] ) mPProfileByMult[m]->Write();
        mFout->Close();
        mFout = nullptr;
        printf( "  wrote %s\n", mOutFile.Data() );
    }

    return kStOk;
}
