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
  mNClustersFound( 0 ),
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
void StFttChargeSharingMaker::clusterAndFill( std::vector<StripHit>& hits, std::vector<float>* outSums )
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
        if ( outSums ) outSums->push_back( sumAdc );

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
        if ( !readyFor( uuid ) ) continue;

        int t = (int)mh->dbcid() - (int)anchorFor( uuid );
        if ( t < mTimeCutLo || t > mTimeCutHi ) continue;
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
    std::map<long long, std::vector<float> > clusterSumsByPQO;
    for ( std::map<long long, std::vector<StripHit> >::iterator it = groups.begin(); it != groups.end(); ++it ) {
        long long key = it->first;
        int orientation = (int)( key & 0xF );
        int quad        = (int)( ( key >> 12 ) & 0xF );
        int plane       = (int)( key >> 16 );
        long long pqoKey = ( (long long)plane * 4 + quad ) * 4 + orientation;
        clusterAndFill( it->second, &clusterSumsByPQO[pqoKey] );
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

            std::vector<float>& xList = clusterSumsByPQO[ base + kFttVertical ];
            std::vector<float>& uList = clusterSumsByPQO[ base + kFttDiagonalV ];
            std::vector<float>& yList = clusterSumsByPQO[ base + kFttHorizontal ];
            std::vector<float>& vList = clusterSumsByPQO[ base + kFttDiagonalH ];
            for ( size_t ix = 0; ix < xList.size(); ix++ ) {
                for ( size_t iu = 0; iu < uList.size(); iu++ ) {
                    mH2ClusterXvsU->Fill( xList[ix], uList[iu] );
                    if ( xList[ix] > 0 ) mHRatioUX->Fill( uList[iu] / xList[ix] );
                }
            }
            for ( size_t iy = 0; iy < yList.size(); iy++ ) {
                for ( size_t iv = 0; iv < vList.size(); iv++ ) {
                    mH2ClusterYvsV->Fill( yList[iy], vList[iv] );
                    if ( yList[iy] > 0 ) mHRatioVY->Fill( vList[iv] / yList[iy] );
                }
            }
            if ( xList.size() == 1 && uList.size() == 1 ) {
                mH2ClusterXvsU_excl->Fill( xList[0], uList[0] );
                if ( xList[0] > 0 ) mHRatioUX_excl->Fill( uList[0] / xList[0] );
            }
            if ( yList.size() == 1 && vList.size() == 1 ) {
                mH2ClusterYvsV_excl->Fill( yList[0], vList[0] );
                if ( yList[0] > 0 ) mHRatioVY_excl->Fill( vList[0] / yList[0] );
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
        mHRatioVY_excl->Write();
        mPProfileAll->Write();
        for ( int m = 2; m <= 5; m++ ) if ( mPProfileByMult[m] ) mPProfileByMult[m]->Write();
        mFout->Close();
        mFout = nullptr;
        printf( "  wrote %s\n", mOutFile.Data() );
    }

    return kStOk;
}
