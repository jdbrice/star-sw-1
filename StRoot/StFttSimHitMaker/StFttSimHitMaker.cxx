#include <cmath>
#include <map>
#include <utility>

#include "StFttSimHitMaker.h"

#include "StEvent.h"
#include "StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFttRawHit.h"

#include "StFttDbMaker/StFttDb.h"

#include "tables/St_g2t_fts_hit_Table.h"

#include "TRandom.h"

ClassImp( StFttSimHitMaker )

//_____________________________________________________________
StFttSimHitMaker::StFttSimHitMaker( const char* name )
: StMaker( name ),
  mEvent( 0 ),
  mFttCollection( 0 ),
  mFttDb( 0 )
{ }

//_____________________________________________________________
StFttSimHitMaker::~StFttSimHitMaker()
{ }

//_____________________________________________________________
Int_t StFttSimHitMaker::Init()
{
    return kStOk;
}

//_____________________________________________________________
Int_t StFttSimHitMaker::Make()
{
    mEvent = (StEvent*)GetInputDS("StEvent");
    if ( !mEvent ) {
        mEvent = new StEvent();
        AddData( mEvent );
        LOG_INFO << "StFttSimHitMaker::Make() - Creating StEvent" << endm;
    }

    if ( mEvent->fttCollection() == nullptr ) {
        LOG_INFO << "StFttSimHitMaker::Make() - Creating FttCollection" << endm;
        mEvent->setFttCollection( new StFttCollection() );
    }
    mFttCollection = mEvent->fttCollection();
    mFttCollection->rawHits().clear();

    mFttDb = static_cast<StFttDb*>( GetDataSet( "fttDb" ) );
    if ( !mFttDb ) {
        LOG_ERROR << "StFttSimHitMaker::Make() - fttDb dataset not found (is StFttDbMaker in the chain, before this maker?), cannot continue" << endm;
        return kStErr;
    }

    St_g2t_fts_hit *geantFtt = (St_g2t_fts_hit*)GetDataSet( "geant/g2t_stg_hit" );
    if ( !geantFtt ) {
        LOG_WARN << "StFttSimHitMaker::Make() - geant/g2t_stg_hit is empty" << endm;
        return kStOk;
    }

    // one truth hit -> one synthetic 3-strip cluster (per orientation, see
    // header comment); dedupe multiple GEANT steps on the same track+volume,
    // same convention as StFttClusterPointMaker::MakeGeantPoints()
    std::map<std::pair<int,int>,int> track_vol_count;

    // Accumulate ADC per physical channel across ALL truth hits this event
    // before creating any StFttRawHit -- two nearby tracks (or one track's
    // center strip landing on another's neighbor spill) can legitimately
    // land on the same channel, and a real DAQ channel reports one summed
    // ADC value per event, not multiple independent records. See
    // status_ftt_sim_maker.txt item 9.
    std::map<long long, ChannelAccum> accum;

    int nTruthHits = 0;

    for ( int i = 0; i < geantFtt->GetNRows(); i++ ) {
        g2t_fts_hit_st *git = (g2t_fts_hit_st*)geantFtt->At(i);
        if ( !git ) continue;

        int track_id  = git->track_p;
        int volume_id = git->volume_id;

        // Dedupe by (track, volume): front/back are genuinely different
        // physical foils (StFttDb::reverseHardwareMap() now disambiguates
        // both orientations at a given row/strip -- see status_ftt_sim_maker.txt
        // item 3), so each gets its own synthetic cluster again.
        if ( ++track_vol_count[ std::make_pair(track_id, volume_id) ] > 1 ) continue;

        int plane_id = (volume_id - 1) / 100;
        // front (even volume_id) = vertical/X strips, back (odd) = horizontal/Y strips --
        // same convention as StFttClusterPointMaker::MakeGeantPoints() / StFttFastSimMaker
        UChar_t intendedOrientation = (volume_id % 2 == 0) ? kFttVertical : kFttHorizontal;

        if ( plane_id < 0 || plane_id >= (int)StFttDb::nPlane ) continue;

        // GEANT truth global position, cm
        double gx = git->x[0];
        double gy = git->x[1];

        // Brute force the quadrant: try all 4, invert
        // StFttClusterPointMaker::MakeGlobalPoints()'s forward transform
        // ( global = ((local*s)+d)/10 ) with the SAME live StFttDb call for each,
        // and keep whichever gives a local (x,y) inside the physical strip range
        // (0..~560mm). This makes the whole chain a self-consistent round trip
        // through the real local<->global formula (whatever it currently is --
        // bug or not, see status_ftt_sim_maker.txt), not an independently
        // re-derived geometry. sx/sy are always +-1, so dividing == multiplying.
        const double kLocalMax = 560.0; // mm, a bit above the largest strip-group edge (548.7mm)
        int quadrant_id = -1;
        double local_x = 0, local_y = 0;
        for ( int q = 0; q < (int)StFttDb::nQuadPerPlane; q++ ) {
            float dx, sx, dy, sy, dz, sz;
            mFttDb->getGloablOffset_ClusterPoint( (UChar_t)plane_id, (UChar_t)q, dx, sx, dy, sy, dz, sz );
            double lx = (gx * 10.0 - dx) * sx;
            double ly = (gy * 10.0 - dy) * sy;
            if ( lx >= 0 && lx <= kLocalMax && ly >= 0 && ly <= kLocalMax ) {
                quadrant_id = q;
                local_x = lx;
                local_y = ly;
                break;
            }
        }
        if ( quadrant_id < 0 ) continue; // didn't land cleanly in any quadrant's valid local range

        // Diagonal strips are a SECOND strip layer etched on the SAME foil as
        // the primary (Vertical or Horizontal) one -- see StEnumerations.h's
        // "diagonal strips on the vertical/horizontal chamber" comment and
        // status_ftt_sim_maker.txt item 5 -- not a separate spatial zone, so
        // try both from the same truth hit and inject whichever is valid.
        UChar_t diagonalOrientation = ( kFttVertical == intendedOrientation ) ? kFttDiagonalV : kFttDiagonalH;
        const UChar_t attempts[2] = { intendedOrientation, diagonalOrientation };

        // Drawn ONCE per truth hit: intendedOrientation and diagonalOrientation
        // are the same physical foil's two strip layers (Vertical+DiagonalV =
        // front foil, Horizontal+DiagonalH = back foil), reading the same
        // charge deposit -- see header comment on setMeanTotalAdc().
        double totalAdc = -mMeanTotalAdc * std::log( gRandom->Rndm() );
        const double kPitch = StFttDb::stripPitch;

        bool anyBuilt = false;
        for ( int a = 0; a < 2; a++ ) {
            UChar_t orientation = attempts[a];
            int row = -1, centerStrip = -1;
            double subStripOffset = 0; // true hit position within the center strip, relative to its own center, mm
            if ( !computeRowStrip( orientation, local_x, local_y, row, centerStrip, subStripOffset ) ) continue;

            if ( mDebug ) {
                printf( "StFttSimHitMaker DEBUG: track=%d plane=%d quad=%d orient=%d g=(%.3f,%.3f) local=(%.2f,%.2f) row=%d strip=%d off=%.2f\n",
                        track_id, plane_id, quadrant_id, (int)orientation, gx, gy, local_x, local_y, row, centerStrip, subStripOffset );
            }

            anyBuilt = true;

            // Real-data-based charge sharing (see header comment): split the
            // shared total charge over center strip +-1 via a Gaussian(sigma)
            // response integrated over each strip's bounds relative to the
            // true sub-strip position -- an off-center hit naturally gets an
            // asymmetric split, unlike a fixed ratio. Threshold is applied
            // once, at final emission below, after summing all contributions
            // per channel (center strip is always kept -- see header comment).
            for ( int d = -1; d <= 1; d++ ) {
                int strip = centerStrip + d;
                if ( strip < 0 ) continue;
                double lo = ( d - 0.5 ) * kPitch - subStripOffset;
                double hi = ( d + 0.5 ) * kPitch - subStripOffset;
                double frac = 0.5 * ( std::erf( hi / ( mChargeShareSigma * std::sqrt(2.0) ) )
                                     - std::erf( lo / ( mChargeShareSigma * std::sqrt(2.0) ) ) );
                float adc = (float)( totalAdc * frac );
                accumulateStrip( accum, plane_id, quadrant_id, row, strip, adc, (UShort_t) track_id, orientation, ( 0 == d ) );
            }
        }
        if ( anyBuilt ) nTruthHits++;
    } // loop on geant hits

    // Now that every truth hit's contribution to every channel has been
    // summed (a channel near two nearby tracks' edges can clear threshold
    // only once both contributions are added -- real hardware behavior),
    // emit one StFttRawHit per channel that clears the hardware-like
    // zero-suppression threshold -- EXCEPT a center (d==0) strip, which is
    // always emitted (floored up to threshold if the draw came in low): see
    // header comment on setStripAdcThreshold() for why (keeps efficiency at
    // 100%, as this project's TODO deliberately still wants for now).
    const float kMaxAdc = 1023.0f; // 10-bit ADC, matches StFttDb::maxADC-2
    int nStripsAdded = 0;
    for ( std::map<long long, ChannelAccum>::iterator it = accum.begin(); it != accum.end(); ++it ) {
        ChannelAccum &c = it->second;
        float adc = c.adcSum;
        if ( adc < mStripAdcThreshold ) {
            if ( !c.isCenter ) continue;
            adc = mStripAdcThreshold;
        }
        if ( adc > kMaxAdc ) adc = kMaxAdc; // saturate, matches real hardware
        size_t before = mFttCollection->rawHits().size();
        addStrip( c.plane, c.quad, c.row, c.strip, adc, c.idTruth, c.orientation );
        if ( mFttCollection->rawHits().size() > before ) nStripsAdded++;
    }

    LOG_INFO << "StFttSimHitMaker made " << nStripsAdded << " raw hits ("
              << accum.size() << " unique channels) from "
              << nTruthHits << " truth hits this event" << endm;

    return kStOk;
}

//_____________________________________________________________
// Invert the local-coordinate formula for the given orientation to get a
// (row, strip) address, self-consistently with whatever the real code
// computes forward:
//   - kFttVertical/kFttHorizontal (rows 0-2): StFttClusterMaker::
//     CalculateClusterInfo()'s x = strip*stripPitch - stripPitch/2 gives the
//     measured axis; StFttDb::YX_StripGroupEdge bands the row axis into rows
//     0/1/2.
//   - kFttDiagonalV/kFttDiagonalH (rows 3-4): StFttClusterPointMaker::
//     MakeLocalPoints()'s rotation x=a+(root2/2)(x'-y'), y=a+(root2/2)(x'+y')
//     (a = StFttDb::D_StripGroupEdge[0]) inverts to x'=(x+y-2a)/root2,
//     y'=(y-x)/root2; row 3 vs 4 is just the sign of y' (MakeLocalPoints
//     negates y' for row==3, both DiagonalV and DiagonalH); x' feeds the same
//     strip*pitch-pitch/2 inversion as X/Y.
// Returns false if the position falls outside this orientation's coverage
// (row 0-2 case) -- strip-range validity for either case is enforced later
// by reverseHardwareMap() itself.
bool StFttSimHitMaker::computeRowStrip( UChar_t orientation, double local_x, double local_y, int &row, int &centerStrip, double &subStripOffset )
{
    if ( kFttVertical == orientation || kFttHorizontal == orientation ) {
        double measuredAxis = ( kFttVertical == orientation ) ? local_x : local_y;
        double rowAxis      = ( kFttVertical == orientation ) ? local_y : local_x;

        row = -1;
        if      ( rowAxis >= StFttDb::YX_StripGroupEdge[2] ) row = 2;
        else if ( rowAxis >= StFttDb::YX_StripGroupEdge[1] ) row = 1;
        else if ( rowAxis >= StFttDb::YX_StripGroupEdge[0] ) row = 0;
        if ( row < 0 ) return false; // falls inside the inner (beampipe-side) cutout

        centerStrip = (int) std::lround( (measuredAxis + StFttDb::stripPitch / 2.0) / StFttDb::stripPitch );
        // true hit position relative to centerStrip's own center (see header
        // comment on setChargeShareSigma): centerStrip's center sits at
        // (centerStrip-0.5)*stripPitch in this same measuredAxis coordinate.
        subStripOffset = measuredAxis - ( centerStrip - 0.5 ) * StFttDb::stripPitch;
        return centerStrip >= 0;
    }

    if ( kFttDiagonalV == orientation || kFttDiagonalH == orientation ) {
        double a = StFttDb::D_StripGroupEdge[0];
        double x_prime = (local_x + local_y - 2.0 * a) / std::sqrt(2.0);
        double y_prime = (local_y - local_x) / std::sqrt(2.0);

        // MakeLocalPoints() uses OPPOSITE row<->sign conventions for the two
        // diagonal orientations: kFttDiagonalV negates y_prime on row==3 (so
        // y_prime>=0 means row==4); kFttDiagonalH negates on row==4 (so
        // y_prime>=0 means row==3). Getting this backwards still finds A
        // valid (row,strip) via reverseHardwareMap() (row 3 and row 4 are
        // both real, distinct hardware channels), just the WRONG one.
        if ( kFttDiagonalV == orientation ) row = ( y_prime >= 0 ) ? 4 : 3;
        else                                row = ( y_prime >= 0 ) ? 3 : 4;

        centerStrip = (int) std::lround( (x_prime + StFttDb::stripPitch / 2.0) / StFttDb::stripPitch );
        subStripOffset = x_prime - ( centerStrip - 0.5 ) * StFttDb::stripPitch;
        return centerStrip >= 0;
    }

    return false;
}

//_____________________________________________________________
// (plane,quad,row,strip) already uniquely identifies one real hardware
// channel -- XY rows (0-2) and diagonal rows (3-4) are disjoint ranges, so
// orientation doesn't need to be part of the key. Plain bit-packed (not
// hashed) so it's trivially injective: plane needs 2 bits, quad 2, row 3,
// strip 8 -- generous shifts below, no overlap.
long long StFttSimHitMaker::packChannelKey( int plane, int quad, int row, int strip )
{
    return ( (long long)plane << 24 )
         | ( (long long)quad  << 20 )
         | ( (long long)row   << 12 )
         | ( (long long)strip );
}

//_____________________________________________________________
void StFttSimHitMaker::accumulateStrip( std::map<long long, ChannelAccum> &accum, int plane, int quad, int row, int strip, float adc, UShort_t idTruth, UChar_t orientation, bool isCenter )
{
    long long key = packChannelKey( plane, quad, row, strip );
    ChannelAccum &c = accum[key]; // default-constructs on first touch of this channel
    c.plane = plane; c.quad = quad; c.row = row; c.strip = strip; c.orientation = orientation;
    c.adcSum += adc;
    if ( isCenter ) c.isCenter = true; // sticky: stays true even if a later, non-center contribution touches this channel too
    if ( adc > c.maxAdc ) { // dominant contributor's truth id represents the merged hit
        c.maxAdc = adc;
        c.idTruth = idTruth;
    }
}

//_____________________________________________________________
void StFttSimHitMaker::addStrip( int plane, int quad, int row, int strip, float adc, UShort_t idTruth, UChar_t expectOrientation )
{
    int rob = -1, feb = -1, vmm = -1, ch = -1;
    // Pre-set orientation to the wanted value: reverseHardwareMap() now
    // disambiguates the (row,strip) collision (both an H and a V -- or
    // DiagonalH/DiagonalV -- channel exist there) using this as a hint.
    UChar_t orientation = expectOrientation;
    bool ok = mFttDb->reverseHardwareMap( rob, feb, vmm, ch, plane, quad, row, strip, orientation );
    if ( !ok ) {
        if ( mDebug ) {
            LOG_INFO << "StFttSimHitMaker::addStrip - no hardware address for plane=" << plane
                      << " quad=" << quad << " row=" << row << " strip=" << strip
                      << " orientation=" << (int)expectOrientation << ", dropping" << endm;
        }
        return;
    }

    // reverseHardwareMap() returns feb/vmm in the table's 1-based convention;
    // StFttRawHit stores them 0-based (StFttDb::hardwareMap() adds the +1 back
    // when packing the key from a raw hit -- see StFttDb.cxx)
    StFttRawHit *hit = new StFttRawHit( (UChar_t)(plane + 1), (UChar_t)(quad + 1),
                                         (UChar_t)(feb - 1), (UChar_t)(vmm - 1), (UChar_t)ch,
                                         (UShort_t)adc, /*bcid*/ 0, /*tb*/ 0, /*bcidDelta*/ 0 );
    hit->setIdTruth( idTruth );
    mFttCollection->addRawHit( hit );
}
