#include <cmath>
#include <map>
#include <utility>

#include "StFttSlowSimMaker.h"

#include "StEvent.h"
#include "StEnumerations.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFttRawHit.h"

#include "StFttDbMaker/StFttDb.h"

#include "tables/St_g2t_fts_hit_Table.h"

#include "TRandom.h"

// One particle crossing one gas layer: GEANT can report it as several steps
// (~6% of crossings, up to 17), which are one physical ionisation trail
// through one foil. They are one charge deposit, so the energy is summed and
// the position is their dE-weighted centroid -- the charge-weighted centre is
// what the strips actually see. Measured against taking the first step
// instead, the centroid moves the hit by less than 0.4mm (an eighth of a
// strip pitch) for 97% of crossings, but by more than a pitch for the ~1%
// that travel along the gas or throw a delta ray.
struct FttGeantCrossing {
    double deSum;      // GeV, summed over the steps of this crossing
    double gxWeighted; // sum of de*x over the steps, cm*GeV
    double gyWeighted; // sum of de*y
    int    trackId;
    int    volumeId;
};

ClassImp( StFttSlowSimMaker )

//_____________________________________________________________
StFttSlowSimMaker::StFttSlowSimMaker( const char* name )
: StMaker( name ),
  mEvent( 0 ),
  mFttCollection( 0 ),
  mFttDb( 0 )
{ }

//_____________________________________________________________
StFttSlowSimMaker::~StFttSlowSimMaker()
{ }

//_____________________________________________________________
Int_t StFttSlowSimMaker::Init()
{
    return kStOk;
}

//_____________________________________________________________
Int_t StFttSlowSimMaker::Make()
{
    mEvent = (StEvent*)GetInputDS("StEvent");
    if ( !mEvent ) {
        mEvent = new StEvent();
        AddData( mEvent );
        LOG_INFO << "StFttSlowSimMaker::Make() - Creating StEvent" << endm;
    }

    if ( mEvent->fttCollection() == nullptr ) {
        LOG_INFO << "StFttSlowSimMaker::Make() - Creating FttCollection" << endm;
        mEvent->setFttCollection( new StFttCollection() );
    }
    mFttCollection = mEvent->fttCollection();
    mFttCollection->rawHits().clear();

    mFttDb = static_cast<StFttDb*>( GetDataSet( "fttDb" ) );
    if ( !mFttDb ) {
        LOG_ERROR << "StFttSlowSimMaker::Make() - fttDb dataset not found (is StFttDbMaker in the chain, before this maker?), cannot continue" << endm;
        return kStErr;
    }

    St_g2t_fts_hit *geantFtt = (St_g2t_fts_hit*)GetDataSet( "geant/g2t_stg_hit" );
    if ( !geantFtt ) {
        LOG_WARN << "StFttSlowSimMaker::Make() - geant/g2t_stg_hit is empty" << endm;
        return kStOk;
    }

    // First pass: collapse the GEANT step list into one entry per crossing
    // (one track through one gas layer), summing the energy loss. Everything
    // below works from these, one crossing -> one synthetic cluster per
    // orientation. Front and back are different physical foils with
    // different volume_ids, so each still gets its own cluster
    // (StFttDb::reverseHardwareMap() disambiguates the two orientations at a
    // given row/strip -- see status_ftt_sim_maker.txt item 3).
    std::map< std::pair<int,int>, FttGeantCrossing > crossings;
    for ( int ig = 0; ig < geantFtt->GetNRows(); ig++ ) {
        g2t_fts_hit_st *step = (g2t_fts_hit_st*)geantFtt->At( ig );
        if ( !step ) continue;
        // operator[] value-initialises a new entry, so the three sums start
        // at zero and the first visit needs no special case.
        FttGeantCrossing &c = crossings[ std::make_pair( step->track_p, step->volume_id ) ];
        c.trackId     = step->track_p;
        c.volumeId    = step->volume_id;
        c.deSum      += step->de;
        c.gxWeighted += step->de * step->x[0];
        c.gyWeighted += step->de * step->x[1];
    }

    // Accumulate ADC per physical channel across ALL truth hits this event
    // before creating any StFttRawHit -- two nearby tracks (or one track's
    // center strip landing on another's neighbor spill) can legitimately
    // land on the same channel, and a real DAQ channel reports one summed
    // ADC value per event, not multiple independent records. See
    // status_ftt_sim_maker.txt item 9.
    std::map<long long, ChannelAccum> accum;

    int nTruthHits = 0, nBadVolumeId = 0, nNoStrip = 0;

    std::map< std::pair<int,int>, FttGeantCrossing >::const_iterator ic;
    for ( ic = crossings.begin(); ic != crossings.end(); ++ic ) {
        const FttGeantCrossing &crossing = ic->second;

        int track_id  = crossing.trackId;
        int volume_id = crossing.volumeId;

        // volume_id = 100*station + 10*pentagon + layer, where layer is 1
        // (front foil) or 2 (back foil): each station in
        // StarVMC/Geometry/StgmGeo/StgmGeo1.xml places four STFM pentagons in
        // turn inside STGM, and each pentagon holds two STMG gas layers.
        // Both station and pentagon are 0-based, from g2t_stg_volume_id():
        //     station = (numbv(1)-1)/4 ; plane = MOD(numbv(1)-1,4)
        // Careful with the names: what g2t calls "plane" is the PENTAGON, and
        // its "station" is the hundreds digit. This is the station, so call it
        // that -- StFttDb's "plane" is the same thing, the pentagon is not.
        // No -1 before the divide: 10*pentagon + layer <= 32 can never carry
        // into the hundreds, so it would change nothing for a real volume_id,
        // while hiding a bad one -- 400 (station 4) would decode as 3 and slip
        // past the range check below instead of being reported.
        int station_id = volume_id / 100;
        int agmlQuad   = (volume_id / 10) % 10;
        // front (even volume_id) = vertical/X strips, back (odd) = horizontal/Y strips --
        // same convention as StFttClusterPointMaker::MakeGeantPoints() / StFttFastSimMaker
        UChar_t intendedOrientation = (volume_id % 2 == 0) ? kFttVertical : kFttHorizontal;

        if ( volume_id < 1 ||
             station_id < 0 || station_id >= (int)StFttDb::nPlane ||
             agmlQuad   < 0 || agmlQuad   >= (int)StFttDb::nQuadPerPlane ) {
            // Not a position we can place: the sTGC geometry cannot produce
            // such a volume_id, so this means the geometry and this decoding
            // have gone out of step -- say so rather than dropping it quietly.
            if ( nBadVolumeId < 10 ) {
                LOG_WARN << "StFttSlowSimMaker::Make() - g2t_stg_hit volume_id " << volume_id
                         << " decodes to station " << station_id << ", pentagon " << agmlQuad
                         << ", outside the sTGC geometry; is StgmGeo's volume numbering still "
                         << "100*station + 10*pentagon + layer?" << endm;
            }
            nBadVolumeId++;
            continue;
        }

        // AGML places the four pentagons of a station in rotation order
        // (A, D, C, B) while StFttDb numbers quadrants A, B, C, D, so the
        // GEANT pentagon copy number and the StFttDb quadrant index differ by
        // this permutation. Measured on a particle-gun sample: of the truth
        // hits whose local position lands on strips, 100% agree with this
        // table and none contradict it (see status_ftt_sim_maker.txt).
        static const int kAgmlToDbQuadrant[4] = { 0, 3, 2, 1 };
        int quadrant_id = kAgmlToDbQuadrant[ agmlQuad ];

        // GEANT truth global position, cm: the dE-weighted centroid of this
        // crossing's steps. deSum > 0 for every crossing we built, but guard
        // anyway -- a zero-energy step list would divide by zero.
        if ( crossing.deSum <= 0 ) continue;
        double gx = crossing.gxWeighted / crossing.deSum;
        double gy = crossing.gyWeighted / crossing.deSum;

        // Invert StFttClusterPointMaker::MakeGlobalPoints()'s forward
        // transform ( global = ((local*s)+d)/10 ) with the SAME live StFttDb
        // call, so the whole chain is a self-consistent round trip through the
        // real local<->global formula (whatever it currently is -- bug or not,
        // see status_ftt_sim_maker.txt) rather than an independently
        // re-derived geometry. sx/sy are always +-1, so dividing == multiplying.
        float dx, sx, dy, sy, dz, sz;
        mFttDb->getGloablOffset_ClusterPoint( (UChar_t)station_id, (UChar_t)quadrant_id, dx, sx, dy, sy, dz, sz );
        double local_x = (gx * 10.0 - dx) * sx;
        double local_y = (gy * 10.0 - dy) * sy;

        // Diagonal strips are a SECOND strip layer etched on the SAME foil as
        // the primary (Vertical or Horizontal) one -- see StEnumerations.h's
        // "diagonal strips on the vertical/horizontal chamber" comment and
        // status_ftt_sim_maker.txt item 5 -- not a separate spatial zone, so
        // try both from the same truth hit and inject whichever is valid.
        UChar_t diagonalOrientation = ( kFttVertical == intendedOrientation ) ? kFttDiagonalV : kFttDiagonalH;
        const UChar_t attempts[2] = { intendedOrientation, diagonalOrientation };

        // Charge for this crossing, from the energy GEANT actually deposited.
        // Drawn ONCE: intendedOrientation and diagonalOrientation are the same
        // physical foil's two strip layers (Vertical+DiagonalV = front foil,
        // Horizontal+DiagonalH = back foil) reading the same charge deposit --
        // see the header comment on setAdcPerKeV().
        double totalAdc = mAdcPerKeV * ( crossing.deSum * 1e6 /* GeV -> keV */ ) * gainFluctuation();
        const double kPitch = StFttDb::stripPitch;

        bool anyBuilt = false;
        for ( int a = 0; a < 2; a++ ) {
            UChar_t orientation = attempts[a];
            int row = -1, centerStrip = -1;
            double subStripOffset = 0; // true hit position within the center strip, relative to its own center, mm
            if ( !computeRowStrip( orientation, local_x, local_y, row, centerStrip, subStripOffset ) ) continue;

            if ( mDebug ) {
                printf( "StFttSlowSimMaker DEBUG: track=%d plane=%d quad=%d orient=%d g=(%.3f,%.3f) local=(%.2f,%.2f) row=%d strip=%d off=%.2f\n",
                        track_id, station_id, quadrant_id, (int)orientation, gx, gy, local_x, local_y, row, centerStrip, subStripOffset );
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
                accumulateStrip( accum, station_id, quadrant_id, row, strip, adc, (UShort_t) track_id, orientation, ( 0 == d ) );
            }
        }
        if ( anyBuilt ) nTruthHits++;
        else             nNoStrip++;
    } // loop on crossings

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

    LOG_INFO << "StFttSlowSimMaker made " << nStripsAdded << " raw hits ("
              << accum.size() << " unique channels) from "
              << nTruthHits << " truth hits this event" << endm;

    // Not an error, and rare: GEANT records a hit anywhere in the STMG gas
    // volume, including the narrow band inside the first strip-group edge
    // (YX_StripGroupEdge[0] = 11.49mm) where the foil carries no strips yet.
    // Measured on an electron gun: 73 of 11404 crossings, 0.64%, about two
    // thirds of them in that band. (A hit is kept if EITHER the X/Y or the
    // diagonal orientation finds a strip, and the diagonal rotation covers
    // most of what the X/Y rows miss, which is why this is so small.)
    // Reported so the number is visible instead of hidden in a silent skip.
    if ( nNoStrip > 0 ) {
        LOG_INFO << "StFttSlowSimMaker - " << nNoStrip
                  << " truth hits fell outside the strip-covered area of their foil "
                  << "(no strips there), no raw hit made for them" << endm;
    }
    if ( nBadVolumeId > 0 ) {
        LOG_WARN << "StFttSlowSimMaker - " << nBadVolumeId
                  << " truth hits had a volume_id this maker could not decode" << endm;
    }

    return kStOk;
}

//_____________________________________________________________
// Gas amplification spread: a Polya/Gamma factor with mean 1 and relative
// width mGainSpread, drawn as the mean of 1/mGainSpread^2 unit exponentials
// (so the achievable widths are 1/sqrt(integer); 0.3 gives 0.302). Small for
// this detector by design: ~0.78keV in the 3mm gap at W~25-30eV is ~26-31
// primary electrons, each avalanching independently, so the spread averages
// down to sqrt(f/N) ~ 13-20%. Real data prefers ~30% -- the same ballpark,
// a little more, as TGC saturated-mode gain plus VMM shaping would give.
// mGainSpread = 0 switches it off entirely (pure scaling of de).
double StFttSlowSimMaker::gainFluctuation()
{
    if ( mGainSpread <= 0 ) return 1.0;
    int nStage = (int) ( 1.0 / ( mGainSpread * mGainSpread ) );
    if ( nStage < 1 ) nStage = 1;
    double g = 0;
    for ( int is = 0; is < nStage; is++ ) g += gRandom->Exp( 1.0 );
    return g / (double)nStage;
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
bool StFttSlowSimMaker::computeRowStrip( UChar_t orientation, double local_x, double local_y, int &row, int &centerStrip, double &subStripOffset )
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
long long StFttSlowSimMaker::packChannelKey( int plane, int quad, int row, int strip )
{
    return ( (long long)plane << 24 )
         | ( (long long)quad  << 20 )
         | ( (long long)row   << 12 )
         | ( (long long)strip );
}

//_____________________________________________________________
void StFttSlowSimMaker::accumulateStrip( std::map<long long, ChannelAccum> &accum, int plane, int quad, int row, int strip, float adc, UShort_t idTruth, UChar_t orientation, bool isCenter )
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
void StFttSlowSimMaker::addStrip( int plane, int quad, int row, int strip, float adc, UShort_t idTruth, UChar_t expectOrientation )
{
    int rob = -1, feb = -1, vmm = -1, ch = -1;
    // Pre-set orientation to the wanted value: reverseHardwareMap() now
    // disambiguates the (row,strip) collision (both an H and a V -- or
    // DiagonalH/DiagonalV -- channel exist there) using this as a hint.
    UChar_t orientation = expectOrientation;
    bool ok = mFttDb->reverseHardwareMap( rob, feb, vmm, ch, plane, quad, row, strip, orientation );
    if ( !ok ) {
        if ( mDebug ) {
            LOG_INFO << "StFttSlowSimMaker::addStrip - no hardware address for plane=" << plane
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
