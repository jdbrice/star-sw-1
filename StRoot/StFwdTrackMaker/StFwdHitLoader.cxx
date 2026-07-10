#include "StFwdTrackMaker/StFwdTrackMaker.h"
#include "StFwdTrackMaker/StFwdHitLoader.h"

// For StEvent
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StFstHitCollection.h"
#include "StEvent/StFstHit.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StFcsHit.h"
#include "StEpdUtil/StEpdGeom.h"
#include "StFcsDbMaker/StFcsDb.h"

// For GEANT
#include "include/Tracker/FwdHit.h"
#include "tables/St_g2t_fts_hit_Table.h"

// For MuDst
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFstCollection.h"
#include "StMuDSTMaker/COMMON/StMuFstHit.h"

// See FstWedgeAligner in StFwdHitLoader.h for what this is and where the
// numbers came from (script/computeWedgeCorrection.C). Recomputed 2026-07-07
// after fixing a +15 deg wedge-phase bug in wedgePhiBin() -- these values
// supersede the ones this array held before that fix.
const double FstWedgeAligner::kWedgePhiCorrection[3][12] = {
    {-0.007271, +0.001096, -0.001476, +0.002857, -0.001234, +0.006262, -0.005856, +0.001923, -0.004920, +0.000666, -0.004712, +0.002365}, // disk 0
    {+0.002934, -0.005147, +0.004865, -0.002575, +0.003776, -0.002012, +0.002006, -0.002858, +0.001207, -0.001799, +0.001198, -0.003505}, // disk 1
    {+0.001411, +0.004961, -0.003647, +0.002564, -0.003000, -0.001889, +0.001911, +0.004294, +0.001591, +0.001295, +0.000639, +0.002373}, // disk 2
};

//#define LOG_DEBUG if(true) std::cerr
//#define LOG_INFO if(true) std::cerr

// Issue #5: default rSize was 3.0, not kFstStripPitchR (2.875 cm) -- the actual
// FST strip pitch in r. phiSize stays numeric (= kFstStripPitchPhi).
TMatrixDSym makeFstCovMat(TVector3 hit, float rSize = kFstStripPitchR, float phiSize = 0.0040906154) {
    // we can calculate the CovMat since we know the det info, but in future we should probably keep this info in the hit itself
    // measurements on a plane only need 2x2
    // for Si geom we need to convert from cylindrical to cartesian coords
    TMatrixDSym cm(2);
    TMatrixD T(2, 2);
    TMatrixD J(2, 2);
    const float x = hit.X();
    const float y = hit.Y();
    const float R = sqrt(x * x + y * y);
    const float cosphi = x / R;
    const float sinphi = y / R;
    const float sqrt12 = sqrt(12.);

    //AAA fix: σ_total = sqrt( (pitch/sqrt12)² + δ_align² )
    //         δ_align_r=1.5 cm, δ_align_phi=0.00409 rad absorb FST alignment offsets
    //         in quadrature with the intrinsic strip resolution (pitch/sqrt12).
    //         Previously (pitch only):  dr=0.83 cm  dphi=0.00118 rad
    //         Now (pitch + alignment):  dr=1.71 cm  dphi=0.00426 rad
    //const float dr = rSize;                     // original: full pitch (no /sqrt12)
    //const float dphi = phiSize;
    //const float dr = rSize / sqrt12;            // fix A: pitch/sqrt12 only
    //const float dphi = phiSize / sqrt12;
    const float delta_align_r   = 1.5f;           //AAA alignment systematic in r (cm)
    const float delta_align_phi = 0.00409f;        //AAA alignment systematic in phi (rad) ≈ kFstStripPitchPhi
    const float dr   = rSize;    // full pitch — keeps hit search window wide in r
    const float dphi = phiSize;  // full pitch — phi pitch/sqrt12 would be correct for weighting
                                  // but tightening phi shrinks the Kalman search window at subsequent
                                  // FST disks, rejecting valid hits (same problem as r/sqrt12).
                                  // FIX NEEDED: separate hit search window from measurement sigma
                                  // in GenFit/FwdTracker — use wide search + tight measurement weight.

    // Setup the Transposed and normal Jacobian transform matrix;
    // note, the si fast sim did this wrong
    // row col
    T(0, 0) = cosphi;
    T(0, 1) = -R * sinphi;
    T(1, 0) = sinphi;
    T(1, 1) = R * cosphi;

    J(0, 0) = cosphi;
    J(0, 1) = sinphi;
    J(1, 0) = -R * sinphi;
    J(1, 1) = R * cosphi;

    TMatrixD cmcyl(2, 2);
    cmcyl(0, 0) = dr * dr;
    cmcyl(1, 1) = dphi * dphi;

    TMatrixD r = T * cmcyl * J;

    // note: float sigmaX = sqrt(r(0, 0));
    // note: float sigmaY = sqrt(r(1, 1));

    cm(0, 0) = r(0, 0);
    cm(1, 1) = r(1, 1);
    cm(0, 1) = r(0, 1);
    cm(1, 0) = r(1, 0);

    TMatrixDSym tamvoc(3);
    tamvoc( 0, 0 ) = cm(0, 0); tamvoc( 0, 1 ) = cm(0, 1); tamvoc( 0, 2 ) = 0.0;
    tamvoc( 1, 0 ) = cm(1, 0); tamvoc( 1, 1 ) = cm(1, 1); tamvoc( 1, 2 ) = 0.0;
    tamvoc( 2, 0 ) = 0.0;      tamvoc( 2, 1 ) = 0.0; tamvoc( 2, 2 )      = 0.1*0.1;


    return tamvoc;
}

int StFwdHitLoader::loadFttHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
    // cout << "log4cxx NDC depth: " << log4cxx::NDC::getDepth() << endl;
    if ( mFttDataSource == StFwdHitLoader::DataSource::IGNORE ){
        // this is a warning because it should not be set during production
        // but it is useful for testing
      if ( kLogLevel >= kLogVerbose ) LOG_WARN << "FTT Data Source is set to IGNORE, not loading FTT hits" << endm;
        return 0;
    }
    if ( kLogLevel >= kLogVerbose ) {
        LOG_DEBUG << "Loading FTT Hits" << endm;
        LOG_DEBUG << "Ftt From Source: " << mFttDataSource << endm;
    }
    if ( StFwdHitLoader::DataSource::GEANT == mFttDataSource ){
        if ( kLogLevel >= kLogVerbose ) {LOG_DEBUG << "Loading FTT Hits from GEANT" << endm;}
        return loadFttPointsFromGEANT( mcTrackMap, hitMap );
    } else if ( StFwdHitLoader::DataSource::STEVENT == mFttDataSource ){
        return loadFttPointsFromStEvent( mcTrackMap, hitMap );
    } else if ( StFwdHitLoader::DataSource::MUDST == mFttDataSource ){
        if ( kLogLevel >= kLogVerbose ) {LOG_DEBUG << "Loading FTT Hits from MuDst (Not Yet Implemented)" << endm;}
        // loadFttHitsFromMuDst( mcTrackMap, hitMap, count );
        return 0;
    }
    return 0;
} // loadFttHits

int StFwdHitLoader::loadFttPointsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap ){
    if ( kLogLevel >= kLogVerbose ) {LOG_DEBUG << "Loading FTT Hits from StEvent" << endm;}
    if ( !mStEvent ) {
        LOG_ERROR << "StEvent is NULL, cannot load Ftt hits" << endm;
        return 0;
    }
    StFttCollection *col = mStEvent->fttCollection();
    size_t numFwdHitsPrior = mFwdHitsFtt.size();
    int count = 0;
    if ( col && col->numberOfPoints() > 0 ){
        if ( kLogLevel >= kLogVerbose ) {LOG_DEBUG << "The Ftt Collection has " << col->numberOfPoints() << " points" << endm;}
        TMatrixDSym hitCov3(3);
        const double sigXY = 0.2; //
        hitCov3(0, 0) = sigXY * sigXY;
        hitCov3(1, 1) = sigXY * sigXY;
        hitCov3(2, 2) = 4; // unused since they are loaded as points on plane
        for ( auto point : col->points() ){

            float xcm = point->xyz().x();
            float ycm = point->xyz().y();
            float zcm = point->xyz().z();

            hitCov3(0, 0) = (double)point->cov()[0][0];
            hitCov3(0, 1) = (double)point->cov()[0][1];
            hitCov3(1, 0) = (double)point->cov()[1][0];
            hitCov3(1, 1) = (double)point->cov()[1][1];

            if ( kLogLevel >= kLogVerbose ){
	      LOG_INFO << "Loaded FTT point with local x: " << xcm << " y: " << ycm << " z: " << zcm << ", disk = " << ((int)point->plane()) << endm;
	      stringstream sstr;
	      sstr << "\t dx = " << sqrt(hitCov3(0, 0)) << " dy = " << sqrt(hitCov3(1, 1));
	      if ( sqrt(hitCov3(0, 0)) > sqrt(hitCov3(1, 1)) ){
                sstr << " (horizontal strip)";
	      } else if ( sqrt(hitCov3(1, 1)) > sqrt(hitCov3(0, 0)) ){
                sstr << " (vertical strip)";
	      } else {
                sstr << " (unknown orientation)";
	      }
	      sstr << endm;
	      LOG_INFO << sstr.str() << endm;
	    }

            if ( sqrt(hitCov3(0, 0)) != sqrt(hitCov3(0, 0)) ){
                if ( kLogLevel >= kLogVerbose ) {LOG_INFO << "Covariance matrix has NaN entries, skipping this hit" << endm;}
                continue;
            }
            if ( sqrt(hitCov3(1, 1)) != sqrt(hitCov3(1, 1)) ){
                if ( kLogLevel >= kLogVerbose ) {LOG_INFO << "Covariance matrix has NaN entries, skipping this hit" << endm;}
                continue;
            }

            // for now we skip diagonal strips
            // We need to develop the matching algorithm in FwdTracker to be able to use these hits, 
            // but for now we want to make sure they are not causing issues in the fitter
            if ( sqrt(hitCov3(0, 0)) == sqrt(hitCov3(1, 1)) ){
	      if ( kLogLevel >= kLogVerbose )
                LOG_INFO << "Skipping FTT point with equal covariance (diagonal/combined): "
                         << "dx=" << sqrt(hitCov3(0, 0)) << " dy=" << sqrt(hitCov3(1, 1))
                         << " disk=" << ((int)point->plane()) << endm;
                continue;
            }


            // get the track id
            int track_id = point->idTruth();
            shared_ptr<McTrack> mcTrack = nullptr;
            if ( mcTrackMap.count(track_id) ) {
                mcTrack = mcTrackMap[track_id];
                if ( kLogLevel >= kLogVerbose ) {LOG_DEBUG << "Adding McTrack to FTT hit: " << track_id << endm;}
            }

            mFwdHitsFtt.push_back(FwdHit(count++, // id
                xcm, ycm, zcm,
                -point->plane(), // volume id
                kFttId, // detid
                track_id, // track id
                hitCov3, // covariance matrix
                mcTrack) // mcTrack
                );
            mSpacepointsFtt.push_back( TVector3( xcm, ycm, zcm)  );
        } // end of loop over points
    } else {
        if ( kLogLevel >= kLogVerbose ) {LOG_DEBUG << "The Ftt Collection is EMPTY points" << endm;}
    }

    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFtt.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFtt[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);
        // add to MC track map
        if ( hit->getMcTrack() ){
            hit->getMcTrack()->addFttHit(hit);
        }
    }

    if ( numFwdHitsPost != numFwdHitsPrior ){
        if ( kLogLevel >= kLogInfo ){LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FTT hits from StEvent" << endm;}
    }
    return count;
}

int StFwdHitLoader::loadFttPointsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
    /************************************************************/
    // STGC Hits
    // St_g2t_fts_hit *mGeantFtt = (St_g2t_fts_hit *)GetDataSet("geant/g2t_stg_hit");
    if ( kLogLevel >= kLogInfo ){LOG_INFO << "Loading FTT Hits from GEANT" << endm;}
    size_t numFwdHitsPrior = mFwdHitsFtt.size();
    if (!mGeantFtt){
      if ( kLogLevel >= kLogVerbose ) LOG_WARN << "geant/g2t_stg_hit is empty" << endm;
        return 0;
    }
    int count = 0;
    // make the Covariance Matrix once and then reuse
    TMatrixDSym hitCov3(3);
    const double sigXY = 0.1;
    hitCov3(0, 0) = sigXY * sigXY;
    hitCov3(1, 1) = sigXY * sigXY;
    hitCov3(2, 2) = 0.1; // unused since they are loaded as points on plane

    int nstg = mGeantFtt->GetNRows();
    if ( kLogLevel >= kLogVerbose ) {LOG_DEBUG << "This event has " << nstg << " stg hits in geant/g2t_stg_hit " << endm;}

    for (int i = 0; i < nstg; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)mGeantFtt->At(i);
        if (!git)
            continue; // invalid geant hit
        int track_id = git->track_p;
        int volume_id = git->volume_id;
        int plane_id = (volume_id - 1) / 100;           // from 1 - 16. four chambers per station

        // only use the hits on the front modules
        if ( volume_id % 2 ==0 )
            continue;

        float x = git->x[0];// + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float y = git->x[1];// + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float z = git->x[2];

        if (plane_id < 0 || plane_id >= 4) {
            continue;
        }
        
        std::shared_ptr<McTrack> mcTrack = nullptr;
        if ( mcTrackMap.count(track_id) ) {
            mcTrack = mcTrackMap[track_id];
            LOG_DEBUG << "Adding McTrack to FTT hit: " << track_id << endm;
        } 
        mFwdHitsFtt.push_back(
            FwdHit(
                count++, // id
                x, y, z, // position
                -plane_id, // volume id
                kFttId, // detid
                track_id, // track id
                hitCov3, // covariance matrix
                mcTrack
                )
            );
        mSpacepointsFtt.push_back( TVector3( x, y, z )  );
    } // loop on hits

    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFtt.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFtt[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);

        if ( dynamic_cast<FwdHit*>(hit)->_mcTrack ){
            dynamic_cast<FwdHit*>(hit)->_mcTrack->addFttHit(hit);
        }
    }

    if ( numFwdHitsPost != numFwdHitsPrior ){
        if ( kLogLevel >= kLogInfo ){LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FTT hits from GEANT" << endm;}
    }
    return count;
} // loadFttHits

/**
 * @brief Loads FST hits from various sources into the hitmap and McTrackMap (if availabale)
 *
 * Order of precedence:
 * MuDst StMuFstCollection (Data)
 * StEvent StFstHitCollection (Data or slowsim)
 * StEvent StRndHitCollection (fast sim)
 * GEANT St_g2t_fts_hit (starsim only) - note if rasterizer is active this takes priority over FastSim
 *
 * @param mcTrackMap : MC track map if running sim
 * @param hitMap : FST hitmap to populate
 * @param count  : number of hits loaded
 */
int StFwdHitLoader::loadFstHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap ){
    if ( mFstDataSource == StFwdHitLoader::DataSource::IGNORE ){
        // this is a warning because it should not be set during production
        // but it is useful for testing
      if ( kLogLevel >= kLogVerbose ) LOG_WARN << "FST Data Source is set to IGNORE, not loading FST hits" << endm;
        return 0;
    }
    if ( kLogLevel >= kLogVerbose ){
      if ( kLogLevel >= kLogVerbose ) LOG_DEBUG << "Loading FST Hits" << endm;
      if ( kLogLevel >= kLogVerbose ) LOG_DEBUG << "Fst From Source: " << mFstDataSource << endm;
    }
    if ( StFwdHitLoader::DataSource::GEANT == mFstDataSource ){
        if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "Loading FST Hits from GEANT" << endm;}
        return loadFstHitsFromGEANT( mcTrackMap, hitMap );
    } else if ( StFwdHitLoader::DataSource::STEVENT == mFstDataSource ){
        if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "Loading FST Hits from StEvent" << endm;}
        return loadFstHitsFromStEvent( mcTrackMap, hitMap );
    } else if ( StFwdHitLoader::DataSource::MUDST == mFstDataSource ){
        if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "Loading FST Hits from MuDst" << endm;}
        return loadFstHitsFromMuDst( mcTrackMap, hitMap );
    }
    return 0;
} // loadFstHits

int StFwdHitLoader::loadFstHitsFromMuDst( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
    int count = 0;

    if ( kLogLevel >= kLogVerbose ) {
      if ( kLogLevel >= kLogVerbose ) LOG_INFO << "Loading FST hits from MuDst" << endm;
    }
    if(!mMuDstMaker) {
      if ( kLogLevel >= kLogVerbose ) LOG_WARN << " No MuDstMaker ... bye-bye" << endm;
        return 0;
    }
    StMuDst *mMuDst = mMuDstMaker->muDst();
    if(!mMuDst) {
      if ( kLogLevel >= kLogVerbose ) LOG_WARN << " No MuDst ... bye-bye" << endm;
        return 0;
    }

    StMuFstCollection * fst = mMuDst->muFstCollection();
    if (!fst) {
      if ( kLogLevel >= kLogVerbose ) LOG_WARN << "No StMuFstCollection ... bye-bye" << endm;
        return 0;
    }

    size_t numFwdHitsPrior = mFwdHitsFst.size();
    if ( kLogLevel >= kLogInfo ){LOG_INFO << "Loading " << fst->numberOfHits() << " StMuFstHits" << endm;}
    TMatrixDSym hitCov3(3);
    for ( unsigned int index = 0; index < fst->numberOfHits(); index++){
        StMuFstHit * muFstHit = fst->getHit( index );

        float vR = muFstHit->localPosition(0);
        float vPhi = muFstHit->localPosition(1);
        float vZ = muFstHit->localPosition(2);

        int diskIndex   = muFstHit->getDisk() - 1;                                    // getDisk() is 1-indexed; convert to 0-indexed
        int wedgeIndex  = (muFstHit->getWedge() - 1) % kFstNumWedgePerDisk;           // getWedge() is global 1-indexed (1-36); convert to per-disk 0-indexed (0-11)
        int sensorIndex = muFstHit->getSensor();                                       // getSensor() is already 0-indexed
        int globalIndex = FwdHit::fstGlobalSensorIndex( diskIndex, wedgeIndex, sensorIndex );
        if ( kLogLevel >= kLogVerbose ) {LOG_INFO << "wedgeIndex: " << wedgeIndex << ", sensorIndex: " << sensorIndex << ", diskIndex: " << diskIndex << ", globalIndex: " << globalIndex << endm;}

        // See bugreport_StFstHitMaker.txt / FstWedgeAligner (StFwdHitLoader.h):
        // the official reconstruction never applies the real per-sensor DB
        // alignment, leaving a measured ~30-degree-periodic (per-wedge) phi
        // offset in every FST hit. Corrected here, at MuDst-load time, rather
        // than at the source (raw DAQ reprocessing is expensive and this is
        // downstream of production anyway) -- off by default.
        if ( mApplyFstWedgeAlignment ) {
            vPhi = FstWedgeAligner::correctPhi( diskIndex, vPhi );
        }

        float x0 = vR * cos( vPhi );
        float y0 = vR * sin( vPhi );
        hitCov3 = makeFstCovMat( TVector3( x0, y0, vZ ) );
        mSpacepointsFst.push_back( TVector3( x0, y0, vZ)  );
        if ( kLogLevel >= kLogVerbose ) {LOG_INFO << TString::Format("FST local position: %f %f %f, global position: %f %f %f", vR, vPhi, vZ, x0, y0, vZ) << endm;}

        // we use d+4 so that both FTT and FST start at 4
        mFwdHitsFst.push_back(
            FwdHit(
                count++, // id
                x0, y0, vZ, // position
                diskIndex+4, // volume id
                kFstId, // detid
                0, // track id
                hitCov3, // covariance matrix
                nullptr // mcTrack
            )
        );
        mFwdHitsFst.back()._genfit_plane_index = globalIndex;
        // Store strip-native local position (no global azimuthal rotation).
        {
            int rStrip = (int)muFstHit->getMeanRStrip();
            mFwdHitsFst.back()._localPosition[0] = kFstrStart[rStrip] + 0.5f * kFstStripPitchR;
            mFwdHitsFst.back()._localPosition[1] = muFstHit->getMeanPhiStrip() * kFstStripPitchPhi;
        }
    } // index

    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFst.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFst[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);
    }

    if ( numFwdHitsPost != numFwdHitsPrior ){
        if ( kLogLevel >= kLogInfo ){LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FST hits from MuDst" << endm;}
    }

    // TODO add to hitmap
    return count;
} // loadFstHitsFromMuDst

int StFwdHitLoader::loadFstHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
    int count = 0;
    //printf("AAA STEVENT loadFstHitsFromStEvent called, mStEvent=%p\n", (void*)mStEvent);
    if (!mStEvent) {
      if ( kLogLevel >= kLogVerbose ) LOG_WARN << "No StEvent, cannot load FST hits from StEvent StFstHitCollection" << endm;
        return 0;
    }
    if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "Got StEvent, loading Fst Hits" << endm;}
    StFstHitCollection *fstHitCollection = mStEvent->fstHitCollection();
    size_t numFwdHitsPrior = mFwdHitsFst.size();
    //printf("AAA STEVENT StFstHitCollection=%p  nHits=%d\n",
    //       (void*)fstHitCollection,
    //       fstHitCollection ? (int)fstHitCollection->numberOfHits() : -1);

    if ( fstHitCollection && fstHitCollection->numberOfHits() > 0){
        // reuse this to store cov mat
        TMatrixDSym hitCov3(3);
        if ( kLogLevel >= kLogVerbose ) {LOG_DEBUG << "StFstHitCollection is NOT NULL, loading FST hits from StEvent" << endm;}
        for ( unsigned int iw = 0; iw < kFstNumWedges; iw++ ){
            StFstWedgeHitCollection * wc = fstHitCollection->wedge( iw );
            if ( !wc ) continue;
            for ( unsigned int is = 0; is < kFstNumSensorsPerWedge; is++ ){
                StFstSensorHitCollection * sc = wc->sensor( is );
                if ( !sc ) continue;
                StSPtrVecFstHit fsthits = sc->hits();
                for ( unsigned int ih = 0; ih < fsthits.size(); ih++ ){
                    float vR   = fsthits[ih]->localPosition(0);
                    float vPhi = fsthits[ih]->localPosition(1);
                    float vZ   = fsthits[ih]->localPosition(2);
                    if ( kLogLevel >= kLogInfo ){LOG_INFO << TString::Format("FST local position: %f %f %f", vR, vPhi, vZ) << endm;}
                    
                    int wedgeIndex  = iw % kFstNumWedgePerDisk;
                    int sensorIndex = is % kFstNumSensorsPerWedge;
                    int diskIndex   = iw / kFstNumWedgePerDisk;
                    int globalIndex = FwdHit::fstGlobalSensorIndex( diskIndex, wedgeIndex, sensorIndex );

                    if ( kLogLevel >= kLogVerbose ) {LOG_INFO << "diskIndex = " << diskIndex << ", wedgeIndex = " << wedgeIndex << ", sensorIndex = " << sensorIndex << ", globalIndex = " << globalIndex << endm;}
                    float x0 = vR * cos( vPhi );
                    float y0 = vR * sin( vPhi );
                    //AAA: use positionError() from StFstFastSimMaker (= dr/sqrt12 = 0.830 cm)
                    //     instead of makeFstCovMat which used full pitch (2.875 cm) as sigma,
                    //     down-weighting FST hits 12x and degrading momentum resolution.
                    StThreeVectorF pe = fsthits[ih]->positionError();
                    float er   = pe.x(); // r error [cm]   (set by StFstFastSimMaker: dr/sqrt12)
                    float ephi = pe.y(); // phi error [rad] (set by StFstFastSimMaker: dphi/sqrt12)
                    if ( er <= 0 || ephi <= 0 ) {
                        // fallback: correct formula if positionError not set
                        er   = kFstStripPitchR   / sqrt(12.f);
                        ephi = kFstStripPitchPhi / sqrt(12.f);
                    }
                    //printf("AAA STEVENT FST disk=%d wedge=%d r=%.4f phi=%.5f z=%.2f er=%.4f ephi=%.6f\n",
                    //       diskIndex, wedgeIndex, sqrt(x0*x0+y0*y0), vPhi, vZ, er, ephi);
                    hitCov3 = makeFstCovMat( TVector3( x0, y0, vZ ), er, ephi );

                       
                    mSpacepointsFst.push_back( TVector3( x0, y0, vZ)  );
                    int track_id = fsthits[ih]->idTruth();
                    shared_ptr<McTrack> mcTrack = nullptr;
                    if ( mcTrackMap.count(track_id) ) {
                        mcTrack = mcTrackMap[track_id];
                        if ( kLogLevel >= kLogVerbose ) {
			  if ( kLogLevel >= kLogVerbose ) LOG_DEBUG << "FST Hit: idTruth = " << track_id << endm;
                          if ( kLogLevel >= kLogVerbose )  LOG_DEBUG << "Adding McTrack to FST hit" << endm;
                        }
                    }

                    // we use d+4 so that both FTT and FST start at 4 for historical reasons
                    mFwdHitsFst.push_back(
                        FwdHit(
                            count++, // id
                            x0, y0, vZ, // position
                            diskIndex+4, // volume id
                            kFstId, // detid
                            track_id, // mc track id
                            hitCov3, // covariance matrix
                            mcTrack // mcTrack
                        )
                    );
                    // store a pointer to the original StFstHit
                    mFwdHitsFst.back()._hit = fsthits[ih];
                    mFwdHitsFst.back()._genfit_plane_index = globalIndex;
                    // Store strip-native local position (no global azimuthal rotation).
                    // kFstrStart[] gives the inner edge of each r-strip in cm; the half-pitch
                    // centers the position within the strip.
                    {
                        int rStrip = (int)fsthits[ih]->getMeanRStrip();
                        mFwdHitsFst.back()._localPosition[0] = kFstrStart[rStrip] + 0.5f * kFstStripPitchR;
                        mFwdHitsFst.back()._localPosition[1] = fsthits[ih]->getMeanPhiStrip() * kFstStripPitchPhi;
                    }
                }
            } // loop is
        } // loop iw
        if ( kLogLevel >= kLogVerbose ) {LOG_DEBUG << " FOUND " << mSpacepointsFst.size() << " FST HITS in StFstHitCollection" << endm;}
    } // fstHitCollection
    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFst.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFst[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);
        // add to MC track map
        if ( hit->getMcTrack() ){
            hit->getMcTrack()->addFstHit(hit);
        }
    }
    if ( numFwdHitsPost != numFwdHitsPrior ){
        if ( kLogLevel >= kLogInfo ){LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FST hits from StEvent" << endm;}
    }
    return count;
} //loadFstHitsFromStEvent

int StFwdHitLoader::loadFstHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap ){
    int count = 0;
    if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "Looking for FST hits in geant struct" << endm;}
    /************************************************************/
    // Load FSI Hits from GEANT
    // St_g2t_fts_hit *mGeantFst = (St_g2t_fts_hit *)GetDataSet("geant/g2t_fsi_hit");

    if ( !mGeantFst ){
        if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "No g2t_fts_hits, cannot load FST hits from GEANT" << endm;}
        return 0;
    }

    int nfsi = mGeantFst->GetNRows();
    size_t numFwdHitsPrior = mFwdHitsFst.size();

    // reuse this to store cov mat
    TMatrixDSym hitCov3(3);

    for (int i = 0; i < nfsi; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)mGeantFst->At(i);

        if (0 == git)
            continue; // geant hit

        int track_id = git->track_p;
        int volume_id = git->volume_id;  // 4, 5, 6
        int d = volume_id / 1000;        // disk id

        int plane_id = d - 4;
        float x = git->x[0];
        float y = git->x[1];
        float z = git->x[2];

        if (plane_id > 3 || plane_id < 0) {
            continue;
        }

        //AAA NOTE: investigated correct strip-center snapping (kFstrStart[], 1536-strip phi)
        //  with both full-pitch and pitch/sqrt12 sigma — ALL variants gave worse results
        //  than the 3.0 cm rasterizer + full pitch sigma (REVERT A baseline).
        //  The tracking seeding/search is tuned to 3.0 cm r-positions; correct geometry
        //  disrupts it. Keep FstRasterizer for GEANT path until seeding is retuned.
        TVector3 rastered = mFstRasterizer.raster(TVector3(git->x[0], git->x[1], git->x[2]));
        x = rastered.X();
        y = rastered.Y();
        hitCov3 = makeFstCovMat( TVector3( x, y, z ) );
        std::shared_ptr<McTrack> mcTrack = nullptr;
        if ( mcTrackMap.count(track_id) ) {
            mcTrack = mcTrackMap[track_id];
            if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "Adding McTrack to FTT hit: " << track_id << endm;}
        }
        mFwdHitsFst.push_back(
            FwdHit(
                count++, // id
                x, y, z, // position
                d, // volume id
                kFstId, // detid
                track_id, // mc track id
                hitCov3, // covariance matrix
                mcTrack // mcTrack
            )
        );
        mSpacepointsFst.push_back( TVector3( x, y, z )  );
    }

    // this has to be done AFTER because the vector reallocates mem when expanding, changing addresses
    size_t numFwdHitsPost = mFwdHitsFst.size();
    for ( size_t iFwdHit = numFwdHitsPrior; iFwdHit < numFwdHitsPost; iFwdHit++){
        FwdHit *hit = &(mFwdHitsFst[ iFwdHit ]);
        // Add the hit to the hit map
        if ( hit->getLayer() >= 0 )
            hitMap[hit->getSector()].push_back(hit);

        // add to MC track map
        if ( hit->getMcTrack() )
            hit->getMcTrack()->addFstHit(hit);
    }
    if ( numFwdHitsPost != numFwdHitsPrior ){
        if ( kLogLevel >= kLogInfo ){LOG_INFO << "Loaded " << numFwdHitsPost - numFwdHitsPrior << " FST hits from GEANT" << endm;}
    }

    return count;
} // loadFstHitsFromGEANT

int StFwdHitLoader::loadEpdHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, StFcsDb *fcsDb  ){

    if ( mEpdDataSource == StFwdHitLoader::DataSource::GEANT ){
        if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "Loading EPD Hits from GEANT (Not Yet Implemented)" << endm;}
        // return loadEpdHitsFromGEANT( mcTrackMap, hitMap );
    } else if ( mEpdDataSource == StFwdHitLoader::DataSource::STEVENT ){
        if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "Loading EPD Hits from StEvent" << endm;}
        return loadEpdHitsFromStEvent( mcTrackMap, hitMap, fcsDb );
    } else if ( mEpdDataSource == StFwdHitLoader::DataSource::MUDST ){
        if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << "Loading EPD Hits from MuDst (Not Yet Implemented)" << endm;}
        // loadEpdHitsFromMuDst( mcTrackMap, hitMap, count );
        return 0;
    }

    return 0;
}
int StFwdHitLoader::loadEpdHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, StFcsDb *fcsDb  ){

    if ( !mStEvent ) {
      if ( kLogLevel >= kLogVerbose ) LOG_ERROR << "StEvent is NULL, cannot load EPD hits" << endm;
        return 0;
    }
    if ( !fcsDb ) {
      if ( kLogLevel >= kLogVerbose ) LOG_ERROR << "StFcsDb is NULL, cannot load EPD hits" << endm;
        return 0;
    }
    StFcsCollection *fcsCol = mStEvent->fcsCollection();

    if (!fcsCol) {
        if ( kLogLevel >= kLogInfo ){LOG_INFO << "StFcsCollection is NULL, cannot load EPD hits" << endm;}
        return 0;
    }

    StEpdGeom epdgeo;
    int count = 0; // counts hits above threshold
    // LOAD PRESHOWER HITS (EPD)
    for ( int det = kFcsPresNorthDetId; det < kFcsPresSouthDetId + 1; det ++ ) {

        StSPtrVecFcsHit& hits = mStEvent->fcsCollection()->hits(det);
        int nh=fcsCol->numberOfHits(det);
        for ( int i = 0; i < nh; i++ ){
            StFcsHit* hit=hits[i];

            if(det!=kFcsPresNorthDetId && det!=kFcsPresSouthDetId) continue;
            if ( hit->energy() < mEpdThreshold ) continue;
            double zepd=375.0;
            int pp,tt,n;
            double x[5],y[5];

            fcsDb->getEPDfromId(det,hit->id(),pp,tt);
            epdgeo.GetCorners(100*pp+tt,&n,x,y);
            double x0 = (x[0] + x[1] + x[2] + x[3]) / 4.0;
            double y0 = (y[0] + y[1] + y[2] + y[3]) / 4.0;
            mSpacepointsEpd.push_back( TVector3( x0, y0, zepd ) );

            // Covariance from tile corners: for a uniform distribution over a
            // quadrilateral, C = (1/3) * mean(corner outer-products) = sum/12.
            // This naturally gives the anisotropic, correlated (x,y) uncertainty
            // from the actual tile shape without any hardcoded approximation.
            TMatrixDSym hitCov3(3);
            double cxx=0, cyy=0, cxy=0;
            for(int k=0;k<4;k++){
                double dx=x[k]-x0, dy=y[k]-y0;
                cxx+=dx*dx; cyy+=dy*dy; cxy+=dx*dy;
            }
            hitCov3(0, 0) = cxx / 12.0;
            hitCov3(1, 1) = cyy / 12.0;
            hitCov3(0, 1) = cxy / 12.0;
            hitCov3(1, 0) = cxy / 12.0;
            hitCov3(2, 2) = 1.0; // σ_z² = 1 cm² (tile thickness ≈ 2 cm / √12 ≈ 0.58 cm; conservative)

            const vector<pair<unsigned int, float>> gt = hit->getGeantTracks();
            int track_id = -1;
            shared_ptr<McTrack> mcTrack = nullptr;
            if ( gt.size() > 0 ){
                track_id = gt[0].first;
                
                if ( mcTrackMap.count(track_id) ) {
                    mcTrack = mcTrackMap[track_id];
                    LOG_DEBUG << "Adding McTrack to EPD hit: " << track_id << endm;
                }
            }
            // TODO: Make FwdHits
            mFwdHitsEpd.push_back(
                FwdHit(
                    count++, // id
                    x0, y0, zepd, // position
                    -7, // volume id
                    kFcsPresId, // detid
                    0, // track id
                    hitCov3, // covariance matrix
                    mcTrack // mcTrack
                )
            );
            if ( kLogLevel >= kLogVerbose ){LOG_DEBUG << TString::Format("Adding Epd as FwdHit: %f %f %f", x0, y0, zepd) << endm;}
        } // for i
    } // for det

    for ( size_t iFwdHit = 0; iFwdHit < mFwdHitsEpd.size(); iFwdHit++){
        FwdHit *hit = &(mFwdHitsEpd[ iFwdHit ]);
        // EPD is 7 (3 FST, 4 FTT, starting at 0)
        hitMap[7].push_back(hit);

        // add to MC track map
        // if ( hit->getMcTrack() )
        //     hit->getMcTrack()->addFstHit(hit);
    }

    return mFwdHitsEpd.size();
}
