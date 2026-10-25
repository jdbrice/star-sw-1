#include "Stypes.h"
#include "TNamed.h"
#include "TGeoMatrix.h"

#include "StFstHitMaker.h"
#include "StFstUtil/StFstCollection.h"
#include "StFstUtil/StFstCluster.h"
#include "StFstUtil/StFstClusterCollection.h"
#include "StFstHit.h"
#include "StFstHitCollection.h"

#include "St_base/StMessMgr.h"
#include "StEvent.h"
#include "StEventTypes.h"
#include "StContainers.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StFstConsts.h"

#include "StFstDbMaker/StFstDb.h"
#include "tables/St_fstControl_Table.h"

ClassImp(StFstHitMaker);


StFstHitMaker::StFstHitMaker( const char *name ) : StMaker(name), mSensorTransforms(0), mApplySurvey(kFALSE)
{}


Int_t StFstHitMaker::InitRun(Int_t runnumber)
{
	TObjectSet *fstDbDataSet = (TObjectSet *) GetDataSet("fst_db");
	StFstDb    *fstDb = 0;

	if (fstDbDataSet) {
		fstDb = (StFstDb *) fstDbDataSet->GetObject();
		assert(fstDb);
	}
	else {
		LOG_ERROR << "InitRun : no fstDb" << endm;
		return kStErr;
	}

	// geometry Db tables
	mSensorTransforms = fstDb->getRotations();

	return kStOk;
}


/**
 * Takes the StFstClusterCollection from StFstCollection that is normally
 * created by the StFstClusterMaker and fills StEvent's StFstHitCollection which
 * is used in tracking.
 */
Int_t StFstHitMaker::Make()
{
	// Obtain hit collection
	StEvent *eventPtr = (StEvent *) GetDataSet("StEvent");

	if (!eventPtr) {
		LOG_ERROR << "Make() - No StEvent found in the chain. Cannot proceed" << endm;
		return kStErr;
	}

	//input clusters info.
	TObjectSet *fstDataSet = (TObjectSet *) GetDataSet("fstRawHitAndCluster");

	if (!fstDataSet) {
		LOG_WARN << "Make() - fstRawHitAndCluster dataset not found. No FST hits will be available for tracking" << endm;
		return kStWarn;
	}

	StFstCollection *fstCollectionPtr = (StFstCollection *) fstDataSet->GetObject();

	if ( !fstCollectionPtr ) {
		LOG_WARN << "Make() - StFstCollection not found. No FST hits will be available for tracking" << endm;
		return kStWarn;
	}

	// Get pointer to an existing StFstHitCollection if any
	StFstHitCollection *fstHitCollection = eventPtr->fstHitCollection();

	// If no fst hit collection, create one
	if (!fstHitCollection) {
		fstHitCollection = new StFstHitCollection();
		eventPtr->setFstHitCollection(fstHitCollection);
		LOG_DEBUG << "Make() - Added new StFstHitCollection to this StEvent" << endm;
	}

	unsigned char  nClusteringType = -1;

	for (unsigned char wedgeIdx = 0; wedgeIdx < kFstNumWedges; ++wedgeIdx) {
		//add new hits from clusters

		StFstClusterCollection *clusterCollectionPtr = fstCollectionPtr->getClusterCollection(wedgeIdx );

		if ( clusterCollectionPtr ) {
			unsigned int numClusters = clusterCollectionPtr->getNumClusters();
			LOG_DEBUG << "Make() - Number of clusters found in wedge " << (int)(wedgeIdx + 1) << ": " << numClusters << endm;

			unsigned short idTruth = 0;
			unsigned char  nRawHits = -1, nRawHitsR = -1, nRawHitsPhi = -1;
			unsigned char  disk = -1, wedge = -1, sensor = -1, apv = -1;
			int  meanRStrip = -1, meanPhiStrip = -1;
                        float  charge = 0., chargeErr = 0.;
			unsigned char  maxTb = -1;
			int	 key = -1;

			for (std::vector< StFstCluster * >::iterator clusterIter = clusterCollectionPtr->getClusterVec().begin(); clusterIter != clusterCollectionPtr->getClusterVec().end(); ++clusterIter) {
				idTruth         = (*clusterIter)->getIdTruth();
				key             = (*clusterIter)->getKey();
				disk            = (*clusterIter)->getDisk();
				wedge           = (*clusterIter)->getWedge();
				sensor          = (*clusterIter)->getSensor();
				apv             = (*clusterIter)->getApv();
				meanRStrip      = (*clusterIter)->getMeanRStrip();
				meanPhiStrip    = (*clusterIter)->getMeanPhiStrip();
				maxTb           = (*clusterIter)->getMaxTimeBin();
				charge          = (*clusterIter)->getTotCharge();
				chargeErr       = (*clusterIter)->getTotChargeErr();
				nRawHits        = (*clusterIter)->getNRawHits();
				nRawHitsR       = (*clusterIter)->getNRawHitsR();
				nRawHitsPhi     = (*clusterIter)->getNRawHitsPhi();
				nClusteringType = (*clusterIter)->getClusteringType();

				StFstHit *newHit = new StFstHit(disk, wedge, sensor, apv, charge, chargeErr, maxTb, meanRStrip, meanPhiStrip, nRawHits, nRawHitsR, nRawHitsPhi);
				newHit->setId(key);
				newHit->setIdTruth(idTruth);

				// strip -> nominal position, from the ONE place that owns the
				// convention. This block used to hold it inline: phiInner and
				// phiOuter from kFstphiStart/kFstphiStop, the per-wedge and
				// per-disk flips from kFstzFilp x kFstzDirct, and the
				// outer-sensor kFstStripGapPhi -- all of which
				// StFstSlowSimMaker restates inverted and StFwdHitLoader
				// patches afterwards. StFstDb::stripPosition() is that algebra
				// moved verbatim, so the producer, the simulation and anything
				// rebuilding a hit from its strips cannot drift apart. The z it
				// returns is sensorZ(), so it follows fstmConfig(): the silicon
				// the strip is on, on the face the wedge is really mounted.
				Double_t g[3] = {0, 0, 0};
				if (!StFstDb::stripPosition(disk, (int)newHit->getWedge(), sensor,
				                            meanRStrip, meanPhiStrip, g)) {
					LOG_WARN << "Make() - no strip position for disk " << (int)disk
					         << " wedge " << (int)newHit->getWedge() << " sensor "
					         << (int)sensor << " rStrip " << meanRStrip
					         << "; hit dropped" << endm;
					delete newHit;
					continue;
				}
				newHit->setLocalPosition(g[0], g[1], g[2]); //nominal-frame position
				fstHitCollection->addHit(newHit);
			} //cluster loop over
		}//end clusterCollectionPtr

		//set global position
		StFstWedgeHitCollection *wedgeHitCollection = fstHitCollection->wedge(wedgeIdx);
		for (int sensorIdx = 0; sensorIdx < kFstNumSensorsPerWedge; sensorIdx++) {
			StFstSensorHitCollection *sensorHitCollection = wedgeHitCollection->sensor(sensorIdx);

			for (int idx = 0; idx < (int) sensorHitCollection->hits().size(); idx++ ) {
				StFstHit *newHit = sensorHitCollection->hits()[idx];
				double local[3];
				double global[3];
				local[0] = newHit->localPosition(0);
				local[1] = newHit->localPosition(1);
				local[2] = newHit->localPosition(2);

				// nominal -> global. The survey matrix is applied ONLY when asked
				// for: see SetApplySurvey() in the header for the two reasons
				// behind that default. Both are about the simulated geometry and
				// the measured residuals; the frame question that used to be the
				// first of them is settled -- StFstDb no longer puts tpcOnGlobal
				// in the product, so this matrix is exactly the displacement AGML
				// gives the TGeo nodes the GenFit planes are built from.
				//
				// Both of the things that used to live here hid the survey rather
				// than applying it. Before 9e70bfd0a9 the LocalToMaster result was
				// overwritten by three lines converting (local[0], local[1]) as
				// polar (r, phi), which was consistent with a producer that stored
				// the global (r, phi, z) but discarded the matrix. 9e70bfd0a9 moved
				// the producer to cartesian via MasterToLocal, which made this an
				// exact round trip -- so the matrix cancelled either way, and no
				// FST hit has ever carried the survey. For the Run 22 tables that
				// is (+0.7126, -3.0102) mm for the x>0 half and (+2.2626, -3.3302)
				// mm for the x<0 half -- which is, to the bit, how far the
				// misaligned AGML geometry moves the sensors of those same two
				// groups (d4378cbde9; measured off the fGeom.root cache, 54 nodes
				// each).
				if (mApplySurvey) {
					int sensorId = StFstDb::sensorSurveyId((int)newHit->getWedge(), (int)newHit->getSensor());
					TGeoHMatrix *geoMSensorOnGlobal = (TGeoHMatrix *) mSensorTransforms->FindObject(Form("R%04i", sensorId));

					if (geoMSensorOnGlobal) {
						geoMSensorOnGlobal->LocalToMaster(local, global);
					}
					else {
						LOG_WARN << "StFstHitMaker: no survey matrix R" << sensorId
						         << ", leaving the hit at its nominal position" << endm;
						global[0] = local[0]; global[1] = local[1]; global[2] = local[2];
					}
				}
				else {
					global[0] = local[0];
					global[1] = local[1];
					global[2] = local[2];
				}

				StThreeVectorF vecGlobal(global);
				newHit->setPosition(vecGlobal); //set global position
			}//end sensor hit collection
		}//end wedge hit collection
	} //wedge loop over

	fstHitCollection->setClusteringType(nClusteringType);

	return kStOk;
}
