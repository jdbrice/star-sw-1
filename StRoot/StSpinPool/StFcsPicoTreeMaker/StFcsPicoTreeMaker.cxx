#include "StFcsPicoTreeMaker.h"

#include "StEnumerations.h"
#include "StMessMgr.h"
#include "Stypes.h"
#include "StEventTypes.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StThreeVectorF.hh"

#include "StFcsDbMaker/StFcsDb.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsCollection.h"
#include "StFcsHit.h"
#include "StFcsCluster.h"
#include "StFcsPoint.h"

StFcsPicoTreeMaker::StFcsPicoTreeMaker(const char* name):StMaker(name)
{
}

StFcsPicoTreeMaker::~StFcsPicoTreeMaker()
{
  delete mDataTree;
  if( mFile!=0 ){mFile->Close();}
  delete mFile;
}

StFcsPicoTree* StFcsPicoTreeMaker::setTree(const char* treename, const char* title,Int_t splitlevel)
{
  if( mDataTree==0 ){ mDataTree = new StFcsPicoTree(treename,title,splitlevel); }
  return mDataTree;
}

void StFcsPicoTreeMaker::turnOffHits()
{
  mDoHits = false;
}
void StFcsPicoTreeMaker::turnOffClusters()
{
  mDoClus = false;
}
void StFcsPicoTreeMaker::turnOffPoints()
{
  mDoPoints = false;
}

Int_t StFcsPicoTreeMaker::Init()
{
  if( !mFcsDb ){ mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb")); }
  //mFcsDb->setDbAccess(0);
  if (!mFcsDb) {
    LOG_ERROR << "StFcsPicoTreeMaker::InitRun Failed to get StFcsDb" << endm;
    return kStFatal;
  }
  if( !mDataTree ){ mDataTree = new StFcsPicoTree("FcsPicoTree",""); }  //Don't use exact class name to avoid confusion
  return kStOk;
}

Int_t StFcsPicoTreeMaker::Finish()
{
  if( mFileName.Length()==0 ){ return kStOk; }
  mFile = new TFile(mFileName.Data(), "RECREATE");
  if( mDataTree ){ mDataTree->Write(); }
  return kStOk;
}

Int_t StFcsPicoTreeMaker::Make()
{
  StEvent* event = (StEvent*)GetInputDS("StEvent");
  if (!event) {
    LOG_ERROR << "StFcsPicoTreeMaker::Make did not find StEvent" << endm;
    return kStErr;
  }
  mFcsColl = event->fcsCollection();
  if (!mFcsColl) {
    LOG_ERROR << "StFcsPicoTreeMaker::Make did not find StEvent->StFcsCollection" << endm;
    return kStErr;
  }
  
  if( !mDoHits && !mDoClus && !mDoPoints ){ LOG_WARN << "StFcsPicoTreeMaker::Make() All tree data turned off skipping Make()" << endm; return kStOk; }

  int totalhits = 0;   //Running counter for all hits for all detectors
  int totalclus = 0;   //Running counter for all clusters for all detectors
  int totalpoints = 0; //Running counter for all points for all detectors
  for( int det=0; det<kFcsNDet; det++ ){
    
    StSPtrVecFcsHit& hits = mFcsColl->hits(det);
    int nh = mFcsColl->numberOfHits(det);
    for( int ihit=0; ihit<nh && mDoHits; ++ihit ){
      StFcsHit* hit=hits[ihit];
      StFcsPicoHit* picohit = mDataTree->ConstructedHit(totalhits++);
      picohit->mDetId  = hit->detectorId();
      picohit->mChId   = hit->id();
      picohit->mAdcSum = hit->adcSum();
      picohit->mEnergy = hit->energy();

      if( det<=kFcsHcalSouthDetId ){ //ECAL and HCAL
	StThreeVectorF xyz = mFcsDb->getStarXYZ(hit);
	picohit->mXstar = xyz.x();
	picohit->mYstar = xyz.y();
	picohit->mZstar = xyz.z();
      }
      else if(det==kFcsPresNorthDetId || det==kFcsPresSouthDetId){//EPD as Pres
	picohit->mZstar = 375.0; //Taken from StFcsEventDisplay for zepd this is in cm
	}
      //std::cout << "|i:"<<ihit<<"|th:"<<totalhits<<"("<<picohit->mXstar<<","<<picohit->mYstar<<","<<picohit->mEnergy<<")"<<std::endl;
    }
    
    StSPtrVecFcsCluster& clusters = mFcsColl->clusters(det);
    int nc = mFcsColl->numberOfClusters(det);
    for( int iclus=0; iclus<nc && mDoClus; ++iclus ){
      StFcsCluster* cluster = clusters[iclus];
      
      StFcsPicoCluster* picoclus = mDataTree->ConstructedCluster(totalclus++);
      picoclus->mId             = cluster->id();
      picoclus->mDetId          = cluster->detectorId();
      picoclus->mCategory       = cluster->category();
      picoclus->mNTowers        = cluster->nTowers();
      picoclus->mNNeighbor      = cluster->nNeighbor();
      picoclus->mNPoints        = cluster->nPoints();
      picoclus->mEnergy         = cluster->energy();
      picoclus->mX              = cluster->x();
      picoclus->mY              = cluster->y();
      picoclus->mSigmaMin       = cluster->sigmaMin();
      picoclus->mSigmaMax       = cluster->sigmaMax();
      picoclus->mTheta          = cluster->theta();
      picoclus->mChi2Ndf1Photon = cluster->chi2Ndf1Photon();
      picoclus->mChi2Ndf2Phoron = cluster->chi2Ndf2Photon();
      //Lorentz 4 momentum of cluster
      StLorentzVectorD clusp = cluster->fourMomentum();
      picoclus->mPx = clusp.px();
      picoclus->mPy = clusp.py();
      picoclus->mPz = clusp.pz();
      picoclus->mE  = clusp.e();
    }

    StSPtrVecFcsPoint& points = mFcsColl->points(det);
    int np = mFcsColl->numberOfPoints(det);
    for( int ipoint=0; ipoint<np && mDoPoints; ++ipoint ){
      StFcsPoint* point=points[ipoint];
      
      StFcsPicoPoint* picopoint = mDataTree->ConstructedPoint(totalpoints++);
      picopoint->mDetId                 = point->detectorId();
      picopoint->mEnergy                = point->energy();
      picopoint->mXlocal                = point->x();
      picopoint->mYlocal                = point->y();
      picopoint->mParentClusterId       = point->parentClusterId();
      picopoint->mNParentClusterPhotons = point->nParentClusterPhotons();
      //STAR xyx
      StThreeVectorF pointxyz = point->xyz();
      picopoint->mXstar = pointxyz.x();
      picopoint->mYstar = pointxyz.y();
      picopoint->mZstar = pointxyz.z();
      //Lorentz 4 momentum of point
      StLorentzVectorD pointp = point->fourMomentum();
      picopoint->mPx = pointp.px();
      picopoint->mPy = pointp.py();
      picopoint->mPz = pointp.pz();
      picopoint->mE  = pointp.e();
    }
  }
  mDataTree->Fill();
  mDataTree->ClearAll();

  return kStOk;
}
