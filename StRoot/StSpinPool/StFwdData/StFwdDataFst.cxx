#include "StFwdDataFst.h"

ClassImp(StFstRawHitInfo)

StFstRawHitInfo::StFstRawHitInfo()
{
  memset( mCharge,0,sizeof(mCharge) );
  memset( mChargeErr,0,sizeof(mChargeErr) );
}

StFstRawHitInfo::~StFstRawHitInfo()
{}

void StFstRawHitInfo::Copy(TObject& object) const
{
  ((StFstRawHitInfo&)object).mChannelId   = mChannelId;
  ((StFstRawHitInfo&)object).mGeoId       = mGeoId;
  ((StFstRawHitInfo&)object).mWedge       = mWedge;
  ((StFstRawHitInfo&)object).mSensor      = mSensor;
  ((StFstRawHitInfo&)object).mPhiStrip    = mPhiStrip;
  ((StFstRawHitInfo&)object).mRStrip      = mRStrip;
  ((StFstRawHitInfo&)object).mRdo         = mRdo;
  ((StFstRawHitInfo&)object).mArm         = mArm;
  ((StFstRawHitInfo&)object).mApv         = mApv;
  ((StFstRawHitInfo&)object).mChannel     = mChannel;
  ((StFstRawHitInfo&)object).mSeedHitFlag = mSeedHitFlag;

  ((StFstRawHitInfo&)object).mMaxTimeBin  = mMaxTimeBin;
  for( int i=0; i<kFstNumTimeBins; ++i ){
    (((StFstRawHitInfo&)object)).mCharge[i] = mCharge[i];
    (((StFstRawHitInfo&)object)).mChargeErr[i] = mChargeErr[i];
  }
}

void StFstRawHitInfo::Clear(Option_t* opt)
{
  mChannelId   = -1;
  mGeoId       = -1;
  mWedge       = -1;
  mSensor      = -1;
  mPhiStrip    = -1;
  mRStrip      = -1;
  mRdo         = -1;
  mArm         = -1;
  mApv         = -1;
  mChannel     = -1;
  mSeedHitFlag = -1;

  mMaxTimeBin  = -1;
  
  memset( mCharge,0,sizeof(mCharge) );
  memset( mChargeErr,0,sizeof(mChargeErr) );
}

void StFstRawHitInfo::Print(Option_t* opt) const
{
  std::cout << "|ChId:"<<mChannelId
	    << "|GeoId:"<<mGeoId
	    << "|Wedge:"<<mWedge
	    << "|Sensor:"<<mSensor
	    << "|PhiStrip:"<<mPhiStrip
	    << "|RStrip:"<<mRStrip
	    << "|Rdo:"<<mRdo
	    << "|Arm:"<<mArm
	    << "|Apv:"<<mApv
	    << "|Ch:"<<mChannel
	    << "|SeedFlag:"<<mSeedHitFlag
	    << "|MaxTb:"<<mMaxTimeBin
	    << std::endl;
  TString option(opt);
  option.ToLower();
  if( option.Contains("charge") ){
    std::cout << "  + ";
    for( int i=0; i<kFstNumTimeBins; ++i ){
      std::cout << i << ":"<<mCharge[i] <<"pm"<<mChargeErr[i] <<"|";
    }
    std::cout << std::endl;
  }
}

ClassImp(StFstHitInfo)

StFstHitInfo::StFstHitInfo()
{}

StFstHitInfo::~StFstHitInfo()
{}

void StFstHitInfo::Copy(TObject& object) const
{
  ((StFstHitInfo&)object).mHitId          = mHitId;
  ((StFstHitInfo&)object).mWedge          = mWedge;
  ((StFstHitInfo&)object).mSensor         = mSensor;
  ((StFstHitInfo&)object).mApv            = mApv;
  ((StFstHitInfo&)object).mMaxTimeBin     = mMaxTimeBin;
  ((StFstHitInfo&)object).mClusteringType = mClusteringType;
  ((StFstHitInfo&)object).mNRawHits       = mNRawHits;
  ((StFstHitInfo&)object).mNRawHitsR      = mNRawHitsR;
  ((StFstHitInfo&)object).mNRawHitsPhi    = mNRawHitsPhi;

  ((StFstHitInfo&)object).mMeanPhiStrip   = mMeanPhiStrip;
  ((StFstHitInfo&)object).mMeanRStrip     = mMeanRStrip;
  ((StFstHitInfo&)object).mLocalR         = mLocalR;
  ((StFstHitInfo&)object).mLocalPhi       = mLocalPhi;
  ((StFstHitInfo&)object).mLocalZ         = mLocalZ;

  ((StFstHitInfo&)object).mX              = mX;
  ((StFstHitInfo&)object).mY              = mY;
  ((StFstHitInfo&)object).mZ              = mZ;

  ((StFstHitInfo&)object).mTotCharge      = mTotCharge;
  ((StFstHitInfo&)object).mTotChargeErr   = mTotChargeErr;
}

void StFstHitInfo::Clear(Option_t* opt)
{
  mHitId          = -1;
  mWedge          = -1;
  mSensor         = -1;
  mApv            = -1;
  mMaxTimeBin     = -1;
  mClusteringType = -1;
  mNRawHits       = -1;
  mNRawHitsR      = -1;
  mNRawHitsPhi    = -1;

  mMeanPhiStrip   = -1;
  mMeanRStrip     = -1;
  mLocalR         = -1;
  mLocalPhi       = -1;
  mLocalZ         = -1;

  mX              = 0;
  mY              = 0;
  mZ              = 0;

  mTotCharge      = 0;
  mTotChargeErr   = 0;
}

void StFstHitInfo::Print(Option_t* opt) const
{
  std::cout << "|HitId:"<<mHitId
	    << "|Wedge:"<<mWedge
	    << "|Sensor:"<<mSensor
	    << "|Apv:"<<mApv
	    << "|MaxTb:"<<mMaxTimeBin
	    << "|ClusType:"<<mClusteringType
	    << "|NRawHits:"<<mNRawHits
	    << "|NRH_R:"<<mNRawHitsR
	    << "|NRH_Phi:"<<mNRawHitsPhi
	    << "|MeanPhiStrip:"<<mMeanPhiStrip
	    << "|MeanRStrip:"<<mMeanRStrip
	    << "|LocalR:"<<mLocalR
	    << "|LocalPhi:"<<mLocalPhi
	    << "|LocalZ:"<<mLocalZ
	    << "|("<<mX <<","<<mY<<","<<mZ<<")"
	    << "|TotCharge:"<<mTotCharge
	    << "|TotChargeErr:"<<mTotChargeErr
	    << std::endl;
}

