#include "StFwdDataFtt.h"

ClassImp(StFttRawHitInfo)

StFttRawHitInfo::StFttRawHitInfo()
{}

StFttRawHitInfo::~StFttRawHitInfo()
{}

void StFttRawHitInfo::Copy(TObject& object) const
{
  ((StFttRawHitInfo&)object).mQuadrant    = mQuadrant;
  ((StFttRawHitInfo&)object).mPlane       = mPlane;
  ((StFttRawHitInfo&)object).mRow         = mRow;
  ((StFttRawHitInfo&)object).mStrip       = mStrip;
  ((StFttRawHitInfo&)object).mOrientation = mOrientation;
  ((StFttRawHitInfo&)object).mAdc         = mAdc;
  ((StFttRawHitInfo&)object).mBcid        = mBcid;
  ((StFttRawHitInfo&)object).mBcidDelta   = mBcidDelta;
  ((StFttRawHitInfo&)object).mTb          = mTb;
  ((StFttRawHitInfo&)object).mTime        = mTime;
}

void StFttRawHitInfo::Clear(Option_t* opt)
{
  mQuadrant    = kFttUnknownQuadrant;
  mPlane       = 255;
  mRow         = 255;
  mStrip       = 255;
  mOrientation = kFttUnknownOrientation;
  mAdc        = 0;
  mBcid       = 0;
  mBcidDelta   = -32000;
  mTb          = -32000;
  mTime        = -32000;
}

void StFttRawHitInfo::Print(Option_t* opt) const
{
  std::cout << "|Quad:"<<mQuadrant
	    << "|Plane:"<<mPlane
	    << "|Row:"<<mRow
	    << "|Strip:"<<mStrip
	    << "|Orientation:"<<mOrientation
	    << "|Adc:"<<mAdc
	    << "|Bcid:"<<mBcid
	    << "|BcidDelta:"<<mBcidDelta
	    << "|Tb:"<<mTb
	    << "|Time:"<<mTime
	    << std::endl;
}


