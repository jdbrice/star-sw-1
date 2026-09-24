/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  Classes related to hold FST information

  DESCRIPTION
  Contains the classes #StFstRawHitInfo and #StFstHitInfo which are essentially mirror copies of #StMuFstRawHit and #StMuFstHit. In this way they can be invoked to fill a TTree in #StFwdAnaData. The FST raw hits denote the single strip level information and bare information. The 'hits' are actually clusters of raw hits and the only difference between clusters and hits is that hits contain information about the cluster in STAR coordinate space.

  LOG
  @[July 30, 2026] > First instance
*/

#ifndef STFWDDATA_STFWDDATAFST_HH
#define STFWDDATA_STFWDDATAFST_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TObject.h"
#include "TString.h"

//STAR Headers
#include "StEvent/StFstConsts.h"

class StFstRawHitInfo : public TObject
{
  StFstRawHitInfo();
  ~StFstRawHitInfo();
  
  Int_t mChannelId   = -1;
  Int_t mGeoId       = -1;
  Int_t mWedge       = -1;
  Int_t mSensor      = -1;
  Int_t mPhiStrip    = -1;
  Int_t mRStrip      = -1;
  Int_t mRdo         = -1;
  Int_t mArm         = -1;
  Int_t mApv         = -1;
  Int_t mChannel     = -1;
  Int_t mSeedHitFlag = -1;

  Int_t mMaxTimeBin  = -1;
  Float_t mCharge[kFstNumTimeBins];
  Float_t mChargeErr[kFstNumTimeBins];

  virtual void Copy(TObject& object) const;           ///< Copy this candidate to object
  virtual void Clear(Option_t* opt="");               ///< Resets all variables to defaults
  virtual void Print(Option_t* opt="") const;         ///< Print all variables, use option "charge" to print out the charge arrays

  ClassDef(StFstRawHitInfo,1);
};

class StFstHitInfo : public TObject
{
  StFstHitInfo();
  ~StFstHitInfo();
  
  Int_t mHitId          = -1;
  Int_t mWedge          = -1;
  Int_t mSensor         = -1;
  Int_t mApv            = -1;
  Int_t mMaxTimeBin     = -1;
  Int_t mClusteringType = -1;
  Int_t mNRawHits       = -1;
  Int_t mNRawHitsR      = -1;
  Int_t mNRawHitsPhi    = -1;
  
  Float_t mMeanPhiStrip = -1;
  Float_t mMeanRStrip   = -1;
  Float_t mLocalR       = -1;
  Float_t mLocalPhi     = -1;
  Float_t mLocalZ       = -1;

  Float_t mX            = 0;
  Float_t mY            = 0;
  Float_t mZ            = 0;

  Float_t mTotCharge    = 0;
  Float_t mTotChargeErr = 0;

  virtual void Copy(TObject& object) const;           ///< Copy this candidate to object
  virtual void Clear(Option_t* opt="");               ///< Resets all variables to defaults
  virtual void Print(Option_t* opt="") const;         ///< Print all variables no options

  ClassDef(StFstHitInfo,1);
};

#endif

