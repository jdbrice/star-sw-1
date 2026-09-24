/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  Classes related to hold FTT (sTGC (small-strip Thin Gap Chamber)) information

  DESCRIPTION
  Contains copies of most variables in StFttRawHit

  LOG
  @[August 4, 2026] > First instance
*/

#ifndef STFWDDATA_STFWDDATAFTT_HH
#define STFWDDATA_STFWDDATAFTT_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TObject.h"
#include "TString.h"

//STAR Headers
#include "StEvent/StEnumerations.h"

class StFttRawHitInfo : public TObject
{
  StFttRawHitInfo();
  ~StFttRawHitInfo();

  UChar_t mQuadrant    = kFttUnknownQuadrant;
  UChar_t mPlane       = 255;
  UChar_t mRow         = 255;
  UChar_t mStrip       = 255;
  UChar_t mOrientation = kFttUnknownOrientation;
  UShort_t mAdc        = 0;
  UShort_t mBcid       = 0;                          ///< time from VMM clock which is an independent 40 MHz clock on VMM boards (high resolution clock)
  Short_t mBcidDelta   = -32000;                     ///< Delta BCID is coming from DAQ
  Short_t mTb          = -32000;                     ///< Digitized time from STAR (low resolution clock)
  Short_t mTime        = -32000;                     ///< Difference between DeltaBCID and the calibrated Delta BCID; ideally, after calibration -3<time<3 should be the triggered bunch crossing

  virtual void Copy(TObject& object) const;           ///< Copy this candidate to object
  virtual void Clear(Option_t* opt="");               ///< Resets all variables to defaults
  virtual void Print(Option_t* opt="") const;         ///< Print all variables, use option "charge" to print out the charge arrays

  ClassDef(StFttRawHitInfo,1);
};
  //StFttRawHit::tb is coming from DAQ low resolution clock
  //BCID is digitized time coming from VMM 40 MHz clock
  //Delta BCID is coming from DAQ
  //StFttRawHit::time is difference between Delta BCID and tb


#endif

