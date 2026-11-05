#ifndef StFstRawHitMaker_hh
#define StFstRawHitMaker_hh

#ifndef __CINT__
#include <array>
#include <iostream>
#endif

#include <vector>

#include "St_base/Stypes.h"
#include "StChain/StRTSBaseMaker.h"
#include "StEvent/StFstConsts.h"


class StEvent;
class StFstCollection;
class StFstEvtCollection;


/**
 * Reads/un-packs a DAQ or SFS file; Accesses calibration DBs;
 * Marks out bad channels/chips; Creates and fills the FST containers.
 * Two working modes included:
 * 1) calibration mode: save raw ADC value only
 * 2) non-calibration (physics) mode: access calibration datasets, signal-like
 *    raw hit decision, full raw hit info., ...
 *
 * \Author Shenghui Zhang
 * \date Aug 2021
 */
class StFstRawHitMaker : public StRTSBaseMaker
{
public:
   StFstRawHitMaker( const char *name = "fst_raw_hit" );
   ~StFstRawHitMaker();
   Int_t Init();
   Int_t InitRun(Int_t runNumber);
   Int_t Make();
   Int_t Finish();
   void Clear( Option_t *opts = "" );

   void setIsCalibrationMode( bool isCaliMode = false )  { mIsCaliMode = isCaliMode; }
   void setMinHitCut(double minhitCut = 2.5)             { mMinHitCut = minhitCut; }
   void setMedHitCut(double medhitCut = 3.5)             { mMedHitCut = medhitCut; }
   void setMaxHitCut(double maxhitCut = 4.0)             { mMaxHitCut = maxhitCut; }
   void setCmnCorrection( bool doCmnCorrection = true )	 { mDoCmnCorrection = doCmnCorrection; }
   void setCmnCut(double cmnCut = 3.)                    { mCmnCut = cmnCut; }
   /// 0 - All data; 1 - non-ZS data; 2 - ZS data; 3 - ZS first data
   void setDataType(int nDataType = 0)                   { mDataType = nDataType; };
   void setDoEmbedding(Bool_t doIt)                      { mDoEmbedding = doIt; }
   UChar_t getDataType() {return mDataType;} // 0 - All data; 1 - non-ZS data; 2 - ZS data; 3 - ZS first data

   /// Occupancy gate for high-multiplicity AuAu. Above maxRawHits raw hits this
   /// maker drops its own collection and switches OFF, for that event only, the
   /// FST cluster/hit makers, the whole FTT chain, forward tracking and the FCS
   /// track match. The event is NOT aborted: FCS and everything else still run
   /// and the event is still written, just with no forward hits or tracks.
   ///
   /// n <= 0 (the default -1) means NO cut. Zero is deliberately not the "off"
   /// value -- it is a legitimate-looking threshold, and a mis-set 0 would
   /// silently kill FST for a whole production.
   ///
   /// The count is the real raw-hit count, tested after each DAQ element, so
   /// the unpacking stops as soon as the threshold is crossed and no more than
   /// one APV's worth of hits is ever built beyond it.
   void setMaxRawHits(int n = -1)                        { mMaxRawHits = n; }
   int  getMaxRawHits() const                            { return mMaxRawHits; }
   /// true if the current event was gated (for the makers downstream and QA)
   Bool_t isOccupancyGated() const                       { return mGated; }

protected:
   Bool_t mIsCaliMode;
   Bool_t mDoEmbedding;
   Bool_t mDoCmnCorrection;
   //control paramters
   double mMinHitCut, mMedHitCut, mMaxHitCut, mCmnCut, mChanMinRmsNoiseLevel, mChanMaxRmsNoiseLevel, mApvMaxCmNoiseLevel;
   Int_t mALLdata, mADCdata, mZSdata, mDefaultTimeBin, mCurrentTimeBinNum;
   Int_t mMinNumOfRawHits, mMaxNumOfRawHits;

   /// Main output container with either corrected or uncorrected ADC values per
   /// channel. Usually the values come from real data DAQ records but can be
   /// merged with simulated hits or entirely provided by StFstSlowSimMaker
   StFstCollection *mFstCollectionPtr;

   /// Pointer to an FST hit container with simulated ADC values usually provided
   /// by StFstSlowSimMaker
   StFstCollection *mFstCollectionSimuPtr;

   std::vector< std::vector< std::vector< float > > > mCmnVec;   ///< APV chip geom. index, common mode (CM) noise
   std::vector< std::vector< float > > mPedVec;   ///< Channel elec. index, pedestal
   std::vector< std::vector< float > > mTotRmsVec;   ///< Channel elec. index, Total RMS noise
   std::vector< std::vector< float > > mRanRmsVec;   ///< Channel elec. index, Random RMS noise
   std::vector< float > mGainVec;  ///< Channel elec. index, gain
   std::vector< int > mMappingVec; ///< Channel elec. index to geometry ID mapping
   std::vector< int > mConfigVec; ///< APV chip configuration status indexed by geom. id

private:

   int FillRawHitCollectionFromAPVData(unsigned char dataFlag, int ntimebin, int counterAdcPerEvent[][kFstNumTimeBins], double sumAdcPerEvent[][kFstNumTimeBins], int apvElecId,
      std::array< std::array<double, kFstNumTimeBins>, kFstNumApvChannels > &signalUnCorrected,
      std::array< std::array<double, kFstNumTimeBins>, kFstNumApvChannels > &signalCorrected,
      std::array< std::array<int, kFstNumTimeBins>, kFstNumApvChannels > &seedFlag,
      std::array<int, kFstNumApvChannels> &idTruth);

   int FillRawHitCollectionFromSimData();

   StEvent* mEvent;
   StFstEvtCollection* mFstEvtCollection;

   Int_t mDataType; ///<  0=all, 1=adc only, 2=zs only

   // ---- occupancy gate ----------------------------------------------------
   enum { kNGateMakers = 8 };
   Int_t    mMaxRawHits;              ///< <=0 means no cut
   Bool_t   mGated;                   ///< was THIS event gated
   Bool_t   mGateResolved;            ///< makers looked up yet
   StMaker *mGateMaker[kNGateMakers]; ///< resolved once in Init(), never owned
   Int_t    mNGateMaker;
   Long64_t mNEventsSeen, mNEventsGated;
   void     resolveGateMakers();              ///< Init(): find them, say what is missing
   void     setGateMakersActive(Bool_t on);   ///< flip the whole controlled set

   ClassDef(StFstRawHitMaker, 1);
};

#endif
