/***************************************************************************
 *
 * StFstSlowSimMaker.h
 *
 * Slow (raw-hit level) simulation maker for the FST.
 *
 * Adapted from StRoot/StIstSimMaker/StIstSlowSimMaker.{h,cxx} (L. Kosarzewski,
 * 2014), which is the maker StFstRawHitMaker was clearly written to pair with:
 * StFstRawHitMaker already consumes a "fstRawAdcSimu" dataset and already
 * implements embedding, but nothing in the tree ever produced that dataset.
 *
 * Chain position:
 *   g2t_fsi_hit -> StFstSlowSimMaker -> "fstRawAdcSimu" (StFstCollection)
 *               -> StFstRawHitMaker  -> StFstClusterMaker -> StFstHitMaker
 *
 * so MC runs through exactly the same clustering and hit-position code as real
 * data, instead of StFstFastSimMaker's parallel reimplementation.
 *
 ****************************************************************************/

#ifndef STAR_StFstSlowSimMaker
#define STAR_StFstSlowSimMaker

#include "StMaker.h"
#include <vector>

class StFstDb;
class StFstCollection;
class TF1;
class TRandom3;

/**
 * Turns GEANT FST hits (g2t_fsi_hit) into StFstRawHits and publishes them as
 * the "fstRawAdcSimu" dataset for StFstRawHitMaker to pick up.
 *
 * SCOPE OF THIS VERSION (v0, minimal):
 *  - one raw hit per GEANT hit, i.e. every cluster is single-strip by
 *    construction. Real data is 74.6% single-strip with mean size 1.39
 *    (measured on run 23081008), so cluster-size distributions will NOT match until
 *    a charge-sharing / crosstalk model is added. Fine for geometry closure,
 *    wrong for anything about clustering efficiency or resolution.
 *  - charge is a flat dE -> ADC scale with a Landau time-bin shape, using the
 *    IST numbers. These need FST-specific calibration.
 *  - no noise, no dead channels, no pedestal.
 *
 * KNOWN LIMITATION, read before trusting a closure test.
 * getStripsFromConsts() inverts the SAME kFstphiStart/kFstrStart formula that
 * StFstHitMaker uses in the forward direction, so the strip <-> position part
 * of the round trip is a tautology and cannot validate the strip formula
 * itself. What it DOES validate, and what has actually been biting us, is
 * everything around it: wedge numbering, disk assignment, the sensor/phiStrip
 * pairing, the elecId<->geoId mapping round trip, the outer-sensor gap sign,
 * and the whole clustering chain. A genuinely independent strip lookup has to
 * come from TGeo (find the FTUS node, MasterToLocal, divide by pitch) -- that is
 * the next step, hooked here as setStripSource(kStripFromTGeo) but not yet
 * implemented.
 */
class StFstSlowSimMaker : public StMaker
{
public:

   enum StripSource { kStripFromConsts = 0, kStripFromTGeo = 1 };

   StFstSlowSimMaker(const char *name = "fstSlowSim");
   virtual ~StFstSlowSimMaker();

   Int_t  Init();
   Int_t  InitRun(Int_t runNumber);
   Int_t  Make();
   void   Clear(Option_t *opts = "");

   /// Flat single-hit efficiency, applied per GEANT hit. Default 1.0.
   void setHitEff(Float_t eff)          { mHitEff = eff; }
   /// Where the (rStrip, phiStrip) lookup comes from. kStripFromTGeo is not
   /// implemented yet and falls back to kStripFromConsts with a warning.
   void setStripSource(StripSource s)   { mStripSource = s; }
   /// dE -> ADC scale factor, ADC counts per GeV deposited.
   void setAdcPerGev(Double_t f)        { mAdcPerGev = f; }

   virtual const char *GetCVS() const
   {
      static const char cvs[] = "Tag $Name: $ built " __DATE__ " " __TIME__;
      return cvs;
   }

protected:

   StFstDb         *mFstDb;
   StFstCollection *mFstCollectionPtr;

   TF1      *fAdc;              ///< Landau shape used for the time-bin profile
   TRandom3 *mRndGen;

   UChar_t   mDefaultTimeBin;
   UChar_t   mCurrentTimeBinNum;

   /// geoId -> elecId, the inverse of Calibrations/fst/fstMapping.
   /// Both ids are 0-based for the FST (the IST table is 1-based -- do not
   /// copy the "-1" from StIstSlowSimMaker).
   std::vector<Int_t> mMappingGeomVec;

   /// sector index (0-11, from phi/30deg) -> moduleIdx (1-12), built in InitRun
   /// from kFstphiStart/kFstphiStop so the GEANT volume_id wedge numbering is
   /// never relied on.
   Int_t     mSectorToModule[12];

   Float_t      mHitEff;
   Double_t     mAdcPerGev;
   StripSource  mStripSource;

   Long64_t  mNHitsSeen, mNHitsUsed, mNHitsBadGeo;

private:

   /// global (x,y,z) of a GEANT hit -> (wedge 1-36, sensor 0-2, rStrip 0-7,
   /// phiStrip 0-127). Returns false if the hit is outside the active area.
   Bool_t getStripsFromConsts(Double_t x, Double_t y, Double_t z, Int_t disk,
                              Int_t &wedge, Int_t &sensor,
                              Int_t &rStrip, Int_t &phiStrip) const;

   /// fill one StFstRawHit for this GEANT hit
   void generateRawHit(Int_t wedge, Int_t sensor, Int_t rStrip, Int_t phiStrip,
                       Double_t de, Int_t idTruth);

   ClassDef(StFstSlowSimMaker, 0)
};

#endif
