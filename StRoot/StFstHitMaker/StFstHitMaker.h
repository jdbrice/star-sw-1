#ifndef StFstHitMaker_hh
#define StFstHitMaker_hh

#include "StMaker.h"

class THashList;


/**
 * Calculates hit global position, and writes FST hits to StFstHitCollection.
 *
 * \author: Shenghui Zhang
 * \date Oct. 2021
 */
class StFstHitMaker : public StMaker
{
public:

   StFstHitMaker( const char *name = "fst_hit" );
   Int_t InitRun(Int_t runnumber);
   Int_t Make();

   /// Apply the FST survey (Geometry/fst/fstOnTpc, hssOnFst, fstWedgeOnHss,
   /// fstSensorOnWedge) to the hit's global position. DEFAULT OFF -- see the
   /// long comment in Make() for the two reasons, both of which have to be
   /// settled before this can become the default:
   ///   - the AGML geometry every production tag uses carries no <Misalign>
   ///     tags, so the planes a hit is fitted to are ideal; moving the hit
   ///     alone would put it 3 mm off its own plane;
   ///   - the survey's ~3 mm y is quoted in the IP frame, and the STAR frame
   ///     the TPC vertex comes in already carries a similar ~2 mm y. Whether
   ///     those are the same offset is an open question, so applying the
   ///     table may double-count it.
   void SetApplySurvey(Bool_t v = kTRUE) { mApplySurvey = v; }
   Bool_t GetApplySurvey() const         { return mApplySurvey; }

protected:

   THashList *mSensorTransforms; ///< A list of TGeo transformations for each FST sensor
   Bool_t     mApplySurvey;      ///< see SetApplySurvey(); default false

   ClassDef(StFstHitMaker, 0);
};

#endif
