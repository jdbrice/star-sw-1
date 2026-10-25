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
   /// fstSensorOnWedge) to the hit's global position. With tpcOnGlobal out of
   /// StFstDb's product the matrix now agrees to the bit with the displacement
   /// AGML gives the TGeo nodes that become the GenFit planes, so a surveyed
   /// hit lands exactly on its own plane. DEFAULT OFF all the same, for two
   /// reasons that are about the data and the simulation, not the frame:
   ///   - the geometry starsim builds from is still ideal (measured: hits
   ///     placed at the pure nominal geometry match the GEANT hit to better
   ///     than 0.4 mm, where the survey would show 3 mm), so in MC this would
   ///     move the hit onto a misaligned plane and away from the particle that
   ///     made it. Fix the simulated geometry first;
   ///   - on data the 117M-hit residuals show no rigid 3 mm offset, and
   ///     FstWedgeAligner (StFwdTrackMaker/StFwdHitLoader.h) already applies a
   ///     measured per-wedge phi correction. Turning this on without taking
   ///     that off double-corrects, and every sTGC pentOnStation iteration
   ///     fitted against FST hits would need a refit.
   void SetApplySurvey(Bool_t v = kTRUE) { mApplySurvey = v; }
   Bool_t GetApplySurvey() const         { return mApplySurvey; }

protected:

   THashList *mSensorTransforms; ///< A list of TGeo transformations for each FST sensor
   Bool_t     mApplySurvey;      ///< see SetApplySurvey(); default false

   ClassDef(StFstHitMaker, 0);
};

#endif
