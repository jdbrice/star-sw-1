/* $Id: StFstDb.h */

#ifndef StFstDb_hh
#define StFstDb_hh

#include "StObject.h"
#include "THashList.h"
#include "TGeoMatrix.h"
#include "StEvent/StFstConsts.h"

class Survey_st;
class fstPedNoise_st;
class fstControl_st;
class fstGain_st;
class fstMapping_st;
class fstChipConfig_st;


/**
 * FST calibration/geometry database access proxy.
 *
 * relation within STAR frame
 * FstOnGlobal = Tpc2Magnet * Fst2Tpc *    Hss2Fst     * Wedge2Hss * Sensor2Wedge * PS

 * Naming of roatation matrices in this maker :
 * positionGlobal  = tpc2Global * fst2Tpc * hss2Fst * wedge2Hss * sensor2Wedge * positionOnSensor

 * numbering
 * Id  = 1000 + (wedge-1)*3 + sensor
 * 1<= wedge <= 36
 * 0<= sensor <= 2
 *
 *
 * \author Shenghui Zhang
 * \date Oct 2021
 */
class StFstDb : public StObject
{

public:
   StFstDb();
   THashList *getRotations() const                                                {return mgRotList; }
   /// The TPC's own alignment. Kept for reference and NOT part of the composed
   /// sensor-on-global matrix: AGML misaligns FSTM in CAVE with Geometry/fst/*
   /// only, so applying it would move the hits out of the frame of the GenFit
   /// planes built from those same TGeo nodes. See setGeoHMatrices().
   const TGeoHMatrix *getGeoHMatrixTpcOnGlobal() const                            {return mGeoHMatrixTpcOnGlobal; }
   const TGeoHMatrix *getGeoHMatrixFstOnTpc() const                               {return &mGeoHMatrixFstOnTpc; }
   /// hssOnFst is per half (halfId(wedge)) and fstWedgeOnHss per wedge, so both
   /// take an index now; the old no-argument form silently returned row 0.
   const TGeoHMatrix *getGeoHMatrixHssOnFst(Int_t half) const                     {return &mGeoHMatrixHssOnFst[half]; }
   const TGeoHMatrix *getGeoHMatrixWedgeOnHss(Int_t wedge) const                   {return &mGeoHMatrixWedgeOnHss[wedge - 1]; }
   const TGeoHMatrix *getGeoHMatrixSensorOnWedge(Int_t wedge, Int_t sensor) const {return &mGeoHMatrixSensorOnWedge[wedge - 1][sensor]; }

   /// GEANT g2t volume_id -> reco numbering. One place, so simulation and
   /// reconstruction cannot drift. All three return -1 if the id is not FST.
   /// See StFstDb.cxx for the measured GEANT-wedge -> phi-sector map.
   Int_t getModuleFromG2TVolumeId(Int_t volumeId) const;   ///< module 1-12 within the disk
   Int_t getWedgeFromG2TVolumeId(Int_t volumeId) const;    ///< reco wedge 1-36
   Int_t getSensorFromG2TVolumeId(Int_t volumeId) const;   ///< 0 inner, 1/2 outer halves
   Int_t getDiskFromG2TVolumeId(Int_t volumeId) const;     ///< disk 1-3 (the id digit is 4-6)
   /// The survey DELTA for a sensor, in the global frame -- NOT a placement,
   /// despite what its old name (getHMatrixSensorOnGlobal, to 2026-10-30) said.
   /// Identity rotation plus a few mm, the same for every wedge of a half. For
   /// where a sensor IS, see sensorPlacement() below.
   static const TGeoHMatrix *getSurveyDeltaSensorOnGlobal(int wedge, int sensor);

   // ---- nominal disk z, and where the silicon actually sits ------------------
   // A wedge is mounted on the FRONT or the BACK of its disk, alternating with
   // wedge parity, and the silicon is offset from the disk's nominal z by a
   // different amount for the inner and the outer sensors. FstmGeoSM.xml:
   //
   //   zdisk              = {151.750, 165.248, 178.781}
   //   innersilicon_z     = ibase2 + HybridThickness + 0.03/2        = 1.7419
   //   outersilicon_z     = sRingThickness - HybridThickness - 0.03/2 = 0.3469
   //   modulus2           = {0,1,0,1,0,1,0,1,0,1,0,1}
   //
   // and the back wedges (which also carry <Rotation alphax="180">) are the ones
   // with modulus2 == 0 on disks 1 and 3, but modulus2 == 1 on disk 2 -- the
   // parity is INVERTED on disk 2, which is why a single formula gets one disk
   // wrong. Measured per (disk, wedge, sensor) with script/zOffsetMap.C: inner
   // -+1.7419, outer -+0.3469, rms 0.003 cm, and the disk-2 flip confirmed on
   // all 12 wedges.
   //
   // StFstHitMaker used to assign the bare nominal disk z to every hit, i.e. the
   // front/back mean, leaving every hit 1.742 cm (inner) or 0.347 cm (outer) from
   // the silicon it was measured on. 0 of 6812 hits in a 2000-event closure run
   // landed within 1.5 mm of the right z.
   static double diskNominalZ(int disk);               ///< disk 1-3 -> cm, 0 if out of range
   /// +1 front (+z) / -1 back, from kFstzFilp and kFstzDirct; 0 if out of range.
   /// Which way round depends on fstmConfig(), below.
   static int    wedgeOrientation(int disk, int wedge);

   // ---- which way the wedges are mounted ----------------------------------
   // Same meaning and the same two values as FSTG_fstmconfig in FstmGeoM.xml /
   // FstmGeoSM.xml, so one word covers the XML and the reco:
   //
   //   kFstmAsInstalled (2)  DEFAULT. The mounting parity of the real detector,
   //                         seen on it and on the installation photographs;
   //                         Flemming's mirror of 2026-08-10, AGML wholewedge
   //                         version 2.0, modulus2 {1,0,1,0,...}. This is what
   //                         y2022c / ft2022c build, and what DATA needs.
   //   kFstmPreMirror   (1)  The parity every geometry carried before the
   //                         mirror: AGML wholewedge version 1.0, modulus2
   //                         {0,1,0,1,...}. Needed to reconstruct anything
   //                         generated on FSTMv1, FSTMv1s, FSTMv1m or FSTMv1sm,
   //                         i.e. on y2022..y2025, dev2021, dev2022, dev2022m,
   //                         dev2022sm or dev2026 -- every sample that existed
   //                         before 2026-10-25.
   //
   // The two differ by 2 x 1.7419 cm on an inner sensor and 2 x 0.3469 cm on an
   // outer one, at every one of the 36 wedges, and by nothing in x or y.
   // Reconstructing a sample with the wrong one puts the hit on the wrong face
   // of its disk. The default is deliberately the real detector and not the
   // old samples, exactly as Flemming set the module default in the XML.
   enum EFstmConfig { kFstmPreMirror = 1, kFstmAsInstalled = 2 };
   static void   setFstmConfig(Int_t config);   ///< kFstmPreMirror or kFstmAsInstalled
   static Int_t  fstmConfig();                  ///< default kFstmAsInstalled

   // ---- the flip changes z, and must NOT change phi ------------------------
   // Confirmed on the detector with Flemming, 2026-10-29: the mirrored mounting
   // is real, and so is the z change that comes with it. What does NOT change is
   // the strip-to-phi relation: a given phi strip sits at the same global azimuth
   // as it always did, which is why production's phi convention survives the
   // 117M-hit residuals (per-wedge corrections of 1-3 mrad, where a mirrored phi
   // would be 15 deg).
   //
   // AGML cannot express that split. It implements the flip as a 180 deg rotation
   // about the wedge's local x axis, which moves the silicon to the other face
   // (right) AND mirrors the sensor's local phi about the wedge bisector (wrong).
   // So a geometry-driven decode -- anything that turns a strip into a position
   // using the FTUS node matrix and shape, as the planar track fit does -- picks
   // up a phi mirror that the detector does not have, up to 15 deg for an outer
   // sensor and 30 deg across the inner one. Measured in the live chain: the
   // stored hit and the geometry decode disagree by up to 25 deg in phi while
   // agreeing in R to the micron.
   //
   // Hence the two senses of this one switch:
   //
   //   kFstmAsInstalled (2, DEFAULT)  z  = the installed face (wedgeOrientation)
   //                                  phi= FLIPPED BACK, i.e. the same phi as
   //                                       before the XML wedge flip
   //   kFstmPreMirror   (1)           z  and phi both as they were before the
   //                                       flip -- for reproducing samples built
   //                                       on the pre-mirror geometry
   //
   // stripPhi() below is deliberately INDEPENDENT of this setting: it already
   // produces the as-installed phi in both configurations, because the detector's
   // strip-to-phi relation never changed. The flip is only ever applied to phi
   // taken FROM THE GEOMETRY, through the two helpers here.
   static bool   mirrorAgmlLocalPhi();                  ///< true for kFstmAsInstalled
   /// Mirror a sensor-shape local azimuth about the WEDGE bisector (local 360 deg),
   /// which is what undoes AGML's 180 deg x-rotation. Mirroring about each sensor's
   /// own bisector would be wrong: the wedge mirror also exchanges the two outer
   /// halves, and this reproduces that. Returns phiLocalDeg unchanged when
   /// mirrorAgmlLocalPhi() is false, so a caller can apply it unconditionally.
   static double agmlLocalPhi(double phiLocalDeg);
   static double sensorZOffset(int disk, int wedge, int sensor); ///< signed cm from the disk nominal
   static double sensorZ(int disk, int wedge, int sensor);       ///< diskNominalZ + sensorZOffset

   // ---- where a sensor IS: the placement matrix, the way AGML builds it -----
   // CAREFUL: getSurveyDeltaSensorOnGlobal() above is NOT this. It was called
   // getHMatrixSensorOnGlobal() until 2026-10-30, and under that name it was
   // read as a placement and put all 108 FST GenFit planes at the origin.
   // The name came from StIstDb/StSstDbMaker, where the survey tables
   // ARE full placements. The FST tables are misalignment DELTAS, so that matrix
   // comes out as identity rotation plus a few mm of translation, carrying no
   // wedge azimuth, no disk z and no face offset. Measured at sdt20211112.000004
   // it is (+0.0713, -0.3010, 0.0000) cm for every wedge of half 0 and
   // (+0.2263, -0.3330, 0.0000) for every wedge of half 1 -- identical across a
   // half, which is what proves it cannot be a placement.
   //
   // This builds the real thing, as FstmGeoSM.xml places FTUS in FSTM:
   //
   //     <Placement block="FTUS" in="FSTM" x="0" y="0" z="..._silicon_motherz">
   //       <Rotation alphax="wedge_xrot"/>     0, or 180 when the wedge is
   //       <Rotation alphaz="wedge_zrot"/>     mounted on the far face
   //
   // so, for a local vector, global = Rz(azimuth) * Rx(flip) * local + T(0,0,z),
   // with z = sensorZ() and the flip taken from wedgeOrientation(). FSTM itself
   // sits at x=y=0 in CAVE, and the XML's z is relative to zcenter, which is why
   // sensorZ() -- diskNominalZ + sensorZOffset -- is already the global z.
   //
   // All three sensors of a wedge share one rotation. The two outer sensors
   // differ only in the TUBS phi extent of the shape, which is not part of the
   // placement, so there is no per-sensor azimuth convention to choose here.
   //
   // FTUS is a tube segment placed on the beam axis, so its LOCAL ORIGIN IS ON
   // THAT AXIS, not at the sensor centroid: local r equals global r, and local
   // x,y are the global x,y turned by -azimuth. Verified against the simulation's
   // own AGML local coordinates, which average r = 13.1 cm over the 6-27 cm span.
   //
   // applySurvey left-multiplies the survey delta, which is exactly what
   // getSurveyDeltaSensorOnGlobal() already holds, composed from the same tables
   // StFstHitMaker applies. Plane and measurement therefore come from one
   // source. NOTE it includes fstOnTpc, which the XML's <Misalign> chain does
   // NOT -- the XML applies only fstSensorOnWedge, fstWedgeOnHss and hssOnFst --
   // so ours and AGML's differ by fstOnTpc, (+0.1488, -0.3170, 0.0000) cm at
   // sdt20211112. That is a real inconsistency between the database and the
   // geometry, not a bug here; it is kept because the hit maker applies it, and
   // a plane that dropped it would sit in a different frame from its own hits.
   static double      wedgeAzimuth(int wedge);   ///< deg; FstmGeoSM wedge_zrot
   static TGeoHMatrix sensorPlacement(int disk, int wedge, int sensor,
                                      bool applySurvey = true);

   // ---- strip -> position, the ONE implementation ---------------------------
   // Moved here verbatim from StFstHitMaker so that the producer, the slow
   // simulation and anything that has to rebuild a hit from its strips all read
   // the convention from one place. Before this, StFstHitMaker held it,
   // StFstSlowSimMaker held the same algebra inverted, and StFwdHitLoader
   // patched its output -- three copies of kFstzFilp x kFstzDirct and two of
   // kFstStripGapPhi, which is why a change to one of them could not be made
   // safely. TrackFitter.h keeps its own GEOMETRY-driven version on purpose
   // (the fit must work from the AGML shape, and that one follows the wedge
   // mirror automatically); script/checkStripMap.C is what keeps the two models
   // from drifting apart.
   //
   // phiStrip is a double on purpose: clustering makes it a charge-weighted
   // mean, and truncating it to int costs a third of a pitch.
   static double stripR  (int rStrip);                                    ///< cm, 0 out of range
   static double stripPhi(int disk, int wedge, int sensor,
                          int rStrip, double phiStrip);                   ///< rad, 0 out of range
   /// Nominal-frame CARTESIAN position of a strip, cm. Returns false (and leaves
   /// xyz untouched) if any index is out of range. z is sensorZ(), so it follows
   /// fstmConfig(); x,y are r and phi of the strip in the ideal detector, which
   /// is already a coordinate in the frame the survey matrices call "local".
   static bool   stripPosition(int disk, int wedge, int sensor,
                               int rStrip, double phiStrip, double xyz[3]);

   // ---- survey-table index convention ------------------------------------
   // fstSensorOnWedge rows are addressed by the Id COLUMN, not by their
   // position in the table: the reader below takes Id and derives (wedge,
   // sensor) from it, so row order is irrelevant and a row carrying the wrong
   // Id lands on the wrong sensor however it is ordered. (The FTT is the other
   // way round -- its rows are positional, see StFttDb::surveyRow.)
   //
   // Anyone WRITING a table needs sensorSurveyId(); the reverse is what the
   // reader uses, and what a check macro needs to say which sensor a row is.
   // wedge is 1..kFstNumWedges and sensor 0..kFstNumSensorsPerWedge-1, the
   // human-friendly numbering the rest of this class uses.
   static int  sensorSurveyId(int wedge, int sensor);
   static bool decodeSensorSurveyId(int id, int &wedge, int &sensor);

   /// The half of the support structure a wedge is mounted on: the row of
   /// Geometry/fst/hssOnFst that applies to it. 0 for the x>0 half, 1 for x<0;
   /// -1 if the wedge index is out of range. See StFstDb.cxx for why it is
   /// derived from the wedge azimuth rather than from a wedge-number cut.
   static int  halfId(int wedge);

   // ---- nominal global <-> surveyed global ---------------------------------
   // Apply the survey delta to a point that is ALREADY global at its nominal
   // position, and take it back, so a caller does not have to fetch the matrix
   // and remember which of LocalToMaster/MasterToLocal it wanted. This is the
   // step StFstHitMaker takes when SetApplySurvey() is on.
   //
   // These are NOT sensor-local <-> global, which is what their names said
   // until 2026-10-30 (sensorLocalToGlobal/globalToSensorLocal, documented here
   // as the fully composed tpc2Global * fst2Tpc * hss2Fst * wedge2Hss *
   // sensor2Wedge). The matrix is fstOnTpc * hssOnFst[half] *
   // fstWedgeOnHss[wedge] * fstSensorOnWedge -- a misalignment delta, with no
   // wedge azimuth, no disk z and no mounting flip in it, and tpcOnGlobal is
   // not in the product either. Feed one a sensor-local point and it silently
   // returns a wrong answer, which is exactly the trap the old name set. For
   // local <-> global use sensorPlacement() with TGeo's own
   // LocalToMaster()/MasterToLocal().
   //
   // wedge 1..36, sensor 0..2. Units are cm, as TGeo uses. False means the
   // index was out of range or the geometry has not been loaded.
   static bool sensorNominalToSurveyed(int wedge, int sensor, const double nominal[3], double surveyed[3]);
   static bool surveyedToSensorNominal(int wedge, int sensor, const double surveyed[3], double nominal[3]);

   const fstPedNoise_st *getPedNoise() const 		{return mFstPedNoise;}
   const fstGain_st *getGain() const     		{return mFstGain;    }
   const fstMapping_st *getMapping() const  		{return mFstMapping; }
   const fstControl_st *getControl() const  		{return mFstControl; }
   const fstChipConfig_st *getChipStatus() const 	{return mFstChipStatus; }
      
   /// nRows, when given, is the row count of each of the four survey tables in
   /// the same order, so a short table is caught instead of read past the end.
   Int_t setGeoHMatrices(Survey_st **tables, const Int_t *nRows = 0);
   void setPedNoise(fstPedNoise_st *pedNoise) 	    { mFstPedNoise   = pedNoise;   }
   void setGain(fstGain_st *gain)		    { mFstGain       = gain;       }
   void setMapping(fstMapping_st *mapping)    	    { mFstMapping    = mapping;    }
   void setControl(fstControl_st *control)    	    { mFstControl    = control;    }
   void setChipStatus(fstChipConfig_st *chipStatus) { mFstChipStatus = chipStatus; }
   
   virtual void Print(Option_t *opt = "") const;

private:
   Int_t moduleFromSector(Int_t sector) const;  ///< phi/30deg sector -> module 1-12, from kFstphiStart/Stop

   enum { kFstNumHss = 2 };     ///< halves of the support structure = rows of hssOnFst

   static THashList 	*mgRotList; ///< A list of TGeoHMatrix transormations for each FST sensor
   static Int_t          mgFstmConfig; ///< see setFstmConfig(); default kFstmAsInstalled
   TGeoHMatrix *mGeoHMatrixTpcOnGlobal;
   TGeoHMatrix mGeoHMatrixFstOnTpc;
   /// hssOnFst has one row per HALF and fstWedgeOnHss one row per WEDGE, so both
   /// are arrays. They used to be single matrices filled from row 0 and applied
   /// to all 36 wedges -- see the comment on halfId() in StFstDb.cxx.
   TGeoHMatrix mGeoHMatrixHssOnFst[kFstNumHss];
   TGeoHMatrix mGeoHMatrixWedgeOnHss[kFstNumWedges];
   TGeoHMatrix mGeoHMatrixSensorOnWedge[kFstNumWedges][kFstNumSensorsPerWedge];

   fstPedNoise_st 	*mFstPedNoise;
   fstGain_st 		*mFstGain;
   fstMapping_st 	*mFstMapping;
   fstControl_st 	*mFstControl;
   fstChipConfig_st 	*mFstChipStatus;
   
   ClassDef(StFstDb, 1)
};

#endif
