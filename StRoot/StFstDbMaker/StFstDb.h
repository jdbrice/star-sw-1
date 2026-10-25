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
   static const TGeoHMatrix *getHMatrixSensorOnGlobal(int wedge, int sensor);

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
   /// +1 front (+z) / -1 back, from kFstzFilp and kFstzDirct; 0 if out of range
   static int    wedgeOrientation(int disk, int wedge);
   static double sensorZOffset(int disk, int wedge, int sensor); ///< signed cm from the disk nominal
   static double sensorZ(int disk, int wedge, int sensor);       ///< diskNominalZ + sensorZOffset

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

   // ---- sensor-local <-> global ------------------------------------------
   // The same convenience StFttDb::localToGlobal()/globalToLocal() give for an
   // sTGC quadrant, so a caller does not have to fetch the matrix and remember
   // which of LocalToMaster/MasterToLocal it wanted. The matrix is the fully
   // composed tpc2Global * fst2Tpc * hss2Fst * wedge2Hss * sensor2Wedge, so
   // these move with the survey tables -- which is the point: a hit and the
   // GenFit plane it is fitted on have to come from the same place.
   //
   // wedge 1..36, sensor 0..2. Units are cm, as TGeo uses. False means the
   // index was out of range or the geometry has not been loaded.
   static bool sensorLocalToGlobal(int wedge, int sensor, const double l[3], double g[3]);
   static bool globalToSensorLocal(int wedge, int sensor, const double g[3], double l[3]);

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
