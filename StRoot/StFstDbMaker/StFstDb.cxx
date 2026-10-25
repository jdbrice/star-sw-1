/* $Id: StFstDb.cxx,v */

#include <assert.h>
#include "StFstDb.h"
#include "StMessMgr.h"
#include "StTpcDb/StTpcDb.h"
#include "St_db_Maker/St_db_Maker.h"

#include "tables/St_Survey_Table.h"
#include "tables/St_fstPedNoise_Table.h"
#include "tables/St_fstGain_Table.h"
#include "tables/St_fstMapping_Table.h"
#include "tables/St_fstControl_Table.h"
#include "tables/St_fstChipConfig_Table.h"

THashList *StFstDb::mgRotList = 0;

ClassImp(StFstDb)


/**
 * \author Shenghui Zhang
 * \date Oct 2021
 */
StFstDb::StFstDb() : StObject()
{
  mGeoHMatrixTpcOnGlobal = NULL;
  mFstPedNoise           = NULL;
  mFstGain               = NULL;
  mFstMapping            = NULL;
  mFstControl            = NULL;
  mFstChipStatus         = NULL;
}
//_____________________________________________________________________________
// GEANT g2t volume_id -> reco numbering.
//
//   volume_id = diskDigit*1000 + wedgeDigit*10 + sensorDigit
//
// with diskDigit 4, 5, 6 -- matching FSTD_4/5/6 in the TGeo tree
// (HALL/CAVE_1/FSTM_1/FSTD_n/FSTW_m/FTUS_k), NOT 1, 2, 3 -- wedgeDigit the GEANT
// wedge 1-12 and sensorDigit 1-3.
//
// The GEANT wedge runs the OTHER WAY ROUND from the phi sector (sector =
// floor(phi/30deg)), with a 3-sector offset, identically on all three disks:
//
//   gWedge   1  2  3  4  5  6  7  8  9 10 11 12
//   sector   2  1  0 11 10  9  8  7  6  5  4  3      = (3 - gWedge) mod 12
//
// Measured on a 2pi electron-gun sample (script/fstWedgeMap.C, 1566 hits,
// dev2022sm / sdt20211112.000004): each of the 36 (disk, gWedge) cells puts
// 82-100% of its hits in that sector and the three disks agree exactly. The
// remainder falls in the neighbouring sector, which is expected -- the sector is
// the IDEAL 30 degree bin while the wedges are surveyed -- so a caller that
// cares about wedge edges should trust this decode over a phi guess.
//
// kFstphiStart/kFstphiStop then map that sector back to a module, and THAT map is
// a reflection too, so the two cancel: the GEANT wedge digit currently equals the
// reco module index. It is still written as the composition rather than as
// "return gWedge", so that if the phi constants are ever renumbered this follows
// them instead of silently becoming wrong.
Int_t StFstDb::getModuleFromG2TVolumeId(Int_t volumeId) const
{
   const Int_t diskDigit  = volumeId / 1000;
   const Int_t gWedge     = (volumeId % 1000) / 10;

   if (diskDigit < 4 || diskDigit > 6)             return -1;
   if (gWedge < 1 || gWedge > kFstNumWedgePerDisk) return -1;

   const Int_t sector = (3 - gWedge + kFstNumWedgePerDisk) % kFstNumWedgePerDisk;

   return moduleFromSector(sector);
}
//_____________________________________________________________________________
Int_t StFstDb::getDiskFromG2TVolumeId(Int_t volumeId) const
{
   const Int_t diskDigit = volumeId / 1000;

   if (diskDigit < 4 || diskDigit > 6) return -1;

   return diskDigit - 3;   // -> 1, 2, 3, the numbering StFstConsts uses
}
//_____________________________________________________________________________
double StFstDb::diskNominalZ(int disk)
{
   // FstmGeoSM.xml, <var name="zdisk" value="{151.750, 165.248, 178.781}" />.
   // These are the same three numbers StFstHitMaker.cxx used to hardcode; the
   // nominal z is unavoidably duplicated between the XML and here, so keep them
   // in one place on this side at least.
   static const double kZDisk[kFstNumDisk] = { 151.750, 165.248, 178.781 };

   if (disk < 1 || disk > kFstNumDisk) return 0.0;

   return kZDisk[disk - 1];
}
//_____________________________________________________________________________
int StFstDb::wedgeOrientation(int disk, int wedge)
{
   if (disk < 1 || disk > kFstNumDisk) return 0;
   if (wedge < 1 || wedge > kFstNumWedges) return 0;

   const int iwedgem = (wedge - 1) % kFstNumWedgePerDisk;   // 0-based within the disk

   // The two arrays this needs already exist in StFstConsts.h and are the ones
   // StFstHitMaker uses for the phi logic, so take the front/back sense from them
   // rather than restating it -- then the z and the phi cannot drift apart:
   //     kFstzDirct[iwedgem] = {1,-1,1,-1,...}  the alternation with wedge parity
   //     kFstzFilp[disk-1]   = {1,-1,1}         the inversion on disk 2
   // FstmGeoSM.xml is the authority behind them: a wedge goes to the BACK (-z,
   // plus <Rotation alphax="180">) when modulus2 == 0 on disks 1 and 3 but when
   // modulus2 == 1 on disk 2.
   //     <If expr="idisk.ne.1"> <If expr="...modulus2....eq.0"> ... </If> </If>
   //     <If expr="idisk.eq.1"> <If expr="...modulus2....eq.1"> ... </If> </If>
   // The leading minus is the measured sign: disk 1 wedge 1 has both array
   // entries +1 and sits at -1.7408 cm. Verified against all 36 wedges.
   return -kFstzFilp[disk - 1] * kFstzDirct[iwedgem];
}
//_____________________________________________________________________________
double StFstDb::sensorZOffset(int disk, int wedge, int sensor)
{
   // FstmGeoSM.xml:
   //   innersilicon_z = ibase2 + HybridThickness + 0.03/2
   //                  = (sRingThickness + oBaseThickness + ioTubeLength
   //                     + iBaseThicknessAB + iBaseThickness) + 0.035 + 0.015
   //                  = (0.3969 + 0.20 + 0.835 + 0.06 + 0.20) + 0.05 = 1.7419
   //   outersilicon_z = sRingThickness - HybridThickness - 0.03/2
   //                  = 0.3969 - 0.035 - 0.015                        = 0.3469
   // Sensor 0 is the inner one; 1 and 2 are the two halves of the outer, which
   // share a z.
   static const double kInnerZ = 1.7419;   // cm
   static const double kOuterZ = 0.3469;   // cm

   if (sensor < 0 || sensor >= kFstNumSensorsPerWedge) return 0.0;

   const int orientation = wedgeOrientation(disk, wedge);
   if (0 == orientation) return 0.0;

   return orientation * ((0 == sensor) ? kInnerZ : kOuterZ);
}
//_____________________________________________________________________________
double StFstDb::sensorZ(int disk, int wedge, int sensor)
{
   const double z0 = diskNominalZ(disk);
   if (0.0 == z0) return 0.0;

   return z0 + sensorZOffset(disk, wedge, sensor);
}
//_____________________________________________________________________________
Int_t StFstDb::getWedgeFromG2TVolumeId(Int_t volumeId) const
{
   const Int_t disk   = getDiskFromG2TVolumeId(volumeId);
   const Int_t module = getModuleFromG2TVolumeId(volumeId);

   if (disk < 1 || module < 1) return -1;

   return module + (disk - 1) * kFstNumWedgePerDisk;   // reco wedge 1-36
}
//_____________________________________________________________________________
// The sensor digit is 1 for the inner sensor and 2, 3 for the two halves of the
// outer one, which StFstConsts counts as 0, 1, 2.
Int_t StFstDb::getSensorFromG2TVolumeId(Int_t volumeId) const
{
   const Int_t diskDigit = volumeId / 1000;

   if (diskDigit < 4 || diskDigit > 6) return -1;

   const Int_t sensorDigit = volumeId % 10;

   if (sensorDigit < 1 || sensorDigit > kFstNumSensorsPerWedge) return -1;

   return sensorDigit - 1;
}
//_____________________________________________________________________________
Int_t StFstDb::moduleFromSector(Int_t sector) const
{
   if (sector < 0 || sector >= kFstNumWedgePerDisk) return -1;

   // kFstphiStart/kFstphiStop are in units of pi/6, one entry per module, and a
   // module spans the sector given by the lower of the two.
   for (Int_t m = 0; m < kFstNumWedgePerDisk; m++) {
      const Int_t a = (Int_t) kFstphiStart[m];
      const Int_t b = (Int_t) kFstphiStop[m];

      if (((a < b) ? a : b) % kFstNumWedgePerDisk == sector) return m + 1;
   }

   return -1;
}
//_____________________________________________________________________________
Int_t StFstDb::setGeoHMatrices(Survey_st **tables)
{
   //using namespace StFstConsts;

   SafeDelete(mgRotList);
   mgRotList = new THashList(kFstNumWedges * kFstNumSensorsPerWedge, 0);
   mgRotList->SetOwner(kFALSE);

   //get TPC positionement relative to STAR
   if (gStTpcDb) {
      mGeoHMatrixTpcOnGlobal = (TGeoHMatrix *)&gStTpcDb->Tpc2GlobalMatrix();
   }
   else {
      if (mGeoHMatrixTpcOnGlobal) delete mGeoHMatrixTpcOnGlobal;

      mGeoHMatrixTpcOnGlobal = new TGeoHMatrix("tpcOnGlobal");
      LOG_WARN << "No gStTpcDb, use null transformation for tpc on global" << endm;
   }

   //obtain FST geomery tables
   Survey_st *fstOnTpc        = tables[0];
   Survey_st *hssOnFst        = tables[1];
   Survey_st *wedgeOnHss      = tables[2];
   Survey_st *sensorsOnWedges = tables[3];

   mGeoHMatrixFstOnTpc.SetName("fstOnTpc");
   mGeoHMatrixFstOnTpc.SetRotation(&fstOnTpc->r00);
   mGeoHMatrixFstOnTpc.SetTranslation(&fstOnTpc->t0);

   mGeoHMatrixHssOnFst.SetName("hssOnFst");
   mGeoHMatrixHssOnFst.SetRotation(&hssOnFst->r00);
   mGeoHMatrixHssOnFst.SetTranslation(&hssOnFst->t0);

   mGeoHMatrixWedgeOnHss.SetName("wedgeOnHss");
   mGeoHMatrixWedgeOnHss.SetRotation(&wedgeOnHss->r00);
   mGeoHMatrixWedgeOnHss.SetTranslation(&wedgeOnHss->t0);
   for (int i = 0; i < kFstNumSensors; i++, sensorsOnWedges++) {
      int id = sensorsOnWedges->Id;
      TGeoHMatrix *comb = (TGeoHMatrix *) mgRotList->FindObject(Form("R%04i", id));

      if (comb) continue;

      comb = new TGeoHMatrix(Form("R%04i", id));
      int wedge = 0, sensor = 0;
      if (!decodeSensorSurveyId(id, wedge, sensor)) continue;  // 1<=wedge<=36, 0<=sensor<=2

      if (wedge <= 0 || wedge > kFstNumWedges) {
         LOG_WARN << "Wedge ID is out of range (1 - 36)!" << endm;
         continue;
      }

      if (sensor < 0 || sensor >= kFstNumSensorsPerWedge) {
         LOG_WARN << "Sensor ID is out of range (0 - 2)!" << endm;
         continue;
      }

      //setting rotation/translation for sensor geometry matrix
      mGeoHMatrixSensorOnWedge[wedge - 1][sensor].SetName(Form("sensorOnWedge%4i%4i", wedge, sensor));
      mGeoHMatrixSensorOnWedge[wedge - 1][sensor].SetRotation(&sensorsOnWedges->r00);
      mGeoHMatrixSensorOnWedge[wedge - 1][sensor].SetTranslation(&sensorsOnWedges->t0);

      TGeoHMatrix *sensorLocal = (TGeoHMatrix *) mgRotList->FindObject(Form("sensorLocal%04i", id));

      if (!sensorLocal) {
         sensorLocal = new  TGeoHMatrix(Form("sensorLocal%04i", id));
         sensorLocal->SetRotation(mGeoHMatrixSensorOnWedge[wedge - 1][sensor].GetRotationMatrix());
         sensorLocal->SetTranslation(mGeoHMatrixSensorOnWedge[wedge - 1][sensor].GetTranslation());
         mgRotList->Add(sensorLocal);
      }

      //calculate sensor global position
      TGeoHMatrix sensorGlobal = (*mGeoHMatrixTpcOnGlobal) * mGeoHMatrixFstOnTpc * mGeoHMatrixHssOnFst * mGeoHMatrixWedgeOnHss * mGeoHMatrixSensorOnWedge[wedge - 1][sensor];

      comb->SetRotation(sensorGlobal.GetRotationMatrix());
      comb->SetTranslation(sensorGlobal.GetTranslation());
      mgRotList->Add(comb);
   }

   return kStOk;
}


/**
 * Returns TGeoHMatrix with complete set of transformations from the sensor
 * local coordinate system to the global one. The wedge and sensor id-s are
 * expected to follow the human friendly numbering scheme.
 */
const TGeoHMatrix *StFstDb::getHMatrixSensorOnGlobal(int wedge, int sensor)
{
   if (wedge < 1 || wedge > kFstNumWedges || sensor < 0 || sensor >= kFstNumSensorsPerWedge)
      return 0;

   int id = sensorSurveyId(wedge, sensor);
   return mgRotList ? (const TGeoHMatrix *) mgRotList->FindObject(Form("R%04i", id)) : 0;
}


/**
 * The fstSensorOnWedge Id convention, in one place and both directions.
 * Id = 1000 + (wedge-1)*kFstNumSensorsPerWedge + sensor, so the 108 sensors
 * occupy 1000..1107. Verified to round trip over all of them.
 */
int StFstDb::sensorSurveyId(int wedge, int sensor)
{
   if (wedge < 1 || wedge > kFstNumWedges) return -1;
   if (sensor < 0 || sensor >= kFstNumSensorsPerWedge) return -1;
   return 1000 + (wedge - 1) * kFstNumSensorsPerWedge + sensor;
}


bool StFstDb::decodeSensorSurveyId(int id, int &wedge, int &sensor)
{
   wedge  = (id - 1000) / kFstNumSensorsPerWedge + 1;
   sensor = (id - 1000) % kFstNumSensorsPerWedge;
   if (id < 1000 || wedge < 1 || wedge > kFstNumWedges ||
       sensor < 0 || sensor >= kFstNumSensorsPerWedge) {
      LOG_WARN << "StFstDb: fstSensorOnWedge Id " << id << " is out of range, "
               << "expected 1000.." << 1000 + kFstNumSensors - 1 << endm;
      return false;
   }
   return true;
}


bool StFstDb::sensorLocalToGlobal(int wedge, int sensor, const double l[3], double g[3])
{
   const TGeoHMatrix *m = getHMatrixSensorOnGlobal(wedge, sensor);
   if (!m) return false;
   m->LocalToMaster(l, g);
   return true;
}


bool StFstDb::globalToSensorLocal(int wedge, int sensor, const double g[3], double l[3])
{
   const TGeoHMatrix *m = getHMatrixSensorOnGlobal(wedge, sensor);
   if (!m) return false;
   m->MasterToLocal(g, l);
   return true;
}


void StFstDb::Print(Option_t *opt) const
{
   mGeoHMatrixTpcOnGlobal->Print();
   mGeoHMatrixFstOnTpc.Print();
   mGeoHMatrixHssOnFst.Print();
   mGeoHMatrixWedgeOnHss.Print();

   for (Int_t iL = 0; iL < kFstNumWedges; iL++) {
      for (Int_t iS = 0; iS < kFstNumSensorsPerWedge; iS++) {
         mGeoHMatrixSensorOnWedge[iL][iS].Print();
      }
   }

   for (Int_t iS = 0; iS < kFstNumSensors; iS++) {
      TGeoHMatrix *sensorOnGlobal = (TGeoHMatrix *) mgRotList->FindObject(Form("R%04i", iS + 1000));
      sensorOnGlobal->Print();
   }
}
