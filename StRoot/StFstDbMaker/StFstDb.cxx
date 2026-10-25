/* $Id: StFstDb.cxx,v */

#include <assert.h>
#include "TMath.h"
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
Int_t      StFstDb::mgFstmConfig = StFstDb::kFstmAsInstalled;

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
void StFstDb::setFstmConfig(Int_t config)
{
   if (kFstmPreMirror != config && kFstmAsInstalled != config) {
      LOG_WARN << "StFstDb::setFstmConfig(" << config << "): only "
               << (int) kFstmPreMirror << " (pre-mirror) and "
               << (int) kFstmAsInstalled << " (as installed) exist; ignored" << endm;
      return;
   }

   mgFstmConfig = config;
   LOG_INFO << "StFstDb: FST wedge mounting set to "
            << (kFstmAsInstalled == config ? "AS INSTALLED (mirrored, AGML wholewedge 2.0)"
                                           : "PRE-MIRROR (AGML wholewedge 1.0)") << endm;
}
//_____________________________________________________________________________
Int_t StFstDb::fstmConfig()
{
   return mgFstmConfig;
}
//_____________________________________________________________________________
int StFstDb::wedgeOrientation(int disk, int wedge)
{
   if (disk < 1 || disk > kFstNumDisk) return 0;
   if (wedge < 1 || wedge > kFstNumWedges) return 0;

   const int iwedgem = (wedge - 1) % kFstNumWedgePerDisk;   // 0-based within the disk

   // The two arrays this needs already exist in StFstConsts.h and are the ones
   // StFstHitMaker (and StFstSlowSimMaker, line ~362) use for the strip-to-phi
   // direction, so take the front/back sense from them rather than restating it:
   //
   // CAREFUL: fstmConfig() flips the Z sense only, NOT the phi sense. That is
   // deliberate. Physically, turning a wedge over reverses both -- but our phi
   // convention is already the one the INSTALLED detector has, because it was
   // judged on data: FstWedgeAligner's 36 measured per-wedge corrections are
   // 1-3 mrad, where a reversed strip direction would be hundreds (half a wedge
   // is ~300 mrad). MC cannot judge it at all, since the slow sim and the hit
   // maker apply the same sign and the error cancels. So before 2026-10-25 the
   // phi sense was as-installed while this z sense was pre-mirror: the two
   // disagreed with each other, and kFstmAsInstalled is what makes them agree.
   //     kFstzDirct[iwedgem] = {1,-1,1,-1,...}  the alternation with wedge parity
   //     kFstzFilp[disk-1]   = {1,-1,1}         the inversion on disk 2
   // FstmGeoSM.xml is the authority behind them: a wedge goes to the BACK (-z,
   // plus <Rotation alphax="180">) when modulus2 == 0 on disks 1 and 3 but when
   // modulus2 == 1 on disk 2.
   //     <If expr="idisk.ne.1"> <If expr="...modulus2....eq.0"> ... </If> </If>
   //     <If expr="idisk.eq.1"> <If expr="...modulus2....eq.1"> ... </If> </If>
   // The sign in front is the whole question. kFstzFilp x kFstzDirct reproduces
   // AGML's wholewedge version 1.0 (modulus2 {0,1,0,1,...}) with a leading
   // MINUS: disk 1 wedge 1 has both array entries +1 and sits at -1.7408 cm in
   // every geometry built before 2026-10-25, measured on all 36 wedges.
   //
   // But that parity is not the detector. The wedges are mounted mirrored --
   // seen on the detector itself and on the installation photographs -- which is
   // why Flemming added AGML wholewedge version 2.0, modulus2 {1,0,1,0,...}, and
   // made it the module default. So kFstmAsInstalled drops the minus, and that
   // is this function's default; kFstmPreMirror restores it for reconstructing
   // samples built on the old geometry. Measured, by azimuth, on the built
   // caches: y2022c/ft2022c (wholewedge 2.0) have the opposite sign to
   // y2022b/dev2022sm (wholewedge 1.0) at all 36 wedges, and nothing else about
   // the sensor positions differs.
   const int mirror = (kFstmPreMirror == mgFstmConfig) ? -1 : +1;

   return mirror * kFstzFilp[disk - 1] * kFstzDirct[iwedgem];
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
double StFstDb::stripR(int rStrip)
{
   if (rStrip < 0 || rStrip >= kFstNumRStripsPerWedge) return 0.0;

   return kFstrStart[rStrip] + 0.5 * kFstStripPitchR;
}
//_____________________________________________________________________________
double StFstDb::stripPhi(int disk, int wedge, int sensor, int rStrip, double phiStrip)
{
   if (disk < 1 || disk > kFstNumDisk) return 0.0;
   if (wedge < 1 || wedge > kFstNumWedges) return 0.0;
   if (sensor < 0 || sensor >= kFstNumSensorsPerWedge) return 0.0;
   if (rStrip < 0 || rStrip >= kFstNumRStripsPerWedge) return 0.0;

   const int    m     = 1 + (wedge - 1) % kFstNumWedgePerDisk;   // module 1-12 in the disk
   const double pitch = kFstStripPitchPhi;
   const double dirct = kFstzDirct[m - 1];
   const double sgn   = kFstzFilp[disk - 1] * dirct;

   // The wedge's two phi edges. Disk 2 is mounted the other way round, so its
   // start and stop swap roles -- the same inversion kFstzFilp carries for z.
   double phiInner = 0.0, phiOuter = 0.0;

   if (2 == disk) {
      phiInner = kFstphiStop [m - 1] * TMath::Pi() / 6.0 - 0.5 * dirct * pitch;
      phiOuter = kFstphiStart[m - 1] * TMath::Pi() / 6.0 + 0.5 * dirct * pitch;
   }
   else {
      phiInner = kFstphiStart[m - 1] * TMath::Pi() / 6.0 + 0.5 * dirct * pitch;
      phiOuter = kFstphiStop [m - 1] * TMath::Pi() / 6.0 - 0.5 * dirct * pitch;
   }

   // INNER vs OUTER is decided by rStrip, not by sensor, which is how
   // StFstHitMaker has always done it and is reproduced here exactly. The two
   // disagree for a handful of hits (36 of 6812 in a 2000-event closure land at
   // +-1.390 cm, i.e. an inner rStrip on an outer sensor id), and those hits get
   // the inner phi with the outer z. Pre-existing; fixing it is a separate
   // change with its own before/after, not something to slip in here.
   if (rStrip < kFstNumRStripsPerWedge / 2) return phiInner + sgn * phiStrip * pitch;

   // Outer: the two halves sit either side of the 1 degree gap on the wedge
   // centreline, sensor 1 before it and sensor 2 after it.
   const double gap = (1 == sensor) ? -sgn * 0.5 * kFstStripGapPhi
                                    : (2 == sensor ?  sgn * 0.5 * kFstStripGapPhi : 0.0);

   return phiOuter - sgn * phiStrip * pitch + gap;
}
//_____________________________________________________________________________
bool StFstDb::stripPosition(int disk, int wedge, int sensor, int rStrip,
                            double phiStrip, double xyz[3])
{
   if (disk < 1 || disk > kFstNumDisk) return false;
   if (wedge < 1 || wedge > kFstNumWedges) return false;
   if (sensor < 0 || sensor >= kFstNumSensorsPerWedge) return false;
   if (rStrip < 0 || rStrip >= kFstNumRStripsPerWedge) return false;

   // An OUTER rStrip carried on sensor id 0 has no phi in this convention: the
   // gap term is defined for the two outer halves only. StFstHitMaker left
   // local[0] = local[1] = 0 in that case, i.e. put the hit on the beam line at
   // the sensor's z, and that is reproduced exactly rather than quietly improved
   // -- it is the same rStrip-vs-sensor mismatch described in stripPhi() and it
   // belongs in its own change.
   if (rStrip >= kFstNumRStripsPerWedge / 2 && 0 == sensor) {
      LOG_WARN << "StFstDb::stripPosition: outer rStrip " << rStrip
               << " on sensor id 0 (disk " << disk << " wedge " << wedge
               << "); no phi is defined, placing it on the beam line as"
               << " StFstHitMaker always has" << endm;
      xyz[0] = 0.0;
      xyz[1] = 0.0;
      xyz[2] = sensorZ(disk, wedge, sensor);
      return true;
   }

   const double r   = stripR(rStrip);
   const double phi = stripPhi(disk, wedge, sensor, rStrip, phiStrip);

   xyz[0] = r * TMath::Cos(phi);
   xyz[1] = r * TMath::Sin(phi);
   xyz[2] = sensorZ(disk, wedge, sensor);

   return true;
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
Int_t StFstDb::setGeoHMatrices(Survey_st **tables, const Int_t *nRows)
{
   //using namespace StFstConsts;

   SafeDelete(mgRotList);
   mgRotList = new THashList(kFstNumWedges * kFstNumSensorsPerWedge, 0);
   mgRotList->SetOwner(kFALSE);

   // TPC positioning relative to STAR. Loaded for the record only -- see the
   // composition below: it is NOT applied to the FST, because AGML does not.
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

   // hssOnFst has one row PER HALF and fstWedgeOnHss one row PER WEDGE. Both used
   // to be read as a single matrix from row 0 and then applied to all 36 wedges,
   // so every wedge got half 0's survey and wedge 1's survey. Harmless only while
   // the tables are identity (fstWedgeOnHss is, in every lane we carry) or while
   // the caller cancels the matrix again, which StFstHitMaker did: with the Run 22
   // tables it would have handed the x<0 half the x>0 half's offset, wrong by
   // 2 x (0.0775, -0.016) cm. FstmGeoSM.xml indexes them the same way we do here:
   //     <Misalign table="Geometry/fst/fstWedgeOnHss" row="iwedgem+12*idisk" />
   //     <Misalign table="Geometry/fst/hssOnFst"      row="halfid" />
   const Int_t nHss   = (nRows && nRows[1] > 0) ? nRows[1] : kFstNumHss;
   const Int_t nWedge = (nRows && nRows[2] > 0) ? nRows[2] : kFstNumWedges;

   if (nHss < kFstNumHss || nWedge < kFstNumWedges) {
      LOG_WARN << "StFstDb: hssOnFst has " << nHss << " rows (expected " << kFstNumHss
               << ") and fstWedgeOnHss " << nWedge << " (expected " << kFstNumWedges
               << "); missing rows are left as the identity" << endm;
   }

   for (Int_t half = 0; half < kFstNumHss; half++) {
      mGeoHMatrixHssOnFst[half].SetName(Form("hssOnFst%i", half));
      if (half >= nHss) continue;
      mGeoHMatrixHssOnFst[half].SetRotation(&hssOnFst[half].r00);
      mGeoHMatrixHssOnFst[half].SetTranslation(&hssOnFst[half].t0);
   }

   for (Int_t iw = 0; iw < kFstNumWedges; iw++) {
      mGeoHMatrixWedgeOnHss[iw].SetName(Form("wedgeOnHss%02i", iw + 1));
      if (iw >= nWedge) continue;
      mGeoHMatrixWedgeOnHss[iw].SetRotation(&wedgeOnHss[iw].r00);
      mGeoHMatrixWedgeOnHss[iw].SetTranslation(&wedgeOnHss[iw].t0);
   }
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
      const int half = halfId(wedge);
      if (half < 0) continue;

      // tpcOnGlobal is deliberately NOT in this product. It is the TPC's own
      // alignment, not part of the FST survey, and AGML does not apply it to the
      // FST: FstmGeoSM.xml places FSTM straight into CAVE and misaligns it with
      // Geometry/fst/* only, so the TGeo nodes -- which are what
      // FwdGeomUtils::buildFstSensorMap hands GenFit as the sensor planes -- are
      // displaced by fstOnTpc * hssOnFst[half] alone. Including tpcOnGlobal here
      // put every FST hit in a different frame from the plane it is fitted on:
      // 2.4/1.7/2.0 mm of translation plus 0.8 mm from the -0.48 mrad R02 acting
      // on the FST's 165 cm z lever arm, which no table of ours accounts for.
      // Without it the composed matrix equals the AGML node displacement exactly,
      // so hit and plane agree to the bit. Put it back when, and only when, the
      // XML starts misaligning the TPC-to-global step as well.
      TGeoHMatrix sensorGlobal = mGeoHMatrixFstOnTpc * mGeoHMatrixHssOnFst[half] * mGeoHMatrixWedgeOnHss[wedge - 1] * mGeoHMatrixSensorOnWedge[wedge - 1][sensor];

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


/**
 * Which half of the support structure a wedge sits on, i.e. which of the two
 * hssOnFst rows applies to it.
 *
 * FstmGeoSM.xml decides this with a cut on the wedge number,
 *
 *     halfid = 0
 *     <If expr="iwedgem.ge.6"> halfid = 1 </If>           (0-based iwedgem)
 *
 * and we deliberately do NOT restate that cut. A wedge-number cut is only as
 * good as the assumption that reco wedge numbering and AGML's iwedgem run the
 * same way round, which is exactly the kind of assumption that produced the
 * FTT positional-row bug and the FST module reflection. The halves are a
 * PHYSICAL split of the disk, so take it from the wedge's own azimuth:
 *
 *     AGML:  wedge_zrot      = 75.0 - 30.0*iwedgem        iwedgem = 0..11
 *     reco:  centre azimuth  = 15*(kFstphiStart + kFstphiStop)  degrees
 *                            = 75, 45, 15, -15, -45, -75, -105, -135, -165,
 *                              165, 135, 105     for module 1..12
 *
 * The two agree term by term, so module i IS iwedgem i-1 -- but the point is
 * that the six wedges with halfid 0 are precisely those with |azimuth| < 90,
 * i.e. the x>0 half. cos(centre) > 0 is therefore the same answer, derived
 * from geometry instead of from numbering, and it stays right if either
 * numbering is ever changed.
 *
 * NOTE this reproduces AGML's FIXED form, iwedgem.ge.6, which splits the disk
 * 6/6. The pre-generated geometry the release ships still has the original
 * .ge.5 and splits it 5/7 (see d4378cbde9: 45 and 63 sensors move, not 54 and
 * 54). Where that difference is live -- a misaligned geometry built from the
 * shipped sources, with reco applying the survey -- it shows up as exactly the
 * wedges at azimuth -75 taking the wrong half's row, 1.55 mm in x.
 */
int StFstDb::halfId(int wedge)
{
   if (wedge < 1 || wedge > kFstNumWedges) return -1;

   const int m = (wedge - 1) % kFstNumWedgePerDisk;   // module, 0-based
   const double centre = 0.5 * (kFstphiStart[m] + kFstphiStop[m]) * TMath::Pi() / 6.0;

   return (TMath::Cos(centre) > 0.0) ? 0 : 1;
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

   for (Int_t half = 0; half < kFstNumHss; half++) mGeoHMatrixHssOnFst[half].Print();

   for (Int_t iw = 0; iw < kFstNumWedges; iw++) mGeoHMatrixWedgeOnHss[iw].Print();

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
