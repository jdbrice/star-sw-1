/***************************************************************************
 *
 * StFstSlowSimMaker.cxx
 *
 * See header file for scope and known limitations.
 * Adapted from StRoot/StIstSimMaker/StIstSlowSimMaker.cxx.
 *
 ****************************************************************************/

#include "StFstSlowSimMaker.h"

#include "StMessMgr.h"
#include "TDataSetIter.h"
#include "TObjectSet.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"

#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_fstMapping_Table.h"
#include "tables/St_fstControl_Table.h"

#include "StFstDbMaker/StFstDb.h"
#include "StFstUtil/StFstCollection.h"
#include "StFstUtil/StFstRawHitCollection.h"
#include "StEvent/StFstRawHit.h"
#include "StEvent/StFstConsts.h"

ClassImp(StFstSlowSimMaker)

//______________________________________________________________________________
StFstSlowSimMaker::StFstSlowSimMaker(const char *name)
   : StMaker(name),
     mFstDb(0), mFstCollectionPtr(0), fAdc(0), mRndGen(0),
     mDefaultTimeBin(4), mCurrentTimeBinNum(kFstNumTimeBins),
     mHitEff(1.0), mAdcPerGev(0.), mStripSource(kStripFromConsts),
     mNHitsSeen(0), mNHitsUsed(0), mNHitsBadGeo(0)
{
   mMappingGeomVec.resize(kFstNumElecIds, -1);

   for (int i = 0; i < 12; i++) mSectorToModule[i] = -1;

   // Landau used only for the relative weight of each ADC time bin, same shape
   // as StIstSlowSimMaker. Not an FST measurement.
   fAdc = new TF1("fAdcFstSlowSim", "landau", 0, kFstNumTimeBins);
   fAdc->SetParameters(1.14288, 3.8, 1.168);

   // Silicon dE/dx 3.87 MeV/cm over a 300 um sensor, MIP = 441 ADC counts
   // (IST numbers -- FST needs its own calibration).
   const Double_t sidEdx          = 3.87;   // MeV/cm
   const Double_t sensorThickness = 0.03;   // cm
   const Double_t mip             = 441.0;  // ADC counts for a MIP
   mAdcPerGev = mip / (sidEdx * 1e-3 * sensorThickness);
}

//______________________________________________________________________________
StFstSlowSimMaker::~StFstSlowSimMaker()
{
   if (fAdc) { delete fAdc; fAdc = 0; }
}

//______________________________________________________________________________
Int_t StFstSlowSimMaker::Init()
{
   LOG_INFO << "StFstSlowSimMaker::Init()" << endm;

   // This is the whole point of the maker: StFstRawHitMaker::Make() looks for
   // exactly this dataset name (StFstRawHitMaker.cxx:195).
   mFstCollectionPtr = new StFstCollection();

   if (!mFstCollectionPtr) {
      LOG_ERROR << "StFstSlowSimMaker::Init() - failed to construct StFstCollection" << endm;
      return kStErr;
   }

   ToWhiteConst("fstRawAdcSimu", mFstCollectionPtr);

   return kStOk;
}

//______________________________________________________________________________
Int_t StFstSlowSimMaker::InitRun(Int_t runnumber)
{
   LOG_INFO << "StFstSlowSimMaker::InitRun() " << runnumber << endm;

   Int_t ierr = kStOk;

   // NOTE: StFstDb inherits StObject, not TDataSet, so ToWhiteConst wraps it --
   // GetDataSet() returns the wrapper and a direct cast silently yields garbage.
   // Go through GetObject().
   TObjectSet *fstDbDataSet = (TObjectSet *) GetDataSet("fst_db");

   if (!fstDbDataSet) {
      LOG_ERROR << "StFstSlowSimMaker::InitRun() - no 'fst_db' dataset. "
                << "Is StFstDbMaker in the chain, before this maker?" << endm;
      return kStErr;
   }

   mFstDb = (StFstDb *) fstDbDataSet->GetObject();

   if (!mFstDb) {
      LOG_ERROR << "StFstSlowSimMaker::InitRun() - 'fst_db' holds no StFstDb" << endm;
      return kStErr;
   }

   // ---- control table: number of ADC time bins -----------------------------
   const fstControl_st *fstControlTable = mFstDb->getControl();

   if (!fstControlTable) {
      LOG_WARN << "StFstSlowSimMaker::InitRun() - no FST control table, "
               << "using defaults (" << (int) kFstNumTimeBins << " time bins)" << endm;
   }
   else {
      mDefaultTimeBin    = fstControlTable[0].kFstDefaultTimeBin;
      mCurrentTimeBinNum = fstControlTable[0].kFstCurrentTimeBinNum;
   }

   if (mCurrentTimeBinNum < 1 || mCurrentTimeBinNum > kFstNumTimeBins) {
      LOG_WARN << "StFstSlowSimMaker::InitRun() - bad time bin count "
               << (int) mCurrentTimeBinNum << ", clamping to " << (int) kFstNumTimeBins << endm;
      mCurrentTimeBinNum = kFstNumTimeBins;
   }

   // ---- mapping table: invert elecId -> geoId ------------------------------
   // Inheriting the real electronics<->geometry relation from the DB is what
   // makes the sensor/phiStrip pairing (sensor 1 <- phiStrip 0-63,
   // sensor 2 <- 64-127) come out automatically instead of being assumed.
   const fstMapping_st *gM = mFstDb->getMapping();

   if (!gM) {
      LOG_ERROR << "StFstSlowSimMaker::InitRun() - FST mapping table is null" << endm;
      return kStErr;
   }

   for (Int_t i = 0; i < kFstNumElecIds; i++) mMappingGeomVec[i] = -1;

   Int_t nBadMap = 0;

   for (Int_t elecId = 0; elecId < kFstNumElecIds; elecId++) {
      Int_t geoId = (Int_t) gM[0].mapping[elecId];   // FST ids are 0-based

      if (geoId < 0 || geoId >= kFstNumElecIds) { nBadMap++; continue; }

      mMappingGeomVec[geoId] = elecId;
   }

   if (nBadMap)
      LOG_WARN << "StFstSlowSimMaker::InitRun() - " << nBadMap
               << " entries of fstMapping out of range" << endm;

   // ---- sector (phi/30deg) -> moduleIdx ------------------------------------
   // Built from kFstphiStart/kFstphiStop so we never depend on the GEANT
   // volume_id wedge numbering, which is not the reco wedge numbering.
   for (Int_t m = 0; m < 12; m++) {
      Int_t a   = (Int_t) kFstphiStart[m];
      Int_t b   = (Int_t) kFstphiStop[m];
      Int_t sec = (a < b) ? a : b;
      sec %= 12;

      if (sec < 0 || sec > 11) {
         LOG_ERROR << "StFstSlowSimMaker::InitRun() - bad sector " << sec
                   << " for module " << m + 1 << endm;
         ierr = kStErr;
         continue;
      }

      if (mSectorToModule[sec] > 0)
         LOG_WARN << "StFstSlowSimMaker::InitRun() - sector " << sec
                  << " claimed twice (modules " << mSectorToModule[sec]
                  << " and " << m + 1 << ")" << endm;

      mSectorToModule[sec] = m + 1;
   }

   for (Int_t s = 0; s < 12; s++) {
      if (mSectorToModule[s] < 0) {
         LOG_ERROR << "StFstSlowSimMaker::InitRun() - sector " << s << " has no module" << endm;
         ierr = kStErr;
      }
   }

   mRndGen = (TRandom3 *) gRandom;

   if (mStripSource == kStripFromTGeo) {
      LOG_WARN << "StFstSlowSimMaker::InitRun() - kStripFromTGeo is not implemented yet, "
               << "falling back to kStripFromConsts (see header, KNOWN LIMITATION)" << endm;
      mStripSource = kStripFromConsts;
   }

   LOG_INFO << "StFstSlowSimMaker::InitRun() - " << (int) mCurrentTimeBinNum
            << " time bins, hit efficiency " << mHitEff
            << ", " << mAdcPerGev << " ADC/GeV" << endm;

   return ierr;
}

//______________________________________________________________________________
Int_t StFstSlowSimMaker::Make()
{
   if (!mFstCollectionPtr) {
      LOG_ERROR << "StFstSlowSimMaker::Make() - no output collection" << endm;
      return kStErr;
   }

   mFstCollectionPtr->setNumTimeBins(mCurrentTimeBinNum);

   // FST GEANT hits live in the g2t table, NOT in StMcEvent -- there is no
   // StMcFstHit. (This is the main structural difference from IST, whose slow
   // sim reads mcEvent->istHitCollection().) Same dataset name the fast sim
   // uses, StFstFastSimMaker.cxx:167.
   St_g2t_fts_hit *hitTable = static_cast<St_g2t_fts_hit *>(GetDataSet("g2t_fsi_hit"));

   if (!hitTable) {
      LOG_INFO << "StFstSlowSimMaker::Make() - no g2t_fsi_hit table" << endm;
      return kStWarn;
   }

   const Int_t nHits = hitTable->GetNRows();

   if (nHits <= 0) {
      LOG_INFO << "StFstSlowSimMaker::Make() - g2t_fsi_hit table is empty" << endm;
      return kStWarn;
   }

   St_g2t_track *trkTable = static_cast<St_g2t_track *>(GetDataSet("g2t_track"));

   Int_t nUsed = 0, nBadGeo = 0;

   for (Int_t i = 0; i < nHits; i++) {

      const g2t_fts_hit_st *hit = (const g2t_fts_hit_st *) hitTable->At(i);

      if (!hit) continue;

      // volume_id encoding, cf. StFstFastSimMaker.cxx:203-206.
      //
      // CAREFUL: the FST disks are numbered 4, 5, 6 in the GEANT volume_id --
      // matching FSTD_4/5/6 in the TGeo tree (HALL/CAVE_1/FSTM_1/FSTD_n/FSTW_m/
      // FTUS_k) -- NOT 1, 2, 3. StFstFastSimMaker keeps that convention and
      // indexes its RMIN/RMAX[6] arrays with disk-1 = 3,4,5. Here we convert to
      // the reco disk numbering 1-3 that StFstHitMaker and StFstConsts use.
      //
      // The wedge and sensor digits are deliberately NOT used: they are the
      // GEANT numbering, which is not the reco wedge numbering, and getting
      // that wrong is one of the classic FST failure modes. Both are derived
      // from the hit position instead, in getStripsFromConsts().
      const Int_t volume_id = hit->volume_id;
      const Int_t diskRaw   = volume_id / 1000;

      if (diskRaw < 4 || diskRaw > 6) continue;   // not an FST disk

      const Int_t disk = diskRaw - 3;             // -> 1, 2, 3

      mNHitsSeen++;

      if (mHitEff < 1.0 && mRndGen && mRndGen->Rndm() > mHitEff) continue;

      Int_t wedge = -1, sensor = -1, rStrip = -1, phiStrip = -1;

      if (!getStripsFromConsts(hit->x[0], hit->x[1], hit->x[2], disk,
                               wedge, sensor, rStrip, phiStrip)) {
         nBadGeo++;
         continue;
      }

      // idTruth is the g2t track key, which is what StFstRawHitMaker propagates
      // through to the cluster and hit for embedding.
      Int_t idTruth = 0;

      if (trkTable && hit->track_p >= 0 && hit->track_p < trkTable->GetNRows()) {
         const g2t_track_st *trk = (const g2t_track_st *) trkTable->At(hit->track_p);
         if (trk) idTruth = trk->id;
      }

      generateRawHit(wedge, sensor, rStrip, phiStrip, hit->de, idTruth);
      nUsed++;
   }

   mNHitsUsed   += nUsed;
   mNHitsBadGeo += nBadGeo;

   LOG_INFO << "StFstSlowSimMaker::Make() - " << nUsed << "/" << nHits
            << " g2t hits -> " << mFstCollectionPtr->getNumRawHits()
            << " raw hits (" << nBadGeo << " outside active area)" << endm;

   return kStOk;
}

//______________________________________________________________________________
/**
 * Invert StFstHitMaker.cxx:130-165.
 *
 * NOTE: the outer-sensor gap term below is the exact inverse of
 * StFstHitMaker.cxx:159,163 AS THEY STAND TODAY, so that MC and data agree. If
 * those signs are ever changed, this must be flipped with them or MC and data
 * will disagree by one full kFstStripGapPhi. Having to say that is itself the
 * argument for doing the strip lookup from TGeo instead.
 */
Bool_t StFstSlowSimMaker::getStripsFromConsts(Double_t x, Double_t y, Double_t z,
                                              Int_t disk, Int_t &wedge, Int_t &sensor,
                                              Int_t &rStrip, Int_t &phiStrip) const
{
   const Double_t r = TMath::Sqrt(x * x + y * y);

   // ---- radial strip -------------------------------------------------------
   if (r < kFstrStart[0]) return kFALSE;

   rStrip = (Int_t) ((r - kFstrStart[0]) / kFstStripPitchR);

   if (rStrip < 0 || rStrip >= kFstNumRStripsPerWedge) return kFALSE;

   // ---- which wedge --------------------------------------------------------
   Double_t phi = TMath::ATan2(y, x);

   while (phi <  0.0)             phi += TMath::TwoPi();
   while (phi >= TMath::TwoPi())  phi -= TMath::TwoPi();

   Int_t sec = (Int_t) (phi / (TMath::Pi() / 6.0));

   if (sec < 0)  sec = 0;
   if (sec > 11) sec = 11;

   const Int_t moduleIdx = mSectorToModule[sec];

   if (moduleIdx < 1 || moduleIdx > 12) return kFALSE;

   wedge = moduleIdx + (disk - 1) * 12;

   // ---- phi strip ----------------------------------------------------------
   const Double_t pitch = kFstStripPitchPhi;
   const Double_t sgn   = kFstzFilp[disk - 1] * kFstzDirct[moduleIdx - 1];

   Double_t phiInner = 0., phiOuter = 0.;

   if (disk == 1 || disk == 3) {
      phiInner = kFstphiStart[moduleIdx - 1] * TMath::Pi() / 6.0
               + 0.5 * kFstzDirct[moduleIdx - 1] * pitch;
      phiOuter = kFstphiStop[moduleIdx - 1] * TMath::Pi() / 6.0
               - 0.5 * kFstzDirct[moduleIdx - 1] * pitch;
   }
   else {   // disk 2 is mounted the other way round
      phiInner = kFstphiStop[moduleIdx - 1] * TMath::Pi() / 6.0
               - 0.5 * kFstzDirct[moduleIdx - 1] * pitch;
      phiOuter = kFstphiStart[moduleIdx - 1] * TMath::Pi() / 6.0
               + 0.5 * kFstzDirct[moduleIdx - 1] * pitch;
   }

   // put phi on the same branch as phiInner/phiOuter before differencing
   Double_t phiRef = 0.5 * (phiInner + phiOuter);

   while (phi - phiRef >  TMath::Pi()) phi -= TMath::TwoPi();
   while (phi - phiRef < -TMath::Pi()) phi += TMath::TwoPi();

   if (rStrip < kFstNumRStripsPerWedge / 2) {
      // inner sensor: phi = phiInner + sgn*phiStrip*pitch
      sensor = 0;

      const Double_t f = (phi - phiInner) / (sgn * pitch);

      phiStrip = (Int_t) TMath::Floor(f + 0.5);
   }
   else {
      // outer: phi = phiOuter - sgn*phiStrip*pitch -/+ sgn*0.5*gap
      // The gap term depends on which outer sensor we are on, and which sensor
      // we are on depends on phiStrip, so resolve it in two passes. The gap is
      // ~2 strips wide, so only strips right at the boundary are ambiguous.
      const Double_t f0 = (phiOuter - phi) / (sgn * pitch);

      Int_t guess = (Int_t) TMath::Floor(f0 + 0.5);

      // sensor 1 <- phiStrip 0-63, sensor 2 <- phiStrip 64-127 (from fstMapping)
      sensor = (guess < kFstNumPhiSegPerWedge / 2) ? 1 : 2;

      const Double_t gapTerm = (sensor == 1)
                             ? -sgn * 0.5 * kFstStripGapPhi
                             : +sgn * 0.5 * kFstStripGapPhi;

      const Double_t f = (phiOuter + gapTerm - phi) / (sgn * pitch);

      phiStrip = (Int_t) TMath::Floor(f + 0.5);
   }

   if (phiStrip < 0 || phiStrip >= kFstNumPhiSegPerWedge) return kFALSE;

   return kTRUE;
}

//______________________________________________________________________________
void StFstSlowSimMaker::generateRawHit(Int_t wedge, Int_t sensor, Int_t rStrip,
                                       Int_t phiStrip, Double_t de, Int_t idTruth)
{
   // geoId layout, from StFstRawHit's accessors:
   //   geoId % 1024 -> sg,  phiStrip = sg % 128,  rStrip = sg / 128
   //   wedge        =  1 + geoId / 1024
   const Int_t geoId = (wedge - 1) * (kFstApvsPerWedge * kFstNumApvChannels)
                     + rStrip * kFstNumPhiSegPerWedge + phiStrip;

   if (geoId < 0 || geoId >= kFstNumElecIds) {
      LOG_WARN << "StFstSlowSimMaker - geoId " << geoId << " out of range" << endm;
      return;
   }

   const Int_t elecId = mMappingGeomVec[geoId];

   if (elecId < 0 || elecId >= kFstNumElecIds) {
      LOG_WARN << "StFstSlowSimMaker - geoId " << geoId
               << " has no elecId in fstMapping" << endm;
      return;
   }

   StFstRawHitCollection *rawHitCollectionPtr = mFstCollectionPtr->getRawHitCollection(wedge - 1);

   if (!rawHitCollectionPtr) {
      LOG_WARN << "StFstSlowSimMaker - no raw hit collection for wedge " << wedge << endm;
      return;
   }

   StFstRawHit *rawHit = rawHitCollectionPtr->getRawHit(elecId);

   if (!rawHit) {
      LOG_WARN << "StFstSlowSimMaker - no raw hit slot for elecId " << elecId << endm;
      return;
   }

   // ---- distribute the charge over the ADC time bins ------------------------
   const Int_t    maxTB = mCurrentTimeBinNum / 2;
   const Double_t mean  = fAdc->GetParameter(1);
   const Double_t sum   = fAdc->Integral(mean - 0.5 - maxTB,
                                         mean - 0.5 + (mCurrentTimeBinNum - maxTB));

   if (TMath::Abs(sum) < 1e-6) {
      LOG_WARN << "StFstSlowSimMaker - ADC time-bin shape integrates to zero" << endm;
      return;
   }

   const Double_t peakFrac = fAdc->Integral(mean - 0.5, mean + 0.5) / sum;

   if (TMath::Abs(peakFrac) < 1e-9) return;

   Float_t chargeMax = 0.;

   for (UChar_t t = 0; t < mCurrentTimeBinNum; t++) {

      // accumulate on top of whatever is already in this channel, so two GEANT
      // hits landing on the same strip add up instead of overwriting
      Float_t adcSum = (rawHit->getChannelId() >= 0) ? rawHit->getCharge(t) : 0.;

      const Double_t frac = fAdc->Integral(mean - 0.5 - maxTB + t,
                                           mean - 0.5 - maxTB + t + 1) / sum;

      const Float_t charge = de * mAdcPerGev * frac / peakFrac;

      adcSum += charge;

      rawHit->setCharge(adcSum, t);
      rawHit->setChargeErr(0., t);

      if (charge > chargeMax) chargeMax = charge;
   }

   rawHit->setChannelId(elecId);
   rawHit->setGeoId(geoId);
   rawHit->setMaxTimeBin(maxTB);
   rawHit->setIdTruth(idTruth);

   // CRITICAL: StFstScanRadiusClusterAlgo only creates a cluster when at least
   // one of its raw hits has the seed flag set (nToSeedhit > 0). Without this
   // the chain runs cleanly and produces zero clusters.
   rawHit->setSeedhitflag(1);
}

//______________________________________________________________________________
void StFstSlowSimMaker::Clear(Option_t *opts)
{
   if (mFstCollectionPtr) {
      for (UChar_t i = 0; i < kFstNumWedges; ++i) {
         StFstRawHitCollection *c = mFstCollectionPtr->getRawHitCollection(i);
         if (c) c->Clear(opts);
      }
   }

   return StMaker::Clear();
}
