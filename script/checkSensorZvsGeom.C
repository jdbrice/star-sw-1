// checkSensorZvsGeom.C -- does StFstDb::sensorZ() agree with the built geometry?
//
//   root4star -l -b -q 'script/checkSensorZvsGeom.C("fGeom.root", 2)'
//
// For every FST silicon node in a TGeo cache this works out which (disk, wedge,
// sensor) the reco calls it, then compares the node's z with
// StFstDb::sensorZ(disk, wedge, sensor). A whole-cache agreement means the reco
// puts hits on the same face of the disk as the geometry they are fitted to.
//
// The second argument is the mounting parity to test, with the same meaning as
// FSTG_fstmconfig in FstmGeoM/SM.xml and StFstDb::EFstmConfig:
//     2  as installed (mirrored)  -- the default, what y2022c / ft2022c build
//     1  pre-mirror               -- y2022..y2025, dev2021/2022/2022m/2022sm/2026
// So a cache built on wholewedge 1.0 should agree at 108/108 with 1 and at
// 0/108 with 2, and the other way round for a wholewedge 2.0 cache.
//
// Azimuth comes from the shape, not from the node translation: FTUS is a TUBS
// placed at x=y=0 with the angle carried by the rotation, so atan2(ty,tx) is 0.
// CINT: fixed arrays, unique loop-variable names, no local std::map.
void checkSensorZvsGeom(const char *cache = "fGeom.root", int fstmConfig = 2)
{
   // libStFstDbMaker needs StObject (StEvent) and the DB table wrappers; loading
   // it bare gives "undefined symbol: _ZTI8StObject".
   gSystem->Load("libStarClassLibrary.so");
   gSystem->Load("libSt_base.so");
   gSystem->Load("libStChain.so");
   gSystem->Load("libStUtilities.so");
   gSystem->Load("libStEvent.so");
   gSystem->Load("libStDb_Tables.so");
   gSystem->Load("libStDbLib.so");
   gSystem->Load("libStTpcDb.so");
   gSystem->Load("libStFstDbMaker.so");

   StFstDb::setFstmConfig(fstmConfig);
   printf("checkSensorZvsGeom: %s, fstmConfig = %d (%s)\n", cache, fstmConfig,
          fstmConfig == 1 ? "pre-mirror" : "as installed");

   TGeoManager *gm = TGeoManager::Import(cache);
   if (!gm) { printf("checkSensorZvsGeom: cannot import %s\n", cache); return; }

   double zdisk[3];
   zdisk[0] = StFstDb::diskNominalZ(1);
   zdisk[1] = StFstDb::diskNominalZ(2);
   zdisk[2] = StFstDb::diskNominalZ(3);

   int nAgree = 0, nWrongFace = 0, nOther = 0;
   double worst = 0.0;

   TGeoIterator it(gm->GetTopVolume());
   TGeoNode *nd = 0;
   while ((nd = it.Next())) {
      TString vn = nd->GetVolume()->GetName();
      if (!vn.BeginsWith("FTUS")) continue;

      TGeoTubeSeg *ts = (TGeoTubeSeg *) nd->GetVolume()->GetShape();
      if (!ts) continue;

      double rmid = 0.5 * (ts->GetRmin() + ts->GetRmax());
      double pmid = 0.5 * (ts->GetPhi1() + ts->GetPhi2());
      double lp[3], gp[3];
      lp[0] = rmid * TMath::Cos(pmid * TMath::DegToRad());
      lp[1] = rmid * TMath::Sin(pmid * TMath::DegToRad());
      lp[2] = 0.0;
      it.GetCurrentMatrix()->LocalToMaster(lp, gp);

      // disk from z, then the offset from the disk nominal
      int idisk = 0;
      for (int jd = 1; jd < 3; jd++)
         if (TMath::Abs(TMath::Abs(gp[2]) - zdisk[jd]) < TMath::Abs(TMath::Abs(gp[2]) - zdisk[idisk])) idisk = jd;
      double dz = TMath::Abs(gp[2]) - zdisk[idisk];

      // sensor class from |dz|: 1.7419 inner (reco sensor 0), 0.3469 outer (1 or 2,
      // which share a z, so 1 stands for both here)
      int sensor = (TMath::Abs(TMath::Abs(dz) - 1.7419) < TMath::Abs(TMath::Abs(dz) - 0.3469)) ? 0 : 1;

      // module 1-12 from the azimuth: module m is centred at 75 - 30*(m-1) deg
      double phi = TMath::ATan2(gp[1], gp[0]) * TMath::RadToDeg();
      int module = TMath::Nint((75.0 - phi) / 30.0) + 1;
      while (module < 1)  module += 12;
      while (module > 12) module -= 12;
      int wedge = module + 12 * idisk;

      double zreco = StFstDb::sensorZ(idisk + 1, wedge, sensor);
      double diff  = zreco - gp[2];
      if (TMath::Abs(diff) > worst) worst = TMath::Abs(diff);

      if (TMath::Abs(diff) < 0.01) nAgree++;
      else if (TMath::Abs(TMath::Abs(diff) - 2 * TMath::Abs(dz)) < 0.01) nWrongFace++;
      else {
         nOther++;
         if (nOther < 6)
            printf("   mismatch: node %s disk %d wedge %2d sensor %d phi %+7.2f"
                   "  geom z %+9.4f  reco z %+9.4f  diff %+8.4f\n",
                   nd->GetName(), idisk + 1, wedge, sensor, phi, gp[2], zreco, diff);
      }
   }

   printf("checkSensorZvsGeom: %d agree (<0.1 mm), %d on the WRONG FACE, %d other;"
          " worst |diff| %.4f cm\n", nAgree, nWrongFace, nOther, worst);
}
