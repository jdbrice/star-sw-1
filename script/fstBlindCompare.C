// fstBlindCompare.C -- which FST strip-to-phi convention does a blind track sit on?
//
//   root4star -l -b -q 'script/fstBlindCompare.C("fwd_blind_diag.root","fGeom.root")'
//
// Input is the fstBlind tree written by FwdTracker.h when TrackFitter:fstNoAdd is
// set (see fwd_afterburner_db.C's trailing fstNoAdd argument). With that flag the
// FST hits are found but never put on the track, and seeding is forced to FTT, so
// projX/projY come from a fit determined by the four sTGC stations alone -- it has
// never seen the silicon. That is what makes this comparison non-circular; every
// FST residual we had before was measured on tracks that contained the hit, via a
// seed the hit helped define (proposal-PR.txt section 25).
//
// For each recorded candidate it rebuilds BOTH decodes of the SAME strip:
//   STORED  production's StFstHitMaker position, as it came out of the MuDst
//   AGML    the sensor shape + node matrix, as TrackFitter does it in planar mode
// and reports which lands closer to the blind projection.
//
// CINT: fixed arrays, unique loop-variable names, no local std::map.
// trackTypeSel: -1 all, else only that StFwdTrack type -- 0 Global (FTT only,
// the FST is extrapolated, projection ~cm), 2 Primary / 4 BLCVtx (FTT + the TPC
// vertex, so the FST is INTERPOLATED and the projection is far tighter; still
// FST-blind under fstNoAdd because the vertex is the TPC's).
void fstBlindCompare(const char *fin = "fwd_blind_diag.root",
                     const char *cache = "fGeom.root", int trackTypeSel = -1)
{
   TFile *f = TFile::Open(fin);
   if (!f) { printf("fstBlindCompare: cannot open %s\n", fin); return; }
   TTree *t = (TTree *) f->Get("fstBlind");
   if (!t) { printf("fstBlindCompare: no fstBlind tree in %s\n", fin); return; }

   Int_t   disk, planeIdx, trackType = -1, nSeed = 0;
   Float_t locR, locPhi, hitX, hitY, hitZ, projX, projY, projZ;
   t->SetBranchAddress("disk", &disk);       t->SetBranchAddress("planeIdx", &planeIdx);
   t->SetBranchAddress("locR", &locR);       t->SetBranchAddress("locPhi", &locPhi);
   t->SetBranchAddress("hitX", &hitX);       t->SetBranchAddress("hitY", &hitY);
   t->SetBranchAddress("hitZ", &hitZ);
   t->SetBranchAddress("projX", &projX);     t->SetBranchAddress("projY", &projY);
   t->SetBranchAddress("projZ", &projZ);
   if (t->GetBranch("trackType")) t->SetBranchAddress("trackType", &trackType);
   if (t->GetBranch("nSeed"))     t->SetBranchAddress("nSeed", &nSeed);
   printf("fstBlindCompare: %lld candidate rows, trackTypeSel = %d\n",
          t->GetEntries(), trackTypeSel);

   // the geometry side, same walk as script/checkStripMap.C
   TGeoManager *gm = TGeoManager::Import(cache);
   if (!gm) { printf("fstBlindCompare: cannot import %s\n", cache); return; }
   TGeoHMatrix mat[108];
   double sp1[108], sp2[108], srn[108], srx[108];
   int    ok[108];
   for (int iz = 0; iz < 108; iz++) { ok[iz] = 0; }
   TGeoIterator it(gm->GetTopVolume());
   TGeoNode *nd = 0;
   while ((nd = it.Next())) {
      TString vn = nd->GetVolume()->GetName();
      if (!vn.BeginsWith("FTUS")) continue;
      TGeoTubeSeg *ts = (TGeoTubeSeg *) nd->GetVolume()->GetShape();
      if (!ts) continue;
      const TGeoHMatrix *here = it.GetCurrentMatrix();
      const Double_t *tr = here->GetTranslation();
      const Double_t *rr = here->GetRotationMatrix();
      int kdisk = (tr[2] < 158.0) ? 0 : ((tr[2] < 172.0) ? 1 : 2);
      double phiW = TMath::ATan2(rr[3], rr[0]) * TMath::RadToDeg();
      if (phiW < 0) phiW += 360.0;
      int kwedge = TMath::Nint((75.0 - phiW) / 30.0);
      kwedge = ((kwedge % 12) + 12) % 12;
      int ksens = (ts->GetRmin() < 10.0) ? 0 : ((ts->GetPhi1() > 180.0) ? 1 : 2);
      int kidx = kdisk * 36 + kwedge * 3 + ksens;
      if (kidx < 0 || kidx >= 108) continue;
      mat[kidx] = *here; sp1[kidx] = ts->GetPhi1(); sp2[kidx] = ts->GetPhi2();
      srn[kidx] = ts->GetRmin(); srx[kidx] = ts->GetRmax(); ok[kidx] = 1;
   }

   const double kPitchPhi = TMath::Pi() * 30.0 / 180.0 / 128.0;   // rad
   long nUse = 0, nStoredCloser = 0, nAgmlCloser = 0;
   double sumS = 0, sumA = 0, sumS2 = 0, sumA2 = 0;
   // The real discriminator at FST occupancy: a 40 deg x 8 cm window holds ~5
   // candidates of which at most ONE is the track's own hit, so the mean and rms
   // are background. The correct convention must show a PEAK AT ZERO; the wrong
   // one spreads its single true hit somewhere else and leaves only background
   // near zero. So count the core.
   const int    kNCore = 5;
   double coreCut[5] = {0.2, 0.5, 1.0, 2.0, 5.0};
   long   coreS[5], coreA[5];
   for (int cz = 0; cz < kNCore; cz++) { coreS[cz] = 0; coreA[cz] = 0; }

   for (Long64_t ie = 0; ie < t->GetEntries(); ie++) {
      t->GetEntry(ie);
      if (trackTypeSel >= 0 && trackType != trackTypeSel) continue;
      if (planeIdx < 0 || planeIdx >= 108 || !ok[planeIdx]) continue;
      int sensor = planeIdx % 3;
      // recover the strip indices the loader put into _localPosition
      int rStrip   = TMath::Nint((locR - 6.4375) / 2.875);
      double phiStrip = locPhi / kPitchPhi;
      if (rStrip < 0 || rStrip > 7) continue;

      // AGML decode of the same strip
      int    nStripSens = (0 == sensor) ? 128 : 64;
      double pitchDeg   = (sp2[planeIdx] - sp1[planeIdx]) / nStripSens;
      double phiLocDeg;
      if (0 == sensor) phiLocDeg = sp2[planeIdx] - (phiStrip + 0.5) * pitchDeg;
      else {
         double sLoc = (2 == sensor) ? (phiStrip - 64.0) : phiStrip;
         phiLocDeg = sp1[planeIdx] + (sLoc + 0.5) * pitchDeg;
      }
      double pitchR = (srx[planeIdx] - srn[planeIdx]) / 4.0;
      int    iRin   = (0 == sensor) ? rStrip : (rStrip - 4);
      double rGeo   = srn[planeIdx] + (iRin + 0.5) * pitchR;
      double lp[3], gp[3];
      lp[0] = rGeo * TMath::Cos(phiLocDeg * TMath::DegToRad());
      lp[1] = rGeo * TMath::Sin(phiLocDeg * TMath::DegToRad());
      lp[2] = 0.0;
      mat[planeIdx].LocalToMaster(lp, gp);

      // both against the blind projection, in r*dphi
      double pphi = TMath::ATan2(projY, projX);
      double prj  = TMath::Sqrt(projX * projX + projY * projY);
      double dS = TMath::ATan2(hitY, hitX) - pphi;
      double dA = TMath::ATan2(gp[1], gp[0]) - pphi;
      while (dS >  TMath::Pi()) dS -= TMath::TwoPi();
      while (dS < -TMath::Pi()) dS += TMath::TwoPi();
      while (dA >  TMath::Pi()) dA -= TMath::TwoPi();
      while (dA < -TMath::Pi()) dA += TMath::TwoPi();
      dS *= prj; dA *= prj;

      nUse++;
      sumS += dS; sumS2 += dS * dS;
      sumA += dA; sumA2 += dA * dA;
      if (TMath::Abs(dS) < TMath::Abs(dA)) nStoredCloser++; else nAgmlCloser++;
      for (int cy = 0; cy < kNCore; cy++) {
         if (TMath::Abs(dS) < coreCut[cy]) coreS[cy]++;
         if (TMath::Abs(dA) < coreCut[cy]) coreA[cy]++;
      }
   }

   if (!nUse) { printf("fstBlindCompare: no usable rows\n"); return; }
   double mS = sumS / nUse, mA = sumA / nUse;
   printf("\n  r*dphi against the BLIND projection, %ld candidates\n", nUse);
   printf("    STORED (production StFstHitMaker) : mean %+8.4f  rms %7.4f cm\n",
          mS, TMath::Sqrt(TMath::Max(0.0, sumS2 / nUse - mS * mS)));
   printf("    AGML   (sensor shape + node)      : mean %+8.4f  rms %7.4f cm\n",
          mA, TMath::Sqrt(TMath::Max(0.0, sumA2 / nUse - mA * mA)));
   printf("    closer to the blind track: STORED %ld (%.1f%%), AGML %ld (%.1f%%)\n",
          nStoredCloser, 100.0 * nStoredCloser / nUse, nAgmlCloser, 100.0 * nAgmlCloser / nUse);
   printf("\n  core counts -- this is the discriminator, not the mean:\n");
   printf("    |r*dphi| <      STORED      AGML     ratio S/A\n");
   for (int cx = 0; cx < kNCore; cx++)
      printf("    %6.2f cm %10ld %9ld %12.2f\n", coreCut[cx], coreS[cx], coreA[cx],
             (coreA[cx] ? (double)coreS[cx] / coreA[cx] : -1.0));
   printf("    (flat background would give both the same count in a given window;\n"
          "     an excess at small |r*dphi| is the track's own hit landing where\n"
          "     that convention says it is)\n");

   printf("\n  NOTE the selection: the record uses a 40 deg window, wider than the\n"
          "  ~25 deg the two decodes differ by, so neither is pre-selected. But the\n"
          "  candidates are all FST hits near the track in SOME frame, so read the\n"
          "  rms and the split together, not the mean alone.\n");
}
