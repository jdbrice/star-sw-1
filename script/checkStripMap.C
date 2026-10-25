// checkStripMap.C -- hold the two strip-to-position models to each other.
//
//   root4star -l -b -q 'script/checkStripMap.C("fGeom.root")'
//
// There are two legitimate implementations of "where is this strip?" and they
// must not drift apart:
//
//   CONSTANTS  StFstDb::stripPosition(), i.e. the algebra StFstHitMaker has
//              always used: kFstphiStart/kFstphiStop for the wedge edge,
//              kFstzFilp x kFstzDirct for the per-wedge and per-disk flip,
//              kFstStripPitchPhi for the step and kFstStripGapPhi for the
//              1 degree outer gap. Works with no geometry loaded, which is why
//              the producer and the slow simulation use it.
//   GEOMETRY   what TrackFitter.h:507-538 does for the fit: the sensor's own
//              AGML shape supplies phi1/phi2 and therefore the pitch and the
//              gap, and the node matrix carries the local point to global. It
//              follows the wedge mirror automatically.
//
// This macro walks the FTUS nodes exactly as FwdGeomUtils::walkFstSensors does
// (disk from z, wedge from the local +x azimuth, sensor from the shape), then
// compares the two models for every sensor and every strip:
//   sensor 0 inner : rStrip 0-3, phiStrip 0-127
//   sensor 1 outer : rStrip 4-7, phiStrip 0-63
//   sensor 2 outer : rStrip 4-7, phiStrip 64-127    <- measured on 300k
//                                                      production hits
// STATUS 2026-10-25, READ BEFORE TRUSTING A NUMBER FROM THIS MACRO.
// The CONSTANTS side is the real StFstDb API. The GEOMETRY side is a
// REIMPLEMENTATION here of what TrackFitter.h does, not the code itself, and it
// does not yet agree with it: this macro reports a mean |dr| of 2.83 cm (inner)
// and 5.8 cm (outer), i.e. a mirror of the strip direction about the wedge
// centreline, uniform across all 36 wedges, with the two outer halves swapping
// on disk 2.
//
// That disagreement is almost certainly in THIS macro, not in the chain, and the
// argument is a measurement: in the live data chain the residual compares the
// stored MuDst position (the constants model) against a track fitted with
// TrackFitter's decode, and that residual has a millimetre-scale core (dx rms
// 5.23 mm). A genuine mirror between the two decodes would put them up to 15 deg
// apart, i.e. ~4 cm at r = 15 cm, which the residual excludes outright.
//
// Nor can the usual tests adjudicate: the MC closure is blind to the convention
// because StFstSlowSimMaker turns a GEANT hit into a strip with the same algebra
// StFstHitMaker turns it back (so the round trip cancels), and
// script/fstLocalMap.C -- the "measured from data" reference quoted in
// TrackFitter.h -- derived the geometry-side convention FROM production's own
// stored hit positions, so it cannot referee the two. It also reads
// /HALL_1/CAVE_1/FSTM_1/FSTD_%d/FSTW_%d/FTUS_%d, a hierarchy that does not exist
// in any M-variant cache (no FSTD node), so it cannot even run on fGeom.root.
//
// THE WAY TO SETTLE IT is the chain's own comparison, not this macro:
// TrackFitter.h:548-556 already prints FwdHit vs the GenFit/AGML global position
// per hit. Its deltaPhi was self-cancelling (it subtracted one expression from
// itself) and is fixed; enable it with kVerbose >= 1 for a handful of events and
// read the real pairs. Then calibrate the geometry side of this macro against
// them, and only then is this a guard worth running in anger.
//
// CINT: fixed arrays, unique loop-variable names, no local std::map.
void checkStripMap(const char *cache = "fGeom.root", int fstmConfig = 2, int verbose = 0)
{
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

   TGeoManager *gm = TGeoManager::Import(cache);
   if (!gm) { printf("checkStripMap: cannot import %s\n", cache); return; }

   // ---- the geometry model: one matrix and one shape phi range per sensor ----
   TGeoHMatrix mat[108];
   double sphi1[108], sphi2[108];
   int    ok[108];
   for (int iz = 0; iz < 108; iz++) { ok[iz] = 0; sphi1[iz] = 0; sphi2[iz] = 0; }

   TGeoIterator it(gm->GetTopVolume());
   TGeoNode *nd = 0;
   int nFound = 0;
   while ((nd = it.Next())) {
      TString vn = nd->GetVolume()->GetName();
      if (!vn.BeginsWith("FTUS")) continue;
      TGeoTubeSeg *seg = (TGeoTubeSeg *) nd->GetVolume()->GetShape();
      if (!seg) continue;

      const TGeoHMatrix *here = it.GetCurrentMatrix();
      const Double_t *tr = here->GetTranslation();
      const Double_t *rr = here->GetRotationMatrix();

      int idisk = (tr[2] < 158.0) ? 0 : ((tr[2] < 172.0) ? 1 : 2);
      double phiW = TMath::ATan2(rr[3], rr[0]) * TMath::RadToDeg();
      if (phiW < 0) phiW += 360.0;
      int iwedge = TMath::Nint((75.0 - phiW) / 30.0);
      iwedge = ((iwedge % 12) + 12) % 12;
      int isens = (seg->GetRmin() < 10.0) ? 0 : ((seg->GetPhi1() > 180.0) ? 1 : 2);
      int idx = idisk * 36 + iwedge * 3 + isens;
      if (idx < 0 || idx >= 108) continue;

      mat[idx] = *here;
      sphi1[idx] = seg->GetPhi1();
      sphi2[idx] = seg->GetPhi2();
      ok[idx] = 1;
      nFound++;
   }
   printf("checkStripMap: %s, fstmConfig %d, %d FTUS nodes mapped\n",
          cache, fstmConfig, nFound);

   // ---- compare, strip by strip ---------------------------------------------
   double wedgeMean[3][12][3], wedgeSum[3][12][3];
   int    wedgeN[3][12][3];
   for (int az = 0; az < 3; az++) for (int bz = 0; bz < 12; bz++) for (int cz = 0; cz < 3; cz++) {
      wedgeMean[az][bz][cz] = 0; wedgeSum[az][bz][cz] = 0; wedgeN[az][bz][cz] = 0;
   }
   double worst[3];     // per sensor class
   double sum[3];
   int    n[3], nBad[3];
   for (int jz = 0; jz < 3; jz++) { worst[jz] = 0; sum[jz] = 0; n[jz] = 0; nBad[jz] = 0; }
   int    worstIdx = -1, worstR = -1;
   double worstPhiStrip = -1;

   for (int kd = 0; kd < 3; kd++) {
      for (int kw = 0; kw < 12; kw++) {
         for (int ks = 0; ks < 3; ks++) {
            int gi = kd * 36 + kw * 3 + ks;
            if (!ok[gi]) { printf("   sensor index %d missing from the geometry\n", gi); continue; }

            const int rLo = (0 == ks) ? 0 : 4;
            const int rHi = (0 == ks) ? 3 : 7;
            const int pLo = (2 == ks) ? 64 : 0;
            const int pHi = (0 == ks) ? 127 : ((1 == ks) ? 63 : 127);
            const int nStripSens = (0 == ks) ? 128 : 64;
            const double pitchDeg = (sphi2[gi] - sphi1[gi]) / nStripSens;

            for (int kr = rLo; kr <= rHi; kr++) {
               for (int kp = pLo; kp <= pHi; kp++) {
                  // CONSTANTS model, through the one API
                  double xyzDb[3];
                  int recoWedge = (kw + 1) + 12 * kd;
                  if (!StFstDb::stripPosition(kd + 1, recoWedge, ks, kr, (double) kp, xyzDb)) continue;

                  // GEOMETRY model, as the fit does it
                  double phiLocDeg;
                  if (0 == ks) phiLocDeg = sphi2[gi] - (kp + 0.5) * pitchDeg;
                  else {
                     double sLocal = (2 == ks) ? (kp - 64) : kp;
                     phiLocDeg = sphi1[gi] + (sLocal + 0.5) * pitchDeg;
                  }
                  double rr2 = StFstDb::stripR(kr);
                  double lp[3], gp[3];
                  lp[0] = rr2 * TMath::Cos(phiLocDeg * TMath::DegToRad());
                  lp[1] = rr2 * TMath::Sin(phiLocDeg * TMath::DegToRad());
                  lp[2] = 0.0;
                  mat[gi].LocalToMaster(lp, gp);

                  double dx = xyzDb[0] - gp[0];
                  double dy = xyzDb[1] - gp[1];
                  double d  = TMath::Sqrt(dx * dx + dy * dy);
                  sum[ks] += d; n[ks]++;
                  wedgeSum[kd][kw][ks] += d; wedgeN[kd][kw][ks]++;
                  if (d > 0.01) nBad[ks]++;
                  if (d > worst[ks]) worst[ks] = d;
                  if (d > ((worstIdx >= 0) ? worst[0] : -1) && d > worst[1] && d > worst[2]) {
                     worstIdx = gi; worstR = kr; worstPhiStrip = kp;
                  }
                  if (verbose && d > 0.01 && nBad[ks] < 4)
                     printf("   disk %d wedge %2d sensor %d rStrip %d phiStrip %3d :"
                            " db (%+8.4f,%+8.4f) geom (%+8.4f,%+8.4f)  d = %7.4f cm\n",
                            kd + 1, kw + 1, ks, kr, kp, xyzDb[0], xyzDb[1], gp[0], gp[1], d);
               }
            }
         }
      }
   }

   for (int qd = 0; qd < 3; qd++) for (int qw = 0; qw < 12; qw++) for (int qs = 0; qs < 3; qs++)
      if (wedgeN[qd][qw][qs]) wedgeMean[qd][qw][qs] = wedgeSum[qd][qw][qs] / wedgeN[qd][qw][qs];

   printf("\n  mean |dr| [cm] per disk and wedge, by sensor class\n");
   for (int ps = 0; ps < 3; ps++) {
      printf("   sensor %d :\n", ps);
      for (int pd = 0; pd < 3; pd++) {
         printf("     disk %d  ", pd + 1);
         for (int pw = 0; pw < 12; pw++) printf("%7.3f", wedgeMean[pd][pw][ps]);
         printf("\n");
      }
   }

   printf("\n  sensor   strips      mean |dr|      worst |dr|    over 0.1 mm\n");
   for (int mz = 0; mz < 3; mz++)
      printf("  %5d %8d %13.5f cm %11.5f cm %10d\n",
             mz, n[mz], (n[mz] ? sum[mz] / n[mz] : 0.0), worst[mz], nBad[mz]);
   printf("\n  A clean result is mean and worst both ~0 and zero strips over 0.1 mm.\n"
          "  The two models then agree and either can be changed with the other as\n"
          "  the reference. Any disagreement is a real convention difference: the\n"
          "  1 degree gap, the pitch, the wedge phase or the mounting parity.\n");
}
