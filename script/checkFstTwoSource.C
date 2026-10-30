// checkFstTwoSource.C -- STEP 4 of proposal_geometry_20261019.txt, FST half:
// "the two-source geometry check". For every GEANT truth hit, compare what WE
// say about the sensor it landed on against what the SIMULATION says, from its
// own AGML geometry. Two independent routes to one number.
//
//   root4star -b -q 'script/checkFstTwoSource.C("ele2pi_np1.y2022c.vz0.run7.fzd", 20000)'
//
// Three tests, and they are not equally finished. Read the summary, not the title.
//
//   TEST 1  MOUNTING PARITY, the open hinge.  global z of the truth hit
//           vs StFstDb::sensorZ(disk, wedge, sensor) = diskNominalZ + sensorZOffset,
//           for BOTH fstmConfig values, all 108 sensors. A wedge on the wrong
//           face is off by 2 x 1.7419 cm (inner) or 2 x 0.3469 cm (outer), so
//           this is an all-or-nothing count: 108/108 for the right parity and
//           0/108 for the wrong one. Measured on the GEANT tracking geometry
//           itself rather than on an fGeom cache, which is what makes it
//           independent of checkSensorZvsGeom.C.
//
//   TEST 2  Is AGML's local z inside the silicon?  |z_local| must be under the
//           0.015 cm half-thickness. A sanity check on route B, and on the
//           g2t_fts.F fix below: if xlocal were garbage this is where it shows.
//
//   TEST 3  THE FULL LOCAL x,y COMPARISON -- NOT POSSIBLE YET, and why that
//           matters more than the number would have. See "THE GAP" below.
//
// THE GAP, measured here 2026-10-30.
// StFstDb has no function that returns where a sensor IS. The name
// getSurveyDeltaSensorOnGlobal() promises one, and in StIstDb/StSstDbMaker -- where
// the name comes from -- it delivers one, because those detectors' survey tables
// are full placements. The FST tables are not: they are MISALIGNMENT DELTAS.
// Composed, R%04i = fstOnTpc * hssOnFst[half] * wedgeOnHss[w] * sensorOnWedge[w][s]
// comes out, at sdt20211112.000004, as
//       identity rotation, translation (+0.0713, -0.3010, 0.0000) cm on half 0
//       identity rotation, translation (+0.2263, -0.3330, 0.0000) cm on half 1
// -- the same for all 18 wedges of a half, so it cannot be a placement: a
// placement could not be identical for wedges 1 to 6. No wedge azimuth, no disk
// z, no face offset. MasterToLocal() with it therefore does not enter the sensor
// frame, it only removes a few mm of survey, which is why route A below returns
// the hit's GLOBAL position and the residual it prints is the entire nominal
// placement (dz for disk 1 wedge 1 sensor 0 comes out at 153.4919 cm, which is
// exactly diskNominalZ(1) + 1.7419).
// Where the placement actually lives today: scattered across the analytic
// helpers (diskNominalZ, sensorZOffset, wedgeOrientation, agmlLocalPhi, stripPhi,
// stripR) and StFstConsts.h -- never composed into a matrix.
// So Step 4's route A needs a NEW StFstDb function that builds the nominal
// placement the way AGML builds it and then applies the survey delta. Until it
// exists, test 3 has nothing to compare against, and this macro prints the size
// of what is missing instead of pretending to measure it.
// CONSEQUENCE BEYOND THIS MACRO: TrackFitter::createAllFstPlanes() on the
// fwd-fst-db-hitmaker branch (PR #39) builds its 108 GenFit planes from
// matrixOrigin/matrixLocalX/matrixLocalY of this same matrix, so every FST plane
// is at (0.07, -0.30, 0.00) cm with u,v = global x,y. Its own comment says the
// code it replaced "failed there and every plane sat at the origin"; this is the
// same symptom from the other direction.
//
// WHY THIS NEEDED A SIMULATION FIX FIRST. STAR PR #36 added xlocal[3] to the g2t
// hit tables and filled it in g2t_fst.F and g2t_stg.F. But FST hits do not take
// the g2t_fst.F path: St_geant_Maker routes FTOS/FTIS/FTUS through
//     AddHits<St_g2t_fts_hit>("FSTH", {"FTOS","FTIS","FTUS"}, "g2t_fsi_hit", g2t_fts)
// and the g2t_fst AddHits line next to it is commented out -- g2t_fst.F is never
// called at all, and the FDSW volume it names exists in no current FST geometry.
// So xlocal was identically zero for every FST hit. g2t_fts.F now fills it the
// same way g2t_fst.F does; without that commit this macro measures nothing.
//
// Inputs: any .fzd whose geometry declares the local hit set for FTUS -- all of
// FstmGeo/GeoM/GeoSM do, <Hit meas="x"/"y"/"z"> at 0.00001 cm bins -- so existing
// samples work and nothing has to be re-simulated. Only the g2t translation
// changed, and that runs here, at reco time.
//
// CINT: unique loop variable names throughout, fixed-size C arrays, no local
// std::map, filename matches the function name.
void checkFstTwoSource(const char* fzd = "ele2pi_np1.y2022c.vz0.run7.fzd",
                       int nev = 20000,
                       const char* tag = "y2022c",
                       const char* sdt = "20211112.000004",
                       int misalign = 1,
                       const char* out = "checkFstTwoSource.root")
{
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");

    // No reconstruction at all: this compares truth against StFstDb, so the
    // chain needs only the geometry, the geant reader and StFstDbMaker.
    // "misalign" is what makes AgPosition::SetReal() apply the <Misalign>
    // tables to the AGML placement, so it decides whether route B is misaligned.
    TString ch = Form("fzin sdt%s %s agml usexgeom %s geant gstar fstDb",
                      sdt, tag, misalign ? "misalign" : "");
    printf("checkFstTwoSource: chain = [%s]\n", ch.Data());
    printf("checkFstTwoSource: fzd = %s, nev = %d, misalign = %d\n", fzd, nev, misalign);
    bfc(-1, ch, fzd);

    StChain* chain = (StChain*)StMaker::GetChain();
    if (chain->Init() != kStOK) { printf("checkFstTwoSource: chain Init FAILED\n"); return; }

    TFile* fout = new TFile(out, "RECREATE");

    // TEST 1: truth z minus StFstDb::sensorZ. The right parity sits in the
    // 0.03 cm silicon; the wrong one is 3.48 cm (inner) or 0.69 cm (outer) away,
    // so bin wide enough to hold both and finely enough to see the silicon.
    TH1F* hDzCfg2 = new TH1F("hDzCfg2", "truth z - StFstDb::sensorZ, fstmConfig 2 "
                             "(as installed);dz [cm];hits", 1000, -5, 5);
    TH1F* hDzCfg1 = new TH1F("hDzCfg1", "truth z - StFstDb::sensorZ, fstmConfig 1 "
                             "(pre-mirror);dz [cm];hits", 1000, -5, 5);
    TH2F* hDzCfg2VsWedge = new TH2F("hDzCfg2VsWedge", "dz (cfg 2) vs wedge;wedge;dz [cm]",
                                    36, 0.5, 36.5, 400, -5, 5);
    TH2F* hDzCfg1VsWedge = new TH2F("hDzCfg1VsWedge", "dz (cfg 1) vs wedge;wedge;dz [cm]",
                                    36, 0.5, 36.5, 400, -5, 5);

    // TEST 2: AGML's own local z, which must be inside the silicon.
    TH1F* hZloc = new TH1F("hZlocGeant", "AGML local z of the hit;z_{local} [cm];hits",
                           400, -0.05, 0.05);
    TH1F* hRloc = new TH1F("hRlocGeant", "AGML local r of the hit;r_{local} [cm];hits",
                           400, 0, 30);

    // TEST 3: the real two-source comparison, through
    // StFstDb::sensorPlacement(). Both variants every event, because which one
    // should match depends on whether AGML was misaligned:
    //   S = placement WITH the survey delta (what StFstHitMaker puts hits on)
    //   N = nominal placement only
    // 2 um bins over +-0.1 cm resolve the GEANT quantisation (0.1 um) and still
    // contain the fstOnTpc-sized disagreement we expect; the wide pair catches
    // anything structural.
    TH1F* hSx = new TH1F("hSdx", "surveyed route - AGML, local x;dx [cm];hits", 1000, -0.1, 0.1);
    TH1F* hSy = new TH1F("hSdy", "surveyed route - AGML, local y;dy [cm];hits", 1000, -0.1, 0.1);
    TH1F* hSz = new TH1F("hSdz", "surveyed route - AGML, local z;dz [cm];hits", 1000, -0.1, 0.1);
    TH1F* hNx = new TH1F("hNdx", "nominal route - AGML, local x;dx [cm];hits", 1000, -0.1, 0.1);
    TH1F* hNy = new TH1F("hNdy", "nominal route - AGML, local y;dy [cm];hits", 1000, -0.1, 0.1);
    TH1F* hNz = new TH1F("hNdz", "nominal route - AGML, local z;dz [cm];hits", 1000, -0.1, 0.1);
    TH1F* hSxW = new TH1F("hSdxWide", "surveyed route - AGML, local x, wide;dx [cm];hits",
                          800, -20, 20);
    TH1F* hSyW = new TH1F("hSdyWide", "surveyed route - AGML, local y, wide;dy [cm];hits",
                          800, -20, 20);
    TH1F* hSzW = new TH1F("hSdzWide", "surveyed route - AGML, local z, wide;dz [cm];hits",
                          800, -20, 20);
    TH1F* hS3  = new TH1F("hS3", "|surveyed route - AGML|;distance [cm];hits", 1000, 0, 5);
    TH2F* hS3VsIndex = new TH2F("hS3VsIndex", "|surveyed route - AGML| vs flat sensor "
                                "index;index 0-107;distance [cm]", 108, -0.5, 107.5, 500, 0, 5);
    TH2F* hSxVsWedge = new TH2F("hSdxVsWedge", "surveyed dx vs wedge;wedge;dx [cm]",
                                36, 0.5, 36.5, 400, -1, 1);
    TH2F* hSyVsWedge = new TH2F("hSdyVsWedge", "surveyed dy vs wedge;wedge;dy [cm]",
                                36, 0.5, 36.5, 400, -1, 1);

    // Per-sensor accumulators, 108 = 36 wedges x 3 sensors in StFstDb's flat
    // index order. Fixed arrays because CINT dislikes local maps.
    int nSen[108], nOk2[108], nOk1[108];
    double sumDz2[108], sumDz1[108];
    double sumSx[108], sumSy[108], sumSz[108], maxS3[108];
    for (int ia = 0; ia < 108; ia++) {
        nSen[ia] = 0; nOk2[ia] = 0; nOk1[ia] = 0; sumDz2[ia] = 0; sumDz1[ia] = 0;
        sumSx[ia] = 0; sumSy[ia] = 0; sumSz[ia] = 0; maxS3[ia] = 0;
    }

    int nEvUsed = 0, nHit = 0, nNoMatrix = 0, nBadId = 0, nZeroLocal = 0;
    int nThick = 0;
    int nS10 = 0, nS100 = 0, nS1mm = 0, nSworse = 0;

    // The silicon is 0.03 cm thick, so anything inside 0.05 cm of the nominal
    // face is "on it" and the two faces are 3.48 cm apart -- no ambiguity.
    const double kTolZ = 0.05;

    for (int iev = 0; iev < nev; iev++) {
        if (chain->Make() != kStOK) break;
        nEvUsed++;

        // geant/g2t_fsi_hit holds the FST: it is a St_g2t_fts_hit (rows are
        // g2t_fts_hit_st), x[3] is GLOBAL cm and xlocal[3] is FTUS-local cm.
        St_g2t_fts_hit* fsi = (St_g2t_fts_hit*)chain->GetDataSet("geant/g2t_fsi_hit");
        if (!fsi) { chain->Clear(); continue; }

        for (int ig = 0; ig < fsi->GetNRows(); ig++) {
            g2t_fts_hit_st* gh = (g2t_fts_hit_st*)fsi->At(ig);

            // volume_id = 1000*disk + 10*wedge + sensor, disk 4-6 "for historical
            // reasons", wedge 1-12 WITHIN the disk, sensor 1-3 with 1 = inner
            // (g2t_volume_id.g, g2t_fst_volume_id).
            int vid    = gh->volume_id;
            int gdisk  = vid / 1000;
            int gwedge = (vid % 1000) / 10;
            int gsens  = vid % 10;
            if (gdisk < 4 || gdisk > 6 || gwedge < 1 || gwedge > 12 ||
                gsens < 1 || gsens > 3) { nBadId++; continue; }

            // to StFstDb's convention: disk 1-3, wedge 1-36 global disk-major,
            // sensor 0-2 with 0 = inner
            int disk   = gdisk - 3;
            int wedge  = (gdisk - 4) * 12 + gwedge;
            int sensor = gsens - 1;
            // The flat 0-107 index. Spelled out rather than calling
            // StFstDb::sensorIndex(), which is PR #39 and not in this tree yet;
            // it is the same flattening, and the same as sensorSurveyId() - 1000.
            int idx    = (wedge - 1) * 3 + sensor;
            if (idx < 0 || idx > 107) { nBadId++; continue; }

            // route B: AGML's own local position of the step point
            double lb[3];
            lb[0] = gh->xlocal[0]; lb[1] = gh->xlocal[1]; lb[2] = gh->xlocal[2];

            // All three exactly zero means xlocal was never written -- the
            // pre-fix state of g2t_fts.F. Count them rather than compare to 0.
            if (0.0 == lb[0] && 0.0 == lb[1] && 0.0 == lb[2]) { nZeroLocal++; continue; }

            double gp[3];
            gp[0] = gh->x[0]; gp[1] = gh->x[1]; gp[2] = gh->x[2];

            // ---- TEST 1: the mounting parity, both ways -------------------
            StFstDb::setFstmConfig(2);
            double zn2 = StFstDb::sensorZ(disk, wedge, sensor);
            StFstDb::setFstmConfig(1);
            double zn1 = StFstDb::sensorZ(disk, wedge, sensor);
            StFstDb::setFstmConfig(2);   // leave it on the default

            double dz2 = gp[2] - zn2;
            double dz1 = gp[2] - zn1;
            hDzCfg2->Fill(dz2); hDzCfg1->Fill(dz1);
            hDzCfg2VsWedge->Fill(wedge, dz2);
            hDzCfg1VsWedge->Fill(wedge, dz1);
            sumDz2[idx] += dz2; sumDz1[idx] += dz1;
            if (TMath::Abs(dz2) < kTolZ) nOk2[idx]++;
            if (TMath::Abs(dz1) < kTolZ) nOk1[idx]++;

            // ---- TEST 2: AGML's local z inside the silicon ----------------
            hZloc->Fill(lb[2]);
            hRloc->Fill(TMath::Sqrt(lb[0]*lb[0] + lb[1]*lb[1]));
            if (TMath::Abs(lb[2]) > 0.015) nThick++;

            // ---- TEST 3: the two-source comparison ------------------------
            // StFstDb::sensorPlacement() is local->global, so MasterToLocal()
            // takes the truth hit's global position into the sensor frame, which
            // is the frame AGML reported xlocal in.
            TGeoHMatrix ms = StFstDb::sensorPlacement(disk, wedge, sensor, true);
            TGeoHMatrix mn = StFstDb::sensorPlacement(disk, wedge, sensor, false);
            double ls[3], ln[3];
            ms.MasterToLocal(gp, ls);
            mn.MasterToLocal(gp, ln);

            double sdx = ls[0] - lb[0], sdy = ls[1] - lb[1], sdz = ls[2] - lb[2];
            double ndx = ln[0] - lb[0], ndy = ln[1] - lb[1], ndz = ln[2] - lb[2];
            double sd3 = TMath::Sqrt(sdx*sdx + sdy*sdy + sdz*sdz);

            hSx->Fill(sdx); hSy->Fill(sdy); hSz->Fill(sdz);
            hNx->Fill(ndx); hNy->Fill(ndy); hNz->Fill(ndz);
            hSxW->Fill(sdx); hSyW->Fill(sdy); hSzW->Fill(sdz);
            hS3->Fill(sd3);
            hS3VsIndex->Fill(idx, sd3);
            hSxVsWedge->Fill(wedge, sdx);
            hSyVsWedge->Fill(wedge, sdy);

            sumSx[idx] += sdx; sumSy[idx] += sdy; sumSz[idx] += sdz;
            if (sd3 > maxS3[idx]) maxS3[idx] = sd3;

            if      (sd3 < 0.001) nS10++;
            else if (sd3 < 0.010) nS100++;
            else if (sd3 < 0.100) nS1mm++;
            else                  nSworse++;

            nSen[idx]++;
            nHit++;
        }
        chain->Clear();
    }

    // ---- summary ------------------------------------------------------------
    printf("\n");
    printf("checkFstTwoSource: %d events, %d compared FST hits\n", nEvUsed, nHit);
    printf("  skipped: %d unusable volume_id, %d no survey matrix, %d xlocal all-zero\n",
           nBadId, nNoMatrix, nZeroLocal);
    if (nZeroLocal > 0 && nHit == 0) {
        printf("\n  xlocal is identically zero for every FST hit. Either this .fzd's\n");
        printf("  geometry did not declare the local hit set for FTUS, or g2t_fts.F in\n");
        printf("  this build does not fill xlocal yet (see the header comment).\n");
    }
    if (nHit == 0) { fout->Write(); fout->Close(); return; }

    // TEST 1
    int nSeen = 0, nAll2 = 0, nAll1 = 0, nSplit = 0;
    for (int ib = 0; ib < 108; ib++) {
        if (nSen[ib] == 0) continue;
        nSeen++;
        if (nOk2[ib] == nSen[ib])      nAll2++;
        else if (nOk1[ib] == nSen[ib]) nAll1++;
        else                           nSplit++;
    }
    printf("\n  TEST 1  mounting parity, from the GEANT tracking geometry\n");
    printf("    sensors seen (of 108)                         %d\n", nSeen);
    printf("    all hits agree with fstmConfig 2 (as installed) %d\n", nAll2);
    printf("    all hits agree with fstmConfig 1 (pre-mirror)   %d\n", nAll1);
    printf("    neither, or mixed                               %d\n", nSplit);
    printf("    dz rms: cfg2 %8.5f cm   cfg1 %8.5f cm\n",
           hDzCfg2->GetRMS(), hDzCfg1->GetRMS());

    // TEST 2
    printf("\n  TEST 2  AGML local z of the hit (silicon half-thickness 0.015 cm)\n");
    printf("    mean %9.6f  rms %9.6f cm, %d of %d outside the silicon\n",
           hZloc->GetMean(), hZloc->GetRMS(), nThick, nHit);
    printf("    local r  mean %7.3f  rms %7.3f cm\n", hRloc->GetMean(), hRloc->GetRMS());

    // TEST 3
    printf("\n  TEST 3  two-source comparison, StFstDb::sensorPlacement vs AGML\n");
    printf("    AGML was built %s\n", misalign ? "MISALIGNED (survey applied)"
                                               : "IDEAL (no survey)");
    printf("    local residual, mean +- rms [um]        dx              dy              dz\n");
    printf("      surveyed placement   %+10.2f %6.2f %+10.2f %6.2f %+10.2f %6.2f\n",
           hSx->GetMean()*1e4, hSx->GetRMS()*1e4,
           hSy->GetMean()*1e4, hSy->GetRMS()*1e4,
           hSz->GetMean()*1e4, hSz->GetRMS()*1e4);
    printf("      nominal placement    %+10.2f %6.2f %+10.2f %6.2f %+10.2f %6.2f\n",
           hNx->GetMean()*1e4, hNx->GetRMS()*1e4,
           hNy->GetMean()*1e4, hNy->GetRMS()*1e4,
           hNz->GetMean()*1e4, hNz->GetRMS()*1e4);
    printf("    (the fine histograms are +-0.1 cm; a mean of 0.00 with rms 0.00 means\n");
    printf("     everything fell in the overflow -- read the wide ones in the file)\n");

    printf("\n    |surveyed route - AGML|      hits     fraction\n");
    printf("      < 10 um              %8d     %6.2f %%\n", nS10,    100.0*nS10   /nHit);
    printf("      10 - 100 um          %8d     %6.2f %%\n", nS100,   100.0*nS100  /nHit);
    printf("      0.1 - 1 mm           %8d     %6.2f %%\n", nS1mm,   100.0*nS1mm  /nHit);
    printf("      > 1 mm               %8d     %6.2f %%\n", nSworse, 100.0*nSworse/nHit);

    int nS3ok = 0, nS3seen = 0;
    double worstS3 = 0.0;
    for (int id = 0; id < 108; id++) {
        if (nSen[id] == 0) continue;
        nS3seen++;
        if (maxS3[id] < 0.010) nS3ok++;
        if (maxS3[id] > worstS3) worstS3 = maxS3[id];
    }
    printf("    per sensor: %d of %d agree to better than 100 um, worst %8.4f cm\n",
           nS3ok, nS3seen, worstS3);

    // Split by mounting orientation. A wrong ROTATION ORDER in the placement
    // (Rz*Rx versus Rx*Rz) is invisible on unflipped wedges and wrong on flipped
    // ones, so this separates it from anything that affects all 108 equally.
    int    nFront = 0, nBack = 0;
    double sFrontX = 0, sFrontY = 0, sBackX = 0, sBackY = 0;
    double mFront = 0, mBack = 0;
    for (int ie = 0; ie < 108; ie++) {
        if (nSen[ie] == 0) continue;
        int ew = ie / 3 + 1;
        int ed = 1 + (ew - 1) / 12;
        if (StFstDb::wedgeOrientation(ed, ew) > 0) {
            nFront++; sFrontX += sumSx[ie]/nSen[ie]; sFrontY += sumSy[ie]/nSen[ie];
            if (maxS3[ie] > mFront) mFront = maxS3[ie];
        } else {
            nBack++;  sBackX  += sumSx[ie]/nSen[ie]; sBackY  += sumSy[ie]/nSen[ie];
            if (maxS3[ie] > mBack) mBack = maxS3[ie];
        }
    }
    printf("\n    split by mounting face (orientation from wedgeOrientation)\n");
    printf("      face   sensors   <dx> [um]   <dy> [um]   worst |d| [cm]\n");
    if (nFront) printf("      front  %7d %11.1f %11.1f %14.4f\n",
                       nFront, 1e4*sFrontX/nFront, 1e4*sFrontY/nFront, mFront);
    if (nBack)  printf("      back   %7d %11.1f %11.1f %14.4f\n",
                       nBack, 1e4*sBackX/nBack, 1e4*sBackY/nBack, mBack);

    printf("\n    per-sensor surveyed residual, first 12 wedges of disk 1\n");
    printf("      wedge sens  azim  orient    <dx> [um]   <dy> [um]   max|d| [cm]\n");
    for (int ig = 0; ig < 36; ig++) {
        if (nSen[ig] == 0) continue;
        int gw = ig / 3 + 1;
        int gs = ig % 3;
        printf("      %5d %4d %5.0f %7d %12.1f %11.1f %13.4f\n",
               gw, gs, StFstDb::wedgeAzimuth(gw), StFstDb::wedgeOrientation(1, gw),
               1e4*sumSx[ig]/nSen[ig], 1e4*sumSy[ig]/nSen[ig], maxS3[ig]);
    }

    // Per-sensor detail only where the parity verdict is not unanimous.
    if (nSplit > 0 || (nAll2 > 0 && nAll1 > 0)) {
        printf("\n  sensors whose parity is not clean:\n");
        printf("    idx  disk wedge sens    hits   ok(cfg2) ok(cfg1)   <dz2>     <dz1>\n");
        for (int ic = 0; ic < 108; ic++) {
            if (nSen[ic] == 0) continue;
            if (nOk2[ic] == nSen[ic] || nOk1[ic] == nSen[ic]) continue;
            int cw = ic / 3 + 1;   // inverse of (wedge-1)*3 + sensor
            int cs = ic % 3;
            printf("    %3d  %4d %5d %4d %7d %9d %8d %9.4f %9.4f\n",
                   ic, 1 + (cw - 1) / 12, cw, cs, nSen[ic], nOk2[ic], nOk1[ic],
                   sumDz2[ic]/nSen[ic], sumDz1[ic]/nSen[ic]);
        }
    }

    fout->Write();
    fout->Close();
    printf("\ncheckFstTwoSource: wrote %s\n", out);
}
