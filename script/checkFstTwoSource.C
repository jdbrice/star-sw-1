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
// getHMatrixSensorOnGlobal() promises one, and in StIstDb/StSstDbMaker -- where
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

    // TEST 3: the size of what route A is missing, so the gap is quantified
    // rather than asserted. la is what globalToSensorLocal() returns, which is
    // the global position minus the survey delta, NOT a local coordinate.
    TH1F* hGapX = new TH1F("hGapX", "route A minus AGML local x (= the missing "
                           "placement);dx [cm];hits", 800, -40, 40);
    TH1F* hGapY = new TH1F("hGapY", "route A minus AGML local y (= the missing "
                           "placement);dy [cm];hits", 800, -40, 40);
    TH1F* hGapZ = new TH1F("hGapZ", "route A minus AGML local z (= the missing "
                           "placement);dz [cm];hits", 400, 140, 190);
    // And the one thing route A IS good for: the survey delta it removes.
    TH1F* hSurvX = new TH1F("hSurvX", "survey delta removed, x;dx [cm];hits", 400, -1, 1);
    TH1F* hSurvY = new TH1F("hSurvY", "survey delta removed, y;dy [cm];hits", 400, -1, 1);

    // Per-sensor accumulators, 108 = 36 wedges x 3 sensors in StFstDb's flat
    // index order. Fixed arrays because CINT dislikes local maps.
    int nSen[108], nOk2[108], nOk1[108];
    double sumDz2[108], sumDz1[108];
    for (int ia = 0; ia < 108; ia++) {
        nSen[ia] = 0; nOk2[ia] = 0; nOk1[ia] = 0; sumDz2[ia] = 0; sumDz1[ia] = 0;
    }

    int nEvUsed = 0, nHit = 0, nNoMatrix = 0, nBadId = 0, nZeroLocal = 0;
    int nThick = 0;

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

            // ---- TEST 3: quantify what route A is missing -----------------
            double la[3];
            if (!StFstDb::globalToSensorLocal(wedge, sensor, gp, la)) { nNoMatrix++; continue; }
            hGapX->Fill(la[0] - lb[0]);
            hGapY->Fill(la[1] - lb[1]);
            hGapZ->Fill(la[2] - lb[2]);
            hSurvX->Fill(la[0] - gp[0]);
            hSurvY->Fill(la[1] - gp[1]);

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
    printf("\n  TEST 3  BLOCKED -- StFstDb has no sensor placement matrix.\n");
    printf("    globalToSensorLocal() removes the survey delta and nothing else:\n");
    printf("      survey removed      dx %+8.5f   dy %+8.5f cm\n",
           hSurvX->GetMean(), hSurvY->GetMean());
    printf("    so route A minus AGML local is the whole missing placement:\n");
    printf("      dx %+9.4f   dy %+9.4f   dz %+9.4f cm\n",
           hGapX->GetMean(), hGapY->GetMean(), hGapZ->GetMean());
    printf("    A new StFstDb placement function has to close this before the\n");
    printf("    local x,y halves of step 4 can be measured. See the header.\n");

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
