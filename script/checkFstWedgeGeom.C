// checkFstWedgeGeom.C
//
// Test A of proposal_residual_mc.txt: geometry closure check for the FST,
// with no simulation and no events.
//
// Compares, for all 108 FST sensors:
//   (1) where the AGML/GEANT geometry actually puts the silicon -- read from
//       the cached TGeoManager in fGeom.root, via each FTUS node's global
//       transformation matrix and its own TGeoTubeSeg shape; against
//   (2) where StFstHitMaker::Make() thinks the strips are -- recomputed here
//       with the exact same formulas and the same StFstConsts.h constants
//       (kFstphiStart/kFstphiStop/kFstzDirct/kFstzFilp/kFstStripPitchPhi/
//       kFstStripGapPhi).
//
// (2) is the code path every REAL-data FST hit position goes through. (1) is
// the placement GEANT hits are generated in. A wedge-dependent disagreement
// between them is, by construction, a candidate explanation for the 12-fold
// (30 deg) step in mean r*dphi seen in residual/index.html -- and one that
// needs no MC digitizer to expose.
//
// A UNIFORM offset (same for all 12 wedges of a disk) is not interesting:
// it is a global rotation, absorbed by tracking. What matters is the spread
// ACROSS wedges, which is what the summary at the end reports.
//
// Node tree in fGeom.root (verified 2026-07-29):
//   HALL/CAVE_1/FSTM_1/FSTD_{4,5,6}/FSTW_{1..12}/FTUS_{1..3}
//   FSTD copy 4,5,6 = disk 1,2,3     (same numbering GEANT volume_id uses)
//   FTUS copy 3 = inner sensor (r 5-16.5), copies 1,2 = outer (r 16.5-28)
//
// Usage: root4star -b -q 'script/checkFstWedgeGeom.C'
//        root4star -b -q 'script/checkFstWedgeGeom.C("fGeom.root")'

#include "TMath.h"

// --- StFstConsts.h values, mirrored here so this macro runs under CINT ---
// Keep in sync with StRoot/StEvent/StFstConsts.h.
const int   kNDisk       = 3;
const int   kNWedge      = 12;
const int   kNSensor     = 3;
const int   kNPhiSeg     = 128;   // kFstNumPhiSegPerWedge
const double kPitchPhi   = TMath::Pi()*30.0/180.0/128.0;  // kFstStripPitchPhi
const double kGapPhi     = TMath::Pi()*1.0/180.0;         // kFstStripGapPhi
const int   kzFilp[3]    = {1,-1,1};                       // kFstzFilp
const int   kzDirct[12]  = {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1};// kFstzDirct
const double kphiStart[12] = {2.0, 2.0, 0.0, 12.0, 10.0, 10.0, 8.0, 8.0, 6.0, 6.0, 4.0, 4.0};
const double kphiStop [12] = {3.0, 1.0, 1.0, 11.0, 11.0,  9.0, 9.0, 7.0, 7.0, 5.0, 5.0, 3.0};

const double kR2D = 180.0/TMath::Pi();

// fold an angle into [0,360)
double fold360(double a){
    while (a <    0.0) a += 360.0;
    while (a >= 360.0) a -= 360.0;
    return a;
}
// signed difference a-b folded into (-180,180]
double dphi180(double a, double b){
    double d = a - b;
    while (d <= -180.0) d += 360.0;
    while (d >   180.0) d -= 360.0;
    return d;
}

void checkFstWedgeGeom(const char* geomFile = "fGeom.root"){

    // ------------------------------------------------------------------
    // 1. Read the true sensor placement out of the geometry
    // ------------------------------------------------------------------
    TFile* fg = TFile::Open(geomFile);
    if (!fg || fg->IsZombie()){ printf("ERROR: cannot open %s\n", geomFile); return; }
    TGeoManager* gm = (TGeoManager*) fg->Get("dyson");   // NOTE: key is "dyson"
    if (!gm){ printf("ERROR: no TGeoManager named 'dyson' in %s\n", geomFile); return; }

    // geant[disk][wedgeCopy][sensorCopy] : global phi of the two angular
    // edges (deg, folded), mid radius, and global z
    // CINT cannot handle function-static multi-dimensional arrays -- keep
    // everything flat and index by SIDX(disk,wedge,sensor) / WIDX(disk,wedge).
    double gLo[108], gHi[108], gRmid[108], gZ[108];
    int    gFilled[108];
    for (int ia=0; ia<108; ia++) gFilled[ia]=0;

    TGeoIterator it(gm->GetTopVolume());
    TGeoNode* nd;
    int nSeen = 0;
    while ((nd = it.Next())){
        if (strcmp(nd->GetVolume()->GetName(), "FTUS") != 0) continue;
        TString path; it.GetPath(path);

        // parse FSTD_<n> / FSTW_<n> / FTUS_<n> out of the path
        int dCopy = -1, wCopy = -1, sCopy = -1;
        {
            TString sD = path; sD.Remove(0, sD.Index("FSTD_")+5); dCopy = atoi(sD.Data());
            TString sW = path; sW.Remove(0, sW.Index("FSTW_")+5); wCopy = atoi(sW.Data());
            TString sS = path; sS.Remove(0, sS.Index("FTUS_")+5); sCopy = atoi(sS.Data());
        }
        int di = dCopy - 4;      // FSTD copy 4,5,6 -> disk index 0,1,2
        int wi = wCopy - 1;
        int si = sCopy - 1;
        if (di<0 || di>2 || wi<0 || wi>11 || si<0 || si>2){
            printf("WARN: unexpected path %s\n", path.Data());
            continue;
        }

        TGeoTubeSeg* ts = (TGeoTubeSeg*) nd->GetVolume()->GetShape();
        double rmid = 0.5*(ts->GetRmin() + ts->GetRmax());
        double p1   = ts->GetPhi1();
        double p2   = ts->GetPhi2();

        const TGeoMatrix* m = it.GetCurrentMatrix();
        double loc[3], mas1[3], mas2[3];
        loc[0] = rmid*cos(p1/kR2D); loc[1] = rmid*sin(p1/kR2D); loc[2] = 0;
        m->LocalToMaster(loc, mas1);
        loc[0] = rmid*cos(p2/kR2D); loc[1] = rmid*sin(p2/kR2D); loc[2] = 0;
        m->LocalToMaster(loc, mas2);

        double a1 = fold360(atan2(mas1[1], mas1[0])*kR2D);
        double a2 = fold360(atan2(mas2[1], mas2[0])*kR2D);

        int sidx = (di*12 + wi)*3 + si;
        gLo[sidx]   = a1;   // image of local phi1
        gHi[sidx]   = a2;   // image of local phi2
        gRmid[sidx] = rmid;
        gZ[sidx]    = mas1[2];
        gFilled[sidx] = 1;
        nSeen++;
    }
    printf("Read %d FTUS sensor nodes from %s\n\n", nSeen, geomFile);
    if (nSeen != 108) printf("WARNING: expected 108 sensors, found %d\n\n", nSeen);

    // ------------------------------------------------------------------
    // 2. Where the geometry puts each wedge (from its INNER sensor, which
    //    is one contiguous 30 deg piece -- the cleanest sector marker)
    // ------------------------------------------------------------------
    printf("=========================================================================\n");
    printf(" GEANT wedge sectors, from the inner sensor (FTUS_3, r 5-16.5 cm)\n");
    printf("=========================================================================\n");
    printf(" disk FSTWcopy   sector [deg]      centre     z [cm]\n");
    double gWedgeCen[36];
    for (int d1=0; d1<3; d1++){
        for (int w1=0; w1<12; w1++){
            int i1 = (d1*12 + w1)*3 + 2;
            if (!gFilled[i1]) continue;
            double e1 = gLo[i1], e2 = gHi[i1];
            double cen = fold360(e1 + 0.5*dphi180(e2, e1));
            gWedgeCen[d1*12 + w1] = cen;
            printf("   %d     %2d     %7.2f -> %7.2f  %7.2f   %8.3f\n",
                   d1+1, w1+1, e1, e2, cen, gZ[i1]);
        }
        printf("\n");
    }

    // ------------------------------------------------------------------
    // 3. Where StFstHitMaker puts the strips, for each reco wedge index.
    //    Formulas copied from StFstHitMaker.cxx:130-165.
    // ------------------------------------------------------------------
    // For reco wedge m (1..12) on disk dd (1..3):
    //   disk 1,3: phiInner = start*pi/6 + 0.5*dirct*pitch
    //             phiOuter = stop *pi/6 - 0.5*dirct*pitch
    //   disk 2:   phiInner/phiOuter swapped
    //   inner (sensor 0), strip k=0..127:
    //             phi = phiInner + filp*dirct*k*pitch
    //   outer sensor 1, strip k:  phi = phiOuter - filp*dirct*k*pitch - filp*dirct*0.5*gap
    //   outer sensor 2, strip k:  phi = phiOuter - filp*dirct*k*pitch + filp*dirct*0.5*gap
    printf("=========================================================================\n");
    printf(" StFstHitMaker predicted strip-centre spans, per reco wedge index\n");
    printf("=========================================================================\n");
    printf(" disk  m   inner span [deg]      centre    outer-s1 span      outer-s2 span\n");
    double rInnerLo[36], rInnerHi[36], rWedgeCen[36];
    double rS1Lo[36], rS1Hi[36], rS2Lo[36], rS2Hi[36];
    for (int d2=0; d2<3; d2++){
        for (int m2=0; m2<12; m2++){
            int    dir  = kzDirct[m2];
            int    filp = kzFilp[d2];
            double phiInner, phiOuter;
            if (d2 == 1){ // disk 2 -- swapped, StFstHitMaker.cxx
                phiInner = kphiStop [m2]*TMath::Pi()/6.0 - 0.5*dir*kPitchPhi;
                phiOuter = kphiStart[m2]*TMath::Pi()/6.0 + 0.5*dir*kPitchPhi;
            } else {
                phiInner = kphiStart[m2]*TMath::Pi()/6.0 + 0.5*dir*kPitchPhi;
                phiOuter = kphiStop [m2]*TMath::Pi()/6.0 - 0.5*dir*kPitchPhi;
            }
            double iA = (phiInner + filp*dir*  0*kPitchPhi)*kR2D;
            double iB = (phiInner + filp*dir*127*kPitchPhi)*kR2D;
            // outer: meanPhiStrip is a wedge-level index 0..127 (see
            // StFstRawHit::getPhiStrip(), = strip % 128), and the two outer
            // sensors take half of it each. WHICH half each takes is NOT
            // derivable from the code -- getSensor() is computed from
            // mChannelId (electronics) while getPhiStrip() comes from mGeoId
            // (geometry), and only StFstDb's mapping table links them.
            // Assumed here: sensor 1 <- 0..63, sensor 2 <- 64..127. The
            // opposite assignment makes the outer sensors match AGML EXACTLY.
            // Real data settles it: script/checkFstInnerOuter.C finds no
            // 0.39 cm outer-sensor structure (measured 0.0015-0.0044 cm), so
            // the opposite assignment is the true one and THERE IS NO OUTER-
            // SENSOR BUG. The +-1.000 deg printed below is the artefact of
            // this assumption, kept only to document the check.
            double s1A = (phiOuter - filp*dir*  0*kPitchPhi - filp*dir*0.5*kGapPhi)*kR2D;
            double s1B = (phiOuter - filp*dir* 63*kPitchPhi - filp*dir*0.5*kGapPhi)*kR2D;
            double s2A = (phiOuter - filp*dir* 64*kPitchPhi + filp*dir*0.5*kGapPhi)*kR2D;
            double s2B = (phiOuter - filp*dir*127*kPitchPhi + filp*dir*0.5*kGapPhi)*kR2D;

            int j2 = d2*12 + m2;
            rInnerLo[j2] = fold360(iA);  rInnerHi[j2] = fold360(iB);
            rS1Lo[j2] = fold360(s1A);    rS1Hi[j2] = fold360(s1B);
            rS2Lo[j2] = fold360(s2A);    rS2Hi[j2] = fold360(s2B);
            rWedgeCen[j2] = fold360(fold360(iA) + 0.5*dphi180(fold360(iB), fold360(iA)));

            printf("   %d   %2d  %7.2f ->%7.2f  %7.2f  %7.2f ->%7.2f  %7.2f ->%7.2f\n",
                   d2+1, m2+1, fold360(iA), fold360(iB), rWedgeCen[j2],
                   fold360(s1A), fold360(s1B), fold360(s2A), fold360(s2B));
        }
        printf("\n");
    }

    // ------------------------------------------------------------------
    // 4. Match reco wedge index <-> GEANT FSTW copy by sector centre, and
    //    report the disagreement
    // ------------------------------------------------------------------
    printf("=========================================================================\n");
    printf(" Matching, and inner-sensor disagreement (the wedge-step candidate)\n");
    printf("=========================================================================\n");
    printf(" Inner sensor: GEANT active edges vs StFstHitMaker strip centres.\n");
    printf(" The strip centres should sit half a pitch (%.3f deg) inside the\n", 0.5*kPitchPhi*kR2D);
    printf(" active edge, so a perfect match gives dLo=dHi=0 after that inset.\n\n");
    printf(" disk  m  FSTWcopy   dCentre   dLo      dHi    | dCentre @r=10cm [cm]\n");

    double dCen[36];
    int    matchCopy[36];
    for (int d3=0; d3<3; d3++){
        for (int m3=0; m3<12; m3++){
            // nearest GEANT wedge by sector centre
            int best = -1; double bestd = 1e9;
            int j3 = d3*12 + m3;
            for (int w3=0; w3<12; w3++){
                if (!gFilled[(d3*12 + w3)*3 + 2]) continue;
                double dd = fabs(dphi180(rWedgeCen[j3], gWedgeCen[d3*12 + w3]));
                if (dd < bestd){ bestd = dd; best = w3; }
            }
            matchCopy[j3] = best;
            if (best < 0) continue;

            // GEANT inner active edges, oriented to match the reco span
            double ge1 = gLo[(d3*12 + best)*3 + 2], ge2 = gHi[(d3*12 + best)*3 + 2];
            // reco strip centres are inset half a pitch from the active edge;
            // pair each reco end with whichever GEANT edge it is nearer to
            double rlo = rInnerLo[j3], rhi = rInnerHi[j3];
            double dA = dphi180(rlo, ge1), dB = dphi180(rlo, ge2);
            double dLo, dHi;
            if (fabs(dA) < fabs(dB)){ dLo = dphi180(rlo, ge1); dHi = dphi180(rhi, ge2); }
            else                    { dLo = dphi180(rlo, ge2); dHi = dphi180(rhi, ge1); }

            dCen[j3] = dphi180(rWedgeCen[j3], gWedgeCen[d3*12 + best]);
            printf("   %d   %2d     %2d     %+7.3f  %+7.3f  %+7.3f  |  %+7.4f\n",
                   d3+1, m3+1, best+1, dCen[j3], dLo, dHi,
                   dCen[j3]/kR2D*10.0);
        }
        printf("\n");
    }

    // ------------------------------------------------------------------
    // 5. Outer sensors -- same comparison
    // ------------------------------------------------------------------
    printf("=========================================================================\n");
    printf(" Outer sensors: GEANT active spans vs StFstHitMaker strip centres\n");
    printf("=========================================================================\n");
    printf(" disk  m  copy |  GEANT FTUS_1        FTUS_2      |  reco s1           s2\n");
    for (int d4=0; d4<3; d4++){
        for (int m4=0; m4<12; m4++){
            int j4 = d4*12 + m4;
            int w4 = matchCopy[j4];
            if (w4 < 0) continue;
            printf("   %d   %2d   %2d  | %6.2f-%6.2f  %6.2f-%6.2f | %6.2f-%6.2f  %6.2f-%6.2f\n",
                   d4+1, m4+1, w4+1,
                   gLo[(d4*12 + w4)*3 + 0], gHi[(d4*12 + w4)*3 + 0],
                   gLo[(d4*12 + w4)*3 + 1], gHi[(d4*12 + w4)*3 + 1],
                   rS1Lo[j4], rS1Hi[j4],
                   rS2Lo[j4], rS2Hi[j4]);
        }
        printf("\n");
    }

    // ------------------------------------------------------------------
    // 5b. Outer-sensor offset, quantified
    // ------------------------------------------------------------------
    // Pair each reco outer span with the GEANT sensor it overlaps, and
    // report the displacement of the reco strip centres from where the
    // silicon actually is. Reco strip centres should sit half a pitch
    // inside the active edge, so that inset is removed first.
    printf("=========================================================================\n");
    printf(" Outer-sensor displacement (reco strip centres - GEANT silicon)\n");
    printf("=========================================================================\n");
    printf(" half-pitch inset (%.3f deg) already removed; 0.000 = perfect.\n", 0.5*kPitchPhi*kR2D);
    printf(" r*dphi columns evaluated at the outer sensor mid-radius r=22.25 cm.\n\n");
    printf(" disk  m  dirct | ds1 [deg]  ds2 [deg] |  ds1 [cm]  ds2 [cm]\n");
    double dOut1[36], dOut2[36];
    for (int d6=0; d6<3; d6++){
        for (int m6=0; m6<12; m6++){
            int j6 = d6*12 + m6;
            int w6 = matchCopy[j6];
            if (w6 < 0) continue;
            double halfp = 0.5*kPitchPhi*kR2D;
            // GEANT edge nearest to each reco span start, inset by half a pitch
            double g1a = gLo[(d6*12 + w6)*3 + 0], g1b = gHi[(d6*12 + w6)*3 + 0];
            double g2a = gLo[(d6*12 + w6)*3 + 1], g2b = gHi[(d6*12 + w6)*3 + 1];
            // reco s1 starts at the wedge's outer edge; that edge is whichever
            // GEANT outer edge it is closest to
            double c1 = fabs(dphi180(rS1Lo[j6], g1a)) < fabs(dphi180(rS1Lo[j6], g1b)) ? g1a : g1b;
            double c2 = fabs(dphi180(rS2Lo[j6], g2a)) < fabs(dphi180(rS2Lo[j6], g2b)) ? g2a : g2b;
            // sign of the half-pitch inset follows the direction the strips run
            // CINT leaks loop-scope names across the whole function -- every
            // local in this macro must have a name used nowhere else.
            double pinset1 = dphi180(rS1Hi[j6], rS1Lo[j6]) > 0 ? +halfp : -halfp;
            double pinset2 = dphi180(rS2Hi[j6], rS2Lo[j6]) > 0 ? +halfp : -halfp;
            dOut1[j6] = dphi180(rS1Lo[j6], c1 + pinset1);
            dOut2[j6] = dphi180(rS2Lo[j6], c2 + pinset2);
            printf("   %d   %2d   %+d   | %+8.3f  %+8.3f  | %+8.4f %+8.4f\n",
                   d6+1, m6+1, kzDirct[m6], dOut1[j6], dOut2[j6],
                   dOut1[j6]/kR2D*22.25, dOut2[j6]/kR2D*22.25);
        }
        printf("\n");
    }
    printf(" NOTE: kFstStripGapPhi = %.3f deg -- the columns above reproduce it\n", kGapPhi*kR2D);
    printf(" exactly, which is the signature of the phiStrip-half assumption\n");
    printf(" above being backwards, NOT of a bug. Real data (checkFstInnerOuter.C)\n");
    printf(" shows no such structure. Treat this table as documentation of a\n");
    printf(" check that came out NEGATIVE.\n\n");

    // ------------------------------------------------------------------
    // 5c. z placement: GEANT vs the single per-disk z StFstHitMaker uses
    // ------------------------------------------------------------------
    // StFstHitMaker.cxx:167-169 assigns ONE z per disk:
    //   disk 1 -> 151.750, disk 2 -> 165.248, disk 3 -> 178.781
    printf("=========================================================================\n");
    printf(" z placement: GEANT sensor z vs the single per-disk z in StFstHitMaker\n");
    printf("=========================================================================\n");
    double recoZ[3]; recoZ[0]=151.750; recoZ[1]=165.248; recoZ[2]=178.781;
    printf(" disk  m  copy |  reco z   inner z   dz     | outer z   dz\n");
    for (int d7=0; d7<3; d7++){
        for (int m7=0; m7<12; m7++){
            int j7 = d7*12 + m7;
            int w7 = matchCopy[j7];
            if (w7 < 0) continue;
            double zin  = gZ[(d7*12 + w7)*3 + 2];
            double zout = gZ[(d7*12 + w7)*3 + 0];
            printf("   %d   %2d   %2d  | %8.3f %8.3f %+7.3f | %8.3f %+7.3f\n",
                   d7+1, m7+1, w7+1, recoZ[d7], zin, zin-recoZ[d7], zout, zout-recoZ[d7]);
        }
        printf("\n");
    }

    // ------------------------------------------------------------------
    // 6. Summary -- the number that actually matters
    // ------------------------------------------------------------------
    printf("=========================================================================\n");
    printf(" SUMMARY: wedge-to-wedge SPREAD of the inner-sensor offset\n");
    printf("=========================================================================\n");
    printf(" A uniform per-disk offset is a global rotation and is absorbed by\n");
    printf(" tracking. Only the spread ACROSS the 12 wedges can make a 12-fold\n");
    printf(" step. For scale, the step seen in data is 0.03-0.07 cm, which at\n");
    printf(" r=10 cm is %.3f-%.3f deg.\n\n", 0.03/10.0*kR2D, 0.07/10.0*kR2D);
    printf(" disk   mean [deg]   rms [deg]   max-min [deg]   max-min @r=10cm [cm]\n");
    for (int d5=0; d5<3; d5++){
        double s=0, s2=0, mn=1e9, mx=-1e9; int nn=0;
        for (int m5=0; m5<12; m5++){
            if (matchCopy[d5*12 + m5] < 0) continue;
            double v = dCen[d5*12 + m5];
            s += v; s2 += v*v; nn++;
            if (v<mn) mn=v;
            if (v>mx) mx=v;
        }
        if (nn<2) continue;
        double mean = s/nn;
        double rms  = sqrt(fabs(s2/nn - mean*mean));
        printf("   %d    %+8.4f   %8.4f    %8.4f        %8.5f\n",
               d5+1, mean, rms, mx-mn, (mx-mn)/kR2D*10.0);
    }

    printf("\n Same, for the OUTER sensors (r*dphi at r=22.25 cm):\n");
    printf(" disk   mean s1 [cm]  mean s2 [cm]   max-min s1 [cm]  max-min s2 [cm]\n");
    for (int d8=0; d8<3; d8++){
        double a1=0,a2=0,n1=0, l1=1e9,h1=-1e9,l2=1e9,h2=-1e9;
        for (int m8=0; m8<12; m8++){
            int j8 = d8*12 + m8;
            if (matchCopy[j8] < 0) continue;
            double v1 = dOut1[j8]/kR2D*22.25, v2 = dOut2[j8]/kR2D*22.25;
            a1 += v1; a2 += v2; n1++;
            if (v1<l1) l1=v1;  if (v1>h1) h1=v1;
            if (v2<l2) l2=v2;  if (v2>h2) h2=v2;
        }
        if (n1<2) continue;
        printf("   %d    %+9.4f    %+9.4f      %9.4f       %9.4f\n",
               d8+1, a1/n1, a2/n1, h1-l1, h2-l2);
    }
    printf("\nDone.\n");
}
