// checkFstStripPhi.C
//
// The concrete version of the outer-sensor question: print the wedge-LOCAL phi
// that AGML and StFstHitMaker each assign to phiStrip = 0, 63, 64, 127.
//
// FRAME. Both columns are in the FSTW (wedge container) local frame: wedge
// centreline at local phi = 0, wedge spanning -15..+15 deg. That is the frame
// AGML defines the FTUS sensor shapes in -- the inner sensor's shape is
// literally phi1=345, phi2=375.
//
// The FSTW placement is a rotation about z, with or without the alphax=180
// mirror, so local -> global is always
//     phi_global = C + sgn * phi_local,   sgn = +-1
// and C, sgn follow from the two shape edges and their images under the node's
// own matrix. That inverts analytically, which is how the StFstHitMaker global
// angles are brought into the local frame here -- no assumed convention.
//
// THE AMBIGUITY BEING RESOLVED. StFstHitMaker's outer branch is selected by
// the SENSOR index (sensor==1 vs sensor==2), but meanPhiStrip is a wedge-level
// 0..127 index (StFstRawHit::getPhiStrip() = strip % 128). Which sensor holds
// which half of that range is set by StFstDb's electronics->geometry mapping
// and is invisible in the code, so BOTH branches are evaluated at all four
// strip numbers and the comparison against AGML picks the assignment.
//
// Usage: root4star -b -q 'script/checkFstStripPhi.C'
//        root4star -b -q 'script/checkFstStripPhi.C("fGeom.root", 1, 1)'

const double kR2Dg     = 57.29577951308232;   // 180/pi
const double kPitchDeg = 30.0/128.0;          // kFstStripPitchPhi, deg
const double kGapDeg   = 1.0;                 // kFstStripGapPhi, deg

double fold180b(double a){
    while (a <= -180.0) a += 360.0;
    while (a >   180.0) a -= 360.0;
    return a;
}

void checkFstStripPhi(const char* geomFile = "fGeom.root",
                      int detailDisk = 1, int detailCopy = 1){

    const int kzD[12] = {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1};   // kFstzDirct
    const int kzF[3]  = {1,-1,1};                           // kFstzFilp
    const double kSt[12] = {2.0,2.0,0.0,12.0,10.0,10.0,8.0,8.0,6.0,6.0,4.0,4.0}; // kFstphiStart
    const double kSp[12] = {3.0,1.0,1.0,11.0,11.0, 9.0,9.0,7.0,7.0,5.0,5.0,3.0}; // kFstphiStop

    TFile* fq = TFile::Open(geomFile);
    if (!fq){ printf("ERROR: cannot open %s\n", geomFile); return; }
    TGeoManager* gq = (TGeoManager*) fq->Get("dyson");
    if (!gq){ printf("ERROR: no TGeoManager 'dyson'\n"); return; }
    printf("geometry loaded\n"); fflush(stdout);

    // Pull out just the three FTUS sensors of the requested wedge.
    // sLo/sHi  = shape local phi edges;  gLo/gHi = their global images
    double sLo[3], sHi[3], sRm[3], gLo[3], gHi[3];
    int    got[3]; got[0]=0; got[1]=0; got[2]=0;

    TString want = Form("FSTD_%d/FSTW_%d/FTUS_", detailDisk+3, detailCopy);
    TGeoIterator it(gq->GetTopVolume());
    TGeoNode* nq;
    int nfound = 0;
    while ((nq = it.Next())){
        if (strcmp(nq->GetVolume()->GetName(), "FTUS") != 0) continue;
        TString pq; it.GetPath(pq);
        if (!pq.Contains(want)) continue;
        TString sq = pq; sq.Remove(0, sq.Index("FTUS_")+5);
        int si = atoi(sq.Data()) - 1;
        if (si < 0 || si > 2) continue;

        TGeoTubeSeg* tq = (TGeoTubeSeg*) nq->GetVolume()->GetShape();
        sLo[si] = tq->GetPhi1();
        sHi[si] = tq->GetPhi2();
        sRm[si] = 0.5*(tq->GetRmin() + tq->GetRmax());

        const TGeoMatrix* mq = it.GetCurrentMatrix();
        double la[3], ma[3], lb[3], mb[3];
        la[0] = sRm[si]*cos(sLo[si]/kR2Dg); la[1] = sRm[si]*sin(sLo[si]/kR2Dg); la[2] = 0;
        lb[0] = sRm[si]*cos(sHi[si]/kR2Dg); lb[1] = sRm[si]*sin(sHi[si]/kR2Dg); lb[2] = 0;
        mq->LocalToMaster(la, ma);
        mq->LocalToMaster(lb, mb);
        gLo[si] = atan2(ma[1], ma[0])*kR2Dg;
        gHi[si] = atan2(mb[1], mb[0])*kR2Dg;
        got[si] = 1;
        nfound++;
        if (nfound == 3) break;      // stop early; no need to walk 3185 volumes
    }
    if (!got[0] || !got[1] || !got[2]){ printf("ERROR: wedge not found (%s)\n", want.Data()); return; }
    printf("found 3 FTUS sensors for %s\n\n", want.Data()); fflush(stdout);

    // local -> global map, taken from the inner sensor (all three FTUS share
    // the FSTW rotation; they differ only by a z translation)
    double dLoc = fold180b(sHi[2] - sLo[2]);
    double dGlo = fold180b(gHi[2] - gLo[2]);
    double sgn  = (dLoc*dGlo > 0) ? +1.0 : -1.0;
    double cOff = fold180b(gLo[2] - sgn*fold180b(sLo[2]));
    printf("local -> global map:  phi_global = %+.4f %+.1f * phi_local\n\n", cOff, sgn);

    // which electronic wedge index m is this FSTW copy?
    double gcen = fold180b(gLo[2] + 0.5*fold180b(gHi[2] - gLo[2]));
    int mS = -1; double mB = 1e9;
    for (int mm = 0; mm < 12; mm++){
        double c1 = kSt[mm]*30.0, c2 = kSp[mm]*30.0;
        double cc = fold180b(c1 + 0.5*fold180b(c2 - c1));
        double dd = fabs(fold180b(cc - gcen));
        if (dd < mB){ mB = dd; mS = mm; }
    }

    int dI = detailDisk - 1;
    int dr = kzD[mS];
    int fl = kzF[dI];
    double phIn, phOut;
    if (dI == 1){                                   // disk 2: swapped
        phIn  = kSp[mS]*30.0 - 0.5*dr*kPitchDeg;
        phOut = kSt[mS]*30.0 + 0.5*dr*kPitchDeg;
    } else {
        phIn  = kSt[mS]*30.0 + 0.5*dr*kPitchDeg;
        phOut = kSp[mS]*30.0 - 0.5*dr*kPitchDeg;
    }

    printf("=========================================================================\n");
    printf(" disk %d, FSTW copy %d  ->  electronic wedge m=%d\n", detailDisk, detailCopy, mS+1);
    printf(" kFstphiStart[m]=%.0f (=%.1f deg)  kFstphiStop[m]=%.0f (=%.1f deg)\n",
           kSt[mS], kSt[mS]*30.0, kSp[mS], kSp[mS]*30.0);
    printf(" kFstzDirct[m]=%+d   kFstzFilp[disk]=%+d\n", dr, fl);
    printf(" phiInner=%+.4f deg   phiOuter=%+.4f deg   (both global)\n", phIn, phOut);
    printf("=========================================================================\n");

    printf("\nAGML active silicon, from the FTUS shapes, in LOCAL phi [deg]:\n");
    printf("  sensor        r range      active edges          strip centres\n");
    printf("  FTUS_3 inner  %4.1f-%4.1f   %+8.3f .. %+8.3f   %+8.3f .. %+8.3f\n",
           5.0, 16.5, fold180b(sLo[2]), fold180b(sHi[2]),
           fold180b(sLo[2])+0.5*kPitchDeg, fold180b(sHi[2])-0.5*kPitchDeg);
    printf("  FTUS_1 outer  %4.1f-%4.1f   %+8.3f .. %+8.3f   %+8.3f .. %+8.3f\n",
           16.5, 28.0, fold180b(sLo[0]), fold180b(sHi[0]),
           fold180b(sLo[0])+0.5*kPitchDeg, fold180b(sHi[0])-0.5*kPitchDeg);
    printf("  FTUS_2 outer  %4.1f-%4.1f   %+8.3f .. %+8.3f   %+8.3f .. %+8.3f\n",
           16.5, 28.0, fold180b(sLo[1]), fold180b(sHi[1]),
           fold180b(sLo[1])+0.5*kPitchDeg, fold180b(sHi[1])-0.5*kPitchDeg);

    int kk[4]; kk[0]=0; kk[1]=63; kk[2]=64; kk[3]=127;

    printf("\nStFstHitMaker INNER sensor\n");
    printf("  local[1] = phiInner + filp*dirct*k*pitch\n");
    printf("  phiStrip    global phi    LOCAL phi\n");
    for (int ka = 0; ka < 4; ka++){
        double gg = phIn + fl*dr*kk[ka]*kPitchDeg;
        printf("  %8d    %+9.4f    %+9.4f\n", kk[ka], fold180b(gg),
               fold180b(sgn*fold180b(gg - cOff)));
    }

    printf("\nStFstHitMaker OUTER, sensor==1 branch\n");
    printf("  local[1] = phiOuter - filp*dirct*k*pitch - filp*dirct*0.5*gap\n");
    printf("  phiStrip    global phi    LOCAL phi\n");
    for (int kb = 0; kb < 4; kb++){
        double gg2 = phOut - fl*dr*kk[kb]*kPitchDeg - fl*dr*0.5*kGapDeg;
        printf("  %8d    %+9.4f    %+9.4f\n", kk[kb], fold180b(gg2),
               fold180b(sgn*fold180b(gg2 - cOff)));
    }

    printf("\nStFstHitMaker OUTER, sensor==2 branch\n");
    printf("  local[1] = phiOuter - filp*dirct*k*pitch + filp*dirct*0.5*gap\n");
    printf("  phiStrip    global phi    LOCAL phi\n");
    for (int kc = 0; kc < 4; kc++){
        double gg3 = phOut - fl*dr*kk[kc]*kPitchDeg + fl*dr*0.5*kGapDeg;
        printf("  %8d    %+9.4f    %+9.4f\n", kk[kc], fold180b(gg3),
               fold180b(sgn*fold180b(gg3 - cOff)));
    }

    printf("\n-------------------------------------------------------------------------\n");
    printf("How to read it: a (branch, strip-range) pair is correct if its two\n");
    printf("endpoints land on an FTUS strip-centre pair listed above. Compare\n");
    printf("k=0..63 and k=64..127 separately against FTUS_1 and FTUS_2.\n");
    printf("-------------------------------------------------------------------------\n");
}
