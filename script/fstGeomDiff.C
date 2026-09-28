// fstGeomDiff.C -- mean FST sensor position in a geometry file, per disk and per half.
//
// Purpose: check whether two fGeom caches differ ONLY by the misalign tables, or
// also by the FST geometry TAG.  The blind-test baseline used fGeom.root (built
// 2026-06-27) while the closure run used dev2022sm, and the sensor hierarchy
// differs between tags (FSTD/FSTW nesting vs 108 FTUS lifted into FSTM), so the
// volume path cannot be assumed -- walk the tree and pick up FTUS* wherever it is.
//
// Reports the mean sensor centre, which is what a rigid misalign shifts, split by
// half (x>0 / x<0) because hssOnFst is half-antisymmetric.
//
// CINT: globals instead of reference args, unique loop variable names.

int    gN;
double gSx, gSy, gSz;                  // sums over all sensors
int    gNL, gNR;
double gSxL, gSyL, gSxR, gSyR;         // per half
int    gNd[3];
double gSxd[3], gSyd[3], gSzd[3];      // per disk

void fgdWalk(TGeoNode* node, TGeoHMatrix mat, int depth) {
    TGeoHMatrix here = mat;
    here.Multiply(node->GetMatrix());
    TString vname = node->GetVolume()->GetName();
    if (vname.BeginsWith("FTUS")) {
        const double* t = here.GetTranslation();
        gN++; gSx += t[0]; gSy += t[1]; gSz += t[2];
        if (t[0] >= 0) { gNR++; gSxR += t[0]; gSyR += t[1]; }
        else           { gNL++; gSxL += t[0]; gSyL += t[1]; }
        int d = -1;
        if      (t[2] < 158) d = 0;
        else if (t[2] < 172) d = 1;
        else                 d = 2;
        if (d >= 0) { gNd[d]++; gSxd[d] += t[0]; gSyd[d] += t[1]; gSzd[d] += t[2]; }
        return;                         // do not descend into the sensor
    }
    for (int i = 0; i < node->GetNdaughters(); i++)
        fgdWalk(node->GetDaughter(i), here, depth + 1);
}

void fstGeomDiff(const char* file) {
    gN = 0; gSx = gSy = gSz = 0;
    gNL = gNR = 0; gSxL = gSyL = gSxR = gSyR = 0;
    for (int k = 0; k < 3; k++) { gNd[k] = 0; gSxd[k] = gSyd[k] = gSzd[k] = 0; }

    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }
    TGeoManager* gm = (TGeoManager*)f->Get("dyson");
    if (!gm) { printf("no TGeoManager 'dyson' in %s\n", file); return; }

    // Start from FSTM, not the top node: walking the whole STAR tree takes minutes.
    // The path below FSTM differs between tags, so only the entry point is fixed.
    TGeoNavigator* nav = gm->GetCurrentNavigator();
    if (!nav) nav = gm->AddNavigator();
    if (!nav->cd("/HALL_1/CAVE_1/FSTM_1")) { printf("  no /HALL_1/CAVE_1/FSTM_1\n"); return; }
    TGeoHMatrix start(*nav->GetCurrentMatrix());
    TGeoNode* fstm = nav->GetCurrentNode();
    for (int i = 0; i < fstm->GetNdaughters(); i++)
        fgdWalk(fstm->GetDaughter(i), start, 0);

    printf("\n%s\n", file);
    if (gN == 0) { printf("  no FTUS sensors found\n"); return; }
    printf("  sensors found : %d\n", gN);
    printf("  mean centre   : x = %+9.5f  y = %+9.5f  z = %9.5f   (cm)\n",
           gSx/gN, gSy/gN, gSz/gN);
    if (gNR) printf("  right half    : x = %+9.5f  y = %+9.5f   (n=%d)\n", gSxR/gNR, gSyR/gNR, gNR);
    if (gNL) printf("  left  half    : x = %+9.5f  y = %+9.5f   (n=%d)\n", gSxL/gNL, gSyL/gNL, gNL);
    for (int kd = 0; kd < 3; kd++)
        if (gNd[kd]) printf("  disk %d        : x = %+9.5f  y = %+9.5f  z = %9.5f  (n=%d)\n",
                            kd, gSxd[kd]/gNd[kd], gSyd[kd]/gNd[kd], gSzd[kd]/gNd[kd], gNd[kd]);
}
