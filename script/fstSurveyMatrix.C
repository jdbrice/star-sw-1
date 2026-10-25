// fstSurveyMatrix.C -- what does the FST survey matrix that StFstHitMaker would apply
// actually contain, at a given timestamp?
//
//   root4star -l -b -q 'script/fstSurveyMatrix.C("20211112.000004")'
//
// StFstDb composes, per sensor,
//     R%04i = tpcOnGlobal * fstOnTpc * hssOnFst[half] * fstWedgeOnHss[wedge] * fstSensorOnWedge[id]
// and StFstHitMaker::SetApplySurvey(true) applies exactly that to the nominal hit
// position. The FST tables are translations only, so anything else in R -- a rotation,
// or a translation that varies smoothly with azimuth -- comes from tpcOnGlobal, which
// is the TPC's own alignment and has nothing to do with the FST survey.
//
// Prints each factor separately, then the composed matrix for one sensor per wedge, so
// a per-wedge pattern is visible as a pattern rather than as a single mean.
//
// CINT: unique loop-variable names, no local std::map.
void fstSurveyMatrix(const char* sdt = "20211112.000004",
                     const char* tag = "dev2022sm")
{
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");
    TString ch = Form("sdt%s %s agml usexgeom misalign db tpcDb fstDb", sdt, tag);
    printf("fstSurveyMatrix: chain = [%s]\n", ch.Data());
    bfc(-1, ch, 0);

    StChain* chain = (StChain*)StMaker::GetChain();
    if (chain->Init() != kStOK) { printf("fstSurveyMatrix: Init FAILED\n"); return; }
    chain->Make();

    StFstDb* db = 0;
    TDataSet* ds = chain->GetDataSet("fst_db");
    if (ds) db = (StFstDb*)ds->GetObject();
    if (!db) {
        TObjectSet* os = (TObjectSet*)chain->GetDataSet("fstDb");
        if (os) db = (StFstDb*)os->GetObject();
    }
    if (!db) { printf("fstSurveyMatrix: no StFstDb in the chain\n"); return; }

    printf("\n--- tpcOnGlobal (gStTpcDb, NOT part of the FST survey) ---\n");
    const TGeoHMatrix* mt = db->getGeoHMatrixTpcOnGlobal();
    if (mt) {
        const Double_t* tt = mt->GetTranslation();
        const Double_t* rr = mt->GetRotationMatrix();
        printf("  t = (%+10.6f, %+10.6f, %+10.6f) cm\n", tt[0], tt[1], tt[2]);
        printf("  R = [%+10.7f %+10.7f %+10.7f]\n", rr[0], rr[1], rr[2]);
        printf("      [%+10.7f %+10.7f %+10.7f]\n", rr[3], rr[4], rr[5]);
        printf("      [%+10.7f %+10.7f %+10.7f]\n", rr[6], rr[7], rr[8]);
        printf("  off-diagonal (a z-rotation would sit in R01/R10): %+.7f %+.7f\n", rr[1], rr[3]);
    }

    printf("\n--- the four FST survey tables as StFstDb loaded them ---\n");
    const TGeoHMatrix* mf = db->getGeoHMatrixFstOnTpc();
    printf("  fstOnTpc        t = (%+9.5f, %+9.5f, %+9.5f)\n",
           mf->GetTranslation()[0], mf->GetTranslation()[1], mf->GetTranslation()[2]);
    for (int ih = 0; ih < 2; ih++) {
        const TGeoHMatrix* mh = db->getGeoHMatrixHssOnFst(ih);
        printf("  hssOnFst[%d]     t = (%+9.5f, %+9.5f, %+9.5f)\n", ih,
               mh->GetTranslation()[0], mh->GetTranslation()[1], mh->GetTranslation()[2]);
    }
    int nWedgeNonId = 0;
    for (int iw = 1; iw <= 36; iw++) {
        const TGeoHMatrix* mw = db->getGeoHMatrixWedgeOnHss(iw);
        const Double_t* tw = mw->GetTranslation();
        if (tw[0] || tw[1] || tw[2]) {
            printf("  fstWedgeOnHss[%2d] t = (%+9.5f, %+9.5f, %+9.5f)\n", iw, tw[0], tw[1], tw[2]);
            nWedgeNonId++;
        }
    }
    printf("  fstWedgeOnHss: %d of 36 rows non-zero\n", nWedgeNonId);
    int nSensNonId = 0;
    for (int jw = 1; jw <= 36; jw++) {
        for (int js = 0; js < 3; js++) {
            const TGeoHMatrix* ms = db->getGeoHMatrixSensorOnWedge(jw, js);
            const Double_t* ts = ms->GetTranslation();
            if (ts[0] || ts[1] || ts[2]) nSensNonId++;
        }
    }
    printf("  fstSensorOnWedge: %d of 108 rows non-zero\n", nSensNonId);

    printf("\n--- composed R%%04i, one sensor (inner) per wedge ---\n");
    printf("  wedge half   t_x        t_y        t_z     | R01       R10      (z-rot)\n");
    for (int kw = 1; kw <= 36; kw++) {
        const TGeoHMatrix* mc = StFstDb::getHMatrixSensorOnGlobal(kw, 0);
        if (!mc) { printf("  %5d   -- no matrix\n", kw); continue; }
        const Double_t* tc = mc->GetTranslation();
        const Double_t* rc = mc->GetRotationMatrix();
        printf("  %5d %4d %+10.6f %+10.6f %+10.6f | %+9.6f %+9.6f\n",
               kw, StFstDb::halfId(kw), tc[0], tc[1], tc[2], rc[1], rc[3]);
    }

    printf("\n--- what applying it does to a point at the nominal hit radius ---\n");
    printf("  A translation moves every wedge the same way; a rotation moves them\n");
    printf("  tangentially, so dx/dy rotates with the wedge azimuth.\n");
    printf("  wedge   phi[deg]   nominal(x,y)        ->  applied(x,y)        dx       dy\n");
    for (int lw = 1; lw <= 12; lw++) {
        double phd = 75.0 - 30.0 * (lw - 1);
        double rr2 = 15.0;
        double l3[3], g3[3];
        l3[0] = rr2 * TMath::Cos(phd * TMath::DegToRad());
        l3[1] = rr2 * TMath::Sin(phd * TMath::DegToRad());
        l3[2] = StFstDb::sensorZ(1, lw, 0);
        if (!StFstDb::sensorLocalToGlobal(lw, 0, l3, g3)) continue;
        printf("  %5d  %+8.1f  (%+8.4f,%+8.4f) -> (%+8.4f,%+8.4f)  %+8.5f %+8.5f\n",
               lw, phd, l3[0], l3[1], g3[0], g3[1], g3[0] - l3[0], g3[1] - l3[1]);
    }
}
