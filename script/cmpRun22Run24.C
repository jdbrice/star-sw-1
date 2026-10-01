// cmpRun22Run24.C -- side-by-side Run22 vs Run24 for the two quantities the Run24
// campaign exists to measure.
//
// Both years are run with the SAME configuration on purpose: calibrated time cut
// (mode 2) with per-run anchors, aligned geometry, DIAG type=2 (Primary). Run24 starts
// from the Run22 final constants, so:
//
//   ALIGNMENT  the Run24 per-quadrant residual IS the Run22 -> Run24 movement, read
//              directly. Expect mm, not cm: the sTGC went back on the same rails.
//              A cm-level number means suspect the carry-over or the quadrant
//              convention before concluding the detector moved.
//   PICKUP     measured on aligned geometry in both years, so a difference is
//              detector//electronics, not geometry. Run24 sTGC is much hotter
//              (FttRawHit 1239/ev vs 326, FttPoint 107/ev vs 2.5).
//
// Usage:
//   root -l -b -q 'cmpRun22Run24.C("run22_Primary.root","run24_Primary.root")'
//
// CINT: unique loop variable names, no std::map (see CLAUDE.md).

double c22x[4], c22y[4], c24x[4], c24y[4];

void cmpPick(const char* fn, double* px, double* py, double& ntr) {
    for (int i = 0; i < 4; i++) { px[i] = 0; py[i] = 0; }
    ntr = 0;
    TFile* f = TFile::Open(fn);
    if (!f || f->IsZombie()) { printf("   cannot open %s\n", fn); return; }
    TH1F* h = (TH1F*)f->Get("PlaneUsage/hPlaneUsage_Primary");
    if (!h) { printf("   no hPlaneUsage_Primary in %s\n", fn); return; }
    ntr = h->GetBinContent(20);
    if (ntr <= 0) return;
    for (int ip = 0; ip < 4; ip++) {
        px[ip] = h->GetBinContent(4 + 3*ip + 1) / ntr;
        py[ip] = h->GetBinContent(4 + 3*ip + 2) / ntr;
    }
}

void cmpRun22Run24(const char* f22 = "fttadd_Primary.root",
                   const char* f24 = "run24pickup_Primary.root") {
    double n22 = 0, n24 = 0;
    cmpPick(f22, c22x, c22y, n22);
    cmpPick(f24, c24x, c24y, n24);

    printf("\n=== FTT pickup probability, Primary tracks (raw, all-track denominator) ===\n");
    printf("  tracks: Run22 %.0f   Run24 %.0f\n", n22, n24);
    printf("  %-8s %18s %18s %10s\n", "station", "Run22 x / y", "Run24 x / y", "ratio x/y");
    double s22 = 0, s24 = 0;
    for (int is = 0; is < 4; is++) {
        s22 += c22x[is] + c22y[is];
        s24 += c24x[is] + c24y[is];
        printf("    %d    %8.4f /%8.4f  %8.4f /%8.4f   %.2f /%.2f\n", is + 1,
               c22x[is], c22y[is], c24x[is], c24y[is],
               (c22x[is] > 0) ? c24x[is]/c22x[is] : 0,
               (c22y[is] > 0) ? c24y[is]/c22y[is] : 0);
    }
    printf("  hits/track over 4 stations x 2 orientations:  Run22 %.3f   Run24 %.3f   ratio %.2f\n",
           s22, s24, (s22 > 0) ? s24/s22 : 0);
    printf("\n  NOTE: this denominator is ALL tracks, so it includes the sTGC design\n");
    printf("  acceptance gap below the beampipe. For the acceptance-corrected number run\n");
    printf("  script/pickupVsAcceptance.C on each year -- it needs h2FttProjAllXY, which\n");
    printf("  exists in the Run24 campaign but in Run22 only in the short acc_Primary run.\n");
    printf("\n  For the alignment side run script/stgcQuadOffsets.C on the Run24 blinddiag:\n");
    printf("  with Run24 seeded from the Run22 finals, those residuals ARE the movement.\n");
}
