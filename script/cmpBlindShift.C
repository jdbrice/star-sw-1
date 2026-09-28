// cmpBlindShift.C -- per-station shift of the blind FST->sTGC residual between two
// geometries, to test the lever-arm prediction.
//
// The FST misalign table at 20211112 moves FST by dy = -3.170 mm, dx = +1.488 mm
// (fstOnTpc t1/t0).  The track is anchored at the TPC vertex, so at the sTGC the
// projection shift is amplified by z_stgc / z_eff with
//     z_eff = sum(z^2)/sum(z) over the three FST disks (151.75,165.25,178.78) = 166.0 cm
// Residual is hit - projection (FwdTracker.h:2366) and the FTT hit does NOT move
// (StFttDb is still hardcoded), so d(residual) = -d(projection).
//
// Fit is the same one that produced the tables (script/stgcQuadOffsets.C): gaus+pol1
// on a +-1.2 cm window around the tallest 3-bin run, fitted BEFORE any scaling.
//
// CINT: globals instead of reference args, unique loop variable names.

double gMu = 0, gSigma = 0, gSignif = 0, gEnt = 0;

void cbsFit(TH1* h) {
    gMu = 0; gSigma = 0; gSignif = 0; gEnt = 0;
    if (!h) return;
    gEnt = h->GetEntries();
    if (gEnt < 200) return;
    int b0 = h->GetXaxis()->FindBin(-2.9), b1 = h->GetXaxis()->FindBin(2.9);
    double best = -1; int bb = b0;
    for (int ib = b0 + 2; ib <= b1 - 2; ib++) {
        double v = h->GetBinContent(ib-1) + h->GetBinContent(ib) + h->GetBinContent(ib+1);
        if (v > best) { best = v; bb = ib; }
    }
    double mu0 = h->GetXaxis()->GetBinCenter(bb), lo = mu0 - 1.2, hi = mu0 + 1.2;
    TF1* fn = new TF1("cbsF", "gaus(0)+pol1(3)", lo, hi);
    double base = 0.5 * (h->GetBinContent(h->GetXaxis()->FindBin(lo))
                       + h->GetBinContent(h->GetXaxis()->FindBin(hi)));
    fn->SetParameters(TMath::Max(1.0, h->GetBinContent(bb) - base), mu0, 0.4, base, 0);
    fn->SetParLimits(1, lo, hi);
    fn->SetParLimits(2, 0.05, 0.9);
    h->Fit(fn, "QNR");
    double A = fn->GetParameter(0), eA = fn->GetParError(0);
    gMu = fn->GetParameter(1);
    gSigma = fn->GetParameter(2);
    gSignif = (eA > 0) ? A / eA : 0;
}

// filled by the two passes so the difference can be printed together
double muI[4][2], muM[4][2], sgI[4][2], sgM[4][2], snI[4][2], snM[4][2];

void cbsPass(const char* file, int which) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }
    for (int id = 0; id < 4; id++) {
        // coord 0 = dx from V strips (x precise), 1 = dy from H strips (y precise)
        for (int ic = 0; ic < 2; ic++) {
            // Use the GATED 1D "All" histograms.  Projecting the ungated 2D
            // ones instead was tried and is useless: the unconditioned landscape
            // is dominated by wrong-quadrant/wrong-row combinatorics (the "X" in
            // the 2D plots) and the fit finds no peak -- significance drops from
            // ~40 to <1.  The conditioning on the OTHER coordinate is what makes
            // these usable, and it leaves the shown coordinate unbiased
            // (FwdTracker.h:2330).  The gates are loose -- 3 sigma of the
            // imprecise coordinate, or 7.5 cm -- so a few-mm shift does not
            // meaningfully re-select.
            const char* nm = (ic == 0) ? Form("hBlindDxAll_V_disk%d", id)
                                       : Form("hBlindDyAll_H_disk%d", id);
            cbsFit((TH1*)f->Get(nm));
            if (which == 0) { muI[id][ic] = gMu; sgI[id][ic] = gSigma; snI[id][ic] = gSignif; }
            else            { muM[id][ic] = gMu; sgM[id][ic] = gSigma; snM[id][ic] = gSignif; }
        }
    }
}

void cmpBlindShift(const char* ideal = "blinddiag_ideal_20260925.root",
                   const char* mis   = "blinddiag_stgcmis_20260925.root") {
    cbsPass(ideal, 0);
    cbsPass(mis, 1);

    // prediction
    const double zFst = 166.0;                                   // sum(z^2)/sum(z)
    const double zSt[4] = {312.8385, 330.8385, 346.8385, 364.8385};
    const double dyFst = -0.31702, dxFst = 0.14876;              // cm, fstOnTpc

    printf("\n=== blind FST->sTGC residual, ideal vs misalign(20211112) ===\n");
    printf("fit: gaus+pol1, centres in cm; prediction = -dFst * z_station / %.1f\n\n", zFst);

    for (int ic = 0; ic < 2; ic++) {
        const double dFst = (ic == 0) ? dxFst : dyFst;
        printf("--- %s ---\n", (ic == 0) ? "dx  (V strips, x precise)" : "dy  (H strips, y precise)");
        printf("%-8s %9s %9s %9s %10s %9s  %7s %7s\n",
               "station", "ideal", "misalign", "measured", "predicted", "diff",
               "sig_id", "sig_mis");
        for (int id = 0; id < 4; id++) {
            double meas = muM[id][ic] - muI[id][ic];
            double pred = -dFst * zSt[id] / zFst;
            printf("%-8d %9.4f %9.4f %9.4f %10.4f %9.4f  %7.1f %7.1f\n",
                   id + 1, muI[id][ic], muM[id][ic], meas, pred, meas - pred,
                   snI[id][ic], snM[id][ic]);
        }
        printf("\n");
    }

    printf("--- widths (cm), should be unchanged: a rigid shift does not broaden ---\n");
    printf("%-8s %8s %8s   %8s %8s\n", "station", "sig_dx_i", "sig_dx_m", "sig_dy_i", "sig_dy_m");
    for (int id = 0; id < 4; id++)
        printf("%-8d %8.4f %8.4f   %8.4f %8.4f\n",
               id + 1, sgI[id][0], sgM[id][0], sgI[id][1], sgM[id][1]);
}
