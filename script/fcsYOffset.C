// fcsYOffset.C -- is the north/south r*dPhi antisymmetry a relative y offset?
//
// Geometry of the test. If the FCS clusters sit at y + dy relative to where the tracks
// project, then at phi~0 (North, +x) the azimuthal direction is +y so r*dPhi ~ +dy, and
// at phi~pi (South, -x) it is -y so r*dPhi ~ -dy. The Cartesian dy residual, by contrast,
// is the same +dy in BOTH halves, and dx is untouched. So:
//     r*dPhi  antisymmetric   north = -south = +dy   (sign per the above)
//     dy      COMMON shift    top  =  bottom = +dy
//     dx      no shift
// Anything else is not a pure y offset.
//
// Second question: is it FCS-side or track-side? Run it per track type inside ONE file.
// A detector offset is common to all types; a track-side bias follows the vertex
// constraint and differs between them.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

double yoMean(TFile* f, const char* panel, const char* xy, const char* type, int rebin) {
    TH1* hs = (TH1*)f->Get(Form("%s%sSame Event_%s", panel, xy, type));
    TH1* hm = (TH1*)f->Get(Form("%s%sMixed Event_%s", panel, xy, type));
    if (!hs || !hm) return -999;
    TH1D* d = (TH1D*)hs->Clone(Form("yo_%s%s_%s", panel, xy, type)); d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("yom_%s%s_%s", panel, xy, type)); m->SetDirectory(0);
    if (rebin > 1) { d->Rebin(rebin); m->Rebin(rebin); }
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= d->GetNbinsX(); ib++) {
        double x = fabs(d->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 160) continue;
        aS += d->GetBinContent(ib); aM += m->GetBinContent(ib);
    }
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    d->GetXaxis()->SetRangeUser(-30, 30);
    return d->GetMean();
}

void fcsYOffset(const char* file, const char* label, int rebin = 8) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }
    const char* tt[6] = {"Global", "Beamline", "Primary", "FwdVtx", "BLCVtx", "FCSTRK"};
    printf("\n=== %s ===\n", label);
    printf("  %-9s | %17s | %17s | %17s | %8s\n", "type",
           "r*dPhi N / S", "dy Top / Bottom", "dx North / South", "dy from");
    printf("  %-9s | %17s | %17s | %17s | %8s\n", "", "(antisym => y off)",
           "(common => y off)", "(should be flat)", "r*dPhi");
    for (int it = 0; it < 6; it++) {
        double pn = yoMean(f, "EcalNorth", "dp", tt[it], rebin);
        double ps = yoMean(f, "EcalSouth", "dp", tt[it], rebin);
        double yt = yoMean(f, "EcalTop",    "dy", tt[it], rebin);
        double yb = yoMean(f, "EcalBottom", "dy", tt[it], rebin);
        double xn = yoMean(f, "EcalNorth", "dx", tt[it], rebin);
        double xs = yoMean(f, "EcalSouth", "dx", tt[it], rebin);
        if (pn < -900 || yt < -900) { printf("  %-9s | %17s\n", tt[it], "missing"); continue; }
        double dyFromPhi = 0.5 * (pn - ps);   // north minus south over two
        printf("  %-9s | %7.2f %8.2f | %7.2f %8.2f | %7.2f %8.2f | %8.2f\n",
               tt[it], pn, ps, yt, yb, xn, xs, dyFromPhi);
    }
    printf("  all in cm, mixed-event subtracted, |d|<30 cm\n");
}
