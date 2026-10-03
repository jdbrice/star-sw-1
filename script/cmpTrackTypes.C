// cmpTrackTypes.C -- the FCS residual for every track type inside ONE picoMatch file,
// so the track-type change is isolated with no other variable moving.
//
// Reports a truncated RMS over |d| < 30 cm, which is comparable across files with
// different binning (old output is 4 cm, new is 0.5 cm), and a core+tail core sigma,
// which is only meaningful where the bins resolve it.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

double ttRms, ttCore, ttCoreE, ttInt, ttBw;

TH1D* ttSub(TFile* f, const char* panel, const char* xy, const char* type) {
    TH1* hs = (TH1*)f->Get(Form("%s%sSame Event_%s", panel, xy, type));
    TH1* hm = (TH1*)f->Get(Form("%s%sMixed Event_%s", panel, xy, type));
    if (!hs || !hm) return 0;
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= hs->GetNbinsX(); ib++) {
        double x = fabs(hs->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 200) continue;
        aS += hs->GetBinContent(ib); aM += hm->GetBinContent(ib);
    }
    TH1D* d = (TH1D*)hs->Clone(Form("tt_%s%s_%s", panel, xy, type));
    d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("ttm_%s%s_%s", panel, xy, type));
    m->SetDirectory(0);
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    ttBw = d->GetXaxis()->GetBinWidth(1);
    return d;
}

void ttMeasure(TH1D* d) {
    ttRms = 0; ttCore = 0; ttCoreE = 0; ttInt = 0;
    if (!d) return;
    ttInt = d->Integral();
    TH1D* q = (TH1D*)d->Clone("ttq"); q->SetDirectory(0);
    q->GetXaxis()->SetRangeUser(-30, 30);
    ttRms = q->GetRMS();
    int pb = d->GetMaximumBin();
    double pk = d->GetXaxis()->GetBinCenter(pb);
    TF1* fn = new TF1("ttf", "gaus(0)+gaus(3)+pol1(6)", pk - 25, pk + 25);
    fn->SetParameters(d->GetBinContent(pb) * 0.7, pk, 1.5,
                      d->GetBinContent(pb) * 0.3, pk, 8.0, 0, 0);
    fn->SetParLimits(1, pk - 5, pk + 5);
    // wide enough that the broad track types do not rail; a value sitting exactly
    // on a limit is flagged with * below and should be read as "fit failed, use rms"
    fn->SetParLimits(2, 0.3, 8.0);
    fn->SetParLimits(4, pk - 10, pk + 10);
    fn->SetParLimits(5, 8.0, 60.0);
    d->Fit(fn, "QNR");
    ttCore = fabs(fn->GetParameter(2)); ttCoreE = fn->GetParError(2);
}

void cmpTrackTypes(const char* file, const char* label) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }
    const char* tt[6] = {"Global", "Beamline", "Primary", "FwdVtx", "BLCVtx", "FCSTRK"};
    const char* pn[4] = {"EcalNorth", "EcalSouth", "EcalTop", "EcalBottom"};
    const char* vx[4] = {"dx", "dx", "dy", "dy"};
    printf("\n=== %s ===\n", label);
    printf("  %-10s %22s %22s %22s %22s\n", "type",
           "North dx rms/core", "South dx rms/core", "Top dy rms/core", "Bottom dy rms/core");
    for (int it = 0; it < 6; it++) {
        printf("  %-10s", tt[it]);
        for (int ip = 0; ip < 4; ip++) {
            TH1D* d = ttSub(f, pn[ip], vx[ip], tt[it]);
            if (!d || d->Integral() <= 0) { printf("%23s", "--"); continue; }
            ttMeasure(d);
            bool rail = (ttCore < 0.31 || ttCore > 7.99);
            printf("   %8.2f /%7.2f%s", ttRms, ttCore, rail ? "*" : " ");
        }
        printf("\n");
    }
    printf("  rms over |d|<30 cm; core from a core+tail fit (bin width %.2f cm)\n", ttBw);
}
