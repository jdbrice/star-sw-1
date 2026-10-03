// fcsResidCore.C -- the FCS residual core, now that picoMatch bins dx/dy at 0.5 cm.
//
// Two jobs:
//  (1) scan the fit window, so we can see whether sigma has stopped tracking the window
//      (at the old 4 cm binning it never did: 9.4 cm at w=60 down to 3.6 cm at w=8);
//  (2) compare the two arms at fixed windows -- FTT hits off vs aligned FTT hits on.
//
// Mixed event is normalised to same over |d| in 80-200 cm and subtracted, as in
// cmpFcsResid.C. CINT: unique loop names, fixed arrays (see CLAUDE.md).

double fcSig, fcMu, fcErrSig, fcSignif, fcInt;

void fcFit(TH1* h, double w) {
    fcSig = 0; fcMu = 0; fcErrSig = 0; fcSignif = 0;
    if (!h) return;
    int pb = h->GetMaximumBin();
    double pk = h->GetXaxis()->GetBinCenter(pb);
    TF1* fn = new TF1(Form("fc_%s_%g", h->GetName(), w), "gaus(0)+pol1(3)", pk - w, pk + w);
    fn->SetParameters(h->GetBinContent(pb), pk, w / 3.0, 0, 0);
    fn->SetParLimits(2, 0.05, 80.0);
    h->Fit(fn, "QNR");
    fcSig = fabs(fn->GetParameter(2));
    fcMu = fn->GetParameter(1);
    fcErrSig = fn->GetParError(2);
    fcSignif = (fn->GetParError(0) > 0) ? fn->GetParameter(0) / fn->GetParError(0) : 0;
}

TH1D* fcSub(TFile* f, const char* panel, const char* xy) {
    TH1* hs = (TH1*)f->Get(Form("%s%sSame Event_Primary", panel, xy));
    TH1* hm = (TH1*)f->Get(Form("%s%sMixed Event_Primary", panel, xy));
    if (!hs || !hm) return 0;
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= hs->GetNbinsX(); ib++) {
        double x = fabs(hs->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 200) continue;
        aS += hs->GetBinContent(ib); aM += hm->GetBinContent(ib);
    }
    TH1D* d = (TH1D*)hs->Clone(Form("sub_%s%s_%p", panel, xy, (void*)f));
    d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("scl_%s%s_%p", panel, xy, (void*)f));
    m->SetDirectory(0);
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    fcInt = d->Integral();
    return d;
}

void fcsResidCore(const char* fOff, const char* fOn) {
    TFile* a = TFile::Open(fOff);
    TFile* b = TFile::Open(fOn);
    if (!a || a->IsZombie() || !b || b->IsZombie()) { printf("cannot open inputs\n"); return; }
    const char* pn[4] = {"EcalNorth", "EcalSouth", "EcalTop", "EcalBottom"};
    const char* vx[4] = {"dx", "dx", "dy", "dy"};
    double ws[6] = {30, 20, 10, 5, 3, 2};

    TH1D* probe = fcSub(a, pn[0], vx[0]);
    printf("\n=== FCS residual core, Primary tracks, Run24 iter-1, bin width %.2f cm ===\n",
           probe ? probe->GetXaxis()->GetBinWidth(1) : 0.0);

    printf("\n  (1) fit-window scan, arm A (FTT hits OFF) -- sigma in cm\n");
    printf("  %-12s %8s %8s %8s %8s %8s %8s\n", "panel", "w=30", "w=20", "w=10", "w=5", "w=3", "w=2");
    for (int ip = 0; ip < 4; ip++) {
        TH1D* d = fcSub(a, pn[ip], vx[ip]);
        printf("  %-12s", pn[ip]);
        for (int jw = 0; jw < 6; jw++) { fcFit(d, ws[jw]); printf(" %8.3f", fcSig); }
        printf("\n");
    }

    printf("\n  (2) arm A vs arm B (aligned FTT hits ON)\n");
    printf("  %-12s %16s %16s %8s %7s\n", "panel", "A sigma", "B sigma", "B/A", "n_sig");
    double wUse[3] = {10, 5, 3};
    for (int kw = 0; kw < 3; kw++) {
        printf("   -- fit window +-%.0f cm --\n", wUse[kw]);
        for (int iq = 0; iq < 4; iq++) {
            TH1D* da = fcSub(a, pn[iq], vx[iq]); fcFit(da, wUse[kw]);
            double sa = fcSig, ea = fcErrSig, sga = fcSignif;
            TH1D* db = fcSub(b, pn[iq], vx[iq]); fcFit(db, wUse[kw]);
            double sb = fcSig, eb = fcErrSig;
            double de = sqrt(ea * ea + eb * eb);
            if (sa > 0 && sb > 0 && sga > 3)
                printf("  %-12s %8.3f +-%5.3f %8.3f +-%5.3f %8.3f %7.1f\n",
                       pn[iq], sa, ea, sb, eb, sb / sa, de > 0 ? (sb - sa) / de : 0);
            else
                printf("  %-12s  weak peak (%.1f sigma)\n", pn[iq], sga);
        }
    }
    printf("\n  A stable sigma across windows means the core is now resolved.\n");
}

// ---------------------------------------------------------------------------
// A single Gaussian never stops sliding with the window (see the scan above), so
// the shape is a narrow core on a wide tail. Fitting BOTH over one wide window
// gives a core sigma that does not depend on where we cut -- use this one for
// tracking alignment changes, not the single-Gaussian number.
double dgC, dgCe, dgT, dgF, dgMu;

void dgFit(TH1* h, double w) {
    dgC = 0; dgCe = 0; dgT = 0; dgF = 0; dgMu = 0;
    if (!h) return;
    int pb = h->GetMaximumBin();
    double pk = h->GetXaxis()->GetBinCenter(pb);
    TF1* fn = new TF1(Form("dg_%s", h->GetName()),
                      "gaus(0)+gaus(3)+pol1(6)", pk - w, pk + w);
    fn->SetParameters(h->GetBinContent(pb) * 0.7, pk, 1.5,
                      h->GetBinContent(pb) * 0.3, pk, 8.0, 0, 0);
    fn->SetParLimits(1, pk - 5, pk + 5);
    fn->SetParLimits(2, 0.3, 4.0);        // core
    fn->SetParLimits(4, pk - 8, pk + 8);
    fn->SetParLimits(5, 4.0, 40.0);       // tail
    h->Fit(fn, "QNR");
    dgC = fabs(fn->GetParameter(2)); dgCe = fn->GetParError(2);
    dgT = fabs(fn->GetParameter(5));
    dgMu = fn->GetParameter(1);
    double aC = fn->GetParameter(0) * dgC, aT = fn->GetParameter(3) * dgT;
    dgF = (aC + aT > 0) ? aC / (aC + aT) : 0;
}

TH1D* dgSub(TFile* f, const char* panel, const char* xy) {
    TH1* hs = (TH1*)f->Get(Form("%s%sSame Event_Primary", panel, xy));
    TH1* hm = (TH1*)f->Get(Form("%s%sMixed Event_Primary", panel, xy));
    if (!hs || !hm) return 0;
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= hs->GetNbinsX(); ib++) {
        double x = fabs(hs->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 200) continue;
        aS += hs->GetBinContent(ib); aM += hm->GetBinContent(ib);
    }
    TH1D* d = (TH1D*)hs->Clone(Form("dsub_%s%s_%p", panel, xy, (void*)f));
    d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("dscl_%s%s_%p", panel, xy, (void*)f));
    m->SetDirectory(0);
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    return d;
}

void fcsResidCoreTail(const char* fa, const char* fb) {
    TFile* a = TFile::Open(fa); TFile* b = TFile::Open(fb);
    const char* pn[4] = {"EcalNorth", "EcalSouth", "EcalTop", "EcalBottom"};
    const char* vx[4] = {"dx", "dx", "dy", "dy"};
    printf("\n  core+tail fit over +-25 cm, Primary, 0.5 cm bins\n");
    printf("  %-12s %18s %7s %7s | %18s %7s %7s | %6s\n",
           "panel", "A core", "A tail", "A frac", "B core", "B tail", "B frac", "B/A");
    for (int ip = 0; ip < 4; ip++) {
        dgFit(dgSub(a, pn[ip], vx[ip]), 25.0);
        double ca = dgC, cae = dgCe, ta = dgT, fa2 = dgF;
        dgFit(dgSub(b, pn[ip], vx[ip]), 25.0);
        double cb = dgC, cbe = dgCe, tb = dgT, fb2 = dgF;
        printf("  %-12s %8.3f +-%6.3f %7.2f %7.3f | %8.3f +-%6.3f %7.2f %7.3f | %6.3f\n",
               pn[ip], ca, cae, ta, fa2, cb, cbe, tb, fb2, ca > 0 ? cb / ca : 0);
    }
    printf("\n  core/tail sigma in cm; frac = core area fraction\n");
}
