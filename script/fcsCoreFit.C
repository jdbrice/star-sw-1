// fcsCoreFit.C -- FCS residual core width WITHOUT mixed-event subtraction.
//
// Why: the mixed event is only a valid background model where same/mixed is flat. It is
// for Run22 (1.22 across 40-140 cm) but NOT for Run24 (1.18 -> 1.08, a 9% slope), so a
// single scale factor over/under-subtracts and the subtracted width becomes a property
// of the normalisation rather than of the data. Cross-dataset comparisons through
// subtraction are therefore unsafe; comparisons within one dataset are fine.
//
// Instead the same-event distribution is fitted directly with
//     narrow gaussian (the match) + broad gaussian (physics tail + combinatorics) + pol1
// so the background is measured per dataset instead of imported from another.
//
// fcCoreRaw()  fits the raw same-event histogram   -- use this across datasets
// fcCoreSub()  fits the mixed-subtracted histogram -- kept to cross-check where the
//              subtraction is valid; the two must agree there or the method is wrong.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

double fcCore, fcCoreE, fcBroad, fcMu, fcChi2;

void fcFitH(TH1* h, double win) {
    fcCore = 0; fcCoreE = 0; fcBroad = 0; fcMu = 0; fcChi2 = 0;
    if (!h) return;
    TH1D* d = (TH1D*)h->Clone(Form("fcw_%s", h->GetName()));
    d->SetDirectory(0);
    int pb = d->GetMaximumBin();
    double pk = d->GetXaxis()->GetBinCenter(pb);
    double amp = d->GetBinContent(pb);
    TF1* fn = new TF1(Form("fcf_%s", h->GetName()),
                      "gaus(0)+gaus(3)+pol1(6)", pk - win, pk + win);
    fn->SetParameters(amp * 0.25, pk, 2.5, amp * 0.75, pk, 20.0, amp * 0.1, 0);
    fn->SetParLimits(0, 0, amp * 5);
    fn->SetParLimits(1, pk - 6, pk + 6);
    fn->SetParLimits(2, 0.5, 8.0);        // the match
    fn->SetParLimits(3, 0, amp * 5);
    fn->SetParLimits(4, pk - 20, pk + 20);
    fn->SetParLimits(5, 9.0, 120.0);      // tail + combinatorics
    d->Fit(fn, "QNR");
    fcCore = fabs(fn->GetParameter(2)); fcCoreE = fn->GetParError(2);
    fcBroad = fabs(fn->GetParameter(5));
    fcMu = fn->GetParameter(1);
    fcChi2 = (fn->GetNDF() > 0) ? fn->GetChisquare() / fn->GetNDF() : 0;
}

TH1* fcGet(TFile* f, const char* panel, const char* xy, const char* type, int sub) {
    TH1* hs = (TH1*)f->Get(Form("%s%sSame Event_%s", panel, xy, type));
    if (!hs || !sub) return hs;
    TH1* hm = (TH1*)f->Get(Form("%s%sMixed Event_%s", panel, xy, type));
    if (!hm) return hs;
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= hs->GetNbinsX(); ib++) {
        double x = fabs(hs->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 200) continue;
        aS += hs->GetBinContent(ib); aM += hm->GetBinContent(ib);
    }
    TH1D* d = (TH1D*)hs->Clone(Form("fcs_%s%s_%s", panel, xy, type));
    d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("fcm_%s%s_%s", panel, xy, type));
    m->SetDirectory(0);
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    return d;
}

// one file, one track type, all four panels, both methods
void fcsCoreFit(const char* file, const char* label, const char* type = "Primary",
                double win = 50.0) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }
    const char* pn[4] = {"EcalNorth", "EcalSouth", "EcalTop", "EcalBottom"};
    const char* vx[4] = {"dx", "dx", "dy", "dy"};
    printf("\n=== %s  [%s] ===\n", label, type);
    printf("  %-12s %18s %8s | %18s %8s\n", "panel",
           "core, raw fit", "chi2/ndf", "core, subtracted", "chi2/ndf");
    for (int ip = 0; ip < 4; ip++) {
        fcFitH(fcGet(f, pn[ip], vx[ip], type, 0), win);
        double cr = fcCore, cre = fcCoreE, x2r = fcChi2;
        fcFitH(fcGet(f, pn[ip], vx[ip], type, 1), win);
        double cs = fcCore, cse = fcCoreE, x2s = fcChi2;
        printf("  %-12s %8.3f +-%6.3f %8.1f | %8.3f +-%6.3f %8.1f\n",
               pn[ip], cr, cre, x2r, cs, cse, x2s);
    }
}
