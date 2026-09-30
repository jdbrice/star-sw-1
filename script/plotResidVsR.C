// plotResidVsR.C -- the test that separates a per-quadrant OFFSET from a radial SCALE.
//
// A translation leaves a residual FLAT in hit radius. A scale leaves one PROPORTIONAL
// to it. The measured sTGC pattern pointed inward in all four quadrants, which sixteen
// translations can absorb without describing -- so this decides whether the constants
// we are fitting are alignment or curve-fitting.
//
// Run it on iteration-2 data, where the per-quadrant means are already ~0: any residual
// r-dependence left is the part the translations could not represent.
//
// CINT: globals instead of reference args, unique loop variable names (see CLAUDE.md).

double rvMu, rvSig, rvSignif, rvEnt;

void rvFit(TH1* h) {
    rvMu = 0; rvSig = 0; rvSignif = 0; rvEnt = 0;
    if (!h) return;
    rvEnt = h->GetEntries();
    if (rvEnt < 300) return;
    int b0 = h->GetXaxis()->FindBin(-2.9), b1 = h->GetXaxis()->FindBin(2.9);
    double best = -1; int bb = b0;
    for (int ib = b0 + 2; ib <= b1 - 2; ib++) {
        double v = h->GetBinContent(ib-1) + h->GetBinContent(ib) + h->GetBinContent(ib+1);
        if (v > best) { best = v; bb = ib; }
    }
    double mu0 = h->GetXaxis()->GetBinCenter(bb), lo = mu0 - 1.2, hi = mu0 + 1.2;
    TF1* fn = new TF1("rvF", "gaus(0)+pol1(3)", lo, hi);
    double base = 0.5 * (h->GetBinContent(h->GetXaxis()->FindBin(lo))
                       + h->GetBinContent(h->GetXaxis()->FindBin(hi)));
    fn->SetParameters(TMath::Max(1.0, h->GetBinContent(bb) - base), mu0, 0.4, base, 0);
    fn->SetParLimits(1, lo, hi);
    fn->SetParLimits(2, 0.05, 0.9);
    h->Fit(fn, "QNR");
    double A = fn->GetParameter(0), eA = fn->GetParError(0);
    rvMu = fn->GetParameter(1);
    rvSig = fn->GetParameter(2);
    rvSignif = (eA > 0) ? A / eA : 0;
}

// slice the 2D in radius, fit each slice, and fit a straight line to centre vs r.
// The slope is the answer: consistent with zero = offset, significantly non-zero = scale.
void rvPanel(TFile* f, const char* base, const char* title, const char* out, int isDx) {
    const int kNSlice = 6;
    const double rlo = 15, rhi = 63;                 // where the sTGC actually has hits
    double rC[kNSlice], rE[kNSlice], mC[kNSlice], mE[kNSlice];
    int n = 0;

    for (int is = 0; is < kNSlice; is++) {
        double r0 = rlo + (rhi - rlo) * is / kNSlice;
        double r1 = rlo + (rhi - rlo) * (is + 1) / kNSlice;
        TH1D* acc = 0;
        for (int id = 0; id < 4; id++) {             // sum the four stations
            TH2* h2 = (TH2*)f->Get(Form("%s_disk%d", base, id));
            if (!h2) continue;
            int bx0 = h2->GetXaxis()->FindBin(r0 + 1e-6);
            int bx1 = h2->GetXaxis()->FindBin(r1 - 1e-6);
            TH1D* p = h2->ProjectionY(Form("rv_%d_%d_%d", isDx, is, id), bx0, bx1);
            if (!acc) { acc = (TH1D*)p->Clone(Form("rvacc_%d_%d", isDx, is)); acc->SetDirectory(0); }
            else acc->Add(p);
        }
        rvFit(acc);
        if (rvSignif < 5) {
            printf("   r %4.1f-%4.1f : weak fit (%.0f entries, %.1f sigma) -- skipped\n",
                   r0, r1, rvEnt, rvSignif);
            continue;
        }
        rC[n] = 0.5 * (r0 + r1); rE[n] = 0.5 * (r1 - r0);
        mC[n] = rvMu;
        mE[n] = rvSig / sqrt(TMath::Max(1.0, rvEnt));   // error on the mean
        printf("   r %4.1f-%4.1f : centre %+7.4f +- %.4f cm   (%.0f entries, %.0f sigma)\n",
               r0, r1, mC[n], mE[n], rvEnt, rvSignif);
        n++;
    }
    if (n < 3) { printf("   too few usable slices to fit a slope\n"); return; }

    TGraphErrors* g = new TGraphErrors(n, rC, mC, rE, mE);
    TF1* lin = new TF1("lin", "pol1", rlo, rhi);
    g->Fit(lin, "QN");
    double slope = lin->GetParameter(1), eslope = lin->GetParError(1);
    double inter = lin->GetParameter(0);
    printf("   --> slope %+.5f +- %.5f cm per cm of radius  (%.1f sigma from zero)\n",
           slope, eslope, eslope > 0 ? fabs(slope)/eslope : 0);
    printf("   --> a pure SCALE of this size would be %+.3f%%\n", 100.0*slope);

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas(Form("crv%d", isDx), "crv", 620, 480);
    gPad->SetGridy();
    g->SetTitle(Form("%s;hit radius r [cm];fitted residual centre [cm]", title));
    g->SetMarkerStyle(20); g->SetMarkerSize(1.1);
    g->SetMarkerColor(isDx ? kAzure+2 : kRed+1);
    g->SetLineColor(isDx ? kAzure+2 : kRed+1);
    g->Draw("AP");
    lin->SetLineColor(kBlack); lin->SetLineStyle(2); lin->Draw("same");
    TLine* z = new TLine(rlo, 0, rhi, 0); z->SetLineColor(kGray+1); z->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.033);
    t.DrawLatex(0.15, 0.86, Form("slope = %+.4f #pm %.4f cm/cm  (%.1f#sigma)",
                                 slope, eslope, eslope > 0 ? fabs(slope)/eslope : 0));
    t.SetTextSize(0.029);
    t.DrawLatex(0.15, 0.81, (eslope > 0 && fabs(slope)/eslope > 3)
        ? "sloped: a radial SCALE, which per-quadrant offsets cannot represent"
        : "flat within errors: consistent with a genuine per-quadrant OFFSET");
    c->SaveAs(out);
}

void plotResidVsR(const char* file = "blinddiag_iter2.root") {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }
    printf("\n=== residual vs hit radius, %s ===\n", file);
    printf("\n dy (H strips):\n");
    rvPanel(f, "hBlindDyVsR_H", "dy vs radius", "FstFttFlipTest/plots/resid_vs_r_dy.png", 0);
    printf("\n dx (V strips):\n");
    rvPanel(f, "hBlindDxVsR_V", "dx vs radius", "FstFttFlipTest/plots/resid_vs_r_dx.png", 1);
}
