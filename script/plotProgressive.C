// plotProgressive.C -- the three-stage progression asked for on slide 2 of
// 20260929_fwd_tracking.pdf:
//
//    1. spacepoint mode          kUseSpacePoints = true, StFttDb hardcoded
//    2. planar mode              kUseSpacePoints = false, StFttDb hardcoded
//    3. planar + FTT tables in   kUseSpacePoints = false, StFttDb from the DB survey
//
// All three are the SAME 2500 events of the same file, so the comparison carries no
// statistics or sample difference -- only the thing being changed. Stage 1 needed a
// rebuild because kUseSpacePoints is a compile-time constant.
//
// Left panel dx from V strips (x is the precise coordinate), right panel dy from H
// strips. Curves are sideband-normalised so three samples with different yields can be
// overlaid; the fitted peak centre is what to read, not the height.
//
// CINT: no reference args, globals instead, unique loop variable names (see CLAUDE.md).

double ppMu[3], ppSig[3], ppSignif[3], ppEnt[3];

// gaus+pol1 on a +-1.2 cm window around the tallest 3-bin run; fit BEFORE any scaling,
// since Integral() after Scale() is not an entry count (that bug cost us a day once).
void ppFit(TH1* h, int islot) {
    ppMu[islot] = 0; ppSig[islot] = 0; ppSignif[islot] = 0; ppEnt[islot] = 0;
    if (!h) return;
    ppEnt[islot] = h->GetEntries();
    if (ppEnt[islot] < 200) return;
    int b0 = h->GetXaxis()->FindBin(-2.9), b1 = h->GetXaxis()->FindBin(2.9);
    double best = -1; int bb = b0;
    for (int ib = b0 + 2; ib <= b1 - 2; ib++) {
        double v = h->GetBinContent(ib-1) + h->GetBinContent(ib) + h->GetBinContent(ib+1);
        if (v > best) { best = v; bb = ib; }
    }
    double mu0 = h->GetXaxis()->GetBinCenter(bb), lo = mu0 - 1.2, hi = mu0 + 1.2;
    TF1* fn = new TF1(Form("ppF%d", islot), "gaus(0)+pol1(3)", lo, hi);
    double base = 0.5 * (h->GetBinContent(h->GetXaxis()->FindBin(lo))
                       + h->GetBinContent(h->GetXaxis()->FindBin(hi)));
    fn->SetParameters(TMath::Max(1.0, h->GetBinContent(bb) - base), mu0, 0.4, base, 0);
    fn->SetParLimits(1, lo, hi);
    fn->SetParLimits(2, 0.05, 0.9);
    h->Fit(fn, "QNR");
    double A = fn->GetParameter(0), eA = fn->GetParError(0);
    ppMu[islot] = fn->GetParameter(1);
    ppSig[islot] = fn->GetParameter(2);
    ppSignif[islot] = (eA > 0) ? A / eA : 0;
}

// sum the four stations, then normalise on the sidebands so the three overlay
TH1F* ppSum(TFile* f, const char* base, const char* tag) {
    if (!f) return 0;
    TH1F* out = 0;
    for (int id = 0; id < 4; id++) {
        TH1F* h = (TH1F*)f->Get(Form("%s_disk%d", base, id));
        if (!h) continue;
        if (!out) { out = (TH1F*)h->Clone(Form("pp_%s_%s", base, tag)); out->SetDirectory(0); }
        else out->Add(h);
    }
    return out;
}

void ppPanel(TFile** f, const char** lab, const char* base, const char* title) {
    const int col[3] = { kGray+2, kAzure+2, kRed+1 };
    TH1F* h[3]; double ymax = 0;
    for (int k = 0; k < 3; k++) {
        h[k] = ppSum(f[k], base, Form("s%d", k));
        ppFit(h[k], k);                         // fit first
        if (!h[k]) continue;
        // sideband normalisation: mean content outside |d| > 3 cm
        double sb = 0; int nsb = 0;
        for (int ib = 1; ib <= h[k]->GetNbinsX(); ib++) {
            double x = h[k]->GetXaxis()->GetBinCenter(ib);
            if (fabs(x) > 3.0) { sb += h[k]->GetBinContent(ib); nsb++; }
        }
        if (nsb && sb > 0) h[k]->Scale(nsb / sb);
        h[k]->SetLineColor(col[k]); h[k]->SetLineWidth(2);
        h[k]->SetTitle(Form("%s;d [cm];sideband-normalised", title));
        if (h[k]->GetMaximum() > ymax) ymax = h[k]->GetMaximum();
    }
    for (int k = 0; k < 3; k++) {
        if (!h[k]) continue;
        h[k]->GetYaxis()->SetRangeUser(0, ymax * 1.45);
        h[k]->GetXaxis()->SetRangeUser(-5, 5);   // the peaks live here; +-15 wastes the frame
        h[k]->Draw(k ? "hist same" : "hist");
    }
    TLegend* lg = new TLegend(0.13, 0.72, 0.62, 0.89);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.031);
    for (int k = 0; k < 3; k++)
        if (h[k]) lg->AddEntry(h[k], Form("%s  %+.2f cm (%.0f#sigma)",
                                          lab[k], ppMu[k], ppSignif[k]), "l");
    lg->Draw();
    TLine* z = new TLine(0, 0, 0, ymax * 1.45);
    z->SetLineStyle(2); z->SetLineColor(kBlack); z->Draw();
}

void plotProgressive(const char* fsp  = "blinddiag_AB_spacepoint.root",
                     const char* fpl  = "blinddiag_AB_dbOFF.root",
                     const char* fdb  = "blinddiag_AB_dbON.root",
                     const char* out  = "FstFttFlipTest/plots/progressive.png") {
    gStyle->SetOptStat(0);
    TFile* f[3];
    f[0] = TFile::Open(fsp);
    f[1] = TFile::Open(fpl);
    f[2] = TFile::Open(fdb);
    const char* lab[3] = { "1 spacepoint", "2 planar", "3 planar + FTT tables" };

    TCanvas* c = new TCanvas("cpp", "cpp", 1100, 480);
    c->Divide(2, 1);
    c->cd(1); gPad->SetGridx();
    ppPanel(f, lab, "hBlindDxAll_V", "V strips:  dx = x_{hit} - x_{proj}");
    double mux[3]; for (int i = 0; i < 3; i++) mux[i] = ppMu[i];
    c->cd(2); gPad->SetGridx();
    ppPanel(f, lab, "hBlindDyAll_H", "H strips:  dy = y_{hit} - y_{proj}");

    c->cd(1);
    TLatex t; t.SetNDC(); t.SetTextSize(0.028);
    t.DrawLatex(0.13, 0.66, "same 2500 events in all three");
    c->SaveAs(out);

    printf("\n=== three-stage progression, fitted peak centres [cm] ===\n");
    printf("  %-24s %10s %10s %10s %10s\n", "stage", "dx", "dy", "sig_dy", "entries_dy");
    // refit cleanly for the printout
    for (int k = 0; k < 3; k++) {
        TH1F* hx = ppSum(f[k], "hBlindDxAll_V", Form("px%d", k)); ppFit(hx, 0);
        double dx = ppMu[0], sx = ppSignif[0];
        TH1F* hy = ppSum(f[k], "hBlindDyAll_H", Form("py%d", k)); ppFit(hy, 0);
        printf("  %-24s %10.3f %10.3f %10.1f %10.0f   (dx %.0f#sigma)\n",
               lab[k], dx, ppMu[0], ppSignif[0], ppEnt[0], sx);
    }
}
