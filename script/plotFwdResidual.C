// plotFwdResidual.C
// Plot FST/FTT per-plane tracking-detector alignment residuals produced by
// StFwdResidualMaker (see StRoot/StFwdResidualMaker, macro/mudst/fwdDetResidual.C).
//
// For each FST disk (0-2):
//   fstResD{d}.png            r*dPhi residual, 1D
//   fstResD{d}_vsX/Y/R/Rphi.png   r*dPhi residual vs x,y,r,r*phi, 2D
//
// For each FTT plane (0-3), separately for H-strip (dy) and V-strip (dx):
//   fttResDyP{p}.png / fttResDxP{p}.png             1D
//   fttResDyP{p}_vsX/Y/R/Rphi.png / fttResDxP{p}_...  2D
//
// Usage: root4star -b -q 'plotFwdResidual.C("fwdDetResidual.root","residual")'

void plot1D(TFile* f, const char* hname, const char* png, const char* xtitle) {
    TH1F* h = (TH1F*)f->Get(hname);
    if (!h) { printf("WARN: missing %s\n", hname); return; }
    TCanvas* c = new TCanvas(Form("c_%s", png), "", 500, 400);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    h->SetLineColor(kBlue+1); h->SetLineWidth(2); h->SetFillColor(kBlue-9);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle("Hits");
    h->Draw("hist");
    gStyle->SetOptStat(1110);
    if (h->GetEntries() > 20) {
        h->Fit("gaus", "QN");
        TF1* fit = h->GetFunction("gaus");
        if (fit) { fit->SetLineColor(kRed); fit->SetLineWidth(2); fit->Draw("same"); }
    }
    c->Print(png);
    delete c;
}

void plot2D(TFile* f, const char* hname, const char* png, const char* xtitle, const char* ytitle) {
    TH2F* h = (TH2F*)f->Get(hname);
    if (!h) { printf("WARN: missing %s\n", hname); return; }
    TCanvas* c = new TCanvas(Form("c_%s", png), "", 500, 400);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.14);
    gStyle->SetOptStat(0);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->Draw("colz");

    // Overlay a linear fit to the per-x-bin mean (TProfile), so a tilted
    // residual band (misalignment) shows up as a nonzero slope.
    TProfile* prof = 0;
    if (h->GetEntries() > 50) {
        prof = h->ProfileX(Form("%s_pfx", h->GetName()));
        prof->SetDirectory(0);
        TFitResultPtr r = prof->Fit("pol1", "Q0S");
        if (r.Get() && r->IsValid()) {
            double slope = r->Parameter(1), slopeErr = r->ParError(1);
            TF1* fit = prof->GetFunction("pol1");
            if (fit) { fit->SetLineColor(kRed); fit->SetLineWidth(3); fit->Draw("same"); }
            prof->SetMarkerStyle(24); prof->SetMarkerSize(1.0); prof->SetMarkerColor(kBlack);
            prof->SetLineColor(kBlack); prof->SetLineWidth(1);
            prof->Draw("same p");
            double sig = (slopeErr > 0) ? fabs(slope)/slopeErr : 0;
            TLatex lat; lat.SetNDC(); lat.SetTextSize(0.042);
            lat.SetTextColor(sig > 3 ? kRed+1 : kBlack);
            lat.DrawLatex(0.15, 0.85, Form("slope = %.4f #pm %.4f", slope, slopeErr));
        }
    }
    c->Print(png);
    // Deleted only now: this TProfile (and its fitted TF1) are live primitives
    // on the pad until Print() rasterizes it -- freeing them earlier silently
    // wipes the overlay from the canvas before the image is written.
    if (prof) delete prof;
    delete c;
}

void plotPlaneUsage(TFile* f, const char* typeName, const char* png) {
    TH1F* h = (TH1F*)f->Get(Form("PlaneUsage/hPlaneUsage_%s", typeName));
    if (!h) { printf("WARN: missing PlaneUsage/hPlaneUsage_%s\n", typeName); return; }
    double all = h->GetBinContent(20);
    TH1F* hf = (TH1F*)h->Clone(Form("%s_frac", h->GetName()));
    hf->Sumw2();
    if (all > 0) hf->Scale(1.0/all);
    TCanvas* c = new TCanvas(Form("c_%s", png), "", 700, 400);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.16); gPad->SetGridy();
    hf->SetFillColor(kAzure-4); hf->SetLineColor(kAzure+2);
    hf->GetYaxis()->SetTitle(Form("fraction of tracks (N=%.0f)", all));
    hf->GetYaxis()->SetRangeUser(0, 1.05);
    hf->GetXaxis()->LabelsOption("h");
    hf->Draw("hist");
    c->Print(png);
    delete c;
}

void plotPlaneUsageAll(TFile* f, const char* png) {
    const char* typeName[6] = {"Global", "BLC", "Primary", "FwdVtx", "BLCVtx", "FCSConstrained"};
    int col[6] = {kBlack, kBlue+1, kGreen+2, kOrange+1, kMagenta+1, kRed+1};
    TCanvas* c = new TCanvas("c_planeUsageAll", "", 900, 500);
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.09); gPad->SetBottomMargin(0.14); gPad->SetGridy();
    gPad->SetRightMargin(0.20);
    TLegend* leg = new TLegend(0.81, 0.55, 0.995, 0.90);
    leg->SetBorderSize(0); leg->SetTextSize(0.035);
    bool first = true;
    for (int t = 0; t < 6; t++) {
        TH1F* h = (TH1F*)f->Get(Form("PlaneUsage/hPlaneUsage_%s", typeName[t]));
        if (!h) continue;
        double all = h->GetBinContent(20);
        TH1F* hf = (TH1F*)h->Clone(Form("frac_%s", typeName[t]));
        if (all > 0) hf->Scale(1.0/all);
        hf->SetLineColor(col[t]); hf->SetLineWidth(2);
        hf->GetYaxis()->SetTitle("fraction of tracks");
        hf->GetYaxis()->SetRangeUser(0, 1.05);
        hf->GetXaxis()->LabelsOption("h");
        hf->Draw(first ? "hist" : "hist same");
        first = false;
        leg->AddEntry(hf, Form("%s (N=%.0f)", typeName[t], all), "l");
    }
    leg->Draw();
    c->Print(png);
    delete c;
}

void plotFwdResidual(const char* fname = "fwdDetResidual.root", const char* outdir = "residual") {
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", fname); return; }
    printf("Reading %s -> %s/\n", fname, outdir);
    gSystem->mkdir(outdir, kTRUE);
    gErrorIgnoreLevel = kWarning;

    // ── FST ──────────────────────────────────────────────────────────────
    for (int d = 0; d < 3; d++) {
        plot1D(f, Form("FST/disk%d/h_fst_d%d_rdphi", d, d),
               Form("%s/fstResD%d.png", outdir, d), "r#upoint#Delta#phi [cm]");
        plot2D(f, Form("FST/disk%d/h2_fst_d%d_rdphi_vs_x", d, d),
               Form("%s/fstResD%d_vsX.png", outdir, d), "x [cm]", "r#upoint#Delta#phi [cm]");
        plot2D(f, Form("FST/disk%d/h2_fst_d%d_rdphi_vs_y", d, d),
               Form("%s/fstResD%d_vsY.png", outdir, d), "y [cm]", "r#upoint#Delta#phi [cm]");
        plot2D(f, Form("FST/disk%d/h2_fst_d%d_rdphi_vs_r", d, d),
               Form("%s/fstResD%d_vsR.png", outdir, d), "r [cm]", "r#upoint#Delta#phi [cm]");
        plot2D(f, Form("FST/disk%d/h2_fst_d%d_rdphi_vs_rphi", d, d),
               Form("%s/fstResD%d_vsRphi.png", outdir, d), "r#upoint#phi [cm]", "r#upoint#Delta#phi [cm]");
    }

    // ── FTT (H-strip -> dy, V-strip -> dx) ─────────────────────────────────
    for (int p = 0; p < 4; p++) {
        plot1D(f, Form("FTT/plane%d/h_ftt_p%d_dy", p, p),
               Form("%s/fttResDyP%d.png", outdir, p), "#Delta y [cm]");
        plot2D(f, Form("FTT/plane%d/h2_ftt_p%d_dy_vs_x", p, p),
               Form("%s/fttResDyP%d_vsX.png", outdir, p), "x [cm]", "#Delta y [cm]");
        plot2D(f, Form("FTT/plane%d/h2_ftt_p%d_dy_vs_y", p, p),
               Form("%s/fttResDyP%d_vsY.png", outdir, p), "y [cm]", "#Delta y [cm]");
        plot2D(f, Form("FTT/plane%d/h2_ftt_p%d_dy_vs_r", p, p),
               Form("%s/fttResDyP%d_vsR.png", outdir, p), "r [cm]", "#Delta y [cm]");
        plot2D(f, Form("FTT/plane%d/h2_ftt_p%d_dy_vs_rphi", p, p),
               Form("%s/fttResDyP%d_vsRphi.png", outdir, p), "r#upoint#phi [cm]", "#Delta y [cm]");

        plot1D(f, Form("FTT/plane%d/h_ftt_p%d_dx", p, p),
               Form("%s/fttResDxP%d.png", outdir, p), "#Delta x [cm]");
        plot2D(f, Form("FTT/plane%d/h2_ftt_p%d_dx_vs_x", p, p),
               Form("%s/fttResDxP%d_vsX.png", outdir, p), "x [cm]", "#Delta x [cm]");
        plot2D(f, Form("FTT/plane%d/h2_ftt_p%d_dx_vs_y", p, p),
               Form("%s/fttResDxP%d_vsY.png", outdir, p), "y [cm]", "#Delta x [cm]");
        plot2D(f, Form("FTT/plane%d/h2_ftt_p%d_dx_vs_r", p, p),
               Form("%s/fttResDxP%d_vsR.png", outdir, p), "r [cm]", "#Delta x [cm]");
        plot2D(f, Form("FTT/plane%d/h2_ftt_p%d_dx_vs_rphi", p, p),
               Form("%s/fttResDxP%d_vsRphi.png", outdir, p), "r#upoint#phi [cm]", "#Delta x [cm]");
    }

    // ── Plane usage (fraction of tracks with a hit on each plane) ──────────
    const char* typeName[6] = {"Global", "BLC", "Primary", "FwdVtx", "BLCVtx", "FCSConstrained"};
    for (int t = 0; t < 6; t++) {
        plotPlaneUsage(f, typeName[t], Form("%s/planeUsage_%s.png", outdir, typeName[t]));
    }
    plotPlaneUsageAll(f, Form("%s/planeUsage_All.png", outdir));

    printf("Done.\n");
}
