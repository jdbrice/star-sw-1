// plotDyVsY.C -- the figure for "a scale, not an offset": dy against the cluster's own y,
// with the straight-line fit drawn. Flat would mean a translation; a slope is a scale.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).
void plotDyVsY(const char* f, const char* type, const char* lab, const char* out) {
    TFile* fp = TFile::Open(f);
    if (!fp || fp->IsZombie()) { printf("cannot open %s\n", f); return; }
    TH2* hs = (TH2*)fp->Get(Form("ydyEcalSame Event_%s", type));
    TH2* hm = (TH2*)fp->Get(Form("ydyEcalMixed Event_%s", type));
    if (!hs || !hm) { printf("missing ydyEcal_%s\n", type); return; }
    const int kN = 6;
    double ylo[kN] = {-90, -60, -35, 10, 35, 60};
    double yhi[kN] = {-60, -35, -10, 35, 60, 90};
    double xc[kN], yv[kN], xe[kN], ye[kN]; int n = 0;
    for (int i = 0; i < kN; i++) {
        int b0 = hs->GetXaxis()->FindBin(ylo[i] + 1e-6), b1 = hs->GetXaxis()->FindBin(yhi[i] - 1e-6);
        TH1D* ps = hs->ProjectionY(Form("pys%d", i), b0, b1); ps->SetDirectory(0);
        TH1D* pm = hm->ProjectionY(Form("pym%d", i), b0, b1); pm->SetDirectory(0);
        double aS = 0, aM = 0;
        for (int jb = 1; jb <= ps->GetNbinsX(); jb++) {
            double x = fabs(ps->GetXaxis()->GetBinCenter(jb));
            if (x < 60 || x > 95) continue;
            aS += ps->GetBinContent(jb); aM += pm->GetBinContent(jb);
        }
        pm->Scale(aM > 0 ? aS / aM : 1.0); ps->Add(pm, -1.0);
        if (ps->Integral() < 400) continue;
        int pb = ps->GetMaximumBin(); double pk = ps->GetXaxis()->GetBinCenter(pb);
        TF1* g = new TF1(Form("pyg%d", i), "gaus", pk - 8, pk + 8);
        ps->Fit(g, "QNR");
        xc[n] = 0.5 * (ylo[i] + yhi[i]); xe[n] = 0.5 * (yhi[i] - ylo[i]);
        yv[n] = g->GetParameter(1);      ye[n] = g->GetParError(1);
        n++;
    }
    if (n < 3) { printf("too few slices\n"); return; }
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("cdy", "cdy", 760, 540);
    gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    TGraphErrors* gr = new TGraphErrors(n, xc, yv, xe, ye);
    gr->SetTitle(Form("%s, ECAL, %s;FCS cluster y [cm];dy = FCS #minus track [cm]", lab, type));
    gr->SetMarkerStyle(20); gr->SetMarkerSize(1.8); gr->SetMarkerColor(kAzure+2);
    gr->SetLineColor(kAzure+2); gr->SetLineWidth(2);
    gr->GetYaxis()->SetTitleOffset(1.3);
    gr->Draw("AP");
    TF1* ln = new TF1("pyln", "pol1", -95, 95);
    gr->Fit(ln, "QN");
    ln->SetLineColor(kRed+1); ln->SetLineWidth(3); ln->SetLineStyle(2); ln->Draw("same");
    TLine* z = new TLine(-95, 0, 95, 0); z->SetLineColor(kGray+2); z->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.042);
    t.DrawLatex(0.17, 0.85, Form("slope = %+.4f cm per cm of y", ln->GetParameter(1)));
    t.DrawLatex(0.17, 0.79, Form("intercept = %+.2f cm", ln->GetParameter(0)));
    t.SetTextSize(0.036); t.SetTextColor(kRed+1);
    t.DrawLatex(0.17, 0.72, "a slope means a SCALE, not an offset");
    c->SaveAs(out);
    printf("  wrote %s  (slope %+.5f, intercept %+.2f)\n", out, ln->GetParameter(1), ln->GetParameter(0));
}
