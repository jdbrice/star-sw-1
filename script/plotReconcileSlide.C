// plotReconcileSlide.C -- slide version of the reconciliation decomposition.
// Same numbers as plotReconcileFactors.C, but sized for a projector: larger canvas,
// big axis and label fonts, big markers, and the explanatory text cut to one line
// since the speaker says the rest out loud.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

void plotReconcileSlide(const char* out = "FstFttFlipTest/plots/reconcile_factors_slide.png") {
    const int kNP = 4;
    const char* pn[kNP] = {"ECAL North dx", "ECAL South dx", "ECAL Top dy", "ECAL Bottom dy"};
    double code [kNP] = {1.063, 1.056, 1.372, 1.153};
    double field[kNP] = {1.072, 1.082, 1.157, 1.114};
    double geo  [kNP] = {0.841, 0.817, 0.846, 0.804};
    int col[kNP] = {kAzure+2, kAzure-3, kRed+1, kOrange+7};
    int mk [kNP] = {20, 21, 22, 33};

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("crs", "crs", 1280, 820);
    gPad->SetGridy();
    gPad->SetLeftMargin(0.155); gPad->SetBottomMargin(0.135);
    gPad->SetTopMargin(0.105);  gPad->SetRightMargin(0.04);

    TH2F* fr = new TH2F("rsfr", ";;FCS residual ratio", 3, 0.4, 3.6, 10, 0.70, 1.58);
    fr->GetXaxis()->SetBinLabel(1, "tracking code");
    fr->GetXaxis()->SetBinLabel(2, "real B field");
    fr->GetXaxis()->SetBinLabel(3, "geometry switches");
    fr->GetXaxis()->SetLabelSize(0.062);
    fr->GetYaxis()->SetLabelSize(0.050);
    fr->GetYaxis()->SetTitleSize(0.058);
    fr->GetYaxis()->SetTitleOffset(1.22);
    fr->Draw();

    TLine* one = new TLine(0.4, 1.0, 3.6, 1.0);
    one->SetLineWidth(3); one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    TLegend* lg = new TLegend(0.60, 0.62, 0.96, 0.87);
    lg->SetFillColor(0); lg->SetBorderSize(0); lg->SetTextSize(0.042);
    for (int ip = 0; ip < kNP; ip++) {
        double x[3] = {1.0 + (ip - 1.5) * 0.14, 2.0 + (ip - 1.5) * 0.14, 3.0 + (ip - 1.5) * 0.14};
        double y[3] = {code[ip], field[ip], geo[ip]};
        TGraph* g = new TGraph(3, x, y);
        g->SetMarkerStyle(mk[ip]); g->SetMarkerSize(ip == 3 ? 3.4 : 2.8);
        g->SetMarkerColor(col[ip]); g->SetLineColor(col[ip]);
        g->Draw("P");
        lg->AddEntry(g, pn[ip], "p");
    }
    lg->Draw();

    TLatex t; t.SetNDC();
    t.SetTextSize(0.055); t.SetTextColor(kBlack);
    t.DrawLatex(0.155, 0.925, "What each change did, same data, one at a time");
    t.SetTextSize(0.040); t.SetTextColor(kGray+3);
    t.DrawLatex(0.175, 0.185, "below 1 = better");
    // arrow making "better" unambiguous without words
    TArrow* a = new TArrow(0.60, 1.03, 0.60, 0.86, 0.022, "|>");
    a->SetLineWidth(3); a->SetLineColor(kGray+2); a->SetFillColor(kGray+2); a->Draw();
    c->SaveAs(out);
    printf("  wrote %s\n", out);
}
