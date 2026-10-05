// plotReconcileFactors.C -- one figure for a general audience: what each change bought.
//
// From the fwd_stream reconciliation (tracking.html section 13): three arms of 74 jobs
// on the SAME Run22 day081 data, each changing one thing, plus the 2026-07-06 result.
// Each factor is a ratio of FCS residual RMS, so 1.0 = no change and BELOW 1 = better.
// The three multiply to the net, which is why they can be shown side by side.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

void plotReconcileFactors(const char* out = "FstFttFlipTest/plots/reconcile_factors.png") {
    const int kNP = 4;
    const char* pn[kNP] = {"ECAL North dx", "ECAL South dx", "ECAL Top dy", "ECAL Bottom dy"};
    // measured, see tracking.html section 13
    double code [kNP] = {1.063, 1.056, 1.372, 1.153};
    double field[kNP] = {1.072, 1.082, 1.157, 1.114};
    double geo  [kNP] = {0.841, 0.817, 0.846, 0.804};
    int col[kNP] = {kAzure+2, kAzure-3, kRed+1, kOrange+7};
    int mk [kNP] = {20, 21, 22, 33};

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("crf", "crf", 820, 560);
    gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
    TH2F* fr = new TH2F("rffr",
        "What each change did to the FCS#minustrack residual;;residual ratio  (below 1 = better)",
        3, 0.4, 3.6, 10, 0.70, 1.62);
    fr->GetXaxis()->SetBinLabel(1, "tracking code");
    fr->GetXaxis()->SetBinLabel(2, "real B field");
    fr->GetXaxis()->SetBinLabel(3, "geometry switches");
    fr->GetXaxis()->SetLabelSize(0.050);
    fr->GetYaxis()->SetTitleSize(0.040);
    fr->GetYaxis()->SetTitleOffset(1.45);
    fr->Draw();

    TLine* one = new TLine(0.4, 1.0, 3.6, 1.0);
    one->SetLineWidth(2); one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();

    TLegend* lg = new TLegend(0.60, 0.68, 0.98, 0.89);
    lg->SetFillColor(0); lg->SetBorderSize(0); lg->SetTextSize(0.033);
    for (int ip = 0; ip < kNP; ip++) {
        double x[3] = {1.0 + (ip - 1.5) * 0.13, 2.0 + (ip - 1.5) * 0.13, 3.0 + (ip - 1.5) * 0.13};
        double y[3] = {code[ip], field[ip], geo[ip]};
        TGraph* g = new TGraph(3, x, y);
        g->SetMarkerStyle(mk[ip]); g->SetMarkerSize(ip == 3 ? 2.0 : 1.6);
        g->SetMarkerColor(col[ip]); g->SetLineColor(col[ip]);
        g->Draw("P");
        lg->AddEntry(g, pn[ip], "p");
    }
    lg->Draw();

    TLatex t; t.SetNDC(); t.SetTextSize(0.030);
    t.DrawLatex(0.155, 0.865, "same data, same files, one change at a time");
    t.SetTextSize(0.027);
    t.DrawLatex(0.15, 0.20, "#bf{geometry switches}: the blind A/B settings, 15#minus20% better");
    t.DrawLatex(0.15, 0.165, "#bf{real B field}: costs 7#minus16%, and is the price of being correct");
    t.DrawLatex(0.15, 0.13, "#bf{tracking code}: 6#minus37%, concentrated in dy #minus still under study");
    c->SaveAs(out);
    printf("  wrote %s\n", out);
}
