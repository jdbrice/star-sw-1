// plotWhatWeGained.C -- July -> October on the SAME 74 forward-stream files, both
// panels the same comparison.  The point of the figure is that the headline gain is in
// HOW MANY tracks match, which a width-only plot hides completely.
//   left  : matched clusters per event -- the yield
//   right : residual width ratio        -- the cost, such as it is
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

void plotWhatWeGained(const char* out = "FstFttFlipTest/plots/what_we_gained.png") {
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("cwg", "cwg", 1320, 560);
    c->Divide(2, 1);

    // ---- left: match rate per event, BLCVtx, |r*dPhi| < 30 cm ----
    c->cd(1);
    gPad->SetGridy(); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.11);
    double jul[2] = {0.97293, 1.03253};
    double oct[2] = {1.39733, 1.41206};
    TH2F* f1 = new TH2F("wgf1", "How many tracks match;;matched clusters per event",
                        2, 0.4, 2.6, 10, 0.0, 1.85);
    f1->GetXaxis()->SetBinLabel(1, "ECAL North");
    f1->GetXaxis()->SetBinLabel(2, "ECAL South");
    f1->GetXaxis()->SetLabelSize(0.065);
    f1->GetYaxis()->SetLabelSize(0.048); f1->GetYaxis()->SetTitleSize(0.052);
    f1->GetYaxis()->SetTitleOffset(1.42);
    f1->SetTitleSize(0.065, "t");
    f1->Draw();
    double xj[2] = {0.85, 1.85}, xo[2] = {1.15, 2.15}, ze[2] = {0, 0};
    TGraph* gj = new TGraph(2, xj, jul); TGraph* go = new TGraph(2, xo, oct);
    gj->SetMarkerStyle(24); gj->SetMarkerSize(3.0); gj->SetMarkerColor(kAzure+2);
    go->SetMarkerStyle(20); go->SetMarkerSize(3.0); go->SetMarkerColor(kRed+1);
    gj->Draw("P"); go->Draw("P");
    for (int i = 0; i < 2; i++) {
        TArrow* a = new TArrow(xj[i] + 0.06, jul[i] + 0.04, xo[i] - 0.06, oct[i] - 0.04, 0.025, "|>");
        a->SetLineWidth(3); a->SetLineColor(kGray+2); a->SetFillColor(kGray+2); a->Draw();
        TLatex* L = new TLatex(xo[i] + 0.06, oct[i], Form("#times%.2f", oct[i]/jul[i]));
        L->SetTextSize(0.062); L->SetTextColor(kRed+1); L->Draw();
    }
    TLegend* l1 = new TLegend(0.21, 0.17, 0.63, 0.34);
    l1->SetFillColor(0); l1->SetBorderSize(0); l1->SetTextSize(0.050);
    l1->AddEntry(gj, "July (old code)", "p");
    l1->AddEntry(go, "October (now)", "p");
    l1->Draw();
    TLatex t1; t1.SetNDC(); t1.SetTextSize(0.050); t1.SetTextColor(kRed+1);
    t1.DrawLatex(0.21, 0.40, "#bf{~40% more matches}");

    // ---- right: residual width, same comparison ----
    c->cd(2);
    gPad->SetGridy(); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.11);
    const int kNP = 4;
    double net[kNP] = {0.958, 0.933, 1.342, 1.033};
    const char* pn[kNP] = {"North dx", "South dx", "Top dy", "Bottom dy"};
    TH2F* f2 = new TH2F("wgf2", "How wide the match is;;width ratio, Oct / July",
                        4, 0.4, 4.6, 10, 0.80, 1.45);
    for (int ip = 0; ip < kNP; ip++) f2->GetXaxis()->SetBinLabel(ip + 1, pn[ip]);
    f2->GetXaxis()->SetLabelSize(0.058);
    f2->GetYaxis()->SetLabelSize(0.048); f2->GetYaxis()->SetTitleSize(0.052);
    f2->GetYaxis()->SetTitleOffset(1.42);
    f2->SetTitleSize(0.065, "t");
    f2->Draw();
    TLine* one = new TLine(0.4, 1.0, 4.6, 1.0);
    one->SetLineWidth(3); one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();
    double xw[kNP] = {1, 2, 3, 4};
    TGraph* gw = new TGraph(kNP, xw, net);
    gw->SetMarkerStyle(21); gw->SetMarkerSize(2.8); gw->SetMarkerColor(kAzure+2);
    gw->Draw("P");
    TLatex t2; t2.SetNDC(); t2.SetTextSize(0.044); t2.SetTextColor(kGray+3);
    t2.DrawLatex(0.19, 0.84, "below 1 = narrower");
    t2.SetTextSize(0.042); t2.SetTextColor(kRed+1);
    t2.DrawLatex(0.19, 0.20, "dy still under study");
    c->SaveAs(out);
    printf("  wrote %s\n", out);
}
