// plotReconcile.C -- the reconciliation ladder on one axis.
// Subtracted truncated RMS (|d|<30) per panel, every configuration side by side,
// all on the old 4 cm axis so binning is never a variable.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

double prRms(const char* file, const char* nm, int rebin) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) return -1;
    TH1* hs = (TH1*)f->Get(Form("%sSame Event_Primary", nm));
    TH1* hm = (TH1*)f->Get(Form("%sMixed Event_Primary", nm));
    if (!hs || !hm) { f->Close(); return -1; }
    TH1D* d = (TH1D*)hs->Clone("prd"); d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone("prm"); m->SetDirectory(0);
    if (rebin > 1) { d->Rebin(rebin); m->Rebin(rebin); }
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= d->GetNbinsX(); ib++) {
        double x = fabs(d->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 160) continue;
        aS += d->GetBinContent(ib); aM += m->GetBinContent(ib);
    }
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    d->GetXaxis()->SetRangeUser(-30, 30);
    double v = d->GetRMS(); f->Close(); return v;
}

void plotReconcile(const char* fOld, const char* fOldGeo, const char* fCur,
                   const char* out) {
    const char* nm[4] = {"EcalNorthdx", "EcalSouthdx", "EcalTopdy", "EcalBottomdy"};
    const char* lb[4] = {"North dx", "South dx", "Top dy", "Bottom dy"};
    double vOld[4], vGeo[4], vCur[4];
    for (int ip = 0; ip < 4; ip++) {
        vOld[ip] = prRms(fOld,    nm[ip], 1);
        vGeo[ip] = prRms(fOldGeo, nm[ip], 8);
        vCur[ip] = prRms(fCur,    nm[ip], 8);
    }
    gStyle->SetOptStat(0);
    TGaxis::SetMaxDigits(3);
    TCanvas* c = new TCanvas("cpr", "cpr", 760, 520);
    gPad->SetGridy(); gPad->SetLeftMargin(0.12);
    TH2F* fr = new TH2F("prfr", "fwd_stream, same data, Primary tracks;;FCS residual RMS [cm]",
                        4, 0.3, 4.7, 10, 0, 24);
    for (int jp = 0; jp < 4; jp++) fr->GetXaxis()->SetBinLabel(jp + 1, lb[jp]);
    fr->GetXaxis()->SetLabelSize(0.045);
    fr->Draw();
    double x1[4] = {0.75, 1.75, 2.75, 3.75};
    double x2[4] = {1.00, 2.00, 3.00, 4.00};
    double x3[4] = {1.25, 2.25, 3.25, 4.25};
    double ze[4] = {0, 0, 0, 0};
    TGraphErrors* g1 = new TGraphErrors(4, x1, vOld, ze, ze);
    TGraphErrors* g2 = new TGraphErrors(4, x2, vGeo, ze, ze);
    TGraphErrors* g3 = new TGraphErrors(4, x3, vCur, ze, ze);
    g1->SetMarkerStyle(20); g1->SetMarkerSize(1.6); g1->SetMarkerColor(kBlack);
    g2->SetMarkerStyle(21); g2->SetMarkerSize(1.6); g2->SetMarkerColor(kOrange + 7);
    g3->SetMarkerStyle(22); g3->SetMarkerSize(1.9); g3->SetMarkerColor(kAzure + 2);
    g1->Draw("P"); g2->Draw("P"); g3->Draw("P");
    TLegend* lg = new TLegend(0.13, 0.74, 0.74, 0.90);
    lg->SetFillColor(0); lg->SetBorderSize(0); lg->SetTextSize(0.034);
    lg->AddEntry(g1, "2026-07-06 production (old code, GEOFIX 1,1,1)", "p");
    lg->AddEntry(g2, "current code, GEOFIX 1,1,1", "p");
    lg->AddEntry(g3, "current code, GEOFIX 1,0,0", "p");
    lg->Draw();
    c->SaveAs(out);
    printf("\n  panel        old    cur1,1,1   cur1,0,0\n");
    for (int kp = 0; kp < 4; kp++)
        printf("  %-12s %7.3f %9.3f %9.3f\n", lb[kp], vOld[kp], vGeo[kp], vCur[kp]);
}
