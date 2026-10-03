// plotFcsResid.C -- FCS residual, FTT hits off vs aligned FTT hits on.
// Two plots: the four ECAL panels overlaid per arm, and the core sigma summary.
// Mixed event normalised over |d| 80-200 cm and subtracted (picoMatch books the
// "Same-Mixed" slots but never fills them).
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

double pcCore, pcCoreE, pcTail, pcFrac;
TF1* pcFn = 0;

TH1D* pcSub(TFile* f, const char* panel, const char* xy, const char* tag) {
    TH1* hs = (TH1*)f->Get(Form("%s%sSame Event_Primary", panel, xy));
    TH1* hm = (TH1*)f->Get(Form("%s%sMixed Event_Primary", panel, xy));
    if (!hs || !hm) return 0;
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= hs->GetNbinsX(); ib++) {
        double x = fabs(hs->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 200) continue;
        aS += hs->GetBinContent(ib); aM += hm->GetBinContent(ib);
    }
    TH1D* d = (TH1D*)hs->Clone(Form("psub_%s%s_%s", panel, xy, tag));
    d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("pscl_%s%s_%s", panel, xy, tag));
    m->SetDirectory(0);
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    return d;
}

void pcFit(TH1* h, const char* tag) {
    pcCore = 0; pcCoreE = 0; pcTail = 0; pcFrac = 0; pcFn = 0;
    if (!h) return;
    int pb = h->GetMaximumBin();
    double pk = h->GetXaxis()->GetBinCenter(pb);
    TF1* fn = new TF1(Form("pcf_%s", tag), "gaus(0)+gaus(3)+pol1(6)", pk - 25, pk + 25);
    fn->SetParameters(h->GetBinContent(pb) * 0.7, pk, 1.5,
                      h->GetBinContent(pb) * 0.3, pk, 8.0, 0, 0);
    fn->SetParLimits(1, pk - 5, pk + 5);
    fn->SetParLimits(2, 0.3, 4.0);
    fn->SetParLimits(4, pk - 8, pk + 8);
    fn->SetParLimits(5, 4.0, 40.0);
    h->Fit(fn, "QNR");
    pcCore = fabs(fn->GetParameter(2)); pcCoreE = fn->GetParError(2);
    pcTail = fabs(fn->GetParameter(5));
    double aC = fn->GetParameter(0) * pcCore, aT = fn->GetParameter(3) * pcTail;
    pcFrac = (aC + aT > 0) ? aC / (aC + aT) : 0;
    pcFn = fn;
}

void plotFcsResid(const char* fa, const char* fb, const char* outdir) {
    TFile* a = TFile::Open(fa);
    TFile* b = TFile::Open(fb);
    if (!a || a->IsZombie() || !b || b->IsZombie()) { printf("cannot open inputs\n"); return; }
    const char* pn[4] = {"EcalNorth", "EcalSouth", "EcalTop", "EcalBottom"};
    const char* vx[4] = {"dx", "dx", "dy", "dy"};
    const char* lb[4] = {"ECAL North, dx", "ECAL South, dx", "ECAL Top, dy", "ECAL Bottom, dy"};
    double cA[4], cAe[4], cB[4], cBe[4];

    gStyle->SetOptStat(0);
    TGaxis::SetMaxDigits(3);
    TCanvas* c1 = new TCanvas("cfr", "cfr", 1000, 760);
    c1->Divide(2, 2);
    for (int ip = 0; ip < 4; ip++) {
        c1->cd(ip + 1);
        gPad->SetGridx(); gPad->SetGridy();
        gPad->SetLeftMargin(0.13);
        TH1D* ha = pcSub(a, pn[ip], vx[ip], "A");
        TH1D* hb = pcSub(b, pn[ip], vx[ip], "B");
        if (!ha || !hb) continue;
        pcFit(ha, Form("a%d", ip)); cA[ip] = pcCore; cAe[ip] = pcCoreE;
        TF1* fA = pcFn;
        pcFit(hb, Form("b%d", ip)); cB[ip] = pcCore; cBe[ip] = pcCoreE;
        TF1* fB = pcFn;
        ha->SetTitle(Form("%s;%s = FCS - track [cm];counts / 0.5 cm", lb[ip], vx[ip]));
        ha->GetYaxis()->SetTitleOffset(1.45);
        ha->GetXaxis()->SetRangeUser(-30, 30);
        ha->SetLineColor(kAzure + 2); ha->SetLineWidth(2);
        hb->SetLineColor(kRed + 1);   hb->SetLineWidth(2);
        ha->Draw("hist");
        hb->Draw("hist same");
        if (fA) { fA->SetLineColor(kAzure + 2); fA->SetLineStyle(2); fA->Draw("same"); }
        if (fB) { fB->SetLineColor(kRed + 1);   fB->SetLineStyle(2); fB->Draw("same"); }
        TLegend* lg = new TLegend(0.58, 0.72, 0.98, 0.90);
        lg->SetFillColor(0); lg->SetBorderSize(0); lg->SetTextSize(0.034);
        lg->AddEntry(ha, Form("FTT off: core %.2f cm", cA[ip]), "l");
        lg->AddEntry(hb, Form("FTT on:  core %.2f cm", cB[ip]), "l");
        lg->Draw();
    }
    c1->SaveAs(Form("%s/fcsresid_panels.png", outdir));

    TCanvas* c2 = new TCanvas("cfs", "cfs", 680, 500);
    gPad->SetGridy();
    double xx[4] = {0.8, 1.8, 2.8, 3.8}, xx2[4] = {1.2, 2.2, 3.2, 4.2}, ex[4] = {0, 0, 0, 0};
    TGraphErrors* gA = new TGraphErrors(4, xx, cA, ex, cAe);
    TGraphErrors* gB = new TGraphErrors(4, xx2, cB, ex, cBe);
    TH2F* fr = new TH2F("fr", "FCS residual core, Primary tracks;;core #sigma [cm]",
                        4, 0.3, 4.7, 10, 0, 4.0);
    for (int jp = 0; jp < 4; jp++) fr->GetXaxis()->SetBinLabel(jp + 1, lb[jp]);
    fr->GetXaxis()->SetLabelSize(0.038);
    fr->Draw();
    gA->SetMarkerStyle(20); gA->SetMarkerSize(1.5); gA->SetMarkerColor(kAzure + 2);
    gA->SetLineColor(kAzure + 2); gA->SetLineWidth(2); gA->Draw("P");
    gB->SetMarkerStyle(21); gB->SetMarkerSize(1.5); gB->SetMarkerColor(kRed + 1);
    gB->SetLineColor(kRed + 1); gB->SetLineWidth(2); gB->Draw("P");
    TLegend* l2 = new TLegend(0.14, 0.75, 0.62, 0.89);
    l2->SetFillColor(0); l2->SetBorderSize(0); l2->SetTextSize(0.036);
    l2->AddEntry(gA, "FTT hits not on the track", "p");
    l2->AddEntry(gB, "aligned FTT hits on the track", "p");
    l2->Draw();
    c2->SaveAs(Form("%s/fcsresid_core.png", outdir));

    printf("\n  panel            A core      B core     B/A\n");
    for (int kp = 0; kp < 4; kp++)
        printf("  %-16s %7.3f     %7.3f   %6.3f\n", lb[kp], cA[kp], cB[kp],
               cA[kp] > 0 ? cB[kp] / cA[kp] : 0);
}
