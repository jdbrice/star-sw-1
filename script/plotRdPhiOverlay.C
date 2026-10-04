// plotRdPhiOverlay.C -- r*dPhi only, ECAL only, mixed-event subtracted only.
// North in the left panel, South in the right, three datasets overlaid, so the
// north/south shift can be read directly.
// Everything is put on the July file's 4 cm axis and area normalised, since the three
// datasets differ in statistics by orders of magnitude.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

TH1D* roSub(const char* file, const char* type, const char* half, int rebin, const char* tag) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) return 0;
    TH1* hs = (TH1*)f->Get(Form("Ecal%sdpSame Event_%s", half, type));
    TH1* hm = (TH1*)f->Get(Form("Ecal%sdpMixed Event_%s", half, type));
    if (!hs || !hm) { printf("  missing Ecal%sdp_%s\n", half, type); return 0; }
    TH1D* d = (TH1D*)hs->Clone(Form("ro_%s_%s", half, tag)); d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("rom_%s_%s", half, tag)); m->SetDirectory(0);
    if (rebin > 1) { d->Rebin(rebin); m->Rebin(rebin); }
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= d->GetNbinsX(); ib++) {
        double x = fabs(d->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 160) continue;
        aS += d->GetBinContent(ib); aM += m->GetBinContent(ib);
    }
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    double it = d->Integral();
    if (it > 0) d->Scale(1.0 / it);
    return d;
}

void plotRdPhiOverlay(const char* fZf, const char* fJul, const char* fOct, const char* out) {
    const char* half[2] = {"North", "South"};
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("rodp", "rodp", 1000, 440);
    c->Divide(2, 1);
    printf("\n  ECAL r*dPhi, mixed-subtracted, area normalised\n");
    printf("  %-8s %22s %22s %22s\n", "half", "ZF22 Primary", "fwd July BLCVtx", "fwd Oct BLCVtx");
    printf("  %-8s %10s %11s %10s %11s %10s %11s\n", "", "peak", "mean", "peak", "mean", "peak", "mean");
    for (int ih = 0; ih < 2; ih++) {
        c->cd(ih + 1); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.13);
        TH1D* a = roSub(fZf,  "Primary", half[ih], 8, "zf");
        TH1D* b = roSub(fJul, "BLCVtx",  half[ih], 1, "jul");
        TH1D* d = roSub(fOct, "BLCVtx",  half[ih], 8, "oct");
        if (!a || !b || !d) continue;
        a->SetTitle(Form("ECAL %s, r*dPhi;fcsR*dPhi [cm];area normalised", half[ih]));
        a->GetXaxis()->SetRangeUser(-60, 60);
        double mx = TMath::Max(a->GetMaximum(), TMath::Max(b->GetMaximum(), d->GetMaximum()));
        a->SetMaximum(mx * 1.3); a->SetMinimum(0);
        a->SetLineColor(kBlack);   a->SetLineWidth(2);
        b->SetLineColor(kAzure+2); b->SetLineWidth(2);
        d->SetLineColor(kRed+1);   d->SetLineWidth(2);
        a->Draw("hist"); b->Draw("hist same"); d->Draw("hist same");
        TLegend* lg = new TLegend(0.58, 0.70, 0.99, 0.90);
        lg->SetFillColor(0); lg->SetBorderSize(0); lg->SetTextSize(0.042);
        lg->AddEntry(a, "ZF22, Primary", "l");
        lg->AddEntry(b, "fwd stream July", "l");
        lg->AddEntry(d, "fwd stream Oct", "l");
        lg->Draw();
        TH1D* qa=(TH1D*)a->Clone("qa"); qa->SetDirectory(0); qa->GetXaxis()->SetRangeUser(-30,30);
        TH1D* qb=(TH1D*)b->Clone("qb"); qb->SetDirectory(0); qb->GetXaxis()->SetRangeUser(-30,30);
        TH1D* qd=(TH1D*)d->Clone("qd"); qd->SetDirectory(0); qd->GetXaxis()->SetRangeUser(-30,30);
        printf("  %-8s %10.2f %11.3f %10.2f %11.3f %10.2f %11.3f\n", half[ih],
               a->GetXaxis()->GetBinCenter(a->GetMaximumBin()), qa->GetMean(),
               b->GetXaxis()->GetBinCenter(b->GetMaximumBin()), qb->GetMean(),
               d->GetXaxis()->GetBinCenter(d->GetMaximumBin()), qd->GetMean());
    }
    c->SaveAs(out);
}
