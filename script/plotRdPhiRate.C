// plotRdPhiRate.C -- r*dPhi, ECAL, mixed-subtracted, normalised PER EVENT rather than
// per area, so the yield is kept instead of thrown away.
//
// Both forward-stream inputs are CURRENT picoMatch output (0.5 cm bins) -- the July one
// is the old production's picos re-analysed, so the tracking code differs but the
// analysis does not. Both are rebinned by 8 onto a 4 cm axis.
// The two forward-stream curves are the SAME 74 input files; July ran them uncapped
// (674938 events) and October at 5000/file (348984), so October is essentially the first
// part of each file and a per-event rate is the honest comparison between them.
// ZF 22 is different data entirely -- its rate is drawn for scale, not for comparison.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

TH1D* rrSub(const char* file, const char* type, const char* half, int rebin,
            double nev, const char* tag) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) return 0;
    TH1* hs = (TH1*)f->Get(Form("Ecal%sdpSame Event_%s", half, type));
    TH1* hm = (TH1*)f->Get(Form("Ecal%sdpMixed Event_%s", half, type));
    if (!hs || !hm) return 0;
    TH1D* d = (TH1D*)hs->Clone(Form("rr_%s_%s", half, tag)); d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("rrm_%s_%s", half, tag)); m->SetDirectory(0);
    if (rebin > 1) { d->Rebin(rebin); m->Rebin(rebin); }
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= d->GetNbinsX(); ib++) {
        double x = fabs(d->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 160) continue;
        aS += d->GetBinContent(ib); aM += m->GetBinContent(ib);
    }
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    if (nev > 0) d->Scale(1.0 / nev);
    return d;
}

void plotRdPhiRate(const char* fZf, double nZf, const char* fJul, double nJul,
                   const char* fOct, double nOct, const char* out) {
    const char* half[2] = {"North", "South"};
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("rrdp", "rrdp", 1000, 440);
    c->Divide(2, 1);
    printf("\n  ECAL r*dPhi, mixed-subtracted, PER EVENT (same 74 files for the two fwd arms)\n");
    printf("  %-7s %13s %13s %9s %13s\n", "half", "July /evt", "Oct /evt", "Oct/July", "ZF22 /evt");
    for (int ih = 0; ih < 2; ih++) {
        c->cd(ih + 1); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.15);
        TH1D* a = rrSub(fZf,  "Primary", half[ih], 8, nZf,  "zf");
        TH1D* b = rrSub(fJul, "BLCVtx",  half[ih], 8, nJul, "jul");
        TH1D* d = rrSub(fOct, "BLCVtx",  half[ih], 8, nOct, "oct");
        if (!a || !b || !d) continue;
        b->SetTitle(Form("ECAL %s, r*dPhi;fcsR*dPhi [cm];matches per event / 4 cm", half[ih]));
        b->GetXaxis()->SetRangeUser(-60, 60);
        double mx = TMath::Max(a->GetMaximum(), TMath::Max(b->GetMaximum(), d->GetMaximum()));
        b->SetMaximum(mx * 1.35); b->SetMinimum(0);
        b->SetLineColor(kAzure+2); b->SetLineWidth(2);
        d->SetLineColor(kRed+1);   d->SetLineWidth(2);
        a->SetLineColor(kGray+2);  a->SetLineWidth(2); a->SetLineStyle(2);
        b->Draw("hist"); d->Draw("hist same"); a->Draw("hist same");
        TLegend* lg = new TLegend(0.55, 0.68, 0.99, 0.90);
        lg->SetFillColor(0); lg->SetBorderSize(0); lg->SetTextSize(0.040);
        lg->AddEntry(b, "fwd stream July", "l");
        lg->AddEntry(d, "fwd stream Oct", "l");
        lg->AddEntry(a, "ZF22 Primary (other data)", "l");
        lg->Draw();
        double ij = 0, io = 0, iz = 0;
        for (int jb = 1; jb <= b->GetNbinsX(); jb++) {
            if (fabs(b->GetXaxis()->GetBinCenter(jb)) <= 30) ij += b->GetBinContent(jb);
        }
        for (int kb = 1; kb <= d->GetNbinsX(); kb++) {
            if (fabs(d->GetXaxis()->GetBinCenter(kb)) <= 30) io += d->GetBinContent(kb);
        }
        for (int lb = 1; lb <= a->GetNbinsX(); lb++) {
            if (fabs(a->GetXaxis()->GetBinCenter(lb)) <= 30) iz += a->GetBinContent(lb);
        }
        printf("  %-7s %13.5f %13.5f %9.3f %13.5f\n", half[ih], ij, io, ij > 0 ? io / ij : 0, iz);
    }
    c->SaveAs(out);
}
