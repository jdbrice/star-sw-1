// cmpOldNew.C -- put an old picoMatch output and a new one on the same axes.
//
// Nothing is un-fixed to do this. The old files kept the same histogram names and the
// same per-track-type layout, so the only thing in the way is binning: old is 4 cm,
// new is 0.5 cm. The new one is rebinned DOWN to the old axis, both are mixed-event
// subtracted the same way, and both are area-normalised because the statistics differ
// by orders of magnitude.
//
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

TH1D* onSub(TFile* f, const char* panel, const char* xy, const char* type, const char* tag) {
    TH1* hs = (TH1*)f->Get(Form("%s%sSame Event_%s", panel, xy, type));
    TH1* hm = (TH1*)f->Get(Form("%s%sMixed Event_%s", panel, xy, type));
    if (!hs || !hm) { printf("   missing %s%s_%s in %s\n", panel, xy, type, tag); return 0; }
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= hs->GetNbinsX(); ib++) {
        double x = fabs(hs->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 200) continue;
        aS += hs->GetBinContent(ib); aM += hm->GetBinContent(ib);
    }
    TH1D* d = (TH1D*)hs->Clone(Form("on_%s%s_%s_%s", panel, xy, type, tag));
    d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("onm_%s%s_%s_%s", panel, xy, type, tag));
    m->SetDirectory(0);
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    return d;
}

void cmpOldNew(const char* fOld, const char* fNew, const char* type,
               const char* labOld, const char* labNew, const char* out) {
    TFile* a = TFile::Open(fOld);
    TFile* b = TFile::Open(fNew);
    if (!a || a->IsZombie() || !b || b->IsZombie()) { printf("cannot open inputs\n"); return; }
    const char* pn[4] = {"EcalNorth", "EcalSouth", "EcalTop", "EcalBottom"};
    const char* vx[4] = {"dx", "dx", "dy", "dy"};
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("con", "con", 1000, 760);
    c->Divide(2, 2);
    printf("\n  %-12s %10s %10s   %10s %10s   %s\n", "panel",
           "old RMS", "new RMS", "old peak", "new peak", "(common 4 cm bins, |d|<30)");
    for (int ip = 0; ip < 4; ip++) {
        c->cd(ip + 1); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.13);
        TH1D* ho = onSub(a, pn[ip], vx[ip], type, "o");
        TH1D* hn = onSub(b, pn[ip], vx[ip], type, "n");
        if (!ho || !hn) continue;
        int rb = (int)(ho->GetXaxis()->GetBinWidth(1) / hn->GetXaxis()->GetBinWidth(1) + 0.5);
        if (rb > 1) hn->Rebin(rb);                       // new 0.5 cm -> old 4 cm
        double io = ho->Integral(), in = hn->Integral();
        if (io > 0) ho->Scale(1.0 / io);
        if (in > 0) hn->Scale(1.0 / in);
        ho->GetXaxis()->SetRangeUser(-60, 60);
        ho->SetTitle(Form("%s %s, %s;%s = FCS - track [cm];area-normalised",
                          pn[ip], vx[ip], type, vx[ip]));
        ho->SetLineColor(kBlack); ho->SetLineWidth(2);
        hn->SetLineColor(kRed + 1); hn->SetLineWidth(2);
        double mx = TMath::Max(ho->GetMaximum(), hn->GetMaximum());
        ho->SetMaximum(mx * 1.25); ho->SetMinimum(0);
        ho->Draw("hist"); hn->Draw("hist same");
        TLegend* lg = new TLegend(0.55, 0.74, 0.98, 0.90);
        lg->SetFillColor(0); lg->SetBorderSize(0); lg->SetTextSize(0.033);
        lg->AddEntry(ho, labOld, "l"); lg->AddEntry(hn, labNew, "l");
        lg->Draw();
        // measure on clones: setting a range on the drawn histograms would clip the plot
        TH1D* qo = (TH1D*)ho->Clone(Form("qo%d", ip)); qo->SetDirectory(0);
        TH1D* qn = (TH1D*)hn->Clone(Form("qn%d", ip)); qn->SetDirectory(0);
        qo->GetXaxis()->SetRangeUser(-30, 30);
        qn->GetXaxis()->SetRangeUser(-30, 30);
        printf("  %-12s %10.3f %10.3f   %10.3f %10.3f\n", pn[ip],
               qo->GetRMS(), qn->GetRMS(),
               qo->GetXaxis()->GetBinCenter(qo->GetMaximumBin()),
               qn->GetXaxis()->GetBinCenter(qn->GetMaximumBin()));
    }
    c->SaveAs(out);
}
