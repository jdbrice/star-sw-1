// plotTrackMatchExtra.C -- the three remaining panels of the 2026-07-06 layout:
//   fcsTrkxy.png       XY of FCS cluster and of the track projection
//   fcsTrkEt.png       FCS and track E and ET
//   fcsTrkQuality.png  track q/pT, eta, phi, DCAz
// Same 2x2 arrangement and the histograms' own ranges, so they match the old page.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

void plotTrackMatchExtra(const char* file, const char* type, const char* outdir,
                         const char* label) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }
    gStyle->SetOptStat(1111);
    printf("\n=== %s [%s] extra panels ===\n", label, type);

    // --- xy ---
    TCanvas* c1 = new TCanvas(Form("xy_%s", type), "xy", 1000, 760);
    c1->Divide(2, 2);
    const char* xyn[4] = {"xyEcal", "xyHcal", "", ""};
    for (int i = 0; i < 4; i++) {
        c1->cd(i + 1);
        TH2* hxy = 0;
        if (i < 2) hxy = (TH2*)f->Get(xyn[i]);
        else       hxy = (TH2*)f->Get(Form(i == 2 ? "xyETrk%s" : "xyHTrk%s", type));
        if (hxy) hxy->Draw("colz");
    }
    c1->SaveAs(Form("%s/fcsTrkxy.png", outdir));

    // --- E and ET ---
    TCanvas* c2 = new TCanvas(Form("et_%s", type), "et", 1000, 760);
    c2->Divide(2, 2);
    const char* en[2] = {"EEcal", "EHcal"};
    const char* tn[2] = {"ETEcal", "ETHcal"};
    for (int j = 0; j < 2; j++) {
        c2->cd(j + 1);   gPad->SetLogy();
        TH1* hee = (TH1*)f->Get(en[j]); if (hee) hee->Draw("hist");
        c2->cd(j + 3);   gPad->SetLogy();
        TH1* het = (TH1*)f->Get(tn[j]); if (het) het->Draw("hist");
    }
    c2->SaveAs(Form("%s/fcsTrkEt.png", outdir));

    // --- track quality ---
    TCanvas* c3 = new TCanvas(Form("q_%s", type), "q", 1000, 760);
    c3->Divide(2, 2);
    const char* qn[4] = {"trkQPt%s", "trkEta%s", "trkPhi%s", "trkDcaZ%s"};
    for (int k = 0; k < 4; k++) {
        c3->cd(k + 1);
        TH1* hq = (TH1*)f->Get(Form(qn[k], type));
        if (hq) { hq->SetLineWidth(2); hq->Draw("hist"); printf("  %-16s entries %.4g\n",
                 hq->GetName(), hq->GetEntries()); }
        else printf("  %s : missing\n", Form(qn[k], type));
    }
    c3->SaveAs(Form("%s/fcsTrkQuality.png", outdir));
    gStyle->SetOptStat(0);
}
