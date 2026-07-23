// plotChargeSharing.C -- render fttChargeSharing.root histograms to PNGs
// for the ftt_sim_hit_maker status page (charge-sharing subpage).
void plotChargeSharing(const char* file = "fttChargeSharing.root",
                        const char* outdir = "/direct/star+u/akio/fcstrk11/star-sw-fwd/ftt_sim_hit_maker/plots_cs") {
    gStyle->SetOptStat(1110);
    gStyle->SetPadLeftMargin(0.12);

    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }

    TH1F* hN = (TH1F*)f->Get("hNStrips");
    TH1F* hSum = (TH1F*)f->Get("hSumAdc");
    TH1F* hSumN1 = (TH1F*)f->Get("hSumAdcN1");
    TH1F* hSumN2p = (TH1F*)f->Get("hSumAdcN2p");
    TH1F* hPeak = (TH1F*)f->Get("hPeakAdc");
    TH1F* hRawInTime = (TH1F*)f->Get("hRawAdcInTime");
    TProfile* pAll = (TProfile*)f->Get("pProfileAll");

    TCanvas c("c", "c", 500, 400);

    c.cd(); hN->GetXaxis()->SetRangeUser(0, 10); hN->Draw();
    c.SetLogy(1); c.SaveAs(Form("%s/hNStrips.png", outdir)); c.SetLogy(0);

    c.cd(); hSum->GetXaxis()->SetRangeUser(0, 2000); hSum->Draw();
    c.SetLogy(1); c.SaveAs(Form("%s/hSumAdc.png", outdir)); c.SetLogy(0);

    c.cd();
    hSumN2p->GetXaxis()->SetRangeUser(0, 2000);
    hSumN2p->SetLineColor(kBlue); hSumN2p->SetTitle("cluster total ADC: nStrips==1 (red) vs nStrips>=2 (blue)");
    hSumN2p->Draw();
    hSumN1->SetLineColor(kRed); hSumN1->Draw("same");
    c.SetLogy(1); c.SaveAs(Form("%s/hSumAdcByMult.png", outdir)); c.SetLogy(0);

    c.cd(); hPeak->Draw();
    c.SaveAs(Form("%s/hPeakAdc.png", outdir));

    c.cd(); hRawInTime->Draw();
    c.SetLogy(1); c.SaveAs(Form("%s/hRawAdcInTime.png", outdir)); c.SetLogy(0);

    c.cd(); pAll->SetMinimum(0); pAll->SetMarkerStyle(20); pAll->Draw();
    c.SaveAs(Form("%s/pProfileAll.png", outdir));

    const char* multTag[6] = { "", "", "2", "3", "4", "5p" };
    for (int m = 2; m <= 5; m++) {
        TProfile* p = (TProfile*)f->Get(Form("pProfile_nStrips%s", multTag[m]));
        if (!p) continue;
        c.cd(); p->SetMinimum(0); p->SetMarkerStyle(20); p->Draw();
        c.SaveAs(Form("%s/pProfile_nStrips%s.png", outdir, multTag[m]));
    }

    printf("wrote plots to %s\n", outdir);
}
