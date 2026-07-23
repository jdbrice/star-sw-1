// inspectChargeSharing.C -- summarize fttChargeSharing.root (produced by
// runFttChargeSharing.C / StFttChargeSharingMaker) so we can decide on a
// charge-sharing model (see status_ftt_sim_maker.txt item 10) before
// touching StFttSimHitMaker's flat 1:8:1 model.
void inspectChargeSharing(const char* file = "fttChargeSharing.root") {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }

    TH1F* hN = (TH1F*)f->Get("hNStrips");
    TH1F* hSum = (TH1F*)f->Get("hSumAdc");
    TH1F* hPeak = (TH1F*)f->Get("hPeakAdc");
    TProfile* pAll = (TProfile*)f->Get("pProfileAll");

    double total = hN->Integral(0, hN->GetNbinsX() + 1);
    printf("\n=== cluster multiplicity (nStrips) -- total clusters = %.0f ===\n", total);
    for (int n = 1; n <= 10; n++) {
        int bin = hN->GetXaxis()->FindBin(n);
        double c = hN->GetBinContent(bin);
        printf("  nStrips=%2d: %8.0f  (%.2f%%)\n", n, c, 100.0 * c / total);
    }

    printf("\n=== sumAdc / peakAdc ===\n");
    printf("  sumAdc:  mean=%.1f rms=%.1f\n", hSum->GetMean(), hSum->GetRMS());
    printf("  peakAdc: mean=%.1f rms=%.1f\n", hPeak->GetMean(), hPeak->GetRMS());

    printf("\n=== mean ADC fraction vs strip offset from peak (all clusters) ===\n");
    for (int b = 1; b <= pAll->GetNbinsX(); b++) {
        double x = pAll->GetXaxis()->GetBinCenter(b);
        double ne = pAll->GetBinEntries(b);
        if (ne <= 0) continue;
        printf("  offset=%+.0f: <frac>=%.4f  entries=%.0f\n", x, pAll->GetBinContent(b), ne);
    }

    const char* multTag[6] = { "", "", "2", "3", "4", "5p" };
    for (int m = 2; m <= 5; m++) {
        TProfile* p = (TProfile*)f->Get(Form("pProfile_nStrips%s", multTag[m]));
        if (!p) continue;
        printf("\n=== mean ADC fraction vs offset, nStrips=%s ===\n", multTag[m]);
        for (int b = 1; b <= p->GetNbinsX(); b++) {
            double x = p->GetXaxis()->GetBinCenter(b);
            double ne = p->GetBinEntries(b);
            if (ne <= 0) continue;
            printf("  offset=%+.0f: <frac>=%.4f  entries=%.0f\n", x, p->GetBinContent(b), ne);
        }
    }
}
