void dumpSumAdcBins(const char* file = "fttChargeSharing.root") {
    TFile f(file);
    TH1F* h = (TH1F*)f.Get("hSumAdc");
    TH1F* h1 = (TH1F*)f.Get("hSumAdcN1");
    TH1F* h2 = (TH1F*)f.Get("hSumAdcN2p");
    TH1F* praw = (TH1F*)f.Get("hRawAdcInTime");
    printf("bin  center   all    n1     n2p\n");
    for (int b=1;b<=30;b++){
        printf("%3d  %6.1f  %6.0f %6.0f %6.0f\n", b, h->GetXaxis()->GetBinCenter(b), h->GetBinContent(b), h1->GetBinContent(b), h2->GetBinContent(b));
    }
    printf("\nraw in-time ADC, first 30 bins (bin width=%.2f)\n", praw->GetXaxis()->GetBinWidth(1));
    for (int b=1;b<=30;b++){
        printf("%3d  center=%6.1f  n=%8.0f\n", b, praw->GetXaxis()->GetBinCenter(b), praw->GetBinContent(b));
    }
}
