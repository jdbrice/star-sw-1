void dumpSumAdcBins2(const char* file = "fttChargeSharing.root") {
    TFile f(file);
    TH1F* h1 = (TH1F*)f.Get("hSumAdcN1");
    TH1F* h2 = (TH1F*)f.Get("hSumAdcN2p");
    printf("bin  center     n1      n2p\n");
    for (int b=1;b<=110;b++){
        double c = h1->GetXaxis()->GetBinCenter(b);
        if (c > 2000) break;
        printf("%3d  %6.1f  %7.0f  %7.0f\n", b, c, h1->GetBinContent(b), h2->GetBinContent(b));
    }
}
