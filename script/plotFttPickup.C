// plotFttPickup.C -- FTT pickup probability per station, from StFwdResidualMaker's
// hPlaneUsage.
//
// The control arm (fttNoAdd=1) reads EXACTLY 0.0000 in every FTT cell, verified on
// pico_iter2_20260930: with no sTGC hit ever put on a track there is nothing to use.
// So the numbers here are the pickup probability directly -- no baseline subtraction,
// and any non-zero cell in a NOADD=0 run is new information.
//
// Binning of hPlaneUsage (see StFwdResidualMaker::fillPlaneUsage):
//   bin 0        vertex used (all types except Global)
//   bins 1-3     FST disk 0,1,2
//   4 + 3*p + 0  FTT station p, V strip -> measures x
//   4 + 3*p + 1  FTT station p, H strip -> measures y
//   4 + 3*p + 2  FTT station p, diagonal u/v -- NOT used in the fit (the loader skips
//                equal-covariance points), so this should stay 0 and is a useful check
//   17,18        Ecal, Hcal
//   19           every track -> the denominator
// TH1 bins are 1-based, hence GetBinContent(n+1) below.
//
// Usage:
//   root -l -b -q 'plotFttPickup.C("fttadd_Primary.root")'
//   root -l -b -q 'plotFttPickup.C("fttadd_Primary.root","Primary","plots/pickup.png")'
//
// CINT: fixed-size arrays, unique loop variable names, no std::map (see CLAUDE.md).

void plotFttPickup(const char* fname = "fttadd_Primary.root",
                   const char* type  = "Primary",
                   const char* out   = "FstFttFlipTest/plots/ftt_pickup.png")
{
    TFile* fpk = TFile::Open(fname);
    if (!fpk || fpk->IsZombie()) { printf("cannot open %s\n", fname); return; }

    TH1F* hu = (TH1F*)fpk->Get(Form("PlaneUsage/hPlaneUsage_%s", type));
    if (!hu) { printf("missing PlaneUsage/hPlaneUsage_%s in %s\n", type, fname); return; }

    double ntr = hu->GetBinContent(20);          // bin 19 = every track
    if (ntr <= 0) { printf("no tracks of type %s\n", type); return; }

    printf("\n=== FTT pickup probability, %s tracks, %s ===\n", type, fname);
    printf("  tracks = %.0f\n", ntr);
    printf("  FST disks : %.4f  %.4f  %.4f\n",
           hu->GetBinContent(2)/ntr, hu->GetBinContent(3)/ntr, hu->GetBinContent(4)/ntr);

    double px[4], py[4], pd[4], ex[4], ey[4], sx[4];
    for (int ip = 0; ip < 4; ip++) {
        double cx = hu->GetBinContent(4 + 3*ip + 1);
        double cy = hu->GetBinContent(4 + 3*ip + 2);
        double cd = hu->GetBinContent(4 + 3*ip + 3);
        px[ip] = cx/ntr; py[ip] = cy/ntr; pd[ip] = cd/ntr;
        // binomial error on a probability
        ex[ip] = sqrt(TMath::Max(1e-12, px[ip]*(1-px[ip])/ntr));
        ey[ip] = sqrt(TMath::Max(1e-12, py[ip]*(1-py[ip])/ntr));
        sx[ip] = ip + 1;
        printf("  station %d : x(V) %.4f +- %.4f   y(H) %.4f +- %.4f   u/v %.4f%s\n",
               ip+1, px[ip], ex[ip], py[ip], ey[ip], pd[ip],
               (pd[ip] > 1e-6) ? "   <-- NONZERO: diagonal strips should be skipped" : "");
    }

    double anyx = 0, anyy = 0;
    for (int iq = 0; iq < 4; iq++) { anyx += px[iq]; anyy += py[iq]; }
    printf("  mean over stations: x %.4f   y %.4f\n", anyx/4, anyy/4);
    if (anyx + anyy < 1e-6)
        printf("  ALL ZERO -- this is a fttNoAdd=1 (control) file, or hits never reached the fit\n");

    gStyle->SetOptStat(0);
    TCanvas* cpk = new TCanvas("cpk", "cpk", 700, 500);
    gPad->SetGridy();
    TGraphErrors* gx = new TGraphErrors(4, sx, px, 0, ex);
    TGraphErrors* gy = new TGraphErrors(4, sx, py, 0, ey);
    double ymx = 0;
    for (int im = 0; im < 4; im++) {
        if (px[im] > ymx) ymx = px[im];
        if (py[im] > ymx) ymx = py[im];
    }
    if (ymx <= 0) ymx = 1;
    gx->SetTitle(Form("FTT pickup probability, %s tracks;sTGC station;hits used / track", type));
    gx->SetMarkerStyle(20); gx->SetMarkerColor(kAzure+2); gx->SetLineColor(kAzure+2);
    gy->SetMarkerStyle(21); gy->SetMarkerColor(kRed+1);   gy->SetLineColor(kRed+1);
    gx->GetYaxis()->SetRangeUser(0, ymx*1.35);
    gx->GetXaxis()->SetLimits(0.5, 4.5);
    gx->Draw("AP");
    gy->Draw("P same");
    TLegend* lg = new TLegend(0.15, 0.75, 0.55, 0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(gx, "V strips (measure x)", "p");
    lg->AddEntry(gy, "H strips (measure y)", "p");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
    tx.DrawLatex(0.15, 0.70, Form("%.0f %s tracks; control arm reads 0.0000 everywhere", ntr, type));
    cpk->SaveAs(out);
    printf("  wrote %s\n", out);
}
