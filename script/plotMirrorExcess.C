// plotMirrorExcess.C
//
// Visualizes the X<->Y mirror test (2026-07-22, see FwdTracker.h comments
// near hBlindDxAll_V_mirror/hBlindDyAll_H_mirror): for each disk, plots the
// ORIGINAL precise-coordinate "All candidates" histogram side by side with
// its MIRRORED counterpart (local x<->y swapped, then the correct
// per-quadrant/per-axis offset+sign re-applied -- NOT a naive global-XY
// swap, which would be physically incoherent given the asymmetric
// per-quadrant calibration).
//
// Usage: root4star -l -b -q 'plotMirrorExcess.C("fwd_blind_diag.root","outdir")'

double avgRangeME(TH1F *h, double lo, double hi) {
    double s = 0; int c = 0;
    for (int ib = 1; ib <= h->GetNbinsX(); ib++) {
        double x = h->GetXaxis()->GetBinCenter(ib);
        if (x >= lo && x < hi) { s += h->GetBinContent(ib); c++; }
    }
    return s / c;
}

void drawOnePanel(TH1F *h, const char* label) {
    double peakVal = avgRangeME(h, -0.6, 0.6);
    double loVal   = avgRangeME(h, -3.0, -0.6);
    double hiVal   = avgRangeME(h, 0.6, 3.0);
    double bgVal   = (loVal + hiVal) / 2.0;
    double loSum = 0, hiSum = 0; int loN = 0, hiN = 0;
    for (int ib = 1; ib <= h->GetNbinsX(); ib++) {
        double x = h->GetXaxis()->GetBinCenter(ib);
        if (x >= -3.0 && x < -0.6) { loSum += h->GetBinContent(ib); loN++; }
        if (x >= 0.6 && x < 3.0)   { hiSum += h->GetBinContent(ib); hiN++; }
    }
    double bgErr = sqrt(loSum + hiSum) / (loN + hiN);
    double excess = peakVal / bgVal;
    double sigma = (bgErr > 0) ? (peakVal - bgVal) / bgErr : 0;

    h->GetXaxis()->SetRangeUser(-5, 5);
    h->SetFillColor(kGray);
    h->SetLineColor(kGray+2);
    double ymin = TMath::Min(bgVal * 0.85, peakVal * 0.85);
    double ymax = TMath::Max(TMath::Max(peakVal, hiVal), bgVal) * 1.15;
    h->SetMinimum(ymin);
    h->SetMaximum(ymax);
    h->SetTitle(Form("%s: excess = %.3fx (%.1f#sigma);dx or dy [cm];strips", label, excess, sigma));
    h->Draw("hist");

    TLine *bgLine = new TLine(-5, bgVal, 5, bgVal);
    bgLine->SetLineColor(kBlue+1);
    bgLine->SetLineStyle(2);
    bgLine->SetLineWidth(2);
    bgLine->Draw();
}

void plotMirrorExcess(const char* infile = "fwd_blind_diag.root", const char* outdir = "20260722_mirrorexcess") {
    gSystem->mkdir(outdir, true);
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", infile); return; }

    for (int d = 0; d < 4; d++) {
        const char* names[4]  = {Form("hBlindDxAll_V_disk%d", d), Form("hBlindDxAll_V_mirror_disk%d", d),
                                  Form("hBlindDyAll_H_disk%d", d), Form("hBlindDyAll_H_mirror_disk%d", d)};
        const char* labels[4] = {Form("disk%d V (ORIGINAL)", d), Form("disk%d V (MIRRORED)", d),
                                  Form("disk%d H (ORIGINAL)", d), Form("disk%d H (MIRRORED)", d)};
        const char* fnames[4] = {"origV", "mirV", "origH", "mirH"};

        for (int i = 0; i < 4; i++) {
            TH1F *h = (TH1F*)f->Get(names[i]);
            if (!h) { printf("Missing %s\n", names[i]); continue; }
            TCanvas *c = new TCanvas(Form("c_%d_%d", d, i), "c", 700, 550);
            drawOnePanel(h, labels[i]);
            c->SaveAs(Form("%s/%s_disk%d.png", outdir, fnames[i], d));
            delete c;
        }
    }
    printf("Done -- wrote plots to %s/\n", outdir);
}
