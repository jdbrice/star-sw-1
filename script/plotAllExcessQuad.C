// plotAllExcessQuad.C
//
// Quadrant-split version of plotAllExcess.C -- plots disk1 (the only disk
// with an independently-confirmed online-QA cross-check), all 4 quadrants,
// both orientations, zoomed to the peak region with local-sideband
// background reference -- same style as plotAllExcess.C.
//
// Usage: root4star -l -b -q 'plotAllExcessQuad.C("fwd_blind_diag.root","outdir",1)'

double avgRangeAEQ(TH1F *h, double lo, double hi) {
    double s = 0; int c = 0;
    for (int ib = 1; ib <= h->GetNbinsX(); ib++) {
        double x = h->GetXaxis()->GetBinCenter(ib);
        if (x >= lo && x < hi) { s += h->GetBinContent(ib); c++; }
    }
    return s / c;
}

void plotAllExcessQuad(const char* infile = "fwd_blind_diag.root", const char* outdir = "20260721_allexcess_quad", int disk = 1) {
    gSystem->mkdir(outdir, true);
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", infile); return; }

    const char* qName[4] = {"A", "B", "C", "D"};
    const char* oriName[2] = {"V", "H"};

    for (int q = 0; q < 4; q++) {
        for (int io = 0; io < 2; io++) {
            TString histName = (io == 0) ? Form("hBlindDxAll_V_disk%d_quad%s", disk, qName[q])
                                          : Form("hBlindDyAll_H_disk%d_quad%s", disk, qName[q]);
            TH1F *h = (TH1F*)f->Get(histName);
            if (!h) { printf("Missing %s\n", histName.Data()); continue; }

            double peakVal = avgRangeAEQ(h, -0.6, 0.6);
            double loVal   = avgRangeAEQ(h, -3.0, -0.6);
            double hiVal   = avgRangeAEQ(h, 0.6, 3.0);
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

            TCanvas *c = new TCanvas(Form("c_%d_%d", q, io), "c", 700, 550);

            h->GetXaxis()->SetRangeUser(-5, 5);
            h->SetFillColor(kGray);
            h->SetLineColor(kGray+2);
            double ymin = TMath::Max(0.0, bgVal * 0.7);
            double ymax = TMath::Max(peakVal, hiVal) * 1.20;
            h->SetMinimum(ymin);
            h->SetMaximum(ymax);
            h->SetTitle(Form("disk%d quad%s %s (precise coord): excess = %.2fx (%.1f#sigma);"
                              "%s [cm];strips",
                              disk, qName[q], oriName[io], excess, sigma,
                              (io == 0) ? "dx" : "dy"));
            h->Draw("hist");

            TLine *bgLine = new TLine(-5, bgVal, 5, bgVal);
            bgLine->SetLineColor(kBlue+1);
            bgLine->SetLineStyle(2);
            bgLine->SetLineWidth(2);
            bgLine->Draw();

            TBox *loBox = new TBox(-3.0, ymin, -0.6, ymax);
            loBox->SetFillColorAlpha(kBlue, 0.08);
            loBox->SetLineColor(kBlue);
            loBox->SetLineStyle(3);
            loBox->Draw();
            TBox *hiBox = new TBox(0.6, ymin, 3.0, ymax);
            hiBox->SetFillColorAlpha(kBlue, 0.08);
            hiBox->SetLineColor(kBlue);
            hiBox->SetLineStyle(3);
            hiBox->Draw();
            TBox *peakBox = new TBox(-0.6, ymin, 0.6, ymax);
            peakBox->SetFillColorAlpha(kRed, 0.10);
            peakBox->SetLineColor(kRed);
            peakBox->SetLineStyle(3);
            peakBox->Draw();

            TLatex lt;
            lt.SetNDC();
            lt.SetTextSize(0.032);
            lt.DrawLatex(0.58, 0.83, Form("bg = %.0f #pm %.0f /bin", bgVal, bgErr));
            lt.DrawLatex(0.58, 0.78, Form("peak = %.0f /bin", peakVal));
            lt.DrawLatex(0.58, 0.73, Form("lo=%.0f  hi=%.0f", loVal, hiVal));

            c->SaveAs(Form("%s/excess_%s_disk%d_quad%s.png", outdir, oriName[io], disk, qName[q]));
            delete c;
        }
    }

    printf("Done -- wrote plots to %s/\n", outdir);
}
