// plotAllExcess.C
//
// Visualizes the section-1 check from proposal_next_step_20260721.txt: for
// each disk/orientation, zoom the "All candidates" (precise-coordinate)
// histogram into the peak region, mark the peak window and the two local
// sideband windows used to estimate the local background, and annotate the
// excess ratio + significance -- the same numbers script/checkAllExcess.C
// prints, shown visually.
//
// Usage: root4star -l -b -q 'plotAllExcess.C("fwd_blind_diag.root","outdir")'

double avgRangeAE(TH1F *h, double lo, double hi) {
    double s = 0; int c = 0;
    for (int ib = 1; ib <= h->GetNbinsX(); ib++) {
        double x = h->GetXaxis()->GetBinCenter(ib);
        if (x >= lo && x < hi) { s += h->GetBinContent(ib); c++; }
    }
    return s / c;
}

void plotAllExcess(const char* infile = "fwd_blind_diag.root", const char* outdir = "20260721_allexcess") {
    gSystem->mkdir(outdir, true);
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", infile); return; }

    const int nDisk = 4;
    const char* oriName[2] = {"V", "H"};

    for (int d = 0; d < nDisk; d++) {
        for (int io = 0; io < 2; io++) {
            TString histName = (io == 0) ? Form("hBlindDxAll_V_disk%d", d)
                                          : Form("hBlindDyAll_H_disk%d", d);
            TH1F *h = (TH1F*)f->Get(histName);
            if (!h) { printf("Missing %s\n", histName.Data()); continue; }

            double peakVal = avgRangeAE(h, -0.6, 0.6);
            double loVal   = avgRangeAE(h, -3.0, -0.6);
            double hiVal   = avgRangeAE(h, 0.6, 3.0);
            double bgVal   = (loVal + hiVal) / 2.0;
            // approximate per-bin sideband count for the error (matches checkAllExcess.C)
            double loSum = 0, hiSum = 0; int loN = 0, hiN = 0;
            for (int ib = 1; ib <= h->GetNbinsX(); ib++) {
                double x = h->GetXaxis()->GetBinCenter(ib);
                if (x >= -3.0 && x < -0.6) { loSum += h->GetBinContent(ib); loN++; }
                if (x >= 0.6 && x < 3.0)   { hiSum += h->GetBinContent(ib); hiN++; }
            }
            double bgErr = sqrt(loSum + hiSum) / (loN + hiN);
            double excess = peakVal / bgVal;
            double sigma = (bgErr > 0) ? (peakVal - bgVal) / bgErr : 0;

            TCanvas *c = new TCanvas(Form("c_%d_%d", d, io), "c", 700, 550);

            h->GetXaxis()->SetRangeUser(-5, 5);
            h->SetFillColor(kGray);
            h->SetLineColor(kGray+2);
            double ymin = bgVal * 0.85;
            double ymax = TMath::Max(peakVal, hiVal) * 1.20;
            h->SetMinimum(ymin);
            h->SetMaximum(ymax);
            h->SetTitle(Form("disk%d %s (precise coord): local excess = %.2fx (%.1f#sigma);"
                              "%s [cm];strips",
                              d, oriName[io], excess, sigma,
                              (io == 0) ? "dx" : "dy"));
            h->Draw("hist");

            // reference line at local background level
            TLine *bgLine = new TLine(-5, bgVal, 5, bgVal);
            bgLine->SetLineColor(kBlue+1);
            bgLine->SetLineStyle(2);
            bgLine->SetLineWidth(2);
            bgLine->Draw();

            // shaded/boxed windows: sidebands (blue) and peak (red)
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

            TLegend *leg = new TLegend(0.14, 0.76, 0.55, 0.89);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(bgLine, "local bg level (sideband avg)", "l");
            leg->AddEntry(peakBox, "peak window |x|<0.6cm", "f");
            leg->AddEntry(loBox, "sideband windows 0.6-3cm", "f");
            leg->Draw();

            TLatex lt;
            lt.SetNDC();
            lt.SetTextSize(0.032);
            lt.DrawLatex(0.58, 0.83, Form("bg = %.0f #pm %.0f /bin", bgVal, bgErr));
            lt.DrawLatex(0.58, 0.78, Form("peak = %.0f /bin", peakVal));
            lt.DrawLatex(0.58, 0.73, Form("lo=%.0f  hi=%.0f", loVal, hiVal));

            c->SaveAs(Form("%s/excess_%s_disk%d.png", outdir, oriName[io], d));
            delete c;
        }
    }

    printf("Done -- wrote plots to %s/\n", outdir);
}
