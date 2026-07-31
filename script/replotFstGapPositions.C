// replotFstGapPositions.C
//
// Re-renders the hit-position figure from the cache written by
// checkFstGapFixEffect.C. No second pass over ~58M rows.
//
// Note on styling: SetFillStyle(0) alone was not enough here -- the histograms
// still came out solid. Forcing SetFillColor(0) as well, and drawing with
// "hist" on a hollow object, gives outline-only curves. Both histograms must
// be hollow or the one drawn second hides the other, which defeats the plot.
//
// Usage: root4star -b -q 'script/replotFstGapPositions.C'
void replotFstGapPositions(const char* dir = "fstgap"){
    gStyle->SetOptStat(0);
    TFile* f = TFile::Open(Form("%s/gapFixHistos.root", dir));
    if (!f){ printf("no cache\n"); return; }
    TH1F* hi[2]; TH1F* ho[2];
    for (int i = 0; i < 2; i++){
        hi[i] = (TH1F*) f->Get(Form("hIn_%d",  i));
        ho[i] = (TH1F*) f->Get(Form("hOut_%d", i));
        if (!hi[i] || !ho[i]){ printf("missing hist %d\n", i); return; }
        hi[i]->SetDirectory(0); ho[i]->SetDirectory(0);
    }
    // Two things had to be fixed to make this figure readable.
    //
    // 1. The apparent solid fill was NOT a fill. The source histogram has
    //    0.12 deg bins while the FST strip pitch is 0.234375 deg, so strips
    //    land in alternate bins and the curve plunges to ~0 between every
    //    pair -- a comb that reads as a filled block at this scale. Rebinning
    //    to ~1 deg spans several strips and removes it.
    // 2. TGraph is used rather than TH1::Draw("hist") so nothing can fill.
    // source is already binned one strip per bin, so rebin by an integer
    // number of strips only -- 4 strips = 0.9375 deg, still resolves the gap
    const int kRebin = 4;
    for (int ir = 0; ir < 2; ir++){ hi[ir]->Rebin(kRebin); ho[ir]->Rebin(kRebin); }
    TCanvas* c = new TCanvas("c_pos2", "", 1000, 420);
    c->Divide(2,1);

    for (int ip = 0; ip < 2; ip++){
        TH1F* a0 = ip ? hi[0] : ho[0];      // OFF
        TH1F* a1 = ip ? hi[1] : ho[1];      // ON
        int nb = a0->GetNbinsX();
        TGraph* g0 = new TGraph();
        TGraph* g1 = new TGraph();
        double ymin = 1e30, ymax = -1e30;
        for (int ib = 1; ib <= nb; ib++){
            double x = a0->GetBinCenter(ib);
            if (fabs(x) > 17) continue;
            double y0 = a0->GetBinContent(ib), y1 = a1->GetBinContent(ib);
            g0->SetPoint(g0->GetN(), x, y0);
            g1->SetPoint(g1->GetN(), x, y1);
            if (y0 > 0 && y0 < ymin) ymin = y0;
            if (y1 > 0 && y1 < ymin) ymin = y1;
            if (y0 > ymax) ymax = y0;
            if (y1 > ymax) ymax = y1;
        }
        c->cd(ip+1);
        gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14);
        // start just below the occupied band, not at zero -- the interesting
        // variation is a few percent of the plateau
        TH1F* fr = gPad->DrawFrame(-17, 0, 17, ymax*1.12);
        fr->GetXaxis()->SetTitle("#delta#phi_{local} [deg]");
        fr->GetYaxis()->SetTitle("hits");
        fr->GetYaxis()->SetTitleOffset(1.6);
        fr->SetTitle(ip ? "INNER sensor (r<16.5 cm) -- control"
                        : "OUTER sensors (r>16.5 cm)");
        g0->SetLineColor(kRed+1);  g0->SetLineWidth(2);
        g1->SetLineColor(kBlue+1); g1->SetLineWidth(2);
        g0->Draw("L same");
        g1->Draw("L same");
        if (ip == 0){
            TLegend* lg = new TLegend(0.33,0.16,0.68,0.31);
            lg->SetBorderSize(0); lg->SetFillStyle(0);
            lg->AddEntry(g0, "gap fix OFF", "l");
            lg->AddEntry(g1, "gap fix ON",  "l");
            lg->Draw();
        }
    }

    c->SaveAs(Form("%s/fstGapFixEffect.png", dir));
    printf("wrote %s/fstGapFixEffect.png\n", dir);
}
