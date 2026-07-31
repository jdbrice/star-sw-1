// plotFstGapResidual.C
//
// Residual SHAPES with the outer-sensor gap fix off and on -- does the
// distribution actually change, or only its mean?
//
// chi2/ndf says the fit improves for outer hits (see
// gapfix_test/CONCLUSIONS.txt section 6), but chi2 is a single number. This
// draws the unbiased r*dphi residual itself, outer and inner separately, and
// quantifies the shape three ways:
//   rms   over the plotted range
//   core  Gaussian sigma from a fit restricted to +-0.3 cm
//   tail  fraction of rows beyond +-0.3 cm
// A real geometric correction should narrow the OUTER core and/or cut its
// tail, and leave the INNER control alone.
//
// Also caches the histograms so the figure can be restyled without another
// pass over ~58M rows per side.
//
// Usage:
//   root4star -b -q 'script/plotFstGapResidual.C("<offGlob>","<onGlob>","fstgap")'

void plotFstGapResidual(const char* offGlob, const char* onGlob,
                        const char* outdir = "fstgap"){

    gSystem->mkdir(outdir, kTRUE);
    gStyle->SetOptStat(0);

    // [pass][region] flattened: pass*2 + region, region 0=inner 1=outer
    TH1F* h[4];
    for (int ia = 0; ia < 4; ia++){
        h[ia] = new TH1F(Form("hres_%d", ia), ";r#Delta#phi [cm];rows", 200, -1.0, 1.0);
        h[ia]->SetDirectory(0);
    }
    double nAll[4]; for (int ib = 0; ib < 4; ib++) nAll[ib] = 0;

    double diskZ[3]; diskZ[0] = 151.750; diskZ[1] = 165.248; diskZ[2] = 178.781;

    for (int ip = 0; ip < 2; ip++){
        TChain* tr = new TChain("alignTree");
        if (tr->Add(ip == 0 ? offGlob : onGlob) <= 0){ printf("no files\n"); return; }

        Int_t   bDet = 0;
        Float_t bHx = 0, bHy = 0, bHz = 0, bPx = 0, bPy = 0;
        tr->SetBranchAddress("detType", &bDet);
        tr->SetBranchAddress("hitX",  &bHx);
        tr->SetBranchAddress("hitY",  &bHy);
        tr->SetBranchAddress("hitZ",  &bHz);
        tr->SetBranchAddress("projX", &bPx);
        tr->SetBranchAddress("projY", &bPy);

        Long64_t ne = tr->GetEntries();
        printf("pass %d: %lld rows\n", ip, ne);
        for (Long64_t ie = 0; ie < ne; ie++){
            tr->GetEntry(ie);
            if (bDet != 0) continue;
            int ok = 0;
            for (int db = 0; db < 3; db++) if (fabs(bHz - diskZ[db]) < 5.0) ok = 1;
            if (!ok) continue;

            double rr   = sqrt(bHx*bHx + bHy*bHy);
            double phih = atan2(bHy, bHx);
            double phip = atan2(bPy, bPx);
            double dph  = phih - phip;
            while (dph >  TMath::Pi()) dph -= 2*TMath::Pi();
            while (dph < -TMath::Pi()) dph += 2*TMath::Pi();
            double rdp = rr * dph;

            int reg = (rr > 16.5) ? 1 : 0;
            nAll[ip*2 + reg] += 1;              // counts EVERYTHING, for the tail fraction
            h[ip*2 + reg]->Fill(rdp);
        }
        delete tr;
    }

    printf("\n=================================================================\n");
    printf(" Residual shape, r*dphi, gap fix off vs on\n");
    printf("=================================================================\n");
    printf(" region  pass |    rows      rms     core sigma   tail>0.3cm\n");
    for (int ir = 1; ir >= 0; ir--){
        for (int jp = 0; jp < 2; jp++){
            int ix = jp*2 + ir;
            double in3 = h[ix]->Integral(h[ix]->FindBin(-0.2999), h[ix]->FindBin(0.2999));
            double tot = nAll[ix];
            h[ix]->Fit("gaus", "QN0", "", -0.3, 0.3);
            TF1* g = (TF1*) gROOT->GetFunction("gaus");
            printf("  %s   %s | %9.0f  %7.4f    %7.4f     %6.2f%%\n",
                   ir ? "OUTER" : "inner", jp ? "ON " : "OFF",
                   tot, h[ix]->GetRMS(), g ? g->GetParameter(2) : 0,
                   tot > 0 ? 100.0*(tot - in3)/tot : 0);
        }
    }

    TFile* fo = new TFile(Form("%s/gapResidualHistos.root", outdir), "RECREATE");
    for (int ic = 0; ic < 4; ic++) h[ic]->Write();
    fo->Close();

    // ---- figure: lines only, no fill (a filled ON hides OFF underneath) ----
    TCanvas* cv = new TCanvas("c_gapres", "", 1000, 420);
    cv->Divide(2,1);
    for (int ipad = 0; ipad < 2; ipad++){
        cv->cd(ipad+1);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogy();
        int reg = (ipad == 0) ? 1 : 0;      // outer first
        TH1F* a = h[0*2 + reg];
        TH1F* b = h[1*2 + reg];
        // normalise so the comparison is of SHAPE, not of how many rows survived
        if (a->Integral() > 0) a->Scale(1.0/a->Integral());
        if (b->Integral() > 0) b->Scale(1.0/b->Integral());
        a->SetLineColor(kRed+1);  a->SetLineWidth(2); a->SetFillStyle(0);
        b->SetLineColor(kBlue+1); b->SetLineWidth(2); b->SetFillStyle(0);
        a->SetTitle(reg ? "OUTER sensors (r>16.5 cm)" : "INNER sensor (r<16.5 cm) -- control");
        a->GetYaxis()->SetTitle("fraction of rows / bin");
        a->Draw("hist");
        b->Draw("hist same");
        TLegend* lg = new TLegend(0.16, 0.75, 0.48, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        lg->AddEntry(a, "gap fix OFF", "l");
        lg->AddEntry(b, "gap fix ON",  "l");
        lg->Draw();
    }
    cv->SaveAs(Form("%s/fstGapResidual.png", outdir));
    printf("\nwrote %s/fstGapResidual.png\n", outdir);
}
