// checkFstGapFixEffect.C
//
// Direct verification that the outer-sensor gap-offset fix
// (StFwdHitLoader::setApplyFstGapFix) does what it is supposed to, using hit
// POSITIONS only -- no residuals, no track fit, nothing that could absorb it.
//
// The signature is unmistakable. Reconstructed outer-sensor hit azimuth,
// measured relative to its own wedge centreline (dphiLocal):
//
// The hits are real; only the position assigned to them moves. Without the fix
// each outer sensor's hits are shifted 1 deg TOWARD the centreline:
//     -phi sensor: true [-15.383,-0.617] -> reconstructed [-14.383, +0.383]
//     +phi sensor: true [+0.617,+15.383] -> reconstructed [-0.383, +14.383]
// so |dphiLocal| < 0.383 deg receives hits from BOTH sensors -- a ~2x density
// BUMP, not an empty gap -- while 14.383 < |dphiLocal| < 15.383 is left EMPTY.
// Total hits are conserved: the bump is paid for by the missing edges.
//
//   gap fix OFF: ~2x density spike at |dphiLocal| < 0.383, edges empty,
//                band ends at +-14.383
//   gap fix ON:  |dphiLocal| < 0.617 EMPTY (the real 1 deg gap plus half a
//                pitch), edges populated, band ends at +-15.383, no spike
//
// Inner-sensor hits carry no gap term and must be identical either way, which
// is the control.
//
// Run it on the two alignment ntuples produced with the flag off and on.
//
// Usage:
//   root4star -b -q 'script/checkFstGapFixEffect.C("off.root","on.root","outdir")'

void checkFstGapFixEffect(const char* fOff, const char* fOn,
                          const char* outdir = "residual/data/alignment/gapfix"){

    gSystem->mkdir(outdir, kTRUE);
    gStyle->SetOptStat(0);

    // dphiLocal distributions: [0] = off, [1] = on;  inner and outer
    TH1F* hIn[2];
    TH1F* hOut[2];
    for (int ia = 0; ia < 2; ia++){
        // Bin width MUST be an exact multiple of the FST strip pitch
        // (kFstStripPitchPhi = 30/128 = 0.234375 deg) or the strips-per-bin
        // count beats against the binning and puts periodic spikes in the
        // figure that look like structure. One strip per bin, over an integer
        // number of strips (+-72 strips = +-16.875 deg), covers the full
        // +-15.383 deg occupancy with no aliasing.
        const int   kNStrip = 144;               // 72 each side
        const double kHalf  = 72.0 * 30.0/128.0; // 16.875 deg
        hIn [ia] = new TH1F(Form("hIn_%d",  ia), ";#delta#phi_{local} [deg];hits", kNStrip, -kHalf, kHalf);
        hOut[ia] = new TH1F(Form("hOut_%d", ia), ";#delta#phi_{local} [deg];hits", kNStrip, -kHalf, kHalf);
        hIn [ia]->SetDirectory(0);
        hOut[ia]->SetDirectory(0);
    }

    double diskZ[3]; diskZ[0] = 151.750; diskZ[1] = 165.248; diskZ[2] = 178.781;

    for (int ipass = 0; ipass < 2; ipass++){
        const char* fn = (ipass == 0) ? fOff : fOn;
        TChain* tr = new TChain("alignTree");
        if (tr->Add(fn) <= 0){ printf("no file matched %s\n", fn); return; }

        Int_t   bDet = 0;
        Float_t bHx = 0, bHy = 0, bHz = 0;
        tr->SetBranchAddress("detType", &bDet);
        tr->SetBranchAddress("hitX", &bHx);
        tr->SetBranchAddress("hitY", &bHy);
        tr->SetBranchAddress("hitZ", &bHz);

        Long64_t ne = tr->GetEntries();
        printf("%s : %lld rows\n", fn, ne);
        for (Long64_t ie = 0; ie < ne; ie++){
            tr->GetEntry(ie);
            if (bDet != 0) continue;
            int ok = 0;
            for (int db = 0; db < 3; db++) if (fabs(bHz - diskZ[db]) < 5.0) ok = 1;
            if (!ok) continue;

            double rr = sqrt(bHx*bHx + bHy*bHy);
            double p360 = atan2(bHy, bHx)*180.0/TMath::Pi();
            if (p360 < 0) p360 += 360.0;
            int sec = (int)(p360/30.0); if (sec > 11) sec = 11;
            double dLoc = p360 - (sec*30.0 + 15.0);

            if (rr < 16.5) hIn [ipass]->Fill(dLoc);
            else           hOut[ipass]->Fill(dLoc);
        }
        delete tr;
    }

    // ---- the numbers that matter ----
    printf("\n==============================================================\n");
    printf(" Outer-sensor hits: occupancy of the wedge centreline region\n");
    printf("==============================================================\n");
    printf(" OFF should show a ~2x density SPIKE in the overlap |dd|<0.383 and\n");
    printf(" empty edges; ON should show an empty core |dd|<0.617 and populated\n");
    printf(" edges. Density is per degree, normalised to the plateau 2<|dd|<10.\n\n");
    printf(" sample | overlap |dd|<0.383 | edge 14.4<|dd|<15.4 |  plateau  | total\n");
    printf("        |  n    density/plat |   n    density/plat |  dens/deg |\n");
    for (int ib = 0; ib < 2; ib++){
        double ovl = 0, edge = 0, plat = 0, tot = 0;
        double wOvl = 0, wEdge = 0, wPlat = 0;
        double bw = hOut[ib]->GetBinWidth(1);
        for (int jb = 1; jb <= hOut[ib]->GetNbinsX(); jb++){
            double xc = fabs(hOut[ib]->GetBinCenter(jb));
            double vv = hOut[ib]->GetBinContent(jb);
            tot += vv;
            if (xc < 0.383)                 { ovl  += vv; wOvl  += bw; }
            if (xc > 14.4 && xc < 15.4)     { edge += vv; wEdge += bw; }
            if (xc > 2.0  && xc < 10.0)     { plat += vv; wPlat += bw; }
        }
        double dPlat = (wPlat > 0) ? plat/wPlat : 0;
        double dOvl  = (wOvl  > 0) ? ovl /wOvl  : 0;
        double dEdge = (wEdge > 0) ? edge/wEdge : 0;
        printf("   %s  | %6.0f   %8.2f  | %6.0f  %8.2f   | %9.1f | %7.0f\n",
               ib == 0 ? "OFF" : "ON ",
               ovl,  dPlat > 0 ? dOvl /dPlat : 0,
               edge, dPlat > 0 ? dEdge/dPlat : 0,
               dPlat, tot);
    }

    printf("\n Inner-sensor control (must be unchanged):\n");
    printf(" sample |  core |dd|<0.6  |  total\n");
    for (int ic = 0; ic < 2; ic++){
        double core2 = 0, tot2 = 0;
        for (int jc = 1; jc <= hIn[ic]->GetNbinsX(); jc++){
            double xc2 = fabs(hIn[ic]->GetBinCenter(jc));
            double vv2 = hIn[ic]->GetBinContent(jc);
            tot2 += vv2;
            if (xc2 < 0.6) core2 += vv2;
        }
        printf("   %s  |  %9.0f     | %8.0f\n", ic == 0 ? "OFF" : "ON ", core2, tot2);
    }

    // ---- plot ----
    // cache the histograms so the figure can be restyled without another
    // full pass over the campaign (58M rows per side)
    TFile* fout = new TFile(Form("%s/gapFixHistos.root", outdir), "RECREATE");
    for (int ih = 0; ih < 2; ih++){ hIn[ih]->Write(); hOut[ih]->Write(); }
    fout->Close();

    TCanvas* cv = new TCanvas("c_gapfix", "", 1000, 420);
    cv->Divide(2,1);
    cv->cd(1);
    // no fill -- a filled ON histogram drawn over OFF hides exactly the
    // comparison the plot exists to show
    hOut[0]->SetLineColor(kRed+1);  hOut[0]->SetLineWidth(2); hOut[0]->SetFillStyle(0);
    hOut[1]->SetLineColor(kBlue+1); hOut[1]->SetLineWidth(2); hOut[1]->SetFillStyle(0);
    hOut[0]->SetTitle("OUTER sensors (r>16.5 cm)");
    hOut[0]->Draw("hist");
    hOut[1]->Draw("hist same");
    TLegend* lg = new TLegend(0.35,0.78,0.65,0.9);
    lg->AddEntry(hOut[0], "gap fix OFF", "l");
    lg->AddEntry(hOut[1], "gap fix ON",  "l");
    lg->Draw();
    cv->cd(2);
    hIn[0]->SetLineColor(kRed+1);  hIn[0]->SetLineWidth(2); hIn[0]->SetFillStyle(0);
    hIn[1]->SetLineColor(kBlue+1); hIn[1]->SetLineWidth(2); hIn[1]->SetFillStyle(0);
    hIn[0]->SetTitle("INNER sensor (r<16.5 cm) -- control");
    hIn[0]->Draw("hist");
    hIn[1]->Draw("hist same");
    cv->SaveAs(Form("%s/fstGapFixEffect.png", outdir));
    printf("\nWrote %s/fstGapFixEffect.png\n", outdir);
}
