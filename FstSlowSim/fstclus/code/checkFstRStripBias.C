// checkFstRStripBias.C
//
// Direct, assumption-free demonstration of what StFstScanRadiusClusterAlgo's
// "meanRStrip = maxRStrip" rule does to real data.
//
// The rule: a cluster spanning several r strips is assigned the LARGEST r-strip
// index it contains, not the charge centroid. If that bites, then multi-strip
// clusters must pile up at the OUTERMOST r strip of their sensor relative to
// single-strip clusters -- and single-strip clusters, which cannot be affected
// by the rule, are the control.
//
// Each sensor has exactly 4 r strips (kFstNumRStripsPerSensor), so the in-sensor
// index is 0..3 (inner sensor rStrip 0-3; outer sensors rStrip 4-7 minus 4).
// A cluster with nRawHitsR = 4 is forced to index 3 by construction; the
// interesting cases are 2 and 3.
//
// Usage:
//   root4star -b -q 'script/checkFstRStripBias.C("<glob>", 0, "FstSlowSim/fstclus")'

void checkFstRStripBias(const char* glob, int nEvents = 0,
                        const char* outdir = "FstSlowSim/fstclus"){

    gSystem->mkdir(outdir, kTRUE);
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    gStyle->SetOptStat(0);

    TChain* ch = new TChain("MuDst");
    if (ch->Add(glob) <= 0){ printf("no files matched %s\n", glob); return; }

    TClonesArray* fstArr = 0;
    ch->SetBranchStatus("*", 0);
    ch->SetBranchStatus("FstHit*", 1);
    ch->SetBranchAddress("FstHit", &fstArr);

    // [region][nRawHitsR-1] -> in-sensor r index 0..3
    TH1F* hIdx[2][4];
    for (int ia = 0; ia < 2; ia++){
        for (int ib = 0; ib < 4; ib++){
            hIdx[ia][ib] = new TH1F(Form("hIdx_%d_%d", ia, ib),
                                    ";in-sensor r index;clusters", 4, -0.5, 3.5);
            hIdx[ia][ib]->SetDirectory(0);
        }
    }
    // radial position actually assigned, single vs multi strip
    TH1F* hRpos[2];
    for (int ic = 0; ic < 2; ic++){
        hRpos[ic] = new TH1F(Form("hRpos_%d", ic), ";assigned r [cm];clusters", 80, 4, 30);
        hRpos[ic]->SetDirectory(0);
    }
    double rStart[8];
    rStart[0]=5.000; rStart[1]=7.875; rStart[2]=10.750; rStart[3]=13.625;
    rStart[4]=16.500; rStart[5]=19.375; rStart[6]=22.250; rStart[7]=25.125;

    double cnt[2][4]; double sumIdx[2][4];
    for (int id = 0; id < 2; id++) for (int ie = 0; ie < 4; ie++){ cnt[id][ie]=0; sumIdx[id][ie]=0; }

    Long64_t ne = ch->GetEntries();
    if (nEvents > 0 && nEvents < ne) ne = nEvents;
    printf("events to read: %lld\n", ne);

    for (Long64_t iev = 0; iev < ne; iev++){
        ch->GetEntry(iev);
        if (!fstArr) continue;
        int nh = fstArr->GetEntriesFast();
        for (int ih = 0; ih < nh; ih++){
            StMuFstHit* hit = (StMuFstHit*) fstArr->UncheckedAt(ih);
            if (!hit) continue;
            int sens = (int) hit->getSensor();
            int nrr  = (int) hit->getNRawHitsR();
            int mrs  = (int) hit->getMeanRStrip();
            if (nrr < 1 || nrr > 4) continue;
            if (mrs < 0 || mrs > 7) continue;

            int reg = (sens == 0) ? 0 : 1;
            int idx = (reg == 0) ? mrs : mrs - 4;
            if (idx < 0 || idx > 3) continue;

            hIdx[reg][nrr-1]->Fill(idx);
            cnt[reg][nrr-1]    += 1;
            sumIdx[reg][nrr-1] += idx;

            double rr = rStart[mrs] + 0.5*2.875;
            hRpos[(nrr == 1) ? 0 : 1]->Fill(rr);
        }
    }

    printf("\n=================================================================\n");
    printf(" Where multi-r-strip clusters get placed  (in-sensor r index 0..3)\n");
    printf("=================================================================\n");
    for (int ir = 0; ir < 2; ir++){
        printf("\n  %s sensor%s\n", ir ? "OUTER" : "INNER", ir ? "s" : "");
        printf("   nRawHitsR |   idx0    idx1    idx2    idx3  |  <idx>   clusters\n");
        for (int jn = 0; jn < 4; jn++){
            if (cnt[ir][jn] <= 0) continue;
            printf("       %d     |", jn+1);
            for (int kx = 0; kx < 4; kx++){
                double cc = hIdx[ir][jn]->GetBinContent(kx+1);
                printf(" %6.2f%%", 100.0*cc/cnt[ir][jn]);
            }
            printf("  | %5.3f  %9.0f\n", sumIdx[ir][jn]/cnt[ir][jn], cnt[ir][jn]);
        }
        if (cnt[ir][0] > 0 && cnt[ir][1] > 0)
            printf("   shift of <idx>, nRawHitsR=2 vs 1 : %+.3f strips = %+.3f cm\n",
                   sumIdx[ir][1]/cnt[ir][1] - sumIdx[ir][0]/cnt[ir][0],
                   (sumIdx[ir][1]/cnt[ir][1] - sumIdx[ir][0]/cnt[ir][0])*2.875);
    }

    printf("\n  If the clusters were r-ADJACENT pairs, a centroid rule would put\n");
    printf("  them at idx+0.5 and the <idx> shift vs single-strip would be about\n");
    printf("  +0.5 strips; the outermost rule doubles that to about +1.0.\n");

    TCanvas* cv = new TCanvas("c_rbias", "", 1000, 420);
    cv->Divide(2,1);
    for (int ip = 0; ip < 2; ip++){
        cv->cd(ip+1);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
        int col[4]; col[0]=kBlack; col[1]=kRed+1; col[2]=kBlue+1; col[3]=kGreen+2;
        double mx = 0;
        for (int jq = 0; jq < 4; jq++){
            if (cnt[ip][jq] <= 0) continue;
            hIdx[ip][jq]->Scale(1.0/cnt[ip][jq]);
            if (hIdx[ip][jq]->GetMaximum() > mx) mx = hIdx[ip][jq]->GetMaximum();
        }
        TLegend* lg = new TLegend(0.16,0.68,0.52,0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        int first = 1;
        for (int jr = 0; jr < 4; jr++){
            if (cnt[ip][jr] <= 0) continue;
            hIdx[ip][jr]->SetLineColor(col[jr]);
            hIdx[ip][jr]->SetLineWidth(2);
            hIdx[ip][jr]->SetFillStyle(0);
            hIdx[ip][jr]->SetMaximum(1.15*mx);
            hIdx[ip][jr]->SetMinimum(0);
            hIdx[ip][jr]->SetTitle(ip ? "OUTER sensors" : "INNER sensor");
            hIdx[ip][jr]->GetYaxis()->SetTitle("fraction of clusters");
            hIdx[ip][jr]->Draw(first ? "hist" : "hist same");
            lg->AddEntry(hIdx[ip][jr], Form("nRawHitsR = %d", jr+1), "l");
            first = 0;
        }
        lg->Draw();
    }
    cv->SaveAs(Form("%s/fstRStripBias.png", outdir));

    TCanvas* cv2 = new TCanvas("c_rpos", "", 620, 440);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    for (int iu = 0; iu < 2; iu++) if (hRpos[iu]->Integral() > 0) hRpos[iu]->Scale(1.0/hRpos[iu]->Integral());
    hRpos[0]->SetLineColor(kBlue+1); hRpos[0]->SetLineWidth(2); hRpos[0]->SetFillStyle(0);
    hRpos[1]->SetLineColor(kRed+1);  hRpos[1]->SetLineWidth(2); hRpos[1]->SetFillStyle(0);
    hRpos[0]->SetTitle("assigned radius");
    hRpos[0]->GetYaxis()->SetTitle("fraction of clusters");
    hRpos[0]->Draw("hist"); hRpos[1]->Draw("hist same");
    TLegend* lg2 = new TLegend(0.45,0.72,0.88,0.86);
    lg2->SetBorderSize(0); lg2->SetFillStyle(0);
    lg2->AddEntry(hRpos[0], "nRawHitsR = 1", "l");
    lg2->AddEntry(hRpos[1], "nRawHitsR > 1", "l");
    lg2->Draw();
    cv2->SaveAs(Form("%s/fstRPos.png", outdir));

    printf("\nwrote %s/fstRStripBias.png, %s/fstRPos.png\n", outdir, outdir);
}
