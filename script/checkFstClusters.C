// checkFstClusters.C
//
// FST cluster properties in REAL DATA MuDst -- the FST counterpart of the FTT
// charge-sharing study (ftt_sim_hit_maker/charge_sharing.html), to inform the
// charge model in StFstSlowSimMaker (see FstSlowSim/index.html, "Limitations").
//
// IMPORTANT difference from the FTT study: the FstRawHit branch is EMPTY in
// this production (0 raw hits in 9210 events), so the per-strip ADC profile
// that the FTT study measured CANNOT be measured here. What survives into the
// MuDst is the cluster-level summary on StMuFstHit:
//     nRawHits, nRawHitsR, nRawHitsPhi   cluster size, split by direction
//     charge, chargeErr                  cluster total ADC
//     meanRStrip, meanPhiStrip           the reconstructed centroid
//     maxTimeBin, apv, disk/wedge/sensor
//
// That is still enough for the two things the sim needs: how big clusters are,
// and how much charge they carry. It is also enough to test what
// StFstScanRadiusClusterAlgo actually does to the position, which is not what
// the variable names suggest:
//
//   - meanRStrip is NOT a centroid. Step 1 sets meanRStrip = maxRStrip, the
//     LARGEST r-strip index in the cluster. So a cluster spanning nRawHitsR
//     strips is placed at its outer edge, biased outward by
//     (nRawHitsR-1)/2 * 2.875 cm relative to its own extent.
//   - meanPhiStrip is the phi index of the max-ADC strip (an integer) unless
//     step 2 merged phi-adjacent clusters, in which case it is charge-weighted
//     between them. So a fractional meanPhiStrip is a direct tag for "this
//     cluster went through the phi merge".
//   - step 2's merge condition is |dMeanRStrip| < 3.5, but r strips only span
//     4 values within a sensor, so the condition is ALWAYS true: any two
//     clusters in adjacent phi strips of the same sensor merge unconditionally,
//     and a run of occupied adjacent phi strips chains into one cluster.
//
// The histograms below measure how often each of those actually bites.
//
// Usage:
//   root4star -b -q 'script/checkFstClusters.C("<glob>", 0, "FstSlowSim/fstclus")'
//     nEvents<=0 means all.

void checkFstClusters(const char* glob, int nEvents = 0,
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

    // ---- histograms -------------------------------------------------------
    // region: 0 = inner sensor (sensor 0), 1 = outer sensors (1,2)
    TH1F* hN[2];  TH1F* hNR[2]; TH1F* hNP[2]; TH1F* hQ[2];
    for (int ia = 0; ia < 2; ia++){
        const char* tg = ia ? "outer" : "inner";
        hN [ia] = new TH1F(Form("hN_%d", ia),  Form("%s;cluster size nRawHits;clusters", tg), 21, -0.5, 20.5);
        hNR[ia] = new TH1F(Form("hNR_%d", ia), Form("%s;nRawHitsR;clusters", tg),              9, -0.5,  8.5);
        hNP[ia] = new TH1F(Form("hNP_%d", ia), Form("%s;nRawHitsPhi;clusters", tg),           21, -0.5, 20.5);
        hQ [ia] = new TH1F(Form("hQ_%d", ia),  Form("%s;cluster charge [ADC];clusters", tg), 200, 0, 4000);
        hN[ia]->SetDirectory(0); hNR[ia]->SetDirectory(0);
        hNP[ia]->SetDirectory(0); hQ[ia]->SetDirectory(0);
    }
    TH1F* hFracR = new TH1F("hFracR", ";frac part of meanRStrip;clusters",   110, -0.05, 1.05);
    TH1F* hFracP = new TH1F("hFracP", ";frac part of meanPhiStrip;clusters", 110, -0.05, 1.05);
    TH2F* hNRvNP = new TH2F("hNRvNP", ";nRawHitsR;nRawHitsPhi", 9, -0.5, 8.5, 13, -0.5, 12.5);
    TH1F* hMulti = new TH1F("hMulti", ";FST clusters per event;events", 100, 0, 400);
    TH1F* hTb    = new TH1F("hTb",    ";maxTimeBin;clusters", 12, -0.5, 11.5);
    TH1F* hRStr  = new TH1F("hRStr",  ";meanRStrip;clusters", 9, -0.5, 8.5);
    // charge vs size, to see whether big clusters are one big deposit or pileup
    TProfile* pQvN = new TProfile("pQvN", ";nRawHits;<cluster charge> [ADC]", 21, -0.5, 20.5);
    hFracR->SetDirectory(0); hFracP->SetDirectory(0); hNRvNP->SetDirectory(0);
    hMulti->SetDirectory(0); hTb->SetDirectory(0); hRStr->SetDirectory(0);
    pQvN->SetDirectory(0);

    // per-disk cluster size, disks 1-3
    TH1F* hNdisk[3];
    for (int ib = 0; ib < 3; ib++){
        hNdisk[ib] = new TH1F(Form("hNdisk_%d", ib), ";cluster size nRawHits;clusters", 21, -0.5, 20.5);
        hNdisk[ib]->SetDirectory(0);
    }

    double nClus = 0, nFracP = 0, nNR2 = 0, nNP2 = 0, nBig = 0;
    double sumRbias = 0;

    Long64_t ne = ch->GetEntries();
    if (nEvents > 0 && nEvents < ne) ne = nEvents;
    printf("events to read: %lld\n", ne);

    for (Long64_t iev = 0; iev < ne; iev++){
        ch->GetEntry(iev);
        if (!fstArr) continue;
        int nh = fstArr->GetEntriesFast();
        hMulti->Fill(nh);

        for (int ih = 0; ih < nh; ih++){
            StMuFstHit* hit = (StMuFstHit*) fstArr->UncheckedAt(ih);
            if (!hit) continue;

            int    sens = (int) hit->getSensor();
            int    disk = (int) hit->getDisk();
            int    nraw = (int) hit->getNRawHits();
            int    nrr  = (int) hit->getNRawHitsR();
            int    nrp  = (int) hit->getNRawHitsPhi();
            float  qq   = hit->getCharge();
            float  mrs  = hit->getMeanRStrip();
            float  mps  = hit->getMeanPhiStrip();

            int reg = (sens == 0) ? 0 : 1;
            hN [reg]->Fill(nraw);
            hNR[reg]->Fill(nrr);
            hNP[reg]->Fill(nrp);
            hQ [reg]->Fill(qq);
            if (disk >= 1 && disk <= 3) hNdisk[disk-1]->Fill(nraw);
            hNRvNP->Fill(nrr, nrp);
            hTb->Fill((int) hit->getMaxTimeBin());
            hRStr->Fill(mrs);
            pQvN->Fill(nraw, qq);

            double fr = mrs - floor(mrs);
            double fp = mps - floor(mps);
            hFracR->Fill(fr);
            hFracP->Fill(fp);

            nClus += 1;
            if (fp > 1e-4 && fp < 1 - 1e-4) nFracP += 1;
            if (nrr > 1) { nNR2 += 1; sumRbias += 0.5*(nrr - 1)*2.875; }
            if (nrp > 1) nNP2 += 1;
            if (nraw > 8) nBig += 1;
        }
    }

    // ---- numbers ----------------------------------------------------------
    printf("\n==============================================================\n");
    printf(" FST cluster properties, real data -- %.0f clusters\n", nClus);
    printf("==============================================================\n");
    printf("  clusters/event  mean %.1f   rms %.1f\n", hMulti->GetMean(), hMulti->GetRMS());

    printf("\n--- cluster size nRawHits ---\n");
    printf("  size |    inner            outer            all\n");
    double tIn = hN[0]->Integral(0, hN[0]->GetNbinsX()+1);
    double tOu = hN[1]->Integral(0, hN[1]->GetNbinsX()+1);
    for (int is = 1; is <= 10; is++){
        double ci = hN[0]->GetBinContent(hN[0]->FindBin(is));
        double co = hN[1]->GetBinContent(hN[1]->FindBin(is));
        printf("  %4d | %9.0f %6.2f%%  %9.0f %6.2f%%  %6.2f%%\n", is,
               ci, tIn>0?100.0*ci/tIn:0, co, tOu>0?100.0*co/tOu:0,
               nClus>0?100.0*(ci+co)/nClus:0);
    }
    printf("  mean size: inner %.3f   outer %.3f\n", hN[0]->GetMean(), hN[1]->GetMean());

    printf("\n--- size split by direction ---\n");
    printf("   n   |  nRawHitsR (in/out)    nRawHitsPhi (in/out)\n");
    for (int js = 1; js <= 6; js++){
        double ri = hNR[0]->GetBinContent(hNR[0]->FindBin(js));
        double ro = hNR[1]->GetBinContent(hNR[1]->FindBin(js));
        double pi_ = hNP[0]->GetBinContent(hNP[0]->FindBin(js));
        double po = hNP[1]->GetBinContent(hNP[1]->FindBin(js));
        printf("  %3d  | %6.2f%% %6.2f%%        %6.2f%% %6.2f%%\n", js,
               tIn>0?100.0*ri/tIn:0, tOu>0?100.0*ro/tOu:0,
               tIn>0?100.0*pi_/tIn:0, tOu>0?100.0*po/tOu:0);
    }

    printf("\n--- what the algorithm does to the position ---\n");
    printf("  clusters with nRawHitsR > 1   : %8.0f  (%.2f%%)\n", nNR2, 100.0*nNR2/nClus);
    printf("      these are placed at the OUTERMOST r strip, not the centroid.\n");
    printf("      mean outward bias over those clusters : %.3f cm\n", nNR2>0? sumRbias/nNR2 : 0);
    printf("      mean outward bias over ALL clusters   : %.3f cm\n", nClus>0? sumRbias/nClus : 0);
    printf("  clusters with nRawHitsPhi > 1 : %8.0f  (%.2f%%)   <- went through phi merge\n",
           nNP2, 100.0*nNP2/nClus);
    printf("  clusters with fractional meanPhiStrip: %8.0f (%.2f%%)\n", nFracP, 100.0*nFracP/nClus);
    printf("  clusters with nRawHits > 8    : %8.0f  (%.2f%%)   <- phi chaining\n",
           nBig, 100.0*nBig/nClus);
    printf("  frac part of meanRStrip: mean %.4f (0 => always an integer strip index)\n",
           hFracR->GetMean());

    printf("\n--- charge ---\n");
    printf("  inner: mean %.1f  rms %.1f\n", hQ[0]->GetMean(), hQ[0]->GetRMS());
    printf("  outer: mean %.1f  rms %.1f\n", hQ[1]->GetMean(), hQ[1]->GetRMS());
    printf("  <charge> vs size:\n");
    for (int ks = 1; ks <= 8; ks++){
        int bb = pQvN->FindBin(ks);
        if (pQvN->GetBinEntries(bb) <= 0) continue;
        printf("    n=%d : %8.1f  (%.0f clusters)\n", ks,
               pQvN->GetBinContent(bb), pQvN->GetBinEntries(bb));
    }

    printf("\n--- per disk, mean cluster size ---\n");
    for (int kd = 0; kd < 3; kd++)
        printf("   disk %d : %.3f   (%.0f clusters)\n", kd+1,
               hNdisk[kd]->GetMean(), hNdisk[kd]->Integral());

    // ---- plots ------------------------------------------------------------
    TCanvas* c1 = new TCanvas("c_fstclus", "", 1100, 700);
    c1->Divide(3,2);

    c1->cd(1); gPad->SetLogy();
    hN[0]->SetLineColor(kBlue+1); hN[0]->SetLineWidth(2); hN[0]->SetFillStyle(0);
    hN[1]->SetLineColor(kRed+1);  hN[1]->SetLineWidth(2); hN[1]->SetFillStyle(0);
    if (tIn>0) hN[0]->Scale(1.0/tIn);
    if (tOu>0) hN[1]->Scale(1.0/tOu);
    hN[0]->SetTitle("cluster size");
    hN[0]->GetYaxis()->SetTitle("fraction of clusters");
    hN[0]->Draw("hist"); hN[1]->Draw("hist same");
    TLegend* lg1 = new TLegend(0.5,0.72,0.88,0.86);
    lg1->SetBorderSize(0); lg1->SetFillStyle(0);
    lg1->AddEntry(hN[0], "inner sensor", "l");
    lg1->AddEntry(hN[1], "outer sensors", "l");
    lg1->Draw();

    c1->cd(2); gPad->SetLogy();
    hNR[0]->SetLineColor(kBlue+1); hNR[0]->SetLineWidth(2); hNR[0]->SetFillStyle(0);
    hNR[1]->SetLineColor(kRed+1);  hNR[1]->SetLineWidth(2); hNR[1]->SetFillStyle(0);
    if (tIn>0) hNR[0]->Scale(1.0/tIn);
    if (tOu>0) hNR[1]->Scale(1.0/tOu);
    hNR[0]->SetTitle("size along R");
    hNR[0]->GetYaxis()->SetTitle("fraction of clusters");
    hNR[0]->Draw("hist"); hNR[1]->Draw("hist same");

    c1->cd(3); gPad->SetLogy();
    hNP[0]->SetLineColor(kBlue+1); hNP[0]->SetLineWidth(2); hNP[0]->SetFillStyle(0);
    hNP[1]->SetLineColor(kRed+1);  hNP[1]->SetLineWidth(2); hNP[1]->SetFillStyle(0);
    if (tIn>0) hNP[0]->Scale(1.0/tIn);
    if (tOu>0) hNP[1]->Scale(1.0/tOu);
    hNP[0]->SetTitle("size along phi");
    hNP[0]->GetYaxis()->SetTitle("fraction of clusters");
    hNP[0]->Draw("hist"); hNP[1]->Draw("hist same");

    c1->cd(4); gPad->SetLogy();
    hQ[0]->SetLineColor(kBlue+1); hQ[0]->SetLineWidth(2); hQ[0]->SetFillStyle(0);
    hQ[1]->SetLineColor(kRed+1);  hQ[1]->SetLineWidth(2); hQ[1]->SetFillStyle(0);
    if (hQ[0]->Integral()>0) hQ[0]->Scale(1.0/hQ[0]->Integral());
    if (hQ[1]->Integral()>0) hQ[1]->Scale(1.0/hQ[1]->Integral());
    hQ[0]->SetTitle("cluster charge");
    hQ[0]->GetYaxis()->SetTitle("fraction of clusters");
    hQ[0]->Draw("hist"); hQ[1]->Draw("hist same");

    c1->cd(5);
    hRStr->SetLineColor(kBlack); hRStr->SetLineWidth(2); hRStr->SetFillStyle(0);
    hRStr->SetTitle("meanRStrip (= outermost strip of cluster)");
    hRStr->Draw("hist");

    c1->cd(6); gPad->SetLogz();
    hNRvNP->SetTitle("size R vs size phi");
    hNRvNP->Draw("colz");

    c1->SaveAs(Form("%s/fstClusterProps.png", outdir));

    TCanvas* c2 = new TCanvas("c_fstclus2", "", 1000, 400);
    c2->Divide(2,1);
    c2->cd(1); gPad->SetLogy();
    hMulti->SetLineColor(kBlack); hMulti->SetLineWidth(2); hMulti->SetFillStyle(0);
    hMulti->SetTitle("FST clusters per event");
    hMulti->Draw("hist");
    c2->cd(2); gPad->SetLogy();
    hFracP->SetLineColor(kBlack); hFracP->SetLineWidth(2); hFracP->SetFillStyle(0);
    hFracP->SetTitle("fractional part of meanPhiStrip");
    hFracP->Draw("hist");
    c2->SaveAs(Form("%s/fstClusterMulti.png", outdir));

    TFile* fo = new TFile(Form("%s/fstClusterHistos.root", outdir), "RECREATE");
    for (int iw = 0; iw < 2; iw++){ hN[iw]->Write(); hNR[iw]->Write(); hNP[iw]->Write(); hQ[iw]->Write(); }
    hFracR->Write(); hFracP->Write(); hNRvNP->Write(); hMulti->Write();
    hTb->Write(); hRStr->Write(); pQvN->Write();
    for (int ix = 0; ix < 3; ix++) hNdisk[ix]->Write();
    fo->Close();
    printf("\nwrote %s/fstClusterProps.png, %s/fstClusterMulti.png, %s/fstClusterHistos.root\n",
           outdir, outdir, outdir);
}
