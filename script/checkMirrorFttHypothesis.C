// checkMirrorFttHypothesis.C
//
// Tests the hypothesis that the FST is mirrored within each wedge (reco vs
// reality) and that this is why FTT hit pickup is inefficient.
//
// The hypothesis: a mirror maps local phi -> -local phi, so it barely moves a
// hit at the wedge CENTRE and moves it by up to 2*15 deg at the wedge edge.
// FST-only seeding would survive everywhere, because all disks and wedges are
// mirrored consistently and three mirrored FST hits still look like a track.
// But the FTT is a different detector and is not mirrored, so extrapolating a
// mirrored FST track to the FTT would miss by ~2*delta. Prediction: tracks
// that DO pick up FTT hits should have their FST hits piled up near
// dphi_local = 0, and tracks with no FTT hits should be spread across the
// whole wedge.
//
// Null hypothesis (no mirror): the two distributions are the same shape, both
// roughly uniform across the wedge.
//
// Rows in the alignment tree are grouped per track, so tracks are accumulated
// by watching for a change in (run, event, trackId).
//
// Usage:
//   root4star -b -q 'script/checkMirrorFttHypothesis.C("<glob>", "outdir")'

void checkMirrorFttHypothesis(const char* glob, const char* outdir = "fstgap"){

    gSystem->mkdir(outdir, kTRUE);
    gStyle->SetOptStat(0);

    TChain* tr = new TChain("alignTree");
    if (tr->Add(glob) <= 0){ printf("no files matched\n"); return; }

    Int_t   bRun = 0, bEvt = 0, bTrk = 0, bDet = 0;
    Float_t bHx = 0, bHy = 0, bHz = 0;
    tr->SetBranchAddress("run",     &bRun);
    tr->SetBranchAddress("event",   &bEvt);
    tr->SetBranchAddress("trackId", &bTrk);
    tr->SetBranchAddress("detType", &bDet);
    tr->SetBranchAddress("hitX", &bHx);
    tr->SetBranchAddress("hitY", &bHy);
    tr->SetBranchAddress("hitZ", &bHz);

    double diskZ[3]; diskZ[0] = 151.750; diskZ[1] = 165.248; diskZ[2] = 178.781;

    // one strip per bin over +-72 strips, same convention as the other plots
    const int    kNB   = 144;
    const double kHalf = 72.0*30.0/128.0;
    TH1F* hWith = new TH1F("hWith", ";#delta#phi_{local} of FST hit [deg];FST hits",
                           kNB, -kHalf, kHalf);
    TH1F* hNo   = new TH1F("hNo",   ";#delta#phi_{local} of FST hit [deg];FST hits",
                           kNB, -kHalf, kHalf);
    hWith->SetDirectory(0); hNo->SetDirectory(0);

    // per-track buffer
    const int kMaxFst = 64;
    double bufDphi[64];
    int    nBufFst = 0, nFtt = 0;
    int    curRun = -1, curEvt = -1, curTrk = -1;
    Long64_t nTrkWith = 0, nTrkNo = 0;

    Long64_t ne = tr->GetEntries();
    printf("rows: %lld\n", ne);

    for (Long64_t ie = 0; ie <= ne; ie++){
        int isNew = 0;
        if (ie < ne){
            tr->GetEntry(ie);
            if (bRun != curRun || bEvt != curEvt || bTrk != curTrk) isNew = 1;
        } else {
            isNew = 1;                    // final flush
        }

        if (isNew && curTrk >= 0){
            // flush the previous track
            for (int k = 0; k < nBufFst; k++){
                if (nFtt > 0) hWith->Fill(bufDphi[k]);
                else          hNo  ->Fill(bufDphi[k]);
            }
            if (nBufFst > 0){ if (nFtt > 0) nTrkWith++; else nTrkNo++; }
            nBufFst = 0; nFtt = 0;
        }
        if (ie >= ne) break;

        if (isNew){ curRun = bRun; curEvt = bEvt; curTrk = bTrk; }

        if (bDet == 1){ nFtt++; continue; }
        if (bDet != 0) continue;

        int ok = 0;
        for (int db = 0; db < 3; db++) if (fabs(bHz - diskZ[db]) < 5.0) ok = 1;
        if (!ok) continue;

        double p360 = atan2(bHy, bHx)*180.0/TMath::Pi();
        if (p360 < 0) p360 += 360.0;
        int sec = (int)(p360/30.0); if (sec > 11) sec = 11;
        double dLoc = p360 - (sec*30.0 + 15.0);
        if (nBufFst < kMaxFst) bufDphi[nBufFst++] = dLoc;
    }

    printf("\ntracks with >=1 FTT hit: %lld   tracks with none: %lld\n", nTrkWith, nTrkNo);
    printf("FST hits on those tracks: %.0f / %.0f\n", hWith->Integral(), hNo->Integral());

    // the discriminator: is the WITH-FTT sample piled up at the wedge centre?
    printf("\n fraction of FST hits within |dphi_local| of ...\n");
    printf("   window     with FTT    without FTT   (uniform expectation)\n");
    double totW = hWith->Integral(), totN = hNo->Integral();
    double win[4]; win[0]=2.0; win[1]=5.0; win[2]=10.0; win[3]=15.0;
    for (int iw = 0; iw < 4; iw++){
        double w = win[iw];
        double fw = hWith->Integral(hWith->FindBin(-w+1e-6), hWith->FindBin(w-1e-6));
        double fn = hNo  ->Integral(hNo  ->FindBin(-w+1e-6), hNo  ->FindBin(w-1e-6));
        printf("   %4.1f deg    %6.2f%%      %6.2f%%          %6.2f%%\n",
               w, totW>0?100.0*fw/totW:0, totN>0?100.0*fn/totN:0, 100.0*w/15.0);
    }
    printf("\n mean |dphi_local| :  with FTT %.3f deg,  without %.3f deg\n",
           fabs(hWith->GetMean()) > 0 ? hWith->GetRMS() : 0, hNo->GetRMS());

    TCanvas* c = new TCanvas("c_mir", "", 700, 480);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    if (totW > 0) hWith->Scale(1.0/totW);
    if (totN > 0) hNo  ->Scale(1.0/totN);
    hWith->SetLineColor(kBlue+1); hWith->SetLineWidth(2); hWith->SetFillStyle(0);
    hNo  ->SetLineColor(kRed+1);  hNo  ->SetLineWidth(2); hNo  ->SetFillStyle(0);
    hWith->SetTitle("FST hit position within its wedge");
    hWith->GetYaxis()->SetTitle("fraction of FST hits / bin");
    hWith->SetMinimum(0);
    hWith->Draw("hist");
    hNo->Draw("hist same");
    TLegend* lg = new TLegend(0.30,0.18,0.72,0.33);
    lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(hWith, "track HAS FTT hits", "l");
    lg->AddEntry(hNo,   "track has NO FTT hits", "l");
    lg->Draw();
    c->SaveAs(Form("%s/fstMirrorFttTest.png", outdir));
    printf("\nwrote %s/fstMirrorFttTest.png\n", outdir);
}
