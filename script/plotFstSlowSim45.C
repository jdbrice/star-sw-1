// plotFstSlowSim45.C
//
// The 45-degree single-electron sanity check, run through the NEW FST slow-sim
// chain (see FstSlowSim/index.html):
//
//   g2t_fsi_hit -> StFstSlowSimMaker -> "fstRawAdcSimu"
//               -> StFstRawHitMaker -> StFstClusterMaker -> StFstHitMaker
//
// Same idea as jpsi/electron_45.html did for the FTT: a single isolated track
// with known kinematics (e-, pT 1.5 GeV, eta 2.5-4.0, phi = 45 +- 0.06 deg at
// the vertex) is the cleanest possible geometry probe -- no multi-track overlap,
// so every reconstructed hit can be attributed unambiguously.
//
// NOTE the generated phi is the phi AT THE VERTEX. The electron is charged and
// bends over the 150-180 cm to the FST, so the arrival phi has real tails and
// hits land in several wedges. That is physics, not a defect: the GEANT truth
// shows exactly the same spread. A "must all be in one wedge" test would be
// wrong.
//
// What this adds over the generic closure test in testFstSlowSim.C: that test
// measures residuals in aggregate, which is blind to a whole-wedge mislabel that
// happens to preserve r and phi. Here the wedge index is checked discretely and
// per hit, against the wedge implied by that same hit's truth phi.
//
// Usage:
//   root4star -l -b -q 'script/plotFstSlowSim45.C("ele45.vz0.run1.fzd", 500)'

void plotFstSlowSim45(const char* inFile = "ele45.vz0.run1.fzd", int nEvents = 500,
                      const char* outdir = "FstSlowSim/electron45"){

    gSystem->mkdir(outdir, kTRUE);
    gStyle->SetOptStat(0);

    TString geom = "y2024 agml usexgeom";
    // sdt20211016 (sim.C's default) predates FST installation -> StFstDbMaker
    // returns StFATAL and the chain dies on event 0. Any 2022+ date works.
    TString chainOpt = Form("fzin %s sdt20220401 fstDb fstSlowSim fstRawHit fstCluster fstHit "
                            "MakeEvent StEvent McEvent ReverseField bigbig", geom.Data());
    printf("Chain:\n  %s\n\n", chainOpt.Data());

    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");
    bfc(-1, chainOpt, inFile);
    if (chain->Init() != kStOK){ printf("chain->Init() failed\n"); return; }

    // ---- histograms -------------------------------------------------------
    TH2F* hXYr[3];   // reconstructed StFstHit
    TH2F* hXYt[3];   // GEANT truth
    TH1F* hWedge[3];
    for (int ia = 0; ia < 3; ia++){
        hXYr[ia] = new TH2F(Form("hXYr_%d", ia), Form("reco, disk %d;x [cm];y [cm]", ia+1),
                            160, -30, 30, 160, -30, 30);
        hXYt[ia] = new TH2F(Form("hXYt_%d", ia), Form("GEANT truth, disk %d;x [cm];y [cm]", ia+1),
                            160, -30, 30, 160, -30, 30);
        hWedge[ia] = new TH1F(Form("hWedge_%d", ia), Form("disk %d;wedge;hits", ia+1),
                              37, -0.5, 36.5);
        hXYr[ia]->SetDirectory(0); hXYt[ia]->SetDirectory(0); hWedge[ia]->SetDirectory(0);
    }
    TH1F* hPhiR = new TH1F("hPhiR", ";#phi [deg];hits", 180, 0, 90);
    TH1F* hPhiT = new TH1F("hPhiT", ";#phi [deg];hits", 180, 0, 90);
    TH1F* hRr   = new TH1F("hRr",   ";r [cm];hits", 120, 0, 30);
    TH1F* hRt   = new TH1F("hRt",   ";r [cm];hits", 120, 0, 30);
    TH1F* hDr   = new TH1F("hDr",   ";r_{reco} - r_{truth} [cm];hits", 120, -3, 3);
    TH1F* hRdp  = new TH1F("hRdp",  ";r#Delta#phi (reco-truth) [cm];hits", 120, -0.3, 0.3);
    TH1F* hSens = new TH1F("hSens", ";sensor;hits", 4, -0.5, 3.5);
    hPhiR->SetDirectory(0); hPhiT->SetDirectory(0); hRr->SetDirectory(0); hRt->SetDirectory(0);
    hDr->SetDirectory(0); hRdp->SetDirectory(0); hSens->SetDirectory(0);

    Long64_t nG2t = 0, nHit = 0, nMatched = 0, nWrongWedge = 0;
    // phi = 45 deg -> sector 1 -> moduleIdx 2 -> wedge 2/14/26 on disk 1/2/3.
    // This is only the DOMINANT wedge, not the only one: the electron is charged
    // and bends in phi over the 150-180 cm to the FST, so the arrival phi has
    // real tails (visible identically in the GEANT truth). The meaningful check
    // is therefore reco-wedge vs the wedge implied by the TRUTH phi of the same
    // hit, done per hit below -- not reco vs a fixed number.
    int expectWedge[3]; expectWedge[0] = 2; expectWedge[1] = 14; expectWedge[2] = 26;

    // sector (phi/30deg) -> moduleIdx, from kFstphiStart/kFstphiStop.
    // (StFstConsts.h consts are not visible to CINT, so the map is inlined.)
    int sectorToModule[12];
    sectorToModule[0]=3;  sectorToModule[1]=2;  sectorToModule[2]=1;  sectorToModule[3]=12;
    sectorToModule[4]=11; sectorToModule[5]=10; sectorToModule[6]=9;  sectorToModule[7]=8;
    sectorToModule[8]=7;  sectorToModule[9]=6;  sectorToModule[10]=5; sectorToModule[11]=4;

    for (int iev = 0; iev < nEvents; iev++){
        chain->Clear();
        int ret = chain->Make();
        if (ret == kStEOF || ret == kStFatal) break;

        St_g2t_fts_hit* g2t = (St_g2t_fts_hit*) chain->GetDataSet("g2t_fsi_hit");
        int nrow = g2t ? g2t->GetNRows() : 0;

        for (int ig = 0; ig < nrow; ig++){
            g2t_fts_hit_st* gh = (g2t_fts_hit_st*) g2t->At(ig);
            if (!gh) continue;
            int dk = gh->volume_id / 1000;          // FST disks are 4,5,6
            if (dk < 4 || dk > 6) continue;
            int d = dk - 4;
            double tx = gh->x[0], ty = gh->x[1];
            hXYt[d]->Fill(tx, ty);
            double pt = atan2(ty, tx)*180.0/TMath::Pi();
            if (pt < 0) pt += 360.0;
            hPhiT->Fill(pt);
            hRt->Fill(sqrt(tx*tx + ty*ty));
            nG2t++;
        }

        StEvent* ev = (StEvent*) chain->GetDataSet("StEvent");
        if (!ev) continue;
        StFstHitCollection* hc = ev->fstHitCollection();
        if (!hc) continue;

        // fixed sizes: kFstNumWedges = 36, kFstNumSensorsPerWedge = 3
        for (int iw = 0; iw < 36; iw++){
            StFstWedgeHitCollection* wc = hc->wedge(iw);
            if (!wc) continue;
            for (int is = 0; is < 3; is++){
                StFstSensorHitCollection* sc = wc->sensor(is);
                if (!sc) continue;
                int nh = sc->hits().size();
                for (int ii = 0; ii < nh; ii++){
                    StFstHit* fh = sc->hits()[ii];
                    if (!fh) continue;
                    int d = (int) fh->getDisk() - 1;
                    if (d < 0 || d > 2) continue;

                    double hx = fh->position().x();
                    double hy = fh->position().y();
                    double hz = fh->position().z();
                    double rr = sqrt(hx*hx + hy*hy);
                    double pr = atan2(hy, hx)*180.0/TMath::Pi();
                    if (pr < 0) pr += 360.0;

                    hXYr[d]->Fill(hx, hy);
                    hPhiR->Fill(pr);
                    hRr->Fill(rr);
                    hWedge[d]->Fill((int) fh->getWedge());
                    hSens->Fill((int) fh->getSensor());
                    nHit++;

                    // nearest GEANT hit on the same disk
                    double best = 1e9; int ibest = -1;
                    for (int im = 0; im < nrow; im++){
                        g2t_fts_hit_st* gg = (g2t_fts_hit_st*) g2t->At(im);
                        if (!gg) continue;
                        int dk2 = gg->volume_id / 1000;
                        if (dk2 < 4 || dk2 > 6) continue;
                        if (fabs(gg->x[2] - hz) > 6.0) continue;
                        double d2 = (gg->x[0]-hx)*(gg->x[0]-hx) + (gg->x[1]-hy)*(gg->x[1]-hy);
                        if (d2 < best){ best = d2; ibest = im; }
                    }
                    if (ibest < 0) continue;
                    g2t_fts_hit_st* gb = (g2t_fts_hit_st*) g2t->At(ibest);
                    double rt = sqrt(gb->x[0]*gb->x[0] + gb->x[1]*gb->x[1]);
                    double dp = atan2(hy,hx) - atan2(gb->x[1], gb->x[0]);
                    while (dp >  TMath::Pi()) dp -= TMath::TwoPi();
                    while (dp < -TMath::Pi()) dp += TMath::TwoPi();
                    hDr->Fill(rr - rt);
                    hRdp->Fill(rr * dp);
                    nMatched++;

                    // the discrete check: does the reco wedge match the wedge
                    // implied by this hit's own TRUTH phi?
                    double ptr = atan2(gb->x[1], gb->x[0])*180.0/TMath::Pi();
                    if (ptr < 0) ptr += 360.0;
                    int sec = (int)(ptr/30.0); if (sec < 0) sec = 0; if (sec > 11) sec = 11;
                    int wTruth = sectorToModule[sec] + d*12;
                    if ((int) fh->getWedge() != wTruth) nWrongWedge++;
                }
            }
        }
    }

    // ---- numbers ----------------------------------------------------------
    printf("\n==============================================================\n");
    printf(" 45-degree electron through the FST slow-sim chain\n");
    printf("==============================================================\n");
    printf("  GEANT FST hits   %lld\n", nG2t);
    printf("  StFstHits        %lld\n", nHit);
    printf("  matched          %lld\n", nMatched);
    printf("\n  phi   truth mean %.3f deg   reco mean %.3f deg\n",
           hPhiT->GetMean(), hPhiR->GetMean());
    printf("  r     truth mean %.3f cm    reco mean %.3f cm\n",
           hRt->GetMean(), hRr->GetMean());
    printf("\n  wedge occupancy (dominant wedge is the prediction;\n"
           "   the tails are real phi bending, present in truth too):\n");
    for (int id = 0; id < 3; id++){
        int nOcc = 0, wBest = -1; double cBest = 0;
        for (int ib = 1; ib <= hWedge[id]->GetNbinsX(); ib++){
            double c = hWedge[id]->GetBinContent(ib);
            if (c <= 0) continue;
            nOcc++;
            if (c > cBest){ cBest = c; wBest = (int)hWedge[id]->GetBinCenter(ib); }
        }
        printf("    disk %d : dominant = %d (expect %d) %s ; %d wedges occupied in total\n",
               id+1, wBest, expectWedge[id],
               (wBest == expectWedge[id]) ? "OK" : "<-- CHECK", nOcc);
    }
    printf("\n  reco wedge != wedge implied by the hit's own truth phi :\n");
    printf("    %lld / %lld  (%.2f%%)   <- this is the real wedge-numbering test\n",
           nWrongWedge, nMatched, nMatched > 0 ? 100.0*nWrongWedge/nMatched : 0.0);
    printf("  sensor split: inner %.0f  outer1 %.0f  outer2 %.0f\n",
           hSens->GetBinContent(1), hSens->GetBinContent(2), hSens->GetBinContent(3));
    printf("\n  dr      mean %+.4f  RMS %.4f cm   (quantization RMS 0.830)\n",
           hDr->GetMean(), hDr->GetRMS());
    printf("  r*dphi  mean %+.4f  RMS %.4f cm\n", hRdp->GetMean(), hRdp->GetRMS());

    // ---- plots ------------------------------------------------------------
    TCanvas* c1 = new TCanvas("c45xy", "", 1100, 720);
    c1->Divide(3,2);
    for (int ip = 0; ip < 3; ip++){
        c1->cd(ip+1); gPad->SetRightMargin(0.13);
        hXYt[ip]->Draw("colz");
        c1->cd(ip+4); gPad->SetRightMargin(0.13);
        hXYr[ip]->Draw("colz");
    }
    c1->SaveAs(Form("%s/fst45_xy.png", outdir));

    TCanvas* c2 = new TCanvas("c45qa", "", 1100, 720);
    c2->Divide(3,2);
    c2->cd(1); gPad->SetLogy();
    hPhiT->SetLineColor(kBlack); hPhiT->SetLineWidth(2); hPhiT->SetFillStyle(0);
    hPhiR->SetLineColor(kRed+1); hPhiR->SetLineWidth(2); hPhiR->SetFillStyle(0);
    hPhiT->SetTitle("#phi: truth (black) vs reco (red)");
    hPhiT->Draw("hist"); hPhiR->Draw("hist same");
    c2->cd(2);
    hRt->SetLineColor(kBlack); hRt->SetLineWidth(2); hRt->SetFillStyle(0);
    hRr->SetLineColor(kRed+1); hRr->SetLineWidth(2); hRr->SetFillStyle(0);
    hRt->SetTitle("r: truth (black) vs reco (red)");
    hRt->Draw("hist"); hRr->Draw("hist same");
    c2->cd(3);
    hSens->SetLineColor(kBlue+1); hSens->SetLineWidth(2); hSens->SetFillStyle(0);
    hSens->SetTitle("sensor (0 = inner)");
    hSens->Draw("hist");
    c2->cd(4);
    hWedge[0]->SetLineColor(kBlue+1); hWedge[0]->SetLineWidth(2); hWedge[0]->SetFillStyle(0);
    hWedge[1]->SetLineColor(kRed+1);  hWedge[1]->SetLineWidth(2); hWedge[1]->SetFillStyle(0);
    hWedge[2]->SetLineColor(kGreen+2);hWedge[2]->SetLineWidth(2); hWedge[2]->SetFillStyle(0);
    hWedge[0]->SetTitle("wedge, disk 1 / 2 / 3");
    hWedge[0]->Draw("hist"); hWedge[1]->Draw("hist same"); hWedge[2]->Draw("hist same");
    c2->cd(5);
    hDr->SetLineColor(kBlue+1); hDr->SetLineWidth(2); hDr->SetFillStyle(0);
    hDr->SetTitle("r_{reco} - r_{truth}");
    hDr->Draw("hist");
    c2->cd(6);
    hRdp->SetLineColor(kBlue+1); hRdp->SetLineWidth(2); hRdp->SetFillStyle(0);
    hRdp->SetTitle("r#Delta#phi (reco - truth)");
    hRdp->Draw("hist");
    c2->SaveAs(Form("%s/fst45_qa.png", outdir));

    printf("\nwrote %s/fst45_xy.png, %s/fst45_qa.png\n", outdir, outdir);
}
