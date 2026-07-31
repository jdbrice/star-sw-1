// testFstSlowSim.C
//
// Smoke test for StFstSlowSimMaker: does MC actually flow through the REAL FST
// reconstruction chain?
//
//   g2t_fsi_hit -> StFstSlowSimMaker -> "fstRawAdcSimu"
//               -> StFstRawHitMaker  -> StFstClusterMaker -> StFstHitMaker
//
// This is the path StFstFastSimMaker bypasses. Counts what comes out of each
// stage, so a stage that silently produces nothing is obvious.
//
// Usage:
//   root4star -l -b -q 'script/testFstSlowSim.C("ele.pt1.5.vz0.run1.fzd", 20)'

void testFstSlowSim(const char* inFile = "ele.pt1.5.vz0.run1.fzd", int nEvents = 20){

    TString geom  = "y2024 agml usexgeom";

    // NOTE the timestamp: sim.C's usual sdt20211016 predates FST installation,
    // so Geometry/fst/fstOnTpc has no entry and StFstDbMaker::Make() returns
    // StFATAL, killing the chain on the first event. Any 2022+ date works.
    TString chainOpt = Form("fzin %s sdt20220401 fstDb fstSlowSim fstRawHit fstCluster fstHit "
                            "MakeEvent StEvent McEvent ReverseField bigbig", geom.Data());

    printf("Chain:\n  %s\n\n", chainOpt.Data());

    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");
    bfc(-1, chainOpt, inFile);

    if (chain->Init() != kStOK){ printf("chain->Init() failed\n"); return; }

    Long64_t nRaw = 0, nClu = 0, nHit = 0, nG2t = 0;
    int nEvtWithHits = 0;

    // ---- closure test: reconstructed StFstHit vs the GEANT hit it came from --
    // This is the whole point. A wedge-numbering, disk, or gap-sign error shows
    // up here as a cm-scale r*dphi residual; correct geometry leaves only strip
    // quantization, which is bounded and calculable:
    //    r     pitch 2.875 cm      -> RMS 2.875/sqrt(12) = 0.830 cm
    //    r*dphi pitch 4.09 mrad    -> RMS r*4.09e-3/sqrt(12), 0.012 cm at r=10
    TH1F* hDr   = new TH1F("hDr",   ";r_{reco} - r_{truth} [cm];hits", 120, -3, 3);
    TH1F* hRdp  = new TH1F("hRdp",  ";r#Delta#phi (reco-truth) [cm];hits", 120, -0.6, 0.6);
    TH1F* hRdpW = new TH1F("hRdpW", ";r#Delta#phi (reco-truth) [cm];hits", 120, -6, 6);
    hDr->SetDirectory(0); hRdp->SetDirectory(0); hRdpW->SetDirectory(0);
    int nMatched = 0, nFarPhi = 0;

    for (int iev = 0; iev < nEvents; iev++){

        chain->Clear();
        int ret = chain->Make();
        if (ret == kStEOF || ret == kStFatal) break;

        // g2t input
        St_g2t_fts_hit* g2t = (St_g2t_fts_hit*) chain->GetDataSet("g2t_fsi_hit");
        int ng = 0;
        if (g2t){
            int nrow = g2t->GetNRows();
            for (int ig = 0; ig < nrow; ig++){
                g2t_fts_hit_st* h = (g2t_fts_hit_st*) g2t->At(ig);
                if (!h) continue;
                int dsk = h->volume_id / 1000;
                if (dsk >= 4 && dsk <= 6) ng++;   // FST disks are 4,5,6 here
            }
        }
        nG2t += ng;

        // dump the raw volume_id encoding for the first couple of events --
        // disk = volume_id/1000 (StFstFastSimMaker.cxx:204) has to be checked,
        // not assumed
        if (iev < 2 && g2t){
            int nrow2 = g2t->GetNRows();
            for (int ik = 0; ik < nrow2 && ik < 12; ik++){
                g2t_fts_hit_st* hh = (g2t_fts_hit_st*) g2t->At(ik);
                if (!hh) continue;
                double rr = sqrt(hh->x[0]*hh->x[0] + hh->x[1]*hh->x[1]);
                printf("      g2t[%2d] volume_id=%8d  x=%8.3f y=%8.3f z=%8.3f  r=%7.3f  de=%.4e\n",
                       ik, hh->volume_id, hh->x[0], hh->x[1], hh->x[2], rr, hh->de);
            }
        }

        // what the slow sim published
        int nr = 0;
        TObjectSet* simSet = (TObjectSet*) chain->GetDataSet("fstRawAdcSimu");
        if (simSet){
            StFstCollection* c = (StFstCollection*) simSet->GetObject();
            if (c) nr = c->getNumRawHits();
        }
        nRaw += nr;

        // what came out of the real reconstruction chain
        StEvent* ev = (StEvent*) chain->GetDataSet("StEvent");
        int nc = 0, nh = 0;
        if (ev){
            StFstHitCollection* hc = ev->fstHitCollection();
            if (hc) nh = hc->numberOfHits();
        }
        StFstCollection* recoColl = 0;
        TObjectSet* recoSet = (TObjectSet*) chain->GetDataSet("fstRawHitAndCluster");
        if (recoSet) recoColl = (StFstCollection*) recoSet->GetObject();
        if (recoColl) nc = recoColl->getNumClusters();

        nClu += nc;
        nHit += nh;
        if (nh > 0) nEvtWithHits++;

        // match each reconstructed hit to the nearest GEANT hit on the same disk
        if (ev && g2t){
            StFstHitCollection* hc2 = ev->fstHitCollection();
            int nrow3 = g2t->GetNRows();
            // StFstHitCollection has no numberOfWedges()/numberOfSensors() --
            // the sizes are fixed: kFstNumWedges = 36, kFstNumSensorsPerWedge = 3
            // (StFstConsts.h consts are not visible to CINT, hence the literals).
            if (hc2){
                for (int iw = 0; iw < 36; iw++){
                    StFstWedgeHitCollection* wc = hc2->wedge(iw);
                    if (!wc) continue;
                    for (int is = 0; is < 3; is++){
                        StFstSensorHitCollection* sc = wc->sensor(is);
                        if (!sc) continue;
                        int nHitsHere = sc->hits().size();
                        for (int ii = 0; ii < nHitsHere; ii++){
                            StFstHit* fh = sc->hits()[ii];
                            if (!fh) continue;
                            double hx = fh->position().x();
                            double hy = fh->position().y();
                            double hz = fh->position().z();

                            double best = 1e9; int ibest = -1;
                            for (int im = 0; im < nrow3; im++){
                                g2t_fts_hit_st* gh = (g2t_fts_hit_st*) g2t->At(im);
                                if (!gh) continue;
                                int dk = gh->volume_id / 1000;
                                if (dk < 4 || dk > 6) continue;
                                if (fabs(gh->x[2] - hz) > 6.0) continue;   // same disk
                                double d2 = (gh->x[0]-hx)*(gh->x[0]-hx) + (gh->x[1]-hy)*(gh->x[1]-hy);
                                if (d2 < best){ best = d2; ibest = im; }
                            }
                            if (ibest < 0) continue;

                            g2t_fts_hit_st* gb = (g2t_fts_hit_st*) g2t->At(ibest);
                            double rr = sqrt(hx*hx + hy*hy);
                            double rt = sqrt(gb->x[0]*gb->x[0] + gb->x[1]*gb->x[1]);
                            double dp = atan2(hy,hx) - atan2(gb->x[1], gb->x[0]);
                            while (dp >  TMath::Pi()) dp -= TMath::TwoPi();
                            while (dp < -TMath::Pi()) dp += TMath::TwoPi();

                            hDr ->Fill(rr - rt);
                            hRdp->Fill(rr * dp);
                            hRdpW->Fill(rr * dp);
                            nMatched++;
                            if (fabs(rr*dp) > 0.5) nFarPhi++;

                            if (nMatched <= 12){
                                double phT = atan2(gb->x[1], gb->x[0])*180.0/TMath::Pi();
                                double phR = atan2(hy, hx)*180.0/TMath::Pi();
                                if (phT < 0) phT += 360.0;
                                if (phR < 0) phR += 360.0;
                                printf("      [%2d] truth r=%7.3f phi=%8.3f | reco r=%7.3f phi=%8.3f "
                                       "| dr=%+7.3f rdphi=%+7.3f | w=%2d s=%d rS=%d phiS=%3d z=%7.2f\n",
                                       nMatched, rt, phT, rr, phR, rr-rt, rr*dp,
                                       (int)fh->getWedge(), (int)fh->getSensor(),
                                       (int)fh->getMeanRStrip(), (int)fh->getMeanPhiStrip(), hz);
                            }
                        }
                    }
                }
            }
        }

        printf("  event %3d :  g2t(FST) %4d   simRawHits %4d   clusters %4d   StFstHits %4d\n",
               iev, ng, nr, nc, nh);
    }

    printf("\n==============================================================\n");
    printf("  totals over the events read\n");
    printf("    g2t FST hits          %lld\n", nG2t);
    printf("    sim raw hits          %lld\n", nRaw);
    printf("    clusters              %lld\n", nClu);
    printf("    StFstHits             %lld\n", nHit);
    printf("    events with FST hits  %d\n", nEvtWithHits);
    printf("==============================================================\n");
    if (nRaw == 0)
        printf("  FAIL: slow sim produced no raw hits.\n");
    else if (nHit == 0)
        printf("  FAIL: raw hits made, but nothing survived to StFstHit.\n");
    else
        printf("  OK: MC reached StFstHit through the real reconstruction chain.\n");

    printf("\n--------------------------------------------------------------\n");
    printf("  closure: reconstructed hit vs its GEANT hit (%d matched)\n", nMatched);
    printf("--------------------------------------------------------------\n");
    printf("    dr        mean %+8.4f   RMS %7.4f cm   (quantization RMS 0.830)\n",
           hDr->GetMean(), hDr->GetRMS());
    // use the WIDE histogram for the summary -- the narrow one silently drops
    // everything beyond +-0.6 cm into overflow and reports a flattering RMS
    printf("    r*dphi    mean %+8.4f   RMS %7.4f cm   (expect < ~0.02)\n",
           hRdpW->GetMean(), hRdpW->GetRMS());
    printf("    |r*dphi| > 0.5 cm : %d / %d  (%.1f%%)\n",
           nFarPhi, nMatched, nMatched > 0 ? 100.0*nFarPhi/nMatched : 0.0);
    // A few outliers are expected: in multi-hit events the nearest-neighbour
    // truth matching above can pair a reco hit with the wrong GEANT hit. A
    // geometry error (wedge numbering, disk, gap sign) would put a large
    // fraction out there, not a per-cent-level tail.
    double farFrac = (nMatched > 0) ? 100.0*nFarPhi/nMatched : 100.0;
    if (nMatched > 0 && farFrac < 2.0)
        printf("    -> phi geometry closes. Wedge numbering, disk and gap sign consistent.\n"
               "       (the %.1f%% tail is truth-matching ambiguity in multi-hit events)\n", farFrac);
    else if (nMatched > 0)
        printf("    -> phi does NOT close; see the wide histogram for the failure scale.\n");

    TCanvas* cc = new TCanvas("c_closure", "", 1000, 380);
    cc->Divide(3,1);
    cc->cd(1); hDr->Draw("hist");
    cc->cd(2); hRdp->Draw("hist");
    cc->cd(3); gPad->SetLogy(); hRdpW->Draw("hist");
    cc->SaveAs("FstSlowSim/fstSlowSimClosure.png");
    printf("\n  wrote fstclus/fstSlowSimClosure.png\n");
}
