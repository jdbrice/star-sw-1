// checkFstInnerOuter.C
//
// Follow-up to script/checkFstWedgeGeom.C (Test A of proposal_residual_mc.txt).
//
// That geometry check compared StFstHitMaker's strip positions against the
// AGML silicon placement for all 108 sensors and found:
//
//   - INNER sensors: exact, 0.0000 deg on all 36 wedges.
//   - OUTER sensors: apparently 1.000 deg off -- but that depended on an
//     assumption about which of the two outer sensors carries wedge-level
//     phiStrip 0-63 vs 64-127. That assignment is NOT derivable from code
//     (StFstRawHit::getSensor() comes from mChannelId, getPhiStrip() from
//     mGeoId; only StFstDb's mapping table links them). With the opposite
//     assignment the reconstruction is exact. This macro is what settles it.
//   - z: StFstHitMaker assigns ONE z per disk, the mean of the front/back
//     wedge z, while the silicon alternates +-1.742 cm (inner) / +-0.347 cm
//     (outer) wedge-to-wedge, with the parity FLIPPED on disk 2. That
//     mismatch is real and is 30-deg periodic -- the surviving candidate.
//
// So this macro tests two things at once:
//   (a) outer-sensor phi: a 1.000 deg error would show as a ~0.39 cm
//       alternating structure in mean r*dphi for r>16.5 cm only. Absence
//       of it means the phiStrip/sensor assignment is the other way and
//       there is no outer-sensor bug.
//   (b) the z stagger: predicts an alternating RADIAL residual of
//       ~0.17 cm (inner) and ~0.05 cm (outer), sign flipping between
//       disk 2 and disks 1,3.
//
// FST radii are quantised by strip (kFstrStart + half pitch):
//   inner rStrip 0-3 -> r = 6.44, 9.31, 12.19, 15.06 cm
//   outer rStrip 4-7 -> r = 17.94, 20.81, 23.69, 26.56 cm
// so r = 16.5 cm is a clean, unambiguous inner/outer split.
//
// Usage:
//   root4star -b -q 'script/checkFstInnerOuter.C'
//   root4star -b -q 'script/checkFstInnerOuter.C("<glob>","<outdir>",nFilesMax)'

void checkFstInnerOuter(
        const char* fname  = "/gpfs01/star/pwg_tasks/FwdCalib/akio/pico_20260704/alignment/*.FwdAlignment_BLCVtx.root",
        const char* outdir = "residual/data/alignment/innerouter",
        int nFilesMax = 0 /* 0 = all */ ){

    TChain* tr = new TChain("alignTree");
    int nAdded = tr->Add(fname);
    if (nAdded <= 0){ printf("No files matched %s\n", fname); return; }
    printf("Matched %d files\n", nAdded);
    gSystem->mkdir(outdir, kTRUE);
    gStyle->SetOptStat(0);

    Int_t   bDet = 0;
    Float_t bHx = 0, bHy = 0, bHz = 0, bPx = 0, bPy = 0;
    tr->SetBranchAddress("detType", &bDet);
    tr->SetBranchAddress("hitX",  &bHx);
    tr->SetBranchAddress("hitY",  &bHy);
    tr->SetBranchAddress("hitZ",  &bHz);
    tr->SetBranchAddress("projX", &bPx);
    tr->SetBranchAddress("projY", &bPy);

    // FST disk z (StFstHitMaker.cxx:167-169)
    double diskZ[3]; diskZ[0] = 151.750; diskZ[1] = 165.248; diskZ[2] = 178.781;

    // profiles of mean r*dphi vs phi, per disk, split inner/outer.
    // 120 bins over 360 deg = 3 deg bins, 10 per wedge.
    const int kNPhiBin = 120;
    TProfile* pIn[3];
    TProfile* pOut[3];
    // radial residual: the z-stagger hypothesis predicts a MUCH bigger signal
    // here than in r*dphi. StFstHitMaker assigns one z per disk (the mean of
    // the front/back wedge z), but the silicon alternates by +-1.742 cm
    // (inner) / +-0.347 cm (outer) wedge-to-wedge. A track projected to the
    // assumed z instead of the true one is displaced radially by
    // dr ~ r*dz/z -- about 0.17 cm at r=15 cm, 0.05 cm at r=22 cm.
    TProfile* qIn[3];
    TProfile* qOut[3];
    for (int da = 0; da < 3; da++){
        pIn [da] = new TProfile(Form("pIn_d%d",  da), Form("FST%d inner (r<16.5cm);#phi [deg];<r#Delta#phi> [cm]", da),
                                kNPhiBin, -180, 180, -3, 3);
        pOut[da] = new TProfile(Form("pOut_d%d", da), Form("FST%d outer (r>16.5cm);#phi [deg];<r#Delta#phi> [cm]", da),
                                kNPhiBin, -180, 180, -3, 3);
        pIn [da]->SetDirectory(0);
        pOut[da]->SetDirectory(0);
        qIn [da] = new TProfile(Form("qIn_d%d",  da), Form("FST%d inner dr;#phi [deg];<#Deltar> [cm]", da),
                                kNPhiBin, -180, 180, -5, 5);
        qOut[da] = new TProfile(Form("qOut_d%d", da), Form("FST%d outer dr;#phi [deg];<#Deltar> [cm]", da),
                                kNPhiBin, -180, 180, -5, 5);
        qIn [da]->SetDirectory(0);
        qOut[da]->SetDirectory(0);
    }

    Long64_t nEnt = tr->GetEntries();
    printf("Total rows: %lld\n", nEnt);

    Long64_t nUsedIn = 0, nUsedOut = 0;
    for (Long64_t ie = 0; ie < nEnt; ie++){
        tr->GetEntry(ie);
        if (bDet != 0) continue;                       // FST only

        int dsel = -1;
        for (int db = 0; db < 3; db++)
            if (fabs(bHz - diskZ[db]) < 5.0) dsel = db;
        if (dsel < 0) continue;

        double rr   = sqrt(bHx*bHx + bHy*bHy);
        double phih = atan2(bHy, bHx);
        double phip = atan2(bPy, bPx);
        double dph  = phih - phip;
        while (dph >  TMath::Pi()) dph -= 2*TMath::Pi();
        while (dph < -TMath::Pi()) dph += 2*TMath::Pi();
        double rdp  = rr * dph;

        double rproj = sqrt(bPx*bPx + bPy*bPy);
        double drad  = rr - rproj;
        double phiDeg = phih * 180.0 / TMath::Pi();
        if (rr < 16.5){ pIn [dsel]->Fill(phiDeg, rdp); qIn [dsel]->Fill(phiDeg, drad); nUsedIn++;  }
        else          { pOut[dsel]->Fill(phiDeg, rdp); qOut[dsel]->Fill(phiDeg, drad); nUsedOut++; }

        if ((ie % 20000000) == 0 && ie > 0) printf("  ... %lld / %lld\n", ie, nEnt);
    }
    printf("Filled: inner %lld, outer %lld\n", nUsedIn, nUsedOut);

    // ---- quantify the 30-deg alternation ----
    // Average each profile into 12 wedge bins (30 deg each, edges at multiples
    // of 30 deg -- matching kFstphiStart, whose sectors are exactly [0,30),
    // [30,60), ...), then report the alternating amplitude:
    //   A = 0.5 * ( <even-index wedges> - <odd-index wedges> )
    // which is what a sign that flips every wedge produces.
    printf("\n===============================================================\n");
    printf(" Mean r*dphi per 30-deg wedge sector [cm]\n");
    printf("===============================================================\n");
    for (int dc = 0; dc < 3; dc++){
        double wIn[12], wOut[12];
        for (int wa = 0; wa < 12; wa++){ wIn[wa] = 0; wOut[wa] = 0; }
        double cIn[12], cOut[12];
        for (int wb = 0; wb < 12; wb++){ cIn[wb] = 0; cOut[wb] = 0; }
        double yIn[12], yOut[12], kIn[12], kOut[12];
        for (int wd = 0; wd < 12; wd++){ yIn[wd]=0; yOut[wd]=0; kIn[wd]=0; kOut[wd]=0; }

        for (int ib = 1; ib <= kNPhiBin; ib++){
            double ph = pIn[dc]->GetBinCenter(ib);
            double p360 = ph; if (p360 < 0) p360 += 360.0;
            int iw = (int)(p360 / 30.0); if (iw > 11) iw = 11;
            if (pIn[dc]->GetBinEntries(ib) > 0){
                wIn[iw]  += pIn[dc]->GetBinContent(ib);  cIn[iw]  += 1;
            }
            if (pOut[dc]->GetBinEntries(ib) > 0){
                wOut[iw] += pOut[dc]->GetBinContent(ib); cOut[iw] += 1;
            }
            if (qIn[dc]->GetBinEntries(ib) > 0){
                yIn[iw]  += qIn[dc]->GetBinContent(ib);  kIn[iw]  += 1;
            }
            if (qOut[dc]->GetBinEntries(ib) > 0){
                yOut[iw] += qOut[dc]->GetBinContent(ib); kOut[iw] += 1;
            }
        }

        printf("\n FST%d   sector    r*dphi in   r*dphi out |    dr in     dr out\n", dc);
        double sEvenI = 0, sOddI = 0, sEvenO = 0, sOddO = 0;
        double nEvenI = 0, nOddI = 0, nEvenO = 0, nOddO = 0;
        double tEvenI = 0, tOddI = 0, tEvenO = 0, tOddO = 0;
        for (int wc = 0; wc < 12; wc++){
            double vi = cIn [wc] > 0 ? wIn [wc]/cIn [wc] : 0;
            double vo = cOut[wc] > 0 ? wOut[wc]/cOut[wc] : 0;
            double ui = kIn [wc] > 0 ? yIn [wc]/kIn [wc] : 0;
            double uo = kOut[wc] > 0 ? yOut[wc]/kOut[wc] : 0;
            printf("       [%3d,%3d)  %+9.4f   %+9.4f | %+9.4f  %+9.4f\n",
                   wc*30, wc*30+30, vi, vo, ui, uo);
            if (wc % 2 == 0){ sEvenI += vi; nEvenI++; sEvenO += vo; nEvenO++; tEvenI += ui; tEvenO += uo; }
            else            { sOddI  += vi; nOddI++;  sOddO  += vo; nOddO++;  tOddI  += ui; tOddO  += uo;  }
        }
        printf("   alternating amplitude  r*dphi: inner %+8.4f  outer %+8.4f  [cm]\n",
               0.5*(sEvenI/nEvenI - sOddI/nOddI), 0.5*(sEvenO/nEvenO - sOddO/nOddO));
        printf("   alternating amplitude      dr: inner %+8.4f  outer %+8.4f  [cm]\n",
               0.5*(tEvenI/nEvenI - tOddI/nOddI), 0.5*(tEvenO/nEvenO - tOddO/nOddO));
        printf("   z-stagger prediction       dr: inner ~0.17     outer ~0.05   [cm],\n");
        printf("   with the sign FLIPPING between disk 2 and disks 1,3.\n");
    }

    // ---- plots ----
    for (int dd = 0; dd < 3; dd++){
        TCanvas* cv = new TCanvas(Form("c_io_%d", dd), "", 900, 400);
        cv->Divide(2,1);
        cv->cd(1); pIn [dd]->SetLineColor(kBlue+1);  pIn [dd]->SetMarkerStyle(20);
        pIn [dd]->SetMarkerSize(0.5); pIn [dd]->Draw();
        cv->cd(2); pOut[dd]->SetLineColor(kRed+1);   pOut[dd]->SetMarkerStyle(20);
        pOut[dd]->SetMarkerSize(0.5); pOut[dd]->Draw();
        cv->SaveAs(Form("%s/fst%dInnerOuter_vsPhi.png", outdir, dd));
    }
    printf("\nWrote plots to %s\n", outdir);
}
