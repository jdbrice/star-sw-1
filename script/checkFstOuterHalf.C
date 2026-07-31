// checkFstOuterHalf.C
//
// Retest of the outer-sensor gap-offset question, after reading Xihe's
// alignment branch (~/xihe_align/star-sw-fwd, Alignment_xihe.pdf slide 9 and
// TrackFitter.h:474-495).
//
// WHY THE EARLIER TEST WAS WRONG
//
// script/checkFstWedgeGeom.C found the two outer FST sensors displaced by
// exactly 1.000 deg (= kFstStripGapPhi) from the AGML silicon, in OPPOSITE
// directions, and that flipping the sign of the +-0.5*kFstStripGapPhi term
// closes it to exactly 0.0000.
//
// script/checkFstInnerOuter.C then appeared to refute that, because outer
// hits showed no 30-deg-periodic structure. But that test averaged ALL outer
// hits in each 30 deg sector together -- and since the two outer halves are
// displaced in OPPOSITE directions, their errors cancel in that average. The
// test was blind to the effect it was meant to find.
//
// Xihe's independent reimplementation of the decode (TrackFitter.h) uses the
// same convention:
//     gapOffset = (stripPhi < halfWedgePhi) ? -0.5*gap : +0.5*gap
//     dphi_outer = stripSign * (edgeToCenterPhi - stripPhi + gapOffset)
// and reproduces the same 1 deg offset against AGML. He flags exactly this as
// unverified on slide 5: "The inner/outer gap and sign tables pass the current
// coordinate-closure plots, but there is no deterministic test covering every
// disk, wedge orientation, strip edge, and outer half."
//
// THE TEST
//
// Profile the mean r*dphi residual against dphiLocal -- the hit's azimuth
// relative to its own wedge centreline, in [-15,+15] deg. The predicted
// signature of a swapped gap offset is a STEP at dphiLocal = 0 of about
// 2 * 1.0 deg * r ~ 0.78 cm at r = 22.25 cm, present for OUTER hits only,
// with inner hits smooth across the centreline as a control.
//
// The global sense of the step flips with the wedge's front/back mounting, so
// even and odd sectors are kept separate; averaging them would cancel it
// again, which is the same mistake as before.
//
// Usage:
//   root4star -b -q 'script/checkFstOuterHalf.C'

void checkFstOuterHalf(
        const char* fname  = "/gpfs01/star/pwg_tasks/FwdCalib/akio/pico_20260704/alignment/*.FwdAlignment_BLCVtx.root",
        const char* outdir = "residual/data/alignment/outerhalf" ){

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

    double diskZ[3]; diskZ[0] = 151.750; diskZ[1] = 165.248; diskZ[2] = 178.781;

    // profiles vs dphiLocal, index = (disk*2 + parity)*2 + isOuter
    TProfile* pf[12];
    for (int ka = 0; ka < 12; ka++){
        int dk = ka/4, pa = (ka/2)%2, ou = ka%2;
        pf[ka] = new TProfile(Form("pf_d%d_p%d_o%d", dk, pa, ou),
                              Form("FST%d %s sectors, %s;#delta#phi_{local} [deg];<r#Delta#phi> [cm]",
                                   dk, pa ? "odd" : "even", ou ? "OUTER" : "inner"),
                              60, -15, 15, -3, 3);
        pf[ka]->SetDirectory(0);
    }

    Long64_t nEnt = tr->GetEntries();
    printf("Total rows: %lld\n", nEnt);

    for (Long64_t ie = 0; ie < nEnt; ie++){
        tr->GetEntry(ie);
        if (bDet != 0) continue;

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
        double rdp = rr * dph;

        double p360 = phih * 180.0 / TMath::Pi();
        if (p360 < 0) p360 += 360.0;
        int sec = (int)(p360 / 30.0); if (sec > 11) sec = 11;
        double dLocal = p360 - (sec*30.0 + 15.0);   // -15 .. +15 within the wedge

        int isOuter = (rr > 16.5) ? 1 : 0;
        int parity  = sec % 2;
        pf[(dsel*2 + parity)*2 + isOuter]->Fill(dLocal, rdp);

        if ((ie % 20000000) == 0 && ie > 0) printf("  ... %lld / %lld\n", ie, nEnt);
    }

    // ---- the step across the wedge centreline ----
    printf("\n==========================================================================\n");
    printf(" Step in <r*dphi> across the wedge centreline (dphiLocal = 0)\n");
    printf("==========================================================================\n");
    printf(" left  = mean over dphiLocal in [-14,-1] deg\n");
    printf(" right = mean over dphiLocal in [+1,+14] deg\n");
    printf(" A swapped gap offset predicts |step| ~ 0.78 cm for OUTER, ~0 for inner.\n\n");
    printf(" disk parity  region |   left      right     step\n");
    for (int kb = 0; kb < 12; kb++){
        int dk = kb/4, pa = (kb/2)%2, ou = kb%2;
        double sl = 0, nl = 0, sr = 0, nr = 0;
        for (int ib = 1; ib <= 60; ib++){
            double xc = pf[kb]->GetBinCenter(ib);
            if (pf[kb]->GetBinEntries(ib) < 50) continue;
            double vv = pf[kb]->GetBinContent(ib);
            if (xc < -1.0 && xc > -14.0){ sl += vv; nl += 1; }
            if (xc >  1.0 && xc <  14.0){ sr += vv; nr += 1; }
        }
        if (nl < 3 || nr < 3) continue;
        printf("   %d    %s    %s | %+8.4f  %+8.4f  %+8.4f\n",
               dk, pa ? "odd " : "even", ou ? "OUTER" : "inner",
               sl/nl, sr/nr, sr/nr - sl/nl);
    }

    // ---- plots ----
    for (int dc = 0; dc < 3; dc++){
        TCanvas* cv = new TCanvas(Form("c_oh_%d", dc), "", 1000, 700);
        cv->Divide(2,2);
        for (int kc = 0; kc < 4; kc++){
            cv->cd(kc+1);
            int idx = dc*4 + kc;
            pf[idx]->SetMarkerStyle(20); pf[idx]->SetMarkerSize(0.6);
            pf[idx]->SetLineColor(kBlue+1); pf[idx]->SetMarkerColor(kBlue+1);
            pf[idx]->Draw();
        }
        cv->SaveAs(Form("%s/fst%dOuterHalf.png", outdir, dc));
    }
    printf("\nWrote plots to %s\n", outdir);
}
