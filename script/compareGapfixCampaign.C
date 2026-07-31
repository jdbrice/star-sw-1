// compareGapfixCampaign.C
//
// Side-by-side off-vs-on comparison for the outer-sensor gap-offset campaigns
// (condor/submit_gapfix_{off,on}.txt). Answers the one question the single-file
// test could not, for lack of statistics: does correcting the 1 deg
// outer-sensor displacement actually change the residual?
//
// Reports, for outer and inner hits separately, per disk:
//   step  = the jump in mean r*dphi across the wedge centreline. This is the
//           quantity a swapped gap offset should produce, because the two outer
//           sensors are displaced in OPPOSITE directions -- averaging them
//           together cancels it, which is what made an earlier version of this
//           test blind (see xihe_crosscheck_20260729.txt section 3).
//   rms   = spread of the 12 per-wedge means about their own average, i.e. the
//           size of the 12-fold structure itself.
//
// Both are measured in a constant-ANGLE window (|dphi| < 30 mrad). A fixed
// r*dphi window is a different angular cut at each radius and biases exactly
// the inner-vs-outer comparison (that mistake is recorded in
// fst_wedge_testA_results.txt section 4).
//
// Expectation if the fix matters for tracking: the OUTER step shrinks
// substantially from off to on, while the INNER numbers stay put (control).
// If the outer step is unchanged, the 1 deg displacement is real (the position
// test in gapfix_test/ proved that) but the residual is insensitive to it, and
// the geometry-vs-data tension is genuine rather than a statistics artefact.
//
// Usage:
//   root4star -b -q 'script/compareGapfixCampaign.C'
//   root4star -b -q 'script/compareGapfixCampaign.C("<offGlob>","<onGlob>")'

void gapfixOnePass(const char* glob, int pass,
                   double* stepOut, double* rmsOut, double* nOut);

void compareGapfixCampaign(
    const char* offGlob = "/gpfs01/star/pwg_tasks/FwdCalib/akio/gapfix_off/alignment/*.FwdAlignment_BLCVtx.root",
    const char* onGlob  = "/gpfs01/star/pwg_tasks/FwdCalib/akio/gapfix_on/alignment/*.FwdAlignment_BLCVtx.root"){

    // [pass][disk][region] flattened: idx = ((pass*3)+disk)*2 + region
    //   pass 0 = off, 1 = on;  region 0 = inner, 1 = outer
    double stepA[12], rmsA[12], nA[12];
    for (int ia = 0; ia < 12; ia++){ stepA[ia] = 0; rmsA[ia] = 0; nA[ia] = 0; }

    gapfixOnePass(offGlob, 0, stepA, rmsA, nA);
    gapfixOnePass(onGlob,  1, stepA, rmsA, nA);

    printf("\n=================================================================\n");
    printf(" GAP-FIX CAMPAIGN: off vs on\n");
    printf("=================================================================\n");
    printf(" centreline step in mean r*dphi [cm], constant-angle window\n\n");
    printf(" region  disk |      OFF              ON            change +- err\n");
    for (int ir = 1; ir >= 0; ir--){          // outer first -- it is the point
        for (int id = 0; id < 3; id++){
            int io = ((0*3)+id)*2 + ir;
            int in = ((1*3)+id)*2 + ir;
            double dch = stepA[in]-stepA[io];
            double edc = sqrt(nA[io]*nA[io] + nA[in]*nA[in]);
            printf("  %s   %d  | %+8.4f+-%.4f %+8.4f+-%.4f  %+8.4f+-%.4f%s\n",
                   ir ? "OUTER" : "inner", id,
                   stepA[io], nA[io], stepA[in], nA[in], dch, edc,
                   (edc > 0 && fabs(dch) > 3*edc) ? "  <-- >3sigma" : "");
        }
    }
    printf("\n 12-fold structure, RMS of the 12 wedge means about their mean [cm]\n\n");
    printf(" region  disk |    OFF       ON     change\n");
    for (int jr = 1; jr >= 0; jr--){
        for (int jd = 0; jd < 3; jd++){
            int jo = ((0*3)+jd)*2 + jr;
            int jn = ((1*3)+jd)*2 + jr;
            printf("  %s   %d  | %8.4f %8.4f  %+8.4f\n",
                   jr ? "OUTER" : "inner", jd, rmsA[jo], rmsA[jn], rmsA[jn]-rmsA[jo]);
        }
    }
    printf("\n Read: OUTER step shrinking off->on means the fix matters for\n");
    printf(" tracking. OUTER step unchanged means the displacement is real (the\n");
    printf(" position test already showed that) but the residual does not see it.\n");
    printf(" INNER is the control and must not move.\n");
}

void gapfixOnePass(const char* glob, int pass,
                   double* stepOut, double* rmsOut, double* nOut){

    TChain* tr = new TChain("alignTree");
    int nf = tr->Add(glob);
    if (nf <= 0){ printf("no files matched %s\n", glob); return; }

    // Completeness check. A job still writing, or one that died before
    // StFwdAlignmentMaker::Finish(), leaves a file with no alignTree; TChain
    // skips it silently and the sample is quietly smaller than it looks.
    // Report any such file by name -- with an off/on comparison a missing
    // file on one side only is worse than a smaller sample, it is a biased
    // one.
    int nbad = 0;
    TObjArray* fl = tr->GetListOfFiles();
    if (fl){
        for (int ifl = 0; ifl < fl->GetEntries(); ifl++){
            const char* fn2 = fl->At(ifl)->GetTitle();
            TFile* ftst = TFile::Open(fn2);
            if (!ftst || ftst->IsZombie()){
                printf("   BAD (unopenable): %s\n", fn2); nbad++;
                if (ftst) delete ftst;
                continue;
            }
            if (!ftst->Get("alignTree")){
                printf("   BAD (no alignTree): %s\n", fn2); nbad++;
            }
            delete ftst;
        }
    }

    Long64_t ne = tr->GetEntries();
    printf("pass %d: %d files (%d unusable), %lld rows  (%s)\n",
           pass, nf, nbad, ne, glob);
    if (nbad) printf("   WARNING: %d unusable file(s) -- exclude them from BOTH\n"
                     "   samples before quoting numbers, or the comparison is biased.\n", nbad);

    Int_t   bDet = 0, bNpt = 0;
    Float_t bHx = 0, bHy = 0, bHz = 0, bPx = 0, bPy = 0, bChi2 = 0;
    Bool_t  bConv = false;
    tr->SetBranchAddress("detType", &bDet);
    tr->SetBranchAddress("hitX",  &bHx);
    tr->SetBranchAddress("hitY",  &bHy);
    tr->SetBranchAddress("hitZ",  &bHz);
    tr->SetBranchAddress("projX", &bPx);
    tr->SetBranchAddress("projY", &bPy);
    // fit-quality branches: the centreline step is confounded by refit
    // coupling (the inner control moves), so chi2/ndf and the converged
    // fraction are the cleaner "does this help tracking" observables.
    tr->SetBranchAddress("chi2ndf",     &bChi2);
    tr->SetBranchAddress("nPointsUsed", &bNpt);
    tr->SetBranchAddress("converged",   &bConv);

    double diskZ[3]; diskZ[0] = 151.750; diskZ[1] = 165.248; diskZ[2] = 178.781;
    const double kCutRad = 0.030;    // constant-angle window

    // sums[disk][region][side]  side 0 = left of centreline, 1 = right
    // and per-wedge sums for the RMS.  All flat for CINT.
    double sSum[12], sCnt[12], sSq[12];        // (disk*2+region)*2 + side
    for (int ka = 0; ka < 12; ka++){ sSum[ka] = 0; sCnt[ka] = 0; sSq[ka] = 0; }
    double wSum[72], wCnt[72];                 // ((disk*2+region)*12) + wedge
    for (int kb = 0; kb < 72; kb++){ wSum[kb] = 0; wCnt[kb] = 0; }
    double chiSum = 0, chiCnt = 0, convSum = 0, allCnt = 0, nptSum = 0;
    // chi2 split by whether the REMOVED hit was inner or outer. Only outer
    // hits move when the gap fix is on, so a real improvement must
    // concentrate in the outer rows; the inner rows are a control that needs
    // no event-by-event join between the two samples.
    double chiR[2], chiN[2], chiSq[2];
    for (int kc = 0; kc < 2; kc++){ chiR[kc] = 0; chiN[kc] = 0; chiSq[kc] = 0; }
    // Same, sliced by nPointsUsed. The off/on samples differ slightly in
    // population (the fix changes which hits pass the cuts), so an inclusive
    // chi2 comparison could be a selection effect. Within a fixed
    // nPointsUsed slice the track composition is controlled, so a surviving
    // improvement is a genuine fit improvement.
    // index = region*20 + nPointsUsed, nPointsUsed clipped to [0,19]
    double cS[40], cN[40], cQ[40];
    for (int ke = 0; ke < 40; ke++){ cS[ke] = 0; cN[ke] = 0; cQ[ke] = 0; }

    for (Long64_t ie = 0; ie < ne; ie++){
        tr->GetEntry(ie);
        if (bDet != 0) continue;
        int dsel = -1;
        for (int db = 0; db < 3; db++) if (fabs(bHz - diskZ[db]) < 5.0) dsel = db;
        if (dsel < 0) continue;

        double rr   = sqrt(bHx*bHx + bHy*bHy);
        double phih = atan2(bHy, bHx);
        double phip = atan2(bPy, bPx);
        double dph  = phih - phip;
        while (dph >  TMath::Pi()) dph -= 2*TMath::Pi();
        while (dph < -TMath::Pi()) dph += 2*TMath::Pi();
        if (fabs(dph) > kCutRad) continue;
        double rdp = rr * dph;

        double p360 = phih*180.0/TMath::Pi();
        if (p360 < 0) p360 += 360.0;
        int sec = (int)(p360/30.0); if (sec > 11) sec = 11;
        double dLoc = p360 - (sec*30.0 + 15.0);

        allCnt += 1;
        if (bConv) convSum += 1;
        nptSum += bNpt;
        if (bConv && bChi2 > 0 && bChi2 < 100){ chiSum += bChi2; chiCnt += 1; }

        int reg  = (rr > 16.5) ? 1 : 0;
        if (bConv && bChi2 > 0 && bChi2 < 100){
            chiR[reg]  += bChi2; chiN[reg] += 1; chiSq[reg] += bChi2*bChi2;
            int np = bNpt; if (np < 0) np = 0; if (np > 19) np = 19;
            int ci = reg*20 + np;
            cS[ci] += bChi2; cN[ci] += 1; cQ[ci] += bChi2*bChi2;
        }
        int side = (dLoc > 0) ? 1 : 0;
        if (fabs(dLoc) > 1.0 && fabs(dLoc) < 14.0){
            int si = (dsel*2 + reg)*2 + side;
            sSum[si] += rdp; sCnt[si] += 1; sSq[si] += rdp*rdp;
        }
        int wi = ((dsel*2 + reg)*12) + sec;
        wSum[wi] += rdp; wCnt[wi] += 1;
    }

    printf("   pass %d fit quality: <chi2/ndf>=%.4f over %.0f converged refits, "
           "converged=%.2f%%, <nPointsUsed>=%.3f\n",
           pass, chiCnt > 0 ? chiSum/chiCnt : 0, chiCnt,
           allCnt > 0 ? 100.0*convSum/allCnt : 0,
           allCnt > 0 ? nptSum/allCnt : 0);
    for (int kd = 0; kd < 2; kd++){
        if (chiN[kd] < 100) continue;
        double m  = chiR[kd]/chiN[kd];
        double v  = chiSq[kd]/chiN[kd] - m*m; if (v < 0) v = 0;
        printf("     removed-hit %s: <chi2/ndf>=%.4f +- %.4f  over %.0f rows\n",
               kd ? "OUTER" : "inner", m, sqrt(v/chiN[kd]), chiN[kd]);
    }
    printf("     population-controlled (fixed nPointsUsed):\n");
    for (int np2 = 6; np2 <= 14; np2++){
        double ni = cN[0*20+np2], no = cN[1*20+np2];
        if (ni < 2000 || no < 2000) continue;
        double mi = cS[0*20+np2]/ni, mo = cS[1*20+np2]/no;
        double vi = cQ[0*20+np2]/ni - mi*mi; if (vi < 0) vi = 0;
        double vo = cQ[1*20+np2]/no - mo*mo; if (vo < 0) vo = 0;
        printf("       nPts=%2d  inner %.4f+-%.4f (%.0f)   OUTER %.4f+-%.4f (%.0f)\n",
               np2, mi, sqrt(vi/ni), ni, mo, sqrt(vo/no), no);
    }

    for (int dd = 0; dd < 3; dd++){
        for (int rg = 0; rg < 2; rg++){
            int iL = (dd*2 + rg)*2 + 0;
            int iR = (dd*2 + rg)*2 + 1;
            // step and its statistical error. Without this the off/on
            // differences cannot be told from refit noise -- and the inner
            // control moving by a comparable amount shows that noise is not
            // negligible.
            double step = 0, estep = 0;
            if (sCnt[iL] > 50 && sCnt[iR] > 50){
                double mL = sSum[iL]/sCnt[iL], mR = sSum[iR]/sCnt[iR];
                double vL = sSq[iL]/sCnt[iL] - mL*mL; if (vL < 0) vL = 0;
                double vR = sSq[iR]/sCnt[iR] - mR*mR; if (vR < 0) vR = 0;
                step  = mR - mL;
                estep = sqrt(vL/sCnt[iL] + vR/sCnt[iR]);
            }

            double mu[12]; double avg = 0, nn = 0;
            for (int wq = 0; wq < 12; wq++){
                int wj = ((dd*2 + rg)*12) + wq;
                mu[wq] = (wCnt[wj] > 50) ? wSum[wj]/wCnt[wj] : 0;
                if (wCnt[wj] > 50){ avg += mu[wq]; nn += 1; }
            }
            double rms = 0;
            if (nn >= 6){
                avg /= nn;
                double ss = 0;
                for (int wr = 0; wr < 12; wr++){
                    int wk = ((dd*2 + rg)*12) + wr;
                    if (wCnt[wk] > 50) ss += (mu[wr]-avg)*(mu[wr]-avg);
                }
                rms = sqrt(ss/nn);
            }
            int oi = ((pass*3)+dd)*2 + rg;
            stepOut[oi] = step;
            rmsOut[oi]  = rms;
            nOut[oi]    = estep;   // carry the error out in place of the raw count
        }
    }
    delete tr;
}
