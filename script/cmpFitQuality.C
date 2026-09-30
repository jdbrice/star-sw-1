// cmpFitQuality.C -- A/B fit quality for the FTT-on-track debugging run.
//
// The point of that run is NOT a residual width: it is whether the planar FTT
// measurement path behaves once GenFit actually uses it. So compare, on the SAME
// files and the same events, the control arm (fttNoAdd=1, no sTGC hit ever on a
// track) against the test arm (fttNoAdd=0):
//
//   chi2            a jump means the FTT hits do not agree with the FST+vertex fit
//                   (wrong plane, wrong covariance, or a real misalignment)
//   nSeedPoints     how many points went in -- the pickup shows up here
//   nFitPoints      how many survived the fit; |value|, since it is charge*nHitsFit
//   ntracks         a drop means fits started failing outright
//
// Control baseline measured on pico_iter2_20260930, Primary tracks, one file:
//   chi2 4.335, nSeed 4.271, 45132 tracks  (nSeed ~ 3 FST disks + vertex). The chi2
//   mean depends on the histogram range, so compare only against this macro's own.
//
// Usage (wildcards fine, they go into a TChain):
//   root -l -b -q 'cmpFitQuality.C("ctrl/*.picoDst.root","test/*.picoDst.root")'
//   root -l -b -q 'cmpFitQuality.C("a.root","b.root",2,"plots/fitq.png")'
//
// CINT: globals instead of reference args, unique loop variable names (see CLAUDE.md).

double fqN, fqChi2, fqSeed, fqFit, fqChi2Rms;

void fqRead(const char* glob, int type) {
    fqN = 0; fqChi2 = 0; fqSeed = 0; fqFit = 0; fqChi2Rms = 0;
    TChain* ch = new TChain("PicoDst");
    if (ch->Add(glob) == 0) { printf("  no files match %s\n", glob); return; }
    TString cut = Form("(FwdTracks.mVtxIndex & 7)==%d", type);

    TH1F* hq = new TH1F(Form("fqq%d", type), "chi2", 400, 0, 200);
    TH1F* hs = new TH1F(Form("fqs%d", type), "seed", 40, 0, 40);
    TH1F* hf = new TH1F(Form("fqf%d", type), "fit",  81, -40, 41);
    ch->Draw(Form("FwdTracks.mChi2>>fqq%d", type), cut, "goff");
    ch->Draw(Form("FwdTracks.mNumberOfSeedPoints>>fqs%d", type), cut, "goff");
    ch->Draw(Form("abs(FwdTracks.mNumberOfFitPoints)>>fqf%d", type), cut, "goff");
    fqN       = hq->GetEntries();
    fqChi2    = hq->GetMean();
    fqChi2Rms = hq->GetRMS();
    fqSeed    = hs->GetMean();
    fqFit     = hf->GetMean();
}

void cmpFitQuality(const char* ctrl, const char* test, int type = 2,
                   const char* out = "FstFttFlipTest/plots/fit_quality.png")
{
    const char* tn[6] = {"Global","BLC","Primary","FwdVtx","BLCVtx","FCSConstrained"};
    printf("\n=== fit quality A/B, %s tracks ===\n", (type < 6) ? tn[type] : "?");

    printf("\n control (fttNoAdd=1): %s\n", ctrl);
    fqRead(ctrl, type);
    double cN = fqN, cQ = fqChi2, cQr = fqChi2Rms, cS = fqSeed, cF = fqFit;
    printf("\n test    (fttNoAdd=0): %s\n", test);
    fqRead(test, type);
    double tN = fqN, tQ = fqChi2, tQr = fqChi2Rms, tS = fqSeed, tF = fqFit;

    printf("\n  %-14s %14s %14s %12s\n", "quantity", "control", "test", "change");
    printf("  %-14s %14.0f %14.0f %11.1f%%\n", "tracks", cN, tN,
           (cN > 0) ? 100.0*(tN-cN)/cN : 0);
    printf("  %-14s %14.3f %14.3f %11.1f%%\n", "chi2 (mean)", cQ, tQ,
           (cQ > 0) ? 100.0*(tQ-cQ)/cQ : 0);
    printf("  %-14s %14.3f %14.3f %11.1f%%\n", "chi2 (rms)", cQr, tQr,
           (cQr > 0) ? 100.0*(tQr-cQr)/cQr : 0);
    printf("  %-14s %14.3f %14.3f %+11.3f\n", "nSeedPoints", cS, tS, tS-cS);
    printf("  %-14s %14.3f %14.3f %+11.3f\n", "|nFitPoints|", cF, tF, tF-cF);

    printf("\n  read: nSeedPoints rising by ~the pickup is the change being looked for.\n");
    if (tS - cS < 0.01)
        printf("  nSeedPoints did NOT rise -- the FTT hits never reached the fit. Check\n"
               "  fttNoAdd actually reached the tracker (grep 'fttNoAdd' in the job log).\n");
    if (cQ > 0 && tQ > 2.0*cQ)
        printf("  chi2 more than doubled -- suspect the FTT plane/covariance, NOT alignment:\n"
               "  this is the first run in which the planar FTT measurement is used at all.\n");
    if (cN > 0 && tN < 0.95*cN)
        printf("  track count dropped >5%% -- fits are failing; look for genfit exceptions.\n");
}
