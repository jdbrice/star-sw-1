// checkFttTiming.C
//
// Quick check: does real per-hit FTT timing (bcid/dbcid -> calibrated time)
// actually come through non-trivially from the MuDst we already have, or is
// it always some default/sentinel value (which would mean a real time cut
// needs raw DAQ reprocessing, not just a config change to the existing
// afterburner chain)?
//
// Usage: root4star -l -b -q 'checkFttTiming.C("file.MuDst.root", 30)'

void loadLibsCFT();

void checkFttTiming(const Char_t *fileList = "st_fwd_23081019_raw_3000043.MuDst.root", size_t nEvents = 30) {
    loadLibsCFT();

    StChain *chain = new StChain("StChain");
    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", 1);
    TChain &muDstChain = *muDstMaker.chain();
    printf("MuDst file has %d events available in tree\n", muDstChain.GetEntries());

    St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
    StMuDst2StEventMaker *mu2ev = new StMuDst2StEventMaker();

    StFttDbMaker *fttDbMk = new StFttDbMaker();
    chain->AddMaker(fttDbMk);
    StFttHitCalibMaker *ftthcm = new StFttHitCalibMaker();
    // NOTE: no StFttClusterMaker here -- we want to inspect rawHits() BEFORE
    // any time-cut-gated clustering happens.

    Int_t iInit = chain->Init();
    if (iInit) chain->Fatal(iInit, "on init");

    size_t nEntries = muDstChain.GetEntries();
    if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

    int nPrinted = 0;
    int nTotal = 0, nSentinel = 0, nZero = 0;
    int timeMin = 999999, timeMax = -999999;

    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) break;

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttRawHit &raw = event->fttCollection()->rawHits();
        for (size_t ih = 0; ih < raw.size(); ih++) {
            StFttRawHit *hit = raw[ih];
            nTotal++;
            int t = (int)hit->time();
            if (t == -4097) nSentinel++;
            if (t == 0) nZero++;
            if (t < timeMin) timeMin = t;
            if (t > timeMax) timeMax = t;
            if (nPrinted < 20) {
                printf("hit: bcid=%d dbcid=%d tb=%d time=%d\n",
                       (int)hit->bcid(), (int)hit->dbcid(), (int)hit->tb(), (int)hit->time());
                nPrinted++;
            }
        }
    }

    printf("\n=== summary ===\n");
    printf("nTotal=%d  nSentinel(-4097)=%d (%.1f%%)  nZero=%d (%.1f%%)\n",
           nTotal, nSentinel, 100.0*nSentinel/nTotal, nZero, 100.0*nZero/nTotal);
    printf("time range: [%d, %d]\n", timeMin, timeMax);
}

void loadLibsCFT() {
    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();

    gSystem->Load("StarMagField");
    gSystem->Load("StMagF");
    gSystem->Load("StDetectorDbMaker");
    gSystem->Load("StTpcDb");
    gSystem->Load("StDaqLib");
    gSystem->Load("StDbBroker");
    gSystem->Load("StDbUtilities");
    gSystem->Load("St_db_Maker");

    gSystem->Load("StEvent");
    gSystem->Load("StEventMaker");

    gSystem->Load("St_base.so");
    gSystem->Load("StUtilities.so");
    gSystem->Load("libPhysics.so");

    gSystem->Load("StarClassLibrary");
    gSystem->Load("StStrangeMuDstMaker");
    gSystem->Load("StMuDSTMaker");

    gSystem->Load("StFttDbMaker");
    gSystem->Load("StFttHitCalibMaker");

    gSystem->Load("StStarLogger.so");
}
