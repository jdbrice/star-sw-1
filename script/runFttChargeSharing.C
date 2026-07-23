// runFttChargeSharing.C -- drives the compiled StFttChargeSharingMaker over
// a real production MuDst file to measure the actual FTT charge-sharing
// profile (per-VMM online time calibration + re-clustering), replacing the
// flat 1:8:1 3-strip model used by StFttSimHitMaker. See
// status_ftt_sim_maker.txt item 10.
void runFttChargeSharing(const char* file = "st_fwd_23081015_raw_1000066.MuDst.root", int nevents = 2463) {
    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();

    gSystem->Load("StarMagField");
    gSystem->Load("StMagF");
    gSystem->Load("StDetectorDbMaker");
    gSystem->Load("StTpcDb");
    gSystem->Load("StDbBroker");
    gSystem->Load("StDbUtilities");
    gSystem->Load("St_db_Maker");

    gSystem->Load("St_base.so");
    gSystem->Load("StUtilities.so");
    gSystem->Load("StEvent.so");
    gSystem->Load("StEventMaker.so");
    gSystem->Load("StarClassLibrary");
    gSystem->Load("StStrangeMuDstMaker");
    gSystem->Load("StMuDSTMaker");

    gSystem->Load("libStFttChargeSharingMaker.so");

    StChain *chain = new StChain("chain");
    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", file, "MuDst.root", 1);
    StFttChargeSharingMaker *csMk = new StFttChargeSharingMaker();
    csMk->setOutputFile("fttChargeSharing.root");
    csMk->setWarmupEvents(20);
    csMk->setTimeCutWindow(-65, 100);
    csMk->setMinSamplesForReady(200);

    Int_t iInit = chain->Init();
    cout << "CHAIN INIT DONE? (good==0): " << iInit << endl;
    if (iInit) chain->Fatal(iInit, "on init");

    int nOk = 0;
    for (int i = 0; i < nevents; i++) {
        chain->Clear();
        if (kStOK != chain->Make()) { printf("event %d: Make() returned non-OK, stopping\n", i); break; }
        nOk++;
        if (i > 0 && i % 200 == 0) printf("...processed %d events\n", i);
    }
    printf("processed %d events total\n", nOk);

    chain->Finish();
}
