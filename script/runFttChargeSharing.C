// runFttChargeSharing.C -- drives the compiled StFttChargeSharingMaker over
// a real production MuDst file to measure the actual FTT charge-sharing
// profile (per-VMM online time calibration + re-clustering), replacing the
// flat 1:8:1 3-strip model used by StFttSlowSimMaker. See
// status_ftt_sim_maker.txt item 10.
// useDb=1 takes the per-VMM anchor AND per-VMM window from
// Calibrations/ftt/fttDataWindowsB. NOTE: StFttDbMaker prefers a local
// fttDataWindow/fttDataWindow.<run>.txt over the database, so run this from a
// directory that has no such file or you will silently test the file instead.
void runFttChargeSharing(const char* file = "st_fwd_23081015_raw_1000066.MuDst.root", int nevents = 2463,
                         int tlo = -65, int thi = 100, const char* outfile = "fttChargeSharing.root",
                         int useDb = 0, double crossTol = 1.6) {
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

    gSystem->Load("StFttDbMaker");

    StChain *chain = new StChain("chain");
    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", file, "MuDst.root", 1);

    St_db_Maker *dbMk = 0;
    if (useDb) {
        dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
        new StFttDbMaker();
    }

    StFttChargeSharingMaker *csMk = new StFttChargeSharingMaker();
    if (useDb) csMk->setUseDb(true);
    csMk->setCrossTolerance(crossTol);
    csMk->setOutputFile(outfile);
    csMk->setWarmupEvents(20);
    csMk->setTimeCutWindow(tlo, thi);
    csMk->setMinSamplesForReady(200);

    Int_t iInit = chain->Init();
    cout << "CHAIN INIT DONE? (good==0): " << iInit << endl;
    if (iInit) chain->Fatal(iInit, "on init");

    // Position the database at this run, exactly as fttDataWindow.C does: the
    // per-run stamped entry is only returned for the run it was measured on.
    if (useDb) {
        TChain &muChain = *muDstMaker->chain();
        muChain.GetEntry(0);
        int run = muDstMaker->muDst()->event()->runNumber();
        TTimeStamp ts((time_t)muDstMaker->muDst()->event()->eventInfo()->time(), 0);
        int date = (int)ts.GetDate(kTRUE), itime = (int)ts.GetTime(kTRUE);
        dbMk->SetDateTime(date, itime);
        chain->InitRun(run);
        printf("runFttChargeSharing: run %d, DB at %d %06d\n", run, date, itime);
    }

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
