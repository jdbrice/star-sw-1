// checkPulseSamples.C
//
// Checks whether the DAQ actually provides more than one (adc,dbcid,tb)
// sample per physical channel per event -- i.e. whether a real per-hit
// "pulse shape" (multiple time samples of one channel's signal) exists at
// all, or whether each channel only ever reports a single peak-sensed
// (adc,time) pair per event (typical of VMM3-style peak/TAC readout, no
// waveform digitization). Groups raw hits within each event by their full
// physical identity (sector,rdo,feb,vmm,channel) via a sort (avoids
// std::map, which CINT has choked on elsewhere this session) and reports
// any group with more than 1 hit, printing the full (dbcid,tb,adc) list.
// Also dumps a handful of individual hits' raw fields for a first look.
//
// Usage: root4star -l -b -q 'checkPulseSamples.C("file.MuDst.root", 200)'

void loadLibsCPS();

void checkPulseSamples(const Char_t *fileList = "st_fwd_23081004_raw_6000003.MuDst.root",
                        size_t nEvents = 200) {
    loadLibsCPS();

    StChain *chain = new StChain("StChain");
    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", 1);
    TChain &muDstChain = *muDstMaker.chain();
    printf("MuDst file has %d events available in tree\n", muDstChain.GetEntries());

    St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
    StMuDst2StEventMaker *mu2ev = new StMuDst2StEventMaker();

    StFttDbMaker *fttDbMk = new StFttDbMaker();
    chain->AddMaker(fttDbMk);
    StFttHitCalibMaker *ftthcm = new StFttHitCalibMaker();

    Int_t iInit = chain->Init();
    if (iInit) chain->Fatal(iInit, "on init");

    size_t nEntries = muDstChain.GetEntries();
    if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

    const int maxHits = 30000;
    Long64_t keyArr[maxHits];
    int idxArr[maxHits];

    int nEventsSeen = 0;
    long nTotalHits = 0;
    long nMultiHitChannels = 0;
    int nMultiPrinted = 0;
    int nSingleDumped = 0;

    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) { printf("event %d: Make() non-OK, stopping\n", (int)iev); break; }
        nEventsSeen++;

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttRawHit &raw = event->fttCollection()->rawHits();
        int n = (int)raw.size();
        nTotalHits += n;
        if (n > maxHits) n = maxHits;

        for (int ih = 0; ih < n; ih++) {
            StFttRawHit *hit = raw[ih];
            Long64_t key = (Long64_t)hit->sector()*100000 + (Long64_t)hit->rdo()*10000 +
                           (Long64_t)hit->feb()*1000 + (Long64_t)hit->vmm()*100 + (Long64_t)hit->channel();
            keyArr[ih] = key;
        }
        TMath::Sort(n, keyArr, idxArr, kFALSE);

        int run = 1;
        for (int k = 1; k <= n; k++) {
            bool sameAsPrev = (k < n) && (keyArr[idxArr[k]] == keyArr[idxArr[k-1]]);
            if (sameAsPrev) {
                run++;
            } else {
                if (run > 1) {
                    nMultiHitChannels++;
                    if (nMultiPrinted < 10) {
                        printf("evt=%d key=%lld nHitsThisChannel=%d:\n", (int)iev, keyArr[idxArr[k-1]], run);
                        for (int r = 0; r < run; r++) {
                            StFttRawHit *h = raw[idxArr[k-1-r]];
                            printf("    dbcid=%d tb=%d adc=%d time=%d\n", (int)h->dbcid(), (int)h->tb(), (int)h->adc(), (int)h->time());
                        }
                        nMultiPrinted++;
                    }
                }
                run = 1;
            }
        }

        if (nSingleDumped < 15 && raw.size() > 0) {
            for (size_t ih = 0; ih < raw.size() && nSingleDumped < 15; ih++) {
                StFttRawHit *h = raw[ih];
                printf("evt=%d hit sec=%d rdo=%d feb=%d vmm=%d ch=%d | bcid=%d dbcid=%d tb=%d adc=%d time=%d\n",
                       (int)iev, (int)h->sector(), (int)h->rdo(), (int)h->feb(), (int)h->vmm(), (int)h->channel(),
                       (int)h->bcid(), (int)h->dbcid(), (int)h->tb(), (int)h->adc(), (int)h->time());
                nSingleDumped++;
            }
        }
    }

    printf("\n=== summary ===\n");
    printf("processed %d events, nTotalHits=%ld\n", nEventsSeen, nTotalHits);
    printf("channel-instances with >1 hit in the SAME event: %ld\n", nMultiHitChannels);
}

void loadLibsCPS() {
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
