// checkOrientationMap.C
//
// Cross-check for the online-QA vs offline-map disagreement on
// disk1(plane index 1)/quadrant A(index 0): online monitoring says the
// Horizontal (Y-measuring) plane is dead there, Vertical (X-measuring)
// alive; our offline 2D hit maps show the opposite (X dead, Y alive).
//
// This dumps, for real raw hits landing in plane=1/quad=0 only, a tally of
// counts by (feb, StFttDb::getOrientation() result) -- i.e. which FEBs are
// actually alive in the real data, and what orientation our code assigns
// them. If the FEBs that are actually alive (have real hit counts) are
// labeled kFttHorizontal by our code, that's the smoking gun: our code
// calls "alive" hits Horizontal when online says Vertical is the alive one
// (and vice versa) -- a real label swap, not just a real dead channel.
//
// Minimal chain: MuDst -> StEvent -> FTT raw hit mapping only (no
// clustering/pointmaking/tracking needed for this check).
//
// Usage: root4star -l -b -q 'checkOrientationMap.C("file.MuDst.root", 500)'

void loadLibsCOM();

void checkOrientationMap(const Char_t *fileList = "st_fwd_23081019_raw_6000036.MuDst.root", size_t nEvents = 500) {
    loadLibsCOM();

    StChain *chain = new StChain("StChain");

    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", 1);
    TChain &muDstChain = *muDstMaker.chain();
    printf("MuDst file has %d events available in tree\n", muDstChain.GetEntries());

    St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");

    StMuDst2StEventMaker *mu2ev = new StMuDst2StEventMaker();

    StFttDbMaker *fttDbMk = new StFttDbMaker();
    chain->AddMaker(fttDbMk);
    StFttHitCalibMaker *ftthcm = new StFttHitCalibMaker();
    StFttClusterMaker *fttclu = new StFttClusterMaker();
    fttclu->SetTimeCut(1, -9999, 9999); // accept-all: this check is about raw-hit mapping, not timing

    Int_t iInit = chain->Init();
    if (iInit) chain->Fatal(iInit, "on init");

    // === direct hardware-map-table probe for rob=5 (plane idx1/disk1, quad
    // idx0/A): does the loaded map have ANY entries for feb=2,4,6 (1-based
    // -- the febs missing from the raw-hit tally below), or are they simply
    // absent from the table entirely (a map-completeness gap, distinct from
    // getOrientation()'s H/V labeling logic)?
    StFttDb *fttDb = (StFttDb*)chain->GetDataSet("fttDb");
    if (fttDb) {
        printf("\n=== hardware map table probe, rob=5 (plane idx1/quad idx0=A) ===\n");
        int rob = 5;
        for (int feb1 = 1; feb1 <= 6; feb1++) {
            int nFound = 0;
            int oriSeen[5] = {0,0,0,0,0};
            for (int vmm = 1; vmm <= 4; vmm++) {
                for (int ch = 0; ch < 64; ch++) {
                    int row = -1, strip = -1;
                    UChar_t oriProbe = 4;
                    if (fttDb->hardwareMap(rob, feb1, vmm, ch, row, strip, oriProbe)) {
                        nFound++;
                        if (oriProbe < 5) oriSeen[oriProbe]++;
                    }
                }
            }
            printf("feb(1-based)=%d: map entries found=%d  (H=%d V=%d DiagH=%d DiagV=%d Unk=%d)\n",
                   feb1, nFound, oriSeen[0], oriSeen[1], oriSeen[2], oriSeen[3], oriSeen[4]);
        }
    } else {
        printf("\n=== could not retrieve fttDb dataset for map probe ===\n");
    }

    size_t nEntries = muDstChain.GetEntries();
    if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

    // tally[feb][orientation] -> count.  feb: 0-15 (headroom), orientation: 0=H,1=V,2=DiagH,3=DiagV,4=Unknown
    const int kMaxFeb = 16;
    const int kMaxOri = 5;
    int tally[kMaxFeb][kMaxOri];
    int tallyAllQuads[kMaxFeb][kMaxOri];
    for (int ifeb = 0; ifeb < kMaxFeb; ifeb++) for (int jori = 0; jori < kMaxOri; jori++) { tally[ifeb][jori] = 0; tallyAllQuads[ifeb][jori] = 0; }

    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) break;

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttRawHit &raw = event->fttCollection()->rawHits();
        for (size_t ih = 0; ih < raw.size(); ih++) {
            StFttRawHit *hit = raw[ih];
            if (hit->plane() != 1) continue;
            int feb = (int)hit->feb();
            int ori = (int)hit->orientation();
            if (feb < 0 || feb >= kMaxFeb || ori < 0 || ori >= kMaxOri) continue;
            tallyAllQuads[feb][ori]++;
            if (hit->quadrant() != 0) continue;
            tally[feb][ori]++;
        }
    }

    const char* oriName[5] = {"Horizontal", "Vertical", "DiagonalH", "DiagonalV", "Unknown"};

    printf("\n=== disk1(plane idx1)/quadA(idx0) raw hit tally by FEB and orientation ===\n");
    bool anyHits = false;
    for (int feb = 0; feb < kMaxFeb; feb++) {
        bool febHasHits = false;
        for (int o = 0; o < kMaxOri; o++) if (tally[feb][o] > 0) febHasHits = true;
        if (!febHasHits) continue;
        anyHits = true;
        printf("feb=%d:", feb);
        for (int o = 0; o < kMaxOri; o++) {
            if (tally[feb][o] > 0) printf("  %s=%d", oriName[o], tally[feb][o]);
        }
        printf("\n");
    }
    if (!anyHits) printf("  (no raw hits found in plane1/quad0 -- check event count / file)\n");

    printf("\n=== disk1(plane idx1), ALL quadrants, for context ===\n");
    for (int feb = 0; feb < kMaxFeb; feb++) {
        bool febHasHits = false;
        for (int o = 0; o < kMaxOri; o++) if (tallyAllQuads[feb][o] > 0) febHasHits = true;
        if (!febHasHits) continue;
        printf("feb=%d:", feb);
        for (int o = 0; o < kMaxOri; o++) {
            if (tallyAllQuads[feb][o] > 0) printf("  %s=%d", oriName[o], tallyAllQuads[feb][o]);
        }
        printf("\n");
    }
}

void loadLibsCOM() {
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
    gSystem->Load("StFttClusterMaker");

    gSystem->Load("StStarLogger.so");
}
