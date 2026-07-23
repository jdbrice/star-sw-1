// checkDisk1QuadA.C
//
// Focused check on disk1(plane index 1, hardware sector 2)/quadA(rdo 1):
// online QA reports VMM 5-8,13-16,21-24 (= feb 2,4,6, 1-based) totally
// empty for this plane/quad. Our StFttDb::getOrientation() parity rule
// (rob=quad+plane*nQuadPerPlane+1=5, odd) maps feb 2,4,6 (even) to
// kFttVertical/kFttDiagonalV -- i.e. our code's "X-measurement" strips.
// Directly tally real raw-hit occupancy per feb, and per-hit corrected time
// (sentinel vs valid, and whether it's inside a plausible time-cut window),
// to see (a) whether feb 2,4,6 really are hit-empty in our own data
// (independent hardware-dead confirmation) and (b) whether whatever DOES
// show up on feb 1,3,5 survives typical time cuts or gets cut away.
//
// Usage: root4star -l -b -q 'checkDisk1QuadA.C("file.MuDst.root", 2000)'

void loadLibsCD1QA();

void checkDisk1QuadA(const Char_t *fileList = "st_fwd_23081004_raw_6000003.MuDst.root",
                      size_t nEvents = 2000) {
    loadLibsCD1QA();

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

    StFttDb *fttDb = (StFttDb*)chain->GetDataSet("fttDb");
    if (!fttDb) { printf("ERROR: could not get fttDb dataset\n"); return; }

    size_t nEntries = muDstChain.GetEntries();
    if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

    // per-feb (0-5, 0-based) tallies, this plane/quad only
    long nHitFeb[6]     = {0,0,0,0,0,0};
    long nSentinelFeb[6] = {0,0,0,0,0,0};
    long nInWin4040[6]  = {0,0,0,0,0,0}; // |time| <= 40 (symmetric test window used earlier)
    long nInWin40100[6] = {0,0,0,0,0,0}; // -40..100 (class built-in default / online-QA value)
    int timeMinFeb[6], timeMaxFeb[6];
    for (int i = 0; i < 6; i++) { timeMinFeb[i] = 999999; timeMaxFeb[i] = -999999; }

    const int targetSector = 2; // hardware sector, 1-based -> plane index 1 -> disk1
    const int targetRdo    = 1; // rdo, 1-based -> quad index 0 -> quadA

    int nEventsSeen = 0;
    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) { printf("event %d: Make() non-OK, stopping\n", (int)iev); break; }
        nEventsSeen++;

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttRawHit &raw = event->fttCollection()->rawHits();
        for (size_t ih = 0; ih < raw.size(); ih++) {
            StFttRawHit *hit = raw[ih];
            if (hit->sector() != targetSector) continue;
            if (hit->rdo()    != targetRdo)    continue;

            int feb = hit->feb(); // 0-based, 0-5
            if (feb < 0 || feb > 5) continue;

            nHitFeb[feb]++;
            int t = (int)hit->time();
            if (t == -4097) { nSentinelFeb[feb]++; continue; }
            if (t < timeMinFeb[feb]) timeMinFeb[feb] = t;
            if (t > timeMaxFeb[feb]) timeMaxFeb[feb] = t;
            if (t >= -40 && t <= 40)   nInWin4040[feb]++;
            if (t >= -40 && t <= 100)  nInWin40100[feb]++;
        }
    }
    printf("processed %d events\n", nEventsSeen);

    // orientation per feb from our own code's parity rule (row=0, i.e.
    // rows 0-2 branch -- feb's orientation is row-independent except for
    // the row3/4 diagonal special case, so this correctly labels the
    // straight H/V assignment for this plane/quad/feb)
    int rob = (targetRdo - 1) + (targetSector - 1) * (int)StFttDb::nQuadPerPlane + 1;
    printf("\nplane index=%d (sector=%d), quad index=%d (rdo=%d), rob=%d\n",
           targetSector-1, targetSector, targetRdo-1, targetRdo, rob);

    printf("\n%-4s %-12s %-10s %-10s %-14s %-16s %-16s %-14s\n",
           "feb", "orientation(ours)", "nHit", "nSentinel", "nSentinel(%)", "nIn[-40,40]", "nIn[-40,100]", "time[min,max]");
    for (int feb1 = 1; feb1 <= 6; feb1++) {
        int febIdx = feb1 - 1;
        UChar_t ori = fttDb->getOrientation(rob, feb1, /*vmm*/1, /*row*/0);
        const char* oriName = (ori==0) ? "Horizontal" : (ori==1) ? "Vertical" : (ori==2) ? "DiagonalH" : (ori==3) ? "DiagonalV" : "Unknown";
        double pctSentinel = (nHitFeb[febIdx] > 0) ? 100.0*nSentinelFeb[febIdx]/nHitFeb[febIdx] : 0.0;
        printf("%-4d %-12s %-10ld %-10ld %-14.1f %-16ld %-16ld [%d,%d]\n",
               feb1, oriName, nHitFeb[febIdx], nSentinelFeb[febIdx], pctSentinel,
               nInWin4040[febIdx], nInWin40100[febIdx], timeMinFeb[febIdx], timeMaxFeb[febIdx]);
    }
}

void loadLibsCD1QA() {
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
