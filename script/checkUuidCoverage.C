// checkUuidCoverage.C
//
// Follow-up on checkFttTimeVsUuid.C's plots: (1) which specific uuids are
// active at all, decoded back to (plane,quad,feb) to explain why only ~4
// bands show up out of 384 possible VMMs; (2) true global min/max raw
// dbcid (not clipped to the -100..200 display window) to know how wide a
// range is needed to see everything; (3) a fine, unbinned integer tally of
// CORRECTED time near the peak, to check whether the jaggy/comb look in
// the 1D projection is a real periodicity in the data or a binning/render
// artifact. Uses TH1 for all tallies (raw C-array + printf combo was
// silently dropping output under CINT for reasons not tracked down).
//
// Usage: root4star -l -b -q 'checkUuidCoverage.C("file.MuDst.root", 2000)'

void loadLibsCUC();

void checkUuidCoverage(const Char_t *fileList = "st_fwd_23081004_raw_6000003.MuDst.root",
                        size_t nEvents = 2000) {
    loadLibsCUC();

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

    TH1I *hUuid = new TH1I("hUuid", "hits per uuid;uuid;nHit", 500, 0, 500);
    TH1I *hDbcid = new TH1I("hDbcid", "raw dbcid;dbcid;nHit", 4200, -100, 4100);
    TH1I *hTimeFine = new TH1I("hTimeFine", "corrected time, fine;time;nHit", 61, -30.5, 30.5);

    int nEventsSeen = 0;
    long nTotal = 0;
    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) { printf("event %d: Make() non-OK, stopping\n", (int)iev); break; }
        nEventsSeen++;

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttRawHit &raw = event->fttCollection()->rawHits();
        for (size_t ih = 0; ih < raw.size(); ih++) {
            StFttRawHit *hit = raw[ih];
            nTotal++;

            UShort_t fob  = (UShort_t)fttDb->fob(hit);
            UShort_t uuid = hit->vmm() + (StFttDb::nVMMPerFob * fob);
            hUuid->Fill(uuid);
            hDbcid->Fill(hit->dbcid());

            int t = (int)hit->time();
            if (t != -4097) hTimeFine->Fill(t);
        }
    }
    printf("processed %d events, nTotal raw hits=%ld\n", nEventsSeen, nTotal);

    int firstBin = hDbcid->FindFirstBinAbove(0);
    int lastBin  = hDbcid->FindLastBinAbove(0);
    printf("global raw dbcid range with entries: [%.0f, %.0f]\n",
           hDbcid->GetXaxis()->GetBinLowEdge(firstBin), hDbcid->GetXaxis()->GetBinUpEdge(lastBin));

    int nActive = 0;
    for (int bA = 1; bA <= hUuid->GetNbinsX(); bA++) if (hUuid->GetBinContent(bA) > 0) nActive++;
    printf("active uuids (nHit>0): %d distinct out of 384 possible\n", nActive);

    printf("uuid,nHit,fob,plane,quad,feb,vmm\n");
    for (int bU = 1; bU <= hUuid->GetNbinsX(); bU++) {
        int nHitU = (int)hUuid->GetBinContent(bU);
        if (nHitU <= 0) continue;
        int uuidU = bU - 1; // bin bU covers [bU-1, bU)
        int vmmU = uuidU % 4;
        int fobU = uuidU / 4;
        int fobIdxU = fobU - 1;
        int planeU = fobIdxU / 24;
        int remU = fobIdxU % 24;
        int quadU = remU / 6;
        int febU  = remU % 6;
        printf("%d,%d,%d,%d,%d,%d,%d\n", uuidU, nHitU, fobU, planeU, quadU, febU, vmmU);
    }

    printf("fine unbinned tally of CORRECTED time, -30..30:\n");
    for (int bT = 1; bT <= hTimeFine->GetNbinsX(); bT++) {
        int tT = (int)TMath::Nint(hTimeFine->GetXaxis()->GetBinCenter(bT));
        printf("t=%d,%d\n", tT, (int)hTimeFine->GetBinContent(bT));
    }
}

void loadLibsCUC() {
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
