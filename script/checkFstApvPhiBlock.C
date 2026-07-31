// checkFstApvPhiBlock.C
//
// Explains the sharp spike at exactly nRawHitsPhi = 16 seen in real data
// (fstclus/index.html): is 16 contiguous phi strips a readout-hardware block?
//
// Each APV25 reads 128 channels and there are 8 per wedge (1024 channels =
// 1 inner sensor of 512 + 2 outer sensors of 256). This macro reads
// Calibrations/fst/fstMapping and reports, for one wedge, how many DISTINCT
// phi strips each APV covers and over what range -- i.e. the phi granularity
// of a single chip. If that number is 16, a chip-level noise/crosstalk burst
// produces exactly 16 adjacent phi strips, and the cluster maker's
// unconditional phi merge chains them into one cluster of nRawHitsPhi = 16.
//
// Usage: root4star -l -b -q 'script/checkFstApvPhiBlock.C'

void loadLibsFAP();

void checkFstApvPhiBlock(int date = 20220401, int time = 0) {
    loadLibsFAP();

    StChain *chain = new StChain("StChain");
    St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
    dbMk->SetDateTime(date, time);
    StFstDbMaker *fstDbMk = new StFstDbMaker();
    chain->AddMaker(fstDbMk);

    Int_t iInit = chain->Init();
    if (iInit) { printf("ERROR: chain->Init() returned %d\n", iInit); return; }
    fstDbMk->InitRun(0);

    // StFstDb inherits StObject, so the dataset is a wrapper -- go via GetObject()
    TDataSet *dsFst = chain->GetDataSet("fst_db");
    if (!dsFst) { printf("ERROR: no 'fst_db' dataset\n"); return; }
    StFstDb *fstDb = (StFstDb*) dsFst->GetObject();
    if (!fstDb) { printf("ERROR: wrapper holds no object\n"); return; }
    const fstMapping_st *mp = fstDb->getMapping();
    if (!mp) { printf("ERROR: getMapping() returned null\n"); return; }

    const int kNElec     = 36864;
    const int kStrInner  = 512;
    const int kStrOuter  = 256;
    const int kPhiSegWdg = 128;
    const int kApvChan   = 128;
    const int kBlock     = kStrInner + 2*kStrOuter;   // 1024 per wedge
    const int kApvPerWdg = 8;

    // seen[apv*128 + phiStrip], flat: CINT cannot handle a local 2D array
    int seen[1024];
    int rseen[64];                 // rseen[apv*8 + rStrip]
    for (int iz = 0; iz < 1024; iz++) seen[iz] = 0;
    for (int iy = 0; iy < 64; iy++)   rseen[iy] = 0;

    // wedge 0 only: elecId 0..1023
    for (int elecId = 0; elecId < kBlock; elecId++){
        int geoId = (int) mp[0].mapping[elecId];
        if (geoId < 0 || geoId >= kNElec) continue;

        int apv = (elecId % kBlock) / kApvChan;        // 0..7
        int sg  = geoId % kBlock;
        int phiStrip = sg % kPhiSegWdg;
        int rStrip   = sg / kPhiSegWdg;
        if (apv < 0 || apv >= kApvPerWdg) continue;
        if (phiStrip < 0 || phiStrip >= kPhiSegWdg) continue;
        if (rStrip < 0 || rStrip >= 8) continue;

        seen[apv*kPhiSegWdg + phiStrip] = 1;
        rseen[apv*8 + rStrip] = 1;
    }

    printf("\n================================================================\n");
    printf(" phi-strip coverage of each APV chip  (wedge 0)\n");
    printf("================================================================\n");
    printf("  apv | sensor | #distinct phiStrips | phi range  | rStrips\n");
    for (int ia = 0; ia < kApvPerWdg; ia++){
        int nphi = 0, lo = 9999, hi = -1;
        for (int ip = 0; ip < kPhiSegWdg; ip++){
            if (!seen[ia*kPhiSegWdg + ip]) continue;
            nphi++;
            if (ip < lo) lo = ip;
            if (ip > hi) hi = ip;
        }
        int nr = 0, rlo = 9999, rhi = -1;
        for (int ir = 0; ir < 8; ir++){
            if (!rseen[ia*8 + ir]) continue;
            nr++;
            if (ir < rlo) rlo = ir;
            if (ir > rhi) rhi = ir;
        }
        // sensor from elecId, same rule as StFstRawHit::getSensor()
        int se0 = ia*kApvChan;
        int sensor = (se0 < kStrInner) ? 0 : (se0 / kStrOuter - 1);
        printf("   %d  |   %d    |        %3d          |  %3d-%3d   |  %d-%d (%d)\n",
               ia, sensor, nphi, hi<0?-1:lo, hi, rhi<0?-1:rlo, rhi, nr);
    }
    printf("\n  A chip-wide burst therefore lights %s contiguous phi strips.\n",
           "the number in the '#distinct phiStrips' column");
}

// same order as script/dumpFstMapping.C -- StDbBroker must precede St_db_Maker
void loadLibsFAP() {
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
    gSystem->Load("StMuDSTMaker");

    gSystem->Load("libStFstDbMaker.so");

    gSystem->Load("StStarLogger.so");
}
