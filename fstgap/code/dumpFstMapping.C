// dumpFstMapping.C
//
// Settles the outer-sensor question directly from the DB instead of by
// elimination between two hypotheses.
//
// The pairing that matters is (sensor, phiStrip), and it is NOT derivable from
// the code, because the two indices come from different places:
//
//   StFstRawHitMaker.cxx:520-539    elecId  = apvElecId + iChan
//                                   geoId   = mMappingVec[elecId]
//                                   setChannelId(elecId);  setGeoId(geoId);
//
//   StFstRawHit::getSensor()        from mChannelId  (electronics)
//   StFstRawHit::getPhiStrip()      from mGeoId      (geometry)
//   StFstRawHit::getRStrip()        from mGeoId
//
// so the elecId->geoId table Calibrations/fst/fstMapping (long mapping[36864])
// is the only thing that links them. This macro reads that table and tallies,
// for every one of the 36864 channels, which (sensor, phiStrip-half) pairs
// actually occur.
//
// Expected outcome if hypothesis B of the geometry check is right:
//   sensor 1 <-> phiStrip 64-127,  sensor 2 <-> phiStrip 0-63
// and the inner sensor 0 should cover the full 0-127.
//
// Usage: root4star -l -b -q 'script/dumpFstMapping.C'
//        root4star -l -b -q 'script/dumpFstMapping.C(20220401, 0)'

void loadLibsDFM();

void dumpFstMapping(int date = 20220401, int time = 0) {
    loadLibsDFM();

    StChain *chain = new StChain("StChain");

    St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
    dbMk->SetDateTime(date, time);

    StFstDbMaker *fstDbMk = new StFstDbMaker();
    chain->AddMaker(fstDbMk);

    Int_t iInit = chain->Init();
    if (iInit) { printf("ERROR: chain->Init() returned %d\n", iInit); return; }

    // the tables are filled in InitRun(), not Init()
    fstDbMk->InitRun(0);

    // NOTE: StFstDb inherits StObject, NOT TDataSet, so ToWhiteConst("fst_db")
    // wraps it in a St_ObjectSet -- GetDataSet() returns the WRAPPER and casting
    // that straight to StFstDb* silently yields garbage (getMapping() then comes
    // back null). Go through GetObject(). StFttDb is a TDataSet, which is why
    // the direct cast works in our StFttDb macros and not here.
    TDataSet *dsFst = chain->GetDataSet("fst_db");
    if (!dsFst) { printf("ERROR: no 'fst_db' dataset\n"); return; }
    StFstDb *fstDb = (StFstDb*) dsFst->GetObject();
    if (!fstDb) { printf("ERROR: 'fst_db' wrapper holds no object\n"); return; }

    const fstMapping_st *mp = fstDb->getMapping();
    if (!mp) { printf("ERROR: getMapping() returned null\n"); return; }
    printf("mapping table loaded for %d %d\n\n", date, time);

    // StFstConsts.h's const ints are not visible to CINT here, so mirror the
    // handful we need. Keep in sync with StRoot/StEvent/StFstConsts.h.
    const int kNElec      = 36864;   // kFstNumElecIds
    const int kStrInner   = 512;     // kFstNumStripsPerInnerSensor
    const int kStrOuter   = 256;     // kFstNumStripsPerOuterSensor
    const int kPhiSegWdg  = 128;     // kFstNumPhiSegPerWedge

    // Reproduce the accessors exactly (StFstRawHit.cxx:57-77, 94-99).
    // per-wedge block = 1*512 + 2*256 = 1024 channels
    const int kBlock  = 1*kStrInner + 2*kStrOuter;   // 1024

    // tally[sensor*2 + phiHalf]; flat, because CINT cannot handle a
    // function-local multi-dimensional array (it abandons the macro silently)
    int tally[6];
    int phiLo[3], phiHi[3], rLo[3], rHi[3];
    for (int ia = 0; ia < 3; ia++){
        tally[ia*2] = 0; tally[ia*2+1] = 0;
        phiLo[ia] = 9999; phiHi[ia] = -1;
        rLo[ia]   = 9999; rHi[ia]   = -1;
    }
    int nBad = 0;

    for (int elecId = 0; elecId < kNElec; elecId++){
        int geoId = (int) mp[0].mapping[elecId];
        if (geoId < 0 || geoId >= kNElec){ nBad++; continue; }

        // sensor, from elecId (this is what getSensor() does)
        int se = elecId % kBlock;
        int sensor;
        if (se < kStrInner) sensor = se / kStrInner;        // 0
        else                sensor = se / kStrOuter - 1;   // 1 or 2
        if (sensor < 0 || sensor > 2){ nBad++; continue; }

        // phiStrip and rStrip, from geoId
        int sg       = geoId % kBlock;
        int phiStrip = sg % kPhiSegWdg;
        int rStrip   = sg / kPhiSegWdg;

        tally[sensor*2 + (phiStrip < kPhiSegWdg/2 ? 0 : 1)]++;
        if (phiStrip < phiLo[sensor]) phiLo[sensor] = phiStrip;
        if (phiStrip > phiHi[sensor]) phiHi[sensor] = phiStrip;
        if (rStrip   < rLo[sensor])   rLo[sensor]   = rStrip;
        if (rStrip   > rHi[sensor])   rHi[sensor]   = rStrip;
    }

    printf("=================================================================\n");
    printf(" (sensor, phiStrip) pairing, over all %d channels\n", kNElec);
    printf("=================================================================\n");
    printf(" sensor is derived from elecId, phiStrip/rStrip from geoId = mapping[elecId]\n\n");
    printf(" sensor | channels with     channels with  | phiStrip | rStrip\n");
    printf("        | phiStrip 0-63    phiStrip 64-127 |  range   |  range\n");
    for (int ib = 0; ib < 3; ib++){
        printf("   %d    |   %8d          %8d      |  %3d-%3d |  %d-%d\n",
               ib, tally[ib*2], tally[ib*2+1],
               phiHi[ib] < 0 ? -1 : phiLo[ib], phiHi[ib],
               rHi[ib]   < 0 ? -1 : rLo[ib],   rHi[ib]);
    }
    if (nBad) printf("\n %d channels skipped (mapping out of range / bad sensor)\n", nBad);

    printf("\n Reading: sensor 0 is the inner sensor and should span phiStrip 0-127\n");
    printf(" and rStrip 0-3. Sensors 1 and 2 are the outer pair and should each\n");
    printf(" take ONE half of phiStrip, with rStrip 4-7. Which half each takes is\n");
    printf(" the answer we are after.\n");
}

void loadLibsDFM() {
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
