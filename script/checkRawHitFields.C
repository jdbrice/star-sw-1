// checkRawHitFields.C
//
// Reconciles a contradiction between two diagnostics this session:
// checkDisk1QuadA.C (direct hit->sector()/hit->rdo()/hit->feb() selection)
// found only 3 of 6 febs alive for disk1/quadA; checkUuidCoverage.C
// (decoding plane/quad/feb back out of fttDb->fob(hit)) found all 6 alive
// for the same file/plane/quad. Dumps both the RAW fields and the
// fob()-based decode side by side for the same hits to see where they
// diverge -- specifically whether hit->plane()/hit->quadrant() are already
// non-sentinel (i.e. pre-mapped by the MuDst conversion) before any
// ApplyHardwareMap() runs in this chain, which would make fttDb->plane(hit)
// use a DIFFERENT source than sector()/rdo().
//
// Usage: root4star -l -b -q 'checkRawHitFields.C("file.MuDst.root", 200)'

void loadLibsCRHF();

void checkRawHitFields(const Char_t *fileList = "st_fwd_23081004_raw_6000003.MuDst.root",
                        size_t nEvents = 200) {
    loadLibsCRHF();

    StChain *chain = new StChain("StChain");
    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", 1);
    TChain &muDstChain = *muDstMaker.chain();

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

    int nPrinted = 0;
    TH1I *hFebSectorRdo = new TH1I("hFebSectorRdo", "feb tally via sector==2,rdo==1;feb (hit->feb());nHit", 8, -1, 7);
    TH1I *hFebPlaneQuad = new TH1I("hFebPlaneQuad", "feb tally via fttDb::plane()==1,quadrant()==0;feb (hit->feb());nHit", 8, -1, 7);

    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) { printf("event %d: Make() non-OK, stopping\n", (int)iev); break; }

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttRawHit &raw = event->fttCollection()->rawHits();
        for (size_t ih = 0; ih < raw.size(); ih++) {
            StFttRawHit *hit = raw[ih];

            int rawPlaneField = (int)hit->plane();       // mapped field, sentinel 255 if unset
            int rawQuadField  = (int)hit->quadrant();     // mapped field, sentinel if unset
            int dbPlane = (int)fttDb->plane(hit);         // class method w/ fallback to sector()-1
            int dbQuad  = (int)fttDb->quadrant(hit);      // class method w/ fallback to rdo()-1

            if (nPrinted < 25) {
                printf("evt=%d sector=%d rdo=%d feb=%d | hit.plane()=%d hit.quad()=%d | fttDb.plane()=%d fttDb.quad()=%d\n",
                       (int)iev, (int)hit->sector(), (int)hit->rdo(), (int)hit->feb(),
                       rawPlaneField, rawQuadField, dbPlane, dbQuad);
                nPrinted++;
            }

            if (hit->sector() == 2 && hit->rdo() == 1) hFebSectorRdo->Fill(hit->feb());
            if (dbPlane == 1 && dbQuad == 0) hFebPlaneQuad->Fill(hit->feb());
        }
    }

    printf("\n=== feb tally via RAW sector()==2 && rdo()==1 ===\n");
    for (int b = 0; b <= 5; b++) printf("  feb=%d : %d\n", b, (int)hFebSectorRdo->GetBinContent(hFebSectorRdo->FindBin(b)));

    printf("\n=== feb tally via fttDb::plane(hit)==1 && fttDb::quadrant(hit)==0 ===\n");
    for (int b = 0; b <= 5; b++) printf("  feb=%d : %d\n", b, (int)hFebPlaneQuad->GetBinContent(hFebPlaneQuad->FindBin(b)));
}

void loadLibsCRHF() {
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
