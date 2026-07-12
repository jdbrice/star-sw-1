// Small-scale test: run afterburner on simulation MuDst and produce PicoDst
// Tests that BLC vertex fields are stored correctly
// Usage: root4star -l -b -q test_blcvtx.C

bool runDb       = false;
bool runFttChain = true;
bool runFcsChain = false;
bool runFwdChain = true;
bool runPico     = true;

#include "StMemStat.h"

void loadLibs();

void test_blcvtx(const Char_t *fileList = "ele.pt1.MuDst.root",
                 size_t nEvents = 10, int debug = 0)
{
    cout << "BLC vertex test: " << fileList << " nEvents=" << nEvents << endl;

    // load all libs before creating any makers
    loadLibs();

    StChain *chain = new StChain("StChain");

    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", 1);
    TChain &muDstChain = *muDstMaker->chain();
    printf("MuDst has %d events\n", (int)muDstChain.GetEntries());

    new StMuDst2StEventMaker();

    if (runFttChain) {
        StFttDbMaker *fttDbMk = new StFttDbMaker();
        chain->AddMaker(fttDbMk);
        new StFttHitCalibMaker();
        StFttClusterMaker *fttclu = new StFttClusterMaker();
        fttclu->SetTimeCut(1, -40, 40);
        new StFttClusterPointMaker();
    }

    StFwdTrackMaker *fwdTrack = NULL;
    if (runFwdChain) {
        fwdTrack = new StFwdTrackMaker();
        fwdTrack->SetDebug(debug);
        fwdTrack->setGeoCache("fGeom.root");
        fwdTrack->setSeedFindingWithFst();
        fwdTrack->setTrackRefit(true);
        fwdTrack->setFitDebugLvl(0);
        fwdTrack->setFitMinIterations(40);
        fwdTrack->setFitMaxIterations(100);
        fwdTrack->setFstHitSource(2 /* MUDST */);
        fwdTrack->setFttHitSource(1 /* STEVENT */);
    }

    if (runPico) {
        StPicoDstMaker *picoMk = (StMaker *)(new StPicoDstMaker(
            StPicoDstMaker::IoWrite, fileList, "picoDst"));
        picoMk->setVtxMode(StPicoDstMaker::Vtxless);
    }

    chain->SetDebug(kError + debug);
    Int_t iInit = chain->Init();
    cout << "Chain Init: " << iInit << " (good==0)" << endl;
    if (iInit) chain->Fatal(iInit, "on init");
    chain->PrintInfo();

    size_t nEntries = muDstChain.GetEntries();
    if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

    for (int i = 0; i < (int)nEntries; i++) {
        chain->Clear();
        if (kStOK != chain->Make()) break;
        printf("EVENT #%d DONE\n", i);
    }

    chain->Finish();
    cout << "=== DONE ===" << endl;
}

void loadLibs()
{
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
    gSystem->Load("libGeom");
    gSystem->Load("St_g2t");
    gSystem->Load("libGeom.so");
    gSystem->Load("St_base.so");
    gSystem->Load("StUtilities.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("StarAgmlUtil.so");
    gSystem->Load("StarAgmlLib.so");
    gSystem->Load("libStarGeometry.so");
    gSystem->Load("libGeometry.so");
    gSystem->Load("xgeometry");
    gSystem->Load("St_geant_Maker");
    gSystem->Load("StarClassLibrary");
    gSystem->Load("StStrangeMuDstMaker");
    gSystem->Load("StMuDSTMaker");
    gSystem->Load("StBTofCalibMaker");
    gSystem->Load("StVpdCalibMaker");
    gSystem->Load("StBTofMatchMaker");
    gSystem->Load("StFcsDbMaker");

    gSystem->Load("StFttDbMaker");
    gSystem->Load("StFttHitCalibMaker");
    gSystem->Load("StFttClusterMaker");
    gSystem->Load("StFttClusterPointMaker");
    gSystem->Load("StFttPointMaker");
    gSystem->Load("libStarGeneratorUtil.so");
    gSystem->Load("libgenfit2");
    gSystem->Load("libKiTrack");
    gSystem->Load("libXMLIO.so");
    gSystem->Load("StFwdTrackMaker.so");
    gSystem->Load("StFwdUtils.so");
    gSystem->Load("libStEpdUtil.so");
    gSystem->Load("StStarLogger.so");
    gSystem->Load("libStPicoEvent");
    gSystem->Load("libStPicoDstMaker");
    gSystem->Load("StMuDst2StEventMaker");
}
