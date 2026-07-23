// simFttResidualTest.C
//
// Full-chain test: FST fast-sim + FTT via StFttSimHitMaker (real StFttDb
// transform, not the GEANT-truth blur) + GenFit StFwdTrackMaker + StFwdResidualMaker,
// modeled on ~/fcstrk11/star-sw-fwd/script/sim.C but with the "fttSim" BFC tag
// (StFttFastSimMaker, disconnected geometry) replaced by the real chain:
//   StFttDbMaker -> StFttSimHitMaker -> StFttClusterMaker -> StFttClusterPointMaker
// inserted immediately before "fwdTrack" so GenFit sees real StFttPoints
// built through the actual StFttDb transform, the same way real data does.
//
// Usage:
//   root4star -b -q 'script/simFttResidualTest.C("ele45.vz0.run1.fzd", 20)'

void simFttResidualTest(const char* inFile = "ele45.vz0.run1.fzd", int nevents = 20, int debug = 0) {

    TString _geom = ""; // use the fGeom.root cache
    TString _chain = Form("fzin %s sdt20211016 fstFastSim fwdTrack MakeEvent StEvent McEvent ReverseField bigbig CMuDST tree", _geom.Data());
    printf("Chain: \n%s\n", _chain.Data());

    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile);

    cout << "Using the Geometry cache: fGeom.root (empty _geom tag + cwd fGeom.root)" << endl;

    gSystem->Load("libStFttDbMaker.so");
    gSystem->Load("libStFttSimHitMaker.so");
    gSystem->Load("libStFttClusterMaker.so");
    gSystem->Load("libStFttClusterPointMaker.so");
    gSystem->Load("libStFwdResidualMaker.so");

    // ------------------------------------------------------------------
    // Real FTT chain, inserted right before "fwdTrack" (StFwdTrackMaker)
    // so GenFit sees StFttPoints built through the real StFttDb transform.
    // ------------------------------------------------------------------
    StFttDbMaker* fttDbMk = new StFttDbMaker();
    chain->AddBefore("fwdTrack", fttDbMk);

    StFttSimHitMaker* fttSimHit = new StFttSimHitMaker();
    fttSimHit->SetDebug(debug);
    chain->AddBefore("fwdTrack", fttSimHit);

    StFttClusterMaker* fttClu = new StFttClusterMaker();
    fttClu->SetDebug(debug);
    fttClu->SetTimeCut(1 /*kTimeCutModeAcceptAll*/, -9999, 9999);
    chain->AddBefore("fwdTrack", fttClu);

    StFttClusterPointMaker* fttCP = new StFttClusterPointMaker();
    fttCP->SetDebug(debug);
    // mUseGeantData left at its constructor default (false) -- real
    // MakeLocalPoints/MakeGlobalPoints path, the whole point of this test.
    chain->AddBefore("fwdTrack", fttCP);

    // ------------------------------------------------------------------
    // Configure the forward tracker
    // ------------------------------------------------------------------
    StFwdTrackMaker* fwdTrack = (StFwdTrackMaker*)chain->GetMaker("fwdTrack");
    if (fwdTrack) {
        fwdTrack->SetDebug(debug);
        fwdTrack->setGeoCache("fGeom.root");
        fwdTrack->setSeedFindingWithFst();
        fwdTrack->setTrackRefit(true);

        fwdTrack->setFitDebugLvl(0);
        fwdTrack->setFitMinIterations(10);
        fwdTrack->setFitMaxIterations(20);
        fwdTrack->setDeltaPval(1e-1);
        fwdTrack->setRelChi2Change(1e-6);

        fwdTrack->setFttHitSource(1 /* StFwdHitLoader::STEVENT -- our real StFttClusterPointMaker output */);
        fwdTrack->setFstHitSource(0 /* StFwdHitLoader::GEANT -- fstFastSim's synthetic hits */);

        fwdTrack->setConfigKeyValue("TrackFitter:refit", true);
    } else {
        cout << "WARNING: fwdTrack maker not found!" << endl;
    }

    // ------------------------------------------------------------------
    // Residual analysis -- reads StEvent::fwdTrackCollection() after MuDst
    // is fully populated (same ordering requirement as the original sim.C)
    // ------------------------------------------------------------------
    TString outName(inFile);
    outName.ReplaceAll(".fzd", ".FwdDetResidual_Global.root");
    StFwdResidualMaker* fwdResidual = new StFwdResidualMaker(outName, "fwdResidual_Global");
    fwdResidual->setTrackType(0 /* Global -- unconstrained, least biased for MC */);
    chain->AddAfter("MuDst", fwdResidual);

    Int_t iInit = chain->Init();
    cout << "CHAIN INIT DONE? (good==0): " << iInit << endl;
    if (iInit) chain->Fatal(iInit, "on init");
    chain->PrintInfo();

    for (int i = 0; i < nevents; i++) {
        cout << "--------->START EVENT: " << i << endl;
        chain->Clear();
        if (kStOK != chain->Make()) break;
        cout << "<---------- END EVENT" << endl;
    }

    fwdResidual->Finish();
}
