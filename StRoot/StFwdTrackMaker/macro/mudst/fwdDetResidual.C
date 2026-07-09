// Thin wrapper — all analysis is in compiled StFwdResidualMaker.
// Usage:
//   root4star -b -q 'fwdDetResidual.C("file.MuDst.root", 200)'
//   root4star -b -q 'fwdDetResidual.C("mudst.list", 5000, 10)'
//   root4star -b -q 'fwdDetResidual.C("file.MuDst.root", 200, 1, "fwdDetResidual.root", 1)'  // BLC tracks
//
// trackType (StFwdTrack::StFwdTrackType): 0=Global 1=BLC 2=Primary 3=FwdVtx 4=BLCVtx 5=FCSConstrained
//   Global (default) has no vertex constraint -> least biased for MC residuals.
//   For real data with a good TPC/BLC vertex, or B=0 straight-track runs, BLC or Primary
//   may give better statistics/cleaner residuals -- switchable here.

void loadLibs();

void fwdDetResidual(const char* fileList = "mudst.list",
                    int   nEvents        = 5000,
                    int   nFiles         = 1,
                    const char* outFile  = "fwdDetResidual.root",
                    int   trackType      = 0)
{
    loadLibs();

    StChain* chain = new StChain("StChain");
    StMuDstMaker* muMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", nFiles);
    printf("MuDst tree has %lld events\n", muMaker->chain()->GetEntries());

    StFwdResidualMaker* resMaker = new StFwdResidualMaker(outFile);
    resMaker->setTrackType((UChar_t)trackType);

    chain->Init();

    for (int i = 0; i < nEvents; i++) {
        chain->Clear();
        if (chain->Make() != kStOK) break;
        if (i % 200 == 0) printf("Event %d\n", i);
    }

    chain->Finish();
    delete chain;
}

void loadLibs() {
    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    gSystem->Load("StEvent");
    gSystem->Load("StStrangeMuDstMaker");
    gSystem->Load("StMuDSTMaker");
    gSystem->Load("StFwdResidualMaker.so");
}
