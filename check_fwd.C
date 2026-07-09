// Check for FWD tracks in picoDst
void check_fwd(const Char_t *fname = "pico/st_fwd_23081004_raw_6000003.picoDst.root") {
    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    gSystem->Load("StEvent");
    gSystem->Load("libStPicoEvent");

    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", fname); return; }
    TTree *t = (TTree *)f->Get("PicoDst");
    if (!t) { printf("No PicoDst\n"); return; }

    printf("Branches:\n");
    t->GetListOfBranches()->ls();

    TClonesArray *fwdArr = new TClonesArray("StPicoFwdTrack", 100);
    t->SetBranchAddress("FwdTrack", &fwdArr);

    int nWithTracks = 0, nTotal = (int)t->GetEntries();
    int maxFwd = 0;
    for (int i = 0; i < nTotal; i++) {
        t->GetEntry(i);
        int nfwd = fwdArr->GetEntries();
        if (nfwd > 0) {
            nWithTracks++;
            if (nfwd > maxFwd) maxFwd = nfwd;
            if (nWithTracks <= 3)
                printf("  Event %d: %d FWD tracks\n", i, nfwd);
        }
    }
    printf("Events with FWD tracks: %d / %d  (max per event: %d)\n", nWithTracks, nTotal, maxFwd);
    f->Close();
}
