// Read back PicoDst and print BLC vertex fields
// Usage: root4star -l -b -q read_blcvtx.C

void read_blcvtx(const Char_t *fname = "ele.pt1.picoDst.root") {
    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    gSystem->Load("StEvent");
    gSystem->Load("libStPicoEvent");

    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) { cout << "Cannot open " << fname << endl; return; }

    TTree *t = (TTree *)f->Get("PicoDst");
    if (!t) { cout << "No PicoDst tree" << endl; return; }

    TClonesArray *evArr = new TClonesArray("StPicoEvent", 1);
    t->SetBranchAddress("Event", &evArr);

    int nBLC = 0;
    Long64_t n = t->GetEntries();
    cout << "Events in tree: " << n << endl;
    cout << Form("%-6s  %-8s  %-8s  %-10s  %-10s  %-8s",
                 "Evt", "BLC_x", "BLC_y", "BLC_z", "sigmaZ", "nTrks") << endl;
    cout << string(58, '-') << endl;

    for (Long64_t i = 0; i < n; i++) {
        t->GetEntry(i);
        StPicoEvent *ev = (StPicoEvent *)evArr->At(0);
        if (!ev) continue;

        TVector3 blc = ev->blcVertex();
        Float_t  sz  = ev->blcVertexSigmaZ();
        UShort_t nt  = ev->blcVertexNTracks();

        cout << Form("%-6lld  %-8.3f  %-8.3f  %-10.3f  %-10.4f  %-8d",
                     i, blc.X(), blc.Y(), blc.Z(), sz, (int)nt) << endl;

        if (nt > 0) nBLC++;
    }

    cout << string(58, '-') << endl;
    cout << "Events with BLC vertex (nTracks>0): " << nBLC << " / " << n << endl;
    f->Close();
}
