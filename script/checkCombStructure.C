// checkCombStructure.C
//
// Checks whether the ~25-unit-spaced comb/jaggy substructure seen in the
// aggregate dbcid/corrected-time plots is:
//  (a) present within a SINGLE channel's own dbcid distribution (fine bin
//      width=1, not lumped) -- tests whether it's a real periodicity at
//      all, vs purely an averaging/binning artifact of the aggregate plot.
//  (b) at the SAME absolute dbcid phase across DIFFERENT channels (would
//      indicate a shared, global cause common to the whole detector/beam,
//      e.g. RHIC's actual bunch/fill pattern) vs a different phase per
//      channel (would indicate a per-channel/electronics-local artifact).
//
// Picks a handful of well-populated, non-wraparound uuids automatically
// and dumps each one's own fine-binned dbcid histogram.
//
// Usage: root4star -l -b -q 'checkCombStructure.C("file.MuDst.root", 2000)'

void loadLibsCCS();

void checkCombStructure(const Char_t *fileList = "st_fwd_23081004_raw_6000003.MuDst.root",
                         size_t nEvents = 2000,
                         const char *outdir = "20260722_combStructure") {
    loadLibsCCS();
    gSystem->mkdir(outdir, true);
    gStyle->SetOptStat(0);

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

    // 6 hand-picked, well-separated, non-wraparound uuids (different
    // plane/quad/feb/vmm each) -- fine (bin width 1) dbcid histograms.
    const int nPick = 6;
    int pickUuid[nPick] = {101, 108, 250, 253, 300, 350};
    TH1I *hFine[nPick];
    for (int i = 0; i < nPick; i++) {
        hFine[i] = new TH1I(Form("hFine_%d", pickUuid[i]), Form("uuid=%d raw dbcid, fine bins;dbcid;n", pickUuid[i]), 4096, 0, 4096);
    }

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

            UShort_t fob  = (UShort_t)fttDb->fob(hit);
            UShort_t uuid = hit->vmm() + (StFttDb::nVMMPerFob * fob);

            for (int i = 0; i < nPick; i++) {
                if ((int)uuid == pickUuid[i]) hFine[i]->Fill(hit->dbcid());
            }
        }
    }
    printf("processed %d events\n", nEventsSeen);

    TFile *fout = new TFile(Form("%s/combStructure.root", outdir), "RECREATE");
    for (int i = 0; i < nPick; i++) {
        printf("uuid=%d entries=%.0f\n", pickUuid[i], hFine[i]->GetEntries());
        hFine[i]->Write();

        // zoom to the peak region: find peak bin, show +-150 around it
        int peakBin = hFine[i]->GetMaximumBin();
        double peakX = hFine[i]->GetXaxis()->GetBinCenter(peakBin);
        TCanvas *c = new TCanvas(Form("c_%d", pickUuid[i]), "c", 800, 500);
        hFine[i]->GetXaxis()->SetRangeUser(peakX - 150, peakX + 150);
        hFine[i]->SetFillColor(kGray);
        hFine[i]->Draw("hist");
        c->SaveAs(Form("%s/fine_uuid%d.png", outdir, pickUuid[i]));
    }
    fout->Close();
    printf("Done -- wrote plots to %s/\n", outdir);
}

void loadLibsCCS() {
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
