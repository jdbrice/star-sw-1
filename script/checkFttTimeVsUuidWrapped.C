// checkFttTimeVsUuidWrapped.C
//
// Same as checkFttTimeVsUuid.C's CORRECTED plots, but with the
// HitCalibHelper::time() wraparound bug (dbcid-anchor, no modular
// correction -- see the literal "TODO: handle wrap around?" in
// HitCalibHelper.h) fixed as a post-process: any |naive time| > 2048 gets
// wrapped by +-4096 (dbcid is a 12-bit/4096-wide circular counter).
// Produces just the 2 "corrected, wrap-fixed" plots (2D + 1D) as separate
// files -- does NOT touch/overwrite the existing raw/corrected(unwrapped)
// plots.
//
// Usage: root4star -l -b -q 'checkFttTimeVsUuidWrapped.C("file.MuDst.root", 2000)'

void loadLibsCFTUW();

void checkFttTimeVsUuidWrapped(const Char_t *fileList = "st_fwd_23081004_raw_6000003.MuDst.root",
                                size_t nEvents = 2000,
                                const char *outdir = "20260722_timeVsUuid") {
    loadLibsCFTUW();
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

    const int nUuid = 400;
    TH2F *hCorr2DWrap = new TH2F("hCorr2DWrap", "CORRECTED time vs uuid, WRAP-FIXED;VMM uuid;corrected time",
                                  nUuid, 0, nUuid, 400, -1000, 1000);

    size_t nEntries = muDstChain.GetEntries();
    if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

    int nTotal = 0, nReady = 0;
    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) break;

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttRawHit &raw = event->fttCollection()->rawHits();
        for (size_t ih = 0; ih < raw.size(); ih++) {
            StFttRawHit *hit = raw[ih];
            nTotal++;

            UShort_t fob  = (UShort_t)fttDb->fob(hit);
            UShort_t uuid = hit->vmm() + (StFttDb::nVMMPerFob * fob);

            int t = (int)hit->time();
            if (t != -4097) {
                nReady++;
                int tWrap = t;
                if (tWrap > 2048) tWrap -= 4096;
                if (tWrap < -2048) tWrap += 4096;
                hCorr2DWrap->Fill(uuid, tWrap);
            }
        }
    }

    printf("\n=== summary ===\n");
    printf("nTotal=%d  nReady(corrected, non-sentinel)=%d (%.1f%%)\n",
           nTotal, nReady, 100.0*nReady/nTotal);

    TH1D *hCorrProjWrap = hCorr2DWrap->ProjectionY("hCorrProjWrap");
    hCorrProjWrap->SetTitle("CORRECTED time, WRAP-FIXED, all VMMs;corrected time;hits");

    TFile *fout = new TFile(Form("%s/fttTimeVsUuidWrapped.root", outdir), "RECREATE");
    hCorr2DWrap->Write();
    hCorrProjWrap->Write();
    fout->Close();

    TCanvas *c1 = new TCanvas("c1", "corr2Dwrap", 900, 600);
    c1->SetRightMargin(0.13);
    c1->SetLogz(1);
    hCorr2DWrap->Draw("colz");
    c1->SaveAs(Form("%s/corrected2D_wrapped.png", outdir));

    TCanvas *c2 = new TCanvas("c2", "corrProjWrap", 700, 550);
    hCorrProjWrap->SetFillColor(kGray);
    hCorrProjWrap->Draw("hist");
    c2->SaveAs(Form("%s/correctedProjY_wrapped.png", outdir));

    printf("Done -- wrote 2 wrap-fixed plots to %s/\n", outdir);
}

void loadLibsCFTUW() {
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
