// checkFttTimeVsUuid.C
//
// 2D diagnostic for the FTT self-calibrating time correction
// (StFttHitCalibMaker / HitCalibHelper): x axis = per-(VMM,FOB) uuid
// (0-400, same uuid = vmm + StFttDb::nVMMPerFob*fob used internally by
// StFttHitCalibMaker::Make()), y axis = per-hit time.
//
// "RAW" = hit->dbcid() as unpacked, before any per-channel anchor
// subtraction.
// "CORRECTED" = hit->time() as set by StFttHitCalibMaker, i.e.
// dbcid - dbcidAnchor[uuid], only filled once that channel's helper is
// ready() (>200 samples) so cold-start sentinel (-4097) hits are excluded.
//
// Produces 4 plots: raw 2D, corrected 2D, and their Y-axis projections.
//
// Usage: root4star -l -b -q 'checkFttTimeVsUuid.C("file.MuDst.root", 2000)'

void loadLibsCFTU();

void checkFttTimeVsUuid(const Char_t *fileList = "st_fwd_23081019_raw_3000043.MuDst.root",
                         size_t nEvents = 2000,
                         const char *outdir = "20260722_timeVsUuid") {
    loadLibsCFTU();
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
    // NOTE: no StFttClusterMaker -- we want every raw hit, not just ones
    // that survive some time cut.

    Int_t iInit = chain->Init();
    if (iInit) chain->Fatal(iInit, "on init");

    StFttDb *fttDb = (StFttDb*)chain->GetDataSet("fttDb");
    if (!fttDb) { printf("ERROR: could not get fttDb dataset\n"); return; }

    const int nUuid = 400;
    // dbcid is a 12-bit unsigned counter, measured range [0,4096) -- give it
    // room to show the whole thing, not just an arbitrary -100..200 window.
    TH2F *hRaw2D  = new TH2F("hRaw2D",  "RAW dbcid vs uuid;VMM uuid;raw dbcid",
                              nUuid, 0, nUuid, 410, 0, 4100);
    // corrected time = dbcid - anchor can run both directions; widen past
    // the old +-100/200 window to see the full pedestal shape.
    TH2F *hCorr2D = new TH2F("hCorr2D", "CORRECTED time vs uuid;VMM uuid;corrected time",
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

            hRaw2D->Fill(uuid, hit->dbcid());

            int t = (int)hit->time();
            if (t != -4097) {
                nReady++;
                hCorr2D->Fill(uuid, t);
            }
        }
    }

    printf("\n=== summary ===\n");
    printf("nTotal=%d  nReady(corrected, non-sentinel)=%d (%.1f%%)\n",
           nTotal, nReady, 100.0*nReady/nTotal);

    TH1D *hRawProj  = hRaw2D->ProjectionY("hRawProj");
    TH1D *hCorrProj = hCorr2D->ProjectionY("hCorrProj");
    hRawProj->SetTitle("RAW dbcid, all VMMs;raw dbcid;hits");
    hCorrProj->SetTitle("CORRECTED time, all VMMs;corrected time;hits");

    TFile *fout = new TFile(Form("%s/fttTimeVsUuid.root", outdir), "RECREATE");
    hRaw2D->Write();
    hCorr2D->Write();
    hRawProj->Write();
    hCorrProj->Write();
    fout->Close();

    TCanvas *c1 = new TCanvas("c1", "raw2D", 900, 600);
    c1->SetRightMargin(0.13);
    c1->SetLogz(1);
    hRaw2D->Draw("colz");
    c1->SaveAs(Form("%s/raw2D.png", outdir));

    TCanvas *c2 = new TCanvas("c2", "corr2D", 900, 600);
    c2->SetRightMargin(0.13);
    c2->SetLogz(1);
    hCorr2D->Draw("colz");
    c2->SaveAs(Form("%s/corrected2D.png", outdir));

    TCanvas *c3 = new TCanvas("c3", "rawProj", 700, 550);
    hRawProj->SetFillColor(kGray);
    hRawProj->Draw("hist");
    c3->SaveAs(Form("%s/rawProjY.png", outdir));

    TCanvas *c4 = new TCanvas("c4", "corrProj", 700, 550);
    hCorrProj->SetFillColor(kGray);
    hCorrProj->Draw("hist");
    c4->SaveAs(Form("%s/correctedProjY.png", outdir));

    printf("Done -- wrote 4 plots + fttTimeVsUuid.root to %s/\n", outdir);
}

void loadLibsCFTU() {
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
