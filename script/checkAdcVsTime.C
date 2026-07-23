// checkAdcVsTime.C
//
// Two things:
// 1) Checks the wraparound hypothesis directly: HitCalibHelper::time() has
//    a literal "// TODO: handle wrap around?" and does a naive
//    dbcid-anchor subtraction with no modular correction. For a uuid whose
//    true dbcid distribution straddles the 4095->0 edge (raw dbcid near
//    both ~0 and ~4095), the naive corrected time should show a
//    lopsided/split distribution; applying a standard shortest-circular-
//    distance correction (subtract/add 4096 when |diff|>2048) after the
//    fact should collapse it back into one clean peak. Dumps 1D corrected
//    time (naive vs wrap-fixed) for one specific uuid to check.
// 2) ADC vs corrected time (wrap-fixed), aggregated over all channels --
//    a first empirical look at the pulse shape now that per-VMM anchoring
//    is correct (post StFttDb::quadrant() fix).
//
// Usage: root4star -l -b -q 'checkAdcVsTime.C("file.MuDst.root", 2000, 350)'
// (350 = the uuid to inspect for the wraparound check; pick from
// checkUuidCoverage.C's table, or use raw2D.png to eyeball one that
// straddles the dbcid edge)

void loadLibsCAVT();

void checkAdcVsTime(const Char_t *fileList = "st_fwd_23081004_raw_6000003.MuDst.root",
                     size_t nEvents = 2000,
                     int uuidToInspect = 350,
                     const char *outdir = "20260722_adcVsTime") {
    loadLibsCAVT();
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

    TH1I *hUuidNaive = new TH1I("hUuidNaive", Form("uuid=%d corrected time, NAIVE (no wrap fix);time;n", uuidToInspect), 400, -4100, 4100);
    TH1I *hUuidWrap  = new TH1I("hUuidWrap",  Form("uuid=%d corrected time, WRAP-FIXED;time;n", uuidToInspect), 400, -1000, 1000);
    TH1I *hUuidDbcid = new TH1I("hUuidDbcid", Form("uuid=%d raw dbcid;dbcid;n", uuidToInspect), 410, 0, 4100);

    TH2F *hAdcVsTimeNaive = new TH2F("hAdcVsTimeNaive", "ADC vs corrected time, NAIVE, all channels;corrected time;adc", 400, -1000, 1000, 256, 0, 1024);
    TH2F *hAdcVsTimeWrap  = new TH2F("hAdcVsTimeWrap",  "ADC vs corrected time, WRAP-FIXED, all channels;corrected time;adc", 400, -1000, 1000, 256, 0, 1024);

    int nEventsSeen = 0;
    long nTotal = 0, nReady = 0;
    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) { printf("event %d: Make() non-OK, stopping\n", (int)iev); break; }
        nEventsSeen++;

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttRawHit &raw = event->fttCollection()->rawHits();
        for (size_t ih = 0; ih < raw.size(); ih++) {
            StFttRawHit *hit = raw[ih];
            nTotal++;

            UShort_t fob  = (UShort_t)fttDb->fob(hit);
            UShort_t uuid = hit->vmm() + (StFttDb::nVMMPerFob * fob);

            int t = (int)hit->time();
            if (t == -4097) continue;
            nReady++;

            int tWrap = t;
            if (tWrap > 2048) tWrap -= 4096;
            if (tWrap < -2048) tWrap += 4096;

            hAdcVsTimeNaive->Fill(t, hit->adc());
            hAdcVsTimeWrap->Fill(tWrap, hit->adc());

            if ((int)uuid == uuidToInspect) {
                hUuidNaive->Fill(t);
                hUuidWrap->Fill(tWrap);
                hUuidDbcid->Fill(hit->dbcid());
            }
        }
    }
    printf("processed %d events, nTotal=%ld nReady=%ld\n", nEventsSeen, nTotal, nReady);
    printf("uuid=%d: raw dbcid entries=%.0f, naive-corrected entries=%.0f, wrap-fixed entries=%.0f\n",
           uuidToInspect, hUuidDbcid->GetEntries(), hUuidNaive->GetEntries(), hUuidWrap->GetEntries());

    TFile *fout = new TFile(Form("%s/adcVsTime.root", outdir), "RECREATE");
    hUuidNaive->Write(); hUuidWrap->Write(); hUuidDbcid->Write();
    hAdcVsTimeNaive->Write(); hAdcVsTimeWrap->Write();
    fout->Close();

    TCanvas *c1 = new TCanvas("c1", "uuidDbcid", 700, 550);
    hUuidDbcid->SetFillColor(kGray);
    hUuidDbcid->Draw("hist");
    c1->SaveAs(Form("%s/uuidDbcid.png", outdir));

    TCanvas *c2 = new TCanvas("c2", "uuidNaive", 700, 550);
    hUuidNaive->SetFillColor(kGray);
    hUuidNaive->Draw("hist");
    c2->SaveAs(Form("%s/uuidCorrectedNaive.png", outdir));

    TCanvas *c3 = new TCanvas("c3", "uuidWrap", 700, 550);
    hUuidWrap->SetFillColor(kGray);
    hUuidWrap->Draw("hist");
    c3->SaveAs(Form("%s/uuidCorrectedWrapFixed.png", outdir));

    TCanvas *c4 = new TCanvas("c4", "adcVsTimeNaive", 900, 600);
    c4->SetRightMargin(0.13);
    c4->SetLogz(1);
    hAdcVsTimeNaive->Draw("colz");
    c4->SaveAs(Form("%s/adcVsTimeNaive.png", outdir));

    TCanvas *c5 = new TCanvas("c5", "adcVsTimeWrap", 900, 600);
    c5->SetRightMargin(0.13);
    c5->SetLogz(1);
    hAdcVsTimeWrap->Draw("colz");
    c5->SaveAs(Form("%s/adcVsTimeWrapFixed.png", outdir));

    printf("Done -- wrote plots to %s/\n", outdir);
}

void loadLibsCAVT() {
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
