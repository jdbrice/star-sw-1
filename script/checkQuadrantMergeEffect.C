// checkQuadrantMergeEffect.C
//
// Directly answers: under the StFttDb::quadrant() bug (fixed 2026-07-22),
// StFttHitCalibMaker's self-calibrating anchor (dbcidAnchor[uuid] = mode of
// that uuid's dbcid distribution) was being fed a HISTOGRAM MERGED FROM ALL
// 4 QUADRANTS sharing a given (plane,feb,vmm), since the buggy fob()
// always used quadrant=4 (out of range) instead of the real 0-3. This
// checks, for feb=0/vmm=0 on each plane: what is each TRUE quadrant's own
// dbcid peak (computed directly from hit->rdo()-1, bypassing the bug), and
// does the buggy MERGED peak (sum of all 4) just land on whichever single
// quadrant happens to have the most raw hits? I.e. is the corrupted
// "calibration" actually just "whichever quadrant is loudest wins, the
// other 3 silently get its anchor instead of their own"?
//
// Usage: root4star -l -b -q 'checkQuadrantMergeEffect.C("file.MuDst.root", 2000)'

void loadLibsCQME();

void checkQuadrantMergeEffect(const Char_t *fileList = "st_fwd_23081004_raw_6000003.MuDst.root",
                               size_t nEvents = 2000,
                               int febSel = 0, int vmmSel = 0) {
    loadLibsCQME();

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

    size_t nEntries = muDstChain.GetEntries();
    if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

    // [plane][quad] true per-quadrant dbcid histograms, feb=febSel,vmm=vmmSel only
    TH1I *hQ[4][4];
    TH1I *hMerged[4]; // sum over quad = what the buggy fob() would have fed the calibration
    for (int p = 0; p < 4; p++) {
        hMerged[p] = new TH1I(Form("hMerged_p%d", p), Form("plane%d feb%d vmm%d merged;dbcid;n", p, febSel, vmmSel), 4096, 0, 4096);
        for (int q = 0; q < 4; q++) {
            hQ[p][q] = new TH1I(Form("hQ_p%d_q%d", p, q), Form("plane%d quad%d feb%d vmm%d;dbcid;n", p, q, febSel, vmmSel), 4096, 0, 4096);
        }
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
            if ((int)hit->feb() != febSel) continue;
            if ((int)hit->vmm() != vmmSel) continue;

            int truePlane = (int)hit->sector() - 1;
            int trueQuad  = (int)hit->rdo() - 1;
            if (truePlane < 0 || truePlane > 3) continue;
            if (trueQuad  < 0 || trueQuad  > 3) continue;

            hQ[truePlane][trueQuad]->Fill(hit->dbcid());
            hMerged[truePlane]->Fill(hit->dbcid());
        }
    }
    printf("processed %d events\n", nEventsSeen);

    const char* qName[4] = {"A", "B", "C", "D"};
    printf("plane,quadA_peak,quadA_n,quadB_peak,quadB_n,quadC_peak,quadC_n,quadD_peak,quadD_n,merged_peak,merged_n,winner\n");
    for (int pl = 0; pl < 4; pl++) {
        int mergedPeakBin = hMerged[pl]->GetMaximumBin();
        int mergedPeak = (int)hMerged[pl]->GetXaxis()->GetBinCenter(mergedPeakBin);
        int mergedN = (int)hMerged[pl]->GetEntries();

        int qPeak[4]; int qN[4];
        for (int qd = 0; qd < 4; qd++) {
            int pb = hQ[pl][qd]->GetMaximumBin();
            qPeak[qd] = (int)hQ[pl][qd]->GetXaxis()->GetBinCenter(pb);
            qN[qd] = (int)hQ[pl][qd]->GetEntries();
        }

        // which true quadrant's own peak is closest to the merged peak?
        int bestQ = 0; int bestDiff = TMath::Abs(qPeak[0]-mergedPeak);
        for (int qd = 1; qd < 4; qd++) {
            int diff = TMath::Abs(qPeak[qd]-mergedPeak);
            if (diff < bestDiff) { bestDiff = diff; bestQ = qd; }
        }

        printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s(diff=%d)\n",
               pl, qPeak[0], qN[0], qPeak[1], qN[1], qPeak[2], qN[2], qPeak[3], qN[3],
               mergedPeak, mergedN, qName[bestQ], bestDiff);
    }
}

void loadLibsCQME() {
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
