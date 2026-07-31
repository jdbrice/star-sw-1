// checkFstStripInMuDst.C
//
// Is the FST strip-native measurement recoverable from a production MuDst, or
// is only the (buggy) derived phi available?
//
// StFstHitMaker computes local[1] = phi from the strip indices and the
// kFstphiStart/kFstzDirct/kFstzFilp/kFstStripGapPhi constants, and that value
// is what ends up in StMuFstHit::localPosition(1). If the gap-sign fix
// (StFstHitMaker.cxx:159,163) is needed, every outer-sensor phi in existing
// MuDst files is wrong by 1 deg -- baked in.
//
// BUT StMuFstHit also carries the inputs to that formula:
//     mMeanRStrip, mMeanPhiStrip, and mHardwarePosition (-> disk/wedge/sensor)
// so the position can be recomputed from scratch at afterburner load time,
// with no DAQ reprocessing. This macro checks those fields are actually
// populated and sane in real data, and shows the size of the correction by
// recomputing phi both ways per hit.
//
// (For MC MuDst they are NOT populated -- StFstFastSimMaker leaves them at -1.
//  See item 2 of the handoff note.)
//
// Usage:
//   root4star -l -b -q 'script/checkFstStripInMuDst.C'
//   root4star -l -b -q 'script/checkFstStripInMuDst.C("file.MuDst.root", 20)'

void loadLibsCFSM();

void checkFstStripInMuDst(
    const Char_t *fileList = "root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/production_pp500_2022/ReversedFullField/P24ia/2022/108/23108014/st_fwd_23108014_raw_2000026.MuDst.root",
    size_t nEvents = 20) {

    loadLibsCFSM();

    // StFstConsts.h values (not visible to CINT here)
    const int    kNRStripWdg = 8;
    const int    kNPhiSegWdg = 128;
    const int    kNWedgeDisk = 12;
    const double kPitchPhi   = TMath::Pi()*30.0/180.0/128.0;
    const double kGapPhi     = TMath::Pi()*1.0/180.0;
    const double kPitchR     = 2.875;
    const int    kzF[3]      = {1,-1,1};
    const int    kzD[12]     = {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1};
    const double kSt[12]     = {2.0,2.0,0.0,12.0,10.0,10.0,8.0,8.0,6.0,6.0,4.0,4.0};
    const double kSp[12]     = {3.0,1.0,1.0,11.0,11.0, 9.0,9.0,7.0,7.0,5.0,5.0,3.0};
    const double kRSt[8]     = {5.000,7.875,10.750,13.625,16.500,19.375,22.250,25.125};

    StChain *chain = new StChain("StChain");
    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", 1);
    TChain *muChain = muDstMaker->chain();
    printf("MuDst tree entries: %lld\n", muChain->GetEntries());

    Int_t iInit = chain->Init();
    if (iInit) { printf("ERROR: Init returned %d\n", iInit); return; }

    Long64_t nAvail = muChain->GetEntries();
    if (nAvail > (Long64_t)nEvents && nEvents > 0) nAvail = nEvents;

    int nHit = 0, nBadR = 0, nBadP = 0, nInner = 0, nOuter = 0;
    double rsLo = 1e9, rsHi = -1e9, psLo = 1e9, psHi = -1e9;
    double dMax = 0, dSum = 0; int dN = 0;

    for (Long64_t iev = 0; iev < nAvail; iev++){
        if (chain->Make(iev)) break;
        StMuFstCollection *fstc = StMuDst::muFstCollection();
        if (!fstc) continue;

        for (unsigned int ih = 0; ih < fstc->numberOfHits(); ih++){
            StMuFstHit *h = fstc->getHit(ih);
            if (!h) continue;
            nHit++;

            double rs = h->getMeanRStrip();
            double ps = h->getMeanPhiStrip();
            if (rs < rsLo) rsLo = rs;   if (rs > rsHi) rsHi = rs;
            if (ps < psLo) psLo = ps;   if (ps > psHi) psHi = ps;
            if (rs < 0 || rs >= kNRStripWdg)  nBadR++;
            if (ps < 0 || ps >= kNPhiSegWdg)  nBadP++;
            if (rs < 0 || rs >= kNRStripWdg || ps < 0 || ps >= kNPhiSegWdg) continue;

            int disk   = h->getDisk();                       // 1-3
            int wedge  = h->getWedge();                       // 1-36
            int sensor = h->getSensor();                      // 0-2
            int mIdx   = wedge - kNWedgeDisk*(disk-1);         // 1-12
            if (mIdx < 1 || mIdx > 12) continue;

            int dr = kzD[mIdx-1];
            int fl = kzF[disk-1];
            double phIn, phOut;
            if (disk == 2){
                phIn  = kSp[mIdx-1]*TMath::Pi()/6.0 - 0.5*dr*kPitchPhi;
                phOut = kSt[mIdx-1]*TMath::Pi()/6.0 + 0.5*dr*kPitchPhi;
            } else {
                phIn  = kSt[mIdx-1]*TMath::Pi()/6.0 + 0.5*dr*kPitchPhi;
                phOut = kSp[mIdx-1]*TMath::Pi()/6.0 - 0.5*dr*kPitchPhi;
            }

            int rStrip = (int)rs;
            double phiAsIs, phiFixed;
            if (rStrip < kNRStripWdg/2){          // inner: no gap term, unaffected
                phiAsIs  = phIn + fl*dr*ps*kPitchPhi;
                phiFixed = phiAsIs;
                nInner++;
            } else {                               // outer
                double sgn = (sensor == 1) ? -1.0 : +1.0;   // as in StFstHitMaker
                phiAsIs  = phOut - fl*dr*ps*kPitchPhi + sgn*fl*dr*0.5*kGapPhi;
                phiFixed = phOut - fl*dr*ps*kPitchPhi - sgn*fl*dr*0.5*kGapPhi;
                nOuter++;
            }

            // compare the recomputed as-is phi against what is stored, to prove
            // the recomputation reproduces StFstHitMaker
            double stored = h->localPosition(1);
            double dd = fabs(phiAsIs - stored);
            while (dd > TMath::Pi()) dd = fabs(dd - 2*TMath::Pi());
            if (dd > dMax) dMax = dd;
            dSum += dd; dN++;

            if (nHit <= 8){
                double rr = kRSt[rStrip] + 0.5*kPitchR;
                printf("  hit %d: disk=%d wedge=%2d sensor=%d rStrip=%.0f phiStrip=%5.1f  "
                       "stored=%+8.5f recomputed=%+8.5f  fixed=%+8.5f  shift=%+7.4f rad (%.3f cm at r=%.2f)\n",
                       nHit, disk, wedge, sensor, rs, ps, stored, phiAsIs, phiFixed,
                       phiFixed-phiAsIs, (phiFixed-phiAsIs)*rr, rr);
            }
        }
    }

    printf("\n=================================================================\n");
    printf(" FST hits seen: %d   (inner %d, outer %d)\n", nHit, nInner, nOuter);
    printf(" meanRStrip   range %.1f .. %.1f    out-of-range: %d\n", rsLo, rsHi, nBadR);
    printf(" meanPhiStrip range %.1f .. %.1f    out-of-range: %d\n", psLo, psHi, nBadP);
    if (dN) printf(" |recomputed - stored| phi:  mean %.3e rad, max %.3e rad, over %d hits\n",
                   dSum/dN, dMax, dN);
    printf("=================================================================\n");
    printf(" If meanRStrip/meanPhiStrip are in range and the recomputed phi\n");
    printf(" matches the stored phi, then the strip-native measurement IS\n");
    printf(" preserved in the MuDst and the gap-sign fix can be applied at\n");
    printf(" afterburner load time -- no DAQ reprocessing needed.\n");
}

void loadLibsCFSM() {
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
    gSystem->Load("StMuDSTMaker");
    gSystem->Load("StStarLogger.so");
}
