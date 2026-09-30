// fttMatchSummary.C
//
// One-line-per-campaign summary of how well FST-only (Global) tracks find
// their sTGC hits. Reads the two products every afterburner campaign writes:
//   <tag>_blinddiag.root : merged fwd_blind_diag.root  (excess + purity)
//   <tag>_<type>.root    : merged FwdDetResidual_<type> (FTT usage fraction),
//                          type = Global (default) / BLC / Primary / BLCVtx. Pass it
//                          as the 5th argument. Primary is what the blind residuals
//                          have actually been selecting since 2026-09-21 (fttDiagType=2),
//                          so quoting usage from _Global against a Primary residual
//                          compares two different track samples.
// Both are hadd merges of the per-job files in <outdir>/blinddiag and
// <outdir>/residual.
//
// Metrics, all for the FST-blind search (no FTT hit is on the track yet):
//   excess  = peak/local-sideband in the "All candidates" histogram at the
//             projection point (|d|<0.6cm vs 0.6-3cm), averaged over 4 disks
//   purity  = background-subtracted peak entries / all matched entries,
//             i.e. of the hits the tracker attaches, the fraction that are
//             the genuinely correct one
//   usage   = fraction of Global tracks with a hit used on that plane
//             (from the quadrant x charge histogram, Pos and Neg averaged)
//
// Usage:
//   root4star -l -b -q 'fttMatchSummary.C("merged_zf2022")'          // one
//   root4star -l -b -q 'fttMatchSummary.C("a","b","c","d")'          // up to 4
// CINT: unique loop variable names, fixed-size arrays, no std::vector.

// excess: peak |d|<0.6cm over sidebands 0.6-3cm      (as script/checkAllExcess.C)
// nReal : background-subtracted entries in |d|<1.0cm, sidebands 1-3cm (as script/estimatePurity.C)
void fmsPeakBg(TH1F *h, double &excess, double &nReal) {
    double ps = 0, ls = 0, hs = 0; int pn = 0, ln = 0, hn = 0;
    double psP = 0, lsP = 0, hsP = 0; int pnP = 0, lnP = 0, hnP = 0;
    for (int ib = 1; ib <= h->GetNbinsX(); ib++) {
        double x = h->GetXaxis()->GetBinCenter(ib), c = h->GetBinContent(ib);
        if (fabs(x) < 0.6)          { ps += c; pn++; }
        if (x >= -3.0 && x < -0.6)  { ls += c; ln++; }
        if (x >=  0.6 && x <  3.0)  { hs += c; hn++; }
        if (fabs(x) < 1.0)          { psP += c; pnP++; }
        if (x >= -3.0 && x < -1.0)  { lsP += c; lnP++; }
        if (x >=  1.0 && x <  3.0)  { hsP += c; hnP++; }
    }
    double bg   = (ls + hs) / (ln + hn);
    double bgP  = (lsP + hsP) / (lnP + hnP);
    excess = (bg > 0) ? (ps / pn) / bg : 0;
    nReal  = psP - bgP * pnP;
}

TString gFmsType = "Global";   // set by fttMatchSummary(); see the header note

void fmsOne(const char* tag) {
    TString fb = TString(tag) + "_blinddiag.root";
    TString fg = TString(tag) + "_" + gFmsType + ".root";
    TFile *b = TFile::Open(fb);
    if (!b || b->IsZombie()) { printf(">>> %-22s  MISSING %s\n", tag, fb.Data()); return; }

    double exV = 0, exH = 0, realSum = 0, matchSum = 0;
    for (int d = 0; d < 4; d++) {
        for (int io = 0; io < 2; io++) {
            const char* nAll = (io == 0) ? Form("hBlindDxAll_V_disk%d", d)     : Form("hBlindDyAll_H_disk%d", d);
            const char* nMat = (io == 0) ? Form("hBlindDxMatched_V_disk%d", d) : Form("hBlindDyMatched_H_disk%d", d);
            TH1F *hA = (TH1F*)b->Get(nAll), *hM = (TH1F*)b->Get(nMat);
            if (!hA || !hM) continue;
            double ex = 0, nReal = 0;
            fmsPeakBg(hA, ex, nReal);
            if (io == 0) exV += ex; else exH += ex;
            realSum  += nReal;
            matchSum += hM->Integral();
        }
    }
    double purity = (matchSum > 0) ? realSum / matchSum : 0;

    // FTT usage fraction for Global tracks, from the quadrant x charge histos
    double useTop = 0, useBot = 0; int nTop = 0, nBot = 0;
    TFile *g = TFile::Open(fg);
    if (g && !g->IsZombie()) {
        const char* qn[4] = {"Ntop", "Nbot", "Stop", "Sbot"};
        // 15-plane scheme of plotFwdResidualQuadCharge.C, bin = index + 2:
        // 0-2 FST1-3, then per FTT disk x,y,uv. uv is never filled, so take x,y only.
        const int fttXY[8] = {3, 4, 6, 7, 9, 10, 12, 13};
        for (int q = 0; q < 4; q++) {
            TH1F *hp = (TH1F*)g->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Pos", qn[q]));
            TH1F *hn = (TH1F*)g->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Neg", qn[q]));
            if (!hp || !hn) continue;
            double nP = hp->GetBinContent(20), nN = hn->GetBinContent(20);   // AllTrk
            if (nP + nN <= 0) continue;
            for (int ip = 0; ip < 8; ip++) {
                int bin = fttXY[ip] + 2;
                double f = (hp->GetBinContent(bin) + hn->GetBinContent(bin)) / (nP + nN);
                if (q == 0 || q == 2) { useTop += f; nTop++; } else { useBot += f; nBot++; }
            }
        }
    }
    printf(">>> %-22s excess V %.3f  H %.3f | purity %5.1f%% | FTT usage top %5.1f%% bot %5.1f%%\n",
           tag, exV / 4.0, exH / 4.0, 100 * purity,
           nTop ? 100 * useTop / nTop : 0, nBot ? 100 * useBot / nBot : 0);
}

void fttMatchSummary(const char* t1, const char* t2 = "", const char* t3 = "", const char* t4 = "",
                     const char* type = "Global") {
    gFmsType = type;
    printf("residual track type: %s\n", gFmsType.Data());
    printf(">>> FST-blind sTGC matching summary (Global = FST-only tracks)\n");
    fmsOne(t1);
    if (strlen(t2)) fmsOne(t2);
    if (strlen(t3)) fmsOne(t3);
    if (strlen(t4)) fmsOne(t4);
}
