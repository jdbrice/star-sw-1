// checkAllExcess.C
//
// Section 1 check from proposal_next_step_20260721.txt: does the "All
// candidates" histogram (hBlindDxAll_V / hBlindDyAll_H -- the PRECISE
// coordinate for each orientation) show any excess over its LOCAL
// background right at the blind-projection center, BEFORE any nearest-hit
// selection is applied?
//
// IMPORTANT (found while writing this): the background is NOT flat across
// the full +/-15cm range -- it can have a strong slope/asymmetry (e.g. real
// disk0 H: ~4150 on one far side, ~1200 on the other). Comparing the peak
// against a symmetric average of DISTANT sidebands is invalid -- it picks
// up the slope, not a real peak. Use LOCAL sidebands immediately adjacent
// to the peak region instead, where the background is least likely to have
// moved much, and report the sideband asymmetry itself so a sloped
// background is visible rather than silently biasing the result.
//
// Usage: root4star -l -b -q 'checkAllExcess.C("fwd_blind_diag.root")'

void checkAllExcess(const char* infile = "fwd_blind_diag.root") {
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", infile); return; }

    const char* oriName[2] = {"V", "H"};
    for (int d = 0; d < 4; d++) {
        for (int io = 0; io < 2; io++) {
            const char* histName = (io == 0) ? Form("hBlindDxAll_V_disk%d", d)
                                              : Form("hBlindDyAll_H_disk%d", d);
            TH1F *h = (TH1F*)f->Get(histName);
            if (!h) { printf("Missing %s\n", histName); continue; }

            double peakSum = 0; int peakN = 0;
            double loSum = 0; int loN = 0;   // local sideband, -3 to -0.6
            double hiSum = 0; int hiN = 0;   // local sideband, +0.6 to +3
            for (int ib = 1; ib <= h->GetNbinsX(); ib++) {
                double x = h->GetXaxis()->GetBinCenter(ib);
                if (fabs(x) < 0.6)              { peakSum += h->GetBinContent(ib); peakN++; }
                if (x >= -3.0 && x < -0.6)      { loSum   += h->GetBinContent(ib); loN++; }
                if (x >=  0.6 && x <  3.0)       { hiSum   += h->GetBinContent(ib); hiN++; }
            }
            double peakPerBin = peakSum / peakN;
            double loPerBin = loSum / loN;
            double hiPerBin = hiSum / hiN;
            double bgPerBin = (loSum + hiSum) / (loN + hiN); // local, symmetric, close to peak
            double bgErrPerBin = sqrt(loSum + hiSum) / (loN + hiN);
            double excessSigma = (bgErrPerBin > 0) ? (peakPerBin - bgPerBin) / bgErrPerBin : 0;
            double sidebandAsym = (hiPerBin - loPerBin) / bgPerBin; // slope check

            printf("disk%d %s: local_bg/bin=%.1f+/-%.1f (lo=%.1f hi=%.1f, asym=%.1f%%)  "
                   "peak/bin=%.1f  excess=%.2fx  significance=%.1f sigma\n",
                   d, oriName[io], bgPerBin, bgErrPerBin, loPerBin, hiPerBin,
                   100*sidebandAsym, peakPerBin, peakPerBin/bgPerBin, excessSigma);
        }
    }
}
