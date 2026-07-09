// computeWedgeCorrection.C
// Computes the per-(disk, wedge-phi-bin) mean unbiased dphi correction table
// consumed by FstWedgeAligner (StRoot/StFwdTrackMaker/StFwdHitLoader.h) --
// see bugreport_StFstHitMaker.txt for why this correction exists (official
// StFstHitMaker.cxx computes but discards the real per-sensor DB alignment,
// leaving an uncorrected ~30-degree-periodic phi offset per wedge).
//
// Usage: root4star -b -q 'script/computeWedgeCorrection.C' from the repo
// root, with DataDisk/pico/alignment/*.FwdAlignment_BLCVtx.root populated
// (see proposal_rerun_data_picodst.txt for how that campaign was produced).
// Paste the printed kFstWedgePhiCorrection[][] array into
// FstWedgeAligner::kWedgePhiCorrection to update the correction with more
// stats later.

#include <algorithm>

int findZClusters(TTree* t, const char* cutExpr, double gapCm, double* lo, double* hi, double* mean) {
    t->SetEstimate(t->GetEntries() + 10);
    Long64_t n = t->Draw("hitZ", cutExpr, "goff");
    if (n <= 0) return 0;
    double* z = t->GetV1();
    double* zs = new double[n];
    for (Long64_t i = 0; i < n; i++) zs[i] = z[i];
    std::sort(zs, zs + n);
    int nc = 0; Long64_t start = 0;
    for (Long64_t i = 1; i <= n; i++) {
        if (i == n || (zs[i] - zs[i-1]) > gapCm) {
            lo[nc] = zs[start] - 0.5; hi[nc] = zs[i-1] + 0.5;
            double s = 0; for (Long64_t k = start; k < i; k++) s += zs[k];
            mean[nc] = s / (i - start);
            nc++; start = i;
        }
    }
    delete[] zs;
    return nc;
}

// Fixed 30-degree phi bins, phase-locked to the boundaries found empirically
// (0,+-30,+-60,...,+-180) -- bin 0 covers [-180,-150), bin index increases
// with phi. This sidesteps needing to reproduce kFstphiStart's
// non-monotonic moduleIdx->phi mapping: we just correct as a function of the
// hit's own measured phi, which is equivalent (each of the 12 measured phi
// bins is exactly one physical wedge, whatever its internal index label is).
int wedgeBinOf(double phiDeg) {
    // Boundary phase confirmed empirically 2026-07-07: real step-jump
    // locations cluster at 0, +-30, +-60, ... +-180 deg (see
    // StFwdHitLoader.h's FstWedgeAligner::wedgePhiBin comment) -- no extra
    // +15 deg offset (that was a bug in an earlier version of this function).
    double shifted = phiDeg + 180.0;
    while (shifted < 0) shifted += 360.0;
    while (shifted >= 360.0) shifted -= 360.0;
    return (int)(shifted / 30.0);
}

void computeWedgeCorrection() {
    TChain* t = new TChain("alignTree");
    int nFiles = t->Add("DataDisk/pico/alignment/*.FwdAlignment_BLCVtx.root");
    printf("nFiles=%d\n", nFiles);

    double lo[20], hi[20], zmean[20];
    int nFst = findZClusters(t, "detType==0", 5.0, lo, hi, zmean);
    printf("nFst=%d\n", nFst);

    Float_t hitX, hitY, hitZ, projX, projY;
    Int_t detType;
    t->SetBranchAddress("detType", &detType);
    t->SetBranchAddress("hitX", &hitX); t->SetBranchAddress("hitY", &hitY); t->SetBranchAddress("hitZ", &hitZ);
    t->SetBranchAddress("projX", &projX); t->SetBranchAddress("projY", &projY);

    double sumDphi[3][12] = {{0}};
    long long cnt[3][12] = {{0}};

    Long64_t nAll = t->GetEntries();
    for (Long64_t i = 0; i < nAll; i++) {
        t->GetEntry(i);
        if (detType != 0) continue;
        int d = -1;
        for (int k = 0; k < nFst && k < 3; k++) { if (hitZ >= lo[k] && hitZ <= hi[k]) { d = k; break; } }
        if (d < 0) continue;
        double r = sqrt(hitX*hitX + hitY*hitY);
        double phi = atan2(hitY, hitX);
        double phiProj = atan2(projY, projX);
        double dphi = phi - phiProj;
        while (dphi > TMath::Pi())  dphi -= 2*TMath::Pi();
        while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
        double rdphi = r * dphi;
        if (fabs(rdphi) > 0.3) continue; // core cut
        int wb = wedgeBinOf(phi * 180.0 / TMath::Pi());
        sumDphi[d][wb] += dphi;
        cnt[d][wb]++;
    }

    printf("\n// Per-(disk, wedge-phi-bin) mean dphi [rad], core-cut sample, %lld total FST rows scanned\n", nAll);
    printf("// wedge-phi-bin 0 covers phi in [-180,-150) deg, increasing by 30 deg per bin\n");
    printf("const double kFstWedgePhiCorrection[3][12] = {\n");
    for (int d = 0; d < 3; d++) {
        printf("    {");
        for (int w = 0; w < 12; w++) {
            double mean = (cnt[d][w] > 0) ? sumDphi[d][w]/cnt[d][w] : 0.0;
            printf("%+.6f%s", mean, (w < 11) ? ", " : "");
        }
        printf("}%s // disk %d\n", (d < 2) ? "," : "", d);
    }
    printf("};\n\n// Diagnostics (N per bin):\n");
    for (int d = 0; d < 3; d++) {
        printf("disk %d: ", d);
        for (int w = 0; w < 12; w++) printf("%lld ", cnt[d][w]);
        printf("\n");
    }
}
