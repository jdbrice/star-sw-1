// compareCorrection15.C
// Full 15-file before/after comparison of the FstWedgeAligner per-wedge phi
// correction (see StRoot/StFwdTrackMaker/StFwdHitLoader.h/.cxx,
// bugreport_StFstHitMaker.txt) against the pre-existing uncorrected baseline,
// for the same run-23081015 files used in the earlier 2- and 4-file checks.
// Corrected ntuples are split across two directories: 11 files were
// reprocessed sequentially on an interactive node, 4 via condor.
//
// Usage: root4star -b -q 'script/compareCorrection15.C'

#include <algorithm>

int findZClusters15(TTree* t, const char* cutExpr, double gapCm, double* lo, double* hi, double* mean) {
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

int wedgeBinOf15(double phiDeg) {
    // See StFwdHitLoader.h's FstWedgeAligner::wedgePhiBin -- no +15 deg offset.
    double shifted = phiDeg + 180.0;
    while (shifted < 0) shifted += 360.0;
    while (shifted >= 360.0) shifted -= 360.0;
    return (int)(shifted / 30.0);
}

void analyze15(TChain* t, const char* label) {
    double lo[20], hi[20], zmean[20];
    int nFst = findZClusters15(t, "detType==0", 5.0, lo, hi, zmean);

    Float_t hitX, hitY, hitZ, projX, projY;
    Int_t detType;
    t->SetBranchAddress("detType", &detType);
    t->SetBranchAddress("hitX", &hitX); t->SetBranchAddress("hitY", &hitY); t->SetBranchAddress("hitZ", &hitZ);
    t->SetBranchAddress("projX", &projX); t->SetBranchAddress("projY", &projY);

    double sumRdphi[3][12] = {{0}};
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
        if (fabs(rdphi) > 0.3) continue;
        int wb = wedgeBinOf15(phi * 180.0 / TMath::Pi());
        sumRdphi[d][wb] += rdphi;
        cnt[d][wb]++;
    }

    printf("\n=== %s (%lld rows) ===\n", label, nAll);
    for (int d = 0; d < 3; d++) {
        printf("disk %d wedge means [cm]: ", d);
        double sumAbs = 0; int nb = 0;
        double vmin=1e9, vmax=-1e9;
        for (int w = 0; w < 12; w++) {
            double m = (cnt[d][w] > 0) ? sumRdphi[d][w]/cnt[d][w] : 0.0;
            printf("%+.4f ", m);
            sumAbs += fabs(m); nb++;
            if (m < vmin) vmin = m;
            if (m > vmax) vmax = m;
        }
        printf("  | mean|level|=%.4f  peak-to-peak=%.4f\n", sumAbs/nb, vmax-vmin);
    }
}

void compareCorrection15() {
    // 11 files reprocessed sequentially, 4 via condor -- same run 23081015
    const char* seqFiles[11] = {
        "6000057","3000049","1000066","5500051","2500048",
        "3500052","6500055","5000047","1500040","3000045","7500055"
    };
    const char* condorFiles[4] = {"2000046","4500029","5000049","7500026"};

    const char* seqDir = "/tmp/claude-2546/-direct-star-u-akio-fcstrk11-star-sw-fwd/999e4fd5-b5d9-4cdc-bbed-3d4d48f7c08f/scratchpad/wedgetest/corrected";
    const char* condorDir = "/gpfs01/star/pwg_tasks/FwdCalib/akio/wedgecorr_condor";

    TChain* tc = new TChain("alignTree");
    for (int i = 0; i < 11; i++)
        tc->Add(Form("%s/st_fwd_23081015_raw_%s.FwdAlignment_BLCVtx.root", seqDir, seqFiles[i]));
    for (int i = 0; i < 4; i++)
        tc->Add(Form("%s/st_fwd_23081015_raw_%s.FwdAlignment_BLCVtx.root", condorDir, condorFiles[i]));
    analyze15(tc, "CORRECTED (15 files)");

    TChain* tu = new TChain("alignTree");
    for (int i = 0; i < 11; i++)
        tu->Add(Form("DataDisk/pico/alignment/st_fwd_23081015_raw_%s.FwdAlignment_BLCVtx.root", seqFiles[i]));
    for (int i = 0; i < 4; i++)
        tu->Add(Form("DataDisk/pico/alignment/st_fwd_23081015_raw_%s.FwdAlignment_BLCVtx.root", condorFiles[i]));
    analyze15(tu, "UNCORRECTED (same 15 files)");
}
