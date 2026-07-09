// testFwdAlignmentTranslation.C
// Direct test of the FST1/FST2 translation signal StFwdResidualMaker's
// (biased) sine-fit found in data/BLCVtx -- A=0.024+-0.005cm (FST1, 4.9sigma),
// A=0.015+-0.003cm (FST2, 5.3sigma) -- using StFwdAlignmentMaker's UNBIASED
// ntuple instead. Same physics (see summarizeFwdResidual.C's fitSine(): a
// rigid (x,y) disk shift decomposes into the azimuthal residual as
// A*sin(phi-phi0)), same fit method (linear cos/sin reparametrization + PDG
// chi2/ndf error scale-factor + epsilon floor against the near-zero-error
// TF1 "N"-option artifact -- see summarizeFwdResidual.C's fitSine() comments
// for why), applied to hit-vs-unbiased-projection instead of hit-vs-biased-
// projection.
//
// Usage: root4star -b -q 'testFwdAlignmentTranslation.C("fwdAlignment_BLCVtx.root")'

#include <vector>
#include <algorithm>

const int kMaxClusters = 20;

// SetEstimate() before Draw() is required once n exceeds ROOT's default
// result-buffer size (historically 1e6): without it, Draw()'s return value
// is the true (uncapped) match count, but GetV1() only backs 'estimate'
// rows -- looping to the full return count then reads/corrupts past the
// buffer (silent heap overrun, e.g. a "*** Break *** segmentation
// violation" that doesn't reproduce at small n). Hit this for real once the
// merged multi-file ntuple pushed the FST-hit count into the tens of
// millions.
int findZClusters(TTree* t, const char* cutExpr, double gapCm, double* lo, double* hi, double* mean) {
    t->SetEstimate(t->GetEntries() + 10);
    Long64_t n = t->Draw("hitZ", cutExpr, "goff");
    if (n <= 0) return 0;
    double* z = t->GetV1();
    double* zs = new double[n];
    for (Long64_t i = 0; i < n; i++) zs[i] = z[i];
    std::sort(zs, zs + n); // was an O(n^2) insertion sort -- fine at ~1e5 rows, not at 1e7+
    int nc = 0; Long64_t start = 0;
    for (Long64_t i = 1; i <= n; i++) {
        if (i == n || (zs[i] - zs[i-1]) > gapCm) {
            if (nc >= kMaxClusters) break;
            lo[nc] = zs[start] - 0.5; hi[nc] = zs[i-1] + 0.5;
            double s = 0; for (Long64_t k = start; k < i; k++) s += zs[k];
            mean[nc] = s / (i - start);
            nc++; start = i;
        }
    }
    delete[] zs;
    return nc;
}

// Same profile-with-nonzero-error-bins approach as summarizeFwdResidual.C's
// profileToGraph(), just fed by an in-memory array instead of a TTree
// expression (we need to compute r*dphi ourselves -- it isn't a stored
// column, unlike in StFwdResidualMaker's histograms).
TGraphErrors* profileToGraphArr(double* phi, double* val, Long64_t n, int nbins, int minBinN) {
    TProfile p("p_tmp", "", nbins, -TMath::Pi(), TMath::Pi());
    p.SetDirectory(0);
    for (Long64_t k = 0; k < n; k++) p.Fill(phi[k], val[k]);
    double vx[400], vy[400], vey[400]; int m = 0;
    for (int i = 1; i <= nbins; i++) {
        if (p.GetBinEntries(i) < minBinN) continue;
        if (p.GetBinError(i) <= 0) continue;
        vx[m] = p.GetBinCenter(i); vy[m] = p.GetBinContent(i); vey[m] = p.GetBinError(i);
        m++;
    }
    if (m < 5) return 0;
    return new TGraphErrors(m, vx, vy, 0, vey);
}

// Same fit + safeguards as summarizeFwdResidual.C::fitSine() -- see that
// file's comments for why each one is there (chi2/ndf scale-factor for
// underestimated per-bin errors; 1e-6 floor against the near-zero-error
// artifact from ROOT's linear-TF1 fast-path).
bool fitSine(double* phi, double* val, Long64_t n, double& amp, double& ampErr, double& chi2ndf, int& nPoints) {
    TGraphErrors* g = profileToGraphArr(phi, val, n, 80, 2);
    if (!g) return false;
    nPoints = g->GetN();
    TF1 fsin("fsin", "[0]*cos(x)+[1]*sin(x)+[2]", -TMath::Pi(), TMath::Pi());
    TFitResultPtr r = g->Fit(&fsin, "Q0S");
    bool ok = r.Get() && r->IsValid() && r->Ndf() > 0;
    if (ok) {
        double a = r->Parameter(0), b = r->Parameter(1);
        double ea = r->ParError(0), eb = r->ParError(1);
        if (ea < 1e-6 || eb < 1e-6) { delete g; return false; }
        chi2ndf = r->Chi2() / r->Ndf();
        double scale = (chi2ndf > 1.0) ? sqrt(chi2ndf) : 1.0;
        double A = sqrt(a*a + b*b);
        amp = A;
        ampErr = (A > 0) ? sqrt(a*a*ea*ea + b*b*eb*eb) / A * scale : 0;
    }
    delete g;
    return ok;
}

// fname: either a single ROOT file, or a glob pattern (e.g.
// "DataDisk/pico/alignment/*.FwdAlignment_BLCVtx.root") naming several
// per-file ntuples. Ntuples don't need hadd'ing the way histograms do --
// unlike a histogram (one object whose bin contents must be summed), an
// ntuple's rows from different files are independent and just need to be
// read in sequence, which is exactly what TChain does without copying any
// data into a new merged file on disk.
void testFwdAlignmentTranslation(const char* fname) {
    TChain* t = new TChain("alignTree");
    int nFiles = t->Add(fname);
    if (nFiles <= 0) { printf("No files matched %s\n", fname); return; }

    double lo[kMaxClusters], hi[kMaxClusters], zmean[kMaxClusters];
    int nFst = findZClusters(t, "detType==0", 5.0, lo, hi, zmean);
    printf("Reading %s: %d files, %lld rows, %d FST disks found\n", fname, nFiles, t->GetEntries(), nFst);

    Float_t hitX, hitY, hitZ, projX, projY, projZ;
    Int_t detType;
    t->SetBranchAddress("detType", &detType);
    t->SetBranchAddress("hitX", &hitX); t->SetBranchAddress("hitY", &hitY); t->SetBranchAddress("hitZ", &hitZ);
    t->SetBranchAddress("projX", &projX); t->SetBranchAddress("projY", &projY); t->SetBranchAddress("projZ", &projZ);
    Long64_t nAll = t->GetEntries();

    // Single pass over the tree, bucketed by disk -- was 3 disks x 2 cut
    // passes = 6 full GetEntry() loops over nAll, fine at ~1e5 rows (the
    // single-file test) but prohibitively slow once nAll reached ~1e8 (the
    // 133-file merged campaign ntuple). The "core cut" pass below just
    // re-filters these same per-disk arrays in memory instead of re-reading
    // the tree.
    std::vector<double> phiAll[kMaxClusters], valAll[kMaxClusters];
    for (int d = 0; d < nFst; d++) {
        phiAll[d].reserve(nAll / (nFst > 0 ? nFst : 1) / 4);
        valAll[d].reserve(nAll / (nFst > 0 ? nFst : 1) / 4);
    }
    for (Long64_t i = 0; i < nAll; i++) {
        t->GetEntry(i);
        if (detType != 0) continue;
        int d = -1;
        for (int k = 0; k < nFst; k++) { if (hitZ >= lo[k] && hitZ <= hi[k]) { d = k; break; } }
        if (d < 0) continue;
        double r = sqrt(hitX*hitX + hitY*hitY);
        double phiHit  = atan2(hitY, hitX);
        double phiProj = atan2(projY, projX);
        double dphi = phiHit - phiProj;
        while (dphi > TMath::Pi())  dphi -= 2*TMath::Pi();
        while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
        phiAll[d].push_back(phiHit);
        valAll[d].push_back(r * dphi);
    }

    // coreCm<=0 means "no cut" (full raw sample); otherwise restrict to
    // |r*dphi|<coreCm before fitting. The unbiased residual has the same
    // tight-core/broad-tail structure documented in plotFwdAlignment.C -- a
    // sine fit to the raw (uncut) sample lets the tail's non-Gaussian scatter
    // distort the per-phi-bin means, which is the likely reason a first,
    // uncut pass gave chi2/ndf~5-6 (vs ~1.2-8 for the well-behaved fits in
    // summarizeFwdResidual.C) and amplitudes that didn't track the biased
    // result consistently disk-to-disk. Report both so the difference itself
    // is visible.
    double coreCms[2] = {-1.0, 0.3};
    const char* labels[2] = {"raw (uncut)", "core only (|r*dphi|<0.3cm)"};
    for (int pass = 0; pass < 2; pass++) {
        double coreCm = coreCms[pass];
        printf("\n--- %s ---\n", labels[pass]);
        printf("%-6s %10s %14s %10s %8s\n", "Disk", "N(fit pts)", "Amplitude [cm]", "chi2/ndf", "sig");
        for (int d = 0; d < nFst; d++) {
            Long64_t nd = phiAll[d].size();
            double* phi = new double[nd];
            double* val = new double[nd];
            Long64_t n = 0;
            for (Long64_t i = 0; i < nd; i++) {
                if (coreCm > 0 && fabs(valAll[d][i]) >= coreCm) continue;
                phi[n] = phiAll[d][i];
                val[n] = valAll[d][i];
                n++;
            }
            double amp = 0, ampErr = 0, chi2ndf = 0; int nPoints = 0;
            bool ok = fitSine(phi, val, n, amp, ampErr, chi2ndf, nPoints);
            if (ok) {
                double sig = (ampErr > 0) ? amp/ampErr : 0;
                printf("FST%d   %10d   %6.4f+-%.4f   %8.2f   %6.1f  (z=%.1f, %lld hits)\n",
                       d, nPoints, amp, ampErr, chi2ndf, sig, zmean[d], n);
            } else {
                printf("FST%d   n/a (not enough usable phi bins, %lld hits)\n", d, n);
            }
            delete[] phi; delete[] val;
        }
    }
    printf("\nCompare to StFwdResidualMaker (biased) data/BLCVtx: "
           "FST1(=disk0) A=0.0241+-0.0049 (4.9sig), FST2(=disk1) A=0.0153+-0.0029 (5.3sig), "
           "FST3(=disk2) A=0.0041+-0.0041 (not significant).\n");
}
