// plotFwdAlignment.C
// Plots unbiased (hit-removed) FST/FTT residuals from StFwdAlignmentMaker's
// output ntuple (see proposal_alignment_path.txt and StFwdAlignmentMaker.*).
//
// Hits are grouped into FST disks / FTT planes by z-clustering rather than by
// genfitPlaneIndex: that index was found not to survive into the leave-one-out
// seed reliably (see the NOTE in StFwdAlignmentMaker::processTrack), so hitZ
// -- a physical, directly-trustworthy hit coordinate -- is used instead. FTT's
// x-strip/y-strip sublayers (~1-2cm apart in z) land in the same z-cluster on
// purpose; they're separated by the stripDir branch instead, not by z.
//
// FST is a polar (r,phi) sensor: the physically meaningful residual is
// r*dphi (see analyzeFstDphi), not independent Cartesian dx=hitX-projX/
// dy=hitY-projY -- those aren't separate measurements of anything for a
// disk sensor read out in r and phi. FTT hits are 1D (a single strip
// measures only x or y -- see stripDir), so only the measured axis is shown
// per stripDir (V-strip -> dx, H-strip -> dy); the other axis isn't a real
// measurement and is skipped rather than plotted as if it were one.
//
// CORE/TAIL: the raw residual distribution is not a single Gaussian -- it's a
// tight core (a real, well-constrained refit) sitting on a broad pedestal
// (poorly-constrained refits, most often when the removed hit isn't bracketed
// by the hits that remain, so the fit has to extrapolate rather than
// interpolate through that z). A plain RMS conflates the two and overstates
// the "typical" residual by a lot -- e.g. FST outer disks were 40% core /
// 60% tail in a first look, with the tail alone pulling the full RMS to
// ~15-20x the biased (StFwdResidualMaker) value. So for every plane/axis this
// macro reports BOTH the full RMS (for reference) AND a core/tail breakdown:
// core fraction (with a binomial error) and a Gaussian sigma fit *restricted
// to* |residual|<coreCm (not a truncated-sample RMS, which is itself biased
// low by the truncation -- an actual fit to the core window is the
// statistically correct way to characterize its width).
//
// Usage: root4star -b -q 'plotFwdAlignment.C("fwdAlignment.root","alignment/data")'
// fname may be a glob pattern (e.g. "DataDisk/pico/alignment/*.FwdAlignment_BLCVtx.root")
// naming several per-file ntuples -- loaded via TChain rather than hadd'd
// into one file first: ntuple rows from different files are independent and
// just need reading in sequence, unlike a histogram whose bins must be summed.

#include <algorithm>

const int kMaxClusters = 20;
const double kCoreCm = 0.3; // core/tail split point, cm -- see file header

// Sorts all hitZ values passing cutExpr, splits into clusters wherever a gap
// exceeds gapCm. Returns cluster count; fills lo[]/hi[]/mean[] (caller-owned,
// sized kMaxClusters). No STL containers in the signature -- keeps this
// CINT-safe like the rest of this session's macros.
int findZClusters(TTree* t, const char* cutExpr, double gapCm, double* lo, double* hi, double* mean) {
    // SetEstimate() before Draw() is required once n exceeds ROOT's default
    // result-buffer size (historically 1e6): without it, Draw()'s return
    // value is the true (uncapped) match count, but GetV1() only backs
    // 'estimate' rows -- looping to the full return count then reads past
    // the buffer (silent heap overrun/segfault). Only shows up once fed a
    // multi-file chain instead of the single small test file this was
    // written against.
    t->SetEstimate(t->GetEntries() + 10);
    Long64_t n = t->Draw("hitZ", cutExpr, "goff");
    if (n <= 0) return 0;
    double* z = t->GetV1();
    double* zs = new double[n];
    for (Long64_t i = 0; i < n; i++) zs[i] = z[i];
    std::sort(zs, zs + n); // was an O(n^2) insertion sort -- fine at ~1e3-1e4 rows, not at 1e7+

    int nc = 0;
    Long64_t start = 0;
    for (Long64_t i = 1; i <= n; i++) {
        if (i == n || (zs[i] - zs[i-1]) > gapCm) {
            if (nc >= kMaxClusters) break;
            lo[nc] = zs[start] - 0.5;
            hi[nc] = zs[i-1] + 0.5;
            double s = 0;
            for (Long64_t k = start; k < i; k++) s += zs[k];
            mean[nc] = s / (i - start);
            nc++;
            start = i;
        }
    }
    delete[] zs;
    return nc;
}

// Pulls varExpr (already restricted by cutExpr) into memory once and computes
// everything downstream from that single copy: full N/mean/RMS, core
// fraction (+binomial error), and a Gaussian sigma fit restricted to
// |v|<coreCm (+fit error). Returns false (all outputs zeroed) if there's no
// data -- callers should treat that as "n/a", not "zero residual".
bool coreTailStats(TTree* t, const char* varExpr, const char* cutExpr, double coreCm,
                    Long64_t& n, double& mean, double& rms,
                    Long64_t& nCore, double& coreFrac, double& coreFracErr,
                    double& coreSigma, double& coreSigmaErr) {
    n = 0; mean = 0; rms = 0; nCore = 0; coreFrac = 0; coreFracErr = 0;
    coreSigma = 0; coreSigmaErr = 0;
    t->SetEstimate(t->GetEntries() + 10); // see findZClusters -- same Draw()/GetV1() buffer-cap issue
    Long64_t nn = t->Draw(varExpr, cutExpr, "goff");
    if (nn <= 0) return false;
    n = nn;
    double* v = t->GetV1();
    double s = 0, s2 = 0;
    static int callCount = 0;
    callCount++;
    TH1F hcore(Form("h_coretail_tmp_%d", callCount), "", 100, -coreCm, coreCm);
    hcore.SetDirectory(0);
    for (Long64_t k = 0; k < n; k++) {
        s += v[k]; s2 += v[k]*v[k];
        if (fabs(v[k]) < coreCm) { nCore++; hcore.Fill(v[k]); }
    }
    mean = s / n;
    rms = sqrt(s2 / n - mean*mean);
    coreFrac = (double)nCore / n;
    coreFracErr = sqrt(coreFrac * (1.0 - coreFrac) / n); // binomial error on a fraction
    if (nCore > 20) {
        // "0S", not "N": "N" (as used elsewhere in earlier iterations of this
        // session's macros) suppresses storing the fit function on the
        // histogram, which then makes GetFunction() always return null --
        // silently zeroing every result. Capture via TFitResultPtr instead.
        TFitResultPtr r = hcore.Fit("gaus", "Q0S");
        if (r.Get() && r->IsValid()) { coreSigma = r->Parameter(2); coreSigmaErr = r->ParError(2); }
    }
    return true;
}

// Shared final-drawing step for plot1D() (tree-Draw-filled) and the FST
// r*dphi path (manually-filled, see analyzeFstDphi) -- same core-fit overlay
// and stats box either way.
void plot1DFromHist(TH1F* h, const char* png, const char* xtitle, double coreCm) {
    TCanvas* c = new TCanvas(Form("c_%s", png), "", 500, 400);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    gStyle->SetOptStat(1110);
    h->SetStats(1);
    h->SetLineColor(kBlue+1); h->SetLineWidth(2); h->SetFillColor(kBlue-9);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle("Hits");
    h->Draw("hist");
    // Overlay the core-only Gaussian (green) -- the fit that actually means
    // something for this two-component shape -- rather than a full-range fit,
    // which would just be dragged wide by the tail.
    if (h->GetEntries() > 20) {
        TF1 coreFit(Form("%s_core", h->GetName()), "gaus", -coreCm, coreCm);
        h->Fit(&coreFit, "QN0R");
        coreFit.SetLineColor(kGreen+2); coreFit.SetLineWidth(2);
        coreFit.DrawCopy("same"); // DrawCopy: coreFit is a stack local, out of scope once this returns
    }
    c->Modified(); c->Update();
    c->Print(png);
    delete c;
}

void plot1D(TTree* t, const char* varExpr, const char* cutExpr, const char* hname, const char* png,
            const char* xtitle, double rangeCm, double coreCm) {
    TH1F* h = new TH1F(hname, "", 100, -rangeCm, rangeCm);
    t->Draw(Form("%s>>%s", varExpr, hname), cutExpr, "goff");
    plot1DFromHist(h, png, xtitle, coreCm);
}

// Same stats as coreTailStats() but from an already-computed array (used by
// the FST r*dphi path, which needs hitX/Y and projX/Y combined with proper
// phi-wrap handling that TTree::Draw's string-expression parser can't do
// cleanly -- see analyzeFstDphi).
bool coreTailStatsFromArray(double* v, Long64_t n, double coreCm,
                             double& mean, double& rms,
                             Long64_t& nCore, double& coreFrac, double& coreFracErr,
                             double& coreSigma, double& coreSigmaErr) {
    mean = 0; rms = 0; nCore = 0; coreFrac = 0; coreFracErr = 0;
    coreSigma = 0; coreSigmaErr = 0;
    if (n <= 0) return false;
    double s = 0, s2 = 0;
    static int callCount = 0;
    callCount++;
    TH1F hcore(Form("h_coretail_arr_%d", callCount), "", 100, -coreCm, coreCm);
    hcore.SetDirectory(0);
    for (Long64_t k = 0; k < n; k++) {
        s += v[k]; s2 += v[k]*v[k];
        if (fabs(v[k]) < coreCm) { nCore++; hcore.Fill(v[k]); }
    }
    mean = s / n;
    rms = sqrt(s2 / n - mean*mean);
    coreFrac = (double)nCore / n;
    coreFracErr = sqrt(coreFrac * (1.0 - coreFrac) / n);
    if (nCore > 20) {
        TFitResultPtr r = hcore.Fit("gaus", "Q0S");
        if (r.Get() && r->IsValid()) { coreSigma = r->Parameter(2); coreSigmaErr = r->ParError(2); }
    }
    return true;
}

// FST-specific counterpart to analyzeAxis(): FST is a polar (r,phi) sensor,
// so the physically meaningful residual is r*dphi (same convention as
// plotFwdAlignmentResidual.C and the rest of this study), not independent
// Cartesian dx=hitX-projX/dy=hitY-projY -- those two aren't separate
// measurements of anything for FST the way they are for FTT's Cartesian
// strips. Computed manually (not via a TTree::Draw string expression) so
// the phi difference can be wrapped to (-pi,pi] correctly -- hitX/Y and
// projX/Y are pulled together via one 4-variable Draw, since GetV1..GetV4
// give parallel arrays for a single Draw() call.
void analyzeFstDphi(TTree* t, FILE* fpRaw, FILE* fpCore, const char* planeLabel, double zmean,
                     const char* cutExpr, const char* outdir, double rangeCm) {
    TString planeLower(planeLabel);
    planeLower.ToLower();
    TString png = Form("%s/%sAlign_rdphi.png", outdir, planeLower.Data());
    TString hname = Form("h_%s_rdphi", planeLower.Data());

    t->SetEstimate(t->GetEntries() + 10);
    Long64_t n = t->Draw("hitX:hitY:projX:projY", cutExpr, "goff");
    if (n <= 0) {
        fprintf(fpRaw, "<tr><td>%s</td><td>%.1f</td><td>rdphi</td><td colspan=3>n/a (no data)</td></tr>\n",
                planeLabel, zmean);
        fprintf(fpCore, "<tr><td>%s</td><td>rdphi</td><td colspan=4>n/a (no data)</td></tr>\n", planeLabel);
        return;
    }
    double* hitX = t->GetV1(); double* hitY = t->GetV2();
    double* projX = t->GetV3(); double* projY = t->GetV4();
    double* rdphi = new double[n];
    for (Long64_t k = 0; k < n; k++) {
        double r = sqrt(hitX[k]*hitX[k] + hitY[k]*hitY[k]);
        double dphi = atan2(hitY[k], hitX[k]) - atan2(projY[k], projX[k]);
        while (dphi > TMath::Pi())  dphi -= 2*TMath::Pi();
        while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
        rdphi[k] = r * dphi;
    }

    TH1F* h = new TH1F(hname, "", 100, -rangeCm, rangeCm);
    h->SetDirectory(0);
    for (Long64_t k = 0; k < n; k++) h->Fill(rdphi[k]);
    plot1DFromHist(h, png, "r#upoint#Delta#phi [cm]", kCoreCm);

    double mean, rms, coreFrac, coreFracErr, coreSigma, coreSigmaErr;
    Long64_t nCore;
    bool ok = coreTailStatsFromArray(rdphi, n, kCoreCm, mean, rms,
                                      nCore, coreFrac, coreFracErr, coreSigma, coreSigmaErr);
    delete[] rdphi;

    fprintf(fpRaw, "<tr><td>%s</td><td>%.1f</td><td>rdphi</td><td>%lld</td>"
                   "<td>%.4f</td><td>%.4f</td><td><a href=\"%s\">plot</a></td></tr>\n",
            planeLabel, zmean, n, mean, rms, gSystem->BaseName(png.Data()));

    if (ok) {
        TString sigmaStr;
        if (nCore > 20) sigmaStr = Form("%.4f&plusmn;%.4f", coreSigma, coreSigmaErr);
        else sigmaStr = "n/a";
        fprintf(fpCore, "<tr><td>%s</td><td>rdphi</td><td>%lld</td>"
                        "<td>%.1f&plusmn;%.1f%%</td><td>%s</td>"
                        "<td>%.1f%%</td></tr>\n",
                planeLabel, n, coreFrac*100, coreFracErr*100,
                sigmaStr.Data(), (1.0-coreFrac)*100);
    }
}

// Writes one row of the raw-stats table plus one row of the core/tail table.
// Returns via out params so the caller can also print a running summary.
void analyzeAxis(TTree* t, FILE* fpRaw, FILE* fpCore, const char* planeLabel, double zmean,
                  const char* axisName, const char* varExpr, const char* cutExpr,
                  const char* outdir, double rangeCm) {
    // planeLabel (e.g. "FST0", "FTT2") is unique per call -- use it (lowercased)
    // for plot/histogram names, so they don't collide across planes.
    TString planeLower(planeLabel);
    planeLower.ToLower();
    TString png = Form("%s/%sAlign_%s.png", outdir, planeLower.Data(), axisName);
    TString hname = Form("h_%s_%s", planeLower.Data(), axisName);
    plot1D(t, varExpr, cutExpr, hname, png, Form("#Delta%s [cm]", axisName+1), rangeCm, kCoreCm);

    Long64_t n, nCore;
    double mean, rms, coreFrac, coreFracErr, coreSigma, coreSigmaErr;
    bool ok = coreTailStats(t, varExpr, cutExpr, kCoreCm, n, mean, rms,
                             nCore, coreFrac, coreFracErr, coreSigma, coreSigmaErr);

    fprintf(fpRaw, "<tr><td>%s</td><td>%.1f</td><td>%s</td><td>%lld</td>"
                   "<td>%.4f</td><td>%.4f</td><td><a href=\"%s\">plot</a></td></tr>\n",
            planeLabel, zmean, axisName, n, mean, rms, gSystem->BaseName(png.Data()));

    if (ok) {
        // coreSigma is only fit (see coreTailStats) when nCore>20; below that,
        // report n/a rather than a misleading 0.0000+-0.0000.
        TString sigmaStr;
        if (nCore > 20) sigmaStr = Form("%.4f&plusmn;%.4f", coreSigma, coreSigmaErr);
        else sigmaStr = "n/a";
        fprintf(fpCore, "<tr><td>%s</td><td>%s</td><td>%lld</td>"
                        "<td>%.1f&plusmn;%.1f%%</td><td>%s</td>"
                        "<td>%.1f%%</td></tr>\n",
                planeLabel, axisName, n, coreFrac*100, coreFracErr*100,
                sigmaStr.Data(), (1.0-coreFrac)*100);
    } else {
        fprintf(fpCore, "<tr><td>%s</td><td>%s</td><td colspan=4>n/a (no data)</td></tr>\n",
                planeLabel, axisName);
    }
}

// trackTypeLabel: printed on the page so it's not left implicit in fname's glob
// pattern (e.g. "BLCVtx") -- pass "" to omit.
void plotFwdAlignment(const char* fname = "fwdAlignment.root", const char* outdir = "alignment",
                       const char* trackTypeLabel = "") {
    TChain* t = new TChain("alignTree");
    int nFiles = t->Add(fname);
    if (nFiles <= 0) { printf("No files matched %s\n", fname); return; }
    gSystem->mkdir(outdir, kTRUE);
    gErrorIgnoreLevel = kWarning;

    printf("Reading %s (%d files, %lld rows) -> %s/\n", fname, nFiles, t->GetEntries(), outdir);

    FILE* fp = fopen(Form("%s/summary.html", outdir), "w");
    fprintf(fp, "\n<HTML>\n<HEAD>\n<TITLE> OGAWA Akio Home page </TITLE>\n</HEAD>\n");
    fprintf(fp, "<BODY BGCOLOR=\"#dfdfff\" TEXT=\"black\" LINK=\"blue\" VLINK=\"darkblue\">\n");
    fprintf(fp, "<H1>Unbiased (hit-removed) FST/FTT Residuals</H1><HR>\n");
    fprintf(fp, "<p>Source: <code>%s</code>%s%s, %d files, %lld rows. Auto-generated by "
                "<code>script/plotFwdAlignment.C</code>. See "
                "<code>proposal_alignment_path.txt</code> and "
                "<code>StFwdAlignmentMaker.cxx</code> for the method: each row is one "
                "FST/FTT hit refit out of its track (leave-one-out), so unlike "
                "StFwdResidualMaker's plots this residual is NOT pulled toward the hit "
                "by the fit that produced it.</p>\n", fname,
                (trackTypeLabel[0] ? " (track type: " : ""),
                (trackTypeLabel[0] ? Form("%s)", trackTypeLabel) : ""),
                nFiles, t->GetEntries());

    fprintf(fp, "<H3>Core/Tail Structure</H3>\n");
    fprintf(fp, "<p>The raw residual is not a single Gaussian: a tight core (a well-constrained "
                "refit) sits on a broad pedestal (a poorly-constrained one -- typically when the "
                "removed hit isn't bracketed in z by the hits that remain, so the refit has to "
                "extrapolate rather than interpolate through that point). The plain RMS in the "
                "table below conflates both and overstates the \"typical\" residual. Core "
                "fraction uses a binomial error; core &sigma; is a Gaussian fit restricted to "
                "|&Delta;|&lt;%.1f&nbsp;cm (not a truncated-sample RMS, which is itself biased "
                "low by the cut) -- these two together are a better summary than the raw RMS "
                "alone.</p>\n", kCoreCm);
    fprintf(fp, "<table border=1 cellpadding=4>\n<tr><th>Plane</th><th>axis</th><th>N</th>"
                "<th>core fraction</th><th>core &sigma; [cm]</th><th>tail fraction</th></tr>\n");

    FILE* fpCore = fp; // write into the same file, just a different table -- see below for split
    // (kept as a second FILE* parameter for analyzeAxis() so both tables can be filled from one pass)

    // Raw-stats table is written to a temp buffer HTML fragment file first, then appended after
    // the core/tail table -- avoids needing two passes over the tree.
    TString rawPath = Form("%s/.rawtable.html", outdir);
    FILE* fpRaw = fopen(rawPath.Data(), "w");
    fprintf(fpRaw, "<table border=1 cellpadding=4>\n<tr><th>Plane</th><th>mean z [cm]</th>"
                   "<th>axis</th><th>N</th><th>mean [cm]</th><th>RMS [cm]</th><th>plot</th></tr>\n");

    // ---- FST: polar sensor, r*dphi is the physically meaningful residual (not dx/dy) ----
    double lo[kMaxClusters], hi[kMaxClusters], mean[kMaxClusters];
    int nFst = findZClusters(t, "detType==0", 5.0, lo, hi, mean);
    printf("FST: %d z-clusters found\n", nFst);
    for (int i = 0; i < nFst; i++) {
        TString cut = Form("detType==0 && hitZ>=%f && hitZ<=%f", lo[i], hi[i]);
        analyzeFstDphi(t, fpRaw, fpCore, Form("FST%d", i), mean[i], cut, outdir, 0.5);
    }

    // ---- FTT: 1D strips, only the measured axis per stripDir ----
    int nFtt = findZClusters(t, "detType==1", 10.0, lo, hi, mean);
    printf("FTT: %d z-clusters found\n", nFtt);
    for (int i = 0; i < nFtt; i++) {
        // V-strip (stripDir==1) measures x -> dx; H-strip (stripDir==2) measures y -> dy
        TString cutX = Form("detType==1 && stripDir==1 && hitZ>=%f && hitZ<=%f", lo[i], hi[i]);
        TString cutY = Form("detType==1 && stripDir==2 && hitZ>=%f && hitZ<=%f", lo[i], hi[i]);
        analyzeAxis(t, fpRaw, fpCore, Form("FTT%d", i), mean[i], "dx", "hitX-projX", cutX, outdir, 2.0);
        analyzeAxis(t, fpRaw, fpCore, Form("FTT%d", i), mean[i], "dy", "hitY-projY", cutY, outdir, 2.0);
    }

    fprintf(fp, "</table>\n");

    fclose(fpRaw);
    fpRaw = fopen(rawPath.Data(), "r");
    fprintf(fp, "<H3>Raw Statistics (full-range mean/RMS)</H3>\n");
    fprintf(fp, "<p>Included for reference/continuity -- the core/tail table above is the more "
                "meaningful summary for this two-component shape.</p>\n");
    char buf[4096];
    while (fgets(buf, sizeof(buf), fpRaw)) fputs(buf, fp);
    fprintf(fp, "</table>\n");
    fclose(fpRaw);
    gSystem->Unlink(rawPath.Data());

    fprintf(fp, "<HR>OGAWA, Akio<ADDRESS><A HREF=\"mailto://akio@bnl.gov\">akio@bnl.gov</A></ADDRESS>\n</BODY></HTML>\n");
    fclose(fp);
    printf("Wrote %s/summary.html\n", outdir);
    printf("Done.\n");
}
