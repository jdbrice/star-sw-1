// plotFwdAlignmentResidual.C
// Unbiased (hit-removed) counterpart to StFwdResidualMaker's biased FST
// r*dPhi-vs-x/y/r/r*phi/phi plots (see script/plotFwdResidual.C and
// StFwdResidualMaker.cxx's bookHistos()/fillFst()) -- same quantities, same
// binning, built from StFwdAlignmentMaker's leave-one-out ntuple instead, so
// it isn't pulled toward the hit by the fit that produced it (see
// script/plotFwdAlignment.C for the core/tail 1D summary already on this
// page; this macro adds the 2D distributions and the sine-fit "translation
// test" plots, both previously text/table-only).
//
// Single TChain pass fills every histogram directly (TH1::Fill/TH2::Fill),
// no per-row array buffering and no repeated TTree::Draw() calls -- that's
// what made the previous multi-Draw approach take over an hour on the full
// 139-file dataset. Two short Draw()-based passes (findZClusters) run first
// just to locate the FST/FTT z bands.
//
// Usage: root4star -b -q 'plotFwdAlignmentResidual.C("DataDisk/pico/alignment/*.FwdAlignment_BLCVtx.root","residual/data/alignment","BLCVtx")'

#include <algorithm>

const int kMaxClusters = 20;
const double kCoreCm = 0.3; // core/tail split point, cm -- see plotFwdAlignment.C

int findZClusters(TTree* t, const char* cutExpr, double gapCm, double* lo, double* hi, double* mean) {
    t->SetEstimate(t->GetEntries() + 10); // see plotFwdAlignment.C -- Draw()/GetV1() buffer-cap gotcha
    Long64_t n = t->Draw("hitZ", cutExpr, "goff");
    if (n <= 0) return 0;
    double* z = t->GetV1();
    double* zs = new double[n];
    for (Long64_t i = 0; i < n; i++) zs[i] = z[i];
    std::sort(zs, zs + n);
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

void plot1D(TH1F* h, const char* png, const char* xtitle) {
    TCanvas* c = new TCanvas(Form("c_%s", png), "", 500, 400);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    gStyle->SetOptStat(1110);
    h->SetStats(1); // belt-and-suspenders: gStyle alone sometimes left the stats
                     // box off in this macro's actual multi-canvas run (not
                     // reproduced in isolation) -- force it on the histogram too.
    h->SetLineColor(kBlue+1); h->SetLineWidth(2); h->SetFillColor(kBlue-9);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle("Hits");
    h->Draw("hist");
    if (h->GetEntries() > 20) {
        h->Fit("gaus", "QN");
        TF1* fit = h->GetFunction("gaus");
        if (fit) { fit->SetLineColor(kRed); fit->SetLineWidth(2); fit->Draw("same"); }
    }
    c->Modified(); c->Update();
    c->Print(png);
    delete c;
}

// Same style as script/plotFwdResidual.C's plot2D(): colz map + a linear fit
// to the per-x-bin mean (TProfile), so a tilted band (rotation-type
// misalignment) shows up as a nonzero slope even though this macro's main
// purpose is the translation (sine-vs-phi) test.
void plot2D(TH2F* h, const char* png, const char* xtitle, const char* ytitle) {
    TCanvas* c = new TCanvas(Form("c_%s", png), "", 500, 400);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.14);
    gStyle->SetOptStat(0);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->Draw("colz");
    TProfile* prof = 0;
    if (h->GetEntries() > 50) {
        prof = h->ProfileX(Form("%s_pfx", h->GetName()));
        prof->SetDirectory(0);
        TFitResultPtr r = prof->Fit("pol1", "Q0S");
        if (r.Get() && r->IsValid()) {
            double slope = r->Parameter(1), slopeErr = r->ParError(1);
            TF1* fit = prof->GetFunction("pol1");
            if (fit) { fit->SetLineColor(kRed); fit->SetLineWidth(3); fit->Draw("same"); }
            prof->SetMarkerStyle(24); prof->SetMarkerSize(1.0); prof->SetMarkerColor(kBlack);
            prof->SetLineColor(kBlack); prof->SetLineWidth(1);
            prof->Draw("same p");
            double sig = (slopeErr > 0) ? fabs(slope)/slopeErr : 0;
            TLatex lat; lat.SetNDC(); lat.SetTextSize(0.042);
            lat.SetTextColor(sig > 3 ? kRed+1 : kBlack);
            lat.DrawLatex(0.15, 0.85, Form("slope = %.4f #pm %.4f", slope, slopeErr));
        }
    }
    c->Print(png);
    if (prof) delete prof; // live primitive until Print() rasterizes it -- see plotFwdResidual.C
    delete c;
}

// Bare colz plot: no profile points, no fit line, no stats box -- for
// presentation use where the fit overlay would be a distraction.
void plot2DBare(TH2F* h, const char* png, const char* xtitle, const char* ytitle) {
    TCanvas* c = new TCanvas(Form("c_%s", png), "", 500, 400);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.14);
    gStyle->SetOptStat(0);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->Draw("colz");
    c->Print(png);
    delete c;
}

// Translation-test plot: profile the rdphi-vs-phi 2D histogram, fit
// A*cos(phi)+B*sin(phi)+C (same linear reparametrization as
// testFwdAlignmentTranslation.C/summarizeFwdResidual.C's fitSine()), draw the
// profile + fit curve, and report amplitude/chi2ndf/significance in the
// corner. Returns false (plot still written, just without a fit overlay) if
// the fit doesn't converge or has a near-zero-error artifact.
bool plotSine(TH2F* h2phi, const char* png, const char* title,
              double& amp, double& ampErr, double& chi2ndf, int& nPoints) {
    amp = ampErr = chi2ndf = 0; nPoints = 0;
    TProfile* p = h2phi->ProfileX(Form("%s_pfx", h2phi->GetName()));
    p->SetDirectory(0);
    int nbins = p->GetNbinsX();
    double vx[400], vy[400], vey[400]; int m = 0;
    for (int i = 1; i <= nbins; i++) {
        if (p->GetBinEntries(i) < 2) continue;
        if (p->GetBinError(i) <= 0) continue;
        vx[m] = p->GetBinCenter(i); vy[m] = p->GetBinContent(i); vey[m] = p->GetBinError(i); m++;
    }
    bool ok = false;
    TF1 fsin("fsin", "[0]*cos(x)+[1]*sin(x)+[2]", -TMath::Pi(), TMath::Pi());
    TGraphErrors* g = 0;
    if (m >= 5) {
        g = new TGraphErrors(m, vx, vy, 0, vey);
        nPoints = m;
        TFitResultPtr r = g->Fit(&fsin, "Q0S");
        ok = r.Get() && r->IsValid() && r->Ndf() > 0;
        if (ok) {
            double a = r->Parameter(0), b = r->Parameter(1);
            double ea = r->ParError(0), eb = r->ParError(1);
            if (ea < 1e-6 || eb < 1e-6) ok = false;
            else {
                chi2ndf = r->Chi2() / r->Ndf();
                double scale = (chi2ndf > 1.0) ? sqrt(chi2ndf) : 1.0;
                double A = sqrt(a*a + b*b);
                amp = A; ampErr = sqrt(a*a*ea*ea + b*b*eb*eb) / A * scale;
            }
        }
    }
    TCanvas* c = new TCanvas(Form("c_%s", png), "", 500, 400);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    gStyle->SetOptStat(0);
    p->SetTitle(title);
    p->GetXaxis()->SetTitle("#phi_{hit} [rad]");
    p->GetYaxis()->SetTitle("r#upoint#Delta#phi [cm]");
    p->GetYaxis()->SetRangeUser(-0.15, 0.15);
    p->SetMarkerStyle(20); p->SetMarkerSize(0.8); p->SetMarkerColor(kBlue+1);
    p->SetLineColor(kBlue+1);
    p->Draw();
    if (ok) {
        fsin.SetLineColor(kRed); fsin.SetLineWidth(2);
        fsin.DrawCopy("same");
        double sig = (ampErr > 0) ? amp/ampErr : 0;
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.040);
        lat.SetTextColor(sig > 3 ? kRed+1 : kBlack);
        lat.DrawLatex(0.15, 0.85, Form("A = %.4f #pm %.4f cm (%.1f#sigma)", amp, ampErr, sig));
        lat.SetTextColor(kBlack);
        lat.DrawLatex(0.15, 0.79, Form("#chi^{2}/ndf = %.1f, N=%d", chi2ndf, m));
    }
    c->Print(png);
    delete c;
    if (g) delete g;
    delete p;
    return ok;
}

void plotFwdAlignmentResidual(const char* fname = "DataDisk/pico/alignment/*.FwdAlignment_BLCVtx.root",
                               const char* outdir = "residual/data/alignment",
                               const char* trackTypeLabel = "BLCVtx") {
    TChain* t = new TChain("alignTree");
    int nFiles = t->Add(fname);
    if (nFiles <= 0) { printf("No files matched %s\n", fname); return; }
    gSystem->mkdir(outdir, kTRUE);
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);

    Long64_t nAll = t->GetEntries();
    printf("Reading %s (%d files, %lld rows)\n", fname, nFiles, nAll);

    double fstLo[kMaxClusters], fstHi[kMaxClusters], fstZ[kMaxClusters];
    int nFst = findZClusters(t, "detType==0", 5.0, fstLo, fstHi, fstZ);
    double fttLo[kMaxClusters], fttHi[kMaxClusters], fttZ[kMaxClusters];
    int nFtt = findZClusters(t, "detType==1", 10.0, fttLo, fttHi, fttZ);
    printf("FST disks=%d FTT planes=%d\n", nFst, nFtt);

    const double FST_RPOS = 30.0; // FST hits sit at r<28cm -- +-65 (FTT's range) left the plot mostly empty
    const double FTT_RPOS = 65.0;
    const int NR = 100, NP = 80;
    const double FST_RRES = 0.5, FTT_RRES = 1.0; // match the ranges already used on this page's 1D plots

    // ---- book histograms: FST rdphi (new, mirrors StFwdResidualMaker exactly) ----
    TH1F* hFstRdphi[kMaxClusters];
    TH2F* h2FstX[kMaxClusters]; TH2F* h2FstY[kMaxClusters]; TH2F* h2FstR[kMaxClusters];
    TH2F* h2FstRphi[kMaxClusters]; TH2F* h2FstPhi[kMaxClusters]; TH2F* h2FstPhiCore[kMaxClusters];
    for (int d = 0; d < nFst; d++) {
        hFstRdphi[d] = new TH1F(Form("h_fstU_d%d_rdphi", d), "", NR, -FST_RRES, FST_RRES);
        hFstRdphi[d]->SetDirectory(0);
        h2FstX[d] = new TH2F(Form("h2_fstU_d%d_vs_x", d), "", NP, -FST_RPOS, FST_RPOS, NR, -FST_RRES, FST_RRES);
        h2FstY[d] = new TH2F(Form("h2_fstU_d%d_vs_y", d), "", NP, -FST_RPOS, FST_RPOS, NR, -FST_RRES, FST_RRES);
        h2FstR[d] = new TH2F(Form("h2_fstU_d%d_vs_r", d), "", 60, 0, FST_RPOS, NR, -FST_RRES, FST_RRES);
        h2FstRphi[d] = new TH2F(Form("h2_fstU_d%d_vs_rphi", d), "", NP, -FST_RPOS, FST_RPOS, NR, -FST_RRES, FST_RRES);
        h2FstPhi[d] = new TH2F(Form("h2_fstU_d%d_vs_phi", d), "", NP, -TMath::Pi(), TMath::Pi(), NR, -FST_RRES, FST_RRES);
        h2FstPhiCore[d] = new TH2F(Form("h2_fstU_d%d_vs_phi_core", d), "", NP, -TMath::Pi(), TMath::Pi(), NR, -kCoreCm, kCoreCm);
        h2FstX[d]->SetDirectory(0); h2FstY[d]->SetDirectory(0); h2FstR[d]->SetDirectory(0);
        h2FstRphi[d]->SetDirectory(0); h2FstPhi[d]->SetDirectory(0); h2FstPhiCore[d]->SetDirectory(0);
    }

    // ---- book histograms: FTT dx/dy (existing quantity, new 2D distributions) ----
    TH1F* hFttRes[kMaxClusters][2];
    TH2F* h2FttX[kMaxClusters][2]; TH2F* h2FttY[kMaxClusters][2];
    TH2F* h2FttR[kMaxClusters][2]; TH2F* h2FttRphi[kMaxClusters][2]; TH2F* h2FttPhi[kMaxClusters][2];
    const char* axName[2] = {"dx", "dy"};
    for (int p = 0; p < nFtt; p++) {
        for (int a = 0; a < 2; a++) {
            hFttRes[p][a]   = new TH1F(Form("h_fttU_p%d_%s", p, axName[a]), "", NR, -FTT_RRES, FTT_RRES);
            h2FttX[p][a]    = new TH2F(Form("h2_fttU_p%d_%s_vs_x", p, axName[a]), "", NP, -FTT_RPOS, FTT_RPOS, NR, -FTT_RRES, FTT_RRES);
            h2FttY[p][a]    = new TH2F(Form("h2_fttU_p%d_%s_vs_y", p, axName[a]), "", NP, -FTT_RPOS, FTT_RPOS, NR, -FTT_RRES, FTT_RRES);
            h2FttR[p][a]    = new TH2F(Form("h2_fttU_p%d_%s_vs_r", p, axName[a]), "", 60, 0, FTT_RPOS, NR, -FTT_RRES, FTT_RRES);
            h2FttRphi[p][a] = new TH2F(Form("h2_fttU_p%d_%s_vs_rphi", p, axName[a]), "", NP, -FTT_RPOS, FTT_RPOS, NR, -FTT_RRES, FTT_RRES);
            h2FttPhi[p][a]  = new TH2F(Form("h2_fttU_p%d_%s_vs_phi", p, axName[a]), "", NP, -TMath::Pi(), TMath::Pi(), NR, -FTT_RRES, FTT_RRES);
            hFttRes[p][a]->SetDirectory(0);
            h2FttX[p][a]->SetDirectory(0); h2FttY[p][a]->SetDirectory(0); h2FttR[p][a]->SetDirectory(0);
            h2FttRphi[p][a]->SetDirectory(0); h2FttPhi[p][a]->SetDirectory(0);
        }
    }

    // ---- single pass: fill everything directly, no per-row buffering ----
    Float_t hitX, hitY, hitZ, projX, projY;
    Int_t detType, stripDir;
    t->SetBranchAddress("detType", &detType);
    t->SetBranchAddress("stripDir", &stripDir);
    t->SetBranchAddress("hitX", &hitX); t->SetBranchAddress("hitY", &hitY); t->SetBranchAddress("hitZ", &hitZ);
    t->SetBranchAddress("projX", &projX); t->SetBranchAddress("projY", &projY);

    for (Long64_t i = 0; i < nAll; i++) {
        t->GetEntry(i);
        double r = sqrt(hitX*hitX + hitY*hitY);
        double phi = atan2(hitY, hitX);
        if (detType == 0) {
            int d = -1;
            for (int k = 0; k < nFst; k++) { if (hitZ >= fstLo[k] && hitZ <= fstHi[k]) { d = k; break; } }
            if (d < 0) continue;
            double phiProj = atan2(projY, projX);
            double dphi = phi - phiProj;
            while (dphi > TMath::Pi())  dphi -= 2*TMath::Pi();
            while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
            double rdphi = r * dphi;
            hFstRdphi[d]->Fill(rdphi);
            h2FstX[d]->Fill(hitX, rdphi);
            h2FstY[d]->Fill(hitY, rdphi);
            h2FstR[d]->Fill(r, rdphi);
            h2FstRphi[d]->Fill(r*phi, rdphi);
            h2FstPhi[d]->Fill(phi, rdphi);
            if (fabs(rdphi) < kCoreCm) h2FstPhiCore[d]->Fill(phi, rdphi);
        } else if (detType == 1) {
            int p = -1;
            for (int k = 0; k < nFtt; k++) { if (hitZ >= fttLo[k] && hitZ <= fttHi[k]) { p = k; break; } }
            if (p < 0) continue;
            int a = -1; double res = 0;
            if (stripDir == 1) { a = 0; res = hitX - projX; }      // V-strip -> dx
            else if (stripDir == 2) { a = 1; res = hitY - projY; } // H-strip -> dy
            if (a < 0) continue;
            hFttRes[p][a]->Fill(res);
            h2FttX[p][a]->Fill(hitX, res);
            h2FttY[p][a]->Fill(hitY, res);
            h2FttR[p][a]->Fill(r, res);
            h2FttRphi[p][a]->Fill(r*phi, res);
            h2FttPhi[p][a]->Fill(phi, res);
        }
    }
    printf("Main pass done, filling plots...\n");

    // ---- FST plots ----
    printf("CHECKPOINT: starting FST plots\n");
    for (int d = 0; d < nFst; d++) {
        printf("CHECKPOINT: FST disk %d\n", d);
        plot1D(hFstRdphi[d], Form("%s/fst%dUnbiasedRdphi.png", outdir, d), "r#upoint#Delta#phi [cm]");
        plot2D(h2FstX[d],    Form("%s/fst%dUnbiasedRdphi_vsX.png", outdir, d), "x [cm]", "r#upoint#Delta#phi [cm]");
        plot2D(h2FstY[d],    Form("%s/fst%dUnbiasedRdphi_vsY.png", outdir, d), "y [cm]", "r#upoint#Delta#phi [cm]");
        plot2D(h2FstR[d],    Form("%s/fst%dUnbiasedRdphi_vsR.png", outdir, d), "r [cm]", "r#upoint#Delta#phi [cm]");
        plot2D(h2FstRphi[d], Form("%s/fst%dUnbiasedRdphi_vsRphi.png", outdir, d), "r#upoint#phi [cm]", "r#upoint#Delta#phi [cm]");
        plot2D(h2FstPhi[d],  Form("%s/fst%dUnbiasedRdphi_vsPhi.png", outdir, d), "#phi [rad]", "r#upoint#Delta#phi [cm]");
        plot2DBare(h2FstPhi[d], Form("%s/fst%dUnbiasedRdphi_vsPhi_bare.png", outdir, d), "#phi [rad]", "r#upoint#Delta#phi [cm]");
    }

    printf("CHECKPOINT: FST plots done, starting FTT plots\n");
    // ---- FTT plots ----
    for (int p = 0; p < nFtt; p++) {
        for (int a = 0; a < 2; a++) {
            printf("CHECKPOINT: FTT plane %d axis %d\n", p, a);
            const char* ax = axName[a];
            const char* ytitle = (a == 0) ? "#Deltax [cm]" : "#Deltay [cm]";
            plot1D(hFttRes[p][a], Form("%s/ftt%dUnbiased%s.png", outdir, p, ax), ytitle);
            plot2D(h2FttX[p][a],    Form("%s/ftt%dUnbiased%s_vsX.png", outdir, p, ax), "x [cm]", ytitle);
            plot2D(h2FttY[p][a],    Form("%s/ftt%dUnbiased%s_vsY.png", outdir, p, ax), "y [cm]", ytitle);
            plot2D(h2FttR[p][a],    Form("%s/ftt%dUnbiased%s_vsR.png", outdir, p, ax), "r [cm]", ytitle);
            plot2D(h2FttRphi[p][a], Form("%s/ftt%dUnbiased%s_vsRphi.png", outdir, p, ax), "r#upoint#phi [cm]", ytitle);
            plot2D(h2FttPhi[p][a],  Form("%s/ftt%dUnbiased%s_vsPhi.png", outdir, p, ax), "#phi [rad]", ytitle);
        }
    }

    printf("CHECKPOINT: FTT plots done, starting translation section\n");
    // ---- Translation-test (sine-vs-phi) plots + HTML subpage ----
    TString transDir = Form("%s/translation", outdir);
    gSystem->mkdir(transDir.Data(), kTRUE);
    FILE* fpT = fopen(Form("%s/index.html", transDir.Data()), "w");
    fprintf(fpT, "<HTML><HEAD><TITLE>More Detailed Look at R*dPhi</TITLE></HEAD>\n");
    fprintf(fpT, "<BODY BGCOLOR=\"#dfdfff\" TEXT=\"black\" LINK=\"blue\" VLINK=\"darkblue\">\n");
    fprintf(fpT, "<H1>More Detailed Look at r#middot#Delta#phi vs #phi</H1><HR>\n");
    fprintf(fpT, "<ul>\n<li>Source: <code>%s</code> (track type: %s), %d files, %lld rows</li>\n",
                 fname, trackTypeLabel, nFiles, nAll);
    fprintf(fpT, "<li>Profile (mean per #phi bin, core cut |r#middot#Delta#phi|<0.3cm) of the 2D "
                 "r#middot#Delta#phi-vs-#phi histogram from <a href=\"../distributions.html\">"
                 "distributions.html</a>, plus a sine fit A#middot sin(#phi-#phi<sub>0</sub>) -- the "
                 "functional form a rigid (x,y) disk shift would produce. The fit describes the data "
                 "badly right now (see Conclusion below) because of the 12-fold wedge-step pattern, "
                 "but zoomed in like this you can also see what looks like short straight-line segments "
                 "within individual wedges, not just flat steps -- worth a closer look once the step "
                 "pattern itself is better understood.</li>\n");
    fprintf(fpT, "<li>Scripts: <code>script/plotFwdAlignmentResidual.C</code>, "
                 "<code>script/testFwdAlignmentTranslation.C</code></li>\n</ul>\n");
    fprintf(fpT, "<table border=1 cellpadding=4>\n<tr><th>Disk</th><th>N(fit pts)</th>"
                 "<th>Amplitude [cm]</th><th>&chi;&sup2;/ndf</th><th>sig</th><th>plot</th></tr>\n");
    for (int d = 0; d < nFst && d < 3; d++) {
        double amp2, ampErr2, chi2ndf2; int nPoints2;
        TString pngCore = Form("fst%d_translation_core.png", d);
        bool okCore = plotSine(h2FstPhiCore[d], Form("%s/%s", transDir.Data(), pngCore.Data()),
                                Form("FST%d (|r#upoint#Delta#phi|<0.3cm)", d),
                                amp2, ampErr2, chi2ndf2, nPoints2);
        fprintf(fpT, "<tr><td>FST%d (z=%.1f)</td><td>%d</td>"
                     "<td>%s</td><td>%.1f</td><td>%s</td>"
                     "<td><a href=\"%s\"><img src=\"%s\" height=150></a></td></tr>\n",
                d, fstZ[d], nPoints2,
                okCore ? Form("%.4f&plusmn;%.4f", amp2, ampErr2) : "n/a",
                chi2ndf2, okCore ? Form("%.1f", (ampErr2>0)?amp2/ampErr2:0) : "n/a",
                pngCore.Data(), pngCore.Data());
    }
    fprintf(fpT, "</table>\n");
    fprintf(fpT, "<H3>Conclusion (2026-07-05)</H3>\n<ul>\n");
    fprintf(fpT, "<li>A single global sin(&phi;-&phi;<sub>0</sub>) is the wrong model &mdash; &chi;&sup2;/ndf "
                 "stays enormous (1000-13000) even after the core cut, because the real unbiased residual "
                 "vs &phi; is a clean <b>12-fold (30&deg;-periodic) step pattern</b> (one step per FST "
                 "wedge), not a smooth sinusoid. See <a href=\"../distributions.html\">"
                 "<code>fst0UnbiasedRdphi_vsPhi.png</code></a>.</li>\n");
    fprintf(fpT, "<li>Candidate explanation (not yet confirmed as the actual cause): in official STAR "
                 "software, <code>StFstHitMaker::Make()</code> (unchanged since it was introduced Jan 2022) "
                 "computes a per-sensor DB/survey alignment correction via "
                 "<code>geoMSensorOnGlobal-&gt;LocalToMaster()</code> (reading "
                 "<code>StFstDb::getRotations()</code>, the <code>Survey_st</code> tables) and then does not "
                 "apply it &mdash; the result is immediately overwritten by an idealized formula-only "
                 "position. Confirmed present in the official <code>star-bnl/star-sw</code> GitHub repo (not "
                 "just our local dev fork).</li>\n");
    fprintf(fpT, "<li><b>Open caveats:</b> (1) we have not checked whether <code>Survey_st</code> is "
                 "actually populated with validated survey offsets for FST, or just unpopulated/identity "
                 "placeholder entries &mdash; if the latter, the discarded value would equal the idealized "
                 "one anyway, and there'd be nothing to explain. (2) Even if populated, we haven't "
                 "independently verified those DB values are themselves correct/calibrated. (3) The "
                 "overwrite may be intentional (a deliberate stopgap until real FST alignment data exists), "
                 "not a bug &mdash; the code's own comment above the overwrite says \"simple transformation... "
                 "need to revisit.\"</li>\n");
    fprintf(fpT, "<li>Decision: not reprocessing from raw DAQ to investigate further right now (expensive) "
                 "&mdash; the planned downstream alignment fit needs to operate <b>per-wedge</b> (12 "
                 "units/disk), not just per-disk, to absorb this step pattern regardless of its root "
                 "cause.</li>\n");
    fprintf(fpT, "<li>Implemented and validated: a per-wedge phi correction applied at MuDst-load time "
                 "cuts the wedge-step amplitude by 42-63%% across all 3 disks (15-file check). See "
                 "<a href=\"../wedgecorrection/index.html\">wedgecorrection/index.html</a> for the "
                 "correction table and validation numbers.</li>\n</ul>\n");
    fprintf(fpT, "<hr>\n<H3>To do</H3>\n<ul>\n<li>Independent cross-check with a B=0, low-luminosity run "
                 "and a real TPC-reconstructed vertex &mdash; BLCVtx's own vertex is built from FST/FTT "
                 "hits, so it can't fully rule out a fit-anchoring effect</li>\n</ul>\n");
    fprintf(fpT, "<HR>OGAWA, Akio<ADDRESS><A HREF=\"mailto://akio@bnl.gov\">akio@bnl.gov</A></ADDRESS>\n</BODY></HTML>\n");
    fclose(fpT);
    printf("Wrote %s/index.html\n", transDir.Data());

    // ---- Distributions subpage (FST rdphi + FTT dx/dy, vs x/y/r/rphi/phi) ----
    FILE* fpD = fopen(Form("%s/distributions.html", outdir), "w");
    fprintf(fpD, "<HTML><HEAD><TITLE>Unbiased FST/FTT Residual Distributions</TITLE></HEAD>\n");
    fprintf(fpD, "<BODY BGCOLOR=\"#dfdfff\" TEXT=\"black\" LINK=\"blue\" VLINK=\"darkblue\">\n");
    fprintf(fpD, "<H1>Unbiased Residual vs x, y, r, r#phi, #phi</H1><HR>\n");
    fprintf(fpD, "<ul>\n<li>Source: <code>%s</code> (track type: %s), %d files, %lld rows</li>\n",
                 fname, trackTypeLabel, nFiles, nAll);
    fprintf(fpD, "<li>Same variables/binning as the biased plots in "
                 "<a href=\"../../index.html\">residual/index.html</a>, built from the unbiased "
                 "hit-removed ntuple instead (<code>script/plotFwdAlignmentResidual.C</code>)</li>\n");
    fprintf(fpD, "<li>FST: combined r#middot#Delta#phi azimuthal residual, axis range #plusmn0.5&nbsp;cm. "
                 "FTT: per-strip #Deltax/#Deltay, axis range #plusmn2&nbsp;cm.</li>\n");
    fprintf(fpD, "<li>A tilted band (red line = linear fit to per-bin mean) &rarr; rotation-type effect. "
                 "Nonzero <a href=\"translation/index.html\">sine-vs-#phi amplitude</a> &rarr; "
                 "translation-type effect.</li>\n");
    fprintf(fpD, "<li>FST \"vs #phi\" column: a clean 12-fold (30&deg;) step pattern, one per wedge "
                 "&mdash; candidate explanation (not yet confirmed): <code>StFstHitMaker::Make()</code> "
                 "computes a DB/survey alignment correction but doesn't apply it. See "
                 "<a href=\"translation/index.html\">translation "
                 "test page</a> for details, and <a href=\"wedgecorrection/index.html\">"
                 "wedgecorrection/index.html</a> for the downstream per-wedge phi correction + "
                 "validation. Bare (no-fit) versions of the vs-#phi plots: "
                 "<a href=\"fst0UnbiasedRdphi_vsPhi_bare.png\">FST0</a>, "
                 "<a href=\"fst1UnbiasedRdphi_vsPhi_bare.png\">FST1</a>, "
                 "<a href=\"fst2UnbiasedRdphi_vsPhi_bare.png\">FST2</a>.</li>\n</ul>\n");
    printf("CHECKPOINT: distributions.html intro written\n");
    fprintf(fpD, "<H3>FST (r#middot#Delta#phi)</H3>\n<table border=1 cellpadding=4>\n"
                 "<tr><th>Disk</th><th>1D</th><th>vs x</th><th>vs y</th><th>vs r</th><th>vs r#phi</th><th>vs #phi</th></tr>\n");
    // Built as a single concatenated TString per row + one fputs, rather than one
    // giant fprintf with ~20+ varargs -- CINT's interpreted variadic-call bridge
    // segfaults past a certain argument count (hit this empirically: the analogous
    // FTT row below, at 23 args, crashed; this FST row, at 14, happened not to).
    for (int d = 0; d < nFst; d++) {
        TString row = Form("<tr><td>FST%d (z=%.1f)</td>", d, fstZ[d]);
        const char* suf[6] = {"", "_vsX", "_vsY", "_vsR", "_vsRphi", "_vsPhi"};
        for (int s = 0; s < 6; s++) {
            TString png = Form("fst%dUnbiasedRdphi%s.png", d, suf[s]);
            row += Form("<td><a href=\"%s\"><img src=\"%s\" height=120></a></td>", png.Data(), png.Data());
        }
        row += "</tr>\n";
        fputs(row.Data(), fpD);
    }
    printf("CHECKPOINT: FST table done\n");
    fprintf(fpD, "</table>\n<H3>FTT (#Deltax V-strip, #Deltay H-strip)</H3>\n<table border=1 cellpadding=4>\n"
                 "<tr><th>Plane</th><th>axis</th><th>1D</th><th>vs x</th><th>vs y</th><th>vs r</th><th>vs r#phi</th><th>vs #phi</th></tr>\n");
    for (int p = 0; p < nFtt; p++) {
        for (int a = 0; a < 2; a++) {
            const char* ax = axName[a];
            TString row = Form("<tr><td>FTT%d (z=%.1f)</td><td>%s</td>", p, fttZ[p], ax);
            TString png1D = Form("ftt%dUnbiased%s.png", p, ax);
            row += Form("<td><a href=\"%s\"><img src=\"%s\" height=120></a></td>", png1D.Data(), png1D.Data());
            const char* suf[5] = {"_vsX", "_vsY", "_vsR", "_vsRphi", "_vsPhi"};
            for (int s = 0; s < 5; s++) {
                TString png = Form("ftt%dUnbiased%s%s.png", p, ax, suf[s]);
                row += Form("<td><a href=\"%s\"><img src=\"%s\" height=120></a></td>", png.Data(), png.Data());
            }
            row += "</tr>\n";
            fputs(row.Data(), fpD);
        }
    }
    fprintf(fpD, "</table>\n<HR>OGAWA, Akio<ADDRESS><A HREF=\"mailto://akio@bnl.gov\">akio@bnl.gov</A></ADDRESS>\n</BODY></HTML>\n");
    fclose(fpD);
    printf("Wrote %s/distributions.html\n", outdir);

    // ---- Save all booked histograms to a ROOT file alongside the PNGs ----
    printf("CHECKPOINT: writing histograms to root file\n");
    TFile* fout = new TFile(Form("%s/alignmentResidualHistos.root", outdir), "RECREATE");
    fout->cd();
    for (int d = 0; d < nFst; d++) {
        hFstRdphi[d]->Write();
        h2FstX[d]->Write(); h2FstY[d]->Write(); h2FstR[d]->Write();
        h2FstRphi[d]->Write(); h2FstPhi[d]->Write(); h2FstPhiCore[d]->Write();
    }
    for (int p = 0; p < nFtt; p++) {
        for (int a = 0; a < 2; a++) {
            hFttRes[p][a]->Write();
            h2FttX[p][a]->Write(); h2FttY[p][a]->Write(); h2FttR[p][a]->Write();
            h2FttRphi[p][a]->Write(); h2FttPhi[p][a]->Write();
        }
    }
    fout->Close();
    printf("Wrote %s/alignmentResidualHistos.root\n", outdir);

    printf("Done.\n");
}
