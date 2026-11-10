// fwoPairRdep.C -- the overlap width versus radius, from the pair centroid.
//
//   root -l -b -q 'script/fwoPairRdep.C("fwo_pairR.root")'
//
// This is the measurement the pitch-vs-radial-shift split rests on, re-made
// with the estimator that needs no conventions. W(r) = <s> + 1 with
// s = d_me + d_nb, per rStrip, where the earlier version of the same number
// came from an equivalent width and therefore from a plateau definition, a
// background definition and a decision about strip 0.
//
// The model, as in script/fstOverlapRdep.C:
//     W(r) = G/p + 128*eps - 128*d/r
// G = 1 deg dead gap on the wedge centreline, p = 30/128 deg nominal pitch,
// eps the fractional pitch error, d the radial offset of the point the strip
// fan converges on. Linear in 1/r, so a two-parameter fit needs no prior:
// pitch alone is flat in r, a fan converging off-axis falls as 1/r.
//
// STRIP 0 IS DROPPED throughout. Its occupancy is 12% above the bulk while its
// coincidence fraction is half the plateau -- edge-truncated clusters, pinned
// to the last strip with their true position outside the silicon, so no
// partner can exist. A pair carries that contamination twice, once per side.
// The outer has four more strips of ridge and can afford it.

const int kN = 12;
const int kBkgFrom = 14;
const int kSMax = 8;

void fwoPairRdep(const char *fn = "fwo_pairR.root")
{
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open(fn);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", fn); return; }

    // per-strip efficiency, folded over both wedge ends, from the occupancy.
    // Not split by rStrip: efficiency is a property of the strip in phi, and
    // with strip 0 dropped the surviving correction is only a few percent.
    TH2F *hp = (TH2F *) f->Get("hPhiStripPerSensor");
    double occ[128];
    for (int c = 0; c < 128; c++) {
        double v = 0;
        for (int sv = 1; sv <= 2; sv++)
            v += hp->GetBinContent(c + 3, hp->GetYaxis()->FindBin(sv));
        occ[c] = v;
    }
    double bulk = 0;
    for (int c = 20; c <= 107; c++) bulk += occ[c];
    bulk /= 88.0;
    double eff[kN];
    for (int d = 0; d < kN; d++) eff[d] = (occ[d] + occ[127 - d]) / (2.0 * bulk);

    double r[4]  = { 17.94, 20.81, 23.69, 26.56 };
    double W[4], eW[4], er[4] = {0,0,0,0};
    double spec[4][kSMax + 1], sMean[4];   // the s distributions, for the overlay

    printf("\nOUTER overlap width per rStrip, pair centroid, strip 0 dropped\n");
    printf("  %-8s %6s %9s %9s %9s %10s\n", "rStrip", "<r>", "<s>", "W", "+- stat", "entries");
    for (int ib = 0; ib < 4; ib++) {
        TH2F *h = (TH2F *) f->Get(Form("hPairR%d", ib + 4));
        if (!h) { printf("  hPairR%d absent\n", ib + 4); return; }
        double cor[kN][kN], rawc[kN][kN];
        for (int i = 0; i < kN; i++)
            for (int j = 0; j < kN; j++) {
                rawc[i][j] = h->GetBinContent(i + 1, j + 1);
                cor[i][j]  = rawc[i][j] / (eff[i] * eff[j]);
            }
        double bkg = 0; int nb = 0;
        for (int i = 0; i < kN; i++)
            for (int j = 0; j < kN; j++)
                if (i + j >= kBkgFrom) { bkg += cor[i][j]; nb++; }
        bkg /= nb;

        double num = 0, den = 0, nraw = 0;
        for (int sv = 0; sv <= kSMax; sv++) {
            double rowv = 0;
            for (int i = 1; i < kN; i++) {
                int j = sv - i; if (j < 1 || j >= kN) continue;
                double v = cor[i][j] - bkg;
                rowv += v; num += sv * v; den += v; nraw += rawc[i][j];
            }
            spec[ib][sv] = rowv;
        }
        double m = num / den, var = 0;
        for (int sv = 0; sv <= kSMax; sv++)
            for (int i = 1; i < kN; i++) {
                int j = sv - i; if (j < 1 || j >= kN) continue;
                var += (cor[i][j] - bkg) * (sv - m) * (sv - m);
            }
        double rms = sqrt(var / den);
        W[ib]  = m + 1;
        sMean[ib] = m;
        for (int sv = 0; sv <= kSMax; sv++) spec[ib][sv] /= den;   // unit area
        eW[ib] = rms / sqrt(nraw);         // centroid error, from the raw count
        printf("  %-8d %6.2f %9.3f %9.3f %9.3f %10.0f\n", ib + 4, r[ib], m, W[ib], eW[ib], nraw);
    }

    // ---- free two-parameter fit W = A + B/r ------------------------------
    double ub = 0, wb = 0, sw = 0;
    for (int i = 0; i < 4; i++) { double w = 1.0 / (eW[i] * eW[i]); ub += w / r[i]; wb += w * W[i]; sw += w; }
    ub /= sw; wb /= sw;
    double Suu = 0, Suw = 0;
    for (int i = 0; i < 4; i++) {
        double w = 1.0 / (eW[i] * eW[i]), du = 1.0 / r[i] - ub;
        Suu += w * du * du; Suw += w * du * (W[i] - wb);
    }
    double B = Suw / Suu, A = wb - B * ub;
    const double pitch = 30.0 / 128.0, gapS = 1.0 / pitch;
    double eps = (A - gapS) / 128.0, dsh = -B / 128.0;
    double eB = sqrt(1.0 / Suu), eEps = (eB / 128.0) / 128.0, eD = eB / 128.0;

    // FIT QUALITY, printed because the two-parameter model does NOT describe
    // these four points and a d quoted without that is an overclaim.
    double c2 = 0, c2const = 0, wsum = 0, wmean = 0;
    for (int i = 0; i < 4; i++) { double w = 1.0/(eW[i]*eW[i]); wsum += w; wmean += w*W[i]; }
    wmean /= wsum;
    for (int i = 0; i < 4; i++) {
        double res = (W[i] - (A + B/r[i])) / eW[i];
        c2 += res*res;
        double rc = (W[i] - wmean) / eW[i];
        c2const += rc*rc;
    }
    double scale = sqrt(c2 / 2.0);

    printf("\n  free fit  W(r) = %.4f + %.4f/r\n", A, B);
    printf("     pitch  %+.3f%%  ->  128 strips span %.3f deg\n", eps * 100, 30 * (1 + eps));
    printf("     d      %+.2f +- %.2f mm  (stat only)\n", dsh * 10, eD * 10);
    printf("\n  FIT QUALITY -- read this before quoting d:\n");
    printf("     A + B/r   chi2 = %.2f / 2 dof   chi2/ndf = %.1f\n", c2, c2/2.0);
    printf("     residuals ");
    for (int i = 0; i < 4; i++) printf(" %+.1f", (W[i] - (A + B/r[i])) / eW[i]);
    printf(" sigma\n");
    printf("     constant  chi2 = %.2f / 3 dof   <- pure pitch, excluded outright\n", c2const);
    printf("     measured steps between adjacent radii:");
    for (int i = 0; i < 3; i++) printf(" %+.3f", W[i+1] - W[i]);
    printf("\n     A + B/r predicts:                     ");
    for (int i = 0; i < 3; i++) printf(" %+.3f", (A + B/r[i+1]) - (A + B/r[i]));
    printf("\n     A+B/r is CONVEX, the data is flat then accelerating -- opposite\n");
    printf("     curvature, so the model is dragged through points not on it.\n");
    printf("\n     errors scaled to chi2/ndf = 1 (x%.1f):  d = %+.2f +- %.2f mm, %.1f sigma\n",
           scale, dsh * 10, eD * 10 * scale, fabs(dsh / (eD * scale)));
    printf("     yield-profile method gave  30.139 deg,  -0.68 mm\n");

    // ---- the plot --------------------------------------------------------
    // CINT shares one scope for the whole macro, so "c" collides with the
    // loop variable used for the occupancy above. Unique names throughout.
    TCanvas *cvRd = new TCanvas("cvRd", "W(r) from the pair centroid", 1200, 520);
    cvRd->Divide(2, 1);
    // LEFT: the raw evidence -- the s = d_me + d_nb spectrum at each radius,
    // strip 0 dropped, accidentals subtracted, each normalised to unit area so
    // the SHAPES compare. The shift is only 0.18 strips on an rms of 1.15, so
    // it is not meant to be read off the curves by eye: the dashed verticals
    // are the centroids, and those are what the fit uses.
    cvRd->cd(1);
    gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    {
        TH2F *fs2 = new TH2F("fs2", "the ridge at each radius, strip 0 dropped;s = d_{me} + d_{nb};fraction of the ridge", 10, 1.5, 7.5, 10, 0, 0.78);
        fs2->GetXaxis()->SetTitleSize(0.05); fs2->GetYaxis()->SetTitleSize(0.05);
        fs2->Draw();
        int colR[4] = { kRed + 1, kOrange + 7, kGreen + 2, kBlue + 1 };
        int mstR[4] = { 20, 21, 22, 23 };
        TLegend *lr = new TLegend(0.50, 0.62, 0.97, 0.90);
        lr->SetFillStyle(0); lr->SetBorderSize(0); lr->SetTextSize(0.034);
        for (int ic = 0; ic < 4; ic++) {
            double xx[kSMax + 1], yy[kSMax + 1];
            for (int sv = 0; sv <= kSMax; sv++) { xx[sv] = sv; yy[sv] = spec[ic][sv]; }
            TGraph *gg = new TGraph(kSMax + 1, xx, yy);
            gg->SetMarkerStyle(mstR[ic]); gg->SetMarkerSize(1.3);
            gg->SetMarkerColor(colR[ic]); gg->SetLineColor(colR[ic]); gg->SetLineWidth(2);
            gg->Draw("PL same");
            TLine *lv = new TLine(sMean[ic], 0, sMean[ic], 0.66);
            lv->SetLineColor(colR[ic]); lv->SetLineStyle(2); lv->SetLineWidth(2); lv->Draw();
            lr->AddEntry(gg, Form("rStrip %d, r = %.1f cm:  #LTs#GT = %.3f",
                                  ic + 4, r[ic], sMean[ic]), "pl");
        }
        lr->Draw();
        TLatex tr; tr.SetNDC(); tr.SetTextSize(0.033);
        tr.DrawLatex(0.17, 0.30, "dashed: the centroid of each");
        tr.DrawLatex(0.17, 0.255, "they walk left as r grows:");
        tr.DrawLatex(0.17, 0.21, "4.407 #rightarrow 4.229, a 0.18 strip shift");
        tr.DrawLatex(0.17, 0.155, "pure pitch would put all four together");
    }
    cvRd->cd(2);
    gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    double ymin = 1e9, ymax = -1e9;
    for (int i = 0; i < 4; i++) { ymin = TMath::Min(ymin, W[i]); ymax = TMath::Max(ymax, W[i]); }
    TH2F *fr = new TH2F("fr", "Overlap width vs radius, from the pair centroid (strip 0 dropped);"
                              "mean radius  r  [cm];overlap width  W = #LTs#GT + 1  [strips]",
                        10, 16.5, 28.0, 10, ymin - 0.52, ymax + 0.46);
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045);
    fr->Draw();

    TF1 *fp = new TF1("fp", "[0]", 16.5, 28.0);
    fp->SetParameter(0, (W[0] + W[1] + W[2] + W[3]) / 4.0);
    fp->SetLineColor(kBlue); fp->SetLineStyle(2); fp->SetLineWidth(3); fp->Draw("same");

    TF1 *fs = new TF1("fs", "[0] - 128*[1]/x", 16.5, 28.0);
    fs->SetParameters(gapS + 128 * ((A - gapS) / 128.0 + 0), -0.176);
    // pure shift normalised to pass through the innermost point
    fs->SetParameter(0, W[0] + 128 * (-0.176) / r[0]);
    fs->SetLineColor(kGreen + 2); fs->SetLineStyle(3); fs->SetLineWidth(3); fs->Draw("same");

    TF1 *ff = new TF1("ff", "[0] + [1]/x", 16.5, 28.0);
    ff->SetParameters(A, B);
    ff->SetLineColor(kRed); ff->SetLineWidth(3); ff->Draw("same");

    TGraphErrors *g = new TGraphErrors(4, r, W, er, eW);
    g->SetMarkerStyle(20); g->SetMarkerSize(1.6); g->SetLineWidth(2); g->Draw("P same");

    TLegend *l = new TLegend(0.50, 0.73, 0.97, 0.93);
    l->SetFillStyle(0); l->SetBorderSize(0); l->SetTextSize(0.034);
    l->AddEntry(g,  "pair centroid, 20 files", "pe");
    l->AddEntry(fp, "pure pitch: flat in r", "l");
    l->AddEntry(fs, "pure radial shift: 1/r", "l");
    l->AddEntry(ff, "free fit: A + B/r", "l");
    l->Draw();

    // The text goes in the band opened up below the curves, not over them:
    // the pure-shift curve runs diagonally through the middle of the frame and
    // was being hidden by this block.
    TLatex t; t.SetNDC(); t.SetTextSize(0.032);
    t.DrawLatex(0.17, 0.260, Form("free fit:  d = %+.2f #pm %.2f mm,  #chi^{2}/ndf = %.1f",
                                  dsh*10, eD*10*scale, c2/2.0));
    t.DrawLatex(0.17, 0.215, "ridge position only -- no plateau, no background convention");
    t.DrawLatex(0.17, 0.170, Form("#chi^{2}/ndf = %.1f: the 1/r form does NOT fit, see the page", c2/2.0));
    t.DrawLatex(0.17, 0.125, "yield-profile method: 30.139#circ, -0.68 mm");

    // redraw the curves last so nothing can sit on top of them
    fp->Draw("same"); fs->Draw("same"); ff->Draw("same"); g->Draw("P same");

    cvRd->SaveAs("FstOuterSensorOverlap/plots/fwo_pair_rdep.png");
    printf("wrote FstOuterSensorOverlap/plots/fwo_pair_rdep.png\n");
}
