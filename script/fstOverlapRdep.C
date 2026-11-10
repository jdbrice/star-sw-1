// fstOverlapRdep.C -- the radial dependence of the FST wedge-overlap width, and
// what the inner sensor would add.
//
//   root -l -b -q script/fstOverlapRdep.C
//
// LEFT: the measurement. Overlap width in strips vs mean radius, outer sensor
// rStrip 4-7, from script/fstWedgeOverlap.C (FstOuterSensorOverlap/index.html).
// Three curves on the same axes:
//   pure pitch error      W = gap/pitch + 128*eps          flat in r
//   pure radial shift     W = gap/pitch - 128*d/r          falls as 1/r
//   free 2-parameter fit  W = A + B/r
// The model is W(r) = G/p + 128*eps - 128*d/r, with G the dead gap on the wedge
// centreline (1 deg for the outer, ZERO for the inner), p = 30/128 deg nominal
// pitch, eps the fractional pitch error, d the radial offset of the point the
// strip fan converges on. Pitch widens the overlap by the same amount at every
// radius; a fan that converges off-axis widens it as 1/r. That is what makes
// them separable, and it is the whole content of this plot.
//
// RIGHT: why the inner sensor is the SHARPER test, not the weaker one. The
// inner has no 1 deg gap, so its width IS the signal rather than a 4% ripple on
// a 4.27-strip pedestal. The three hypotheses give 0.99 / 1.92 / 3.50 strips at
// rStrip 0 -- a factor 3.5 spread, against 5% across the whole outer range.
//
// Errors: the published widths are quoted to 0.01 strips and the raw
// coincidence counts behind them were not kept, so the bars here are the
// residual scatter about the fit, 0.017 strips (2 dof). That is an estimate of
// the scatter, NOT a propagated statistical error -- it needs the counts back
// from a re-run to be a real uncertainty, and it is drawn dashed to say so.

void fstOverlapRdep()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(0.055);

    const double pitch  = 30.0 / 128.0;     // nominal, deg
    const double gapStr = 1.0 / pitch;      // 1 deg gap in strips = 4.2667

    // ---- the measurement, outer rStrip 4-7 -------------------------------
    const int    nR = 4;
    double r[nR]  = { 17.94, 20.81, 23.69, 26.56 };
    double w[nR]  = {  5.26,  5.22,  5.14,  5.10 };
    double er[nR] = {  0.00,  0.00,  0.00,  0.00 };
    double ew[nR] = { 0.017, 0.017, 0.017, 0.017 };

    // free fit W = A + B/r, done in python alongside this macro; repeated here
    // so the plot is self-contained and the numbers cannot drift apart
    double ub = 0, wb = 0;
    for (int i = 0; i < nR; i++) { ub += 1.0 / r[i]; wb += w[i]; }
    ub /= nR; wb /= nR;
    double Suu = 0, Suw = 0;
    for (int i = 0; i < nR; i++) {
        double du = 1.0 / r[i] - ub;
        Suu += du * du;
        Suw += du * (w[i] - wb);
    }
    const double B = Suw / Suu;
    const double A = wb - B * ub;
    const double eps = (A - gapStr) / 128.0;        // fractional pitch error
    const double dsh = -B / 128.0;                  // radial offset, cm

    // the two single-hypothesis comparisons quoted in the note
    const double epsPure = 0.0077;                  // pitch alone explains 5.25
    const double dPure   = -0.176;                  // shift alone, cm

    TCanvas *c = new TCanvas("c", "FST overlap width vs radius", 1200, 520);
    c->Divide(2, 1);

    // ================= LEFT: outer, measured ==============================
    c->cd(1);
    gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);

    TH2F *fr1 = new TH2F("fr1", "OUTER sensor: measured, and what each hypothesis needs;"
                                "mean radius  r  [cm];overlap width  [strips]",
                         10, 16.5, 28.0, 10, 4.95, 5.70);
    fr1->GetXaxis()->SetTitleSize(0.05); fr1->GetYaxis()->SetTitleSize(0.05);
    fr1->Draw();

    TF1 *fPitch = new TF1("fPitch", "[0]", 16.5, 28.0);
    fPitch->SetParameter(0, gapStr + 128 * epsPure);
    fPitch->SetLineColor(kBlue); fPitch->SetLineStyle(2); fPitch->SetLineWidth(3);
    fPitch->Draw("same");

    TF1 *fShift = new TF1("fShift", "[0] - 128*[1]/x", 16.5, 28.0);
    fShift->SetParameters(gapStr, dPure);
    fShift->SetLineColor(kGreen + 2); fShift->SetLineStyle(3); fShift->SetLineWidth(3);
    fShift->Draw("same");

    TF1 *fFit = new TF1("fFit", "[0] + [1]/x", 16.5, 28.0);
    fFit->SetParameters(A, B);
    fFit->SetLineColor(kRed); fFit->SetLineWidth(3);
    fFit->Draw("same");

    TGraphErrors *g = new TGraphErrors(nR, r, w, er, ew);
    g->SetMarkerStyle(20); g->SetMarkerSize(1.5); g->SetLineWidth(2);
    g->Draw("P same");

    TLegend *l1 = new TLegend(0.40, 0.60, 0.97, 0.90);
    l1->SetFillStyle(0); l1->SetBorderSize(0); l1->SetTextSize(0.035);
    l1->AddEntry(g, "measured, rStrip 4-7", "pe");
    l1->AddEntry(fPitch, Form("pure pitch  (+%.2f%%): flat", epsPure * 100), "l");
    l1->AddEntry(fShift, Form("pure shift  (d = %.2f cm): 1/r", dPure), "l");
    l1->AddEntry(fFit, Form("free fit: pitch %+.2f%%, d = %+.2f mm", eps * 100, dsh * 10), "l");
    l1->Draw();

    TLatex t; t.SetNDC(); t.SetTextSize(0.034);
    t.DrawLatex(0.16, 0.22, Form("128 strips span %.3f#circ, not 30.000#circ", 30 * (1 + eps)));
    t.DrawLatex(0.16, 0.17, "measured trend is between the two, much closer to flat");

    // ================= RIGHT: inner, predicted ============================
    c->cd(2);
    gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);

    TH2F *fr2 = new TH2F("fr2", "INNER sensor (no 1#circ gap): the same three hypotheses;"
                                "mean radius  r  [cm];overlap width  [strips]",
                         10, 5.0, 17.0, 10, 0.0, 4.0);
    fr2->GetXaxis()->SetTitleSize(0.05); fr2->GetYaxis()->SetTitleSize(0.05);
    fr2->Draw();

    TF1 *iPitch = new TF1("iPitch", "[0]", 5.0, 17.0);
    iPitch->SetParameter(0, 128 * epsPure);
    iPitch->SetLineColor(kBlue); iPitch->SetLineStyle(2); iPitch->SetLineWidth(3);
    iPitch->Draw("same");

    TF1 *iShift = new TF1("iShift", "-128*[0]/x", 5.0, 17.0);
    iShift->SetParameter(0, dPure);
    iShift->SetLineColor(kGreen + 2); iShift->SetLineStyle(3); iShift->SetLineWidth(3);
    iShift->Draw("same");

    TF1 *iFit = new TF1("iFit", "[0] + [1]/x", 5.0, 17.0);
    iFit->SetParameters(A - gapStr, B);
    iFit->SetLineColor(kRed); iFit->SetLineWidth(3);
    iFit->Draw("same");

    // the one inner number we have: 1.0 strip at <r> = 12.81, and the +-0.30
    // is the quoted inner/outer width ratio error, so it is a ~30% band
    double ri[1] = { 12.81 }, wi[1] = { 1.00 }, eri[1] = { 0.0 }, ewi[1] = { 0.30 };
    TGraphErrors *gi = new TGraphErrors(1, ri, wi, eri, ewi);
    gi->SetMarkerStyle(21); gi->SetMarkerSize(1.6); gi->SetLineWidth(2);
    gi->Draw("P same");

    // where the four inner rStrip bins would land
    double rb[4] = { 6.44, 9.31, 12.19, 15.06 };
    for (int i = 0; i < 4; i++) {
        TLine *ln = new TLine(rb[i], 0.0, rb[i], 3.6);
        ln->SetLineColor(kGray + 1); ln->SetLineStyle(7); ln->Draw();
    }

    TLegend *l2 = new TLegend(0.36, 0.62, 0.97, 0.90);
    l2->SetFillStyle(0); l2->SetBorderSize(0); l2->SetTextSize(0.035);
    l2->AddEntry(gi, "measured, inner summed over r", "pe");
    l2->AddEntry(iPitch, "pure pitch: 0.99 at every r", "l");
    l2->AddEntry(iShift, "pure shift: 3.50 #rightarrow 1.50", "l");
    l2->AddEntry(iFit, "this fit: 1.92 #rightarrow 1.10", "l");
    l2->Draw();

    TLatex t2; t2.SetNDC(); t2.SetTextSize(0.034);
    t2.DrawLatex(0.16, 0.25, "dashed grey: rStrip 0,1,2,3 bin centres");
    t2.DrawLatex(0.16, 0.20, "spread between hypotheses at rStrip 0: #times3.5");
    t2.DrawLatex(0.16, 0.15, "against #times1.05 anywhere on the left panel");

    c->SaveAs("FstOuterSensorOverlap/plots/fwo_rdep.png");
    printf("\nfree fit W(r) = %.5f + %.4f/r\n", A, B);
    printf("  pitch error %+.4f%% -> 128 strips span %.3f deg\n", eps * 100, 30 * (1 + eps));
    printf("  radial offset d = %+.4f cm = %+.2f mm\n", dsh, dsh * 10);
    printf("wrote FstOuterSensorOverlap/plots/fwo_rdep.png\n");
}
