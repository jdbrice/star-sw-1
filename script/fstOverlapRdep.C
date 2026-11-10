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
    // Primary points are the 48-job condor sample merged 2026-11-10:
    // 60,376,989 events, 500x the 120k run the first version of this page used.
    // Statistical errors are ~0.003 strips and invisible; the bars drawn are the
    // ESTIMATOR systematic, the spread over reasonable background and plateau
    // conventions, which is what actually limits these numbers.
    //
    // THE ESTIMATOR MATTERS, and getting it wrong looks like a real effect. The
    // d=0 bin is edge-suppressed and its suppression is not even monotonic in r
    // ((n0-tail)/plateau = 0.457, 0.304, 0.307, 0.585 across rStrip 4-7), so a
    // width that uses its measured count comes out 4.81, 4.59, 4.52, 4.78 --
    // non-monotonic, trend destroyed. Counting strip 0 as one full strip of
    // overlap instead, which is what it is, gives the monotonic row below. The
    // background subtracted is the pred row's own far tail (mean of d=7,8,9),
    // NOT the mirrored cell: the mirror is 68.9k per bin where the tail is
    // 329k, because the pred row carries same-event combinatoric background
    // that the mirrored cell does not see.
    const int    nR = 4;
    double r[nR]  = { 17.94, 20.81, 23.69, 26.56 };
    double w[nR]  = { 5.347, 5.281, 5.211, 5.198 };   // 60.4M events
    double er[nR] = {  0.00,  0.00,  0.00,  0.00 };
    double ew[nR] = { 0.090, 0.090, 0.090, 0.090 };   // estimator systematic

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
                         10, 16.5, 28.0, 10, 5.00, 5.75);
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
    l1->AddEntry(g, "measured, 60.4M events", "pe");
    l1->AddEntry(fPitch, Form("pure pitch  (+%.2f%%): flat", epsPure * 100), "l");
    l1->AddEntry(fShift, Form("pure shift  (d = %.2f cm): 1/r", dPure), "l");
    l1->AddEntry(fFit, Form("free fit: pitch %+.2f%%, d = %+.2f mm", eps * 100, dsh * 10), "l");
    l1->Draw();

    TLatex t; t.SetNDC(); t.SetTextSize(0.034);
    t.DrawLatex(0.16, 0.22, Form("128 strips span %.3f#circ, not 30.000#circ", 30 * (1 + eps)));
    t.DrawLatex(0.16, 0.17, "bars are the estimator systematic; statistical error is ~0.003");

    // ================= RIGHT: inner, MEASURED ============================
    // The inner overlap is 1-2 strips, so there is no plateau to normalise an
    // equivalent width to, and the partner-finding efficiency cannot be read off
    // its own profile. A SHAPE ratio needs neither: n(d=2)/n(d=0) with the flat
    // accidental pedestal (the mirrored cell) subtracted. Efficiency multiplies
    // every distance bin alike, so it cancels.
    //
    // The point of the panel: a pitch error gives the same overlap in STRIPS at
    // every radius, hence the same profile shape at every radius, hence a FLAT
    // line here. The data falls by a factor 5.1. No model of the profile shape
    // is needed to say that.
    c->cd(2);
    gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);

    const int nI = 4;
    double ri[nI]  = {  6.44,  9.31, 12.19, 15.06 };
    double qi[nI]  = { 0.765, 0.424, 0.304, 0.150 };
    double eqi[nI] = { 0.023, 0.013, 0.014, 0.010 };
    double eri[nI] = { 0, 0, 0, 0 };

    TH2F *fr2 = new TH2F("fr2", "INNER sensor (no 1#circ gap): measured shape vs radius;"
                                "mean radius  r  [cm];n(d=2) / n(d=0),  accidentals subtracted",
                         10, 5.0, 16.5, 10, 0.0, 1.0);
    fr2->GetXaxis()->SetTitleSize(0.05); fr2->GetYaxis()->SetTitleSize(0.05);
    fr2->Draw();

    TF1 *flat = new TF1("flat", "[0]", 5.0, 16.5);
    flat->SetParameter(0, qi[0]);
    flat->SetLineColor(kBlue); flat->SetLineStyle(2); flat->SetLineWidth(3);
    flat->Draw("same");

    TGraphErrors *gi = new TGraphErrors(nI, ri, qi, eri, eqi);
    gi->SetMarkerStyle(21); gi->SetMarkerSize(1.5); gi->SetLineWidth(2);
    gi->Draw("P same");

    // no TLegend here: its marker sample lands inside the frame and reads as a
    // fifth data point. Label the line and the points directly instead.
    TLatex t2; t2.SetNDC(); t2.SetTextSize(0.034);
    t2.SetTextColor(kBlue);
    t2.DrawLatex(0.42, 0.835, "what pure pitch requires: flat");
    t2.SetTextColor(kBlack);
    t2.DrawLatex(0.42, 0.665, "#bullet  measured, inner rStrip 0-3");
    t2.DrawLatex(0.42, 0.585, "a pitch error gives the same overlap");
    t2.DrawLatex(0.42, 0.540, "in strips at every radius, so the same");
    t2.DrawLatex(0.42, 0.495, "profile shape");
    t2.DrawLatex(0.17, 0.30, Form("falls by #times%.1f  #Rightarrow  pure pitch excluded", qi[0] / qi[3]));
    t2.DrawLatex(0.17, 0.25, "a radial shift of the strip fan is required");
    t2.DrawLatex(0.17, 0.17, "200,956 events; the 48-job rerun adds 500#times");

    c->SaveAs("FstOuterSensorOverlap/plots/fwo_rdep.png");
    printf("\nfree fit W(r) = %.5f + %.4f/r\n", A, B);
    printf("  pitch error %+.4f%% -> 128 strips span %.3f deg\n", eps * 100, 30 * (1 + eps));
    printf("  radial offset d = %+.4f cm = %+.2f mm\n", dsh, dsh * 10);
    printf("wrote FstOuterSensorOverlap/plots/fwo_rdep.png\n");
}
