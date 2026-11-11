// fwoPairInner.C -- the inner sensors through the pair-centroid method, and
// both sensors against what the SENSOR MASK says the overlap should be.
//
// WHY THE INNER IS WORTH PLOTTING THIS WAY
//   The inner sensors have no 1 deg gap: AGML places them phi1 = -15, phi2 =
//   +15 on a 30 deg wedge pitch, so adjacent inner sensors abut exactly and
//   the design overlap is ZERO. Anything measured here is as-built mounting,
//   with no pedestal to sit on. The outer's 1 deg design overlap is
//   1.0 / (30/128) = 4.2667 strips, flat in r.
//
//   The earlier reading was that the inner centroid cannot work because the
//   ridge sits at s = 0, hard against the s >= 0 wall. On the full 60.4M
//   sample that is only true at the LARGE inner radii. At rStrip 0 the
//   spectrum peaks at s = 3 with s = 0 at less than half the peak, so the
//   ridge is off the wall and <s> is usable there. It walks inward as r grows
//   and only then runs into the wall -- which biases <s> UP at large r, so the
//   measured fall is a LOWER limit on the true one.
//
// Usage:
//   root4star -l -b -q 'script/fwoPairInner.C("<hadd-ed fwo file>")'

const int kFPIsMax = 10;

void fwoPairInner(const char *fn = "fwo_pairR.root",
                  const char *sampleLabel = "60.4M events, 48 files")
{
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleFontSize(0.055);

   TFile *fIn = TFile::Open(fn);
   if (!fIn || fIn->IsZombie()) { printf("cannot open %s\n", fn); return; }

   const double pitchDeg = 30.0 / 128.0;        // exact, confirmed from the GDS
   const double outerDesign = 1.0 / pitchDeg;   // 4.2667 strips
   const double innerDesign = 0.0;

   double rIn[4] = {  6.44,  9.31, 12.19, 15.06 };
   double rOu[4] = { 17.94, 20.81, 23.69, 26.56 };
   double wIn[4], wOu[4], eIn[4], eOu[4], zero[4] = {0,0,0,0};
   double wOuD[4], eOuD[4];        // outer with the edge strip dropped
   double specIn[4][kFPIsMax + 1];

   printf("INNER, pair centroid per rStrip. Design overlap is ZERO: adjacent\n");
   printf("inner sensors abut at +-15 deg on a 30 deg wedge pitch.\n\n");
   printf("  rStrip   <r>    peak s    <s>       W     +-stat   ridge entries\n");

   for (int side = 0; side < 2; side++) {
      for (int k = 0; k < 4; k++) {
         TH2F *hh = (TH2F*)fIn->Get(side == 0 ? Form("hPairIR%d", k)
                                              : Form("hPairR%d", k + 4));
         if (!hh) { printf("  missing histogram\n"); return; }

         // Flat accidental floor per cell, from the far corner s >= 14, where
         // no real pair can live: the two strips are further apart than any
         // credible overlap.
         double flo = 0; int nfl = 0;
         for (int i = 0; i < 12; i++) for (int j = 0; j < 12; j++)
            if (i + j >= 14) { flo += hh->GetBinContent(i+1, j+1); nfl++; }
         if (nfl) flo /= nfl;

         double sp[24], n2[24], spD[24], n2D[24];
         for (int s = 0; s < 24; s++) { sp[s] = 0; n2[s] = 0; spD[s] = 0; n2D[s] = 0; }
         for (int i = 0; i < 12; i++) for (int j = 0; j < 12; j++) {
            double c = hh->GetBinContent(i+1, j+1);
            sp[i+j] += c - flo;
            n2[i+j] += c;                        // for the statistical error
            // the published outer convention also drops the edge strip on
            // BOTH sides; for the inner that would remove the ridge itself at
            // the large radii, so it is only computed for the outer.
            if (i >= 1 && j >= 1) { spD[i+j] += c - flo; n2D[i+j] += c; }
         }

         // <s> over 0..6. Beyond that the subtracted spectrum is consistent
         // with zero and only adds noise; for the outer the ridge is at 4 and
         // this window is generous either way.
         double num = 0, den = 0, var = 0;
         for (int s = 0; s <= 6; s++) { num += sp[s]*s; den += sp[s]; }
         double sm = den > 0 ? num/den : 0;
         for (int s = 0; s <= 6; s++) var += n2[s]*(s - sm)*(s - sm);
         double err = den > 0 ? sqrt(var)/den : 0;

         double norm = 0;
         for (int s = 0; s <= kFPIsMax; s++) norm += sp[s];

         if (side == 0) {
            wIn[k] = sm + 1; eIn[k] = err;
            for (int s = 0; s <= kFPIsMax; s++) specIn[k][s] = norm > 0 ? sp[s]/norm : 0;
            int pk = 0; for (int s = 1; s <= 6; s++) if (sp[s] > sp[pk]) pk = s;
            printf("    %d   %6.2f   %4d   %7.3f %7.3f  %7.4f  %10.0f\n",
                   k, rIn[k], pk, sm, sm + 1, err, den);
         } else {
            wOu[k] = sm + 1; eOu[k] = err;
            double numD = 0, denD = 0, varD = 0;
            for (int s = 2; s <= 8; s++) { numD += spD[s]*s; denD += spD[s]; }
            double smD = denD > 0 ? numD/denD : 0;
            for (int s = 2; s <= 8; s++) varD += n2D[s]*(s - smD)*(s - smD);
            wOuD[k] = smD + 1; eOuD[k] = denD > 0 ? sqrt(varD)/denD : 0;
         }
      }
   }

   printf("\n  inner W falls %.3f -> %.3f, a factor %.2f in the EXCESS over a\n",
          wIn[0], wIn[3], (wIn[0]-1)/(wIn[3]-1));
   printf("  zero-overlap design. The large-r points are biased UP by the\n");
   printf("  s >= 0 wall, so the true fall is steeper than this.\n");
   printf("\n  OUTER for contrast (same estimator, same window):\n");
   for (int k = 0; k < 4; k++)
      printf("    rStrip %d  r = %5.2f   W = %6.3f (all d)  %6.3f (drop d=0)"
             "   design %.3f, flat\n",
             k + 4, rOu[k], wOu[k], wOuD[k], outerDesign);

   // ---------------------------------------------------------------- plot
   TCanvas *cvPI = new TCanvas("cvPI", "inner overlap from the pair map", 1800, 780);
   cvPI->Divide(2, 1);

   cvPI->cd(1);
   gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
   TH2F *fi1 = new TH2F("fi1", "INNER: the ridge moves with radius;"
                               "s = d_{me} + d_{nb};fraction of the ridge",
                        10, -0.5, 8.5, 10, 0, 0.42);
   fi1->GetXaxis()->SetTitleSize(0.05); fi1->GetYaxis()->SetTitleSize(0.05);
   fi1->Draw();
   int colI[4] = { kRed + 1, kOrange + 7, kGreen + 2, kBlue + 1 };
   int mstI[4] = { 20, 21, 22, 23 };
   TLegend *li = new TLegend(0.46, 0.63, 0.97, 0.90);
   li->SetFillStyle(0); li->SetBorderSize(0); li->SetTextSize(0.034);
   for (int k = 0; k < 4; k++) {
      double xx[kFPIsMax + 1], yy[kFPIsMax + 1];
      for (int s = 0; s <= kFPIsMax; s++) { xx[s] = s; yy[s] = specIn[k][s]; }
      TGraph *gi = new TGraph(kFPIsMax + 1, xx, yy);
      gi->SetMarkerStyle(mstI[k]); gi->SetMarkerSize(1.3);
      gi->SetMarkerColor(colI[k]); gi->SetLineColor(colI[k]); gi->SetLineWidth(2);
      gi->Draw("PL same");
      li->AddEntry(gi, Form("rStrip %d, r = %.1f cm:  W = %.2f", k, rIn[k], wIn[k]), "pl");
   }
   li->Draw();
   TLatex ti; ti.SetNDC(); ti.SetTextSize(0.031);
   ti.DrawLatex(0.45, 0.555, "design overlap for the inner is ZERO:");
   ti.DrawLatex(0.45, 0.510, "adjacent sensors abut at #pm15#circ exactly");
   ti.DrawLatex(0.45, 0.450, "at r = 15 cm the ridge is against the s #geq 0 wall,");
   ti.DrawLatex(0.45, 0.405, "so those W are upper limits and the fall is steeper");

   cvPI->cd(2);
   gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
   TH2F *fi2 = new TH2F("fi2", "Apparent overlap vs radius, against the sensor mask;"
                               "mean radius  r  [cm];overlap width  W  [strips]",
                        10, 4.0, 29.0, 10, 0.0, 6.6);
   fi2->GetXaxis()->SetTitleSize(0.045); fi2->GetYaxis()->SetTitleSize(0.045);
   fi2->Draw();

   TLine *ldI = new TLine(5.0, innerDesign, 16.5, innerDesign);
   ldI->SetLineColor(kBlack); ldI->SetLineStyle(2); ldI->SetLineWidth(4); ldI->Draw();
   TLine *ldO = new TLine(16.5, outerDesign, 28.0, outerDesign);
   ldO->SetLineColor(kBlack); ldO->SetLineStyle(2); ldO->SetLineWidth(4); ldO->Draw();
   TLine *lsep = new TLine(16.5, 0.0, 16.5, 6.6);
   lsep->SetLineColor(kGray + 2); lsep->SetLineStyle(3); lsep->SetLineWidth(2); lsep->Draw();

   TGraphErrors *gI = new TGraphErrors(4, rIn, wIn, zero, eIn);
   gI->SetMarkerStyle(21); gI->SetMarkerSize(1.6);
   gI->SetMarkerColor(kRed + 1); gI->SetLineColor(kRed + 1); gI->SetLineWidth(2);
   gI->Draw("P same");
   TGraphErrors *gO = new TGraphErrors(4, rOu, wOu, zero, eOu);
   gO->SetMarkerStyle(24); gO->SetMarkerSize(1.6);
   gO->SetMarkerColor(kBlue + 1); gO->SetLineColor(kBlue + 1); gO->SetLineWidth(2);
   gO->Draw("P same");
   TGraphErrors *gOD = new TGraphErrors(4, rOu, wOuD, zero, eOuD);
   gOD->SetMarkerStyle(20); gOD->SetMarkerSize(1.6);
   gOD->SetMarkerColor(kBlue + 1); gOD->SetLineColor(kBlue + 1); gOD->SetLineWidth(2);
   gOD->Draw("P same");

   // top LEFT: the inner points sit low and the outer points sit right, so
   // the empty corner is above the inner ones.
   TLegend *l2 = new TLegend(0.145, 0.70, 0.78, 0.93);
   l2->SetFillStyle(0); l2->SetBorderSize(0); l2->SetTextSize(0.033);
   l2->AddEntry(gI, Form("inner, measured (%s)", sampleLabel), "pe");
   l2->AddEntry(gOD, "outer, measured (edge strip dropped)", "pe");
   l2->AddEntry(gO, "outer, all d -- the convention spread", "pe");
   l2->AddEntry(ldI, "design, from the GDS mask: 0 inner, 4.267 outer", "l");
   l2->Draw();

   TLatex t2; t2.SetNDC(); t2.SetTextSize(0.030);
   t2.DrawLatex(0.155, 0.285, "the mask is exact: pitch 30/128#circ, fan on the axis to <0.1 #mum,");
   t2.DrawLatex(0.155, 0.240, "inner fan 30.000#circ, outer 15.000#circ, r = 5.000-16.5, 16.5-28.0 cm");
   t2.DrawLatex(0.155, 0.185, "so NEITHER a pitch error NOR an off-axis fan is available:");
   t2.DrawLatex(0.155, 0.140, "every strip of this excess is how the wedges are MOUNTED");

   cvPI->SaveAs("FstOuterSensorOverlap/plots/fwo_pair_inner.png");
   printf("\nwrote FstOuterSensorOverlap/plots/fwo_pair_inner.png\n");
}
