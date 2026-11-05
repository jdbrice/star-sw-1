// fwoRead.C -- read a merged fstWedgeOverlap output and print the derived widths.
//
// Why this exists: fstWedgeOverlap.C accumulates its coincidence counters in plain
// C arrays and used to print them and nothing else, so the only way to get more
// statistics was to run it in one process. The counters are now written as TH1D
// (see "the counters, as histograms" at the end of fstWedgeOverlap.C), so a condor
// run merges with plain hadd and this macro reads the sum.
//
//   hadd -f fwo_merged.root fwo_*.root
//   root4star -l -b -q 'script/fwoRead.C("fwo_merged.root")'
//
// THE FLATTENING ORDER HERE MUST MATCH fstWedgeOverlap.C. Both are written out in
// the loop order documented there; change one and you must change the other.
//
// What the larger sample is for: the PER-BOUNDARY widths (section 3 below). The
// single-process result gave an r-slope of -0.17 strips where pure pitch gives 0.00
// and a pure radial shift d = -0.176 cm gives -0.40, i.e. a mixture. That argument
// assumes the widths are not being faked by a per-wedge rotation pattern, and only
// per-boundary widths can show that: geometry is common to all 36 boundaries, a
// rotation pattern is not.
//
// Width definition, stated because the number depends on it: the profile is
// (predicted cell) minus (mirrored cell) as a function of distance from the strip
// edge; the mirrored cell is the accidental floor. Width = sum over distances
// 0..11 of (profile / plateau), i.e. the effective number of strips that overlap.
//
// Plateau = the mean of distances 1, 2 and 3, NOT including distance 0. The
// profile rises before it falls -- distance 0 is the outermost strip, on the
// physical sensor edge, and it sits below the strips just inside it. Normalising
// on it would inflate every width. Both definitions are printed side by side
// (w13 and w01) so the size of that choice is visible, and the raw profile is
// printed with every width so any other definition can be recomputed without
// re-running the jobs.

double fwoWidth(double *prof, int lo, int hi)
{
   double plateau = 0;
   for (int ia0 = lo; ia0 <= hi; ia0++) plateau += prof[ia0];
   plateau /= (hi - lo + 1);
   if (plateau <= 0) return 0;
   double sum = 0;
   for (int ib0 = 0; ib0 < 12; ib0++) sum += prof[ib0] / plateau;
   return sum;
}

void fwoRead(const Char_t *mergedFile = "fwo_merged.root")
{
   TFile *fin = TFile::Open(mergedFile);
   if (!fin || fin->IsZombie()) { printf("ERROR: cannot open %s\n", mergedFile); return; }

   // ---------- 1. how much data ----------
   TH1D *cTot = (TH1D *) fin->Get("cTot");
   if (!cTot) {
      printf("ERROR: no cTot in %s -- was this written by a version of\n"
             "fstWedgeOverlap.C that persists its counters? Re-run with the current one.\n",
             mergedFile);
      return;
   }
   printf("==== %s ====\n", mergedFile);
   printf("events seen %.0f, with FST %.0f, hits %.0f (outer %.0f), raw %.0f, inner/outer pairs %.0f\n",
          cTot->GetBinContent(1), cTot->GetBinContent(2), cTot->GetBinContent(3),
          cTot->GetBinContent(4), cTot->GetBinContent(5), cTot->GetBinContent(6));

   // ---------- 2. outer width by r-strip ----------
   TH1D *cProfR = (TH1D *) fin->Get("cProfR");
   if (cProfR) {
      double rmid[4] = {17.94, 20.81, 23.69, 26.56};
      double wR[4]   = {0, 0, 0, 0};
      printf("\nOUTER WIDTH BY rStrip -- pitch error or radial shift of the strip fan?\n");
      printf("pure pitch predicts 5.25 in every row; a radial shift d = -0.176 cm\n");
      printf("predicts 5.51, 5.34, 5.21, 5.11 falling with r.\n");
      printf("  rStrip   <r>    w13    w01 |  profile (pred-mirr) by distance from edge 0..11\n");
      for (int ja = 0; ja < 4; ja++) {
         double profR[12];
         for (int jb = 0; jb < 12; jb++) {
            int kb = ja * 24 + jb * 2;                       // [r][dist][pred] -> bin kb+1
            profR[jb] = cProfR->GetBinContent(kb + 1) - cProfR->GetBinContent(kb + 2);
         }
         double sumR = fwoWidth(profR, 1, 3);
         wR[ja] = sumR;
         printf("     %d    %5.2f  %6.2f %6.2f |", ja + 4, rmid[ja], sumR, fwoWidth(profR, 0, 1));
         for (int jd = 0; jd < 12; jd++) printf(" %6.0f", profR[jd]);
         printf("\n");
      }
      // slope by least squares in r-strip index, which is what the hypotheses predict
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      for (int je = 0; je < 4; je++) { sx += je; sy += wR[je]; sxx += je * je; sxy += je * wR[je]; }
      double den = 4 * sxx - sx * sx;
      if (den != 0) {
         double slope = (4 * sxy - sx * sy) / den;
         printf("  slope %+.3f strips per rStrip  (pure pitch 0.00, pure radial shift -0.40)\n", slope);
      }
   }

   // ---------- 3. per-boundary width: geometry or a rotation pattern? ----------
   // Each wedge contributes two ends, so every boundary is measured twice, once
   // from each side. Geometry is common to all of them; a per-wedge rotation is not.
   printf("\nPER-BOUNDARY WIDTH -- the point of the larger sample.\n");
   printf("A common geometry term gives the same width everywhere. A per-wedge\n");
   printf("rotation pattern does not, and it is antisymmetric: a wedge rotated into\n");
   printf("its neighbour widens one of its ends and narrows the other, so the two\n");
   printf("ends of the SAME wedge move in opposite directions while geometry moves\n");
   printf("them together. Watch the (end0 - end1) column, not the widths alone.\n");
   printf(" disk wedge  m  end0   end1   mean   end0-end1\n");
   double wAll = 0; int nAll = 0;
   double adiff = 0; int ndiff = 0;
   for (int jf = 1; jf <= 3; jf++) {
      for (int jg = 1; jg <= 12; jg++) {
         TH2F *hwd = (TH2F *) fin->Get(Form("hD%dW%02d", jf, jg));
         if (!hwd) continue;
         int par = (jg % 2) ? 0 : 1;                       // 0 = odd m, 1 = even m
         int fl  = ((par == 1) ? 1 : 0) ^ ((jf == 2) ? 1 : 0);
         double wEnd[2] = {0, 0};
         for (int jh = 0; jh < 2; jh++) {
            int pr = (jh == 0) ? (fl ? 2 : 0) : (fl ? 1 : 3);
            int mr = (jh == 0) ? (fl ? 0 : 2) : (fl ? 3 : 1);
            double profB[12];
            for (int ji = 0; ji < 12; ji++) {
               int strip = (jh == 0) ? ji : 127 - ji;       // distance from this end
               profB[ji] = hwd->GetBinContent(strip + 1, pr + 1)
                         - hwd->GetBinContent(strip + 1, mr + 1);
            }
            wEnd[jh] = fwoWidth(profB, 1, 3);
         }
         double mean = 0.5 * (wEnd[0] + wEnd[1]);
         printf("   %d   %2d   %2d  %5.2f  %5.2f  %5.2f   %+6.2f\n",
                jf, jg, jg, wEnd[0], wEnd[1], mean, wEnd[0] - wEnd[1]);
         if (wEnd[0] > 0 && wEnd[1] > 0) {
            wAll += mean; nAll++;
            adiff += TMath::Abs(wEnd[0] - wEnd[1]); ndiff++;
         }
      }
   }
   if (nAll > 0) {
      double wmean = wAll / nAll;
      double rms = 0;
      // second pass for the spread, so a wide spread is visible as a number
      int nr = 0;
      for (int jk = 1; jk <= 3; jk++) {
         for (int jl = 1; jl <= 12; jl++) {
            TH2F *hwr = (TH2F *) fin->Get(Form("hD%dW%02d", jk, jl));
            if (!hwr) continue;
            int par2 = (jl % 2) ? 0 : 1;
            int fl2  = ((par2 == 1) ? 1 : 0) ^ ((jk == 2) ? 1 : 0);
            double we[2] = {0, 0};
            for (int jm = 0; jm < 2; jm++) {
               int pr2 = (jm == 0) ? (fl2 ? 2 : 0) : (fl2 ? 1 : 3);
               int mr2 = (jm == 0) ? (fl2 ? 0 : 2) : (fl2 ? 3 : 1);
               double pf[12];
               for (int jn = 0; jn < 12; jn++) {
                  int st2 = (jm == 0) ? jn : 127 - jn;
                  pf[jn] = hwr->GetBinContent(st2 + 1, pr2 + 1)
                         - hwr->GetBinContent(st2 + 1, mr2 + 1);
               }
               we[jm] = fwoWidth(pf, 1, 3);
            }
            if (we[0] > 0 && we[1] > 0) {
               double mn2 = 0.5 * (we[0] + we[1]);
               rms += (mn2 - wmean) * (mn2 - wmean); nr++;
            }
         }
      }
      if (nr > 1) rms = TMath::Sqrt(rms / (nr - 1));
      printf("  %d boundaries: mean width %.2f, rms %.2f strips\n", nAll, wmean, rms);
      printf("  mean |end0-end1| = %.2f strips over %d wedges\n",
             (ndiff > 0) ? adiff / ndiff : 0.0, ndiff);
      printf("  A geometry term shows up as a small rms with |end0-end1| ~ 0.\n");
      printf("  A rotation pattern shows up as |end0-end1| comparable to the rms.\n");
   }

   // ---------- 4. inner vs outer: the gap, via the partnered-edge-hit fraction ----------
   // NOT the profile integral of sections 2 and 3. The inner profile has no
   // plateau to normalise on -- the overlap is about one strip, so the profile
   // falls from distance 0 instead of sitting flat and then falling -- and
   // dividing by a mean over distances 1..3 would be meaningless there. The
   // published inner number is the fraction of edge hits that have a partner in
   // the predicted cell (25.0% inner against 49.3% outer), so that is what this
   // reproduces, from the counters rather than from a shape assumption.
   TH1D *cCoin  = (TH1D *) fin->Get("cCoin");
   TH1D *cCoinI = (TH1D *) fin->Get("cCoinI");
   TH1D *cEdge  = (TH1D *) fin->Get("cEdge");
   TH1D *cEdgeI = (TH1D *) fin->Get("cEdgeI");
   if (cCoin && cCoinI && cEdge && cEdgeI) {
      double parO = 0, parI = 0, edgO = 0, edgI = 0;
      for (int jp = 0; jp < 3; jp++) {                      // disk 0..2
         for (int jq = 0; jq < 2; jq++) {                   // parity
            for (int jr = 0; jr < 2; jr++) {                // end
               int fl3 = ((jq == 1) ? 1 : 0) ^ ((jp == 1) ? 1 : 0);
               int pr3 = (jr == 0) ? (fl3 ? 2 : 0) : (fl3 ? 1 : 3);
               // The inner predicted cell is the outer's WITH THE ENDS SWAPPED:
               // the inner strips run the opposite way inside a wedge, so inner
               // s0 pairs with the neighbour that outer s124-127 pairs with.
               // See the comment at "inner s0 end pairs with" in fstWedgeOverlap.C.
               int pi3 = (jr == 0) ? (fl3 ? 0 : 2) : (fl3 ? 3 : 1);
               int ke  = (jp * 2 + jq) * 2 + jr;                  // [disk][par][end]
               int kc  = ke * 4 + pr3;                            // [disk][par][end][y]
               int ki  = ke * 4 + pi3;
               parO += cCoin ->GetBinContent(kc + 1);
               parI += cCoinI->GetBinContent(ki + 1);
               edgO += cEdge ->GetBinContent(ke + 1);
               edgI += cEdgeI->GetBinContent(ke + 1);
            }
         }
      }
      printf("\nPARTNERED EDGE-HIT FRACTION (the published estimator)\n");
      printf("  outer %6.0f / %6.0f = %5.1f%%   (published 49.3%%)\n",
             parO, edgO, (edgO > 0) ? 100.0 * parO / edgO : 0.0);
      printf("  inner %6.0f / %6.0f = %5.1f%%   (published 25.0%%)\n",
             parI, edgI, (edgI > 0) ? 100.0 * parI / edgI : 0.0);
      printf("  Outer ~5.2 strips, inner ~1.0: their difference is the 1 deg gap in\n");
      printf("  units of pitch, nominal 4.27, which is what pins the gap as the\n");
      printf("  innocent one and puts the excess in the strips.\n");
   }

   // The inner profile, raw, for the record -- no width derived from it on purpose.
   TH1D *cProfI = (TH1D *) fin->Get("cProfI");
   if (cProfI) {
      printf("  inner profile (pred-mirr) by distance 0..11 :");
      for (int js = 0; js < 12; js++)
         printf(" %6.0f", cProfI->GetBinContent(js * 2 + 1) - cProfI->GetBinContent(js * 2 + 2));
      printf("\n");
   }

   fin->Close();
}
