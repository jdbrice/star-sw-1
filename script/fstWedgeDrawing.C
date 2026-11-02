// fstWedgeDrawing.C
//
// Draw three adjacent FST wedges with their 9 sensors, the 1 degree gap on each
// wedge centreline, the 0.5 degree overhangs that make the 1 degree overlap at
// each wedge boundary, and the strip numbering and its direction.
//
// Default is disk 1, wedges 12, 1 and 2: the 60-90 deg wedge and its two
// neighbours.
//
// NUMBERING IS 1-BASED THROUGHOUT, because that is what the code and the data
// use: StMuFstHit::getWedge() returns 1-36, getDisk() returns 1-3, and
// StFstDb::stripPhi rejects wedge < 1. Within a disk, m = 1 + (wedge-1)%12,
// so m runs 1-12 and wedge 1 is the 60-90 deg one. The 0-based iwedgem =
// (wedge-1)%12 exists only to index kFstzDirct etc. Talking 0-based flips the
// odd/even parity of m, which is the whole front/back grouping, so do not.
//
// Every phi here is computed with the SAME formula as StFstDb::stripPhi (the
// constants are inlined because StFstConsts.h is often invisible to CINT), so
// the picture is that function's claim, not a redrawing of someone's idea of
// it. Changing stripPhi without changing this macro will make the two disagree.
//
// VIEW: from the IR (z = 0) looking towards +z, which is what you need for
// "which face points at the IR". In that view STAR +x points LEFT, +y is up,
// and STAR phi increases CLOCKWISE. Panel A says so on the picture, because
// this is exactly the handedness that makes mirror questions confusing.
//
// Usage: root4star -l -b -q 'fstWedgeDrawing.C(1,"fstWedgeDrawing.png")'

double fwdPitch()  { return 30.0 / 128.0; }        // deg, kFstStripPitchPhi
double fwdGap()    { return 1.0; }                 // deg, kFstStripGapPhi

// StFstDb::stripPhi, in degrees, for disk 1/3 and disk 2
double fwdStripPhi(int disk, int m, int sensor, int rStrip, double phiStrip)
{
   double phiStart[12] = {2.0, 2.0, 0.0, 12.0, 10.0, 10.0, 8.0, 8.0, 6.0, 6.0, 4.0, 4.0};
   double phiStop [12] = {3.0, 1.0, 1.0, 11.0, 11.0,  9.0, 9.0, 7.0, 7.0, 5.0, 5.0, 3.0};
   int    zDirct  [12] = {1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1};
   int    zFilp   [3]  = {1, -1, 1};

   double pitch = fwdPitch();
   double dirct = zDirct[m - 1];
   double sgn   = zFilp[disk - 1] * dirct;
   double phiInner, phiOuter;

   if (2 == disk) {
      phiInner = phiStop [m - 1] * 30.0 - 0.5 * dirct * pitch;
      phiOuter = phiStart[m - 1] * 30.0 + 0.5 * dirct * pitch;
   } else {
      phiInner = phiStart[m - 1] * 30.0 + 0.5 * dirct * pitch;
      phiOuter = phiStop [m - 1] * 30.0 - 0.5 * dirct * pitch;
   }

   if (rStrip < 4) return phiInner + sgn * phiStrip * pitch;

   double gap = (1 == sensor) ?  sgn * 0.5 * fwdGap()
                             : ((2 == sensor) ? -sgn * 0.5 * fwdGap() : 0.0);
   return phiOuter - sgn * phiStrip * pitch + gap;
}

int fwdOrient(int disk, int m)          // +1 = front (larger z, away from the IR)
{
   int zDirct[12] = {1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1};
   int zFilp [3]  = {1, -1, 1};
   return zFilp[disk - 1] * zDirct[m - 1];      // kFstmAsInstalled, the default
}

void fstWedgeDrawing(int disk = 1, const Char_t *out = "fstWedgeDrawing.png")
{
   const double rIn1 = 5.000, rIn2 = 13.625, rOut2 = 28.000;   // kFstrStart/Stop
   const double zDisk[3] = {151.750, 165.248, 178.781};
   const double kInnerZ = 1.7419, kOuterZ = 0.3469;
   int mm[3] = {12, 1, 2};                     // code numbering, phi descending
   const Char_t *lbl[3] = {"wedge 12", "wedge 1", "wedge 2"};

   gStyle->SetOptStat(0);
   TCanvas *c = new TCanvas("c", "FST wedge geometry", 1700, 1150);

   // ================= panel A : true-scale IR view =================
   c->cd();
   TPad *pA = new TPad("pA", "", 0.00, 0.50, 0.33, 1.00); pA->Draw();
   pA->cd(); pA->Range(-35.0, -35.9, 35.0, 35.9);
   TLatex *tA = new TLatex(); tA->SetTextSize(0.038);
   for (int ia0 = 1; ia0 <= 12; ia0++) {
      int isThree = (ia0 == 12 || ia0 == 1 || ia0 == 2);
      double p0 = fwdStripPhi(disk, ia0, 1, 7, 0.0);
      double p1 = fwdStripPhi(disk, ia0, 2, 7, 127.0);
      double lo = (p0 < p1) ? p0 : p1, hi = (p0 < p1) ? p1 : p0;
      double xs[62], ys[62];
      int np = 0;
      for (int ib = 0; ib <= 20; ib++) {
         double ph = (lo + (hi - lo) * ib / 20.0) * TMath::DegToRad();
         xs[np] = -rOut2 * TMath::Cos(ph); ys[np] = rOut2 * TMath::Sin(ph); np++;
      }
      for (int ic = 20; ic >= 0; ic--) {
         double ph = (lo + (hi - lo) * ic / 20.0) * TMath::DegToRad();
         xs[np] = -rIn1 * TMath::Cos(ph); ys[np] = rIn1 * TMath::Sin(ph); np++;
      }
      xs[np] = xs[0]; ys[np] = ys[0]; np++;
      TPolyLine *pl = new TPolyLine(np, xs, ys);
      if (isThree) pl->SetFillColor((fwdOrient(disk, ia0) > 0) ? kOrange - 4 : kAzure - 9);
      else         pl->SetFillColor((fwdOrient(disk, ia0) > 0) ? 19 : 18);
      pl->SetLineColor(kBlack); pl->SetLineWidth(isThree ? 3 : 1);
      pl->Draw("f"); pl->Draw();
      double pmid = 0.5 * (lo + hi) * TMath::DegToRad();
      tA->SetTextAlign(22); tA->SetTextSize(isThree ? 0.045 : 0.032);
      tA->DrawLatex(-21.0 * TMath::Cos(pmid), 21.0 * TMath::Sin(pmid),
                    Form("%d", ia0));
      if (isThree)
         tA->DrawLatex(-26.0 * TMath::Cos(pmid), 26.0 * TMath::Sin(pmid),
                       (fwdOrient(disk, ia0) > 0) ? "FRONT" : "BACK");
   }
   tA->SetTextAlign(12); tA->SetTextSize(0.036);
   tA->DrawLatex(-34.5, 34.0, Form("A. disk %d (of disks 1-3) from the IR towards +z", disk));
   tA->SetTextSize(0.030);
   tA->DrawLatex(-34.5, 31.3, "wedge 1-12 as in the code, true scale");
   tA->DrawLatex(-34.5, -29.5, "STAR +x points LEFT, +y up:");
   tA->DrawLatex(-34.5, -32.0, "#phi increases CLOCKWISE in this view");
   tA->DrawLatex(-34.5, -34.5, "FRONT = larger z = away from the IR");
   // the two boundaries the other panels talk about
   tA->SetTextAlign(21); tA->SetTextSize(0.032);
   tA->DrawLatex(0.0, 28.9, "#phi=90#circ  (wedge 1 | wedge 12)");
   tA->SetTextAlign(32);
   tA->DrawLatex(-29.2, 0.0, "#phi=0");
   TArrow *axr = new TArrow(0, 0, -9, 0, 0.025, "|>"); axr->SetLineWidth(2); axr->Draw();
   tA->SetTextAlign(22); tA->SetTextSize(0.038);
   tA->DrawLatex(-10.5, 2.2, "+x");
   TArrow *axu = new TArrow(0, 0, 0, 9, 0.025, "|>"); axu->SetLineWidth(2); axu->Draw();
   tA->DrawLatex(2.6, 9.5, "+y");
   // phi runs clockwise on screen: draw it outside the disk, sweeping that way,
   // with the head at the END of the sweep so the sense cannot be misread
   double arcX[21], arcY[21];
   for (int ip0 = 0; ip0 <= 20; ip0++) {
      double aa = (-16.0 - 52.0 * ip0 / 20.0) * TMath::DegToRad();   // -16 -> -68 deg
      arcX[ip0] = 31.5 * TMath::Cos(aa); arcY[ip0] = 31.5 * TMath::Sin(aa);
   }
   TPolyLine *parc = new TPolyLine(21, arcX, arcY);
   parc->SetLineColor(kRed + 1); parc->SetLineWidth(3); parc->Draw();
   TArrow *ahead = new TArrow(arcX[18], arcY[18], arcX[20], arcY[20], 0.030, "|>");
   ahead->SetLineColor(kRed + 1); ahead->SetFillColor(kRed + 1); ahead->SetLineWidth(3); ahead->Draw();
   tA->SetTextColor(kRed + 1); tA->SetTextAlign(22); tA->SetTextSize(0.030);
   tA->DrawLatex(20.5, -20.5, "STAR #phi");
   tA->DrawLatex(20.5, -23.3, "increases");
   tA->SetTextColor(kBlack);

   c->cd();
   TPad *pB = new TPad("pB", "", 0.33, 0.50, 1.00, 1.00); pB->Draw();
   pB->cd(); pB->SetLeftMargin(0.07); pB->SetBottomMargin(0.11);
   TH2F *fr = new TH2F("fr", Form("B. disk %d (of disks 1-3) unrolled: 3 wedges, 9 sensors (phi exaggerated);STAR #phi [deg];r [cm]", disk),
                       100, 26, 124, 100, 2, 31);
   fr->Draw();
   TLatex *tB = new TLatex(); tB->SetTextSize(0.028); tB->SetTextAlign(22);
   for (int id = 0; id < 3; id++) {
      int m = mm[id];
      int orient = fwdOrient(disk, m);
      // wedge 0 gets its own colour and its own height, so the three wedges'
      // arrows and strip numbers never sit on top of each other. Same colours
      // as panel D: orange = wedge 0, blue = wedges 11 and 1.
      int isW0  = (m == 1);
      int colW  = isW0 ? (kOrange + 7) : (kAzure + 2);
      double yAo = isW0 ? 25.6 : 22.9;      // outer arrow
      double yLo2 = isW0 ? 24.4 : 21.7;     // outer strip labels
      double yAi = isW0 ? 11.8 :  8.8;      // inner arrow
      double yLi = isW0 ? 10.6 :  7.6;      // inner strip labels
      const double kIn = 1.6;               // push end labels inside the sensor

      // inner sensor: strips 0 and 127 at rStrip 0-3
      double i0 = fwdStripPhi(disk, m, 0, 0, 0.0);
      double i1 = fwdStripPhi(disk, m, 0, 0, 127.0);
      double iLo = (i0 < i1) ? i0 : i1, iHi = (i0 < i1) ? i1 : i0;
      TBox *bi = new TBox(iLo - 0.5 * fwdPitch(), rIn1, iHi + 0.5 * fwdPitch(), rIn2);
      bi->SetFillColor(kGray + 1); bi->SetLineColor(kBlack); bi->Draw("l");
      // the two outer half sensors
      for (int ie = 1; ie <= 2; ie++) {
         double sLo = (ie == 1) ? 0.0 : 64.0, sHi = (ie == 1) ? 63.0 : 127.0;
         double o0 = fwdStripPhi(disk, m, ie, 7, sLo);
         double o1 = fwdStripPhi(disk, m, ie, 7, sHi);
         double oLo = (o0 < o1) ? o0 : o1, oHi = (o0 < o1) ? o1 : o0;
         TBox *bo = new TBox(oLo - 0.5 * fwdPitch(), rIn2, oHi + 0.5 * fwdPitch(), rOut2);
         bo->SetFillColor((orient > 0) ? kOrange - 4 : kAzure - 9);
         bo->SetLineColor(kBlack); bo->Draw("l");
         tB->SetTextColor(kBlack); tB->SetTextSize(0.026);
         tB->DrawLatex(0.5 * (oLo + oHi), 27.3, Form("outer %d", ie));
         TArrow *ar = new TArrow(o0, yAo, o1, yAo, 0.012, "|>");
         ar->SetLineWidth(2); ar->SetLineColor(colW); ar->SetFillColor(colW); ar->Draw();
         // only the two extreme ends are labelled, pushed inside the sensor so
         // nothing is clipped by the frame; 63 and 64 meet across the gap and
         // are called out in the legend instead
         tB->SetTextColor(colW); tB->SetTextSize(0.028);
         if ((int)sLo == 0)   tB->DrawLatex(o0 + ((o0 > o1) ? -kIn : kIn), yLo2, "s0");
         if ((int)sHi == 127) tB->DrawLatex(o1 + ((o1 > o0) ? -kIn : kIn), yLo2, "s127");
         tB->SetTextColor(kBlack);
      }
      // inner strip direction
      TArrow *ai = new TArrow(i0, yAi, i1, yAi, 0.012, "|>");
      ai->SetLineWidth(2); ai->SetLineColor(colW); ai->SetFillColor(colW); ai->Draw();
      tB->SetTextColor(colW); tB->SetTextSize(0.028);
      tB->DrawLatex(i0 + ((i0 > i1) ? -kIn : kIn), yLi, "s0");
      tB->DrawLatex(i1 + ((i1 > i0) ? -kIn : kIn), yLi, "s127");
      tB->SetTextColor(kBlack); tB->SetTextSize(0.026);
      tB->DrawLatex(0.5 * (iLo + iHi), 12.9, "inner");
      tB->SetTextColor(colW); tB->SetTextSize(0.030);
      tB->DrawLatex(0.5 * (iLo + iHi), 4.0, Form("%s  %s", lbl[id], (orient > 0) ? "FRONT" : "BACK"));
      tB->SetTextColor(kBlack); tB->SetTextSize(0.026);
      // the 1 deg gap on the centreline
      double gLo = fwdStripPhi(disk, m, 1, 7, 63.0), gHi = fwdStripPhi(disk, m, 2, 7, 64.0);
      double gA = (gLo < gHi) ? gLo : gHi, gB = (gLo < gHi) ? gHi : gLo;
      TBox *bg = new TBox(gA + 0.5 * fwdPitch(), rIn2, gB - 0.5 * fwdPitch(), rOut2);
      bg->SetFillColor(kWhite); bg->SetLineColor(kRed); bg->SetLineStyle(2); bg->Draw("l");
      tB->SetTextColor(kRed);
      tB->DrawLatex(0.5 * (gA + gB), 29.3, "1#circ gap");
      tB->SetTextColor(kBlack);
   }
   // the overlap bands at the two wedge boundaries
   for (int ig = 0; ig < 2; ig++) {
      double bnd = (ig == 0) ? 90.0 : 60.0;
      TBox *bb = new TBox(bnd - 0.5, rIn2, bnd + 0.5, rOut2);
      bb->SetFillStyle(3354); bb->SetFillColor(kGreen + 2); bb->Draw();
      tB->SetTextColor(kGreen + 3);
      tB->DrawLatex(bnd, 29.3, "1#circ overlap");
      tB->SetTextColor(kBlack);
   }
   tB->SetTextAlign(12); tB->SetTextSize(0.024);
   tB->DrawLatex(26.8, 19.2, "arrow = direction of");
   tB->DrawLatex(26.8, 17.9, "increasing strip index");
   tB->DrawLatex(26.8, 16.6, "s63 and s64 meet");
   tB->DrawLatex(26.8, 15.3, "across each 1#circ gap");
   tB->SetTextAlign(22);

   // ================= panel C : z cross section =================
   c->cd();
   TPad *pC = new TPad("pC", "", 0.00, 0.00, 0.40, 0.50); pC->Draw();
   pC->cd(); pC->SetLeftMargin(0.11); pC->SetBottomMargin(0.12);
   double z0 = zDisk[disk - 1];
   TH2F *fz = new TH2F("fz", Form("C. z section, disk %d (of disks 1-3);z [cm];r [cm]", disk),
                       100, z0 - 3.2, z0 + 3.2, 100, 2, 31);
   fz->Draw();
   TLatex *tC = new TLatex(); tC->SetTextSize(0.034);
   TLine *ld = new TLine(z0, 2, z0, 31); ld->SetLineStyle(2); ld->SetLineColor(kGray + 2); ld->Draw();
   tC->SetTextAlign(21); tC->DrawLatex(z0, 29.6, Form("disk nominal z = %.3f", z0));
   for (int ih = 0; ih < 2; ih++) {
      int orient = (ih == 0) ? +1 : -1;            // front, then back
      int col = (orient > 0) ? kOrange - 4 : kAzure - 9;
      double zi = z0 + orient * kInnerZ, zo = z0 + orient * kOuterZ;
      TBox *bi2 = new TBox(zi - 0.07, rIn1, zi + 0.07, rIn2);
      bi2->SetFillColor(kGray + 1); bi2->Draw();
      TBox *bo2 = new TBox(zo - 0.07, rIn2, zo + 0.07, rOut2);
      bo2->SetFillColor(col); bo2->Draw();
      tC->SetTextSize(0.030);
      tC->SetTextAlign((orient > 0) ? 12 : 32);
      double tx = zi + ((orient > 0) ? 0.18 : -0.18);
      tC->DrawLatex(tx, 9.2, Form("inner z = %.3f", zi));
      tC->DrawLatex(zo + ((orient > 0) ? 0.18 : -0.18), 26.0, Form("outer z = %.3f", zo));
      tC->SetTextAlign((orient > 0) ? 12 : 32); tC->SetTextSize(0.034);
      tC->DrawLatex(zo + ((orient > 0) ? 0.20 : -0.20), 17.0, (orient > 0) ? "FRONT" : "BACK");
      tC->SetTextSize(0.026);
      tC->DrawLatex(zo + ((orient > 0) ? 0.20 : -0.20), 15.3,
                    (orient > 0) ? "wedge 1" : "wedges 12, 2");
   }
   tC->SetTextAlign(12); tC->SetTextSize(0.026);
   tC->DrawLatex(z0 - 3.0, 23.0, "inner silicon sits 1.7419 cm");
   tC->DrawLatex(z0 - 3.0, 21.6, "from the disk plane, outer");
   tC->DrawLatex(z0 - 3.0, 20.2, "0.3469 cm, SAME side:");
   tC->DrawLatex(z0 - 3.0, 18.8, "the wedge flip moves both.");
   TArrow *air = new TArrow(z0 - 1.0, 4.0, z0 - 2.9, 4.0, 0.025, "|>");
   air->SetLineWidth(3); air->Draw();
   tC->SetTextAlign(12); tC->SetTextSize(0.034);
   tC->DrawLatex(z0 - 2.8, 5.3, "towards the IR");

   // ================= panel D : zoom on the 90 deg boundary =================
   c->cd();
   TPad *pD = new TPad("pD", "", 0.40, 0.00, 1.00, 0.50); pD->Draw();
   pD->cd(); pD->SetLeftMargin(0.07); pD->SetBottomMargin(0.12);
   TH2F *fd = new TH2F("fd", "D. zoom on the wedge 12 / wedge 1 boundary: the 4 overlapping outer strips;STAR #phi [deg];",
                       100, 89.3, 90.7, 100, 0, 10);
   fd->GetYaxis()->SetLabelSize(0); fd->GetYaxis()->SetTickLength(0);
   fd->Draw();
   TLatex *tD = new TLatex(); tD->SetTextSize(0.036); tD->SetTextAlign(22);
   TBox *bov = new TBox(89.5, 0, 90.5, 10);
   bov->SetFillStyle(3003); bov->SetFillColor(kGreen + 1); bov->Draw();
   for (int ie2 = 0; ie2 < 2; ie2++) {
      TLine *lb = new TLine(89.5 + ie2, 0, 89.5 + ie2, 10);
      lb->SetLineColor(kGreen + 3); lb->SetLineStyle(2); lb->SetLineWidth(2); lb->Draw();
   }
   // wedge 0 (m=1) outer strips 0-5 and wedge 11 (m=12) outer strips 0-5
   for (int ii = 0; ii <= 5; ii++) {
      double pw0 = fwdStripPhi(disk, 1, 1, 7, (double)ii);
      TLine *l0 = new TLine(pw0, 5.4, pw0, 8.2); l0->SetLineColor(kOrange + 7); l0->SetLineWidth(3); l0->Draw();
      tD->SetTextColor(kOrange + 7); tD->DrawLatex(pw0, 8.8, Form("%d", ii));
      double pw11 = fwdStripPhi(disk, 12, 1, 7, (double)ii);
      TLine *l1 = new TLine(pw11, 1.8, pw11, 4.6); l1->SetLineColor(kAzure + 2); l1->SetLineWidth(3); l1->Draw();
      tD->SetTextColor(kAzure + 2); tD->DrawLatex(pw11, 1.2, Form("%d", ii));
   }
   tD->SetTextColor(kOrange + 7); tD->SetTextAlign(12); tD->SetTextSize(0.034);
   tD->DrawLatex(89.35, 9.5, "wedge 1 (60-90), FRONT: outer strip index grows towards LOWER phi");
   tD->SetTextColor(kAzure + 2);
   tD->DrawLatex(89.35, 0.5, "wedge 12 (90-120), BACK: outer strip index grows towards HIGHER phi");
   tD->SetTextColor(kGreen + 3); tD->SetTextAlign(22);
   tD->DrawLatex(90.0, 10.0 - 0.0, "");
   tD->SetTextColor(kBlack); tD->SetTextAlign(22); tD->SetTextSize(0.030);
   tD->DrawLatex(90.0, 6.9, "1#circ overlap, 4.27 pitches: both wedges have strips 0-3 here");
   tD->SetTextSize(0.027);
   tD->DrawLatex(90.0, 3.4, "so w1 s0 pairs with w12 s3, s1 with s2 ... (offset 1.0 - 4#times0.234 = 0.0625#circ)");

   c->cd();
   c->SaveAs(out);
   printf("\nwrote %s\n", out);

   // the numbers the picture is drawn from, so they can be checked by eye
   printf("\ndisk %d, as drawn (deg):\n", disk);
   for (int ij = 0; ij < 3; ij++) {
      int m = mm[ij];
      printf("  %-9s (m=%2d) %-5s : inner s0 %7.3f  s127 %7.3f | outer s0 %7.3f  s127 %7.3f\n",
             lbl[ij], m, (fwdOrient(disk, m) > 0) ? "FRONT" : "BACK",
             fwdStripPhi(disk, m, 0, 0, 0.0),   fwdStripPhi(disk, m, 0, 0, 127.0),
             fwdStripPhi(disk, m, 1, 7, 0.0),   fwdStripPhi(disk, m, 2, 7, 127.0));
   }
}
