// fttDataWindowPlot.C -- look at what the window fit saw, from the saved histograms.
//
// Pure analysis: it reads fttDataWindow.<run>.txt and root/fttDataWindow.<run>.root
// and never opens a MuDst. The fit it draws is fitFttWindow() from fttDataWindow.C,
// loaded from there rather than copied, so the picture cannot drift from the code.
//
//   root4star -l -b -q 'fttDataWindowPlot.C("fttDataWindow",23081049)'
//
//   dir        holds fttDataWindow.<run>.txt and root/fttDataWindow.<run>.root
//   outdir     where the png files go
//   vmmA, anchorA, anchorB   optional: draw one VMM twice, against two anchors, to
//              show what an anchor error does. vmmA < 0 skips it.
//   vmmChan    optional: a VMM whose per-channel dbcid was kept by fttDataWindow.C's
//              vmmQA argument. Answers "is this VMM shifted as a whole, or only some
//              of its channels?". vmmChan < 0 skips it.
//   chanAnchor the reference to measure that per-channel view against, normally the
//              anchor of the VMM's FOB mates.
//
// Figures written:
//   fttWindow_floor_<run>.png    the busiest VMMs: the fit asks for less than the
//                                floor and the floor holds the window open
//   fttWindow_wider_<run>.png    VMMs whose own crossing reaches outside the floor,
//                                if the run has any (after the half-period cap and
//                                the plain-argmax anchor, usually it has none)
//   fttWindow_anchor_<run>.png   the two-anchor comparison
//   fttWindow_chan_<run>.png     dbcid vs channel for one VMM
//
// CINT notes: fixed-size arrays only, unique loop variable names.

const char *kDWMacro = "/star/u/akio/fcstrk11/star-sw-fwd/StRoot/StFttDbMaker/macro/FttDataWindows/fttDataWindow.C";
const int kPlotHalf = 14;     // ticks either side of the anchor to draw

void onePad(TFile *F, int vmm, int anchor, int storedLo, int storedHi,
            double bxPeriod, double winKeep, int floorLo, int floorHi,
            int nHitsTotal, const char *why, double inTimeHalf)
{
  TH1 *h = (TH1*) F->Get(Form("hDbcid_%03d", vmm));
  if (!h) { printf("no histogram for VMM %d\n", vmm); return; }
  if (!fitFttWindow(h, anchor, bxPeriod, winKeep, inTimeHalf)) { printf("VMM %d did not fit\n", vmm); return; }
  int fitLo = gFitLo, fitHi = gFitHi;

  TH1F *hc = new TH1F(Form("hc%d_%d", vmm, anchor), "", 2 * kPlotHalf + 1, -kPlotHalf - 0.5, kPlotHalf + 0.5);
  TH1F *hb = new TH1F(Form("hb%d_%d", vmm, anchor), "", 2 * kPlotHalf + 1, -kPlotHalf - 0.5, kPlotHalf + 0.5);
  double tot = 0;
  for (int it = -kPlotHalf; it <= kPlotHalf; it++) {
    double c = h->GetBinContent(((anchor + it) % 4096 + 4096) % 4096 + 1);
    hc->SetBinContent(it + kPlotHalf + 1, c);
    hb->SetBinContent(it + kPlotHalf + 1, gFitBkg[it + 40]);
    tot += c;
  }
  // a semicolon in the title would be read by ROOT as an axis-title separator
  hc->SetTitle(Form("VMM %d -- %s;t = dbcid - anchor  [ticks];hits", vmm, why));
  hc->SetLineColor(kBlack); hc->SetLineWidth(2);
  hc->SetMaximum(hc->GetMaximum() * 1.35); hc->SetMinimum(0);
  hc->Draw("hist");
  hb->SetLineColor(kRed); hb->SetLineStyle(2); hb->SetLineWidth(2); hb->Draw("hist same");

  double y0 = 0, y1 = hc->GetMaximum();
  TBox *bfl = new TBox(floorLo - 0.5, y0, floorHi + 0.5, y1 * 0.97);
  bfl->SetFillColor(17); bfl->SetFillStyle(1001); bfl->Draw();
  hc->Draw("hist same"); hb->Draw("hist same");
  TBox *bf = new TBox(fitLo - 0.5, y0, fitHi + 0.5, y1 * 0.97);
  bf->SetFillColorAlpha(kBlue, 0.12); bf->SetLineColor(kBlue); bf->SetLineWidth(2); bf->SetFillStyle(1001);
  bf->Draw("l");
  TLine *ls1 = new TLine(storedLo - 0.5, y0, storedLo - 0.5, y1 * 0.97);
  TLine *ls2 = new TLine(storedHi + 0.5, y0, storedHi + 0.5, y1 * 0.97);
  ls1->SetLineColor(kGreen + 2); ls1->SetLineWidth(3); ls1->Draw();
  ls2->SetLineColor(kGreen + 2); ls2->SetLineWidth(3); ls2->Draw();

  TLatex *tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.042);
  tx->DrawLatex(0.15, 0.88, Form("%d hits (%.0f in view)", nHitsTotal, tot));
  tx->SetTextColor(kBlue);      tx->DrawLatex(0.15, 0.82, Form("fit    [%d, %d]", fitLo, fitHi));
  tx->SetTextColor(kGreen + 2); tx->DrawLatex(0.15, 0.76, Form("stored [%d, %d]", storedLo, storedHi));
  tx->SetTextColor(14);         tx->DrawLatex(0.15, 0.70, Form("floor  [%d, %d]", floorLo, floorHi));
  tx->SetTextColor(kRed);       tx->DrawLatex(0.15, 0.64, "pileup background");
  tx->SetTextColor(kBlack);
}

void fttDataWindowPlot(const char *dir, int run,
                       const char *outdir = "png",
                       int vmmA = -1, int anchorA = 0, int anchorB = 0,
                       int vmmChan = -1, int chanAnchor = 0,
                       double bxPeriod = 4.263, double winKeep = 0.95,
                       int floorLo = -1, int floorHi = 3)
{
  gROOT->LoadMacro(kDWMacro);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05, "t");

  TFile *F = TFile::Open(Form("%s/root/fttDataWindow.%d.root", dir, run));
  if (!F || F->IsZombie()) { printf("no root file in %s\n", dir); return; }

  int anc[400], slo[400], shi[400], nh[400];
  for (int iz = 0; iz < 400; iz++) { anc[iz] = -1; slo[iz] = 0; shi[iz] = 0; nh[iz] = 0; }
  char line[512];
  FILE *I = fopen(Form("%s/fttDataWindow.%d.txt", dir, run), "r");
  if (!I) { printf("no calibration in %s\n", dir); return; }
  while (fgets(line, sizeof(line), I)) {
    if (line[0] == '#') continue;
    // CINT's sscanf takes at most 12 arguments, so read only what is needed here
    int u, m, lo, hi, a, st; double nn;
    if (sscanf(line, "%d %d %d %d %d %d %lf", &u, &m, &lo, &hi, &a, &st, &nn) < 7) continue;
    if (u < 0 || u >= 384) continue;
    anc[u] = a; slo[u] = lo; shi[u] = hi; nh[u] = (int) nn;
  }
  fclose(I);
  gSystem->mkdir(outdir, kTRUE);

  // ---- VMMs whose own crossing reaches outside the floor -------------------------
  int wide[3], nw = 0;
  for (int iw = 0; iw < 3; iw++) wide[iw] = -1;
  for (int jv = 0; jv < 384 && nw < 3; jv++)
    if (anc[jv] >= 0 && (slo[jv] < floorLo || shi[jv] > floorHi)) { wide[nw] = jv; nw++; }
  if (nw > 0) {
    TCanvas *c1 = new TCanvas("cwide", "wider than the floor", 1500, 420);
    c1->Divide(3, 1);
    for (int pw = 0; pw < nw; pw++) {
      c1->cd(pw + 1); gPad->SetLeftMargin(0.12);
      onePad(F, wide[pw], anc[wide[pw]], slo[wide[pw]], shi[wide[pw]],
             bxPeriod, winKeep, floorLo, floorHi, nh[wide[pw]], "the fit opens it wider", 0);
    }
    c1->SaveAs(Form("%s/fttWindow_wider_%d.png", outdir, run));
  } else {
    printf("fttDataWindowPlot: no VMM is wider than the floor in this run\n");
  }

  // ---- the busiest VMMs sitting exactly on the floor ------------------------------
  int tight[3];
  for (int kv = 0; kv < 3; kv++) tight[kv] = -1;
  for (int lv = 0; lv < 384; lv++) {
    if (anc[lv] < 0 || slo[lv] != floorLo || shi[lv] != floorHi) continue;
    for (int mv = 0; mv < 3; mv++) {
      if (tight[mv] < 0) { tight[mv] = lv; break; }
      if (nh[lv] > nh[tight[mv]]) {
        for (int nv2 = 2; nv2 > mv; nv2--) tight[nv2] = tight[nv2 - 1];
        tight[mv] = lv; break;
      }
    }
  }
  TCanvas *c2 = new TCanvas("ctight", "the floor holds it open", 1500, 420);
  c2->Divide(3, 1);
  for (int pt = 0; pt < 3; pt++) {
    if (tight[pt] < 0) continue;
    c2->cd(pt + 1); gPad->SetLeftMargin(0.12);
    onePad(F, tight[pt], anc[tight[pt]], slo[tight[pt]], shi[tight[pt]],
           bxPeriod, winKeep, floorLo, floorHi, nh[tight[pt]], "the fit asks for less, the floor holds", 0);
  }
  c2->SaveAs(Form("%s/fttWindow_floor_%d.png", outdir, run));

  // ---- one VMM against two anchors ------------------------------------------------
  if (vmmA >= 0) {
    TCanvas *c3 = new TCanvas("canch", "two anchors", 1500, 420);
    c3->Divide(3, 1);
    int mate = -1;
    for (int im = 4 * (vmmA / 4); im < 4 * (vmmA / 4) + 4; im++)
      if (im != vmmA && anc[im] >= 0 && (mate < 0 || nh[im] > nh[mate])) mate = im;
    if (mate >= 0) {
      c3->cd(1); gPad->SetLeftMargin(0.12);
      onePad(F, mate, anc[mate], slo[mate], shi[mate], bxPeriod, winKeep,
             floorLo, floorHi, nh[mate], "healthy VMM on the same FOB", 0);
    }
    c3->cd(2); gPad->SetLeftMargin(0.12);
    onePad(F, vmmA, anchorA, floorLo, floorHi, bxPeriod, winKeep,
           floorLo, floorHi, nh[vmmA], Form("anchor %d", anchorA), 0);
    c3->cd(3); gPad->SetLeftMargin(0.12);
    onePad(F, vmmA, anchorB, slo[vmmA], shi[vmmA], bxPeriod, winKeep,
           floorLo, floorHi, nh[vmmA], Form("same data, anchor %d", anchorB), 0);
    c3->SaveAs(Form("%s/fttWindow_anchor_%d.png", outdir, run));
  }

  // ---- dbcid vs channel for one VMM, if it was kept -------------------------------
  if (vmmChan >= 0) {
    TH2 *h2 = (TH2*) F->Get(Form("hChDbcid_%03d", vmmChan));
    if (!h2) {
      printf("fttDataWindowPlot: no per-channel histogram for VMM %d "
             "(rerun fttDataWindow.C with vmmQA=\"%d\")\n", vmmChan, vmmChan);
    } else {
      int ref = (chanAnchor > 0) ? chanAnchor : anc[vmmChan];
      TH2F *hz = new TH2F("hzchan", Form("VMM %d, reference dbcid %d;channel;t = dbcid - reference",
                                         vmmChan, ref), 64, -0.5, 63.5, 41, -20.5, 20.5);
      for (int ic = 0; ic < 64; ic++)
        for (int it2 = -20; it2 <= 20; it2++)
          hz->SetBinContent(ic + 1, it2 + 21,
                            h2->GetBinContent(ic + 1, ((ref + it2) % 4096 + 4096) % 4096 + 1));
      TCanvas *c4 = new TCanvas("cchan", "dbcid vs channel", 800, 600);
      gPad->SetRightMargin(0.13);
      hz->Draw("colz");
      c4->SaveAs(Form("%s/fttWindow_chan_%d.png", outdir, run));
    }
  }
  printf("fttDataWindowPlot: wrote figures into %s\n", outdir);
}
