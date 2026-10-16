// fttDataWindowRefit.C -- change the stored time WINDOW of an existing calibration
// without re-reading any MuDst.
//
// fttDataWindow.C saves every VMM's raw dbcid histogram in
// root/fttDataWindow.<run>.root. The window is measured from exactly those
// histograms, so a different window criterion is a re-analysis of the saved file,
// not a new pass over the data. The anchors are NOT recomputed: they are read back
// from the existing calibration, so only columns 3 and 4 (min,max) can change.
//
// Usage:
//   root4star -l -b -q 'fttDataWindowRefit.C(23081049)'                  // same settings, round trip
//   root4star -l -b -q 'fttDataWindowRefit.C(23081049,"out",0.90)'       // keep 90% of the crossing
//   root4star -l -b -q 'fttDataWindowRefit.C(23081049,"out",0.95,4.263,-2,4)'  // a wider floor
//
//   run        run number
//   outdir     where to write the new fttDataWindow.<run>.txt. Default "refit".
//              The input is read from indir (default "."), which must hold both
//              fttDataWindow.<run>.txt and root/fttDataWindow.<run>.root.
//   winKeep    fraction of the in-time crossing the window must hold (0.95)
//   bxPeriod   ticks per bunch crossing (4.263)
//   winFloorMin, winFloorMax   the window is always widened to contain these (-1, 3)
//
// The algorithm itself lives in fttDataWindow.C (fitFttWindow) and is loaded from
// there, so there is one copy of it.
//
// CINT notes: fixed-size arrays only, unique loop variable names.

void fttDataWindowRefit(int run,
                        const char *outdir = "refit",
                        double winKeep     = 0.95,
                        double bxPeriod    = 4.263,
                        int    winFloorMin = -1,
                        int    winFloorMax = 3,
                        const char *indir  = ".")
{
  gROOT->LoadMacro(Form("%s/fttDataWindow.C", gSystem->DirName(__FILE__)));

  TString fin  = Form("%s/fttDataWindow.%d.txt", indir, run);
  TString froot = Form("%s/root/fttDataWindow.%d.root", indir, run);
  if (gSystem->AccessPathName(fin.Data()))   { printf("fttDataWindowRefit: no %s\n", fin.Data());   return; }
  if (gSystem->AccessPathName(froot.Data())) { printf("fttDataWindowRefit: no %s\n", froot.Data()); return; }

  // ---- read the existing calibration: anchors, and everything we pass through ----
  int  uuid[400], mode[400], omin[400], omax[400], anch[400];
  char rest[400][256];
  int  nLine = 0;
  TString hdr = "";
  char line[512];
  FILE *I = fopen(fin.Data(), "r");
  while (fgets(line, sizeof(line), I)) {
    if (line[0] == '#') { hdr += line; continue; }
    int u, m, lo, hi, a, nc = 0;
    if (sscanf(line, "%d %d %d %d %d %n", &u, &m, &lo, &hi, &a, &nc) < 5) continue;
    if (nLine >= 400) break;
    uuid[nLine] = u; mode[nLine] = m; omin[nLine] = lo; omax[nLine] = hi; anch[nLine] = a;
    strncpy(rest[nLine], line + nc, 255); rest[nLine][255] = 0;
    // strip the trailing newline so the new line can be assembled cleanly
    for (int ic = 0; ic < 255; ic++) if (rest[nLine][ic] == '\n') { rest[nLine][ic] = 0; break; }
    nLine++;
  }
  fclose(I);
  printf("fttDataWindowRefit: run %d, %d lines from %s\n", run, nLine, fin.Data());

  TFile *F = TFile::Open(froot.Data());

  // ---- refit the window of every calibrated VMM ---------------------------------
  int nFit = 0, nFloor = 0, nSkip = 0, nSame = 0;
  int newLo[400], newHi[400], fitLo[400], fitHi[400];
  for (int jl = 0; jl < nLine; jl++) {
    newLo[jl] = omin[jl]; newHi[jl] = omax[jl]; fitLo[jl] = 0; fitHi[jl] = 0;
    if (uuid[jl] < 0 || uuid[jl] >= 384) continue;    // slot 384 is the run stamp
    if (anch[jl] < 0) { nSkip++; continue; }          // no usable anchor: keep the wide default
    TH1 *h = (TH1*) F->Get(Form("hDbcid_%03d", uuid[jl]));
    if (!h) { nSkip++; continue; }
    if (!fitFttWindow(h, anch[jl], bxPeriod, winKeep)) { nSkip++; continue; }
    fitLo[jl] = gFitLo; fitHi[jl] = gFitHi;
    int lo = gFitLo, hi = gFitHi;
    if (lo > winFloorMin || hi < winFloorMax) nFloor++;
    if (lo > winFloorMin) lo = winFloorMin;
    if (hi < winFloorMax) hi = winFloorMax;
    newLo[jl] = lo; newHi[jl] = hi; nFit++;
    if (lo == omin[jl] && hi == omax[jl]) nSame++;
  }

  // the run's median window, for the VMMs that are calibrated but could not be fitted
  int medLo = winFloorMin, medHi = winFloorMax;
  {
    int vlo[400], vhi[400], nv = 0;
    for (int jm = 0; jm < nLine; jm++)
      if (fitLo[jm] != 0 || fitHi[jm] != 0) { vlo[nv] = newLo[jm]; vhi[nv] = newHi[jm]; nv++; }
    if (nv > 0) { std::sort(vlo, vlo + nv); std::sort(vhi, vhi + nv); medLo = vlo[nv / 2]; medHi = vhi[nv / 2]; }
  }
  for (int jn = 0; jn < nLine; jn++)
    if (uuid[jn] >= 0 && uuid[jn] < 384 && anch[jn] >= 0 && fitLo[jn] == 0 && fitHi[jn] == 0) {
      newLo[jn] = medLo; newHi[jn] = medHi;
    }

  // ---- write ---------------------------------------------------------------------
  gSystem->mkdir(outdir, kTRUE);
  TString fout = Form("%s/fttDataWindow.%d.txt", outdir, run);
  FILE *O = fopen(fout.Data(), "w");
  fprintf(O, "%s", hdr.Data());
  fprintf(O, "# REFIT %s: window recomputed from root/fttDataWindow.%d.root, anchors unchanged.\n",
          TDatime().AsSQLString(), run);
  fprintf(O, "#       winKeep %.2f bxPeriod %.3f floor [%d, %d]; %d VMMs refitted, %d widened to the\n"
             "#       floor, %d unchanged, %d left alone (no anchor or no fit), run median [%d, %d].\n",
          winKeep, bxPeriod, winFloorMin, winFloorMax, nFit, nFloor, nSame, nSkip, medLo, medHi);
  for (int jw = 0; jw < nLine; jw++)
    fprintf(O, "%3d %d %4d %4d %5d  %s\n", uuid[jw], mode[jw], newLo[jw], newHi[jw], anch[jw], rest[jw]);
  fclose(O);
  F->Close();

  printf("fttDataWindowRefit: %d refitted (%d widened to the floor, %d unchanged), "
         "%d left alone, median [%d, %d]\n", nFit, nFloor, nSame, nSkip, medLo, medHi);
  printf("fttDataWindowRefit: wrote %s\n", fout.Data());
}
