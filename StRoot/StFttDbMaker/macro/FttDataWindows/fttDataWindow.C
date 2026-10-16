// fttDataWindow.C -- per-VMM FTT (sTGC) time calibration from MuDst, for one run.
//
// Output (in <outdir>, default ./fttDataWindow/):
//   fttDataWindow.<run>.txt        calibration, one line per VMM, first 5 columns in
//                                  Calibrations/ftt/fttDataWindowsB order
//                                  (uuid mode min max anchor), QA columns after that
//   fttDataWindow.<year>.QA.txt    one summary line APPENDED per run (year from run number)
//   png/fttDataWindow.<run>.png    QA plots
//   root/fttDataWindow.<run>.root  per-VMM dbcid histograms
// Push to DB with script/fttDataWindow_db.C.
//
// Usage:
//   root4star -l -b -q 'fttDataWindow.C("st_fwd_23081049_raw_5500015.MuDst.root")'
//   root4star -l -b -q 'fttDataWindow.C("run23063028.list")'   // text list, one MuDst per line
//
// WHY THIS EXISTS. StFttHitCalibMaker finds each VMM's dbcid anchor on the fly,
// separately in every job. On sparse data (the zero-field physics-stream runs have
// ~31x fewer FTT hits/event than the fwd stream) that never converges, so a
// calibrated time cut throws most hits away. Measuring anchors once per run from
// enough events, and storing them, removes the per-file warm-up.
//
// METHOD. Raw dbcid per VMM is histogrammed (the anchor does not depend on
// StFttHitCalibMaker or on any time cut). Periodically, and at the end:
//   anchor  = most populated dbcid (same definition as HitCalibHelper), found on a
//             4-bin rebinned histogram and refined to the raw bin inside that group
//   t       = dbcid - anchor, wrapped to [-2048, 2047] (12-bit circular counter)
//   bkg     = mean counts/bin in the sidebands 500 <= |t| <= 1500
//   signal  = counts in [winMin, winMax] - bkg * window width
// Because the far sidebands are nearly empty (the out-of-window hits form a broad
// plateau right next to the spike, not a flat background), "signal >= minPeak" alone
// is too weak: a window dropped anywhere on a busy plateau passes it. So a VMM is
// only accepted if its in-window FRACTION is >= fracTol x the run median. A VMM
// that fails is retried at the median anchor of the OTHER accepted VMMs on the same
// FOB (anchors are constant per FOB = 4 VMMs; measured on run 23081049, FOBs within
// one quadrant differ by ~30 ticks and occasionally by ~1500, so the quadrant is the
// wrong grouping). Accepted with status 4 if that passes, otherwise anchor -1.
// BUT a low in-window fraction alone does not mean a wrong anchor -- it also flags
// noisy channels whose anchor is right. So first: if the VMM's own anchor agrees with
// its FOB median to within fobTol ticks, keep it (status 5). Measured: on ZF run
// 23063028, 11 of 17 fraction-failures sat 0-3 ticks from all their FOB neighbours
// (one with 13773 hits); the other 6, and all 10 on field-on run 23081049, were
// 250-2000 ticks away. VMMs 200 and 369 were ~500 off in BOTH runs -- likely a
// chip/configuration problem, not statistics.
// SAFETY FLOOR. The fitted window is only ever allowed to GROW. Whatever the fit
// returns, the stored window is widened to contain [winFloorMin, winFloorMax],
// default [-1, 3] -- the one-crossing window of the summed-over-all-VMM study
// (ftt_sim_hit_maker/timecut/index.html section 5). Rationale: the measurement is
// per VMM and per run, from a few hundred events, so a VMM whose crossing is a
// tick thin for want of statistics must not end up with a cut TIGHTER than the
// number we already trust. A window that is one tick too wide costs a little
// pileup; one that is a tick too narrow throws real hits away, in production,
// silently. The fit therefore only ever tells us where to open the window FURTHER.
// Measured (runs 24192016/17): the fit lands on [-1,2] for about half the VMMs, so
// the floor's upper edge is what is usually stored and the fit mostly decides the
// tail -- which is the intended asymmetry.
//
// Stops when goodFrac of ACTIVE VMMs (>= nActive hits) have signal >= minPeak, and
// at least minEvents were read; or at maxEvents; or at end of input. So dense runs
// stop after a few hundred events and sparse ones run as long as they need.
//
// CONVENTIONS -- following what production actually stores (checked against the
// existing DB entry, 2026-09-14). The idl uuid comment is fixed on Daniel's dev
// (PR #20); the idl mode comment ("0 = timebin, 1 = bcid") is still inverted:
//   uuid = StFttDb::vmmId(), 0-based 0..383
//        = vmm + 4*(feb + 6*((rdo-1) + 4*(sector-1)))
//   mode 0 = cut on CALIBRATED time (StFttClusterMaker kHitCalibratedTime = 0)
// anchor = -1 marks a VMM with no usable calibration (status 1,2,3); consumers must
// fall back to the on-the-fly anchor for it.
//
// CINT notes: fixed-size arrays and histograms only, unique loop variable names.

// ---- constants inlined from StFttDb.h (class statics are not reliably visible) --
const int kNVMM   = 384;
#include <algorithm>

const int kNDBCID = 4096;
// Half-range scanned when fitting a VMM's in-time crossing. Must cover the
// neighbour crossings used for the background: 6 * 4.263 + 1 = 27 ticks.
const int kWinScan = 40;
const int kNOrient = 5;   // kFttHorizontal=0, kFttVertical=1, kFttDiagonalH=2, kFttDiagonalV=3, unknown=4

// ---- the window fit, used by fttDataWindow.C and by fttDataWindowRefit.C ------
// Out-of-time hits arrive in crossings every bxPeriod ticks, so what sits under the
// in-time peak is pileup at the SAME phase: bg(t) = mean of counts at t + k*bxPeriod,
// k = +-1..+-6, interpolated between ticks, skipping neighbours that fall inside the
// in-time crossing itself. A linear slope of the plateau cancels. (sugMin/sugMax in
// the QA columns subtract the FAR SIDEBAND instead, which leaves the whole plateau
// in -- which is why they land around [-100,+225]: that is the plateau, not the
// crossing.) The window then grows outward from the peak, always taking the richer
// side, until it holds winKeep of the crossing's signal -- so the parameter is the
// "signal kept" of the timecut study rather than a threshold on an arbitrary shape.
// One crossing keeps ~83% of the signal at ~8x the purity of [-65,100];
// see ftt_sim_hit_maker/timecut/index.html section 5.
//
// h is a VMM's raw dbcid histogram and aQ its anchor BIN. Returns true and sets
// gFitLo/gFitHi (globals rather than reference arguments: CINT resolves signatures
// at load time and is happier passing nothing by reference across macros).
// gFitSig/gFitBkg hold the per-tick signal and background of the last call, for
// plotting; index t runs -kWinScan..+kWinScan at offset kWinScan.
int gFitLo = 0, gFitHi = 0;
double gFitSig[2 * kWinScan + 1], gFitBkg[2 * kWinScan + 1];
int    gFitPeak = 0;

// inTimeHalf <= 0 means "half a bunch crossing", which is the right answer unless
// you are reproducing old behaviour. Not a default argument: CINT does not honour
// one when the function is called from another loaded macro.
bool fitFttWindow(TH1 *h, int aQ, double bxPeriod, double winKeep, double inTimeHalf)
{
  // inTimeHalf is the half-width of the in-time crossing, in ticks. It is used in two
  // places that must agree: neighbours this close to t are left out of the pileup
  // background, and the signal integral that the window is a fraction of runs over
  // the same range. It must not exceed HALF A BUNCH CROSSING, which is the default --
  // beyond that the range reaches into the neighbouring crossings, and on a VMM with
  // a weak peak the "95% of the crossing" criterion is then met by admitting their
  // shoulders rather than real in-time hits. Measured on run 23081049: the in-time
  // crossing is the three ticks t = -1, 0, +1, ticks +-2 and +-3 are an empty gap,
  // and the next crossings stand up clearly at t = +-4. At 3.0 ticks the fit reached
  // into that gap and out to the neighbours -- VMM 193 came out [-3,3], 7 ticks wide
  // for one 4.263-tick crossing, gaining 7 counts of "signal" at t=+3 while admitting
  // 16 of background. Half a period cannot reach a neighbour, by construction.
  if (inTimeHalf <= 0) inTimeHalf = 0.5 * bxPeriod;
  double cntW[2 * kWinScan + 1];
  for (int tw = -kWinScan; tw <= kWinScan; tw++)
    cntW[tw + kWinScan] = h->GetBinContent(((aQ + tw) % kNDBCID + kNDBCID) % kNDBCID + 1);
  for (int tw2 = -kWinScan; tw2 <= kWinScan; tw2++) {
    double bsum = 0; int nb = 0;
    for (int kk = -6; kk <= 6; kk++) {
      if (kk == 0) continue;
      double x = tw2 + kk * bxPeriod;
      if (fabs(x) <= inTimeHalf) continue;
      int x0 = (int)floor(x);
      double fx = x - x0;
      int i0 = ((aQ + x0)     % kNDBCID + kNDBCID) % kNDBCID;
      int i1 = ((aQ + x0 + 1) % kNDBCID + kNDBCID) % kNDBCID;
      bsum += (1 - fx) * h->GetBinContent(i0 + 1) + fx * h->GetBinContent(i1 + 1);
      nb++;
    }
    gFitBkg[tw2 + kWinScan] = (nb > 0 ? bsum / nb : 0);
    gFitSig[tw2 + kWinScan] = cntW[tw2 + kWinScan] - gFitBkg[tw2 + kWinScan];
  }
  int tPk = 0; double sPk = -1e30;
  for (int tw3 = -3; tw3 <= 5; tw3++) if (gFitSig[tw3 + kWinScan] > sPk) { sPk = gFitSig[tw3 + kWinScan]; tPk = tw3; }
  gFitPeak = tPk;
  int xlo = tPk - (int)inTimeHalf, xhi = tPk + (int)inTimeHalf;
  double totX = 0;
  for (int tw4 = xlo; tw4 <= xhi; tw4++)
    if (gFitSig[tw4 + kWinScan] > 0) totX += gFitSig[tw4 + kWinScan];
  if (sPk <= 0 || totX <= 0) return false;
  int lo = tPk, hi = tPk;
  double got = sPk;
  while (got < winKeep * totX && (lo - 1 >= xlo || hi + 1 <= xhi)) {
    double sl = (lo - 1 >= xlo) ? gFitSig[lo - 1 + kWinScan] : -1e30;
    double sh = (hi + 1 <= xhi) ? gFitSig[hi + 1 + kWinScan] : -1e30;
    if (sl >= sh) { lo--; if (sl > 0) got += sl; }
    else          { hi++; if (sh > 0) got += sh; }
  }
  gFitLo = lo; gFitHi = hi;
  return true;
}

int wrapT(int t)
{
  while (t >= 2048) t -= 4096;
  while (t < -2048) t += 4096;
  return t;
}

// Integral of a dbcid histogram over [d1, d2] in dbcid units, circular.
double circInt(TH1 *h, int d1, int d2)
{
  if (d2 < d1) return 0;
  if (d2 - d1 + 1 >= kNDBCID) return h->Integral(1, kNDBCID);
  int a = d1, b = d2;
  while (a < 0) { a += kNDBCID; b += kNDBCID; }
  while (a >= kNDBCID) { a -= kNDBCID; b -= kNDBCID; }
  if (b < kNDBCID) return h->Integral(a + 1, b + 1);
  return h->Integral(a + 1, kNDBCID) + h->Integral(1, b - kNDBCID + 1);
}

// Anchor: coarse argmax on a 4-bin rebin (robust at low statistics), refined to
// the raw argmax inside that 4-bin group.
int findAnchor(TH1 *h)
{
  // Fullest window of 4 consecutive dbcid, SLIDING and circular, then the fullest
  // raw bin inside it.
  //
  // This used to be Rebin(4), i.e. windows aligned to multiples of 4, and that loses
  // a weak peak that straddles a group boundary. Measured on VMM 54 of run 23081049:
  // its real peak is 60+60 counts at dbcid 1070,1071 over a level of 25, which is
  // 9.9 sigma and certainly not a fluctuation -- but 1071 is the last bin of its
  // aligned group, so that group summed 170 while two neighbouring groups, holding
  // only the broad tail, summed 173 and 171. The anchor came out 1075, five ticks
  // late, and the window fit then chased the peak to [-5,0]. The sliding window has
  // no boundary to straddle: it finds 1069..1072 = 192 and refines to 1070, which is
  // both the raw argmax and the anchor of the other three VMMs on the FOB.
  //
  // Cost is 4 x 4096 adds per call, which is nothing next to reading the events.
  int nb = h->GetNbinsX();
  double bestS = -1;
  int at = 0;
  for (int db = 0; db < nb; db++) {
    double sum4 = 0;
    for (int dk = 0; dk < 4; dk++) sum4 += h->GetBinContent((db + dk) % nb + 1);
    if (sum4 > bestS) { bestS = sum4; at = db; }
  }
  int best = at;
  double bestC = -1;
  for (int dj = 0; dj < 4; dj++) {
    double c = h->GetBinContent((at + dj) % nb + 1);
    if (c > bestC) { bestC = c; best = (at + dj) % nb; }
  }
  return best;
}

void loadLibsDW()
{
  gSystem->Load("libStarClassLibrary.so");
  gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StarMagField");
  gSystem->Load("StMagF");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbBroker");
  gSystem->Load("StDbUtilities");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StEvent");
  gSystem->Load("StEventMaker");
  gSystem->Load("StMuDSTMaker");
  gSystem->Load("StFttDbMaker");
}

void fttDataWindow(const char *input,
                   int    maxEvents  = 50000,
                   int    minPeak    = 200,     // bkg-subtracted in-window hits per VMM
                   double goodFrac   = 0.95,
                   int    winMin     = -65,     // wide fallback, kept by VMMs with no
                   int    winMax     = 100,     // usable anchor; also the DB entry today
                   const char *outdir = "fttDataWindow",
                   int    checkEvery = 200,
                   int    minEvents  = 200,
                   int    nActive    = 50,
                   double fracTol    = 0.7,     // accept if inWin fraction >= fracTol * run median
                   int    fobTol     = 20,      // ...or if own anchor is within this many ticks of FOB median
                   double bxPeriod   = 4.263,   // ticks per bunch crossing (40 MHz / 9.383 MHz)
                   double winKeep    = 0.95,    // window holds this fraction of the in-time crossing
                   int    winFloorMin = -1,     // the stored window always CONTAINS
                   int    winFloorMax = 3)      // [winFloorMin, winFloorMax] -- see below
{
  loadLibsDW();
  gSystem->mkdir(Form("%s/png", outdir), kTRUE);
  gSystem->mkdir(Form("%s/root", outdir), kTRUE);

  // ---- input: one MuDst, or a text list of them -------------------------------
  TString in(input);
  TString first = "";
  int nFiles = 0;
  const int kMaxFiles = 2000;
  TString files[2000];
  if (in.EndsWith(".root")) {
    files[0] = in; nFiles = 1;
  } else {
    ifstream lst(input);
    TString line;
    while (line.ReadLine(lst) && nFiles < kMaxFiles) {
      line = line.Strip(TString::kBoth);
      if (line.Length() == 0 || line.BeginsWith("#")) continue;
      files[nFiles++] = line;
    }
  }
  if (nFiles == 0) { printf("fttDataWindow: no input files in %s\n", input); return; }

  StChain *chain = new StChain("StChain");
  StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", files[0].Data(), "MuDst.root", 1);
  TChain &muChain = *muDstMaker->chain();
  for (int iF = 1; iF < nFiles; iF++) muChain.Add(files[iF].Data());
  // Only the event header and FTT raw hits are needed -- a large speed-up on
  // physics-stream files.
  muDstMaker->SetStatus("*", 0);
  muDstMaker->SetStatus("MuEvent", 1);
  muDstMaker->SetStatus("FttRawHit", 1);

  St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
  StMuDst2StEventMaker *mu2ev = new StMuDst2StEventMaker();
  StFttDbMaker *fttDbMk = new StFttDbMaker();

  if (chain->Init()) { printf("fttDataWindow: chain Init failed\n"); return; }
  if (muChain.GetEntries() <= 0) {
    printf("fttDataWindow: ERROR no events readable from %s (missing file or unavailable disk?)\n", input);
    return;
  }

  muChain.GetEntry(0);
  int run = muDstMaker->muDst()->event()->runNumber();
  // UTC date/time of the first event, for the DB query. TTimeStamp rather than
  // gmtime(): CINT does not reliably know struct tm without extra includes.
  TTimeStamp firstTs((time_t)muDstMaker->muDst()->event()->eventInfo()->time(), 0);
  int date  = (int)firstTs.GetDate(kTRUE);
  int itime = (int)firstTs.GetTime(kTRUE);
  dbMk->SetDateTime(date, itime);
  chain->InitRun(run);
  printf("fttDataWindow: run %d, first event %d %06d, %d file(s), %lld entries\n",
         run, date, itime, nFiles, (Long64_t)muChain.GetEntries());

  // ---- accumulators -------------------------------------------------------------
  TH1I *hV[384];
  for (int iv = 0; iv < kNVMM; iv++) {
    hV[iv] = new TH1I(Form("hDbcid_%03d", iv), Form("VMM %d;dbcid;hits", iv), kNDBCID, 0, kNDBCID);
    hV[iv]->SetDirectory(0);
  }
  double orientCount[384][5];
  for (int io = 0; io < kNVMM; io++) for (int jo = 0; jo < kNOrient; jo++) orientCount[io][jo] = 0;

  int anchor[384], status[384];
  double nHitsV[384], sigV[384], bkgV[384], inWinV[384];

  Long64_t nEntries = muChain.GetEntries();
  if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;
  double nHitsTot = 0, nBadVmm = 0, nBadDbcid = 0, nOtherRun = 0;
  int nRead = 0;
  TString stopReason = "endOfInput";

  for (Long64_t iev = 0; iev < muChain.GetEntries(); iev++) {
    if (maxEvents > 0 && iev >= maxEvents) { stopReason = "maxEvents"; break; }
    chain->Clear();
    int rc = chain->Make();
    if (rc == kStEOF) break;
    if (rc != kStOK) continue;
    nRead++;

    StEvent *ev = (StEvent *)chain->GetDataSet("StEvent");
    if (!ev) continue;
    if (ev->runId() != run) { nOtherRun++; continue; }
    StFttCollection *fc = ev->fttCollection();
    if (!fc) continue;
    StFttDb *db = (StFttDb *)chain->GetDataSet("fttDb");
    if (!db) continue;

    int nh = fc->rawHits().size();
    for (int ih = 0; ih < nh; ih++) {
      StFttRawHit *hit = fc->rawHits()[ih];
      int v = (int)db->vmmId(hit);
      if (v < 0 || v >= kNVMM) { nBadVmm++; continue; }
      int d = hit->dbcid();
      if (d < 0 || d >= kNDBCID) { nBadDbcid++; continue; }
      hV[v]->Fill(d);
      nHitsTot++;
      db->hardwareMap(hit);
      int o = (int)db->orientation(hit);
      if (o < 0 || o > 4) o = 4;
      orientCount[v][o]++;
    }

    // ---- convergence check ----------------------------------------------------
    if ((nRead % checkEvery) == 0 && nRead >= minEvents) {
      int nAct = 0, nGood = 0;
      for (int ic = 0; ic < kNVMM; ic++) {
        double nC = hV[ic]->GetEntries();
        if (nC < nActive) continue;
        nAct++;
        int aC = hV[ic]->GetMaximumBin() - 1;      // cheap raw argmax is enough here
        double bC = (circInt(hV[ic], aC + 500, aC + 1500) + circInt(hV[ic], aC - 1500, aC - 500)) / 2002.0;
        double inC = circInt(hV[ic], aC + winMin, aC + winMax);
        double sC = inC - bC * (winMax - winMin + 1);
        // same fraction requirement as the final pass, with a nominal 0.5 median,
        // so a window parked on the plateau cannot declare convergence early
        if (sC >= minPeak && inC / nC >= fracTol * 0.5) nGood++;
      }
      double frC = nAct > 0 ? (double)nGood / nAct : 0.0;
      printf("fttDataWindow: %6d events  active VMMs %3d  good %3d  (%.3f)\n", nRead, nAct, nGood, frC);
      if (nAct > 0 && frC >= goodFrac) { stopReason = "converged"; break; }
    }
  }

  // ---- final per-VMM calibration ------------------------------------------------
  const int kTZoom = 300;
  TH1D *hAligned = new TH1D("hAligned", "all good VMMs;t = dbcid - anchor;hits", 1001, -500.5, 500.5);
  TH2D *hTvsV = new TH2D("hTvsV", "t vs VMM (good);vmmId;t", kNVMM, 0, kNVMM, 2 * kTZoom + 1, -kTZoom - 0.5, kTZoom + 0.5);
  TH1D *hAnchor = new TH1D("hAnchor", "anchor (good);vmmId;anchor dbcid", kNVMM, 0, kNVMM);
  TH1D *hAnchorBad = new TH1D("hAnchorBad", "anchor (not ok)", kNVMM, 0, kNVMM);
  TH1D *hN = new TH1D("hN", "hits per VMM;vmmId;hits", kNVMM, 0, kNVMM);
  TH1D *hFrac = new TH1D("hFrac", "bkg-subtracted signal fraction in window;vmmId;fraction", kNVMM, 0, kNVMM);

  int nOK = 0, nLow = 0, nNoPeak = 0, nDead = 0;
  double orInWin[5], orTot[5];
  for (int ko = 0; ko < kNOrient; ko++) { orInWin[ko] = 0; orTot[ko] = 0; }
  int dominant[384];
  int sugMin[384], sugMax[384];
  // Per-VMM window actually written to the DB (columns 3-4). Starts as the wide
  // default and is replaced below by the in-time bunch crossing where one is found.
  int winMinV[384], winMaxV[384];
  // What the fit itself returned, before the safety floor -- QA columns 14-15, so
  // that a run's real per-VMM measurement stays visible in the file we push from.
  int fitLoV[384], fitHiV[384];
  int nWinOwn = 0, nWinFloor = 0;

  for (int jv = 0; jv < kNVMM; jv++) {
    nHitsV[jv] = hV[jv]->GetEntries();
    dominant[jv] = 4;
    double bestO = 0;
    for (int lo = 0; lo < kNOrient; lo++) if (orientCount[jv][lo] > bestO) { bestO = orientCount[jv][lo]; dominant[jv] = lo; }
    sugMin[jv] = winMin; sugMax[jv] = winMax;
    winMinV[jv] = winMin; winMaxV[jv] = winMax;
    fitLoV[jv] = 0; fitHiV[jv] = 0;

    if (nHitsV[jv] <= 0) { status[jv] = 3; anchor[jv] = -1; sigV[jv] = 0; bkgV[jv] = 0; inWinV[jv] = 0; nDead++; continue; }

    int a = findAnchor(hV[jv]);
    double b = (circInt(hV[jv], a + 500, a + 1500) + circInt(hV[jv], a - 1500, a - 500)) / 2002.0;
    double inW = circInt(hV[jv], a + winMin, a + winMax);
    double s = inW - b * (winMax - winMin + 1);
    bkgV[jv] = b; inWinV[jv] = inW; sigV[jv] = s;

    if (nHitsV[jv] < nActive)      status[jv] = 1;
    else if (s < minPeak)          status[jv] = 2;
    else                           status[jv] = 0;
    anchor[jv] = a;               // provisional; finalised in the consistency pass below

  }

  // ---- consistency pass: in-window fraction vs run median, quadrant fallback ------
  int anchorRaw[384];
  for (int rv = 0; rv < kNVMM; rv++) anchorRaw[rv] = anchor[rv];
  double fracOwn[384];
  double fracList[384];
  int nFracList = 0;
  for (int cv = 0; cv < kNVMM; cv++) {
    fracOwn[cv] = (nHitsV[cv] > 0) ? inWinV[cv] / nHitsV[cv] : 0.0;
    if (status[cv] == 0) fracList[nFracList++] = fracOwn[cv];
  }
  double fracMed = 0;
  if (nFracList > 0) {
    int idxSort[384];
    TMath::Sort(nFracList, fracList, idxSort, kFALSE);
    fracMed = fracList[idxSort[nFracList / 2]];
  }

  // FOB median anchor (circular), from accepted VMMs on the same FOB (vmmId / 4)
  int fobMed[96];
  for (int fm = 0; fm < 96; fm++) {
    int ref = -1, nrel = 0;
    double rel[4];
    for (int fq = fm * 4; fq < fm * 4 + 4; fq++) {
      if (status[fq] != 0 || fracOwn[fq] < fracTol * fracMed) continue;
      if (ref < 0) ref = anchor[fq];
      rel[nrel++] = wrapT(anchor[fq] - ref);
    }
    fobMed[fm] = -1;
    if (nrel >= 2) {                         // need at least two agreeing neighbours
      int idxF[4];
      TMath::Sort(nrel, rel, idxF, kFALSE);
      fobMed[fm] = ((ref + (int)rel[idxF[nrel / 2]]) % kNDBCID + kNDBCID) % kNDBCID;
    }
  }

  nOK = 0; nLow = 0; nNoPeak = 0;
  int nFromFob = 0, nRejected = 0, nNoisy = 0;
  for (int fv = 0; fv < kNVMM; fv++) {
    if (status[fv] == 3) continue;
    if (status[fv] == 1) { anchor[fv] = -1; nLow++; continue; }
    if (fracOwn[fv] >= fracTol * fracMed && status[fv] == 0) { nOK++; continue; }
    int ff = fv / 4;
    // low fraction but anchor agrees with the FOB: correct anchor, noisy channel
    if (fobMed[ff] >= 0 && abs(wrapT(anchor[fv] - fobMed[ff])) <= fobTol) {
      status[fv] = 5; nNoisy++;
      continue;
    }
    // own anchor unusable -- try the median of the other accepted VMMs on this FOB
    if (fobMed[ff] >= 0) {
      int am = fobMed[ff];
      double inWm = circInt(hV[fv], am + winMin, am + winMax);
      double bm = (circInt(hV[fv], am + 500, am + 1500) + circInt(hV[fv], am - 1500, am - 500)) / 2002.0;
      double sm = inWm - bm * (winMax - winMin + 1);
      if (inWm / nHitsV[fv] >= fracTol * fracMed && sm >= minPeak) {
        anchor[fv] = am; status[fv] = 4; inWinV[fv] = inWm; bkgV[fv] = bm; sigV[fv] = sm;
        nFromFob++;
        continue;
      }
    }
    anchor[fv] = -1; status[fv] = 2; nNoPeak++; nRejected++;
  }
  printf("fttDataWindow: median in-window fraction %.3f; ok %d, noisy-but-FOB-consistent %d, anchor from FOB median %d, rejected %d\n",
         fracMed, nOK, nNoisy, nFromFob, nRejected);

  // ---- QA histograms and suggested windows, for VMMs with a final anchor ----------
  for (int pv = 0; pv < kNVMM; pv++) {
    if (status[pv] == 3) continue;
    double sQ = sigV[pv], bQ = bkgV[pv], inWQ = inWinV[pv];
    int aQ = anchor[pv];
    orInWin[dominant[pv]] += inWQ;
    orTot[dominant[pv]]   += nHitsV[pv];
    hN->SetBinContent(pv + 1, nHitsV[pv]);
    hFrac->SetBinContent(pv + 1, nHitsV[pv] > 0 ? sQ / nHitsV[pv] : 0);

    if (aQ < 0) { hAnchorBad->SetBinContent(pv + 1, anchorRaw[pv]); continue; }
    hAnchor->SetBinContent(pv + 1, aQ);

    // aligned distributions, and a data-derived window suggestion (QA only):
    // central 98% of the background-subtracted signal within [-300, +400]
    double cum[701];
    double tot = 0;
    for (int it = -300; it <= 400; it++) {
      double c = hV[pv]->GetBinContent(((aQ + it) % kNDBCID + kNDBCID) % kNDBCID + 1);
      double sc = c - bQ; if (sc < 0) sc = 0;
      tot += sc; cum[it + 300] = tot;
    }
    for (int it2 = -300; it2 <= 400 && tot > 0; it2++) if (cum[it2 + 300] >= 0.01 * tot) { sugMin[pv] = it2; break; }
    for (int it3 = 400; it3 >= -300 && tot > 0; it3--) if (cum[it3 + 300] <= 0.99 * tot) { sugMax[pv] = it3 + 1; break; }

      // ---- the window actually stored: the in-time bunch crossing -----------------
      if (fitFttWindow(hV[pv], aQ, bxPeriod, winKeep, 0)) {
        fitLoV[pv] = gFitLo; fitHiV[pv] = gFitHi;
        int lo = gFitLo, hi = gFitHi;
        // Grow-only: never store a window tighter than the floor.
        if (lo > winFloorMin || hi < winFloorMax) nWinFloor++;
        if (lo > winFloorMin) lo = winFloorMin;
        if (hi < winFloorMax) hi = winFloorMax;
        winMinV[pv] = lo; winMaxV[pv] = hi; nWinOwn++;
      }

    for (int dd = 0; dd < kNDBCID; dd++) {
      double c2 = hV[pv]->GetBinContent(dd + 1);
      if (c2 <= 0) continue;
      int t = wrapT(dd - aQ);
      if (t >= -500 && t <= 500) hAligned->Fill(t, c2);
      if (t >= -kTZoom && t <= kTZoom) hTvsV->Fill(pv, t, c2);
    }
  }


  // ---- write calibration --------------------------------------------------------
  TString fout = Form("%s/fttDataWindow.%d.txt", outdir, run);
  FILE *F = fopen(fout.Data(), "w");
  fprintf(F, "# fttDataWindowsB calibration, run %d\n", run);
  fprintf(F, "# run %d firstEventDate %d firstEventTime %06d\n", run, date, itime);
  fprintf(F, "# input %s nFiles %d eventsRead %d rawHits %.0f stop %s\n", input, nFiles, nRead, nHitsTot, stopReason.Data());
  // A calibrated VMM whose crossing could not be fitted takes the run's median
  // window rather than the wide default -- its anchor is good, only its shape was
  // too thin to measure. VMMs with anchor -1 keep the wide default: they fall back
  // to the consumer's on-the-fly anchor, which is not validated to the +-1 tick a
  // one-crossing window needs.
  int medLo = winMin, medHi = winMax;
  int nNoAnchor = 0;
  for (int nv0 = 0; nv0 < kNVMM; nv0++) if (anchor[nv0] < 0) nNoAnchor++;
  {
    int vlo[384], vhi[384], nv = 0;
    for (int mv = 0; mv < kNVMM; mv++)
      if (anchor[mv] >= 0 && winMinV[mv] != winMin) { vlo[nv] = winMinV[mv]; vhi[nv] = winMaxV[mv]; nv++; }
    if (nv > 0) {
      std::sort(vlo, vlo + nv); std::sort(vhi, vhi + nv);
      medLo = vlo[nv / 2]; medHi = vhi[nv / 2];
      for (int mv2 = 0; mv2 < kNVMM; mv2++)
        if (anchor[mv2] >= 0 && winMinV[mv2] == winMin) { winMinV[mv2] = medLo; winMaxV[mv2] = medHi; }
    }
  }
  printf("fttDataWindow: per-VMM window fitted for %d VMMs, run median [%d, %d]; "
         "%d widened to the floor [%d, %d]; %d calibrated VMMs took the median, "
         "%d uncalibrated kept [%d, %d]\n",
         nWinOwn, medLo, medHi, nWinFloor, winFloorMin, winFloorMax,
         nOK + nNoisy + nFromFob - nWinOwn, nNoAnchor, winMin, winMax);

  fprintf(F, "# window: per VMM, the in-time bunch crossing (bxPeriod %.3f, holding %.2f of its signal),\n"
             "#         then widened if needed to contain the floor [%d, %d] -- grow-only, a fit\n"
             "#         is never allowed to make the cut tighter than that. %d VMMs were widened.\n"
             "#         run median [%d, %d]; uncalibrated VMMs keep the wide [%d, %d].\n"
             "#         minPeak %d, goodFrac %.3f\n",
          bxPeriod, winKeep, winFloorMin, winFloorMax, nWinFloor,
          medLo, medHi, winMin, winMax, minPeak, goodFrac);
  fprintf(F, "# conventions: uuid = StFttDb::vmmId() 0..383; mode 0 = calibrated time (the idl mode comment says the opposite)\n");
  fprintf(F, "# cols 1-5 -> DB: uuid mode min max anchor   (anchor -1 = no usable calibration, use on-the-fly anchor)\n");
  fprintf(F, "# cols 6+  QA only, ignored by fttDataWindow_db.C: status nHits sig bkgPerBin sigFrac orient sugMin sugMax fitLo fitHi\n");
  fprintf(F, "#           fitLo/fitHi = what the crossing fit returned before the floor was applied (0 0 = not fitted)\n");
  fprintf(F, "# status 0=ok 1=lowstat(<%d hits) 2=rejected(sig<%d or inWin fraction < %.2f x median %.3f) 3=nohits 4=anchor taken from FOB median 5=own anchor kept, agrees with FOB within %d ticks but low in-window fraction (noisy); orient 0=H 1=V 2=dH 3=dV 4=unknown\n", nActive, minPeak, fracTol, fracMed, fobTol);
  for (int wv = 0; wv < kNVMM; wv++) {
    fprintf(F, "%3d %d %4d %4d %5d  %d %9.0f %9.1f %8.3f %6.3f %d %4d %4d %4d %4d\n",
            wv, 0, winMinV[wv], winMaxV[wv], anchor[wv],
            status[wv], nHitsV[wv], sigV[wv], bkgV[wv],
            nHitsV[wv] > 0 ? sigV[wv] / nHitsV[wv] : 0.0, dominant[wv], sugMin[wv], sugMax[wv],
            fitLoV[wv], fitHiV[wv]);
  }
  fclose(F);
  printf("fttDataWindow: wrote %s\n", fout.Data());

  // ---- one QA line per run ------------------------------------------------------
  int runYear = 2000 + run / 1000000 - 1;       // 23xxxxxx -> 2022
  TString fqa = Form("%s/fttDataWindow.%d.QA.txt", outdir, runYear);
  bool newQa = gSystem->AccessPathName(fqa.Data());     // true if it does NOT exist
  FILE *Q = fopen(fqa.Data(), "a");
  // Fixed-width columns so that "cat *.qa" is readable as a table. The header is written
  // as its own fragment, 0_header.qa, which sorts before any run number so it lands first.
  const char *qaHdrFmt = "#%7s %8s %6s %2s %7s %10s %4s %5s %4s %4s %4s %4s %5s %5s %4s %4s %5s %5s %5s %5s %5s %s\n";
  const char *qaRowFmt = "%8d %8d %06d %2d %7d %10.0f %4d %5d %4d %4d %4d %4d %5d %5d %4d %4d %5.3f %5.3f %5.3f %5.3f %5.3f %s\n";
  if (newQa) fprintf(Q, qaHdrFmt, "run", "date", "time", "nF", "events", "rawHits",
                     "nOK", "noisy", "fFob", "lowS", "rej", "dead", "wLo", "wHi", "nFit", "nFlr",
                     "all", "H", "V", "dH", "dV", "stop");
  double allIn = 0, allTot = 0;
  for (int qo = 0; qo < 4; qo++) { allIn += orInWin[qo]; allTot += orTot[qo]; }
  fprintf(Q, qaRowFmt,
          run, date, itime, nFiles, nRead, nHitsTot, nOK, nNoisy, nFromFob, nLow, nNoPeak, nDead,
          medLo, medHi, nWinOwn, nWinFloor,
          allTot > 0 ? allIn / allTot : 0.0,
          orTot[0] > 0 ? orInWin[0] / orTot[0] : 0.0, orTot[1] > 0 ? orInWin[1] / orTot[1] : 0.0,
          orTot[2] > 0 ? orInWin[2] / orTot[2] : 0.0, orTot[3] > 0 ? orInWin[3] / orTot[3] : 0.0,
          stopReason.Data());
  fclose(Q);

  // Standalone header fragment, same widths as the rows: runfttcal copies this into
  // <outdir>/qa/0_header.qa so "cat *.qa" prints a header line first.
  {
    TString fhd = Form("%s/fttDataWindow.%d.header.qa.txt", outdir, runYear);
    FILE *H = fopen(fhd.Data(), "w");
    if (H) {
      fprintf(H, qaHdrFmt, "run", "date", "time", "nF", "events", "rawHits",
              "nOK", "noisy", "fFob", "lowS", "rej", "dead", "wLo", "wHi", "nFit", "nFlr",
              "all", "H", "V", "dH", "dV", "stop");
      fclose(H);
    }
  }

  // ---- QA plots -----------------------------------------------------------------
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1", "fttDataWindow", 1500, 950);
  c1->Divide(3, 2);

  c1->cd(1); gPad->SetLogy();
  hAligned->SetTitle(Form("run %d: all ok VMMs, t = dbcid - anchor", run));
  hAligned->Draw("hist");
  double ymax = hAligned->GetMaximum() * 2, ymin = 0.5;
  // red = the wide default kept by uncalibrated VMMs; green = the run's median stored
  // window, which is what a calibrated VMM actually gets.
  TLine *l1 = new TLine(winMin, ymin, winMin, ymax); l1->SetLineColor(kRed); l1->Draw();
  TLine *l2 = new TLine(winMax, ymin, winMax, ymax); l2->SetLineColor(kRed); l2->Draw();
  TLine *l3 = new TLine(medLo, ymin, medLo, ymax); l3->SetLineColor(kGreen + 2); l3->Draw();
  TLine *l4 = new TLine(medHi, ymin, medHi, ymax); l4->SetLineColor(kGreen + 2); l4->Draw();

  c1->cd(2); gPad->SetLogz();
  hTvsV->Draw("colz");

  c1->cd(3);
  hAnchor->SetMarkerStyle(20); hAnchor->SetMarkerSize(0.5); hAnchor->SetMinimum(0); hAnchor->SetMaximum(4096);
  hAnchor->Draw("p");
  hAnchorBad->SetMarkerStyle(24); hAnchorBad->SetMarkerSize(0.5); hAnchorBad->SetMarkerColor(kRed);
  hAnchorBad->Draw("p same");

  c1->cd(4); gPad->SetLogy();
  hN->Draw("hist");

  c1->cd(5);
  hFrac->SetMinimum(-0.1); hFrac->SetMaximum(1.1);
  hFrac->Draw("hist");

  c1->cd(6);
  TPaveText *pt = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
  pt->SetTextAlign(12); pt->SetFillColor(0);
  pt->AddText(Form("run %d   first event %d %06d", run, date, itime));
  pt->AddText(Form("files %d   events read %d   stop: %s", nFiles, nRead, stopReason.Data()));
  pt->AddText(Form("raw hits %.0f", nHitsTot));
  pt->AddText(Form("VMM ok %d  noisy(FOB-ok) %d  fromFob %d", nOK, nNoisy, nFromFob));
  pt->AddText(Form("lowstat %d  rejected %d  nohits %d", nLow, nNoPeak, nDead));
  pt->AddText(Form("median in-window fraction %.3f (accept >= %.2f x)", fracMed, fracTol));
  pt->AddText(Form("window: median [%d, %d] stored, floor [%d, %d], %d fitted %d widened",
                   medLo, medHi, winFloorMin, winFloorMax, nWinOwn, nWinFloor));
  pt->AddText(Form("uncalibrated keep [%d, %d]   minPeak %d", winMin, winMax, minPeak));
  pt->AddText(Form("in-window (raw) all %.3f", allTot > 0 ? allIn / allTot : 0.0));
  pt->AddText(Form("  H %.3f  V %.3f  dH %.3f  dV %.3f",
                   orTot[0] > 0 ? orInWin[0] / orTot[0] : 0.0, orTot[1] > 0 ? orInWin[1] / orTot[1] : 0.0,
                   orTot[2] > 0 ? orInWin[2] / orTot[2] : 0.0, orTot[3] > 0 ? orInWin[3] / orTot[3] : 0.0));
  if (nBadVmm > 0 || nBadDbcid > 0 || nOtherRun > 0)
    pt->AddText(Form("skipped: badVmm %.0f badDbcid %.0f otherRunEvents %.0f", nBadVmm, nBadDbcid, nOtherRun));
  pt->Draw();

  TString fpng = Form("%s/png/fttDataWindow.%d.png", outdir, run);
  c1->SaveAs(fpng.Data());

  printf("fttDataWindow: run %d events %d stop %s  VMM ok %d noisyFobOK %d fromFob %d lowstat %d rejected %d nohits %d\n",
         run, nRead, stopReason.Data(), nOK, nNoisy, nFromFob, nLow, nNoPeak, nDead);

  // per-VMM dbcid histograms, so anchors/windows can be re-derived without rerunning
  TFile *fRoot = new TFile(Form("%s/root/fttDataWindow.%d.root", outdir, run), "RECREATE");
  for (int sv = 0; sv < kNVMM; sv++) hV[sv]->Write();
  hAligned->Write(); hTvsV->Write(); hAnchor->Write(); hN->Write(); hFrac->Write();
  fRoot->Close();
}
