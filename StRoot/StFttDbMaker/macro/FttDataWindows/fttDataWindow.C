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
// default [-1, 2]. Rationale: the measurement is per VMM and per run, from a few
// hundred events, so a VMM whose crossing is a tick thin for want of statistics
// must not end up with a cut TIGHTER than the number we trust. A window that is
// one tick too wide costs a little pileup; one that is a tick too narrow throws
// real hits away, in production, silently. The fit therefore only ever tells us
// where to open the window FURTHER.
//
// WHY [-1, 2] AND NOT [-1, 3]. One bunch crossing is 4.263 ticks, so [-1,2] is four
// ticks -- narrower than the period, and unable to hold two crossings. [-1,3] is
// five, and must hold part of a neighbour for some events. Measured per tick,
// summed over every calibrated VMM in each VMM's own anchor frame (signal is
// background-subtracted with the same phase-matched model the fit uses):
//
//                 pp 23081049        AuAu 24192017      zero field 25174011
//   t         signal    bkg        signal    bkg        signal     bkg
//   -1         18165   7887         38068    859        368870    7027
//    0         64610  15556         93254    895        805234    9793
//   +1         15342   6821         26784    496        318857    6037
//   +2          2411   1512          7145    333         33000    4600
//   +3          -660   4795          5477    380          6676    4956
//
// Tick +3 is not signal on the pp field-on data the [-1,3] number came from: 4135
// raw counts against 4795 of background, negative after subtraction. On zero field
// it is 0.4% of the crossing. Only AuAu shows a real +3.3%, and there the tail
// keeps going past it (+4 2.0%, +5 1.0%, +6 0.9%), so +3 is an arbitrary stopping
// point rather than an edge. In all three it costs 15-18% more background.
//
// The older [-1,3] came from the summed-over-all-VMM distribution of
// ftt_sim_hit_maker/timecut section 5, which is broadened by the per-VMM anchor
// offsets measured since -- some VMMs sit up to 12 ticks from their FOB mates. Per
// VMM the crossing really is the three ticks -1, 0, +1.
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
  // The most populated dbcid, full stop -- no grouping of any kind.
  //
  // This is what the consumer does: HitCalibHelper::add() keeps a running argmax
  // over single dbcid values ("Anchor = the most populated dbcid"). The macro has
  // always claimed the same definition, but it used to group 4 bins first, and that
  // grouping was the bug behind VMM 54 of run 23081049: its real peak, 60+60 counts
  // at dbcid 1070,1071 over a level of 25 (9.9 sigma, not a fluctuation), straddles
  // a group boundary, so its group summed 170 while two neighbours holding only the
  // broad tail summed 173 and 171. The anchor came out 1075, five ticks late. The
  // raw argmax was 1070 throughout -- the same as all three of its FOB mates.
  //
  // Grouping was meant to help when statistics are thin. Measured, it does not:
  // counting the VMMs whose anchor disagrees with the median of their FOB mates,
  // the plain argmax against a sliding 4-tick sum gives
  //
  //   23081049  pp field on, 4063 hits/VMM    17 / 360  vs  17 / 360
  //   24192017  AuAu,        ~1000 hits/VMM   29 / 354  vs  29 / 354
  //   25174011  zero field,  7206 hits/VMM    88 / 344  vs  88 / 344
  //   25174011  zero field,   522 hits/VMM    92 / 330  vs  93 / 330
  //
  // identical everywhere except the thinnest case, where the plain argmax is
  // slightly better. So there is nothing to buy and a boundary artefact to avoid.
  int nb = h->GetNbinsX();
  int best = 0;
  double bestC = -1;
  for (int db = 0; db < nb; db++) {
    double c = h->GetBinContent(db + 1);
    if (c > bestC) { bestC = c; best = db; }
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
                   int    winFloorMax = 2,      // [winFloorMin, winFloorMax] -- see below
                   const char *vmmQA  = "")     // up to 8 VMMs to also histogram per
                                                // CHANNEL, e.g. "54,172,372"
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

  // ---- redo mode: the input is one of OUR OWN root/fttDataWindow.<run>.root -----
  // Everything the analysis needs per VMM is the raw dbcid histogram, so a saved
  // file can be re-analysed with no MuDst: anchors, acceptance and window are all
  // recomputed. That is how a calibration is redone after the algorithm changes,
  // and how someone without the MuDsts can still produce a calibration text file.
  //
  // Two things are not in the histogram file and are read from the .txt beside it
  // (../fttDataWindow.<run>.txt, as fttDataWindow.C writes them): the run metadata
  // for the header, and each VMM's strip orientation for the QA columns. Without
  // that .txt the run number still comes from the filename and the analysis is
  // unaffected; only the orientation split in the QA line is lost.
  // The histograms can also come from the forward QA job
  // (StFwdRun22Qa_<run>_raw_*.root), which keeps dbcid vs CHANNEL per VMM as
  // H2S_FttRawHit_dbcidVch_<vmm>, 64 x 4096. Summed over channel that IS the
  // hDbcid_<vmm> this macro saves: the same StFttRawHit::dbcid(), the same
  // StFttDb::vmmId() 0..383, and no time cut applied upstream. Verified on run
  // 23081049 -- anchors derived from the QA file match the MuDst calibration for
  // all 360 calibrated VMMs exactly, on 6.4x the statistics (13.4M hits vs 2.0M).
  // 376 of 384 VMMs come out with identical DB columns; the 8 that differ are all
  // statistics (3 VMMs rise above the sig>=200 cut, 1 window closes from +3 to +2,
  // and one FOB is rejected as noisy on the larger sample).
  bool redoMode = false, qaMode = false;
  TFile *fIn = 0;
  if (in.EndsWith(".root")) {
    fIn = TFile::Open(in.Data());
    if (fIn && !fIn->IsZombie() && fIn->Get("hDbcid_000")) redoMode = true;
    else if (fIn && !fIn->IsZombie() && fIn->Get("H2S_FttRawHit_dbcidVch_0")) {
      redoMode = true; qaMode = true;
    }
    else if (fIn) { fIn->Close(); fIn = 0; }
    // back to memory, or the accumulators below would be created inside fIn and
    // collide by name with what is already in it
    gROOT->cd();
  }

  StChain *chain = 0;
  int run = 0, date = 0, itime = 0;
  Long64_t nEntries = 0;
  StMuDstMaker *muDstMaker = 0;

  if (!redoMode) {
    chain = new StChain("StChain");
    muDstMaker = new StMuDstMaker(0, 0, "", files[0].Data(), "MuDst.root", 1);
    TChain &muChain = *muDstMaker->chain();
    for (int iF = 1; iF < nFiles; iF++) muChain.Add(files[iF].Data());
    // Only the event header and FTT raw hits are needed -- a large speed-up on
    // physics-stream files.
    muDstMaker->SetStatus("*", 0);
    muDstMaker->SetStatus("MuEvent", 1);
    muDstMaker->SetStatus("FttRawHit", 1);

    St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
    new StMuDst2StEventMaker();
    new StFttDbMaker();

    if (chain->Init()) { printf("fttDataWindow: chain Init failed\n"); return; }
    if (muChain.GetEntries() <= 0) {
      printf("fttDataWindow: ERROR no events readable from %s (missing file or unavailable disk?)\n", input);
      return;
    }

    muChain.GetEntry(0);
    run = muDstMaker->muDst()->event()->runNumber();
    // UTC date/time of the first event, for the DB query. TTimeStamp rather than
    // gmtime(): CINT does not reliably know struct tm without extra includes.
    TTimeStamp firstTs((time_t)muDstMaker->muDst()->event()->eventInfo()->time(), 0);
    date  = (int)firstTs.GetDate(kTRUE);
    itime = (int)firstTs.GetTime(kTRUE);
    dbMk->SetDateTime(date, itime);
    chain->InitRun(run);
    nEntries = (Long64_t) muChain.GetEntries();
    printf("fttDataWindow: run %d, first event %d %06d, %d file(s), %lld entries\n",
           run, date, itime, nFiles, nEntries);
  }

  // ---- accumulators -------------------------------------------------------------
  TH1I *hV[384];
  for (int iv = 0; iv < kNVMM; iv++) {
    hV[iv] = new TH1I(Form("hDbcid_%03d", iv), Form("VMM %d;dbcid;hits", iv), kNDBCID, 0, kNDBCID);
    hV[iv]->SetDirectory(0);
  }
  double orientCount[384][5];
  for (int io = 0; io < kNVMM; io++) for (int jo = 0; jo < kNOrient; jo++) orientCount[io][jo] = 0;

  // Optional: dbcid per CHANNEL for a few named VMMs, so that "is this VMM shifted as
  // a whole, or only some of its channels?" can be asked without another pass over
  // the MuDsts. 64 channels x 4096 dbcid for every VMM would be far too much to keep,
  // which is why it is opt-in.
  int qaVmm[8], nQa = 0;
  TH2S *hCh[8];
  for (int iq = 0; iq < 8; iq++) { qaVmm[iq] = -1; hCh[iq] = 0; }
  {
    TString qs(vmmQA);
    TObjArray *toks = qs.Tokenize(",");
    for (int it = 0; it < toks->GetEntries() && nQa < 8; it++) {
      TString one = ((TObjString *)toks->At(it))->GetString();
      one.ReplaceAll(" ", "");
      if (one.Length() == 0) continue;
      int vq = one.Atoi();
      if (vq < 0 || vq >= kNVMM) continue;
      qaVmm[nQa] = vq;
      hCh[nQa] = new TH2S(Form("hChDbcid_%03d", vq),
                          Form("VMM %d;channel;dbcid", vq), 64, -0.5, 63.5, kNDBCID, 0, kNDBCID);
      hCh[nQa]->SetDirectory(0);
      nQa++;
    }
    delete toks;
    if (nQa > 0) printf("fttDataWindow: also histogramming dbcid per channel for %d VMM(s)\n", nQa);
  }

  int anchor[384], status[384];
  double nHitsV[384], sigV[384], bkgV[384], inWinV[384];

  // nEntries was set from the chain above; redo mode leaves it 0 and never loops
  if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;
  double nHitsTot = 0, nBadVmm = 0, nBadDbcid = 0, nOtherRun = 0;
  int nRead = 0;
  TString stopReason = "endOfInput";
  double nHitsTotRedo = 0;

  // ---- redo mode: fill the accumulators from the saved file instead of events ----
  if (redoMode) {
    for (int ir = 0; ir < kNVMM; ir++) {
      if (qaMode) {
        TH2 *h2In = (TH2 *) fIn->Get(Form("H2S_FttRawHit_dbcidVch_%d", ir));
        if (!h2In) continue;
        if (h2In->GetNbinsY() != kNDBCID) {
          printf("fttDataWindow: QA VMM %d has %d dbcid bins, expected %d -- skipped\n",
                 ir, h2In->GetNbinsY(), kNDBCID);
          continue;
        }
        double sumV = 0;
        for (int iy = 1; iy <= kNDBCID; iy++) {
          double cy = 0;
          for (int ix = 1; ix <= h2In->GetNbinsX(); ix++) cy += h2In->GetBinContent(ix, iy);
          if (cy > 0) hV[ir]->SetBinContent(iy, (int) (cy + 0.5));
          sumV += cy;
        }
        hV[ir]->SetEntries(sumV);
        nHitsTotRedo += sumV;
        continue;
      }
      TH1 *hIn = (TH1 *) fIn->Get(Form("hDbcid_%03d", ir));
      if (!hIn) continue;
      for (int ib = 1; ib <= kNDBCID; ib++) hV[ir]->SetBinContent(ib, hIn->GetBinContent(ib));
      hV[ir]->SetEntries(hIn->GetEntries());
      nHitsTotRedo += hIn->GetEntries();
    }
    // Run metadata and the per-VMM orientation come from the .txt beside the file,
    // ../fttDataWindow.<run>.txt. Both are QA only -- the calibration itself is
    // recomputed from the histograms either way.
    // The run number: from the .txt header when there is one, otherwise from the
    // file name -- fttDataWindow.<run>.root, or the QA job's
    // StFwdRun22Qa_<run>_raw_<n>.root, where it is an embedded 8-digit field.
    int runFromName = 0;
    TString bn = gSystem->BaseName(in.Data());
    TString b1 = bn; b1.ReplaceAll("fttDataWindow.", ""); b1.ReplaceAll(".root", "");
    if (b1.IsDigit()) runFromName = b1.Atoi();
    int nbn = bn.Length();
    for (int ic = 0; ic + 8 <= nbn && runFromName == 0; ic++) {
      if (ic > 0      && bn[ic-1] >= '0' && bn[ic-1] <= '9') continue;  // mid-number
      if (ic + 8 < nbn && bn[ic+8] >= '0' && bn[ic+8] <= '9') continue;  // longer than 8
      TString cand = bn(ic, 8);
      if (cand.IsDigit() && cand.Atoi() > 10000000) runFromName = cand.Atoi();
    }

    TString dirIn  = gSystem->DirName(in.Data());      // .../root
    TString baseIn = gSystem->BaseName(in.Data());     // fttDataWindow.<run>.root
    baseIn.ReplaceAll(".root", ".txt");
    // A QA file's own name is not fttDataWindow.<run>.txt, so look for the
    // calibration text file by run number instead -- orientation and metadata are
    // picked up if a previous calibration of this run is sitting there.
    if (qaMode && runFromName > 0) baseIn = Form("fttDataWindow.%d.txt", runFromName);
    TString txt = Form("%s/../%s", dirIn.Data(), baseIn.Data());   // the usual layout
    FILE *T = fopen(txt.Data(), "r");
    if (!T) {                                          // or sitting beside the file
      txt = Form("%s/%s", dirIn.Data(), baseIn.Data());
      T = fopen(txt.Data(), "r");
    }
    if (T) {
      char ln[512];
      while (fgets(ln, sizeof(ln), T)) {
        if (ln[0] == '#') {
          int r2 = 0, d2 = 0, t2 = 0;
          if (sscanf(ln, "# run %d firstEventDate %d firstEventTime %d", &r2, &d2, &t2) == 3) {
            run = r2; date = d2; itime = t2;
          }
          continue;
        }
        // cols 1, 7 and 11 are uuid, nHits and orient. The rest are skipped with %*,
        // which costs no argument: CINT's sscanf takes at most 12.
        int u2 = -1, or2 = -1; double nn2 = 0;
        if (sscanf(ln, "%d %*d %*d %*d %*d %*d %lf %*f %*f %*f %d", &u2, &nn2, &or2) == 3)
          if (u2 >= 0 && u2 < kNVMM && or2 >= 0 && or2 < kNOrient) orientCount[u2][or2] = nn2;
      }
      fclose(T);
      printf("fttDataWindow: redo -- metadata and orientation from %s\n", txt.Data());
    } else {
      printf("fttDataWindow: redo -- no .txt beside the histogram file; the QA\n"
             "               orientation split and the run metadata will be blank\n");
    }
    if (run == 0) run = runFromName;   // no .txt header: take it from the filename
    nRead = 0;
    nHitsTot = nHitsTotRedo;
    stopReason = "redoFromHistograms";
    printf("fttDataWindow: REDO run %d from %s%s -- %0.f hits in 384 VMM histograms,\n"
           "               anchors and windows recomputed, no MuDst read\n",
           run, in.Data(),
           qaMode ? " (forward QA file, dbcid summed over channel)" : "",
           nHitsTotRedo);
  }

  for (Long64_t iev = 0; !redoMode && iev < nEntries; iev++) {
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
      for (int jq = 0; jq < nQa; jq++)
        if (v == qaVmm[jq]) hCh[jq]->Fill((int)hit->channel(), d);
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
  for (int sq = 0; sq < nQa; sq++) hCh[sq]->Write();
  hAligned->Write(); hTvsV->Write(); hAnchor->Write(); hN->Write(); hFrac->Write();
  fRoot->Close();
}
