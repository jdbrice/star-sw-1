// fttDataWindowCheck.C -- what does the CONSUMER get for a run?
//
// fttDataWindow_db.C says what the database serves at a run's timestamp and whether
// it matches the file. This says what StFttDbMaker then does with it, which is the
// part that decides whether production uses the anchors:
//
//   stamped for THIS run     -> per-VMM anchors ENABLED
//   stamped for another run  -> anchors IGNORED, windows still applied, and
//                               StFttHitCalibMaker derives anchors on the fly
//   unstamped (legacy entry) -> same, quietly
//
//   root4star -l -b -q 'fttDataWindowCheck.C(25174011)'
//
// WHERE TO RUN IT. Two requirements, both the same as a condor job:
//
//   - NO ./fttDataWindow/ directory. A local fttDataWindow/fttDataWindow.<run>.txt
//     overrides the database, so a tree that has one tests the file, not the DB.
//     The macro checks for it and refuses.
//   - YOUR OWN BUILD, i.e. .$STAR_HOST_SYS reachable. Without it root4star loads
//     the official $STAR StFttDbMaker, which is older (2026-06-03) and has no run
//     stamp at all -- no dataWindowStampRun, no setDataWindowAnchorsValid -- so it
//     applies whatever entry it finds and prints none of the lines this macro is
//     looking for. The check then silently tests the wrong code.
//
//       mkdir -p /tmp/dbcheck && cd /tmp/dbcheck
//       ln -s <tree>/.$STAR_HOST_SYS . && ln -s <tree>/StRoot .
//
//     StRoot itself is optional here: StFttDbMaker::InitRun reads
//     StRoot/StFttDbMaker/Row*.txt first, but a missing one is only a warning
//     ("Strip center file not found") and the data windows load regardless. Link it
//     anyway to keep the log clean.
//
// Measured example of what it is for -- run 23081001 pushed, 23081002 not:
//
//   run 23081001  entry begins 04:08:18, query 04:08:28
//                 per-VMM time anchors ENABLED: entry stamped for run 23081001
//   run 23081002  the SAME entry is served, it is the latest at or before 04:40:20
//                 WARNING ... stamped for run 23081001 but this is run 23081002:
//                 its anchors are IGNORED ... derives anchors on the fly
//
// which is the whole point of the run stamp: without it 23081002 would silently
// reconstruct with 23081001's anchors.
//
// CINT notes: StFttDb may only appear inside a function body -- signatures are
// resolved when the file is loaded, before loadLibs() makes the class known.

void loadLibsChk()
{
  gROOT->LoadMacro(Form("%s/fttDataWindow.C", gSystem->DirName(__FILE__)));
  gROOT->LoadMacro(Form("%s/fttDataWindow_db.C", gSystem->DirName(__FILE__)));
  loadLibsDW();
}

// The run start comes from fttDataWindow_db.C's runStartTime(), loaded rather than
// copied: the RunLog port depends on the run PERIOD, which does not line up with the
// run-number prefix in December, and that rule belongs in one place.

void fttDataWindowCheck(int run, const char *flavor = "ofl")
{
  loadLibsChk();

  TString local = Form("fttDataWindow/fttDataWindow.%d.txt", run);
  if (!gSystem->AccessPathName(local.Data())) {
    printf("\n*** %s EXISTS in this directory. StFttDbMaker prefers it over the\n"
           "*** database, so this would test the FILE, not the DB. Run from elsewhere.\n\n",
           local.Data());
    return;
  }

  int st = runStartTime(run);
  if (st <= 0) { printf("fttDataWindowCheck: no startRunTime for run %d\n", run); return; }
  TTimeStamp ts((time_t)st, 0);
  int date = (int)ts.GetDate(kTRUE), time = (int)ts.GetTime(kTRUE);
  printf("fttDataWindowCheck: run %d, DB timestamp %d %06d (its start)\n", run, date, time);

  StChain *chain = new StChain("chain");
  St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
  dbMk->SetFlavor(flavor);
  dbMk->SetDateTime(date, time);
  StFttDbMaker *fttDbMk = new StFttDbMaker("fttDbMkr");
  if (chain->Init()) { printf("fttDataWindowCheck: chain Init failed\n"); return; }
  // InitRun on the maker itself: StChain::InitRun did not reach it here, and this is
  // the call production makes.
  fttDbMk->InitRun(run);

  // The answer is in StFttDbMaker's two loadDataWindows lines above, which is where
  // production shows it too. They are not read back through an accessor: CINT cannot
  // call StFttDb's inline methods on a pointer it got from the chain.
  printf("\n==== run %d: the two \"StFttDbMaker::loadDataWindows\" lines above say it ====\n", run);
  printf("  ENABLED  -> this run reconstructs with its own measured anchors\n");
  printf("  IGNORED  -> the entry belongs to another run. Correct, and expected for any\n");
  printf("              run with no entry of its own: the windows still apply and\n");
  printf("              StFttHitCalibMaker derives the anchors on the fly.\n");
  printf("  unstamped-> only the old global entry was found for this run\n");
  printf("\n");
}
