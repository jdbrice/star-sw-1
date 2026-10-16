// fttDataWindow_db.C -- push per-run FTT time calibration (from script/fttDataWindow.C)
// into Calibrations/ftt/fttDataWindowsB, and/or read it back.
// Adapted from $STAR/StRoot/StFcsDbMaker/macro/fcsTimedepGainCorr_db.C: kept the
// run -> startRunTime lookup and the StDbManager store / St_db_Maker read-back.
//
// Usage:
//   root4star -l -b -q 'fttDataWindow_db.C(23081049)'                 // read back + compare only
//   root4star -l -b -q 'fttDataWindow_db.C(23081049,"writedb")'       // store, then read back
//   root4star -l -b -q 'fttDataWindow_db.C(0,"writedb")'              // every file in data/
// Options (substrings of opt): writedb  store to DB (needs write access)
//                              force    store even if the calibration looks poor
//
// TABLE CONTENT
//   slots 0..383 : one per VMM, uuid = StFttDb::vmmId() (0-based), mode 0 = calibrated
//                  time, min/max window, anchor (-1 = no usable calibration for that VMM)
//   slot  384    : RUN STAMP -- uuid = -1, mode = 0, min = run / 10000, max = run % 10000,
//                  anchor = stamp format version (1)
// vmmId only spans 0..383, so slot 384 is spare.
//
// WHY A RUN STAMP. STAR DB entries have a begin time only (the end time is left in the
// far future by policy), so St_db_Maker returns the most recent entry at or before the
// query time. A run with no entry of its own would therefore silently get the PREVIOUS
// run's anchors -- and anchors differ run to run, they are not slowly varying. The
// consumer (StFttDbMaker) compares the stamp with the run being reconstructed and uses
// the anchors only on an exact match. Unstamped entries (the pre-existing global one,
// uuid[384] == 384) are treated as window-only.

// Read-only host for the "is it already there?" query. Writes go through StDbLib
// to the load balancer; this is only used to look, so one node is fine.
const char *kReadHost = "star-db01.sdcc.bnl.gov";
const int   kReadPort = 3316;
// The write server of $DB_SERVER_LOCAL_CONFIG. Only used to confirm a row we have
// just written: a replica may not have it yet.
const char *kWriteHost = "robinson.star.bnl.gov";
const int   kWritePort = 3306;

// A run's first entry goes kPreRoll seconds BEFORE the run start, not at it.
// Resolution is "latest beginTime at or before the requested time" and every
// event of the run is at or after the run start, so an entry anywhere in
// [start-kPreRoll, start] serves the WHOLE run -- while an entry after the
// start would leave the first seconds of the run on the previous entry. That
// gives kPreRoll+1 slots to recalibrate a run into before anything has to be
// deactivated, each revision one second later than the last.
//
// 10 s is safe against landing inside the previous run: measured over 2149
// runs of 2022 days 060-099, the smallest gap between one run ending and the
// next starting was 26 s, and none was below 10 s.
const int kPreRoll = 10;

const int kNVMMdb     = 384;
const int kStampSlot  = 384;
const int kStampUuid  = -1;
const int kStampVers  = 1;

// Parsed calibration, held in plain arrays: CINT resolves function SIGNATURES when the
// file is loaded -- before loadDbLibs() makes fttDataWindowsB_st known -- so the struct
// may only appear inside function bodies.
short gUuid[385];
short gMode[385];
short gMin[385];
short gMax[385];
short gAnchor[385];
int   gNCal = 0, gNActive = 0, gFileRun = -1;

// One St_db_Maker reused for every run. Creating one per run costs ~30-60 s, which
// is minutes of MySQL setup over a 56-run pass. Held as TObject* because CINT
// resolves signatures at load time, before loadDbLibs() makes St_db_Maker known.
TObject *gDbMk = 0;

int  dbLoaded = 0;
void loadDbLibs()
{
  if (dbLoaded) return;
  gROOT->Macro("LoadLogger.C");
  gSystem->Load("St_base.so");
  gSystem->Load("libStDb_Tables.so");
  gSystem->Load("StDbLib.so");
  gSystem->Load("StChain.so");
  gSystem->Load("StBFChain");
  gSystem->Load("StUtilities");
  gSystem->Load("StIOMaker");
  gSystem->Load("StarClassLibrary");
  gSystem->Load("St_Tables");
  gSystem->Load("StDbLib");
  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");
  dbLoaded = 1;
}

// run start from RunLog (same host/port rule as fcsTimedepGainCorr_db.C)
int runStartTime(int run)
{
  int year = run / 1000000 - 1;
  int port = 3400 + year - 1;
  TString cmd = Form("mysql -h db04.star.bnl.gov --port=%d -N -s -e \"SELECT startRunTime FROM RunLog.runDescriptor WHERE runNumber=%d LIMIT 1\"",
                     port, run);
  TString st = gSystem->GetFromPipe(cmd.Data());
  return st.Atoi();
}

// UTC "YYYY-MM-DD HH:MM:SS" for a unix time -- the form beginTime is stored in.
TString utcStamp(int t)
{
  return gSystem->GetFromPipe(Form("date -u -d @%d +'%%Y-%%m-%%d %%H:%%M:%%S'", t));
}

// Unix time of a UTC "YYYY-MM-DD HH:MM:SS".
int unixOfStamp(const char *ts)
{
  return gSystem->GetFromPipe(Form("date -u -d '%s UTC' +%%s", ts)).Atoi();
}

// beginTime of the entry that applies at this run -- the latest at or before the
// run start. Empty if there is none. WHICH run it belongs to is not decided here:
// dbCompare reads the run stamp, which is the reliable answer.
// Latest active beginTime at or before starttime. onMaster picks the write server,
// which is the only place a row written seconds ago is certain to be visible.
TString dbLatestBegin(int starttime, const char *flavor, int onMaster = 0)
{
  TString ts = utcStamp(starttime);
  TString cmd = Form("mysql -h %s --port=%d %s -N -s -e "
                     "\"SELECT MAX(beginTime) FROM Calibrations_ftt.fttDataWindowsB "
                     "WHERE beginTime<='%s' AND flavor='%s' AND deactive=0\" 2>/dev/null",
                     onMaster ? kWriteHost : kReadHost,
                     onMaster ? kWritePort : kReadPort,
                     onMaster ? "" : "-u loader", ts.Data(), flavor);
  TString out = gSystem->GetFromPipe(cmd.Data());
  out = out.Strip(TString::kBoth);
  if (out == "NULL") out = "";
  return out;
}

// Parse fttDataWindow.<run>.txt into the g* arrays. Returns VMM lines read, <0 on error.
int readCalibFile(const char *fname, int run)
{
  FILE *F = fopen(fname, "r");
  if (!F) { printf("ERROR: cannot open %s\n", fname); return -1; }

  for (int j = 0; j < kNVMMdb; j++) {          // defaults for any VMM missing from the file
    gUuid[j] = j; gMode[j] = 0; gMin[j] = -65; gMax[j] = 100; gAnchor[j] = -1;
  }
  gUuid[kStampSlot]   = kStampUuid;
  gMode[kStampSlot]   = 0;
  gMin[kStampSlot]    = run / 10000;
  gMax[kStampSlot]    = run % 10000;
  gAnchor[kStampSlot] = kStampVers;

  char line[512];
  int nLines = 0, nBad = 0;
  int seen[384];
  for (int js = 0; js < kNVMMdb; js++) seen[js] = 0;
  gNCal = 0; gNActive = 0; gFileRun = -1;
  while (fgets(line, sizeof(line), F)) {
    if (line[0] == '#') {
      int r = 0;
      if (sscanf(line, "# run %d", &r) == 1) gFileRun = r;
      continue;
    }
    int u, m, lo, hi, an, stt = -1;
    double nh = -1;
    int nf = sscanf(line, "%d %d %d %d %d %d %lf", &u, &m, &lo, &hi, &an, &stt, &nh);
    if (nf < 5) continue;
    if (u < 0 || u >= kNVMMdb || m < 0 || m > 1 || lo >= hi || lo < -32768 || hi > 32767 || an < -1 || an > 4095) {
      printf("WARNING: %s: bad line (skipped): %s", fname, line);
      nBad++;
      continue;
    }
    if (seen[u]) printf("WARNING: %s: uuid %d appears twice, last one wins\n", fname, u);
    seen[u] = 1;
    gUuid[u] = u; gMode[u] = m; gMin[u] = lo; gMax[u] = hi; gAnchor[u] = an;
    nLines++;
    if (an >= 0) gNCal++;
    if (nf >= 7 && nh > 0) gNActive++;        // QA column: hits seen for this VMM
    else if (nf < 7 && an >= 0) gNActive++;   // no QA columns: count calibrated ones
  }
  fclose(F);
  int nMissing = 0;
  for (int jm = 0; jm < kNVMMdb; jm++) if (!seen[jm]) nMissing++;
  if (nMissing) printf("WARNING: %s: %d VMMs missing, filled with defaults and anchor -1\n", fname, nMissing);
  if (nBad)     printf("WARNING: %s: %d malformed lines skipped\n", fname, nBad);
  return nLines;
}

// Read back what St_db_Maker returns at a given time and compare with the parsed
// file. Returns  -2 nothing there, -1 present but not this run's (unstamped, or
// stamped for another run), 0 identical, >0 number of differing slots.
//
// The stamp makes this self-validating: "stamped for run R and every slot equal
// to the file for run R" can only be run R's own entry, so a cached or stale
// read cannot make us skip a write we needed -- it would show the wrong stamp
// and we would write.
int dbCompare(int run, int date, int time, const char *flavor, int verbose)
{
  St_db_Maker *dbMk = (St_db_Maker *)gDbMk;
  if (!dbMk) {
    dbMk = new St_db_Maker("dbFtt", "MySQL:StarDb", "$STAR/StarDb");
    gDbMk = dbMk;
    dbMk->SetFlavor(flavor);
    dbMk->SetDateTime(date, time);
    dbMk->Init();
  } else {
    dbMk->SetFlavor(flavor);
    dbMk->SetDateTime(date, time);
  }
  dbMk->Make();

  TDataSet *DB = dbMk->GetInputDB("Calibrations/ftt");
  if (!DB) { if (verbose) printf("READBACK run %d: no Calibrations/ftt\n", run); return -2; }
  St_fttDataWindowsB *t = (St_fttDataWindowsB *)DB->Find("fttDataWindowsB");
  if (!t) { if (verbose) printf("READBACK run %d: no fttDataWindowsB\n", run); return -2; }
  fttDataWindowsB_st *s = t->GetTable();

  bool stamped  = (s[0].uuid[kStampSlot] == kStampUuid);
  int  stampRun = stamped ? s[0].min[kStampSlot] * 10000 + s[0].max[kStampSlot] : -1;
  int  nAnch    = 0;
  for (int ja = 0; ja < kNVMMdb; ja++) if (s[0].anchor[ja] >= 0) nAnch++;

  if (verbose) {
    if (stamped)
      printf("READBACK run %d at %d %06d: stamped for run %d%s, VMMs with anchor %d/%d\n", run, date, time,
             stampRun, stampRun == run ? " (MATCH)" : " (DIFFERENT RUN -- consumer will ignore its anchors)",
             nAnch, kNVMMdb);
    else
      printf("READBACK run %d at %d %06d: UNSTAMPED legacy entry -- window only, its anchor column is not calibration and is ignored\n",
             run, date, time);
  }
  if (!stamped || stampRun != run) return -1;

  int nDiff = 0;
  for (int jc = 0; jc <= kStampSlot; jc++) {
    if (s[0].uuid[jc] != gUuid[jc] || (short)s[0].mode[jc] != gMode[jc] ||
        s[0].min[jc]  != gMin[jc]  || s[0].max[jc]  != gMax[jc] ||
        s[0].anchor[jc] != gAnchor[jc]) {
      if (verbose && nDiff < 5)
        printf("  DIFF slot %3d  db(u=%d m=%d %d..%d a=%d)  file(u=%d m=%d %d..%d a=%d)\n", jc,
               s[0].uuid[jc], (int)s[0].mode[jc], s[0].min[jc], s[0].max[jc], s[0].anchor[jc],
               gUuid[jc], (int)gMode[jc], gMin[jc], gMax[jc], gAnchor[jc]);
      nDiff++;
    }
  }
  if (verbose)
    printf("READBACK run %d: %s (%d of %d slots differ)\n", run,
           nDiff == 0 ? "DB content IDENTICAL to file" : "DB content DIFFERS from file", nDiff, kStampSlot + 1);
  return nDiff;
}

int doRun(int run, const char *dir, TString option, const char *flavor)
{
  TString fname = Form("%s/fttDataWindow.%d.txt", dir, run);
  int n = readCalibFile(fname.Data(), run);
  if (n < 0) return -1;
  int nCal = gNCal, nActive = gNActive, fileRun = gFileRun;

  int starttime = runStartTime(run);
  if (starttime <= 0) { printf("ERROR: run %d: no startRunTime from RunLog -- skipped\n", run); return -1; }
  int readDate = gSystem->GetFromPipe(Form("date -u -d @%d +%%Y%%m%%d", starttime)).Atoi();
  int readTime = gSystem->GetFromPipe(Form("date -u -d @%d +%%H%%M%%S", starttime)).Atoi();

  printf("run %d: startRunTime %d (%d %06d UTC)  file %s: %d VMM lines, %d calibrated / %d active\n",
         run, starttime, readDate, readTime, fname.Data(), n, nCal, nActive);

  bool ok = true;
  if (fileRun != run) { printf("ERROR: run %d: file header says run %d\n", run, fileRun); ok = false; }
  if (nActive > 0 && nCal < 0.5 * nActive) {
    printf("ERROR: run %d: only %d of %d active VMMs calibrated\n", run, nCal, nActive);
    ok = false;
  }

  if (!option.Contains("writedb")) {           // look only -- say what a push would do
    int cmpL = dbCompare(run, readDate, readTime, flavor, 1);
    if (!ok && !option.Contains("force"))
      printf("  would NOT store run %d: fails the quality gate\n", run);
    else if (cmpL == 0)
      printf("  nothing to do for run %d: already in the DB and identical\n", run);
    else if (cmpL > 0)
      printf("  run %d is in the DB but DIFFERS (%d slots): a push needs \"update\"\n", run, cmpL);
    else
      printf("  would store run %d at %s (run start %s minus %d s)\n", run,
             utcStamp(starttime - kPreRoll).Data(), utcStamp(starttime).Data(), kPreRoll);
    return 0;
  }

  if (!ok && !option.Contains("force")) {
    printf("NOT WRITING run %d (add \"force\" to override)\n", run);
    return -1;
  }

  // Do not insert blindly. The run stamp in slot 384 says whether what the DB
  // already serves for this run is this run's entry, and dbCompare reads it.
  int cmp = dbCompare(run, readDate, readTime, flavor, 0);
  int storeAt = starttime - kPreRoll;        // first entry for a run

  if (cmp == 0) {
    printf("SKIP run %d: already in the DB and identical\n", run);
    return 0;
  }
  if (cmp > 0) {                             // this run's entry, but different numbers
    TString prev = dbLatestBegin(starttime, flavor);
    if (prev == "") { printf("NOT WRITING run %d: cannot find the existing entry\n", run); return -1; }
    if (!option.Contains("update")) {
      printf("SKIP run %d: it is in the DB already and the content DIFFERS in %d of %d slots.\n"
             "  The existing entry begins %s. Nothing is overwritten -- "
             "(nodeID,elementID,beginTime,flavor,deactive)\n"
             "  is the primary key. Add \"update\" to store the new numbers one second later,\n"
             "  which then wins for the whole run. See README.\n",
             run, cmp, kStampSlot + 1, prev.Data());
      dbCompare(run, readDate, readTime, flavor, 1);
      return -1;
    }
    storeAt = unixOfStamp(prev.Data()) + 1;
    if (storeAt > starttime) {
      printf("NOT WRITING run %d: the %d slots from %s to the run start are all used.\n"
             "  Going later than the run start would leave the first seconds of the run on\n"
             "  the old entry, so ask the database team to deactivate one instead.\n",
             run, kPreRoll + 1, utcStamp(starttime - kPreRoll).Data());
      return -1;
    }
    printf("UPDATE run %d: %d slots differ, storing a revision at %s\n", run, cmp, utcStamp(storeAt).Data());
  }

  // DB_ACCESS_MODE is set once in fttDataWindow_db(), before anything connects --
  // see the note there. Setting it here would be too late.
  StDbManager *mgr = StDbManager::Instance();
  StDbConfigNode *node = mgr->initConfig("Calibrations_ftt");
  mgr->setStoreTime(storeAt);
  fttDataWindowsB_st dw;                 // struct only inside a body (see note at top)
  memset(&dw, 0, sizeof(dw));
  for (int jw = 0; jw <= kStampSlot; jw++) {
    dw.uuid[jw] = gUuid[jw]; dw.mode[jw] = (unsigned char)gMode[jw];
    dw.min[jw] = gMin[jw]; dw.max[jw] = gMax[jw]; dw.anchor[jw] = gAnchor[jw];
  }
  StDbTable *table = node->addDbTable("fttDataWindowsB");
  table->SetTable((char *)&dw, 1);
  table->setFlavor(flavor);
  bool stored = mgr->storeDbTable(table);
  if (!stored) {
    printf("\n*** WRITE FAILED for run %d at %s -- NOTHING WAS STORED ***\n"
           "    see the MysqlDb/StDbSql errors above\n", run, utcStamp(storeAt).Data());
    return -2;                           // fatal: stop, do not try the other runs
  }
  printf("STORE run %d at %s flavor %s: OK\n", run, utcStamp(storeAt).Data(), flavor);

  // NOT dbCompare() here. St_db_Maker caches the dataset it read before the store, so
  // asking it again in this process returns the state from BEFORE the write and looks
  // like the write did nothing. Confirm the row directly instead, on the server it was
  // written to, and leave the slot-by-slot comparison to a fresh read-only run.
  TString now = dbLatestBegin(storeAt, flavor, 1);
  TString want = utcStamp(storeAt);
  if (now == want)
    printf("VERIFY run %d: the entry at %s is in %s. Re-run without \"writedb\" for the\n"
           "       slot-by-slot comparison -- this process cannot, its read is cached.\n",
           run, want.Data(), kWriteHost);
  else
    printf("VERIFY run %d: expected an entry at %s but %s reports %s -- CHECK THIS\n",
           run, want.Data(), kWriteHost, now == "" ? "none" : now.Data());
  return 0;
}

void fttDataWindow_db(int run = 0, const char *opt = "", const char *dir = "",
                      const char *flavor = "ofl")
{
  loadDbLibs();
  TString option(opt);

  // Default to the data/ directory beside this macro, so it works from any cwd.
  // Pass dir explicitly to read a calibration you have just produced elsewhere.
  TString sdir(dir);
  if ( sdir == "" ) sdir = Form("%s/data", gSystem->DirName(__FILE__));

  printf("fttDataWindow_db: run=%d opt=\"%s\" dir=%s flavor=%s  (writedb=%d)\n",
         run, opt, sdir.Data(), flavor, (int)option.Contains("writedb"));

  // DB_ACCESS_MODE has to be set BEFORE the first connection, not just before the
  // store. StDbManager opens one connection and reuses it, and this macro reads the
  // existing entry before deciding to write -- so if the mode is set only around
  // storeDbTable, the store goes out on the connection the read already opened, to a
  // replica, and MySQL answers
  //
  //   The MySQL server is running with the --read-only option so it cannot execute
  //   this statement
  //
  // with the macro correctly reporting that nothing was stored. (fstGeometry_db.C
  // never hit this: it only writes, so its first connection is already a write one.)
  //
  // Setting it here means the read-back also goes to the write server named in
  // $DB_SERVER_LOCAL_CONFIG, which is robinson.star.bnl.gov:3306. That is the master,
  // so the "is it already there?" check sees our own writes immediately rather than
  // waiting for replication -- worth more than the trivial extra load of one query
  // per run. A read-only invocation leaves the mode alone and uses the replicas.
  if (option.Contains("writedb")) {
    gSystem->Setenv("DB_ACCESS_MODE", "write");
    printf("fttDataWindow_db: DB_ACCESS_MODE=write (writes and read-backs both go to\n"
           "                  the write server in %s)\n",
           gSystem->Getenv("DB_SERVER_LOCAL_CONFIG") ? gSystem->Getenv("DB_SERVER_LOCAL_CONFIG")
                                                     : "$DB_SERVER_LOCAL_CONFIG");
  }

  if (run > 0) { doRun(run, sdir.Data(), option, flavor); return; }

  TString list = gSystem->GetFromPipe(Form("ls %s/fttDataWindow.*.txt 2>/dev/null | sort", sdir.Data()));
  TObjArray *tok = list.Tokenize("\n");
  int nDone = 0;
  for (int it = 0; it < tok->GetEntries(); it++) {
    TString f = ((TObjString *)tok->At(it))->GetString();
    TString base = gSystem->BaseName(f.Data());
    base.ReplaceAll("fttDataWindow.", "");
    base.ReplaceAll(".txt", "");
    if (!base.IsDigit()) continue;
    if (doRun(base.Atoi(), sdir.Data(), option, flavor) == -2) {
      printf("*** STOPPING -- %d run(s) processed before the failure ***\n", nDone);
      return;
    }
    nDone++;
  }
  printf("fttDataWindow_db: processed %d run(s) from %s\n", nDone, sdir.Data());
}
