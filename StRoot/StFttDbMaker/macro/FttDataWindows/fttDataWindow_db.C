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

// Read back what St_db_Maker returns at a given time; compare with expected if given.
void readBack(int run, int date, int time, const char *flavor, int compare)
{
  St_db_Maker *dbMk = new St_db_Maker(Form("db%d", run), "MySQL:StarDb", "$STAR/StarDb");
  dbMk->SetDateTime(date, time);
  dbMk->SetFlavor(flavor);
  dbMk->Init();
  dbMk->Make();
  TDataSet *DB = dbMk->GetInputDB("Calibrations/ftt");
  if (!DB) { printf("READBACK run %d: no Calibrations/ftt\n", run); return; }
  St_fttDataWindowsB *t = (St_fttDataWindowsB *)DB->Find("fttDataWindowsB");
  if (!t) { printf("READBACK run %d: no fttDataWindowsB\n", run); return; }
  fttDataWindowsB_st *s = t->GetTable();

  bool stamped = (s[0].uuid[kStampSlot] == kStampUuid);
  int stampRun = stamped ? s[0].min[kStampSlot] * 10000 + s[0].max[kStampSlot] : -1;
  int nAnch = 0;
  for (int ja = 0; ja < kNVMMdb; ja++) if (s[0].anchor[ja] >= 0) nAnch++;
  if (stamped)
    printf("READBACK run %d at %d %06d: stamped for run %d%s, VMMs with anchor %d/%d\n", run, date, time,
           stampRun, stampRun == run ? " (MATCH)" : " (DIFFERENT RUN -- consumer will ignore its anchors)",
           nAnch, kNVMMdb);
  else
    printf("READBACK run %d at %d %06d: UNSTAMPED legacy entry -- window only, its anchor column is not calibration and is ignored\n",
           run, date, time);

  if (!compare) return;
  int nDiff = 0;
  for (int jc = 0; jc <= kStampSlot; jc++) {
    if (s[0].uuid[jc] != gUuid[jc] || (short)s[0].mode[jc] != gMode[jc] ||
        s[0].min[jc]  != gMin[jc]  || s[0].max[jc]  != gMax[jc] ||
        s[0].anchor[jc] != gAnchor[jc]) {
      if (nDiff < 5)
        printf("  DIFF slot %3d  db(u=%d m=%d %d..%d a=%d)  file(u=%d m=%d %d..%d a=%d)\n", jc,
               s[0].uuid[jc], (int)s[0].mode[jc], s[0].min[jc], s[0].max[jc], s[0].anchor[jc],
               gUuid[jc], (int)gMode[jc], gMin[jc], gMax[jc], gAnchor[jc]);
      nDiff++;
    }
  }
  printf("READBACK run %d: %s (%d of %d slots differ)\n", run,
         nDiff == 0 ? "DB content IDENTICAL to file" : "DB content DIFFERS from file", nDiff, kStampSlot + 1);
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

  if (option.Contains("writedb")) {
    if (!ok && !option.Contains("force")) {
      printf("NOT WRITING run %d (add \"force\" to override)\n", run);
    } else {
      gSystem->Setenv("DB_ACCESS_MODE", "write");
      StDbManager *mgr = StDbManager::Instance();
      StDbConfigNode *node = mgr->initConfig("Calibrations_ftt");
      mgr->setStoreTime(starttime);
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
      printf("STORE run %d at %d %06d flavor %s: %s\n", run, readDate, readTime, flavor, stored ? "OK" : "FAILED");
    }
  }
  readBack(run, readDate, readTime, flavor, 1);
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
    doRun(base.Atoi(), dir, option, flavor);
    nDone++;
  }
  printf("fttDataWindow_db: processed %d run(s) from %s\n", nDone, dir);
}
