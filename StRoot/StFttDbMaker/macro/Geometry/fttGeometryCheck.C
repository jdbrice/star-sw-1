// fttGeometryCheck.C -- what geometry does the DATABASE ALONE give for a timestamp?
//
//   root4star -l -b -q 'fttGeometryCheck.C(20220322,221026)'   // a Run 22 real-data time
//   root4star -l -b -q 'fttGeometryCheck.C()'                  // same, the default
//
// St_db_Maker is created with MySQL as its ONLY source, which is what the BFC
// option NoCintDb does: no $STAR/StarDb, no ./StarDb. That matters here for the
// same reason it does on the FST side -- a Geometry/stgc or Geometry/ftt directory
// shipped with the release, or sitting in the working directory, wins over the DB
// and the DB is then never exercised. This is the check that says what production
// will get from the database itself.
//
// For each table it prints which entry was served (the validity range), how many
// rows are identity, and the first row's translation, which is where the survey
// shows itself.
//
// Expected, once Geometry_ftt is populated (see ftt/README):
//
//   2021-10-12   pentOnStation row0 t = (+0.80900, +3.63400, 0)   the SURVEY
//   2022-03-22   pentOnStation row0 t = (+2.09845, +4.60535, 0)   Run 22 iter 3
//   2024-02-10   pentOnStation row0 t = (+1.94794, +5.75562, 0)   Run 24 iter 1
//   before 2021-10-10                                             nothing served
//
// All three tables must answer, and with 1 / 4 / 16 rows: StFttDb::loadGeometry()
// chains them, so a timestamp that resolves only some of them is a half-built
// geometry, not a partial one.
//
// CINT note: Survey_st and St_db_Maker may only appear inside a function BODY --
// signatures are resolved when the file is loaded, before loadLibs makes them known.
// That is why this is one function rather than a helper per table.

void loadLibsFttChk()
{
  gSystem->Load("libStarClassLibrary.so");
  gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StarMagField");
  gSystem->Load("StMagF");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbBroker");        // St_db_Maker needs it; missing it is a
  gSystem->Load("StDbUtilities");     // symbol lookup error at load time
  gSystem->Load("St_db_Maker");
  gSystem->Load("libStDb_Tables.so");
}

void fttGeometryCheck(int date = 20220322, int time = 221026)
{
  loadLibsFttChk();
  // MySQL ONLY: no $STAR/StarDb, no ./StarDb -- the BFC NoCintDb case
  St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb");
  dbMk->SetFlavor("ofl");
  dbMk->SetDateTime(date, time);
  dbMk->Init();
  dbMk->Make();

  printf("\nfttGeometryCheck: Geometry/ftt from the DATABASE ALONE at %d %06d\n", date, time);
  TDataSet *DB = dbMk->GetInputDB("Geometry/ftt");
  if (!DB) { printf("  no Geometry/ftt -- nothing is served at this timestamp\n"); return; }

  const char *names[3];
  const int   want[3] = { 1, 4, 16 };
  names[0] = "stgcOnTpc"; names[1] = "stationOnStgc"; names[2] = "pentOnStation";

  int nServed = 0;
  for (int it = 0; it < 3; it++) {
    St_Survey *t = (St_Survey *) DB->Find(names[it]);
    if (!t) { printf("  %-16s NOT FOUND\n", names[it]); continue; }
    TDatime val[2];
    dbMk->GetValidity((TTable *)t, val);
    Survey_st *s = t->GetTable();
    int n = t->GetNRows(), nIdent = 0;
    for (int ir = 0; ir < n; ir++) {
      bool ident = (s[ir].r00 == 1 && s[ir].r11 == 1 && s[ir].r22 == 1 &&
                    s[ir].r01 == 0 && s[ir].r02 == 0 && s[ir].r10 == 0 &&
                    s[ir].r12 == 0 && s[ir].r20 == 0 && s[ir].r21 == 0 &&
                    s[ir].t0  == 0 && s[ir].t1  == 0 && s[ir].t2  == 0);
      if (ident) nIdent++;
    }
    printf("  %-16s %3d rows%s, %3d identity   begin %s   row0 t = (%+.5f, %+.5f, %+.5f)\n",
           names[it], n, n == want[it] ? "    " : " (!)", nIdent,
           val[0].AsSQLString(), s[0].t0, s[0].t1, s[0].t2);
    if (n == want[it]) nServed++;
  }

  if (nServed == 3)
    printf("\n  all three tables served with the right row counts.\n\n");
  else
    printf("\n  *** %d of 3 tables served correctly -- loadGeometry() chains all three,\n"
           "      so this timestamp gives a HALF-BUILT geometry ***\n\n", nServed);
}
