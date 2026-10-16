// fstGeometryCheck.C -- what geometry does the DATABASE ALONE give for a timestamp?
//
//   root4star -l -b -q 'fstGeometryCheck.C(20220322,221026)'   // a Run 22 real-data time
//   root4star -l -b -q 'fstGeometryCheck.C()'                  // same, the default
//
// St_db_Maker is created with MySQL as its ONLY source, which is what the BFC
// option NoCintDb does: no $STAR/StarDb, no ./StarDb. That matters because the
// release still ships Geometry/fst at 20211110.000001, LATER than every DB entry,
// so with the file sources in play the files win and the DB is never exercised.
// This is therefore the check that says what production will get once the release
// no longer carries those files.
//
// For each table it prints which entry was served (the validity range), how many
// rows are identity, and the first row's translation, which is where the survey
// shows itself: hssOnFst row 1 is t0 = -0.0775, t1 = +0.016 in the survey and all
// zeros in an identity entry.
//
// CINT note: Survey_st and St_db_Maker may only appear inside a function BODY --
// signatures are resolved when the file is loaded, before loadLibs makes them known.
// That is why this is one function rather than a helper per table.

void loadLibsFstChk()
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

void fstGeometryCheck(int date = 20220322, int time = 221026)
{
  loadLibsFstChk();
  // MySQL ONLY: no $STAR/StarDb, no ./StarDb -- the BFC NoCintDb case
  St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb");
  dbMk->SetFlavor("ofl");
  dbMk->SetDateTime(date, time);
  dbMk->Init();
  dbMk->Make();

  printf("\nfstGeometryCheck: Geometry/fst from the DATABASE ALONE at %d %06d\n", date, time);
  TDataSet *DB = dbMk->GetInputDB("Geometry/fst");
  if (!DB) { printf("  no Geometry/fst\n"); return; }

  const char *names[4];
  names[0] = "fstOnTpc"; names[1] = "hssOnFst";
  names[2] = "fstWedgeOnHss"; names[3] = "fstSensorOnWedge";
  for (int it = 0; it < 4; it++) {
    St_Survey *t = (St_Survey *) DB->Find(names[it]);
    if (!t) { printf("  %-18s NOT FOUND\n", names[it]); continue; }
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
    printf("  %-18s %3d rows, %3d identity   begin %s   row0 t = (%+.5f, %+.5f, %+.5f)\n",
           names[it], n, nIdent, val[0].AsSQLString(), s[0].t0, s[0].t1, s[0].t2);
  }
  printf("\n  hssOnFst row0 t = (-0.07750, +0.01600, 0) is THE SURVEY; all zeros is identity.\n\n");
}
