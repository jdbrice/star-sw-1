// fttBridgeCheck.C -- agmlQuadLocalToDbLocal: the swap, the round trip, the refusal.
//   root4star -l -b -q 'fttBridgeCheck.C'
// CINT: unique loop variable names (see CLAUDE.md).
void loadLibsBr(){
  gSystem->Load("libStarClassLibrary.so"); gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StarMagField"); gSystem->Load("StMagF");
  gSystem->Load("StDetectorDbMaker"); gSystem->Load("StTpcDb");
  gSystem->Load("StDbBroker"); gSystem->Load("StDbUtilities");
  gSystem->Load("St_db_Maker"); gSystem->Load("libStDb_Tables.so");
  gSystem->Load("StFttDbMaker");
}
void fttBridgeCheck(int date = 20211012, int time = 0){
  loadLibsBr();
  St_db_Maker *dbMk = new St_db_Maker("db","StarDb");
  dbMk->SetFlavor("ofl"); dbMk->SetDateTime(date,time); dbMk->Init(); dbMk->Make();
  StFttDb *db = new StFttDb("fttDb");
  TDataSet *DB = dbMk->GetInputDB("Geometry/stgc");
  St_Survey *t0 = DB ? (St_Survey*)DB->Find("stgcOnTpc") : 0;
  St_Survey *t1 = DB ? (St_Survey*)DB->Find("stationOnStgc") : 0;
  St_Survey *t2 = DB ? (St_Survey*)DB->Find("pentOnStation") : 0;
  db->loadGeometryFromDb(t0,t1,t2);
  printf("\nfttBridgeCheck at %d: geometryFromDb = %d\n", date, (int)db->geometryFromDb());

  const char *qnm[4]; qnm[0]="A"; qnm[1]="B"; qnm[2]="C"; qnm[3]="D";
  double pin[3]; pin[0]=137.0; pin[1]=-42.5; pin[2]=3.0;
  double pout[3], pback[3];
  double worst = 0;
  printf("\n  quad  in (x,y)            -> db local (x,y)        round trip max|d|\n");
  for (int iqb = 0; iqb < 4; iqb++){
    if (!db->agmlQuadLocalToDbLocal(0,(UChar_t)iqb,pin,pout)) { printf("   %s  REFUSED\n",qnm[iqb]); continue; }
    db->dbLocalToAgmlQuadLocal(0,(UChar_t)iqb,pout,pback);
    double dmax = 0;
    for (int ic = 0; ic < 3; ic++){ double dd = fabs(pback[ic]-pin[ic]); if (dd>dmax) dmax=dd; }
    if (dmax > worst) worst = dmax;
    printf("   %s  (%8.3f,%8.3f) -> (%8.3f,%8.3f)    %.3g\n",
           qnm[iqb], pin[0], pin[1], pout[0], pout[1], dmax);
  }
  printf("\n  worst round-trip residual: %.3g  %s\n", worst, worst==0 ? "(exact)" : "");
  printf("  expected: A and C identity, B and D swapped\n");

  // and it must refuse once the DB geometry is dropped
  db->resetGeometryToHardcoded();
  bool got = db->agmlQuadLocalToDbLocal(0,1,pin,pout);
  printf("\n  after resetGeometryToHardcoded(): returns %s  -- %s\n",
         got ? "true" : "false", got ? "WRONG, it should refuse" : "correct, it refuses");
}
