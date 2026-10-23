// fttOriginProbe.C -- where is StFttDb's quadrant-local origin, relative to the
// pentagon the strips actually live on?
//
// If StFttDb's origin sits ON the pentagon corner, strips must start at positive
// local coordinates (scMapXY's 15.95 mm) and the analytic (strip-0.5)*pitch,
// which puts strip 0 at -1.6, is missing a real offset.
// Uses the MISALIGNED cache so the built geometry and StFttDb's matrices come
// from the same survey tables -- that is the only case where they share an origin.
//   root4star -l -b -q 'fttOriginProbe.C'
void loadLibsOP(){
  gSystem->Load("libStarClassLibrary.so"); gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StarMagField"); gSystem->Load("StMagF");
  gSystem->Load("StDetectorDbMaker"); gSystem->Load("StTpcDb");
  gSystem->Load("StDbBroker"); gSystem->Load("StDbUtilities");
  gSystem->Load("St_db_Maker"); gSystem->Load("libStDb_Tables.so");
  gSystem->Load("StFttDbMaker");
}
void fttOriginProbe(const char* gfile = "fGeom_dev2026_misal.root",
                   int date = 20211016, int time = 0){
  loadLibsOP();
  St_db_Maker *dbMk = new St_db_Maker("db","StarDb");
  dbMk->SetFlavor("ofl"); dbMk->SetDateTime(date,time); dbMk->Init(); dbMk->Make();
  StFttDb *db = new StFttDb("fttDb");
  TDataSet *DB = dbMk->GetInputDB("Geometry/stgc");
  St_Survey *s0 = DB ? (St_Survey*)DB->Find("stgcOnTpc") : 0;
  St_Survey *s1 = DB ? (St_Survey*)DB->Find("stationOnStgc") : 0;
  St_Survey *s2 = DB ? (St_Survey*)DB->Find("pentOnStation") : 0;
  db->loadGeometryFromDb(s0,s1,s2);
  printf("\nfttOriginProbe: geometryFromDb=%d\n", (int)db->geometryFromDb());

  TGeoManager *gm = TGeoManager::Import(gfile);
  if(!gm){printf("cannot import %s\n", gfile); return;}

  // quad A station 1 is STFM_1 (pent order A,D,C,B -> row 4*station+pent)
  const char *want[4]; want[0]="STFM_1"; want[1]="STFM_2"; want[2]="STFM_3"; want[3]="STFM_4";
  const char *qlab[4]; qlab[0]="A(pent0)"; qlab[1]="D(pent1)"; qlab[2]="C(pent2)"; qlab[3]="B(pent3)";
  int quadOf[4]; quadOf[0]=0; quadOf[1]=3; quadOf[2]=2; quadOf[3]=1;  // pent -> StFttDb quad

  printf("\n  node     quad        pentagon origin, global (cm)        -> StFttDb local (mm)\n");
  TGeoIterator ip(gm->GetTopVolume());
  TGeoNode *ndp=0;
  while((ndp=ip.Next())){
    TString nn = ndp->GetName();
    for(int iw=0; iw<4; iw++){
      if(nn != want[iw]) continue;
      const TGeoHMatrix *mp = ip.GetCurrentMatrix();
      const Double_t *tp = mp->GetTranslation();
      double g[3]; g[0]=tp[0]; g[1]=tp[1]; g[2]=tp[2];
      double l[3];
      if(!db->globalToLocal((UChar_t)0,(UChar_t)quadOf[iw], g, l)){ printf("  %s: globalToLocal failed\n", want[iw]); continue; }
      printf("  %-8s %-9s (%9.4f,%9.4f,%9.4f)   (%9.3f,%9.3f)\n",
             want[iw], qlab[iw], g[0],g[1],g[2], l[0]*10.0, l[1]*10.0);
    }
  }
  printf("\n  If StFttDb's origin were the pentagon corner these would be ~(0,0) mm.\n");
  printf("  scMapXY puts strip 0 at +15.95 mm; the analytic form puts it at -1.60 mm.\n");
}
