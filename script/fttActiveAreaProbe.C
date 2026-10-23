// fttActiveAreaProbe.C -- where does the sensitive gas actually start, in
// StFttDb quadrant-local mm? If the strips are right, strip 0's left edge
// (14.35 mm for XY) should sit at the edge of the active region.
//   root4star -l -b -q 'fttActiveAreaProbe.C'
void loadLibsAA(){
  gSystem->Load("libStarClassLibrary.so"); gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StarMagField"); gSystem->Load("StMagF");
  gSystem->Load("StDetectorDbMaker"); gSystem->Load("StTpcDb");
  gSystem->Load("StDbBroker"); gSystem->Load("StDbUtilities");
  gSystem->Load("St_db_Maker"); gSystem->Load("libStDb_Tables.so");
  gSystem->Load("StFttDbMaker");
}
void fttActiveAreaProbe(const char* gfile="fGeom_dev2026_misal.root", int date=20211016){
  loadLibsAA();
  St_db_Maker *dbMk = new St_db_Maker("db","StarDb");
  dbMk->SetFlavor("ofl"); dbMk->SetDateTime(date,0); dbMk->Init(); dbMk->Make();
  StFttDb *db = new StFttDb("fttDb");
  TDataSet *DB = dbMk->GetInputDB("Geometry/stgc");
  db->loadGeometryFromDb( DB?(St_Survey*)DB->Find("stgcOnTpc"):0,
                          DB?(St_Survey*)DB->Find("stationOnStgc"):0,
                          DB?(St_Survey*)DB->Find("pentOnStation"):0 );
  TGeoManager *gm = TGeoManager::Import(gfile);
  if(!gm){printf("cannot import\n");return;}

  // walk into STFM_1 (quad A, station 1) and take every sensitive gas sub-volume
  double lxmin=1e9, lxmax=-1e9, lymin=1e9, lymax=-1e9;
  int nfound=0;
  TGeoIterator ia(gm->GetTopVolume());
  TGeoNode *na=0;
  bool inside=false;
  while((na=ia.Next())){
    TString vn = na->GetVolume()->GetName();
    if(!(vn.BeginsWith("STGP")||vn.BeginsWith("STGL")||vn.BeginsWith("STGS"))) continue;
    const TGeoHMatrix *ma = ia.GetCurrentMatrix();
    const Double_t *ta = ma->GetTranslation();
    if( fabs(ta[2]-312.342) > 2.0 ) continue;     // station 1 only
    TGeoBBox *bb = (TGeoBBox*) na->GetVolume()->GetShape();
    double hx=bb->GetDX(), hy=bb->GetDY();
    for(int ic=0; ic<4; ic++){
      double sgx = (ic&1)? +1:-1, sgy=(ic&2)? +1:-1;
      double loc[3]; loc[0]=sgx*hx; loc[1]=sgy*hy; loc[2]=0;
      double glo[3]; ma->LocalToMaster(loc, glo);
      double dbl[3];
      if(!db->globalToLocal((UChar_t)0,(UChar_t)0, glo, dbl)) continue;
      double mx=dbl[0]*10.0, my=dbl[1]*10.0;
      // Station 1 holds all FOUR pentagons at the same z, and globalToLocal for
      // quad A maps the other three to large negative local coordinates. Keep
      // only what lands inside quad A's own pentagon.
      if(mx < -60.0 || mx > 900.0 || my < -60.0 || my > 900.0) continue;
      if(mx<lxmin)lxmin=mx; if(mx>lxmax)lxmax=mx;
      if(my<lymin)lymin=my; if(my>lymax)lymax=my;
    }
    nfound++;
    if(nfound<=0) printf("  %-5s centre global (%8.3f,%8.3f,%8.3f) halfsize (%6.2f,%6.2f) cm\n",
                         vn.Data(), ta[0],ta[1],ta[2], hx, hy);
  }
  printf("\n  %d sensitive sub-volumes at station 1\n", nfound);
  printf("  active region in StFttDb quadrant-local mm:  x [%.2f, %.2f]   y [%.2f, %.2f]\n",
         lxmin,lxmax,lymin,lymax);
  printf("\n  strip 0 left edge, measured table : 14.35 mm   (centre 15.95 - pitch/2)\n");
  printf("  strip 0 left edge, analytic form  : -3.20 mm   (centre -1.60 - pitch/2)\n");
  printf("  last strip right edge, measured   : 548.75 mm\n");
}
