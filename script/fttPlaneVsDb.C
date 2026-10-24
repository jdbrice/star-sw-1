// fttPlaneVsDb.C -- do the FTT tracking PLANE and the FTT hit FRAME describe the
// same geometry?
//
// Planar mode projects an FTT hit as
//     diff = (hit global xyz) - plane->getO();  u = diff.Dot(plane->getU()) ...
// with the plane taken from the TGeo node /HALL_1/CAVE_1/STGM_1/STFM_<n>/ in the
// fGeom cache, and the hit global xyz built by StFttClusterPointMaker from
// StFttDb's offsets. If the cache carries the survey (bfc "misalign") both sides
// are the surveyed geometry and the projection is unbiased; if it does not, every
// FTT measurement is a surveyed hit on a nominal plane -- Xihe's NOTE B.
//
// This prints, for all 16 (station, quadrant):
//   |O_tgeo - O_reco|    plane origin vs the origin the RECO path puts the
//                        quadrant at -- getGloablOffset_ClusterPoint(), which is
//                        a translation plus axis REFLECTIONS, not AGML's rotation
//   |O_tgeo - O_mtx|     the same for quadToGlobal()/localToGlobal(), the proper
//                        survey matrices
//   dU, dV               the in-plane bias from the 17.55 mm strip offset, for
//                        each of the two views -- what the alignment fit absorbed
//
// The two StFttDb paths are NOT interchangeable: the reflection path is what
// StFttClusterPointMaker::MakeGlobalPoints() actually uses, so it is the one that
// has to agree with the plane.
//
//   root4star -l -b -q 'script/fttPlaneVsDb.C("fG.root",20240622,0)'
//
// CINT: fixed arrays, unique loop variable names, no local std::map (CLAUDE.md).
void loadLibsPVD(){
  gSystem->Load("libStarClassLibrary.so"); gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gSystem->Load("StarMagField"); gSystem->Load("StMagF");
  gSystem->Load("StDetectorDbMaker"); gSystem->Load("StTpcDb");
  gSystem->Load("StDbBroker"); gSystem->Load("StDbUtilities");
  gSystem->Load("St_db_Maker"); gSystem->Load("libStDb_Tables.so");
  gSystem->Load("StFttDbMaker");
}

void fttPlaneVsDb(const char* gfile = "fGeom.root", int date = 20240622, int time = 0)
{
  loadLibsPVD();

  St_db_Maker *dbMk = new St_db_Maker("db","StarDb");
  dbMk->SetFlavor("ofl"); dbMk->SetDateTime(date,time); dbMk->Init(); dbMk->Make();
  StFttDb *db = new StFttDb("fttDb");
  TDataSet *DB = dbMk->GetInputDB("Geometry/stgc");
  St_Survey *s0 = DB ? (St_Survey*)DB->Find("stgcOnTpc")     : 0;
  St_Survey *s1 = DB ? (St_Survey*)DB->Find("stationOnStgc") : 0;
  St_Survey *s2 = DB ? (St_Survey*)DB->Find("pentOnStation") : 0;
  db->loadGeometryFromDb(s0,s1,s2);
  printf("\nfttPlaneVsDb: geometry file [%s]  date %d  geometryFromDb=%d\n",
         gfile, date, (int)db->geometryFromDb());

  TGeoManager *gm = TGeoManager::Import(gfile);
  if(!gm){ printf("fttPlaneVsDb: cannot import %s\n", gfile); return; }

  // Collect the 16 STFM node matrices by node name, keyed STFM_1..STFM_16.
  // Node name is the key -- TGeoIterator::GetPath(TString&) kills the loop in CINT.
  double nO[16][3], nU[16][3], nV[16][3];
  int    nSeen[16];
  for(int iz=0; iz<16; iz++){ nSeen[iz]=0;
    for(int jz=0; jz<3; jz++){ nO[iz][jz]=0; nU[iz][jz]=0; nV[iz][jz]=0; } }

  TGeoIterator itp(gm->GetTopVolume());
  TGeoNode *ndp = 0;
  while((ndp = itp.Next())){
    TString nn = ndp->GetName();
    if(!nn.BeginsWith("STFM_")) continue;
    int idx = TString(nn(5, nn.Length()-5)).Atoi() - 1;   // STFM_1 -> 0
    if(idx < 0 || idx > 15) continue;
    const TGeoHMatrix *mp = itp.GetCurrentMatrix();
    const Double_t *tp = mp->GetTranslation();
    const Double_t *rp = mp->GetRotationMatrix();         // row major
    nO[idx][0]=tp[0]; nO[idx][1]=tp[1]; nO[idx][2]=tp[2];
    nU[idx][0]=rp[0]; nU[idx][1]=rp[3]; nU[idx][2]=rp[6]; // column 0
    nV[idx][0]=rp[1]; nV[idx][1]=rp[4]; nV[idx][2]=rp[7]; // column 1
    nSeen[idx]=1;
  }
  int nGot=0; for(int ia=0; ia<16; ia++) nGot += nSeen[ia];
  printf("fttPlaneVsDb: found %d of 16 STFM nodes\n\n", nGot);

  // pent -> StFttDb quadrant, as FwdGeomUtils/createAllFttPlanes measured it
  int quadOfPent[4]; quadOfPent[0]=0; quadOfPent[1]=3; quadOfPent[2]=2; quadOfPent[3]=1;
  const char *qn[4]; qn[0]="A"; qn[1]="B"; qn[2]="C"; qn[3]="D";

  const double kOff = 1.755;   // cm -- the strip offset, 17.55 mm

  printf(" idx st pent q |  plane O (cm)                        | |O_tgeo-O| (mm) |  bias of +17.55 mm, in-plane (mm)\n");
  printf("                |                                     |  reco     mtx   |   H view         V view\n");
  printf("                |                                     |                 |   dU      dV     dU      dV\n");
  printf(" ---------------+-------------------------------------+-----------------+-----------------------------\n");

  double maxDOr=0, maxDOm=0;
  for(int ib=0; ib<16; ib++){
    if(!nSeen[ib]) continue;
    int st   = ib/4;
    int pent = ib%4;
    int q    = quadOfPent[pent];

    // --- the RECO path: translation + reflections, mm, exactly as
    // StFttClusterPointMaker::MakeGlobalPoints() applies it
    float dx=0, sx=1, dy=0, sy=1, dz=0, sz=1;
    db->getGloablOffset_ClusterPoint( (UChar_t)st, (UChar_t)q, dx, sx, dy, sy, dz, sz );
    double Or[3]; Or[0]=dx/10.0; Or[1]=dy/10.0; Or[2]=dz;      // local (0,0,0) -> global, cm

    // --- the proper survey matrices
    double l0[3]; l0[0]=0; l0[1]=0; l0[2]=0;
    double Om[3];
    if(!db->localToGlobal((UChar_t)st,(UChar_t)q,l0,Om)){ printf(" idx %2d: localToGlobal failed\n", ib); continue; }

    double dOr = sqrt( (nO[ib][0]-Or[0])*(nO[ib][0]-Or[0])
                     + (nO[ib][1]-Or[1])*(nO[ib][1]-Or[1]) ) * 10.0;   // mm, x-y only
    double dOm = sqrt( (nO[ib][0]-Om[0])*(nO[ib][0]-Om[0])
                     + (nO[ib][1]-Om[1])*(nO[ib][1]-Om[1]) ) * 10.0;

    // --- the bias. The bug adds +17.55 mm to the quadrant-local PRECISION
    // coordinate, which is the cluster's x for a horizontal view and its y for a
    // vertical one; MakeGlobalPoints then reflects. So in global cm:
    double bh[3]; bh[0]=kOff*sx; bh[1]=0;        bh[2]=0;   // H view
    double bv[3]; bv[0]=0;       bv[1]=kOff*sy;  bv[2]=0;   // V view
    double hU = (bh[0]*nU[ib][0] + bh[1]*nU[ib][1] + bh[2]*nU[ib][2]) * 10.0;
    double hV = (bh[0]*nV[ib][0] + bh[1]*nV[ib][1] + bh[2]*nV[ib][2]) * 10.0;
    double vU = (bv[0]*nU[ib][0] + bv[1]*nU[ib][1] + bv[2]*nU[ib][2]) * 10.0;
    double vV = (bv[0]*nV[ib][0] + bv[1]*nV[ib][1] + bv[2]*nV[ib][2]) * 10.0;

    printf(" %3d %2d  %2d  %s | (%9.4f,%9.4f,%9.3f)        | %7.3f %7.3f | %6.2f %6.2f %6.2f %6.2f\n",
           ib, st, pent, qn[q], nO[ib][0], nO[ib][1], nO[ib][2], dOr, dOm, hU, hV, vU, vV);

    if(dOr > maxDOr) maxDOr = dOr;
    if(dOm > maxDOm) maxDOm = dOm;
  }

  printf("\n  worst |O_tgeo - O_reco| = %.3f mm   (reflection path, what reco uses)\n", maxDOr);
  printf("  worst |O_tgeo - O_mtx|  = %.3f mm   (quadToGlobal survey matrices)\n", maxDOm);
  printf("\n  O agreeing => the plane and the hit share a geometry, and the only\n");
  printf("  remaining FTT position error is the 17.55 mm strip offset, whose\n");
  printf("  in-plane bias is tabulated above.\n");
}
