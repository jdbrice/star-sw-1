// fttQuadConvention.C -- measure the AGML pentagon placement convention and compare it
// against StFttDb's (dx,dy,sx,sy) reflection convention, per station and quadrant.
//
// AGML places the four pentagons of a station by pure ROTATIONS about z (0/90/180/270).
// StFttDb maps local->global by REFLECTIONS:  global = (sx*u + dx, sy*v + dy).
// Rotations and reflections differ in handedness, so for two of the four quadrants the
// two conventions cannot agree -- this macro measures which two, and by what.
//
// For each STFM node we take the global image of the pentagon local origin and of the
// local +x / +y unit vectors.  Axis directions are origin-independent, so they compare
// directly with StFttDb's (sx,sy) even though the two frames use different origins
// (AGML: the shape corner; StFttDb: the pin hole).
//
// CINT: no std::map, no reference args, unique loop variable names (see CLAUDE.md).

const char* qname[4] = {"A","B","C","D"};   // StFttDb quad index 0..3
// StFttDb::getGloablOffset_ClusterPoint, inlined (consts from headers are not visible in CINT)
double sfX[4][4] = { {   8.09,    8.34,    6.62,    7.54},   // A
                     { 112.74,  112.14,  113.30,  113.49},   // B
                     {-107.51, -108.22, -109.87, -108.91},   // C
                     {  -3.69,   -4.96,   -4.36,   -3.75} }; // D
double sfY[4][4] = { {  95.34,   94.33,   96.03,   95.01},
                     {  84.24,   83.37,   83.61,   83.81},
                     {  83.60,   83.42,   84.37,   82.81},
                     {  95.70,   94.40,   95.55,   94.16} };
double sfSX[4] = { 1, 1,-1,-1};
double sfSY[4] = { 1,-1,-1, 1};

// which StFttDb quadrant does a global (x,y) corner pair sit in
int quadOf(double x, double y){
  if (x>0 && y>0) return 0;
  if (x>0 && y<0) return 1;
  if (x<0 && y<0) return 2;
  return 3;
}
const char* axname(double ax, double ay){
  if (ax> 0.9) return "+x";
  if (ax<-0.9) return "-x";
  if (ay> 0.9) return "+y";
  if (ay<-0.9) return "-y";
  return "??";
}

void fttQuadConvention(const char* gfile = "fGeom.root"){
  TFile* fg = TFile::Open(gfile);
  if(!fg){ printf("cannot open %s\n", gfile); return; }
  TGeoManager* gm = (TGeoManager*)fg->Get("dyson");
  if(!gm){ printf("no TGeoManager 'dyson' in %s\n", gfile); return; }
  gm->cd();

  printf("\n=== AGML pentagon placement, from %s ===\n", gfile);
  printf("%-26s %8s %8s %8s   %-4s %-4s  %-6s %5s\n",
         "path","Ox(cm)","Oy(cm)","z(cm)","u->","v->","quad","det");

  // collect per (station,quad)
  double gOx[4][4], gOy[4][4], gOz[4][4];
  int    gUx[4][4], gVx[4][4];        // encoded axis: 0=+x 1=-x 2=+y 3=-y
  int    gDet[4][4];
  int    gHave[4][4];
  for (int is=0; is<4; is++) for (int iq=0; iq<4; iq++) gHave[is][iq]=0;

  if (!gm->cd("/HALL_1/CAVE_1/STGM_1")) { printf("no /HALL_1/CAVE_1/STGM_1\n"); return; }
  TGeoNode* mother = gm->GetCurrentNode();
  for (int ip=0; ip<mother->GetNdaughters(); ip++){
      TGeoNode* pn = mother->GetDaughter(ip);
      TString pnm = pn->GetName();
      if (!pnm.BeginsWith("STFM")) continue;
      TString path = Form("/HALL_1/CAVE_1/STGM_1/%s", pn->GetName());
      if (!gm->cd(path)) { printf("cd failed: %s\n", path.Data()); continue; }
      TGeoHMatrix* hm = gm->GetCurrentMatrix();
      double lo[3] = {0,0,0}, go[3];
      hm->LocalToMaster(lo, go);
      double lu[3] = {1,0,0}, gu[3], lv[3] = {0,1,0}, gv[3];
      hm->LocalToMasterVect(lu, gu);
      hm->LocalToMasterVect(lv, gv);
      double lm[3] = {30,30,0}, gmid[3];
      hm->LocalToMaster(lm, gmid);
      int q  = quadOf(gmid[0], gmid[1]);
      double det = gu[0]*gv[1] - gu[1]*gv[0];
      printf("%-26s %8.3f %8.3f %8.3f   %-4s %-4s  %-6s %5.1f\n",
             path.Data(), go[0], go[1], go[2], axname(gu[0],gu[1]), axname(gv[0],gv[1]),
             qname[q], det);
      int st = -1;
      double zz = go[2];
      if      (zz < 320) st = 0;
      else if (zz < 338) st = 1;
      else if (zz < 356) st = 2;
      else               st = 3;
      gOx[st][q]=go[0]; gOy[st][q]=go[1]; gOz[st][q]=go[2];
      gUx[st][q] = (gu[0]>0.9)?0:((gu[0]<-0.9)?1:((gu[1]>0.9)?2:3));
      gVx[st][q] = (gv[0]>0.9)?0:((gv[0]<-0.9)?1:((gv[1]>0.9)?2:3));
      gDet[st][q]= (det>0)?1:-1;
      gHave[st][q]=1;
  }

  const char* axlbl[4] = {"+x","-x","+y","-y"};
  printf("\n=== convention comparison, per quadrant (station 1 shown; all 4 identical) ===\n");
  printf("%-6s | %-18s | %-18s | %s\n","quad","AGML (rotation)","StFttDb (reflection)","verdict");
  printf("-------+--------------------+--------------------+------------------------\n");
  for (int iq=0; iq<4; iq++){
    if (!gHave[0][iq]) { printf("%-6s | (not found)\n", qname[iq]); continue; }
    int sdbU = (sfSX[iq]>0)?0:1;          // StFttDb: local u -> sx * global x
    int sdbV = (sfSY[iq]>0)?2:3;          //          local v -> sy * global y
    int same = (gUx[0][iq]==sdbU && gVx[0][iq]==sdbV);
    printf("%-6s | u->%-3s  v->%-3s det%+d | u->%-3s  v->%-3s det%+d | %s\n",
           qname[iq], axlbl[gUx[0][iq]], axlbl[gVx[0][iq]], gDet[0][iq],
           axlbl[sdbU], axlbl[sdbV], (sfSX[iq]*sfSY[iq]>0)?1:-1,
           same ? "AGREE" : "DIFFER: local u<->v swap");
  }

  printf("\n=== per-station z and origin of each pentagon ===\n");
  for (int is=0; is<4; is++){
    for (int iq=0; iq<4; iq++){
      if(!gHave[is][iq]) continue;
      // StFttDb origin of the same quadrant, in cm
      double dbx = sfX[iq][is]/10.0, dby = sfY[iq][is]/10.0;
      printf("  station %d quad %s : AGML corner (%8.3f,%8.3f) z=%8.3f   StFttDb pin (%8.3f,%8.3f) z=%8.3f\n",
             is+1, qname[iq], gOx[is][iq], gOy[is][iq], gOz[is][iq], dbx, dby, 0.0);
    }
  }
}

// --------------------------------------------------------------------------
// Closure test for FwdGeomUtils::fttChamberIndex( plane, quadrant ):
// take the centre of each quadrant's StFttDb footprint and ask the geometry
// which STFM pentagon actually contains it.  The answer must be
// 4*plane + (4-quad)%4.
// --------------------------------------------------------------------------
void fttChamberClosure(const char* gfile = "fGeom.root"){
  TFile* fg = TFile::Open(gfile);
  TGeoManager* gm = (TGeoManager*)fg->Get("dyson");
  if(!gm){ printf("no geometry\n"); return; }
  gm->cd();
  int nbad = 0;
  printf("\n%-6s %-6s | %-24s | %-8s %-8s | %s\n",
         "plane","quad","StFttDb footprint centre","expect","found","z(cm)");
  for (int ip=0; ip<4; ip++){
    for (int iq=0; iq<4; iq++){
      // footprint centre: local (275,275) mm, mapped with StFttDb's own convention
      double gx = (275.0*sfSX[iq] + sfX[iq][ip])/10.0;
      double gy = (275.0*sfSY[iq] + sfY[iq][ip])/10.0;
      int expect = 4*ip + ((4-iq)%4);
      // probe at THIS station's z (where the hit really is), not at each
      // candidate's own z -- otherwise every station matches on |dz| and the
      // first one always wins.
      if (!gm->cd(Form("/HALL_1/CAVE_1/STGM_1/STFM_%d", 4*ip+1))) continue;
      double gz = gm->GetCurrentMatrix()->GetTranslation()[2];
      int found = -1; double fz = 0;
      for (int ic=0; ic<16; ic++){
        if (!gm->cd(Form("/HALL_1/CAVE_1/STGM_1/STFM_%d", ic+1))) continue;
        TGeoHMatrix* hm = gm->GetCurrentMatrix();
        double g[3] = {gx, gy, gz}, l[3];
        hm->MasterToLocal(g, l);
        if (l[0]>0 && l[1]>0 && l[0]<61 && l[1]<61 && fabs(l[2])<5.0){ found=ic; fz=hm->GetTranslation()[2]; break; }
      }
      if (found != expect) nbad++;
      printf("%-6d %-6s | (%8.2f,%8.2f)       | %-8d %-8d | %8.3f %s\n",
             ip, qname[iq], gx, gy, expect, found, fz, (found==expect)?"":"  <-- MISMATCH");
    }
  }
  printf("\n%s : %d mismatches out of 16\n", nbad?"FAIL":"PASS", nbad);
}
