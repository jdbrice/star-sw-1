// fstWedgePhase.C -- does fGeom.root (AGML) agree on front/back phase with
// kFstzFilp x kFstzDirct in StFstConsts.h?
// AGML: a wedge placed with the 180 deg x-rotation has zz = -1.
// Wedge index from the measured DB mapping: centre 75 deg -> wedge 0, 45 -> 1, 15 -> 2, ...
int zzOf[3][12];
void pwalk(TGeoNode* nd, TGeoHMatrix m, int depth){
  TGeoHMatrix mm = m; mm.Multiply(nd->GetMatrix());
  TString v = nd->GetVolume()->GetName();
  if (v.BeginsWith("FTUS")){
    const double* t = mm.GetTranslation();
    const double* r = mm.GetRotationMatrix();
    double a = TMath::ATan2(r[3], r[0])*180./TMath::Pi(); if(a<0) a+=360.;
    int d = (t[2] < 158.) ? 0 : ((t[2] < 172.) ? 1 : 2);
    int w = (int)TMath::Nint((75. - a)/30.); w = ((w % 12) + 12) % 12;
    zzOf[d][w] = (r[8] > 0) ? +1 : -1;
    return;
  }
  if (depth > 8) return;
  for (int i=0;i<nd->GetNdaughters();i++) pwalk(nd->GetDaughter(i), mm, depth+1);
}
void fstWedgePhase(){
  const int kFilp[3]  = {1,-1,1};
  const int kDirct[12]= {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1};
  for(int a=0;a<3;a++) for(int b=0;b<12;b++) zzOf[a][b]=0;
  TGeoManager::Import("/star/u/akio/fcstrk11/star-sw-fwd/fGeom.root");
  TGeoHMatrix id;
  pwalk(gGeoManager->GetTopNode(), id, 0);
  printf(">>> wedge  centre |  disk0 AGML  const |  disk1 AGML  const |  disk2 AGML  const\n");
  int nagree=0, ntot=0;
  for (int w=0; w<12; w++){
    double c = 75. - 30.*w; if(c<0) c+=360.;
    printf(">>>   %2d   %5.0f  |", w, c);
    for (int d=0; d<3; d++){
      int prod = kFilp[d]*kDirct[w];
      // observed mapping: AGML zz = -1  <->  product = +1
      int agree = ( (zzOf[d][w] == -1 && prod == +1) || (zzOf[d][w] == +1 && prod == -1) );
      printf("   %-5s  %+2d %s |", zzOf[d][w]<0?"BACK":"FRONT", prod, agree?"ok":"XX");
      nagree += agree; ntot++;
    }
    printf("\n");
  }
  printf(">>> agreement: %d / %d wedges\n", nagree, ntot);
}
