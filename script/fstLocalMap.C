// fstLocalMap.C -- derive the true strip -> sensor-local-phi mapping from data.
// For each FST hit, transform its stored global position into each of the three
// FTUS local frames of its wedge; the correct sensor is the one whose shape
// (rmin..rmax, Phi1..Phi2) contains the hit. Then report, per hit-sensor id,
// which FTUS it really is and how local phi runs against meanPhiStrip.
const int kE2G[12]  = {2, 7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8};
const int kE2G2[12] = {7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8, 2};
TGeoHMatrix gM[3]; double gP1[3], gP2[3], gRmin[3], gRmax[3];
bool loadWedge(int disk, int ew){
  int planeIndex = disk + 4;
  const int* wmap = (planeIndex==5) ? kE2G2 : kE2G;
  int wedgeIndex = wmap[ew];
  for (int s=1; s<=3; s++){
    TString path = Form("/HALL_1/CAVE_1/FSTM_1/FSTD_%d/FSTW_%d/FTUS_%d", planeIndex, wedgeIndex, s);
    if (!gGeoManager->cd(path)) return false;
    gM[s-1] = *(gGeoManager->GetCurrentMatrix());
    TGeoShape* sh = gGeoManager->GetCurrentNode()->GetVolume()->GetShape();
    TGeoTubeSeg* ts = (TGeoTubeSeg*)sh;
    gP1[s-1]=ts->GetPhi1(); gP2[s-1]=ts->GetPhi2(); gRmin[s-1]=ts->GetRmin(); gRmax[s-1]=ts->GetRmax();
  }
  return true;
}
void fstLocalMap(const char* file, int nev=300){
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  TGeoManager::Import("/star/u/akio/fcstrk11/star-sw-fwd/fGeom.root");
  StMuDstMaker* mk = new StMuDstMaker(0,0,"",file,"st:MuDst.root",1);
  // per hit-sensor id: which FTUS matched, and slope/offset of localphi vs strip
  long  cnt[3][3];  double sx[3][3], sy[3][3], sxx[3][3], sxy[3][3];
  for(int a=0;a<3;a++) for(int b=0;b<3;b++){ cnt[a][b]=0; sx[a][b]=sy[a][b]=sxx[a][b]=sxy[a][b]=0; }
  int nread=0;
  for (int iev=0; iev<nev; iev++){
    if (mk->Make()) break;
    StMuDst* d=mk->muDst(); if(!d) continue;
    StMuFstCollection* fc=d->muFstCollection(); if(!fc) continue;
    nread++;
    for (int ih=0; ih<(int)fc->numberOfHits(); ih++){
      StMuFstHit* fh=fc->getHit(ih); if(!fh) continue;
      int disk=(int)fh->getDisk()-1, gw=(int)fh->getWedge()-1, sn=(int)fh->getSensor();
      double ps=fh->getMeanPhiStrip();
      if(disk<0||disk>2||sn<0||sn>2||ps<0) continue;
      if(!loadWedge(disk, gw%12)) continue;
      TVector3 g=fh->xyz();
      double gl[3]={g.X(),g.Y(),g.Z()}, lo[3];
      for (int s=0;s<3;s++){
        gM[s].MasterToLocal(gl,lo);
        double rl=sqrt(lo[0]*lo[0]+lo[1]*lo[1]);
        double pl=TMath::ATan2(lo[1],lo[0])*180./TMath::Pi();
        double p1=gP1[s],p2=gP2[s]; double plw=pl; while(plw<p1) plw+=360.; 
        if (rl>=gRmin[s]-0.3 && rl<=gRmax[s]+0.3 && plw>=p1-0.3 && plw<=p2+0.3){
          cnt[sn][s]++; sx[sn][s]+=ps; sy[sn][s]+=plw; sxx[sn][s]+=ps*ps; sxy[sn][s]+=ps*plw;
        }
      }
    }
  }
  printf(">>> events %d\n", nread);
  printf(">>> hit-sensor  FTUS  hits    localphi range        d(localphi)/d(strip)\n");
  for (int a=0;a<3;a++) for(int b=0;b<3;b++){
    if(cnt[a][b]<50) continue;
    double n=cnt[a][b];
    double slope=(n*sxy[a][b]-sx[a][b]*sy[a][b])/(n*sxx[a][b]-sx[a][b]*sx[a][b]);
    double icpt=(sy[a][b]-slope*sx[a][b])/n;
    printf(">>>     %d       FTUS_%d %6ld   shape[%7.2f,%7.2f]   %+8.5f deg/strip   at strip0: %+8.3f\n",
           a, b+1, cnt[a][b], gP1[b], gP2[b], slope, icpt);
  }
}
