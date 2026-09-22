// Is the wedge flip ALREADY in the reconstructed hits? For each wedge index, correlate the
// hit's phi within its wedge against the strip number. A per-wedge alternating SIGN means
// StFstHitMaker already reverses the strip direction (kFstzFilp * kFstzDirct).
void fstStripDirection(const char* file, int nev=300){
  gSystem->Load("libStarClassLibrary.so"); gSystem->Load("libStarRoot.so");
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries(); gSystem->Load("StEvent"); gSystem->Load("StMuDSTMaker");
  StChain* chain = new StChain("chain");
  StMuDstMaker* mk = new StMuDstMaker(0,0,"",file,"MuDst.root",1);
  chain->Init();
  double sx[36], sy[36], sxy[36], sxx[36]; int cn[36];
  for (int i=0;i<36;i++){ sx[i]=sy[i]=sxy[i]=sxx[i]=0; cn[i]=0; }
  for (int ie=0; ie<nev; ie++){
    if (chain->Make()!=0) break;
    StMuDst* d = mk->muDst(); if(!d) continue;
    StMuFstCollection* fst = d->muFstCollection(); if(!fst) continue;
    for (int ih=0; ih<fst->numberOfHits(); ih++){
      StMuFstHit* hh = fst->getHit(ih); if(!hh) continue;
      int wg = hh->getWedge()-1; if (wg<0||wg>35) continue;
      float ps = hh->getMeanPhiStrip(); if (ps < 0) continue;
      const TVector3 &g = hh->xyz();
      double p = atan2(g.Y(), g.X())*180.0/TMath::Pi(); while(p<0) p+=360;
      double ph = fmod(p, 30.0);                       // phase inside the wedge
      sx[wg]+=ps; sy[wg]+=ph; sxy[wg]+=ps*ph; sxx[wg]+=ps*ps; cn[wg]++;
    }
  }
  printf(">>> wedge(global) disk  nHits  slope d(phi)/d(strip)   direction\n");
  for (int w=0; w<36; w++){
    if (cn[w] < 20) continue;
    double den = cn[w]*sxx[w]-sx[w]*sx[w]; if (fabs(den)<1e-9) continue;
    double slope = (cn[w]*sxy[w]-sx[w]*sy[w])/den;
    printf(">>>   %2d          %d   %5d    %+8.4f            %s\n", w, w/12, cn[w], slope, slope>0?"+":"-");
  }
}
