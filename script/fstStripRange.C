// fstStripRange.C -- is meanPhiStrip numbered per SENSOR (0-63) or per WEDGE ROW (0-127)
// for the two outer FST sensors? Decides whether the two outer sensors land in different
// halves of the wedge or on top of each other.
void fstStripRange(const char* file, int nev=400){
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  StMuDstMaker* mk = new StMuDstMaker(0,0,"",file,"st:MuDst.root",1);
  TH2F* h = new TH2F("h","meanPhiStrip vs sensor;sensor;meanPhiStrip",4,-0.5,3.5,140,-0.5,139.5);
  TH2F* hp= new TH2F("hp","phi within wedge vs sensor;sensor;phi mod 30 [deg]",4,-0.5,3.5,60,0,30);
  double smin[4], smax[4]; int scnt[4];
  for(int i=0;i<4;i++){ smin[i]=1e9; smax[i]=-1e9; scnt[i]=0; }
  int nread=0;
  for (int iev=0; iev<nev; iev++){
    if (mk->Make()) break;
    StMuDst* d = mk->muDst(); if(!d) continue;
    StMuFstCollection* fc = d->muFstCollection(); if(!fc) continue;
    nread++;
    for (int ih=0; ih<(int)fc->numberOfHits(); ih++){
      StMuFstHit* fh = fc->getHit(ih); if(!fh) continue;
      int    sn = (int)fh->getSensor();
      double ps = fh->getMeanPhiStrip();
      double rs = fh->getMeanRStrip();
      if (sn<0||sn>3) continue;
      h->Fill(sn, ps);
      if (ps<smin[sn]) smin[sn]=ps;
      if (ps>smax[sn]) smax[sn]=ps;
      scnt[sn]++;
      TVector3 g = fh->xyz();
      double ph = TMath::ATan2(g.Y(), g.X())*180./TMath::Pi(); if(ph<0) ph+=360.;
      hp->Fill(sn, fmod(ph,30.0));
    }
  }
  printf(">>> events read %d\n", nread);
  for (int s=0; s<4; s++){
    if (!scnt[s]) continue;
    printf(">>> sensor %d : %7d hits, meanPhiStrip range %6.1f .. %6.1f\n", s, scnt[s], smin[s], smax[s]);
  }
  // phi-within-wedge occupancy per sensor, coarse
  for (int s=0; s<4; s++){
    if (!scnt[s]) continue;
    printf(">>> sensor %d phi-in-wedge profile (3 deg bins): ", s);
    for (int b=0; b<10; b++){
      double n=0; for(int k=b*6+1;k<=(b+1)*6;k++) n+=hp->GetBinContent(s+1,k);
      printf("%5.1f%% ", scnt[s]?100.0*n/scnt[s]:0);
    }
    printf("\n");
  }
}
