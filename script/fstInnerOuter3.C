// fstInnerOuter3.C -- select genuine track pairs by POLAR ANGLE (r/z), which is
// independent of phi, then ask where the phi correlation sits.
//   correct inner<->outer mapping => peak at dphi = 0
//   mirrored mapping              => peak at dphi = 2*delta_inner
// Both are histogrammed; whichever peaks wins. Random pairs are flat in both.
void fstInnerOuter3(const char* file, int nev=1500){
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  StMuDstMaker* mk = new StMuDstMaker(0,0,"",file,"st:MuDst.root",1);
  TH1F* hA = new TH1F("hA","inner-outer, r/z matched;dphi [deg];pairs",120,-30,30);
  TH1F* hB = new TH1F("hB","inner-outer, r/z matched;dphi - 2*delta [deg];pairs",120,-30,30);
  int nread=0; const int MAXH=400;
  for (int iev=0; iev<nev; iev++){
    if (mk->Make()) break;
    StMuDst* d=mk->muDst(); if(!d) continue;
    StMuFstCollection* fc=d->muFstCollection(); if(!fc) continue;
    nread++;
    int nh=(int)fc->numberOfHits(); if(nh>MAXH) nh=MAXH;
    double ph[MAXH], toz[MAXH]; int dk[MAXH], inr[MAXH], sec[MAXH]; int n=0;
    for (int ih=0; ih<nh; ih++){
      StMuFstHit* fh=fc->getHit(ih); if(!fh) continue;
      TVector3 g=fh->xyz();
      double p=TMath::ATan2(g.Y(),g.X())*180./TMath::Pi(); if(p<0) p+=360.;
      double r=sqrt(g.X()*g.X()+g.Y()*g.Y());
      double z=g.Z(); if (fabs(z)<1) continue;
      ph[n]=p; toz[n]=r/z; dk[n]=(int)fh->getDisk();
      inr[n]=((int)fh->getSensor()==0)?1:0; sec[n]=(int)(p/30.0); n++;
    }
    for (int i=0;i<n;i++) for (int j=0;j<n;j++){
      if (i==j) continue;
      if (!inr[i] || inr[j]) continue;        // i inner, j outer
      if (dk[i]==dk[j]) continue;
      if (sec[i]!=sec[j]) continue;
      if (fabs(toz[i]-toz[j]) > 0.010) continue;   // same polar angle: a real track
      double dp=ph[i]-ph[j];
      while(dp>180) dp-=360;  while(dp<-180) dp+=360;
      double dl=fmod(ph[i],30.0)-15.0;
      hA->Fill(dp);
      hB->Fill(dp-2.0*dl);
    }
  }
  printf(">>> events %d   r/z-matched inner-outer same-wedge pairs %.0f\n", nread, hA->GetEntries());
  double tot=TMath::Max(1.0,hA->Integral());
  printf(">>> hypothesis A (mapping correct, peak at dphi=0)      : |dphi|<2deg      = %5.1f%%   rms %6.2f\n",
         100.0*hA->Integral(hA->FindBin(-2.0),hA->FindBin(2.0))/tot, hA->GetRMS());
  printf(">>> hypothesis B (mirrored,        peak at dphi=2delta) : |dphi-2del|<2deg = %5.1f%%   rms %6.2f\n",
         100.0*hB->Integral(hB->FindBin(-2.0),hB->FindBin(2.0))/tot, hB->GetRMS());
  printf(">>> A profile (3deg bins, -30..30): ");
  for (int b=1;b<=20;b++) printf("%4.0f ", hA->Integral((b-1)*6+1,b*6));
  printf("\n>>> B profile (3deg bins, -30..30): ");
  for (int b=1;b<=20;b++) printf("%4.0f ", hB->Integral((b-1)*6+1,b*6));
  printf("\n");
}
