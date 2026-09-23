// fstGapSign2.C -- gap sign without the same-wedge truncation bias.
// The inner partner may be in ANY wedge (only |dphi|<4 deg), so there is no asymmetric
// cut at the wedge edge. <dphi> is profiled against the OUTER hit's position in its
// wedge: the gap error is a STEP of 1 deg at delta_out = 0, not a smooth trend.
void fstGapSign2(const char* file, int nev=3000){
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  StMuDstMaker* mk = new StMuDstMaker(0,0,"",file,"st:MuDst.root",1);
  TProfile* pr = new TProfile("pr","<dphi> vs outer-hit position in wedge;delta_{outer} [deg];<dphi> [deg]",30,-15,15,-4,4);
  int nread=0; const int MAXH=400;
  for (int iev=0; iev<nev; iev++){
    if (mk->Make()) break;
    StMuDst* d=mk->muDst(); if(!d) continue;
    StMuFstCollection* fc=d->muFstCollection(); if(!fc) continue;
    nread++;
    int nh=(int)fc->numberOfHits(); if(nh>MAXH) nh=MAXH;
    double ph[MAXH], toz[MAXH]; int dk[MAXH], inr[MAXH]; int n=0;
    for (int ih=0; ih<nh; ih++){
      StMuFstHit* fh=fc->getHit(ih); if(!fh) continue;
      TVector3 g=fh->xyz();
      double p=TMath::ATan2(g.Y(),g.X())*180./TMath::Pi(); if(p<0) p+=360.;
      double r=sqrt(g.X()*g.X()+g.Y()*g.Y()), z=g.Z(); if(fabs(z)<1) continue;
      ph[n]=p; toz[n]=r/z; dk[n]=(int)fh->getDisk();
      inr[n]=(((int)fh->getSensor())==0)?1:0; n++;
    }
    for (int i=0;i<n;i++) for (int j=0;j<n;j++){
      if (i==j) continue;
      if (!inr[i] || inr[j]) continue;             // i inner, j outer
      if (dk[i]==dk[j]) continue;
      if (fabs(toz[i]-toz[j]) > 0.010) continue;
      double dp=ph[i]-ph[j];
      while(dp>180) dp-=360;  while(dp<-180) dp+=360;
      if (fabs(dp) > 4.0) continue;                // no same-wedge cut: no edge truncation
      pr->Fill(fmod(ph[j],30.0)-15.0, dp);
    }
  }
  printf(">>> events %d\n", nread);
  printf(">>> delta_out : ");  for (int b=1;b<=30;b++) printf("%5.1f", pr->GetBinCenter(b));   printf("\n");
  printf(">>> <dphi>    : ");  for (int b=1;b<=30;b++) printf("%+5.2f", pr->GetBinContent(b)); printf("\n");
  printf(">>> entries   : ");  for (int b=1;b<=30;b++) printf("%5.0f", pr->GetBinEntries(b));  printf("\n");
  // step across delta_out = 0, using only |delta| in 1..10 deg (clear of the gap and the edges)
  double sp=0,sm=0; int np=0,nm=0;
  for (int b=1;b<=30;b++){
    double c=pr->GetBinCenter(b), v=pr->GetBinContent(b), e=pr->GetBinEntries(b);
    if (e<20 || fabs(c)<1.0 || fabs(c)>10.0) continue;
    if (c>0){ sp+=v*e; np+=e; } else { sm+=v*e; nm+=e; }
  }
  if (np&&nm) printf(">>> STEP across the wedge centre, |delta| 1-10 deg: %+6.3f deg  (wrong gap sign predicts +1.0, correct predicts 0)\n",
                     sp/np - sm/nm);
}
