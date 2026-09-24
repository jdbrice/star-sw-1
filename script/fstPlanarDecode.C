// fstPlanarDecode.C -- does the planar decode's dphi match the hit's real angular
// offset from its wedge centre?  TrackFitter.h uses
//     dphi = kFstzDirct[electronicWedge] * (localPosition[1] - phi_half)
// with NO kFstzFilp[disk], while StFstHitMaker builds the global position with
// kFstzFilp[disk] * kFstzDirct[wedge].  If the disk flip is really missing, the two
// must be anti-correlated on the disk where kFstzFilp = -1 (the middle disk).
void fstPlanarDecode(const char* file, int nev=800){
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  StMuDstMaker* mk = new StMuDstMaker(0,0,"",file,"st:MuDst.root",1);
  const double pitch = 30.0/128.0;            // deg per phi strip
  const double phi_half = 0.5*128*pitch;      // = 15 deg
  const int kDirct[12] = {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1};
  double sxy[3], sxx[3], syy[3]; int nn[3];
  for (int i=0;i<3;i++){ sxy[i]=sxx[i]=syy[i]=0; nn[i]=0; }
  int nread=0;
  for (int iev=0; iev<nev; iev++){
    if (mk->Make()) break;
    StMuDst* d=mk->muDst(); if(!d) continue;
    StMuFstCollection* fc=d->muFstCollection(); if(!fc) continue;
    nread++;
    for (int ih=0; ih<(int)fc->numberOfHits(); ih++){
      StMuFstHit* fh=fc->getHit(ih); if(!fh) continue;
      int disk = (int)fh->getDisk() - 1;                 // 0,1,2
      int gw   = (int)fh->getWedge() - 1;                // 0-35
      int ew   = gw % 12;                                // electronic wedge 0-11
      double ps = fh->getMeanPhiStrip();
      if (disk<0||disk>2||ps<0) continue;
      TVector3 g = fh->xyz();
      double p = TMath::ATan2(g.Y(),g.X())*180./TMath::Pi(); if(p<0) p+=360.;
      double dphiStored = fmod(p,30.0) - 15.0;           // real offset from wedge centre
      double dphiPlanar = kDirct[ew] * (ps*pitch - phi_half);   // what TrackFitter computes
      sxy[disk]+=dphiStored*dphiPlanar; sxx[disk]+=dphiPlanar*dphiPlanar;
      syy[disk]+=dphiStored*dphiStored; nn[disk]++;
    }
  }
  printf(">>> events %d\n", nread);
  printf(">>> FST disk   hits    slope d(stored)/d(planar)    correlation\n");
  for (int k=0;k<3;k++){
    if(!nn[k]) continue;
    double slope = sxx[k]>0 ? sxy[k]/sxx[k] : 0;
    double corr  = (sxx[k]>0&&syy[k]>0) ? sxy[k]/sqrt(sxx[k]*syy[k]) : 0;
    printf(">>>    %d      %7d          %+7.3f              %+6.3f   %s\n",
           k, nn[k], slope, corr, (slope<0)?"<-- MIRRORED":"ok");
  }
}
