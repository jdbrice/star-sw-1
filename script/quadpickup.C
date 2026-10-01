// per-quadrant, per-station FTT pickup. Quadrant convention from
// StFwdResidualMaker::fillPlaneUsage: 0=North-top 1=North-bottom 2=South-top 3=South-bottom
// (North = x<0). Bins: 4+3*p+0 = V(x), +1 = H(y); bin 19 = every track.
void quadpickup(const char* fn="fttadd_Primary.root"){
  TFile* f=TFile::Open(fn); if(!f||f->IsZombie()){printf("cannot open\n");return;}
  const char* qn[4]={"Ntop","Nbot","Stop","Sbot"};
  const char* cn[2]={"Pos","Neg"};
  printf("\n=== FTT pickup per quadrant (Pos+Neg summed) ===\n");
  printf("  %-7s %-7s", "quad", "tracks");
  for(int ip=0;ip<4;ip++) printf("   s%d x/y      ",ip+1);
  printf("\n");
  for(int iq=0;iq<4;iq++){
    double ntr=0, ux[4]={0,0,0,0}, uy[4]={0,0,0,0};
    for(int ic=0;ic<2;ic++){
      TH1F* h=(TH1F*)f->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_%s",qn[iq],cn[ic]));
      if(!h) continue;
      ntr+=h->GetBinContent(20);
      for(int ip=0;ip<4;ip++){ ux[ip]+=h->GetBinContent(4+3*ip+1); uy[ip]+=h->GetBinContent(4+3*ip+2); }
    }
    if(ntr<=0){printf("  %-7s  no tracks\n",qn[iq]);continue;}
    printf("  %-7s %7.0f",qn[iq],ntr);
    for(int ip=0;ip<4;ip++) printf("   %.3f/%.3f",ux[ip]/ntr,uy[ip]/ntr);
    printf("\n");
  }
}
