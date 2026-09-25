// stgcQuadOffsets.C -- per disk/quadrant FST->sTGC residuals, decomposed into the
// three sTGC survey tables.
//
//   stgcOnTpc     (1 row)   global average        -- moves all of STGM
//   stationOnStgc (4 rows)  per-disk average minus global
//   pentOnStation (16 rows) per-quadrant residual, row = 4*station + k
//
// Residual measured is d = hit - projection, so the CORRECTION to apply is -d.
// Quadrants come from the real StFttDb footprints (see FwdTracker.h).
double gMu[4][4][2], gSig[4][4][2]; int gOK[4][4][2];  // [disk][quad][0=dx,1=dy]
// NOTE: disk1/quadA has ZERO V-strip entries in the data -- only one orientation is
// populated there, which is the online-QA-vs-offline-map H/V swap suspected since
// 2026-07-21. Cells with no usable fit are EXCLUDED from the averages rather than
// counted as zero, which would drag the station term.
TH1F* soSum(TFile* f, const char* base, const char* q, int disk){
  return (TH1F*)f->Get(Form("%s_disk%d_quad%s", base, disk, q));
}
void soFit(TH1F* h, int d, int q, int c){
  gMu[d][q][c]=0; gSig[d][q][c]=0; gOK[d][q][c]=0;
  if(!h || h->GetEntries()<200) return;
  int b0=h->GetXaxis()->FindBin(-2.9), b1=h->GetXaxis()->FindBin(2.9);
  double best=-1; int bb=b0;
  for(int ib=b0+2; ib<=b1-2; ib++){
    double v=h->GetBinContent(ib-1)+h->GetBinContent(ib)+h->GetBinContent(ib+1);
    if(v>best){ best=v; bb=ib; }
  }
  double mu0=h->GetXaxis()->GetBinCenter(bb), lo=mu0-1.2, hi=mu0+1.2;
  TF1* fn=new TF1(Form("so%d%d%d",d,q,c),"gaus(0)+pol1(3)",lo,hi);
  double base=0.5*(h->GetBinContent(h->GetXaxis()->FindBin(lo))+h->GetBinContent(h->GetXaxis()->FindBin(hi)));
  fn->SetParameters(TMath::Max(1.0,h->GetBinContent(bb)-base),mu0,0.4,base,0);
  fn->SetParLimits(1,lo,hi); fn->SetParLimits(2,0.05,0.9);
  h->Fit(fn,"QNR");
  double A=fn->GetParameter(0), eA=fn->GetParError(0);
  gMu[d][q][c]=fn->GetParameter(1); gSig[d][q][c]=(eA>0)?A/eA:0;
  gOK[d][q][c]=(gSig[d][q][c]>0)?1:0;
}
// writeRow: one Survey_st row in the style of the existing StarDb/Geometry/stgc files
// (same 1-byte memcpy for the comment as every sibling file; the readable text is the
// string literal itself, which is what a human reads in the .C)
void writeRow(FILE* fp, int id, double t0, double t1, const char* cmt){
  fprintf(fp,"\n    memset(&row,0,tableSet->GetRowSize());\n");
  fprintf(fp,"        row.Id   = %d;\n", id);
  fprintf(fp,"        row.r00  = 1.0;\n        row.r01  = 0.0;\n        row.r02  = 0.0;\n");
  fprintf(fp,"        row.r10  = 0.0;\n        row.r11  = 1.0;\n        row.r12  = 0.0;\n");
  fprintf(fp,"        row.r20  = 0.0;\n        row.r21  = 0.0;\n        row.r22  = 1.0;\n");
  fprintf(fp,"        row.t0   = %.5f;\n", t0);
  fprintf(fp,"        row.t1   = %.5f;\n", t1);
  fprintf(fp,"        row.t2   = 0.0;\n");
  fprintf(fp,"    memcpy(&row.comment,\"%s\\x00\",1);\n", cmt);
  fprintf(fp,"    tableSet->AddAt(&row);\n");
}
void openTable(FILE* fp, const char* name, int n, const char* hdr){
  fprintf(fp,"// %s\n", hdr);
  fprintf(fp,"TDataSet *CreateTable() {\n");
  fprintf(fp,"    if (!TClass::GetClass(\"St_Survey\")) return 0;\n");
  fprintf(fp,"Survey_st row;\n");
  fprintf(fp,"St_Survey *tableSet = new St_Survey(\"%s\",%d);\n", name, n);
}
void closeTable(FILE* fp){ fprintf(fp,"\nreturn (TDataSet *)tableSet;\n}\n"); }

void stgcQuadOffsets(const char* file, const char* outdir = 0){
  TFile* f=TFile::Open(file);
  const char* qn[4]={"A","B","C","D"};
  for(int d=0;d<4;d++) for(int q=0;q<4;q++){
    soFit(soSum(f,"hBlindDxAll_V",qn[q],d), d,q,0);
    soFit(soSum(f,"hBlindDyAll_H",qn[q],d), d,q,1);
  }
  printf(">>> MEASURED RESIDUALS  d = hit - projection  [cm]\n");
  printf(">>> disk quad      dx    (sig)      dy    (sig)\n");
  for(int d=0;d<4;d++){ for(int q=0;q<4;q++)
      printf(">>>   %d   %s     %s   %s\n", d, qn[q],
             gOK[d][q][0] ? Form("%+6.3f (%5.1f)", gMu[d][q][0], gSig[d][q][0]) : "   --  (no V hits)",
             gOK[d][q][1] ? Form("%+6.3f (%5.1f)", gMu[d][q][1], gSig[d][q][1]) : "   --  (no H hits)");
  }
  // global average
  double gx=0, gy=0; int nx=0, ny=0;
  for(int d=0;d<4;d++) for(int q=0;q<4;q++){
    if(gOK[d][q][0]){ gx+=gMu[d][q][0]; nx++; }
    if(gOK[d][q][1]){ gy+=gMu[d][q][1]; ny++; } }
  if(nx) gx/=nx; if(ny) gy/=ny;
  printf(">>>  (averages use %d/16 dx and %d/16 dy cells; the rest had no hits)\n", nx, ny);
  printf(">>>\n>>> stgcOnTpc  (1 row)  correction = -average\n");
  printf(">>>   t0 = %+8.5f   t1 = %+8.5f   cm   (average residual %+.3f, %+.3f)\n", -gx, -gy, gx, gy);
  // per-disk average minus global
  printf(">>>\n>>> stationOnStgc  (4 rows)  correction = -(disk average - global)\n");
  double dxD[4], dyD[4];
  for(int d=0;d<4;d++){
    double sx=0, sy=0; int mx=0, my=0;
    for(int q=0;q<4;q++){ if(gOK[d][q][0]){ sx+=gMu[d][q][0]; mx++; }
                          if(gOK[d][q][1]){ sy+=gMu[d][q][1]; my++; } }
    dxD[d]=(mx?sx/mx:gx)-gx; dyD[d]=(my?sy/my:gy)-gy;
    printf(">>>   row %d (station %d)  t0 = %+8.5f   t1 = %+8.5f\n", d, d, -dxD[d], -dyD[d]);
  }
  // per-quadrant residual
  printf(">>>\n>>> pentOnStation  (16 rows)  correction = -(quad - global - disk)\n");
  printf(">>>   row = 4*station + k   (k -> quadrant mapping still to be confirmed)\n");
  for(int d=0;d<4;d++) for(int q=0;q<4;q++){
    if(!gOK[d][q][0] || !gOK[d][q][1]){
      printf(">>>   station %d quad %s   -- incomplete (%s missing), leave identity\n",
             d, qn[q], gOK[d][q][0]?"dy":"dx"); continue; }
    double rx = gMu[d][q][0]-gx-dxD[d];
    double ry = gMu[d][q][1]-gy-dyD[d];
    printf(">>>   station %d quad %s   t0 = %+8.5f   t1 = %+8.5f\n", d, qn[q], -rx, -ry);
  }
  if(!outdir) return;

  // ---- emit the three DB tables ----------------------------------------
  // stationOnStgc is left IDENTITY (per-disk terms measured at <=0.08 cm), so
  // pentOnStation must absorb the full per-quadrant term relative to the global:
  //     pentOnStation[d][q] = -(m[d][q] - global)
  // Row order is station-major, and within a station k: 0->A 1->D 2->C 3->B,
  // measured from the geometry (script/pentMap.C), NOT A,B,C,D.
  const int k2q[4]={0,3,2,1};                       // k -> index into qn[] = A,B,C,D
  const char* k2qs[4]={"A,x>0,y>0","D,x<0,y>0","C,x<0,y<0","B,x>0,y<0"};
  char path[512]; FILE* fp;
  sprintf(path,"%s/stgcOnTpc.20211112.000001.C",outdir); fp=fopen(path,"w");
  openTable(fp,"stgcOnTpc",1,"sTGC global offset from FST->sTGC residuals (blind matching, ZF run 23063028). correction = -average residual");
  writeRow(fp,1,-gx,-gy,"whole STGC on TPC");
  closeTable(fp); fclose(fp); printf(">>> wrote %s\n",path);

  sprintf(path,"%s/stationOnStgc.20211112.000001.C",outdir); fp=fopen(path,"w");
  openTable(fp,"stationOnStgc",4,"per-station: IDENTITY. Measured per-disk common terms are <=0.08 cm in x and <=0.04 cm in y, an order of magnitude below the per-quadrant terms, so nothing is put here.");
  for(int d=0;d<4;d++){
    char c[128]; sprintf(c,"Station=%d(Disk=%d) identity",d+1,d);
    writeRow(fp,d+1,0.0,0.0,c);
  }
  closeTable(fp); fclose(fp); printf(">>> wrote %s\n",path);

  sprintf(path,"%s/pentOnStation.20211112.000001.C",outdir); fp=fopen(path,"w");
  openTable(fp,"pentOnStation",16,"per-quadrant offsets from FST->sTGC residuals. row = 4*station + k, k: 0=A 1=D 2=C 3=B (measured from the geometry, script/pentMap.C). correction = -(residual - global)");
  for(int d=0;d<4;d++) for(int k=0;k<4;k++){
    int q=k2q[k]; char c[160];
    double t0=0, t1=0;
    if(gOK[d][q][0] && gOK[d][q][1]){ t0=-(gMu[d][q][0]-gx); t1=-(gMu[d][q][1]-gy); }
    sprintf(c,"Station=%d(Disk=%d),Pent=%d(Quad=%s)%s", d+1,d,k+1,k2qs[k],
            (gOK[d][q][0]&&gOK[d][q][1])?"":" NO V HITS - identity");
    writeRow(fp, 4*d+k+1, t0, t1, c);
  }
  closeTable(fp); fclose(fp); printf(">>> wrote %s\n",path);
}
