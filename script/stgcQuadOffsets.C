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
int gEdge[4][4][2];                                    // peak search or fit ran into a limit
// Half-width of the peak SEARCH range, cm. 2.9 is the historical value and the right
// one while the tables are converged. Widen it when the residual has genuinely moved
// -- iteration 3 needed 5.9 for quadrant C's dx -- and the value used is printed with
// the results, because a residual only means something against the window it was
// found in. The FIT window stays peak+-1.2 either way.
double gSearchHalf = 2.9;
// NOTE: disk1/quadA has ZERO V-strip entries in the data -- only one orientation is
// populated there, which is the online-QA-vs-offline-map H/V swap suspected since
// 2026-07-21. Cells with no usable fit are EXCLUDED from the averages rather than
// counted as zero, which would drag the station term.
TH1F* soSum(TFile* f, const char* base, const char* q, int disk){
  return (TH1F*)f->Get(Form("%s_disk%d_quad%s", base, disk, q));
}
// The peak is searched only over +-2.9 cm and the fit's mean is then limited to
// peak+-1.2. That is plenty while the tables are converged and the residual sits near
// zero, but it CLIPS SILENTLY if a code change moves the residual: the 17.55 mm
// strip-frame fix (825d1fe21e) moved every quadrant's peak to ~1.755 cm in one step,
// which is already more than half the search range. gEdge flags a cell whose peak
// landed against the end of the search range, or whose fitted mean ran into its own
// parameter limit -- in either case the number is a boundary, not a measurement, and
// the window has to be widened before the delta is believed.
void soFit(TH1F* h, int d, int q, int c){
  gMu[d][q][c]=0; gSig[d][q][c]=0; gOK[d][q][c]=0; gEdge[d][q][c]=0;
  if(!h || h->GetEntries()<200) return;
  int b0=h->GetXaxis()->FindBin(-gSearchHalf), b1=h->GetXaxis()->FindBin(gSearchHalf);
  double best=-1; int bb=b0;
  for(int ib=b0+2; ib<=b1-2; ib++){
    double v=h->GetBinContent(ib-1)+h->GetBinContent(ib)+h->GetBinContent(ib+1);
    if(v>best){ best=v; bb=ib; }
  }
  if(bb<=b0+3 || bb>=b1-3) gEdge[d][q][c]=1;        // peak against the search range
  double mu0=h->GetXaxis()->GetBinCenter(bb), lo=mu0-1.2, hi=mu0+1.2;
  TF1* fn=new TF1(Form("so%d%d%d",d,q,c),"gaus(0)+pol1(3)",lo,hi);
  double base=0.5*(h->GetBinContent(h->GetXaxis()->FindBin(lo))+h->GetBinContent(h->GetXaxis()->FindBin(hi)));
  fn->SetParameters(TMath::Max(1.0,h->GetBinContent(bb)-base),mu0,0.4,base,0);
  fn->SetParLimits(1,lo,hi); fn->SetParLimits(2,0.05,0.9);
  h->Fit(fn,"QNR");
  double A=fn->GetParameter(0), eA=fn->GetParError(0);
  gMu[d][q][c]=fn->GetParameter(1); gSig[d][q][c]=(eA>0)?A/eA:0;
  gOK[d][q][c]=(gSig[d][q][c]>0)?1:0;
  // 1% of the half-width: Minuit parks a parameter exactly on its limit, so this is
  // a tolerance against rounding, not a near-miss criterion.
  if(gMu[d][q][c] <= lo+0.012 || gMu[d][q][c] >= hi-0.012) gEdge[d][q][c]=1;
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

// iteration = 1  : residuals were measured with StFttDb in HARDCODED mode, so the
//                  absolute form below is valid and the three tables are emitted.
// iteration >= 2 : residuals were measured with StFttDb in DB mode, i.e. hits already
//                  sat at (AGML placeholder + previous table). The absolute form would
//                  re-inject the hardcoded origin AND discard the previous table's
//                  alignment. Only the delta is written; add it to the previous table.
void stgcQuadOffsets(const char* file, const char* outdir = 0, int iteration = 1,
                     double searchHalf = 2.9){
  gSearchHalf = searchHalf;
  TFile* f=TFile::Open(file);
  const char* qn[4]={"A","B","C","D"};
  for(int d=0;d<4;d++) for(int q=0;q<4;q++){
    soFit(soSum(f,"hBlindDxAll_V",qn[q],d), d,q,0);
    soFit(soSum(f,"hBlindDyAll_H",qn[q],d), d,q,1);
  }
  printf(">>> MEASURED RESIDUALS  d = hit - projection  [cm]   (peak searched over +-%.1f cm)\n", gSearchHalf);
  printf(">>> disk quad      dx    (sig)      dy    (sig)\n");
  // de/qe, not d/q: CINT does not reliably reset a for-variable reused across two
  // loops in one function, and the loop can then run zero times with no error.
  int nEdge=0;
  for(int de=0;de<4;de++){ for(int qe=0;qe<4;qe++){
      const char* mark = "";
      if( (gOK[de][qe][0] && gEdge[de][qe][0]) || (gOK[de][qe][1] && gEdge[de][qe][1]) ){
        mark = "  <- AT WINDOW EDGE"; nEdge++;
      }
      printf(">>>   %d   %s     %s   %s%s\n", de, qn[qe],
             gOK[de][qe][0] ? Form("%+6.3f (%5.1f)", gMu[de][qe][0], gSig[de][qe][0]) : "   --  (no V hits)",
             gOK[de][qe][1] ? Form("%+6.3f (%5.1f)", gMu[de][qe][1], gSig[de][qe][1]) : "   --  (no H hits)",
             mark);
  } }
  if(nEdge)
    printf(">>>  WARNING: %d cell(s) have a peak or a fitted mean sitting on a limit of the\n"
           ">>>  +-%.1f cm search / +-1.2 cm fit window. Those numbers are the boundary, not a\n"
           ">>>  measurement -- rerun with a larger 4th argument before using their deltas.\n",
           nEdge, gSearchHalf);
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

  // ---- guard: the absolute formula below is iteration-1 only -------------------
  if( iteration != 1 ){
    printf("\n>>> iteration=%d: NOT writing absolute tables.\n", iteration);
    printf(">>> Hits already sat at (AGML placeholder + previous table), so the\n");
    printf(">>> absolute form would re-inject StFttDb's hardcoded origin and drop the\n");
    printf(">>> previous table's alignment. Correct update is  new = previous + delta,\n");
    printf(">>> with the delta being -(measured residual):\n\n");
    printf(">>> # disk quad   d_t0(cm)   d_t1(cm)\n");
    // Significance floor. A weak fit produces a large, meaningless centre -- station 3
    // quad C gave dx = -1.86 at 1.0 sigma against +0.71/+0.71/+0.39 elsewhere -- and
    // feeding that back as a delta injects it straight into the table. A cell whose fit
    // is not believable is left UNCHANGED, which is the safe direction: the previous
    // value stays and the next iteration can revisit it.
    const double kMinSig = 5.0;
    for(int dd=0; dd<4; dd++) for(int qq=0; qq<4; qq++){
      if(!gOK[dd][qq][0] || !gOK[dd][qq][1]){
        printf(">>> DELTA %d %s     (incomplete, leave unchanged)\n", dd, qn[qq]);
        continue;
      }
      // A cell sitting on the search/fit boundary gets NO parseable DELTA line at all,
      // the same shape as the incomplete branch above, so stgcTableAddDelta.py cannot
      // consume a boundary value by accident. Widen the window in soFit() and rerun --
      // unlike a low-sigma cell, this is not a weak measurement, it is a broken one.
      if(gEdge[dd][qq][0] || gEdge[dd][qq][1]){
        printf(">>> DELTA %d %s     (at window edge, leave unchanged -- rerun with a larger\n"
               ">>>                    4th argument than %.1f)\n", dd, qn[qq], gSearchHalf);
        continue;
      }
      // Per-coordinate, not per-cell: quad C has excellent dy fits (150+ sigma) and weak
      // dx ones, and dropping the whole cell would throw away a good dy correction to
      // avoid a bad dx one. A failing coordinate gets a zero delta, i.e. left unchanged.
      bool okx = gSig[dd][qq][0] >= kMinSig, oky = gSig[dd][qq][1] >= kMinSig;
      double dx = okx ? -gMu[dd][qq][0] : 0.0;
      double dy = oky ? -gMu[dd][qq][1] : 0.0;
      printf(">>> DELTA %d %s  %+9.5f  %+9.5f%s\n", dd, qn[qq], dx, dy,
             (okx&&oky) ? "" :
             Form("   <- %s held (dx %.1f, dy %.1f sigma)",
                  !okx ? (!oky ? "dx and dy" : "dx") : "dy",
                  gSig[dd][qq][0], gSig[dd][qq][1]));
    }
    printf("\n>>> Apply with script/stgcTableAddDelta.py (previous table + these deltas).\n");
    return;
  }


  // ---- emit the three DB tables ----------------------------------------
  // stationOnStgc is left IDENTITY (per-disk terms measured at <=0.08 cm), so
  // pentOnStation must absorb the full per-quadrant term relative to the global:
  //     pentOnStation[d][q] = -(m[d][q] - global)
  // Row order is station-major, and within a station k: 0->A 1->D 2->C 3->B,
  // measured from the geometry (script/pentMap.C), NOT A,B,C,D.
  const int k2q[4]={0,3,2,1};                       // k -> index into qn[] = A,B,C,D
  const char* k2qs[4]={"A,x>0,y>0","D,x<0,y>0","C,x<0,y<0","B,x>0,y<0"};
  // ITERATION 1 ONLY -- see the guard above and the note at the top of this function.
  // The AGML placement is a symmetric PLACEHOLDER: the quadrant local origins sit at
  // A(0,5.9) D(0,5.9) C(-6.5,5.9) B(+6.5,5.9), the same on every station. StFttDb knows
  // the real positions. So the table has to carry BOTH terms:
  //    pentOnStation = (StFttDb origin - AGML placeholder) + (our measured correction)
  // Then AGML and StFttDb describe the same true position, and StFttDb can later be
  // changed to compute dx,dy = placeholder + table instead of hardcoding.
  const double phx[4]={ 0.0,  0.0, -6.5, +6.5};     // by k: A, D, C, B
  const double phy[4]={ 5.9,  5.9,  5.9,  5.9};
  // StFttDb per-plane origins (cm), indexed [k][station]
  const double dbx[4][4]={{0.809,0.834,0.662,0.754},{-0.369,-0.496,-0.436,-0.375},
                          {-10.751,-10.822,-10.987,-10.891},{11.274,11.214,11.330,11.349}};
  const double dby[4][4]={{9.534,9.433,9.603,9.501},{9.570,9.440,9.555,9.416},
                          {8.360,8.342,8.437,8.281},{8.424,8.337,8.361,8.381}};
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
    int q=k2q[k]; char c[200];
    double geo0 = dbx[k][d]-phx[k];        // placeholder -> StFttDb
    double geo1 = dby[k][d]-phy[k];
    double cor0 = 0, cor1 = 0;             // measured alignment correction
    bool haveCor = (gOK[d][q][0] && gOK[d][q][1]);
    if(haveCor){ cor0=-(gMu[d][q][0]-gx); cor1=-(gMu[d][q][1]-gy); }
    double t0 = geo0+cor0, t1 = geo1+cor1;
    sprintf(c,"Station=%d(Disk=%d),Pent=%d(Quad=%s) geo(%+.3f,%+.3f)%s", d+1,d,k+1,k2qs[k],
            geo0, geo1, haveCor ? Form(" + align(%+.3f,%+.3f)",cor0,cor1)
                                : " ; NO V HITS - geometry term only");
    writeRow(fp, 4*d+k+1, t0, t1, c);
  }
  closeTable(fp); fclose(fp); printf(">>> wrote %s\n",path);
}
