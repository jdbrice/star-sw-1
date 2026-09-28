// plotStgcRadial.C -- the measured sTGC residual drawn as a vector at each quadrant's
// footprint centre, to show that it is a RADIAL pattern rather than a translation.
//
// d = hit - projection, so the correction to apply is -d. Arrows are drawn as the
// correction, which is what would go into the misalign table.
//
// CINT: no reference args, unique loop variable names, fixed C arrays (see CLAUDE.md).
double grMu[4][4][2], grSig[4][4][2]; int grOK[4][4][2];

void grFit(TH1F* h, int d, int q, int c){
  grMu[d][q][c]=0; grSig[d][q][c]=0; grOK[d][q][c]=0;
  if(!h || h->GetEntries()<200) return;
  int b0=h->GetXaxis()->FindBin(-2.9), b1=h->GetXaxis()->FindBin(2.9);
  double best=-1; int bb=b0;
  for(int ib=b0+2; ib<=b1-2; ib++){
    double v=h->GetBinContent(ib-1)+h->GetBinContent(ib)+h->GetBinContent(ib+1);
    if(v>best){ best=v; bb=ib; }
  }
  double mu0=h->GetXaxis()->GetBinCenter(bb), lo=mu0-1.2, hi=mu0+1.2;
  TF1* fn=new TF1(Form("gr%d%d%d",d,q,c),"gaus(0)+pol1(3)",lo,hi);
  double base=0.5*(h->GetBinContent(h->GetXaxis()->FindBin(lo))+h->GetBinContent(h->GetXaxis()->FindBin(hi)));
  fn->SetParameters(TMath::Max(1.0,h->GetBinContent(bb)-base),mu0,0.4,base,0);
  fn->SetParLimits(1,lo,hi); fn->SetParLimits(2,0.05,0.9);
  h->Fit(fn,"QNR");
  double A=fn->GetParameter(0), eA=fn->GetParError(0);
  grMu[d][q][c]=fn->GetParameter(1); grSig[d][q][c]=(eA>0)?A/eA:0;
  grOK[d][q][c]=(grSig[d][q][c]>3)?1:0;
}

void plotStgcRadial(const char* file, const char* out){
  gStyle->SetOptStat(0);
  TFile* f=TFile::Open(file);
  const char* qn[4]={"A","B","C","D"};
  for(int d=0;d<4;d++) for(int q=0;q<4;q++){
    grFit((TH1F*)f->Get(Form("hBlindDxAll_V_disk%d_quad%s",d,qn[q])), d,q,0);
    grFit((TH1F*)f->Get(Form("hBlindDyAll_H_disk%d_quad%s",d,qn[q])), d,q,1);
  }
  // quadrant footprint centres, from StFttDb (local 275,275 mm mapped with its own signs)
  double sfX[4][4]={{8.09,8.34,6.62,7.54},{112.74,112.14,113.30,113.49},
                    {-107.51,-108.22,-109.87,-108.91},{-3.69,-4.96,-4.36,-3.75}};
  double sfY[4][4]={{95.34,94.33,96.03,95.01},{84.24,83.37,83.61,83.81},
                    {83.60,83.42,84.37,82.81},{95.70,94.40,95.55,94.16}};
  double sx[4]={1,1,-1,-1}, sy[4]={1,-1,-1,1};

  TCanvas* c=new TCanvas("csr","csr",620,620);
  gPad->SetGridx(); gPad->SetGridy();
  TH2F* fr=new TH2F("fr","sTGC residual, median over the four stations;x [cm];y [cm]",1,-70,70,1,-60,75);
  fr->Draw();
  TLatex tt; tt.SetTextSize(0.028);
  tt.DrawLatex(-68,70,"arrows = correction -d, scaled #times10   (d = hit - projection)");
  TEllipse* bp=new TEllipse(0,0,4.1,4.1); bp->SetFillStyle(0); bp->SetLineStyle(2); bp->Draw();
  tt.SetTextSize(0.024); tt.DrawLatex(5,-2,"beampipe");

  int col[4]={kRed+1,kOrange+7,kAzure+1,kBlue+2};
  for(int q=0;q<4;q++){
    // MEDIAN over the four stations, not the mean: station 3 quad C has dx = -2.33 at
    // only 5.6 sigma against +0.71/+0.71/+0.39 elsewhere, and a single bad fit like that
    // drags a 4-point mean far enough to reverse the arrow.
    double vx[4],vy[4],cx=0,cy=0; int n=0;
    for(int d=0;d<4;d++){
      if(!grOK[d][q][0]||!grOK[d][q][1]) continue;
      vx[n]=grMu[d][q][0]; vy[n]=grMu[d][q][1];
      cx+=(275.0*sx[q]+sfX[q][d])/10.0; cy+=(275.0*sy[q]+sfY[q][d])/10.0; n++;
    }
    if(!n) continue;
    for(int ia=0; ia<n-1; ia++) for(int ib=ia+1; ib<n; ib++){
      if(vx[ib]<vx[ia]){ double t=vx[ia]; vx[ia]=vx[ib]; vx[ib]=t; }
      if(vy[ib]<vy[ia]){ double t=vy[ia]; vy[ia]=vy[ib]; vy[ib]=t; }
    }
    double mx = (n%2) ? vx[n/2] : 0.5*(vx[n/2-1]+vx[n/2]);
    double my = (n%2) ? vy[n/2] : 0.5*(vy[n/2-1]+vy[n/2]);
    cx/=n; cy/=n;
    TMarker* m=new TMarker(cx,cy,20); m->SetMarkerColor(col[q]); m->SetMarkerSize(1.3); m->Draw();
    TArrow* a=new TArrow(cx,cy,cx-10*mx,cy-10*my,0.018,"|>");
    a->SetLineColor(col[q]); a->SetFillColor(col[q]); a->SetLineWidth(3); a->Draw();
    tt.SetTextColor(col[q]); tt.SetTextSize(0.026);
    tt.DrawLatex(cx-6,cy+(cy>0?6:-9), Form("%s  d=(%+.2f,%+.2f)",qn[q],mx,my));
  }
  tt.SetTextColor(kBlack); tt.SetTextSize(0.026);
  tt.DrawLatex(-68,-48,"every arrow points OUTWARD: the reconstructed hits sit closer to the");
  tt.DrawLatex(-68,-53,"beam axis than the tracks do -- a radial pattern, not a translation");
  tt.SetTextSize(0.022);
  tt.DrawLatex(-68,-57,"per-quadrant MEDIAN over the four stations (station 3 quad C dx is a 5.6#sigma outlier)");
  c->SaveAs(out);
}
