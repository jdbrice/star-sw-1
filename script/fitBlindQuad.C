// fitBlindQuad.C -- is the dx double peak a north/south (x>0 vs x<0) effect?
// Sums the 4 disks per quadrant, fits gaus+pol1 around the local maximum.
TH1F* qsum(TFile* f, const char* base, const char* q){
  TH1F* s = 0;
  for (int dd=0; dd<4; dd++){
    TH1F* h = (TH1F*)f->Get(Form("%s_disk%d_quad%s", base, dd, q));
    if(!h) continue;
    if(!s){ s = (TH1F*)h->Clone(Form("sum_%s_%s", base, q)); s->SetDirectory(0); }
    else s->Add(h);
  }
  return s;
}
void qfit(TH1F* h, const char* label){
  if(!h){ printf(">>> %-22s missing\n", label); return; }
  int b0 = h->GetXaxis()->FindBin(-2.9), b1 = h->GetXaxis()->FindBin(2.9);
  double best=-1; int bb=b0;
  for (int ib=b0+2; ib<=b1-2; ib++){
    double v = h->GetBinContent(ib-1)+h->GetBinContent(ib)+h->GetBinContent(ib+1);
    if (v>best){ best=v; bb=ib; }
  }
  double mu0 = h->GetXaxis()->GetBinCenter(bb);
  double lo = mu0-1.2, hi = mu0+1.2;
  TF1* fn = new TF1(Form("f_%s",label), "gaus(0)+pol1(3)", lo, hi);
  double base = 0.5*(h->GetBinContent(h->GetXaxis()->FindBin(lo))+h->GetBinContent(h->GetXaxis()->FindBin(hi)));
  fn->SetParameters(TMath::Max(1.0,h->GetBinContent(bb)-base), mu0, 0.4, base, 0);
  fn->SetParLimits(1, lo, hi); fn->SetParLimits(2, 0.05, 0.9);
  h->Fit(fn,"QNR");
  double A=fn->GetParameter(0), mu=fn->GetParameter(1), sg=fabs(fn->GetParameter(2)), eA=fn->GetParError(0);
  double nreal = A*sg*sqrt(2*TMath::Pi())/h->GetBinWidth(1);
  printf(">>> %-22s centre %+6.3f  sigma %5.3f  amp %7.0f (%5.1f sig)  peak-entries %8.0f  hist-entries %9.0f\n",
         label, mu, sg, A, eA>0?A/eA:0, nreal, h->Integral());
}
void fitBlindQuad(const char* file, const char* tag){
  TFile* f = TFile::Open(file);
  if(!f||f->IsZombie()){ printf(">>> %s MISSING\n", tag); return; }
  printf(">>> ===== %s =====\n", tag);
  const char* qn[4] = {"A","B","C","D"};
  const char* qd[4] = {"A x>0,y>0","B x>0,y<0","C x<0,y<0","D x<0,y>0"};
  TH1F* dxq[4]; TH1F* dyq[4];
  for (int i=0;i<4;i++){ dxq[i]=qsum(f,"hBlindDxAll_V",qn[i]); dyq[i]=qsum(f,"hBlindDyAll_H",qn[i]); }
  for (int j=0;j<4;j++) qfit(dxq[j], Form("dx quad%s", qd[j]));
  for (int k=0;k<4;k++) qfit(dyq[k], Form("dy quad%s", qd[k]));
  // x>0 (A+B) vs x<0 (C+D)
  TH1F* dxp = (TH1F*)dxq[0]->Clone("dx_xpos"); dxp->SetDirectory(0); dxp->Add(dxq[1]);
  TH1F* dxn = (TH1F*)dxq[2]->Clone("dx_xneg"); dxn->SetDirectory(0); dxn->Add(dxq[3]);
  TH1F* dyp = (TH1F*)dyq[0]->Clone("dy_ypos"); dyp->SetDirectory(0); dyp->Add(dyq[3]); // y>0 = A+D
  TH1F* dyn = (TH1F*)dyq[1]->Clone("dy_yneg"); dyn->SetDirectory(0); dyn->Add(dyq[2]); // y<0 = B+C
  qfit(dxp, "dx  x>0 (A+B)"); qfit(dxn, "dx  x<0 (C+D)");
  qfit(dyp, "dy  y>0 (A+D)"); qfit(dyn, "dy  y<0 (B+C)");
  // double-gaussian on the all-quadrant dx sum, to quantify the two peaks
  TH1F* dxall = (TH1F*)dxp->Clone("dx_all"); dxall->SetDirectory(0); dxall->Add(dxn);
  TF1* g2 = new TF1("g2","gaus(0)+gaus(3)+pol1(6)", -3.0, 3.0);
  g2->SetParameters(dxall->GetMaximum()*0.5, -1.0, 0.4, dxall->GetMaximum()*0.5, 1.0, 0.4, 100, 0);
  g2->SetParLimits(1,-2.5,0.2); g2->SetParLimits(4,-0.2,2.5);
  g2->SetParLimits(2,0.1,0.9);  g2->SetParLimits(5,0.1,0.9);
  dxall->Fit(g2,"QNR");
  printf(">>> dx ALL two-gauss: peak1 %+6.3f (sig %5.3f, A %7.0f)   peak2 %+6.3f (sig %5.3f, A %7.0f)\n",
         g2->GetParameter(1), fabs(g2->GetParameter(2)), g2->GetParameter(0),
         g2->GetParameter(4), fabs(g2->GetParameter(5)), g2->GetParameter(3));
}
