// plotBlindQuadFit.C -- quadrant-split blind FST->sTGC residuals with the peak fitted.
//
// Quadrants come from the REAL sTGC footprints (StFttDb per-plane origins), not from
// the sign of x,y: all four local origins sit ~8.4-9.5 cm above the beampipe, so a
// y=0 split mislabels every hit with 0 < y < ~9 cm.
//
// Two traps this macro avoids, both hit while writing it:
//   - CINT does not propagate reference parameters in an interpreted macro, so results
//     go to globals rather than to double&/TF1&.
//   - the fit must run on the UNSCALED histogram; sideband-normalising first distorts
//     the bin errors (and an Integral-based entry guard then skips every fit).
//     So: fit first, then scale the histogram and the fitted function together.
double gQmu[4], gQsg[4], gQsig[4]; TF1* gQfn[4];
// disk = -1 sums all four sTGC disks; 0-3 selects one
TH1F* qfSum(TFile* f, const char* base, const char* q, int disk){
  TH1F* s=0;
  for (int d=0; d<4; d++){
    if (disk >= 0 && d != disk) continue;
    TH1F* h=(TH1F*)f->Get(Form("%s_disk%d_quad%s", base, d, q));
    if(!h) continue;
    if(!s){ s=(TH1F*)h->Clone(Form("qf_%s_%s_d%d", base, q, disk)); s->SetDirectory(0); }
    else s->Add(h);
  }
  return s;
}
void qfFit(TH1F* h, int k){
  gQmu[k]=0; gQsg[k]=0; gQsig[k]=0; gQfn[k]=0;
  if(!h || h->GetEntries()<200) return;
  int b0=h->GetXaxis()->FindBin(-2.9), b1=h->GetXaxis()->FindBin(2.9);
  double best=-1; int bb=b0;
  for(int ib=b0+2; ib<=b1-2; ib++){
    double v=h->GetBinContent(ib-1)+h->GetBinContent(ib)+h->GetBinContent(ib+1);
    if(v>best){ best=v; bb=ib; }
  }
  double mu0=h->GetXaxis()->GetBinCenter(bb), lo=mu0-1.2, hi=mu0+1.2;
  TF1* fn=new TF1(Form("fq_%s",h->GetName()),"gaus(0)+pol1(3)",lo,hi);
  double base=0.5*(h->GetBinContent(h->GetXaxis()->FindBin(lo))+h->GetBinContent(h->GetXaxis()->FindBin(hi)));
  fn->SetParameters(TMath::Max(1.0,h->GetBinContent(bb)-base),mu0,0.4,base,0);
  fn->SetParLimits(1,lo,hi); fn->SetParLimits(2,0.05,0.9);
  h->Fit(fn,"QNR");
  double A=fn->GetParameter(0), eA=fn->GetParError(0);
  gQmu[k]=fn->GetParameter(1); gQsg[k]=fabs(fn->GetParameter(2));
  gQsig[k]=(eA>0)?A/eA:0; gQfn[k]=fn;
}
void qfPanel(TFile* f, const char* base, const char* title, const char* note, int disk){
  const char* qn[4]={"A","B","C","D"};
  const char* ql[4]={"A  +x,+y (upper right)","B  +x,-y (lower right)",
                     "C  -x,-y (lower left)", "D  -x,+y (upper left)"};
  int col[4]={kRed+1,kOrange+7,kAzure+2,kBlue+2};
  TH1F* h[4]; double mx=0;
  for(int i=0;i<4;i++){
    h[i]=qfSum(f,base,qn[i],disk);
    if(!h[i]) continue;
    qfFit(h[i], i);                               // fit BEFORE scaling
    double sb=h[i]->Integral(h[i]->GetXaxis()->FindBin(2.0),h[i]->GetXaxis()->FindBin(5.0))
             +h[i]->Integral(h[i]->GetXaxis()->FindBin(-5.0),h[i]->GetXaxis()->FindBin(-2.0));
    if(sb>0){
      h[i]->Scale(1.0/sb);
      if(gQfn[i]){                                // scale the fit the same way
        gQfn[i]->SetParameter(0, gQfn[i]->GetParameter(0)/sb);
        gQfn[i]->SetParameter(3, gQfn[i]->GetParameter(3)/sb);
        gQfn[i]->SetParameter(4, gQfn[i]->GetParameter(4)/sb);
      }
    }
    h[i]->SetLineColor(col[i]); h[i]->SetLineWidth(2);
    h[i]->GetXaxis()->SetRangeUser(-5,5);
    if(h[i]->GetMaximum()>mx) mx=h[i]->GetMaximum();
  }
  int first=-1; for(int i=0;i<4;i++) if(h[i]){ first=i; break; }
  if(first<0) return;
  h[first]->SetTitle(Form("%s;d [cm];sideband-normalised", title));
  h[first]->SetMaximum(mx*1.55); h[first]->Draw("hist");
  for(int i=0;i<4;i++) if(h[i]&&i!=first) h[i]->Draw("hist same");
  for(int i=0;i<4;i++) if(gQfn[i]){ gQfn[i]->SetLineColor(col[i]); gQfn[i]->SetLineStyle(2);
                                    gQfn[i]->SetLineWidth(2); gQfn[i]->Draw("same"); }
  TLatex t; t.SetNDC(); t.SetTextSize(0.034);
  for(int i=0;i<4;i++){
    t.SetTextColor(col[i]);
    if(h[i]) t.DrawLatex(0.14,0.88-0.050*i,
        Form("%s  %+.2f cm (%.0f#sigma)", ql[i], gQmu[i], gQsig[i]));
  }
  t.SetTextColor(kBlack); t.SetTextSize(0.030); t.DrawLatex(0.14,0.66,note);
}
// disk = -1 (default) sums all four sTGC disks; 0-3 makes the single-disk version
void plotBlindQuadFit(const char* file, const char* out, int disk = -1){
  gStyle->SetOptStat(0);
  TFile* f=TFile::Open(file);
  TString tag = (disk<0) ? TString("all 4 disks") : TString(Form("sTGC disk %d", disk));
  TCanvas* c=new TCanvas(Form("cqf%d",disk),"cqf",1000,450); c->Divide(2,1);
  c->cd(1); gPad->SetGridx();
  qfPanel(f,"hBlindDxAll_V",Form("V strips: dx = x_{hit} - x_{proj}  (%s)",tag.Data()),
          "dashed = gaus+pol1 fit", disk);
  c->cd(2); gPad->SetGridx();
  qfPanel(f,"hBlindDyAll_H",Form("H strips: dy = y_{hit} - y_{proj}  (%s)",tag.Data()),
          "quadrants from the real StFttDb footprints", disk);
  c->SaveAs(out);
}
