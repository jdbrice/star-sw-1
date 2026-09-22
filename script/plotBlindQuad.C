TH1F* psum(TFile* f, const char* base, const char* q){
  TH1F* s=0;
  for(int dd=0; dd<4; dd++){
    TH1F* h=(TH1F*)f->Get(Form("%s_disk%d_quad%s",base,dd,q));
    if(!h) continue;
    if(!s){ s=(TH1F*)h->Clone(Form("p_%s_%s",base,q)); s->SetDirectory(0);} else s->Add(h);
  }
  return s;
}
void norm(TH1F* h){ // normalise to the sideband 2-5 cm so shapes are comparable
  if(!h) return;
  double sb=h->Integral(h->GetXaxis()->FindBin(2.0),h->GetXaxis()->FindBin(5.0))
           +h->Integral(h->GetXaxis()->FindBin(-5.0),h->GetXaxis()->FindBin(-2.0));
  if(sb>0) h->Scale(1.0/sb);
}
void plotBlindQuad(const char* file, const char* out){
  gStyle->SetOptStat(0);
  TFile* f=TFile::Open(file);
  TCanvas* c=new TCanvas("c","c",1000,450); c->Divide(2,1);
  const char* qn[4]={"A","B","C","D"};
  const char* ql[4]={"A  x>0,y>0","B  x>0,y<0","C  x<0,y<0","D  x<0,y>0"};
  int col[4]={kRed+1,kOrange+7,kAzure+2,kBlue+2};
  c->cd(1); gPad->SetGridx();
  TH1F* dx[4];
  for(int i=0;i<4;i++){ dx[i]=psum(f,"hBlindDxAll_V",qn[i]); norm(dx[i]);
    dx[i]->SetLineColor(col[i]); dx[i]->SetLineWidth(2); }
  dx[3]->GetXaxis()->SetRangeUser(-5,5);
  dx[3]->SetTitle("V strips: dx = x_{hit} - x_{blind proj}, by quadrant;dx [cm];sideband-normalised");
  dx[3]->SetMaximum(dx[3]->GetMaximum()*1.35);
  dx[3]->Draw("hist");
  for(int j=0;j<3;j++) dx[j]->Draw("hist same");
  TLatex t; t.SetNDC(); t.SetTextSize(0.040);
  for(int k=0;k<4;k++){ t.SetTextColor(col[k]); t.DrawLatex(0.14,0.86-0.055*k, ql[k]); }
  t.SetTextColor(kBlack); t.SetTextSize(0.035);
  t.DrawLatex(0.14,0.62,"x>0 at -1.2 cm, x<0 at +1.2 cm");
  c->cd(2); gPad->SetGridx();
  TH1F* dy[4];
  for(int m=0;m<4;m++){ dy[m]=psum(f,"hBlindDyAll_H",qn[m]); norm(dy[m]);
    dy[m]->SetLineColor(col[m]); dy[m]->SetLineWidth(2); }
  dy[3]->GetXaxis()->SetRangeUser(-5,5);
  dy[3]->SetTitle("H strips: dy = y_{hit} - y_{blind proj}, by quadrant;dy [cm];sideband-normalised");
  dy[3]->SetMaximum(dy[3]->GetMaximum()*1.35);
  dy[3]->Draw("hist");
  for(int n=0;n<3;n++) dy[n]->Draw("hist same");
  for(int p=0;p<4;p++){ t.SetTextColor(col[p]); t.SetTextSize(0.040); t.DrawLatex(0.14,0.86-0.055*p, ql[p]); }
  t.SetTextColor(kBlack); t.SetTextSize(0.035);
  t.DrawLatex(0.14,0.62,"all four positive: common-mode");
  c->SaveAs(out);
}
