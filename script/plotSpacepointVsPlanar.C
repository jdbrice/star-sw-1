// plotSpacepointVsPlanar.C -- the blind FST->sTGC peak for the spacepoint path,
// the planar path as it was, and the planar path after the geometry-driven fixes.
TH1F* spvSum(TFile* f, const char* base){
  TH1F* s=0;
  for (int d=0; d<4; d++){
    TH1F* h=(TH1F*)f->Get(Form("%s_disk%d", base, d));
    if(!h) continue;
    if(!s){ s=(TH1F*)h->Clone(Form("s_%s_%p", base, (void*)f)); s->SetDirectory(0); }
    else s->Add(h);
  }
  if (s){ // normalise to its own 2-5 cm sideband so shapes compare
    double sb = s->Integral(s->GetXaxis()->FindBin(2.0), s->GetXaxis()->FindBin(5.0))
              + s->Integral(s->GetXaxis()->FindBin(-5.0), s->GetXaxis()->FindBin(-2.0));
    if (sb>0) s->Scale(1.0/sb);
  }
  return s;
}
void onePanel(TFile* a, TFile* b, TFile* c, const char* base, const char* title){
  TH1F* h1=spvSum(a,base); TH1F* h2=spvSum(b,base); TH1F* h3=spvSum(c,base);
  if(!h1) return;
  h1->SetLineColor(kBlack);   h1->SetLineWidth(2);
  if(h2){ h2->SetLineColor(kRed+1);   h2->SetLineWidth(2); }
  if(h3){ h3->SetLineColor(kAzure+2); h3->SetLineWidth(2); }
  double mx=h1->GetMaximum();
  if(h3&&h3->GetMaximum()>mx) mx=h3->GetMaximum();
  h1->GetXaxis()->SetRangeUser(-5,5);
  h1->SetTitle(Form("%s;d [cm];sideband-normalised", title));
  h1->SetMaximum(mx*1.45); h1->Draw("hist");
  if(h2) h2->Draw("hist same");
  if(h3) h3->Draw("hist same");
  TLatex t; t.SetNDC(); t.SetTextSize(0.037);
  t.SetTextColor(kBlack);   t.DrawLatex(0.15,0.86,"spacepoint  17.7%");
  t.SetTextColor(kRed+1);   t.DrawLatex(0.15,0.81,"planar, as it was  5.4%");
  t.SetTextColor(kAzure+2); t.DrawLatex(0.15,0.76,"planar, geometry-driven  17.3%");
}
void plotSpacepointVsPlanar(const char* fsp, const char* fold, const char* fnew, const char* out){
  gStyle->SetOptStat(0);
  TFile* a=TFile::Open(fsp); TFile* b=TFile::Open(fold); TFile* c=TFile::Open(fnew);
  TCanvas* cv=new TCanvas("csp","csp",1000,450); cv->Divide(2,1);
  cv->cd(1); gPad->SetGridx(); onePanel(a,b,c,"hBlindDyAll_H","H strips: dy (the clean channel)");
  cv->cd(2); gPad->SetGridx(); onePanel(a,b,c,"hBlindDxAll_V","V strips: dx");
  cv->SaveAs(out);
}
