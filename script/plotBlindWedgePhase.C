// plotBlindWedgePhase.C -- overlay the FST->sTGC match peak in the four wedge-phase
// bins, for one configuration. Flat = the wedge convention is right; peak only in
// the first bin = the convention is globally mirrored.
// Usage: root4star -l -b -q 'plotBlindWedgePhase.C("a.root","g100","b.root","g110","out.png")'
TH1F* wpg(TFile* f, const char* base, const char* wp){
  TH1F* s = 0;
  for (int dd = 0; dd < 4; dd++){
    TH1F* h = (TH1F*)f->Get(Form("%s_disk%d_wp%s", base, dd, wp));
    if (!h) continue;
    if (!s){ s = (TH1F*)h->Clone(Form("g_%s_%s_%p", base, wp, (void*)f)); s->SetDirectory(0); }
    else s->Add(h);
  }
  return s;
}
void wpnorm(TH1F* h){
  if (!h) return;
  double sb = h->Integral(h->GetXaxis()->FindBin(2.0), h->GetXaxis()->FindBin(5.0))
            + h->Integral(h->GetXaxis()->FindBin(-5.0), h->GetXaxis()->FindBin(-2.0));
  if (sb > 0) h->Scale(1.0/sb);
}
void onePad(TFile* f, const char* base, const char* title, const char* note){
  const char* wpN[4] = {"0to4","4to8","8to12","12to15"};
  const char* wpL[4] = {"|#delta| 0-4#circ (centreline)","|#delta| 4-8#circ",
                        "|#delta| 8-12#circ","|#delta| 12-15#circ (edge)"};
  int col[4] = {kRed+1, kOrange+7, kAzure+2, kBlue+2};
  TH1F* h[4]; double mx = 0;
  for (int i = 0; i < 4; i++){
    h[i] = wpg(f, base, wpN[i]); wpnorm(h[i]);
    if (!h[i]) continue;
    h[i]->SetLineColor(col[i]); h[i]->SetLineWidth(2);
    h[i]->GetXaxis()->SetRangeUser(-5,5);
    if (h[i]->GetMaximum() > mx) mx = h[i]->GetMaximum();
  }
  if (!h[0]) return;
  h[0]->SetTitle(Form("%s;d [cm];sideband-normalised", title));
  h[0]->SetMaximum(mx*1.35); h[0]->Draw("hist");
  for (int j = 1; j < 4; j++) if (h[j]) h[j]->Draw("hist same");
  TLatex t; t.SetNDC(); t.SetTextSize(0.038);
  for (int k = 0; k < 4; k++){ t.SetTextColor(col[k]); t.DrawLatex(0.14, 0.86-0.052*k, wpL[k]); }
  t.SetTextColor(kBlack); t.SetTextSize(0.034); t.DrawLatex(0.14, 0.63, note);
}
void plotBlindWedgePhase(const char* fa, const char* ta, const char* fb, const char* tb, const char* out){
  gStyle->SetOptStat(0);
  TFile* A = TFile::Open(fa); TFile* B = TFile::Open(fb);
  TCanvas* c = new TCanvas("cwp","cwp",1000,450); c->Divide(2,1);
  c->cd(1); gPad->SetGridx(); onePad(A, "hBlindDxAll_V", Form("%s (mirror OFF): dx by wedge phase", ta), "");
  c->cd(2); gPad->SetGridx(); onePad(B, "hBlindDxAll_V", Form("%s (mirror ON): dx by wedge phase", tb), "positive control");
  c->SaveAs(out);
}
