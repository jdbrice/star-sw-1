// plotWedgePhaseSummary.C -- purity and peak significance vs |delta| from the wedge
// centreline, for mirror-off vs mirror-on. A global wedge-phase error would drive the
// correspondence to zero at the wedge edge; a correct convention keeps it flat.
TH1F* spSum(TFile* f, const char* base, const char* wp){
  TH1F* s = 0;
  for (int dd = 0; dd < 4; dd++){
    TH1F* h = (TH1F*)f->Get(Form("%s_disk%d_wp%s", base, dd, wp));
    if (!h) continue;
    if (!s){ s = (TH1F*)h->Clone(Form("q_%s_%s_%p", base, wp, (void*)f)); s->SetDirectory(0); }
    else s->Add(h);
  }
  return s;
}
double spFit(TH1F* h, double &sig){
  sig = 0; if (!h || h->Integral() < 200) return 0;
  int b0 = h->GetXaxis()->FindBin(-2.9), b1 = h->GetXaxis()->FindBin(2.9);
  double best = -1; int bb = b0;
  for (int ib = b0+2; ib <= b1-2; ib++){
    double v = h->GetBinContent(ib-1)+h->GetBinContent(ib)+h->GetBinContent(ib+1);
    if (v > best){ best = v; bb = ib; }
  }
  double mu0 = h->GetXaxis()->GetBinCenter(bb), lo = mu0-1.2, hi = mu0+1.2;
  TF1* fn = new TF1(Form("sf_%s", h->GetName()), "gaus(0)+pol1(3)", lo, hi);
  double base = 0.5*(h->GetBinContent(h->GetXaxis()->FindBin(lo))+h->GetBinContent(h->GetXaxis()->FindBin(hi)));
  fn->SetParameters(TMath::Max(1.0,h->GetBinContent(bb)-base), mu0, 0.4, base, 0);
  fn->SetParLimits(1, lo, hi); fn->SetParLimits(2, 0.05, 0.9);
  h->Fit(fn,"QNR");
  double A = fn->GetParameter(0), sg = fabs(fn->GetParameter(2)), eA = fn->GetParError(0);
  sig = (eA > 0) ? A/eA : 0;
  return A*sg*sqrt(2*TMath::Pi())/h->GetBinWidth(1);
}
void plotWedgePhaseSummary(const char* f1, const char* f2, const char* f3, const char* out){
  gStyle->SetOptStat(0);
  const char* wpN[4] = {"0to4","4to8","8to12","12to15"};
  double xc[4] = {2.0, 6.0, 10.0, 13.5};
  const char* files[3] = {f1, f2, f3};
  const char* labs[3]  = {"g100  mirror OFF", "g110  mirror ON", "g010  mirror ON, no gap fix"};
  int cols[3] = {kBlack, kRed+1, kOrange+7};
  TGraph* gp[3]; TGraph* gs[3];
  for (int c = 0; c < 3; c++){
    TFile* f = TFile::Open(files[c]);
    gp[c] = new TGraph(4); gs[c] = new TGraph(4);
    for (int w = 0; w < 4; w++){
      double sV, sH;
      TH1F* aV = spSum(f,"hBlindDxAll_V",wpN[w]);      TH1F* aH = spSum(f,"hBlindDyAll_H",wpN[w]);
      TH1F* mV = spSum(f,"hBlindDxMatched_V",wpN[w]);  TH1F* mH = spSum(f,"hBlindDyMatched_H",wpN[w]);
      double pk = spFit(aV,sV) + spFit(aH,sH);
      double mt = (mV?mV->Integral():0) + (mH?mH->Integral():0);
      gp[c]->SetPoint(w, xc[w], mt>0 ? 100.0*pk/mt : 0);
      gs[c]->SetPoint(w, xc[w], sH);
    }
    gp[c]->SetLineColor(cols[c]); gp[c]->SetMarkerColor(cols[c]);
    gp[c]->SetLineWidth(3); gp[c]->SetMarkerStyle(20); gp[c]->SetMarkerSize(1.3);
    gs[c]->SetLineColor(cols[c]); gs[c]->SetMarkerColor(cols[c]);
    gs[c]->SetLineWidth(3); gs[c]->SetMarkerStyle(20); gs[c]->SetMarkerSize(1.3);
  }
  TCanvas* c1 = new TCanvas("cs","cs",1000,450); c1->Divide(2,1);
  c1->cd(1); gPad->SetGridy(); gPad->SetGridx();
  TH2F* fr1 = new TH2F("fr1","FST#rightarrowsTGC match purity vs distance from wedge centreline;|#delta| from wedge centreline [deg];purity [%]",10,0,15.5,10,-1,23);
  fr1->Draw();
  for (int a = 0; a < 3; a++) gp[a]->Draw("PL same");
  TLatex t; t.SetNDC(); t.SetTextSize(0.038);
  for (int k = 0; k < 3; k++){ t.SetTextColor(cols[k]); t.DrawLatex(0.16, 0.36-0.055*k, labs[k]); }
  t.SetTextColor(kBlack); t.SetTextSize(0.032);
  t.DrawLatex(0.16, 0.86, "a wrong wedge phase would fall to zero here #rightarrow");
  c1->cd(2); gPad->SetGridy(); gPad->SetGridx(); gPad->SetLogy();
  TH2F* fr2 = new TH2F("fr2","H-strip peak significance vs distance from wedge centreline;|#delta| from wedge centreline [deg];peak significance [#sigma]",10,0,15.5,10,0.1,200);
  fr2->Draw();
  for (int b = 0; b < 3; b++) gs[b]->Draw("PL same");
  for (int m = 0; m < 3; m++){ t.SetTextColor(cols[m]); t.SetTextSize(0.038); t.DrawLatex(0.16, 0.36-0.055*m, labs[m]); }
  c1->SaveAs(out);
}
