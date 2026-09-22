// fitBlindWedgePhase.C -- FST->sTGC match quality vs |delta| from the wedge centreline.
//
// A global FST wedge-orientation phase error mirrors every hit about its wedge
// centreline, displacing it by 2*|delta|: nothing at the centre, 30 deg at the edge.
// So a WRONG convention keeps the correspondence only in the first bin and loses it
// towards the wedge edges; a RIGHT convention is flat in |delta|.
// Run it on g100 (question) and on g110 (positive control: hits really are mirrored
// there, so g110 MUST show the falling pattern if the binning works at all).
//
// Usage: root4star -l -b -q 'fitBlindWedgePhase.C("fwd_blind_diag.root","g100")'
TH1F* wpSum(TFile* f, const char* base, const char* wp){
  TH1F* s = 0;
  for (int dd = 0; dd < 4; dd++){
    TH1F* h = (TH1F*)f->Get(Form("%s_disk%d_wp%s", base, dd, wp));
    if (!h) continue;
    if (!s){ s = (TH1F*)h->Clone(Form("s_%s_%s", base, wp)); s->SetDirectory(0); }
    else s->Add(h);
  }
  return s;
}
// returns fitted gaussian area (entries); fills centre/sigma/significance
double wpFit(TH1F* h, double &mu, double &sg, double &sig){
  mu = 0; sg = 0; sig = 0;
  if (!h || h->Integral() < 200) return 0;
  int b0 = h->GetXaxis()->FindBin(-2.9), b1 = h->GetXaxis()->FindBin(2.9);
  double best = -1; int bb = b0;
  for (int ib = b0+2; ib <= b1-2; ib++){
    double v = h->GetBinContent(ib-1) + h->GetBinContent(ib) + h->GetBinContent(ib+1);
    if (v > best){ best = v; bb = ib; }
  }
  double mu0 = h->GetXaxis()->GetBinCenter(bb), lo = mu0-1.2, hi = mu0+1.2;
  TF1* fn = new TF1(Form("wf_%s", h->GetName()), "gaus(0)+pol1(3)", lo, hi);
  double base = 0.5*(h->GetBinContent(h->GetXaxis()->FindBin(lo)) + h->GetBinContent(h->GetXaxis()->FindBin(hi)));
  fn->SetParameters(TMath::Max(1.0, h->GetBinContent(bb)-base), mu0, 0.4, base, 0);
  fn->SetParLimits(1, lo, hi); fn->SetParLimits(2, 0.05, 0.9);
  h->Fit(fn, "QNR");
  double A = fn->GetParameter(0), eA = fn->GetParError(0);
  mu = fn->GetParameter(1); sg = fabs(fn->GetParameter(2)); sig = (eA>0) ? A/eA : 0;
  return A*sg*sqrt(2*TMath::Pi())/h->GetBinWidth(1);
}
void fitBlindWedgePhase(const char* file, const char* tag){
  TFile* f = TFile::Open(file);
  if (!f || f->IsZombie()){ printf(">>> %s MISSING\n", tag); return; }
  const char* wpN[4] = {"0to4","4to8","8to12","12to15"};
  const char* wpL[4] = {"0-4 deg (centreline)","4-8 deg","8-12 deg","12-15 deg (edge)"};
  printf(">>> ===== %s : match quality vs |delta| from the wedge centreline =====\n", tag);
  printf(">>> %-22s %10s %10s %8s %7s %6s %8s\n",
         "wedge phase","candidates","matched","peak","purity","sigma","centre");
  double p0 = 0;
  for (int w = 0; w < 4; w++){
    TH1F* aV = wpSum(f,"hBlindDxAll_V",wpN[w]);   TH1F* aH = wpSum(f,"hBlindDyAll_H",wpN[w]);
    TH1F* mV = wpSum(f,"hBlindDxMatched_V",wpN[w]); TH1F* mH = wpSum(f,"hBlindDyMatched_H",wpN[w]);
    double muV,sgV,sigV,muH,sgH,sigH;
    double pkV = wpFit(aV,muV,sgV,sigV), pkH = wpFit(aH,muH,sgH,sigH);
    double cand = (aV?aV->Integral():0) + (aH?aH->Integral():0);
    double match= (mV?mV->Integral():0) + (mH?mH->Integral():0);
    double pk = pkV + pkH;
    double pur = (match>0) ? 100.0*pk/match : 0;
    if (w==0) p0 = pur;
    printf(">>> %-22s %10.0f %10.0f %8.0f %6.2f%% %6.1f  V%+5.2f H%+5.2f  (V %4.1f sig, H %4.1f sig)%s\n",
           wpL[w], cand, match, pk, pur, 0.5*(sgV+sgH), muV, muH, sigV, sigH,
           (w>0 && p0>0) ? Form("  [%.2f x bin0]", pur/p0) : "");
  }
}
