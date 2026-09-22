// fitBlindPeak.C -- fit the FST-blind match peak WITHOUT assuming it sits at zero.
// Gaussian + linear background over |d| < 3 cm; reports the peak centre, width and
// amplitude significance, and a purity recomputed in a window centred on the fit.
// Usage: root4star -l -b -q 'fitBlindPeak.C("b_g000_P_blinddiag.root","g000_P")'
void fbpOne(TFile* f, int d, int io, double &sumReal, double &sumMatch,
            double &sumMean, double &sumSig, double &sumWid, int &n){
  const char* nA = (io==0) ? Form("hBlindDxAll_V_disk%d",d)     : Form("hBlindDyAll_H_disk%d",d);
  const char* nM = (io==0) ? Form("hBlindDxMatched_V_disk%d",d) : Form("hBlindDyMatched_H_disk%d",d);
  TH1F* hA=(TH1F*)f->Get(nA); TH1F* hM=(TH1F*)f->Get(nM);
  if(!hA||!hM) return;
  // locate the peak first: the match is NOT centred at zero (measured offset ~1 cm), and
  // the V distributions sit on a broad shoulder, so fit locally around the maximum with a
  // width limit rather than a single wide fit that lets the gaussian eat the shoulder.
  int b0 = hA->GetXaxis()->FindBin(-2.9), b1 = hA->GetXaxis()->FindBin(2.9);
  double best = -1; int bbest = b0;
  for (int ib = b0+2; ib <= b1-2; ib++){
    double v = hA->GetBinContent(ib-1) + hA->GetBinContent(ib) + hA->GetBinContent(ib+1);
    if (v > best){ best = v; bbest = ib; }
  }
  double mu0 = hA->GetXaxis()->GetBinCenter(bbest);
  double lo = mu0 - 1.2, hi = mu0 + 1.2;
  TF1* fn = new TF1(Form("fn%d%d",d,io),"gaus(0)+pol1(3)", lo, hi);
  double base = 0.5*(hA->GetBinContent(hA->GetXaxis()->FindBin(lo)) + hA->GetBinContent(hA->GetXaxis()->FindBin(hi)));
  fn->SetParameters(TMath::Max(1.0, hA->GetBinContent(bbest) - base), mu0, 0.4, base, 0);
  fn->SetParLimits(1, lo, hi); fn->SetParLimits(2, 0.05, 0.9);
  hA->Fit(fn,"QNR");
  double A=fn->GetParameter(0), mu=fn->GetParameter(1), sg=fabs(fn->GetParameter(2));
  double eA=fn->GetParError(0);
  double binw = hA->GetBinWidth(1);
  double nreal = A*sg*sqrt(2*TMath::Pi())/binw;          // integral of the gaussian, in entries
  printf(">>>   disk%d %s  centre %+6.3f cm  sigma %5.3f  amp %8.0f +- %6.0f (%5.1f sig)  peak entries %8.0f\n",
         d, io==0?"V":"H", mu, sg, A, eA, eA>0?A/eA:0, nreal);
  sumReal += nreal; sumMatch += hM->Integral();
  sumMean += mu; sumSig += (eA>0?A/eA:0); sumWid += sg; n++;
}
void fitBlindPeak(const char* file, const char* tag){
  TFile* f=TFile::Open(file);
  if(!f||f->IsZombie()){ printf(">>> %s MISSING\n", tag); return; }
  printf(">>> %s\n", tag);
  double sr=0, sm=0, smu=0, ssig=0, swid=0; int n=0;
  fbpOne(f,0,0,sr,sm,smu,ssig,swid,n); fbpOne(f,1,0,sr,sm,smu,ssig,swid,n);
  fbpOne(f,2,0,sr,sm,smu,ssig,swid,n); fbpOne(f,3,0,sr,sm,smu,ssig,swid,n);
  fbpOne(f,0,1,sr,sm,smu,ssig,swid,n); fbpOne(f,1,1,sr,sm,smu,ssig,swid,n);
  fbpOne(f,2,1,sr,sm,smu,ssig,swid,n); fbpOne(f,3,1,sr,sm,smu,ssig,swid,n);
  if(n) printf(">>> %-10s SUMMARY  purity %6.2f%%  mean centre %+6.3f cm  mean sigma %5.3f  mean sig %5.1f\n",
               tag, sm>0?100*sr/sm:0, smu/n, swid/n, ssig/n);
}
