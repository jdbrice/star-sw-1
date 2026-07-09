// plot_fcstrk_mom.C
// Standalone FCSTRK momentum resolution and charge-sign QA macro.
// Reads a picoMatch match.root file and produces:
//   fcsTrkMomPerEvt.png    — N BLCVtx / FCSTRK tracks per event
//   fcsTrkMomResolution.png — q/pT residual: BLCVtx vs FCSTRK, double-Gaussian fit
//   fcsTrkMomSignQ.png      — q/pT distribution (charge misID: fraction with q/pT>0 for e-)
//
// Track selection: 1 track per event, highest |pT| (primary electron for ele.e40 sample).
// Histograms used: hNTrkPerEvtBLCVtx, hNTrkPerEvtFCSTRK,
//                  hQPtSelBLCVtx, hQPtSelFCSTRK,
//                  hDqPtSelBLCVtx, hDqPtSelFCSTRK
// These are filled by picoMatch.C with the per-event highest-|pT| selection.
//
// Usage: root4star -b -q 'plot_fcstrk_mom.C("hist_match/all.match.root")'

void fitDoubleGaus(TH1F *h, int col,
                   double &sig1, double &sig2, double &impA, double &impN) {
    // Double Gaussian with shared center: A1*G(x,mu,s1) + A2*G(x,mu,s2), s1 < s2
    // Returns sig1 (narrow core) and sig2 (wide tail).
    // impA = fractional area in core, impN = core amplitude fraction.
    if(!h || h->GetEntries()<10){ sig1=sig2=impA=impN=0; return; }
    double mean = h->GetMean(), rms = h->GetRMS();
    TF1 *f = new TF1(Form("dg_%s",h->GetName()),
                     "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4])",
                     h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    f->SetParameters(h->GetMaximum(), mean, 0.5*rms,
                     0.3*h->GetMaximum(), 2.0*rms);
    double binW = h->GetBinWidth(1);
    f->SetParLimits(2, binW*0.5, rms);
    f->SetParLimits(4, rms*0.5, rms*10.0);
    f->SetLineColor(col); f->SetLineStyle(2); f->SetLineWidth(2);
    h->Fit(f, "RQN");
    sig1 = fabs(f->GetParameter(2));
    sig2 = fabs(f->GetParameter(4));
    if(sig1 > sig2){ double tmp=sig1; sig1=sig2; sig2=tmp; }
    double a1 = f->GetParameter(0)*sig1;
    double a2 = f->GetParameter(3)*sig2;
    impA = (a1+a2>0) ? a1/(a1+a2) : 0;
    impN = (f->GetParameter(0)+f->GetParameter(3)>0) ?
           f->GetParameter(0)/(f->GetParameter(0)+f->GetParameter(3)) : 0;
    f->Draw("same");
}

void plot_fcstrk_mom(const char* fname = "hist_match/all.match.root") {
    TFile *F = TFile::Open(fname);
    if(!F || F->IsZombie()){ printf("Cannot open %s\n",fname); return; }
    printf("Reading %s\n",fname);

    // Load histograms
    TH1F *hNGlob = (TH1F*)F->Get("nTrkPerEvtGlobal");
    TH1F *hNBLC  = (TH1F*)F->Get("nTrkPerEvtBLC");
    TH1F *hNBlc  = (TH1F*)F->Get("nTrkPerEvtBLCVtx");
    TH1F *hNFcs  = (TH1F*)F->Get("nTrkPerEvtFCSTRK");
    TH1F *hQGlob = (TH1F*)F->Get("qPtSelGlobal");
    TH1F *hQBLC  = (TH1F*)F->Get("qPtSelBLC");
    TH1F *hQBlc  = (TH1F*)F->Get("qPtSelBLCVtx");
    TH1F *hQFcs  = (TH1F*)F->Get("qPtSelFCSTRK");
    TH1F *hDBlc  = (TH1F*)F->Get("dqPtSelBLCVtx");
    TH1F *hDFcs  = (TH1F*)F->Get("dqPtSelFCSTRK");

    bool hasSel = hQBlc && hQFcs && hDBlc && hDFcs &&
                  hQBlc->GetEntries()>0 && hQFcs->GetEntries()>0;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // ── Page 1: N tracks per event ──────────────────────────────────────────
    if(hNBlc && hNFcs){
      TCanvas *c1 = new TCanvas("cN","",800,500);
      c1->Divide(2,1);
      for(int pad=0; pad<2; pad++){
        c1->cd(pad+1);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetLogy();
        TH1F *h = (pad==0) ? hNBlc : hNFcs;
        const char *name = (pad==0) ? "BLCVtx (type 4)" : "FCSTRK (type 5)";
        int col = (pad==0) ? kBlue : kRed;
        h->SetLineColor(col); h->SetLineWidth(2); h->SetFillColor(col-9);
        h->GetXaxis()->SetTitle("N tracks / event");
        h->GetYaxis()->SetTitle("Events");
        h->Draw("hist");
        double n0 = h->GetBinContent(1);
        double ntot = h->GetEntries();
        double frac0 = (ntot>0) ? 100.*n0/ntot : 0;
        double frac2p = 0;
        for(int ib=3; ib<=h->GetNbinsX(); ib++) frac2p += h->GetBinContent(ib);
        frac2p = (ntot>0) ? 100.*frac2p/ntot : 0;
        TLatex *lt = new TLatex(); lt->SetNDC(); lt->SetTextSize(0.047); lt->SetTextColor(col);
        lt->DrawLatex(0.55,0.83,name);
        lt->SetTextColor(kBlack); lt->SetTextSize(0.040);
        lt->DrawLatex(0.55,0.74,Form("N events: %.0f",ntot));
        lt->DrawLatex(0.55,0.66,Form("0 tracks: %.1f%%",frac0));
        lt->DrawLatex(0.55,0.58,Form("#geq2 tracks: %.1f%%",frac2p));
        printf("[%s] N events=%.0f  0 tracks=%.1f%%  >=2 tracks=%.1f%%\n",
               name, ntot, frac0, frac2p);
      }
      c1->SaveAs("fcsTrkMomPerEvt.png");
      printf("Saved fcsTrkMomPerEvt.png\n");
    }

    if(!hasSel){ printf("No per-event selection histograms found; exiting.\n"); return; }

    // ── Page 2: q/pT residual — Global / BLC / BLCVtx / FCSTRK ─────────────
    {
      const char *DNAME[4]  = {"dqPtSelGlobal","dqPtSelBLC","dqPtSelBLCVtx","dqPtSelFCSTRK"};
      const char *DLABEL[4] = {"Global","BLC","BLCVtx","FCSTRK"};
      int DCOL[4]; DCOL[0]=kBlack; DCOL[1]=kGreen+2; DCOL[2]=kBlue; DCOL[3]=kRed;

      TH1F *hD[4]; TH1F *hN[4];
      for(int i=0; i<4; i++){ hD[i]=0; hN[i]=0; }
      for(int i=0; i<4; i++){
        hD[i] = (TH1F*)F->Get(DNAME[i]);
        if(!hD[i] || hD[i]->GetEntries()<10) continue;
        hN[i] = (TH1F*)hD[i]->Clone(Form("hDn%d",i));
        if(hN[i]->Integral()>0) hN[i]->Scale(1./hN[i]->Integral());
        hN[i]->GetXaxis()->SetRangeUser(-2.0, 2.0);
        hN[i]->SetLineColor(DCOL[i]); hN[i]->SetLineWidth(2);
      }

      TCanvas *c2 = new TCanvas("cDqPt","",800,600);
      gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogy(1);

      // Find max y over all types for range
      double ymx = 0;
      for(int i=0; i<4; i++) if(hN[i] && hN[i]->GetMaximum()>ymx) ymx=hN[i]->GetMaximum();

      // Draw frame from first available
      TH1F *hFirst = 0;
      for(int i=0; i<4; i++) if(hN[i]){ hFirst=hN[i]; break; }
      if(hFirst){
        hFirst->SetMinimum(5e-5);
        hFirst->SetMaximum(ymx*5.0);
        hFirst->GetXaxis()->SetTitle("q/p_{T}^{reco} - q/p_{T}^{true}  [(GeV/c)^{-1}]");
        hFirst->GetYaxis()->SetTitle("Norm. entries / event");
        hFirst->Draw("hist");
        for(int i=0; i<4; i++) if(hN[i] && hN[i]!=hFirst) hN[i]->Draw("hist same");
      }

      // Fit each and draw colored text label
      double s1[4]; double s2[4]; double aFr[4]; double nFr[4];
      for(int i=0; i<4; i++){ s1[i]=0; s2[i]=0; aFr[i]=0; nFr[i]=0; }
      for(int i=0; i<4; i++) if(hN[i]) fitDoubleGaus(hN[i],DCOL[i],s1[i],s2[i],aFr[i],nFr[i]);

      TLatex *tl = new TLatex();
      tl->SetNDC(); tl->SetTextSize(0.033);
      double yleg = 0.38;
      for(int i=0; i<4; i++){
        if(!hN[i]) continue;
        tl->SetTextColor(DCOL[i]);
        tl->DrawLatex(0.30, yleg,
          Form("%s   #sigma_{core}=%.3f   #sigma_{tail}=%.3f", DLABEL[i], s1[i], s2[i]));
        yleg -= 0.058;
      }

      c2->SaveAs("fcsTrkMomResolution.png");
      printf("Saved fcsTrkMomResolution.png\n");
      printf("Resolution (double-Gaussian core, 1 track/event highest |pT|):\n");
      for(int i=0; i<4; i++)
        if(hD[i] && hD[i]->GetEntries()>0)
          printf("  %s: sigma_core=%.4f  sigma_tail=%.4f\n", DLABEL[i], s1[i], s2[i]);
    }

    // ── Page 3: q/pT absolute — charge sign misID (all 4 track types) ─────────
    // For e- (negative charge), correct sign = q/pT < 0.
    // Fraction with q/pT > 0 is the charge misidentification rate.
    // Zoom to [-1.5, 1.5]; misID fraction computed from full range.
    {
      const char *QLABEL[4] = {"Global","BLC","BLCVtx","FCSTRK"};
      int QCOL[4];  QCOL[0]=kBlack;  QCOL[1]=kGreen+2; QCOL[2]=kBlue;   QCOL[3]=kRed;
      int QFILL[4]; QFILL[0]=kGray;  QFILL[1]=kGreen-9; QFILL[2]=kBlue-9; QFILL[3]=kRed-9;
      TH1F *hQSrc[4]; hQSrc[0]=hQGlob; hQSrc[1]=hQBLC; hQSrc[2]=hQBlc; hQSrc[3]=hQFcs;
      TH1F *hNSrc[4]; hNSrc[0]=hNGlob; hNSrc[1]=hNBLC; hNSrc[2]=hNBlc; hNSrc[3]=hNFcs;

      // Pre-compute misID rates
      double qTot[4]; double qMis[4]; double qRate[4];
      for(int i=0; i<4; i++){ qTot[i]=0; qMis[i]=0; qRate[i]=0; }
      for(int i=0; i<4; i++){
        if(!hQSrc[i] || hQSrc[i]->GetEntries()<10) continue;
        int b0i = hQSrc[i]->FindBin(0.0);
        qTot[i] = hQSrc[i]->GetEntries();
        for(int ib=b0i+1; ib<=hQSrc[i]->GetNbinsX(); ib++) qMis[i] += hQSrc[i]->GetBinContent(ib);
        qRate[i] = (qTot[i]>0) ? 100.*qMis[i]/qTot[i] : 0;
      }

      TCanvas *c3 = new TCanvas("cSignQ","",800,600);
      gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogy();

      // Normalize to events with >=1 track of each type
      TH1F *hQN[4];
      for(int i=0; i<4; i++) hQN[i]=0;
      for(int i=0; i<4; i++){
        if(!hQSrc[i] || hQSrc[i]->GetEntries()<10) continue;
        hQN[i] = (TH1F*)hQSrc[i]->Clone(Form("hQN%d",i));
        double n1 = hNSrc[i] ? (hNSrc[i]->GetEntries() - hNSrc[i]->GetBinContent(1))
                              : hQSrc[i]->GetEntries();
        if(n1>0) hQN[i]->Scale(1./n1);
        hQN[i]->GetXaxis()->SetRangeUser(-1.5, 1.5);
        hQN[i]->SetLineColor(QCOL[i]); hQN[i]->SetLineWidth(2);
      }

      double ymxQ = 0;
      for(int i=0; i<4; i++) if(hQN[i] && hQN[i]->GetMaximum()>ymxQ) ymxQ=hQN[i]->GetMaximum();

      // Draw frame from Global (index 0), fall back to first available
      TH1F *hQFirst = hQN[0];
      if(!hQFirst) for(int i=1; i<4; i++) if(hQN[i]){ hQFirst=hQN[i]; break; }
      if(hQFirst){
        hQFirst->SetMinimum(5e-5);
        hQFirst->SetMaximum(ymxQ*5);
        hQFirst->GetXaxis()->SetTitle("q/p_{T} [(GeV/c)^{-1}]");
        hQFirst->GetYaxis()->SetTitle("Tracks / event");
        hQFirst->Draw("hist");
        for(int i=1; i<4; i++) if(hQN[i]) hQN[i]->Draw("hist same");
      }

      // Filled misID region (q/pT>0): draw Global→BLC→BLCVtx→FCSTRK so narrower on top
      for(int i=0; i<4; i++){
        if(!hQN[i]) continue;
        TH1F *hMis = (TH1F*)hQN[i]->Clone(Form("hMis%d",i));
        for(int ib=1; ib<=hMis->GetNbinsX(); ib++)
          if(hMis->GetBinCenter(ib)<=0) hMis->SetBinContent(ib,0);
        hMis->SetFillColor(QFILL[i]); hMis->SetLineColor(QCOL[i]);
        hMis->Draw("hist same");
      }

      // Redraw lines on top of fills
      for(int i=0; i<4; i++) if(hQN[i]) hQN[i]->Draw("hist same");

      // q/pT = 0 dividing line
      TLine *lz = new TLine(0, 5e-5, 0, ymxQ*3);
      lz->SetLineColor(kBlack); lz->SetLineWidth(2); lz->SetLineStyle(2); lz->Draw();

      // Colored text legend — 1 line per type, no marker
      TLatex *tleg = new TLatex(); tleg->SetNDC(); tleg->SetTextSize(0.031);
      double yleg = 0.42;
      printf("Charge misID rate (q/pT > 0 for e-):\n");
      for(int i=0; i<4; i++){
        if(!hQN[i]) continue;
        tleg->SetTextColor(QCOL[i]);
        tleg->DrawLatex(0.14, yleg,
          Form("%s  N=%.0f  misID: %.1f%%", QLABEL[i], qTot[i], qRate[i]));
        yleg -= 0.055;
        printf("  %s: %.0f / %.0f = %.2f%%\n", QLABEL[i], qMis[i], qTot[i], qRate[i]);
      }

      TLatex *lt3 = new TLatex(); lt3->SetNDC(); lt3->SetTextSize(0.040);
      lt3->SetTextColor(kBlack);
      lt3->DrawLatex(0.52, 0.80, "e^{-} only (q/p_{T}>0 is MisID)");

      c3->SaveAs("fcsTrkMomSignQ.png");
      printf("Saved fcsTrkMomSignQ.png\n");
    }

    F->Close();
}
