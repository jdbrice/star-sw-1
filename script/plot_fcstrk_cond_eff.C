// plot_fcstrk_cond_eff.C
// Conditional efficiency: given BLCVtx track, fraction that become FCSTRK (trkType=5).
// Plots efficiency vs eta, phi, |p|, pT, and 2D eta-phi.
// Uses reconstructed track kinematics (no helix approximation needed).
//
// Usage: root4star -b -q 'plot_fcstrk_cond_eff.C("hist_match/all.match.root")'

void calcEff1D(TH1F *hEff, TH1F *hNum, TH1F *hDen) {
    // Compute bin-by-bin efficiency with binomial errors.
    // Direct calculation avoids ROOT Sumw2 complications with TH1::Divide("B").
    hEff->Reset();
    for(int ib = 1; ib <= hEff->GetNbinsX(); ib++){
        double den = hDen->GetBinContent(ib);
        double num = hNum->GetBinContent(ib);
        if(den > 0){
            double p = (num > den) ? 1.0 : num/den;
            hEff->SetBinContent(ib, p);
            hEff->SetBinError(ib, (den > 1) ? sqrt(p*(1-p)/den) : 0);
        } else {
            hEff->SetBinContent(ib, 0);
            hEff->SetBinError(ib, 0);
        }
    }
}

void plot_fcstrk_cond_eff(const char* fname = "hist_match/all.match.root") {
    TFile *F = TFile::Open(fname);
    if(!F || F->IsZombie()){ printf("Cannot open %s\n", fname); return; }
    printf("Reading %s\n", fname);

    // Load BLCVtx (denominator) and FCSTRK (numerator) 1D histograms
    TH1F *hEtaBLC  = (TH1F*)F->Get("trkEtaBLCVtx");
    TH1F *hEtaFCS  = (TH1F*)F->Get("trkEtaFCSTRK");
    TH1F *hPhiBLC  = (TH1F*)F->Get("trkPhiBLCVtx");
    TH1F *hPhiFCS  = (TH1F*)F->Get("trkPhiFCSTRK");
    TH1F *hPBLC    = (TH1F*)F->Get("trkPBLCVtx");
    TH1F *hPFCS    = (TH1F*)F->Get("trkPFCSTRK");
    TH1F *hPtBLC   = (TH1F*)F->Get("trkPtBLCVtx");
    TH1F *hPtFCS   = (TH1F*)F->Get("trkPtFCSTRK");
    TH1F *hPBLC_ng  = (TH1F*)F->Get("trkP_noNSgapBLCVtx");
    TH1F *hPFCS_ng  = (TH1F*)F->Get("trkP_noNSgapFCSTRK");
    TH1F *hPtBLC_ng = (TH1F*)F->Get("trkPt_noNSgapBLCVtx");
    TH1F *hPtFCS_ng = (TH1F*)F->Get("trkPt_noNSgapFCSTRK");
    TH2F *hEPBLC   = (TH2F*)F->Get("trkEtaPhiBLCVtx");
    TH2F *hEPFCS   = (TH2F*)F->Get("trkEtaPhiFCSTRK");

    if(!hEtaBLC || !hEtaFCS){ printf("Missing eta histograms — re-run picoMatch.C first\n"); return; }

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(55);   // rainbow palette for 2D

    double nBLC = hEtaBLC->GetEntries();
    double nFCS = hEtaFCS->GetEntries();
    printf("BLCVtx tracks: %.0f   FCSTRK tracks: %.0f   overall: %.1f%%\n",
           nBLC, nFCS, nBLC>0 ? 100.*nFCS/nBLC : 0);

    // ── Page 1: efficiency vs eta and phi (1D) ─────────────────────────────
    {
      TCanvas *c1 = new TCanvas("cCondEff1","",1000,450);
      c1->Divide(2,1);

      // eta
      c1->cd(1);
      gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
      TH1F *hEffEta = (TH1F*)hEtaBLC->Clone("hEffEta");
      calcEff1D(hEffEta, hEtaFCS, hEtaBLC);
      hEffEta->SetLineColor(kRed); hEffEta->SetMarkerColor(kRed); hEffEta->SetMarkerStyle(20);
      hEffEta->GetYaxis()->SetRangeUser(0, 1.3);
      hEffEta->GetXaxis()->SetTitle("#eta");
      hEffEta->GetYaxis()->SetTitle("FCSTRK / BLCVtx");
      hEffEta->Draw("E");
      TLine *l1 = new TLine(hEffEta->GetXaxis()->GetXmin(),1.0,hEffEta->GetXaxis()->GetXmax(),1.0);
      l1->SetLineColor(kGray+2); l1->SetLineStyle(2); l1->Draw();
      TLatex *lt = new TLatex(); lt->SetNDC(); lt->SetTextSize(0.045);
      lt->DrawLatex(0.18,0.87,Form("Overall: %.1f%%",nBLC>0?100.*nFCS/nBLC:0));

      // phi
      c1->cd(2);
      gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
      TH1F *hEffPhi = (TH1F*)hPhiBLC->Clone("hEffPhi");
      calcEff1D(hEffPhi, hPhiFCS, hPhiBLC);
      hEffPhi->SetLineColor(kRed); hEffPhi->SetMarkerColor(kRed); hEffPhi->SetMarkerStyle(20);
      hEffPhi->GetYaxis()->SetRangeUser(0, 1.3);
      hEffPhi->GetXaxis()->SetTitle("#phi [rad]");
      hEffPhi->GetYaxis()->SetTitle("FCSTRK / BLCVtx");
      hEffPhi->Draw("E");
      TLine *l2 = new TLine(hEffPhi->GetXaxis()->GetXmin(),1.0,hEffPhi->GetXaxis()->GetXmax(),1.0);
      l2->SetLineColor(kGray+2); l2->SetLineStyle(2); l2->Draw();
      // NS gap at phi = +/-pi/2 (horizontal split between North and South halves)
      TLine *lgap1 = new TLine( TMath::Pi()/2, 0,  TMath::Pi()/2, 1.3);
      TLine *lgap2 = new TLine(-TMath::Pi()/2, 0, -TMath::Pi()/2, 1.3);
      lgap1->SetLineColor(kBlack); lgap1->SetLineStyle(3); lgap1->Draw();
      lgap2->SetLineColor(kBlack); lgap2->SetLineStyle(3); lgap2->Draw();
      TLatex *ltg = new TLatex(); ltg->SetNDC(); ltg->SetTextSize(0.036);
      ltg->DrawLatex(0.50, 0.87, "dotted: #phi=#pm#pi/2 (NS gap)");

      c1->SaveAs("fcsTrkCondEff1D.png");
      printf("Saved fcsTrkCondEff1D.png\n");
    }

    // ── Page 2: efficiency vs |p| and pT (all phi, top) and away NS gap (bottom) ──
    if(hPBLC && hPFCS){
      TCanvas *c2 = new TCanvas("cCondEff2","",1000,900);
      c2->Divide(2,2);

      // helper macro to draw one efficiency pad
      // (ROOT 5 CINT: no lambdas; inline the four times)

      // [1] |p| all phi
      c2->cd(1);
      gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
      TH1F *hEffP = (TH1F*)hPBLC->Clone("hEffP");
      calcEff1D(hEffP, hPFCS, hPBLC);
      hEffP->SetLineColor(kRed); hEffP->SetMarkerColor(kRed); hEffP->SetMarkerStyle(20);
      hEffP->GetYaxis()->SetRangeUser(0, 1.3);
      hEffP->GetXaxis()->SetTitle("|p| [GeV/c]");
      hEffP->GetYaxis()->SetTitle("FCSTRK / BLCVtx");
      hEffP->Draw("E");
      { TLine *l=new TLine(hEffP->GetXaxis()->GetXmin(),1.0,hEffP->GetXaxis()->GetXmax(),1.0);
        l->SetLineColor(kGray+2); l->SetLineStyle(2); l->Draw(); }
      { TLatex *lt=new TLatex(); lt->SetNDC(); lt->SetTextSize(0.042);
        lt->DrawLatex(0.18,0.87,"all #phi"); }

      // [2] pT all phi
      c2->cd(2);
      gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
      TH1F *hEffPt = (TH1F*)hPtBLC->Clone("hEffPt");
      calcEff1D(hEffPt, hPtFCS, hPtBLC);
      hEffPt->SetLineColor(kRed); hEffPt->SetMarkerColor(kRed); hEffPt->SetMarkerStyle(20);
      hEffPt->GetYaxis()->SetRangeUser(0, 1.3);
      hEffPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
      hEffPt->GetYaxis()->SetTitle("FCSTRK / BLCVtx");
      hEffPt->Draw("E");
      { TLine *l=new TLine(hEffPt->GetXaxis()->GetXmin(),1.0,hEffPt->GetXaxis()->GetXmax(),1.0);
        l->SetLineColor(kGray+2); l->SetLineStyle(2); l->Draw(); }
      { TLatex *lt=new TLatex(); lt->SetNDC(); lt->SetTextSize(0.042);
        lt->DrawLatex(0.18,0.87,"all #phi"); }

      // [3] |p| away NS gap
      if(hPBLC_ng && hPFCS_ng){
        c2->cd(3);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
        TH1F *hEffPng = (TH1F*)hPBLC_ng->Clone("hEffPng");
        calcEff1D(hEffPng, hPFCS_ng, hPBLC_ng);
        hEffPng->SetLineColor(kBlue); hEffPng->SetMarkerColor(kBlue); hEffPng->SetMarkerStyle(20);
        hEffPng->GetYaxis()->SetRangeUser(0, 1.3);
        hEffPng->GetXaxis()->SetTitle("|p| [GeV/c]");
        hEffPng->GetYaxis()->SetTitle("FCSTRK / BLCVtx");
        hEffPng->Draw("E");
        { TLine *l=new TLine(hEffPng->GetXaxis()->GetXmin(),1.0,hEffPng->GetXaxis()->GetXmax(),1.0);
          l->SetLineColor(kGray+2); l->SetLineStyle(2); l->Draw(); }
        { TLatex *lt=new TLatex(); lt->SetNDC(); lt->SetTextSize(0.042);
          lt->DrawLatex(0.18,0.87,"|#phi|<0.9 or |#phi|>2.5"); }
      }

      // [4] pT away NS gap
      if(hPtBLC_ng && hPtFCS_ng){
        c2->cd(4);
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
        TH1F *hEffPtng = (TH1F*)hPtBLC_ng->Clone("hEffPtng");
        calcEff1D(hEffPtng, hPtFCS_ng, hPtBLC_ng);
        hEffPtng->SetLineColor(kBlue); hEffPtng->SetMarkerColor(kBlue); hEffPtng->SetMarkerStyle(20);
        hEffPtng->GetYaxis()->SetRangeUser(0, 1.3);
        hEffPtng->GetXaxis()->SetTitle("p_{T} [GeV/c]");
        hEffPtng->GetYaxis()->SetTitle("FCSTRK / BLCVtx");
        hEffPtng->Draw("E");
        { TLine *l=new TLine(hEffPtng->GetXaxis()->GetXmin(),1.0,hEffPtng->GetXaxis()->GetXmax(),1.0);
          l->SetLineColor(kGray+2); l->SetLineStyle(2); l->Draw(); }
        { TLatex *lt=new TLatex(); lt->SetNDC(); lt->SetTextSize(0.042);
          lt->DrawLatex(0.18,0.87,"|#phi|<0.9 or |#phi|>2.5"); }
      }

      c2->SaveAs("fcsTrkCondEffMom.png");
      printf("Saved fcsTrkCondEffMom.png\n");
    }

    // ── Page 3: 2D efficiency eta vs phi ───────────────────────────────────
    if(hEPBLC && hEPFCS){
      TCanvas *c4 = new TCanvas("cCondEff2D","",800,550);
      gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.14); gPad->SetBottomMargin(0.13);

      TH2F *hEff2D = (TH2F*)hEPFCS->Clone("hEff2D");
      hEff2D->Divide(hEPBLC);   // simple bin-by-bin ratio
      // zero out bins where denominator < 3 (unreliable)
      for(int ix=1; ix<=hEff2D->GetNbinsX(); ix++)
        for(int iy=1; iy<=hEff2D->GetNbinsY(); iy++)
          if(hEPBLC->GetBinContent(ix,iy) < 3)
            hEff2D->SetBinContent(ix,iy, 0);
      hEff2D->SetMinimum(0);
      hEff2D->SetMaximum(1);
      hEff2D->SetTitle("FCSTRK / BLCVtx conditional efficiency; #eta; #phi [rad]");
      hEff2D->Draw("colz");
      TLatex *lt4 = new TLatex(); lt4->SetNDC(); lt4->SetTextSize(0.038);
      lt4->DrawLatex(0.12,0.92,"FCSTRK / BLCVtx (blank = <3 BLCVtx tracks)");

      c4->SaveAs("fcsTrkCondEff2D.png");
      printf("Saved fcsTrkCondEff2D.png\n");
    }

    F->Close();
    printf("Done.\n");
}
