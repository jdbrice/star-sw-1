static const int mNCut=7;
const char* nameCut[mNCut] = {"All","ETOT","HTOT","Cone","SigmaMax","TrackMatch","ChargeSign"};
const char* TTYPE[6]={"Global","Beamline","Primary","FwdVtx","BLCVtx","FCSTRK"};
TCanvas* c1;

void plotDilep(int plt=-1, int cut=6, char* data=".", int run=1, int log=1, int trkType=1){
  if(plt==-1) {
    plotDilepX(1,4,data,run,log,trkType);
    plotDilepX(1,6,data,run,log,trkType);
    plotDilepX(2,6,data,run,log,trkType);
    plotDilepX(3,6,data,run,log,trkType);
    plotDilepX(4,6,data,run,log,trkType);
    plotDilepX(5,6,data,run,log,trkType);
    plotDilepX(6,6,data,run,log,trkType);
    plotDilepX(7,6,data,run,log,trkType);
  }else{
    plotDilepX(plt,cut,data,run,log,trkType);
  }
}

void plotDilepX(int plt=0, int cut=6, char* data=".", int run=-1, int log=1, int trkType=1){
  TString Data(data);
  const char* tag = gSystem->BaseName(data);
  const char* ttag = TTYPE[trkType];
  char file[100];
  if(run==-1){
    sprintf(file,"%s/hist_dilep/all.dilep.root",data);
  }else if(run==1){
    sprintf(file,"%s/hist_dilep/1.0.dilep.root",data);
  }else{
    sprintf(file,"%s/hist_dilep/%d.0.dilep.root",data,run);
  }
  printf("Reading %s (trkType=%d=%s)\n",file,trkType,ttag);
  TFile *F = new TFile(file,"old");

  if(c1) delete c1;
  c1 = new TCanvas("c1","FCS DiLepton",50,0,1500,1200);
  gStyle->SetLabelSize(0.1,"xy");
  gStyle->SetPalette(1);
  gStyle->SetStatW(0.4);

  TH1F* h1;
  TH2F* h2;
  char hname[100];
  if(plt==0 || plt==1) {
    c1->Clear();
    c1->Divide(2,2);

    c1->cd(1); gPad->SetLogz(); h2=(TH2F*)F->Get(Form("XFPT_%s_%s",ttag,nameCut[cut]));   if(h2) h2->Draw("colz");
    c1->cd(2); gPad->SetLogz(); h2=(TH2F*)F->Get(Form("ET12_%s_%s",ttag,nameCut[cut]));   if(h2) h2->Draw("colz");
    c1->cd(3); h2=(TH2F*)F->Get(Form("XY_%s_%s",  ttag,nameCut[cut]));   if(h2) h2->Draw("colz");
    c1->cd(4); gPad->SetLogz(); h2=(TH2F*)F->Get(Form("PTET_%s_%s",ttag,nameCut[cut]));   if(h2) h2->Draw("colz");
    c1->cd(0); TText *t = new TText(0.4,0.95,Form("%s Cut=%s",ttag,nameCut[cut])); t->Draw();

    c1->SaveAs(Form("dilep.cut%d.%s.%s.%d.png",cut,ttag,tag,run));
  }

  if(plt==0 || plt==2) {
    c1->Clear();
    c1->Divide(2,3);

    int c=0;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("RETot_%s_%s",   ttag,nameCut[c])); if(h1) h1->Draw();
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("RHTot_%s_%s",   ttag,nameCut[c])); if(h1) h1->Draw();
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("RCone_%s_%s",   ttag,nameCut[c])); if(h1) h1->Draw();
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Sigmax_%s_%s",  ttag,nameCut[c])); if(h1) h1->Draw();
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("EToverPT_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetMinimum(0.2); h1->Draw(); }
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("ChargeSum_%s_%s",ttag,nameCut[c]));if(h1){ h1->SetMinimum(0.2); h1->Draw(); }

    c=4;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("RETot_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("RHTot_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("RCone_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Sigmax_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("EToverPT_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("ChargeSum_%s_%s",ttag,nameCut[c]));if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }

    c=5;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("RETot_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("RHTot_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("RCone_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Sigmax_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("EToverPT_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("ChargeSum_%s_%s",ttag,nameCut[c]));if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }

    c=6;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("RETot_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("RHTot_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("RCone_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Sigmax_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("EToverPT_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("ChargeSum_%s_%s",ttag,nameCut[c]));if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }

    c1->cd(1);
    TText *t;
    t=new TText(0.15, 0.8,"No cut");        t->SetTextColor(kBlack);   t->SetNDC(); t->Draw();
    t=new TText(0.15, 0.7,"FCS cuts");      t->SetTextColor(kBlue);    t->SetNDC(); t->Draw();
    t=new TText(0.15, 0.6,"Track Matched"); t->SetTextColor(kMagenta); t->SetNDC(); t->Draw();
    t=new TText(0.15, 0.5,"Charge Sign");   t->SetTextColor(kRed);     t->SetNDC(); t->Draw();
    c1->SaveAs(Form("dilep1.%s.%s.%d.png",ttag,tag,run));
  }

  if(plt==0 || plt==3) {
    c1->Clear();
    c1->Divide(2,3);

    int c=0;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ET_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetMinimum(0.2); h1->Draw(); }
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("EZ_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetMinimum(0.2); h1->Draw(); }
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("M_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetMinimum(0.2); h1->Draw(); }
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Z_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetMinimum(0.2); h1->Draw(); }
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("CosT_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetMinimum(0.2); h1->Draw(); }
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("Phi_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetMinimum(0.2); h1->Draw(); }

    c=4;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ET_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("EZ_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("M_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Z_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("CosT_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("Phi_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->Draw("same"); }

    c=5;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ET_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("EZ_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("M_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Z_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("CosT_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("Phi_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->Draw("same"); }

    c=6;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ET_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("EZ_%s_%s",  ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("M_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("Z_%s_%s",   ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(5)->SetLogy(log); h1=(TH1F*)F->Get(Form("CosT_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }
    c1->cd(6)->SetLogy(log); h1=(TH1F*)F->Get(Form("Phi_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->Draw("same"); }

    c1->SaveAs(Form("dilep2.%s.%s.%d.png",ttag,tag,run));
  }

  if(plt==0 || plt==4) {
    c1->Clear();
    c1->Divide(2,2);

    int c=0;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ZVTX_%s_%s", ttag,nameCut[c])); if(h1) {h1->SetMinimum(1); h1->Draw();}
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("ZVTXA_%s_%s",ttag,nameCut[c])); if(h1) {h1->SetMinimum(1); h1->Draw();}
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("ZVTXD_%s_%s",ttag,nameCut[c])); if(h1) {h1->SetMinimum(1); h1->Draw();}
    c=4;
    c1->cd(1); h1=(TH1F*)F->Get(Form("ZVTX_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->SetMinimum(1); h1->Draw("same"); }
    c1->cd(2); h1=(TH1F*)F->Get(Form("ZVTXA_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->SetMinimum(1); h1->Draw("same"); }
    c1->cd(3); h1=(TH1F*)F->Get(Form("ZVTXD_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->SetMinimum(1); h1->Draw("same"); }
    c=5;
    c1->cd(1); h1=(TH1F*)F->Get(Form("ZVTX_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->SetMinimum(1); h1->Draw("same"); }
    c1->cd(2); h1=(TH1F*)F->Get(Form("ZVTXA_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->SetMinimum(1); h1->Draw("same"); }
    c1->cd(3); h1=(TH1F*)F->Get(Form("ZVTXD_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->SetMinimum(1); h1->Draw("same"); }
    c=6;
    c1->cd(1); h1=(TH1F*)F->Get(Form("ZVTX_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->SetMinimum(1); h1->Draw("same"); }
    c1->cd(2); h1=(TH1F*)F->Get(Form("ZVTXA_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->SetMinimum(1); h1->Draw("same"); }
    c1->cd(3); h1=(TH1F*)F->Get(Form("ZVTXD_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->SetMinimum(1); h1->Draw("same"); }

    // Panel 4: event-level BLC vertex Z, same cut/color scheme as panels 1-3.
    // Lets you directly compare the fitted BLC vertex Z spread against the
    // per-track DCAZ-based panels above, e.g. to gauge whether a
    // z-vertex-binned mixing scheme would help the mixed-event subtraction.
    c=0;
    c1->cd(4)->SetLogy(log); h1=(TH1F*)F->Get(Form("BLCVtxZ_%s",nameCut[c])); if(h1) {h1->SetMinimum(1); h1->Draw();}
    c=4;
    c1->cd(4); h1=(TH1F*)F->Get(Form("BLCVtxZ_%s",nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->SetMinimum(1); h1->Draw("same"); }
    c=5;
    c1->cd(4); h1=(TH1F*)F->Get(Form("BLCVtxZ_%s",nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->SetMinimum(1); h1->Draw("same"); }
    c=6;
    c1->cd(4); h1=(TH1F*)F->Get(Form("BLCVtxZ_%s",nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->SetMinimum(1); h1->Draw("same"); }

    c1->SaveAs(Form("dilep3.%s.%s.%d.png",ttag,tag,run));
  }

  if(plt==0 || plt==5) {
    //Mixed-event background subtraction, z-vertex x exact-charge-combo
    //binned: picoDilepton.C mixes within (zbin, combo) bins (3 z-vertex bins
    //x 4 exact charge combos -- N+S-, N-S+, N+S+, N-S-). The 2 OS combos
    //(N+S-, N-S+) are the physics channel; each is normalized against its
    //LS counterpart's same-event count (N+S- <-> N+S+, N-S+ <-> N-S-) before
    //summing across z-bins, so bins with different acceptance/statistics
    //each get their own correct normalization instead of one flat factor.
    //Normalizing to the LIKE-SIGN same-event count (not the OS same-event
    //count, and not a mass-window sideband) avoids a real bias: the OS
    //same-event count includes true dilepton signal on top of its own
    //combinatorial background, so normalizing the (purely combinatorial)
    //mixed template to match it over-subtracts everywhere except right at the
    //signal peak. Real dileptons are essentially all opposite-sign, so the
    //like-sign count is a signal-free measure of the combinatorial rate --
    //standard practice for dilepton continuum background (cf. the
    //N_comb~2*sqrt(N++ N--) like-sign method) and doesn't require guessing
    //where the signal isn't.
    //Bins with zero mixed or zero like-sign entries are skipped -- no
    //background estimate is available from them.
    const int mNZBin=3;
    const char* zbinName[mNZBin]={"z0","z1","z2"};
    const char* osCombo[2] = {"NpSm","NmSp"}; //physics-channel combos
    const char* lsCombo[2] = {"NpSp","NmSm"}; //their LS normalization reference
    TH1F* hM = (TH1F*)F->Get(Form("M_%s_%s",ttag,nameCut[6]));
    if(hM){
      hM->Sumw2();
      TH1F* hMixWeighted=(TH1F*)hM->Clone(Form("MmixWeighted_%s",ttag));
      hMixWeighted->Reset();
      hMixWeighted->Sumw2();
      double totalWeightedMix=0;
      for(int zb=0; zb<mNZBin; zb++){
        for(int o=0; o<2; o++){
          TH1F* hMixBin  = (TH1F*)F->Get(Form("MmixBin_%s_%s_%s", ttag,zbinName[zb],osCombo[o]));
          TH1F* hLSameBin= (TH1F*)F->Get(Form("MsameBin_%s_%s_%s",ttag,zbinName[zb],lsCombo[o]));
          if(!hMixBin || !hLSameBin) continue;
          double nMix  = hMixBin->Integral();
          double nLS   = hLSameBin->Integral();
          if(nMix<=0 || nLS<=0) continue;
          double w = nLS/nMix;
          hMixBin->Sumw2();
          TH1F* hMixBinScaled=(TH1F*)hMixBin->Clone(Form("MmixBinScaled_%s_%s_%s",ttag,zbinName[zb],osCombo[o]));
          hMixBinScaled->Scale(w);
          hMixWeighted->Add(hMixBinScaled);
          totalWeightedMix += nLS;
        }
      }

      TH1F* hSub=(TH1F*)hM->Clone(Form("MSub_%s",ttag));
      hSub->Add(hM,hMixWeighted,1.0,-1.0); //signal - weighted mixed, errors propagated via Sumw2

      c1->Clear();
      c1->Divide(1,2);

      c1->cd(1)->SetLogy(log);
      hM->SetMinimum(0.2);
      hM->SetLineColor(kBlack);
      hM->SetTitle(Form("Mass cut=%s %s;M [GeV];Counts",nameCut[6],ttag));
      hM->Draw("E");
      hMixWeighted->SetLineColor(kRed);
      hMixWeighted->Draw("E same");
      TLegend *legMix=new TLegend(0.55,0.7,0.88,0.88);
      legMix->AddEntry(hM,"Same-event","le");
      legMix->AddEntry(hMixWeighted,"Mixed-event bg (z/charge-binned, weighted)","le");
      legMix->Draw();

      c1->cd(2)->SetLogy(0);
      hSub->SetLineColor(kBlack);
      hSub->SetTitle(Form("Mixed-event-subtracted Mass cut=%s %s;M [GeV];Counts",nameCut[6],ttag));
      hSub->Draw("E");
      TLine *lSub0=new TLine(hSub->GetXaxis()->GetXmin(),0,hSub->GetXaxis()->GetXmax(),0);
      lSub0->SetLineStyle(2);
      lSub0->Draw();
    }

    c1->SaveAs(Form("dilep4.%s.%s.%d.png",ttag,tag,run));
  }

  if(plt==0 || plt==6) {
    //Pi0 (gamma-gamma) mass from FCS ECAL cluster pairs (RunPi0() in
    //picoDilepton.C) -- NOT track-type dependent (photons, no tracking
    //involved), so this panel is identical regardless of the trkType
    //argument. Simple inclusive combinatorial pairing, no background
    //subtraction; the true pi0 mass (0.135 GeV) is marked for
    //reference as an independent FCS ECAL calibration check.
    c1->Clear();
    c1->Divide(1,1);
    TH1F* hPi0 = (TH1F*)F->Get("Pi0Mass");
    if(hPi0){
      //c1->cd(1)->SetLogy(log);
      c1->cd(1);
      hPi0->SetMinimum(0.2);
      hPi0->SetLineColor(kBlack);
      hPi0->Draw("E");
      TLine *lPi0=new TLine(0.135,hPi0->GetMinimum(),0.135,hPi0->GetMaximum());
      lPi0->SetLineStyle(2);
      lPi0->SetLineColor(kRed);
      lPi0->Draw();
      TLegend *legPi0=new TLegend(0.5,0.55,0.88,0.68);
      legPi0->AddEntry(hPi0,"Cluster pairs (cuts in title above)","le");
      legPi0->AddEntry(lPi0,"True #pi^{0} mass (0.135 GeV)","l");
      legPi0->Draw();
    }
    c1->SaveAs(Form("pi0.%s.%d.png",tag,run));
  }

  if(plt==0 || plt==7) {
    //Per-(zbin,combo) mixed-event-subtracted mass, NOT summed across bins
    //like plt=5's dilep4 -- shows each of the 12 (3 z-vertex bins x 4 exact
    //charge combos) cells' own subtracted spectrum separately.
    //OS combo columns (N+S-, N-S+) directly test whether the negative-dip/
    //positive-bump structure seen in the summed spectrum (todo.txt #6) is
    //uniform across z-vertex (would argue against a vertex-position-
    //dependent acceptance effect) or concentrated in specific bins (would
    //support the #15 double-arm acceptance-vs-z-vertex hypothesis) --
    //normalized the same cross-referenced way as plt=5.
    //LS combo columns (N+S+, N-S-) are the null test: self-normalized (own
    //same-event count vs own mixed count), so with no real dilepton signal
    //in like-sign, these should come out flat if mixing is behaving. If they
    //show the same shape as the OS columns, that points at a mixing-
    //methodology artifact rather than physics near the OS mass features.
    const int mNZBin=3;
    const int mNCombo=4;
    const char* zbinName[mNZBin]={"z0","z1","z2"};
    const float zBinEdges[mNZBin+1]={-50.0,30.0,70.0,120.0}; //must match picoDilepton.C's mZBinEdges
    const char* comboName[mNCombo]={"NpSm","NmSp","NpSp","NmSm"};
    //normalization reference for each combo: OS combos (0,1) cross-reference
    //their LS counterpart's same-event count; LS combos (2,3) self-reference.
    const char* normCombo[mNCombo]={"NpSp","NmSm","NpSp","NmSm"};

    c1->Clear();
    c1->Divide(4,3);

    printf("=== plt=7 per-(zbin,combo) event counts (%s, cut=%s) ===\n",ttag,nameCut[6]);
    double sumSame[mNCombo]={0,0,0,0};
    for(int zb=0; zb<mNZBin; zb++){
      for(int c=0; c<mNCombo; c++){
        int pad = zb*mNCombo + c + 1;
        TH1F* hSameBin  = (TH1F*)F->Get(Form("MsameBin_%s_%s_%s", ttag,zbinName[zb],comboName[c]));
        TH1F* hMixBin   = (TH1F*)F->Get(Form("MmixBin_%s_%s_%s",  ttag,zbinName[zb],comboName[c]));
        TH1F* hNormBin  = (TH1F*)F->Get(Form("MsameBin_%s_%s_%s", ttag,zbinName[zb],normCombo[c]));
        c1->cd(pad);
        TString label = Form("z=[%.0f,%.0f) %s",zBinEdges[zb],zBinEdges[zb+1],comboName[c]);
        if(!hSameBin || !hMixBin || !hNormBin){
          printf("  %-22s : no histo\n",label.Data());
          TText *tmiss=new TText(0.5,0.5,"no histo"); tmiss->SetNDC(); tmiss->SetTextAlign(22); tmiss->Draw();
          continue;
        }
        double nSame = hSameBin->Integral();
        double nMix  = hMixBin->Integral();
        double nNorm = hNormBin->Integral();
        sumSame[c] += nSame;
        printf("  %-22s : same=%-8.0f mixed=%-8.0f norm(%s)=%-8.0f\n",
               label.Data(),nSame,nMix,normCombo[c],nNorm);
        if(nMix<=0 || nNorm<=0){
          TText *tempty=new TText(0.5,0.5,Form("%s: no bg",label.Data())); tempty->SetNDC(); tempty->SetTextAlign(22); tempty->Draw();
          continue;
        }
        double w = nNorm/nMix;
        hSameBin->Sumw2();
        hMixBin->Sumw2();
        TH1F* hMixBinScaled=(TH1F*)hMixBin->Clone(Form("MmixBinScaled7_%s_%s_%s",ttag,zbinName[zb],comboName[c]));
        hMixBinScaled->Scale(w);
        TH1F* hSubBin=(TH1F*)hSameBin->Clone(Form("MSubBin_%s_%s_%s",ttag,zbinName[zb],comboName[c]));
        hSubBin->Add(hSameBin,hMixBinScaled,1.0,-1.0);
        // Auto stat box (Entries/Mean/RMS) is misleading here: this is a
        // difference of two histograms that's *designed* to sum close to
        // zero net integral, so GetMean()=sum(x*w)/sum(w) divides by a
        // near-zero denominator and blows up to a meaningless value. Bin
        // content/errors are still correct -- just suppress the box and
        // show the (meaningful) same/mixed/norm counts instead.
        hSubBin->SetStats(0);
        hSubBin->SetTitle(Form("%s;M [GeV];Counts",label.Data()));
        hSubBin->SetLineColor(kBlack);
        hSubBin->Draw("E");
        TText *tcounts=new TText(0.15,0.85,Form("same=%.0f mix=%.0f norm=%.0f",nSame,nMix,nNorm));
        tcounts->SetNDC(); tcounts->SetTextSize(0.06); tcounts->Draw();
        TLine *lSub0=new TLine(hSubBin->GetXaxis()->GetXmin(),0,hSubBin->GetXaxis()->GetXmax(),0);
        lSub0->SetLineStyle(2);
        lSub0->Draw();
      }
    }
    printf("=== plt=7 charge-combo event counts, summed over all %d z bins (%s, cut=%s) ===\n",
           mNZBin,ttag,nameCut[6]);
    printf("  %s=%.0f  %s=%.0f  %s=%.0f  %s=%.0f\n",
           comboName[0],sumSame[0],comboName[1],sumSame[1],comboName[2],sumSame[2],comboName[3],sumSame[3]);
    c1->SaveAs(Form("dilep5.%s.%s.%d.png",ttag,tag,run));
  }
}
