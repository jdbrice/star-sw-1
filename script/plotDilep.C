static const int mNCut=7;
const char* nameCut[mNCut] = {"All","ETOT","HTOT","Cone","SigmaMax","TrackMatch","ChargeSign"};
const char* TTYPE[6]={"Global","Beamline","Primary","FwdVtx","BLCVtx","FCSTRK"};
TCanvas* c1;

void plotDilep(int plt=-1, int cut=6, char* data=".", int run=1, int log=1, int trkType=1){
  if(plt==-1) {
    plotDilepX(1,6,data,run,log,trkType);
    plotDilepX(2,6,data,run,log,trkType);
    plotDilepX(3,6,data,run,log,trkType);
    plotDilepX(4,6,data,run,log,trkType);
    plotDilepX(1,4,data,run,log,trkType);
    plotDilepX(5,6,data,run,log,trkType);
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

    c1->cd(1); h2=(TH2F*)F->Get(Form("XFPT_%s_%s",ttag,nameCut[cut]));   if(h2) h2->Draw("colz");
    c1->cd(2); h2=(TH2F*)F->Get(Form("ET12_%s_%s",ttag,nameCut[cut]));   if(h2) h2->Draw("colz");
    c1->cd(3); h2=(TH2F*)F->Get(Form("XY_%s_%s",  ttag,nameCut[cut]));   if(h2) h2->Draw("colz");
    c1->cd(4); h2=(TH2F*)F->Get(Form("PTET_%s_%s",ttag,nameCut[cut]));   if(h2) h2->Draw("colz");
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
    //Mixed-event background subtraction, z-vertex x north-charge-sign binned:
    //picoDilepton.C only mixes an event's candidates with pool entries from
    //the same (zbin,chargebin) bin (6 z-vertex bins x 2 north-charge bins =
    //12 bins), and fills a per-bin mixed histogram (MmixBin) alongside a
    //per-bin same-event LIKE-SIGN histogram (MLSameBin). Each bin's mixed
    //spectrum is normalized independently -- weight =
    //Integral(MLSameBin)/Integral(MmixBin) -- before summing across bins, so
    //bins with different acceptance/statistics (e.g. very-forward vs.
    //very-backward z-vertex) each get their own correct normalization instead
    //of one flat pool-size-based factor.
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
    const int mNZBin=6;
    const int mNChBin=2;
    const char* zbinName[mNZBin]={"z0","z1","z2","z3","z4","z5"};
    const char* chbinName[mNChBin]={"Np","Nm"};
    TH1F* hM = (TH1F*)F->Get(Form("M_%s_%s",ttag,nameCut[6]));
    if(hM){
      hM->Sumw2();
      TH1F* hMixWeighted=(TH1F*)hM->Clone(Form("MmixWeighted_%s",ttag));
      hMixWeighted->Reset();
      hMixWeighted->Sumw2();
      double totalWeightedMix=0;
      for(int zb=0; zb<mNZBin; zb++){
        for(int cb=0; cb<mNChBin; cb++){
          TH1F* hMixBin  = (TH1F*)F->Get(Form("MmixBin_%s_%s_%s",  ttag,zbinName[zb],chbinName[cb]));
          TH1F* hLSameBin= (TH1F*)F->Get(Form("MLSameBin_%s_%s_%s",ttag,zbinName[zb],chbinName[cb]));
          if(!hMixBin || !hLSameBin) continue;
          double nMix  = hMixBin->Integral();
          double nLS   = hLSameBin->Integral();
          if(nMix<=0 || nLS<=0) continue;
          double w = nLS/nMix;
          hMixBin->Sumw2();
          TH1F* hMixBinScaled=(TH1F*)hMixBin->Clone(Form("MmixBinScaled_%s_%s_%s",ttag,zbinName[zb],chbinName[cb]));
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
}
