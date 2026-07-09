static const int mNCut=7;
const char* nameCut[mNCut] = {"All","ETOT","HTOT","Cone","SigmaMax","TrackMatch","ChargeSign"};
const char* TTYPE[5]={"Global","Beamline","Primary","FwdVtx","BLCVtx"};
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
    c1->Divide(1,3);

    int c=0;
    c1->cd(1)->SetLogy(log); h1=(TH1F*)F->Get(Form("ZVTX_%s_%s", ttag,nameCut[c])); if(h1) {h1->SetMinumum(1); h1->Draw();}
    c1->cd(2)->SetLogy(log); h1=(TH1F*)F->Get(Form("ZVTXA_%s_%s",ttag,nameCut[c])); if(h1) {h1->SetMinumum(1); h1->Draw();}
    c1->cd(3)->SetLogy(log); h1=(TH1F*)F->Get(Form("ZVTXD_%s_%s",ttag,nameCut[c])); if(h1) {h1->SetMinumum(1); h1->Draw();}
    c=4;
    c1->cd(1); h1=(TH1F*)F->Get(Form("ZVTX_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->SetMinumum(1); h1->Draw("same"); }
    c1->cd(2); h1=(TH1F*)F->Get(Form("ZVTXA_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->SetMinumum(1); h1->Draw("same"); }
    c1->cd(3); h1=(TH1F*)F->Get(Form("ZVTXD_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kBlue); h1->SetMinumum(1); h1->Draw("same"); }
    c=5;
    c1->cd(1); h1=(TH1F*)F->Get(Form("ZVTX_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->SetMinumum(1); h1->Draw("same"); }
    c1->cd(2); h1=(TH1F*)F->Get(Form("ZVTXA_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->SetMinumum(1); h1->Draw("same"); }
    c1->cd(3); h1=(TH1F*)F->Get(Form("ZVTXD_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kMagenta); h1->SetMinumum(1); h1->Draw("same"); }
    c=6;
    c1->cd(1); h1=(TH1F*)F->Get(Form("ZVTX_%s_%s", ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->SetMinumum(1); h1->Draw("same"); }
    c1->cd(2); h1=(TH1F*)F->Get(Form("ZVTXA_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->SetMinumum(1); h1->Draw("same"); }
    c1->cd(3); h1=(TH1F*)F->Get(Form("ZVTXD_%s_%s",ttag,nameCut[c])); if(h1){ h1->SetLineColor(kRed); h1->SetMinumum(1); h1->Draw("same"); }

    c1->SaveAs(Form("dilep3.%s.%s.%d.png",ttag,tag,run));
  }

  if(plt==0 || plt==5) {
    //Mixed-event background subtraction: Mmix has 2 mixed pairs filled per
    //source event (this-N+prev-S and prev-N+this-S), so normalize by 1/2 to
    //put it on the same per-event footing as the same-event M spectrum.
    TH1F* hM   = (TH1F*)F->Get(Form("M_%s_%s",ttag,nameCut[6]));
    TH1F* hMix = (TH1F*)F->Get(Form("Mmix_%s",ttag));
    if(hM && hMix){
      hM->Sumw2();
      hMix->Sumw2();
      TH1F* hMixScaled=(TH1F*)hMix->Clone(Form("MmixScaled_%s",ttag));
      hMixScaled->Scale(0.5);

      TH1F* hSub=(TH1F*)hM->Clone(Form("MSub_%s",ttag));
      hSub->Add(hM,hMix,1.0,-0.5); //signal - 0.5*mixed, errors propagated via Sumw2

      c1->Clear();
      c1->Divide(1,2);

      c1->cd(1)->SetLogy(log);
      hM->SetMinimum(0.2);
      hM->SetLineColor(kBlack);
      hM->SetTitle(Form("Mass cut=%s %s;M [GeV];Counts",nameCut[6],ttag));
      hM->Draw("E");
      hMixScaled->SetLineColor(kRed);
      hMixScaled->Draw("E same");
      TLegend *legMix=new TLegend(0.55,0.7,0.88,0.88);
      legMix->AddEntry(hM,"Same-event","le");
      legMix->AddEntry(hMixScaled,"Mixed-event bg (#times1/2)","le");
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
