// fstMapPlots.C -- two panels for the FST strip-mapping checks, straight from MuDst hits.
// Left : inner-outer dphi for r/z-matched pairs, against the mirrored hypothesis.
// Right: <dphi> vs the outer hit's position in its wedge -- a STEP means the gap term.
void fstMapPlots(const char* file, const char* out, int nev=4000){
  gROOT->Macro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  gStyle->SetOptStat(0);
  StMuDstMaker* mk = new StMuDstMaker(0,0,"",file,"st:MuDst.root",1);
  TH1F* hA = new TH1F("hA2","inner#minusouter, r/z matched;#Delta#phi [deg];pairs",120,-30,30);
  TH1F* hB = new TH1F("hB2","mirrored hypothesis;#Delta#phi #minus 2#delta [deg];pairs",120,-30,30);
  TProfile* pr = new TProfile("pr2","",30,-15,15,-4,4);
  int nread=0; const int MAXH=400;
  for (int iev=0; iev<nev; iev++){
    if (mk->Make()) break;
    StMuDst* d=mk->muDst(); if(!d) continue;
    StMuFstCollection* fc=d->muFstCollection(); if(!fc) continue;
    nread++;
    int nh=(int)fc->numberOfHits(); if(nh>MAXH) nh=MAXH;
    double ph[MAXH], toz[MAXH]; int dk[MAXH], inr[MAXH], sec[MAXH]; int n=0;
    for (int ih=0; ih<nh; ih++){
      StMuFstHit* fh=fc->getHit(ih); if(!fh) continue;
      TVector3 g=fh->xyz();
      double p=TMath::ATan2(g.Y(),g.X())*180./TMath::Pi(); if(p<0) p+=360.;
      double r=sqrt(g.X()*g.X()+g.Y()*g.Y()), z=g.Z(); if(fabs(z)<1) continue;
      ph[n]=p; toz[n]=r/z; dk[n]=(int)fh->getDisk();
      inr[n]=(((int)fh->getSensor())==0)?1:0; sec[n]=(int)(p/30.0); n++;
    }
    for (int i=0;i<n;i++) for (int j=0;j<n;j++){
      if (i==j) continue;
      if (!inr[i] || inr[j]) continue;
      if (dk[i]==dk[j]) continue;
      if (fabs(toz[i]-toz[j]) > 0.010) continue;
      double dp=ph[i]-ph[j];
      while(dp>180) dp-=360;  while(dp<-180) dp+=360;
      double dlIn = fmod(ph[i],30.0)-15.0, dlOut = fmod(ph[j],30.0)-15.0;
      if (sec[i]==sec[j] && fabs(dp)<35.0){ hA->Fill(dp); hB->Fill(dp-2.0*dlIn); }
      if (fabs(dp)<4.0) pr->Fill(dlOut, dp);
    }
  }
  printf(">>> events %d\n", nread);
  TCanvas* c=new TCanvas("cm","cm",1000,450); c->Divide(2,1);
  c->cd(1); gPad->SetGridx();
  hA->SetLineColor(kBlue+2); hA->SetLineWidth(2);
  hB->SetLineColor(kRed+1);  hB->SetLineWidth(2);
  hA->SetTitle("inner#minusouter FST pairs, matched in r/z;#Delta#phi  or  #Delta#phi#minus2#delta  [deg];pairs");
  hA->Draw("hist"); hB->Draw("hist same");
  TLatex t; t.SetNDC(); t.SetTextSize(0.037);
  t.SetTextColor(kBlue+2); t.DrawLatex(0.14,0.84,"#Delta#phi : mapping correct");
  t.SetTextColor(kRed+1);  t.DrawLatex(0.14,0.79,"#Delta#phi#minus2#delta : mirrored");
  t.SetTextColor(kBlack);  t.SetTextSize(0.032);
  t.DrawLatex(0.14,0.72,"sharp peak at zero #Rightarrow");
  t.DrawLatex(0.14,0.67,"inner#leftrightarrowouter mapping is right");
  c->cd(2); gPad->SetGridx(); gPad->SetGridy();
  pr->SetTitle("outer-sensor gap: a STEP at the wedge centre;#delta of the outer hit in its wedge [deg];#LT#Delta#phi#GT [deg]");
  pr->SetMarkerStyle(20); pr->SetMarkerSize(1.0); pr->SetMarkerColor(kBlack);
  pr->SetLineColor(kBlack); pr->SetLineWidth(2); pr->SetMinimum(-2.0); pr->SetMaximum(2.0);
  pr->Draw("P");
  TLine* l0=new TLine(0,-2,0,2); l0->SetLineStyle(2); l0->SetLineColor(kGray+2); l0->Draw();
  t.SetTextSize(0.032); t.SetTextColor(kBlack);
  t.DrawLatex(0.16,0.86,"flat on both sides, step #approx 1.6#circ");
  t.DrawLatex(0.16,0.81,"a start-#phi or direction error would give");
  t.DrawLatex(0.16,0.76,"an offset or a slope, not a step");
  c->SaveAs(out);
}
