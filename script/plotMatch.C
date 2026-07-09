#include <iostream>

static const int mNType=6;
const char *TTYPE[mNType]={"Global","Beamline","Primary","FwdVtx","BLCVtx","FCSTRK"};
// Colors per track type: black, blue, red, green, magenta, cyan
static const int TCOL[mNType]={kBlack, kBlue, kRed, kGreen+2, kMagenta+1, kCyan+2};

static const int mNCut=3;
const char *EH[2]={"Ecal","Hcal"};
const char *NSTB[4][2]={{"North","South"},{"Top","Bottom"},{"North","South"},{"R<70","R>70"}};
const char *CUT[4][mNCut+1]={{"Same Event","Mixed Event","Same-Mixed North","Same-Mixed South"},
                             {"Same Event","Mixed Event","Same-Mixed Top","Same-Mixed Bottom"},
                             {"Same Event","Mixed Event","Same-Mixed North","Same-Mixed South"},
                             {"Same Event","Mixed Event","Same-Mixed R<70","Same-Mixed R>70"}};

void plotMatch(char* data="202604",int run=0, int set=0, int trkType=4){
  //void plotMatch(char* data=".",int run=1,int set=0){
  const char *ttag = TTYPE[trkType];
  char file[100];
  if(run==0){
    sprintf(file,"%s/hist_match/all.match.root",data);
  }else{
    sprintf(file,"%s/hist_match/%d.%d.match.root",data,run,set);
  }
  printf("Reading %s (trkType=%d=%s)\n",file,trkType,ttag);
  TFile *F = new TFile(file,"old");

  TH1F *hTrkQPt[mNType], *hTrkEta[mNType];
  TH1F *hTrkPhi[mNType], *hTrkDcaZ[mNType];
  // MC truth-matching histograms (only present in simulation match files)
  TH1F *hTrkDqPt[mNType], *hTrkDpTrel[mNType];
  TH2F *hTrkQPtRecovsTru[mNType];

  TH2F *hXY[2];
  TH1F *hE[3];
  TH1F *hET[3];
  TH1F *hdx[2][2][3], *hdy[2][2][3], *hdp[2][2][3], *hdr[2][2][3];
  TH2F *hxdx[2][3], *hydx[2][3], *hpdx[2][3], *hrdx[2][3];
  TH2F *hxdy[2][3], *hydy[2][3], *hpdy[2][3], *hrdy[2][3];
  TH2F *hxdp[2][3], *hydp[2][3], *hpdp[2][3], *hrdp[2][3];
  TH2F *hxdr[2][3], *hydr[2][3], *hpdr[2][3], *hrdr[2][3];

  int CUTCOL[4]={kBlack, kBlue, kRed, kGreen+2};

  // Track-quality histograms by type
  for(int tt=0; tt<mNType; tt++){
    hTrkQPt[tt]  = (TH1F*)F->Get(Form("trkQPt%s",  TTYPE[tt]));
    hTrkEta[tt]  = (TH1F*)F->Get(Form("trkEta%s",  TTYPE[tt]));
    hTrkPhi[tt]  = (TH1F*)F->Get(Form("trkPhi%s",  TTYPE[tt]));
    hTrkDcaZ[tt] = (TH1F*)F->Get(Form("trkDcaZ%s", TTYPE[tt]));
    if(hTrkQPt[tt])  hTrkQPt[tt] ->SetLineColor(TCOL[tt]);
    if(hTrkEta[tt])  hTrkEta[tt] ->SetLineColor(TCOL[tt]);
    if(hTrkPhi[tt])  hTrkPhi[tt] ->SetLineColor(TCOL[tt]);
    if(hTrkDcaZ[tt]) hTrkDcaZ[tt]->SetLineColor(TCOL[tt]);
    hTrkDqPt[tt]         = (TH1F*)F->Get(Form("trkDqPt%s",         TTYPE[tt]));
    hTrkDpTrel[tt]       = (TH1F*)F->Get(Form("trkDpTrel%s",       TTYPE[tt]));
    hTrkQPtRecovsTru[tt] = (TH2F*)F->Get(Form("trkQPtRecovsTru%s", TTYPE[tt]));
    if(hTrkDqPt[tt])         hTrkDqPt[tt]        ->SetLineColor(TCOL[tt]);
    if(hTrkDpTrel[tt])       hTrkDpTrel[tt]      ->SetLineColor(TCOL[tt]);
  }

  // Matching histograms are booked per track type; select trkType's slice via ttag
  hXY[0] = (TH2F*)F->Get("xyEcal");
  hXY[1] = (TH2F*)F->Get("xyHcal");
  TH2F *hXYETrk = (TH2F*)F->Get(Form("xyETrk%s",ttag));
  TH2F *hXYHTrk = (TH2F*)F->Get(Form("xyHTrk%s",ttag));
  hE[2]  = (TH1F*)F->Get(Form("ETrk%s", ttag));
  hET[2] = (TH1F*)F->Get(Form("PTTrk%s",ttag));
  for(int eh=0; eh<2; eh++){
    hE[eh]  = (TH1F*)F->Get(Form("E%s",EH[eh]));
    hET[eh] = (TH1F*)F->Get(Form("ET%s",EH[eh]));
    for(int cut=0; cut<mNCut-1; cut++){
      hxdx[eh][cut] = (TH2F*)F->Get(Form("xdx%s%s_%s",EH[eh],CUT[0][cut],ttag));
      hydx[eh][cut] = (TH2F*)F->Get(Form("ydx%s%s_%s",EH[eh],CUT[0][cut],ttag));
      hpdx[eh][cut] = (TH2F*)F->Get(Form("pdx%s%s_%s",EH[eh],CUT[0][cut],ttag));
      hrdx[eh][cut] = (TH2F*)F->Get(Form("rdx%s%s_%s",EH[eh],CUT[0][cut],ttag));
      hxdy[eh][cut] = (TH2F*)F->Get(Form("xdy%s%s_%s",EH[eh],CUT[1][cut],ttag));
      hydy[eh][cut] = (TH2F*)F->Get(Form("ydy%s%s_%s",EH[eh],CUT[1][cut],ttag));
      hpdy[eh][cut] = (TH2F*)F->Get(Form("pdy%s%s_%s",EH[eh],CUT[1][cut],ttag));
      hrdy[eh][cut] = (TH2F*)F->Get(Form("rdy%s%s_%s",EH[eh],CUT[1][cut],ttag));
      hxdp[eh][cut] = (TH2F*)F->Get(Form("xdp%s%s_%s",EH[eh],CUT[2][cut],ttag));
      hydp[eh][cut] = (TH2F*)F->Get(Form("ydp%s%s_%s",EH[eh],CUT[2][cut],ttag));
      hpdp[eh][cut] = (TH2F*)F->Get(Form("pdp%s%s_%s",EH[eh],CUT[2][cut],ttag));
      hrdp[eh][cut] = (TH2F*)F->Get(Form("rdp%s%s_%s",EH[eh],CUT[2][cut],ttag));
      hxdr[eh][cut] = (TH2F*)F->Get(Form("xdr%s%s_%s",EH[eh],CUT[3][cut],ttag));
      hydr[eh][cut] = (TH2F*)F->Get(Form("ydr%s%s_%s",EH[eh],CUT[3][cut],ttag));
      hpdr[eh][cut] = (TH2F*)F->Get(Form("pdr%s%s_%s",EH[eh],CUT[3][cut],ttag));
      hrdr[eh][cut] = (TH2F*)F->Get(Form("rdr%s%s_%s",EH[eh],CUT[3][cut],ttag));
      for(int nstb=0; nstb<2; nstb++){
        hdx[eh][nstb][cut] = (TH1F*)F->Get(Form("%s%sdx%s_%s",EH[eh],NSTB[0][nstb],CUT[0][cut],ttag));
        hdy[eh][nstb][cut] = (TH1F*)F->Get(Form("%s%sdy%s_%s",EH[eh],NSTB[1][nstb],CUT[1][cut],ttag));
        hdp[eh][nstb][cut] = (TH1F*)F->Get(Form("%s%sdp%s_%s",EH[eh],NSTB[2][nstb],CUT[2][cut],ttag));
        hdr[eh][nstb][cut] = (TH1F*)F->Get(Form("%s%sdr%s_%s",EH[eh],NSTB[3][nstb],CUT[3][cut],ttag));
        hdx[eh][nstb][cut]->SetLineColor(CUTCOL[cut]);
        hdy[eh][nstb][cut]->SetLineColor(CUTCOL[cut]);
        hdp[eh][nstb][cut]->SetLineColor(CUTCOL[cut]);
        hdr[eh][nstb][cut]->SetLineColor(CUTCOL[cut]);
      }
    }
  }
  
  for(int eh=0; eh<2; eh++){
    for(int nstb=0; nstb<2; nstb++){
      hdx[eh][nstb][2]= (TH1F*)hdx[eh][nstb][0]->Clone(); hdx[eh][nstb][2]->Add(hdx[eh][nstb][1],-1); hdx[eh][nstb][2]->SetLineColor(CUTCOL[2+nstb]);
      hdy[eh][nstb][2]= (TH1F*)hdy[eh][nstb][0]->Clone(); hdy[eh][nstb][2]->Add(hdy[eh][nstb][1],-1); hdy[eh][nstb][2]->SetLineColor(CUTCOL[2+nstb]);
      hdp[eh][nstb][2]= (TH1F*)hdp[eh][nstb][0]->Clone(); hdp[eh][nstb][2]->Add(hdp[eh][nstb][1],-1); hdp[eh][nstb][2]->SetLineColor(CUTCOL[2+nstb]);
      hdr[eh][nstb][2]= (TH1F*)hdr[eh][nstb][0]->Clone(); hdr[eh][nstb][2]->Add(hdr[eh][nstb][1],-1); hdr[eh][nstb][2]->SetLineColor(CUTCOL[2+nstb]);
    }
    hdx[eh][0][0]->Add(hdx[eh][1][0]);
    hdy[eh][0][0]->Add(hdy[eh][1][0]);
    hdp[eh][0][0]->Add(hdp[eh][1][0]);
    hdr[eh][0][0]->Add(hdr[eh][1][0]);
    hdx[eh][0][1]->Add(hdx[eh][1][1]);
    hdy[eh][0][1]->Add(hdy[eh][1][1]);
    hdp[eh][0][1]->Add(hdp[eh][1][1]);
    hdr[eh][0][1]->Add(hdr[eh][1][1]);
    hxdx[eh][2]= (TH2F*)hxdx[eh][0]->Clone(); hxdx[eh][2]->Add(hxdx[eh][1],-1);
    hydx[eh][2]= (TH2F*)hydx[eh][0]->Clone(); hydx[eh][2]->Add(hydx[eh][1],-1);
    hpdx[eh][2]= (TH2F*)hpdx[eh][0]->Clone(); hpdx[eh][2]->Add(hpdx[eh][1],-1);
    hrdx[eh][2]= (TH2F*)hrdx[eh][0]->Clone(); hrdx[eh][2]->Add(hrdx[eh][1],-1);
    hxdy[eh][2]= (TH2F*)hxdy[eh][0]->Clone(); hxdy[eh][2]->Add(hxdy[eh][1],-1);
    hydy[eh][2]= (TH2F*)hydy[eh][0]->Clone(); hydy[eh][2]->Add(hydy[eh][1],-1);
    hpdy[eh][2]= (TH2F*)hpdy[eh][0]->Clone(); hpdy[eh][2]->Add(hpdy[eh][1],-1);
    hrdy[eh][2]= (TH2F*)hrdy[eh][0]->Clone(); hrdy[eh][2]->Add(hrdy[eh][1],-1);
    hxdp[eh][2]= (TH2F*)hxdp[eh][0]->Clone(); hxdp[eh][2]->Add(hxdp[eh][1],-1);
    hydp[eh][2]= (TH2F*)hydp[eh][0]->Clone(); hydp[eh][2]->Add(hydp[eh][1],-1);
    hpdp[eh][2]= (TH2F*)hpdp[eh][0]->Clone(); hpdp[eh][2]->Add(hpdp[eh][1],-1);
    hrdp[eh][2]= (TH2F*)hrdp[eh][0]->Clone(); hrdp[eh][2]->Add(hrdp[eh][1],-1);
    hxdr[eh][2]= (TH2F*)hxdr[eh][0]->Clone(); hxdr[eh][2]->Add(hxdr[eh][1],-1);
    hydr[eh][2]= (TH2F*)hydr[eh][0]->Clone(); hydr[eh][2]->Add(hydr[eh][1],-1);
    hpdr[eh][2]= (TH2F*)hpdr[eh][0]->Clone(); hpdr[eh][2]->Add(hpdr[eh][1],-1);
    hrdr[eh][2]= (TH2F*)hrdr[eh][0]->Clone(); hrdr[eh][2]->Add(hrdr[eh][1],-1);
  }

  TCanvas *c1=new TCanvas("FCSTRK", "FCSTRK",50,0,1200,1200);
  TText* t;
  float lx=0.15,ly0=0.8,ldy=0.04, ly;

  c1->Divide(2,2);
  c1->cd(1); hdx[0][0][0]->Draw(); hdx[0][0][1]->Draw("same"); hdx[0][0][2]->Draw("same"); hdx[0][1][2]->Draw("same"); 
  c1->cd(2); hdy[0][0][0]->Draw(); hdy[0][0][1]->Draw("same"); hdy[0][0][2]->Draw("same"); hdy[0][1][2]->Draw("same");
  c1->cd(3); hdp[0][0][0]->Draw(); hdp[0][0][1]->Draw("same"); hdp[0][0][2]->Draw("same"); hdp[0][1][2]->Draw("same"); 
  c1->cd(4); hdr[0][0][0]->Draw(); hdr[0][0][1]->Draw("same"); hdr[0][0][2]->Draw("same"); hdr[0][1][2]->Draw("same");
  for(int j=0; j<4; j++){
    c1->cd(j+1);
    ly=ly0;
    for(int i=0; i<mNCut+1; i++){
      t=new TText(lx,ly,CUT[j][i]); 
      t->SetNDC(); t->SetTextColor(CUTCOL[i]); t->SetTextSize(0.03); 
      t->Draw();
      ly-=ldy;
    }
  }
  c1->SaveAs(Form("fcsTrkMatchEcal.%s.png",ttag));

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hdx[1][0][0]->Draw(); hdx[1][0][1]->Draw("same"); hdx[1][0][2]->Draw("same"); hdx[1][1][2]->Draw("same"); 
  c1->cd(2); hdy[1][0][0]->Draw(); hdy[1][0][1]->Draw("same"); hdy[1][0][2]->Draw("same"); hdy[1][1][2]->Draw("same");
  c1->cd(3); hdp[1][0][0]->Draw(); hdp[1][0][1]->Draw("same"); hdp[1][0][2]->Draw("same"); hdp[1][1][2]->Draw("same"); 
  c1->cd(4); hdr[1][0][0]->Draw(); hdr[1][0][1]->Draw("same"); hdr[1][0][2]->Draw("same"); hdr[1][1][2]->Draw("same");
  ly=ly0;
  for(int j=0; j<4; j++){
    c1->cd(j+1);
    ly=ly0;
    for(int i=0; i<mNCut+1; i++){      
      t=new TText(lx,ly,CUT[j][i]); 
      t->SetNDC(); t->SetTextColor(CUTCOL[i]); t->SetTextSize(0.03); 
      t->Draw();
      ly-=ldy;
    }
  }
  c1->SaveAs(Form("fcsTrkMatchHcal.%s.png",ttag));

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hXY[0]->Draw("colz");
  c1->cd(2); hXY[1]->Draw("colz");
  c1->cd(3); if(hXYETrk) hXYETrk->Draw("colz");
  c1->cd(4); if(hXYHTrk) hXYHTrk->Draw("colz");
  c1->SaveAs(Form("fcsTrkxy.%s.png",ttag));

  c1->Clear();
  c1->Divide(3,2);
  c1->cd(1); hET[0]->Draw();
  c1->cd(2); hET[1]->Draw();
  c1->cd(3); hET[2]->Draw();
  c1->cd(4); hE[0]->Draw();
  c1->cd(5); hE[1]->Draw();
  c1->cd(6); hE[2]->Draw();
  c1->SaveAs(Form("fcsTrkEt.%s.png",ttag));

  // ── Ecal dx/dy 2D with slope measurement via ProfileX ──────────────────
  // Helper: draw 2D colz + ProfileX fitted with linear f(x)=p0+p1*x
  // Returns slope p1 and prints on pad.
  TF1 *flin = new TF1("flin","[0]+[1]*x",-150,150);
  flin->SetLineColor(kRed); flin->SetLineWidth(2);

  // Measure slope of dx vs FcsX and dy vs FcsY for same/mixed/signal.
  // Profile same [0] and mixed [1] separately; signal = same - mixed.
  // Printed to stdout only — not drawn on plots.
  { TF1 *fs=new TF1("fls","[0]+[1]*x",-150,150);
    TF1 *fm=new TF1("flm","[0]+[1]*x",-150,150);
    TProfile *ps, *pm;
    double ss, ms, sig;
    // Ecal dx vs FcsX
    ps=(TProfile*)hxdx[0][0]->ProfileX("_pxs",1,-1,"s"); ps->Fit("fls","Q0");
    pm=(TProfile*)hxdx[0][1]->ProfileX("_pxm",1,-1,"s"); pm->Fit("flm","Q0");
    ss=fs->GetParameter(1); ms=fm->GetParameter(1); sig=ss-ms;
    printf("Ecal dx/dFcsX:  same=%.4f  mixed=%.4f  signal=%.4f\n",ss,ms,sig);
    // Ecal dy vs FcsY
    ps=(TProfile*)hydy[0][0]->ProfileX("_pys",1,-1,"s"); ps->Fit("fls","Q0");
    pm=(TProfile*)hydy[0][1]->ProfileX("_pym",1,-1,"s"); pm->Fit("flm","Q0");
    ss=fs->GetParameter(1); ms=fm->GetParameter(1); sig=ss-ms;
    printf("Ecal dy/dFcsY:  same=%.4f  mixed=%.4f  signal=%.4f\n",ss,ms,sig);
    // Hcal dx vs FcsX
    ps=(TProfile*)hxdx[1][0]->ProfileX("_hpxs",1,-1,"s"); ps->Fit("fls","Q0");
    pm=(TProfile*)hxdx[1][1]->ProfileX("_hpxm",1,-1,"s"); pm->Fit("flm","Q0");
    ss=fs->GetParameter(1); ms=fm->GetParameter(1); sig=ss-ms;
    printf("Hcal dx/dFcsX:  same=%.4f  mixed=%.4f  signal=%.4f\n",ss,ms,sig);
    // Hcal dy vs FcsY
    ps=(TProfile*)hydy[1][0]->ProfileX("_hpys",1,-1,"s"); ps->Fit("fls","Q0");
    pm=(TProfile*)hydy[1][1]->ProfileX("_hpym",1,-1,"s"); pm->Fit("flm","Q0");
    ss=fs->GetParameter(1); ms=fm->GetParameter(1); sig=ss-ms;
    printf("Hcal dy/dFcsY:  same=%.4f  mixed=%.4f  signal=%.4f\n",ss,ms,sig);
    delete fs; delete fm; }

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdx[0][2]->Draw("colz");
  c1->cd(2); hydx[0][2]->Draw("colz");
  c1->cd(3); hpdx[0][2]->Draw("colz");
  c1->cd(4); hrdx[0][2]->Draw("colz");
  c1->SaveAs(Form("fcsTrkEcaldx.%s.png",ttag));
  c1->Clear();

  c1->Divide(2,2);
  c1->cd(1); hxdy[0][2]->Draw("colz");
  c1->cd(2); hydy[0][2]->Draw("colz");
  c1->cd(3); hpdy[0][2]->Draw("colz");
  c1->cd(4); hrdy[0][2]->Draw("colz");
  c1->SaveAs(Form("fcsTrkEcaldy.%s.png",ttag));

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdp[0][2]->Draw("colz");
  c1->cd(2); hydp[0][2]->Draw("colz");
  c1->cd(3); hpdp[0][2]->Draw("colz");
  c1->cd(4); hrdp[0][2]->Draw("colz");
  c1->SaveAs(Form("fcsTrkEcaldp.%s.png",ttag));

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdr[0][2]->Draw("colz");
  c1->cd(2); hydr[0][2]->Draw("colz");
  c1->cd(3); hpdr[0][2]->Draw("colz");
  c1->cd(4); hrdr[0][2]->Draw("colz");
  c1->SaveAs(Form("fcsTrkEcaldr.%s.png",ttag));

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdx[1][2]->Draw("colz");
  c1->cd(2); hydx[1][2]->Draw("colz");
  c1->cd(3); hpdx[1][2]->Draw("colz");
  c1->cd(4); hrdx[1][2]->Draw("colz");
  c1->SaveAs(Form("fcsTrkHcaldx.%s.png",ttag));
  c1->Clear();

  c1->Divide(2,2);
  c1->cd(1); hxdy[1][2]->Draw("colz");
  c1->cd(2); hydy[1][2]->Draw("colz");
  c1->cd(3); hpdy[1][2]->Draw("colz");
  c1->cd(4); hrdy[1][2]->Draw("colz");
  c1->SaveAs(Form("fcsTrkHcaldy.%s.png",ttag));

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdp[1][2]->Draw("colz");
  c1->cd(2); hydp[1][2]->Draw("colz");
  c1->cd(3); hpdp[1][2]->Draw("colz");
  c1->cd(4); hrdp[1][2]->Draw("colz");
  c1->SaveAs(Form("fcsTrkHcaldp.%s.png",ttag));

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdr[1][2]->Draw("colz");
  c1->cd(2); hydr[1][2]->Draw("colz");
  c1->cd(3); hpdr[1][2]->Draw("colz");
  c1->cd(4); hrdr[1][2]->Draw("colz");
  c1->SaveAs(Form("fcsTrkHcaldr.%s.png",ttag));

  // ── Track quality page: all 4 track types overlaid ──────────────────────
  c1->Clear();
  c1->Divide(2,3);

  // Macro to draw 4 types overlaid with auto-scale; repeated for each pad
  #define DRAWTYPES(arr) { \
    double _mx=0; \
    for(int _t=0;_t<mNType;_t++) if(arr[_t]) _mx=TMath::Max(_mx,arr[_t]->GetMaximum()); \
    bool _first=true; \
    for(int _t=0;_t<mNType;_t++){ if(!arr[_t]) continue; \
      arr[_t]->SetMaximum(_mx*1.15); arr[_t]->Draw(_first?"":"same"); _first=false; } }

  c1->cd(1); gPad->SetLogy(); DRAWTYPES(hTrkQPt);
  c1->cd(2); DRAWTYPES(hTrkEta);
  c1->cd(3); DRAWTYPES(hTrkPhi);
  c1->cd(4); gPad->SetLogy(); DRAWTYPES(hTrkDcaZ);
  #undef DRAWTYPES

  // Legend in pad 5+6 (pads 5 and 6 in 2×3 layout, bottom row)
  c1->cd(5);
  float lx=0.2, ly=0.8, ldy=0.12;
  for(int tt=0; tt<mNType; tt++){
    t = new TText(lx, ly-tt*ldy, Form("Type %d: %s",tt,TTYPE[tt]));
    t->SetNDC(); t->SetTextColor(TCOL[tt]); t->SetTextSize(0.08);
    t->Draw();
  }
  c1->SaveAs("fcsTrkQuality.png");

  // ── MC truth-matching resolution page (only filled for simulation) ─────
  // Check if histograms exist and have entries
  bool hasTruth = false;
  for(int tt=0; tt<mNType; tt++)
    if(hTrkDqPt[tt] && hTrkDqPt[tt]->GetEntries()>0) hasTruth=true;

  if(hasTruth){
    c1->Clear();
    c1->Divide(3,2);

    // DqPt = q/pT_reco - q/pT_true, one pad per track type (skip FwdVtx if empty)
    #define DRAWTYPES_TRUTH(arr,title) { \
      bool _first=true; double _mx=0; \
      for(int _t=0;_t<mNType;_t++) if(arr[_t]&&arr[_t]->GetEntries()>0) _mx=TMath::Max(_mx,arr[_t]->GetMaximum()); \
      for(int _t=0;_t<mNType;_t++){ if(!arr[_t]||arr[_t]->GetEntries()==0) continue; \
        arr[_t]->SetMaximum(_mx*1.2); arr[_t]->Draw(_first?"":"same"); _first=false; } \
      TLatex *_la=new TLatex(0.15,0.85,title); _la->SetNDC(); _la->SetTextSize(0.05); _la->Draw(); }

    c1->cd(1); DRAWTYPES_TRUTH(hTrkDqPt,   "q/p_{T}^{reco}-q/p_{T}^{true}");
    c1->cd(2); DRAWTYPES_TRUTH(hTrkDpTrel, "(p_{T}^{reco}-p_{T}^{true})/p_{T}^{true}");

    // 2D reco vs truth for Beamline (most useful for BLC resolution)
    c1->cd(3);
    if(hTrkQPtRecovsTru[1] && hTrkQPtRecovsTru[1]->GetEntries()>0){
      hTrkQPtRecovsTru[1]->Draw("colz");
      TLine *rDiagBLC=new TLine(-5,-5,5,5); rDiagBLC->SetLineColor(kRed); rDiagBLC->Draw();
    }

    // 2D reco vs truth for Primary
    c1->cd(4);
    if(hTrkQPtRecovsTru[2] && hTrkQPtRecovsTru[2]->GetEntries()>0){
      hTrkQPtRecovsTru[2]->Draw("colz");
      TLine *rDiagPri=new TLine(-5,-5,5,5); rDiagPri->SetLineColor(kRed); rDiagPri->Draw();
    }

    // Legend + resolution summary: fit ±0.3 Gaussian for core width
    c1->cd(5);
    TLatex *rLtx = new TLatex();
    rLtx->SetNDC(); rLtx->SetTextSize(0.055);
    rLtx->DrawLatex(0.05,0.92,"MC truth resolution (Gauss fit #pm0.3)");
    float rYpos=0.78;
    for(int tt=0; tt<mNType; tt++){
      if(!hTrkDqPt[tt]||hTrkDqPt[tt]->GetEntries()==0) continue;
      double mean = hTrkDqPt[tt]->GetMean();
      TF1 *gcore = new TF1(Form("gc%d",tt),"gaus", mean-0.3, mean+0.3);
      hTrkDqPt[tt]->Fit(gcore,"QRN");
      float sigRes  = gcore->GetParameter(2);
      float biasRes = gcore->GetParameter(1);
      rLtx->SetTextColor(TCOL[tt]);
      rLtx->DrawLatex(0.05, rYpos, Form("%s: #sigma_{core}=%.3f bias=%+.3f",
                                        TTYPE[tt], sigRes, biasRes));
      rYpos -= 0.13;
    }
    #undef DRAWTYPES_TRUTH

    c1->SaveAs("fcsTrkResolution.png");

    // ── Dedicated BLCVtx vs FCSTRK momentum resolution page ────────────────
    // Shows only types 4 (BLCVtx) and 5 (FCSTRK); fits double Gaussian.
    // Double Gaussian: narrow core (sigma1) + wide tail (sigma2), same center.
    // Core sigma1 is the "good track" resolution; tail captures secondaries/failures.
    // Range of fit: full histogram (tails important for shape), but report sigma1.
    TH1F *hDqBLC = hTrkDqPt[4], *hDqFCS = hTrkDqPt[5];
    if(hDqBLC && hDqBLC->GetEntries()>0 && hDqFCS && hDqFCS->GetEntries()>0){
      TCanvas *cMom = new TCanvas("cMomRes","",800,600);
      gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
      // normalize to unit area for shape comparison
      TH1F *hBnorm = (TH1F*)hDqBLC->Clone("hBnorm");
      TH1F *hFnorm = (TH1F*)hDqFCS->Clone("hFnorm");
      if(hBnorm->Integral()>0) hBnorm->Scale(1.0/hBnorm->Integral());
      if(hFnorm->Integral()>0) hFnorm->Scale(1.0/hFnorm->Integral());
      hBnorm->SetLineColor(kBlue);  hBnorm->SetLineWidth(2);
      hFnorm->SetLineColor(kRed);   hFnorm->SetLineWidth(2);
      double ymx = TMath::Max(hBnorm->GetMaximum(), hFnorm->GetMaximum());
      hBnorm->SetMaximum(ymx*1.4);
      hBnorm->SetTitle("q/p_{T} residual: BLCVtx vs FCSTRK;q/p_{T}^{reco}-q/p_{T}^{true} [(GeV/c)^{-1}];Norm. entries");
      hBnorm->Draw("hist"); hFnorm->Draw("hist same");
      // Double-Gaussian fit: [0]*G(x,mu,s1) + [3]*G(x,mu,s2), s1<s2 (CINT-compatible, no lambda)
      double s1B,s2B,f1B, s1F,s2F,f1F;
      { // BLCVtx fit
        TH1F *h=hBnorm; double mean=h->GetMean(), rms=h->GetRMS();
        TF1 *f=new TF1("dgB","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[1],[4])",
                       h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
        f->SetParameters(h->GetMaximum(),mean,0.5*rms,0.3*h->GetMaximum(),2.0*rms);
        f->SetParLimits(2,1e-4,rms); f->SetParLimits(4,rms*0.5,rms*10.);
        f->SetLineColor(kBlue); f->SetLineStyle(2); f->SetLineWidth(2);
        h->Fit(f,"RQN"); s1B=f->GetParameter(2); s2B=f->GetParameter(4);
        if(s1B>s2B){double t=s1B;s1B=s2B;s2B=t;}
        double a1=f->GetParameter(0),a2=f->GetParameter(3);
        f1B=(a1+a2>0)?a1*s1B/(a1*s1B+a2*s2B):0; f->Draw("same");
      }
      { // FCSTRK fit
        TH1F *h=hFnorm; double mean=h->GetMean(), rms=h->GetRMS();
        TF1 *f=new TF1("dgF","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[1],[4])",
                       h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
        f->SetParameters(h->GetMaximum(),mean,0.5*rms,0.3*h->GetMaximum(),2.0*rms);
        f->SetParLimits(2,1e-4,rms); f->SetParLimits(4,rms*0.5,rms*10.);
        f->SetLineColor(kRed); f->SetLineStyle(2); f->SetLineWidth(2);
        h->Fit(f,"RQN"); s1F=f->GetParameter(2); s2F=f->GetParameter(4);
        if(s1F>s2F){double t=s1F;s1F=s2F;s2F=t;}
        double a1=f->GetParameter(0),a2=f->GetParameter(3);
        f1F=(a1+a2>0)?a1*s1F/(a1*s1F+a2*s2F):0; f->Draw("same");
      }
      double improv = (s1B>0) ? 100.*(1.-s1F/s1B) : 0;
      TLegend *lMom = new TLegend(0.55,0.55,0.88,0.88);
      lMom->SetBorderSize(0);
      lMom->AddEntry(hBnorm,Form("BLCVtx   N=%.0f", hDqBLC->GetEntries()),"l");
      lMom->AddEntry((TObject*)0,Form("  #sigma_{core}=%.3f  #sigma_{tail}=%.3f",s1B,s2B),"");
      lMom->AddEntry(hFnorm,Form("FCSTRK  N=%.0f", hDqFCS->GetEntries()),"l");
      lMom->AddEntry((TObject*)0,Form("  #sigma_{core}=%.3f  #sigma_{tail}=%.3f (%.1f%% better)",s1F,s2F,improv),"");
      lMom->Draw();
      cMom->SaveAs("fcsTrkMomResolution.png");
      printf("Momentum resolution (double-Gaussian core):\n");
      printf("  BLCVtx: sigma_core=%.3f  sigma_tail=%.3f\n",s1B,s2B);
      printf("  FCSTRK: sigma_core=%.3f  sigma_tail=%.3f  improvement=%.1f%%\n",s1F,s2F,improv);
    }
  }

  // ── BLCVtx vertex quality page ──────────────────────────────────────────
  // Analogous to plotDilep.C panel 4: vertex z, avg, diff (ZVTX/ZVTXA/ZVTXD)
  //   Pad 1: BLC track DCA-z (individual, wide → single-track z resolution)
  //   Pad 2: BLCVtx track DCA-z after refit (narrow ~0 → improvement from vertex)
  //   Pad 3: per-event avg BLC DCA-z = estimated z_vtx (peaks at true vz)
  //   Pad 4: BLC dcaZ_i - dcaZ_j for all pairs (vertex resolution proxy, ~2σ wide)
  //   Pad 5: estimated z_vtx - event primary z (= MC truth vz in simulation)
  //   Pad 6: N BLC tracks per event
  TH1F *hBLCVtxN       = (TH1F*)F->Get("blcVtxN");
  TH1F *hBLCVtxZ       = (TH1F*)F->Get("blcVtxZ");
  TH1F *hBLCVtxZres    = (TH1F*)F->Get("blcVtxZres");
  TH1F *hBLCVtxZdiff   = (TH1F*)F->Get("blcVtxZdiff");
  TH1F *hBLCVtxTrkDcaZ = (TH1F*)F->Get("blcVtxTrkDcaZ");

  bool hasBLCVtx = (hBLCVtxN && hBLCVtxN->GetEntries()>0) ||
                   (hBLCVtxTrkDcaZ && hBLCVtxTrkDcaZ->GetEntries()>0);
  if(hasBLCVtx){
    c1->Clear();
    c1->Divide(2,3);

    c1->cd(1); gPad->SetLogy();
    if(hTrkDcaZ[1]){ hTrkDcaZ[1]->SetLineColor(kBlue); hTrkDcaZ[1]->Draw(); }
    TLatex *vl1=new TLatex(0.15,0.85,"BLC DCA-z (individual)"); vl1->SetNDC(); vl1->SetTextSize(0.055); vl1->Draw();

    c1->cd(2); gPad->SetLogy();
    if(hBLCVtxTrkDcaZ){ hBLCVtxTrkDcaZ->SetLineColor(kMagenta+1); hBLCVtxTrkDcaZ->Draw(); }
    TLatex *vl2=new TLatex(0.15,0.85,"BLCVtx DCA-z (after refit)"); vl2->SetNDC(); vl2->SetTextSize(0.055); vl2->Draw();

    c1->cd(3); gPad->SetLogy();
    if(hBLCVtxZ){ hBLCVtxZ->SetLineColor(kRed); hBLCVtxZ->Draw(); }
    TLatex *vl3=new TLatex(0.15,0.85,"avg BLC DCA-z = z_{vtx}"); vl3->SetNDC(); vl3->SetTextSize(0.055); vl3->Draw();

    c1->cd(4); gPad->SetLogy();
    if(hBLCVtxZdiff){ hBLCVtxZdiff->SetLineColor(kBlack); hBLCVtxZdiff->Draw(); }
    TLatex *vl4=new TLatex(0.15,0.85,"BLC dcaZ_{i}-dcaZ_{j} (pairs)"); vl4->SetNDC(); vl4->SetTextSize(0.055); vl4->Draw();

    c1->cd(5); gPad->SetLogy();
    if(hBLCVtxZres){ hBLCVtxZres->SetLineColor(kGreen+2); hBLCVtxZres->Draw(); }
    TLatex *vl5=new TLatex(0.15,0.85,"z_{vtx}^{BLC} - z_{vtx}^{event}"); vl5->SetNDC(); vl5->SetTextSize(0.055); vl5->Draw();

    c1->cd(6);
    if(hBLCVtxN){ hBLCVtxN->SetLineColor(kBlack); hBLCVtxN->Draw(); }
    TLatex *vl6=new TLatex(0.15,0.85,"N BLC tracks/event"); vl6->SetNDC(); vl6->SetTextSize(0.055); vl6->Draw();

    c1->SaveAs("fcsTrkBLCVtx.png");
  }

  // ── ECAL position resolution vs cluster energy ──────────────────────────
  // Requires histograms written by picoMatch.C with mTrackType=2 (Primary).
  // dx = fcsX - Primary_track_ecalProjection.X().
  // Primary track is constrained to vertex; sigma_track_proj << sigma_pos,
  // so sigma(dx) ~ sigma_pos(E) directly.
  TH2F *hFcsEdxVsE = (TH2F*)F->Get("FcsEdxVsE");
  TH2F *hFcsEdyVsE = (TH2F*)F->Get("FcsEdyVsE");
  if(hFcsEdxVsE && hFcsEdxVsE->GetEntries()>0){
    // Group 5 consecutive 1-GeV bins per fit to get ~50-150 entries per group
    // Avoids ROOT5 RebinX() in-place vs return-value ambiguity.
    int ne_orig = hFcsEdxVsE->GetNbinsX(); // 50 bins × 1 GeV
    int ngrp = 5;  // 5 GeV per group
    TGraphErrors *gSigDx = new TGraphErrors();
    TGraphErrors *gSigDy = new TGraphErrors();
    int npt=0;
    // storage for multi-panel slice plot (all bins with enough entries)
    TH1D *slDxArrE[12]; TH1D *slDyArrE[12];
    TF1  *fDxArrE[12];  TF1  *fDyArrE[12];
    double eCenArrE[12], sDxArrE[12], sDyArrE[12];
    int   passArrE[12], nSlE=0;
    for(int i=0;i<12;i++){slDxArrE[i]=0;slDyArrE[i]=0;fDxArrE[i]=0;fDyArrE[i]=0;passArrE[i]=0;}
    printf("\n=== ECAL Position Resolution (sigma(dx) and sigma(dy) vs E) ===\n");
    printf("  E_lo-hi  n_dx  sigma_dx  err_dx  sigma_dy  err_dy  [cm]\n");
    for(int igrp=0; igrp*ngrp < ne_orig; igrp++){
      int ib1 = igrp*ngrp + 1;
      int ib2 = ib1+ngrp-1; if(ib2>ne_orig) ib2=ne_orig;
      double eCenter = 0.5*(hFcsEdxVsE->GetXaxis()->GetBinLowEdge(ib1) +
                            hFcsEdxVsE->GetXaxis()->GetBinUpEdge(ib2));
      TH1D *slDx = hFcsEdxVsE->ProjectionY(Form("slDxG%d",igrp), ib1, ib2);
      TH1D *slDy = hFcsEdyVsE->ProjectionY(Form("slDyG%d",igrp), ib1, ib2);
      if(slDx->GetEntries()<10) continue;
      const double fitWin = 5.0;
      TF1 *fDx = new TF1(Form("fEdxG%d",igrp),"gaus",-20,20);
      double peakDx = slDx->GetBinCenter(slDx->GetMaximumBin());
      fDx->SetParameters(slDx->GetMaximum(), peakDx, 2.0);
      fDx->SetParLimits(2, 0.05, 12.0);
      fDx->SetRange(peakDx-fitWin, peakDx+fitWin);
      slDx->Fit(fDx,"QRBN");
      TF1 *fDy = new TF1(Form("fEdyG%d",igrp),"gaus",-20,20);
      double peakDy = slDy->GetBinCenter(slDy->GetMaximumBin());
      fDy->SetParameters(slDy->GetMaximum(), peakDy, 2.0);
      fDy->SetParLimits(2, 0.05, 12.0);
      fDy->SetRange(peakDy-fitWin, peakDy+fitWin);
      slDy->Fit(fDy,"QRBN");
      double sDx = fDx->GetParameter(2), eDx = fDx->GetParError(2);
      double sDy = fDy->GetParameter(2), eDy = fDy->GetParError(2);
      int pass = (sDx>=0.05 && sDx<12 && sDy>=0.05 && sDy<12);
      if(nSlE<12){ slDxArrE[nSlE]=slDx; slDyArrE[nSlE]=slDy;
                   fDxArrE[nSlE]=fDx;   fDyArrE[nSlE]=fDy;
                   eCenArrE[nSlE]=eCenter; sDxArrE[nSlE]=sDx; sDyArrE[nSlE]=sDy;
                   passArrE[nSlE]=pass; nSlE++; }
      if(!pass) continue;
      printf("  %6.1f  %5.0f   %7.4f  %6.4f  %7.4f  %6.4f\n",
             eCenter,(double)slDx->GetEntries(),sDx,eDx,sDy,eDy);
      gSigDx->SetPoint(npt, eCenter, sDx); gSigDx->SetPointError(npt, 0, eDx);
      gSigDy->SetPoint(npt, eCenter, sDy); gSigDy->SetPointError(npt, 0, eDy);
      npt++;
    }

    TF1 *fResX = new TF1("fResX","[0]/sqrt(x)+[1]", 1.0, 50.0);
    fResX->SetParameters(0.5, 0.2); fResX->SetLineColor(kBlue); fResX->SetLineWidth(2);
    TF1 *fResY = new TF1("fResY","[0]/sqrt(x)+[1]", 1.0, 50.0);
    fResY->SetParameters(0.5, 0.2); fResY->SetLineColor(kRed);  fResY->SetLineWidth(2);

    gSigDx->SetMarkerStyle(20); gSigDx->SetMarkerColor(kBlue); gSigDx->SetLineColor(kBlue);
    gSigDy->SetMarkerStyle(24); gSigDy->SetMarkerColor(kRed);  gSigDy->SetLineColor(kRed);

    c1->Clear(); c1->Divide(2,2);
    c1->cd(1); hFcsEdxVsE->Draw("colz");
    c1->cd(2); hFcsEdyVsE->Draw("colz");

    c1->cd(3);
    gSigDx->SetTitle("#sigma(dx) vs E_{cluster} (ECAL);E_{cluster} [GeV];#sigma(dx) [cm]");
    gSigDx->Draw("AP");
    if(npt>=3){ gSigDx->Fit("fResX","R"); fResX->Draw("same"); }
    TLatex *laDx = new TLatex(); laDx->SetNDC(); laDx->SetTextSize(0.048);
    laDx->DrawLatex(0.18,0.82, Form("dx: A=%.3f B=%.3f cm",
        fResX->GetParameter(0), fResX->GetParameter(1)));
    // draw theory curve for comparison
    TF1 *fTheory = new TF1("fTheory","0.5/sqrt(x)+0.2", 1.0, 50.0);
    fTheory->SetLineColor(kGray+2); fTheory->SetLineStyle(2); fTheory->Draw("same");
    laDx->DrawLatex(0.18,0.75,"-- theory: 0.50/#sqrt{E}+0.20");

    c1->cd(4);
    gSigDy->SetTitle("#sigma(dy) vs E_{cluster} (ECAL);E_{cluster} [GeV];#sigma(dy) [cm]");
    gSigDy->Draw("AP");
    if(npt>=3){ gSigDy->Fit("fResY","R"); fResY->Draw("same"); }
    TLatex *laDy = new TLatex(); laDy->SetNDC(); laDy->SetTextSize(0.048);
    laDy->DrawLatex(0.18,0.82, Form("dy: A=%.3f B=%.3f cm",
        fResY->GetParameter(0), fResY->GetParameter(1)));
    TF1 *fTheoryY = (TF1*)fTheory->Clone("fTheoryY");
    fTheoryY->Draw("same");
    laDy->DrawLatex(0.18,0.75,"-- theory: 0.50/#sqrt{E}+0.20");

    printf("  dx fit: A=%.3f+/-%.3f cm*sqrt(GeV), B=%.3f+/-%.3f cm\n",
        fResX->GetParameter(0),fResX->GetParError(0),
        fResX->GetParameter(1),fResX->GetParError(1));
    printf("  dy fit: A=%.3f+/-%.3f cm*sqrt(GeV), B=%.3f+/-%.3f cm\n",
        fResY->GetParameter(0),fResY->GetParError(0),
        fResY->GetParameter(1),fResY->GetParError(1));

    c1->SaveAs("fcsTrkFcsResolution.png");

    // Multi-panel: one column per energy bin, top row=dx, bottom row=dy
    if(nSlE>0){
      TCanvas *cSlE = new TCanvas("cFcsSlice","",270*nSlE,540);
      cSlE->Divide(nSlE,2);
      TLatex *lxE=new TLatex(); lxE->SetNDC(); lxE->SetTextSize(0.08);
      TLatex *lyE=new TLatex(); lyE->SetNDC(); lyE->SetTextSize(0.08);
      for(int ip=0; ip<nSlE; ip++){
        double cxE, rxE, sxE;
        // dx
        cSlE->cd(ip+1);
        gPad->SetLeftMargin(0.18); gPad->SetBottomMargin(0.18);
        slDxArrE[ip]->SetTitle(Form("E=%.0f GeV;dx [cm];",eCenArrE[ip]));
        slDxArrE[ip]->SetLineColor(kBlue);
        cxE=fDxArrE[ip]->GetParameter(1); sxE=sDxArrE[ip];
        rxE=sxE*6>3.0?sxE*6:3.0; if(rxE>18)rxE=18;
        slDxArrE[ip]->GetXaxis()->SetRangeUser(cxE-rxE,cxE+rxE);
        slDxArrE[ip]->Draw();
        fDxArrE[ip]->SetLineColor(passArrE[ip]?kBlue:kRed);
        fDxArrE[ip]->SetLineWidth(2); fDxArrE[ip]->Draw("same");
        lxE->SetTextColor(passArrE[ip]?kBlue:kRed);
        lxE->DrawLatex(0.22,0.85,Form("#sigma_{x}=%.2f cm",sDxArrE[ip]));
        // dy
        cSlE->cd(ip+1+nSlE);
        gPad->SetLeftMargin(0.18); gPad->SetBottomMargin(0.18);
        slDyArrE[ip]->SetTitle(Form("E=%.0f GeV;dy [cm];",eCenArrE[ip]));
        slDyArrE[ip]->SetLineColor(kRed);
        cxE=fDyArrE[ip]->GetParameter(1); sxE=sDyArrE[ip];
        rxE=sxE*6>3.0?sxE*6:3.0; if(rxE>18)rxE=18;
        slDyArrE[ip]->GetXaxis()->SetRangeUser(cxE-rxE,cxE+rxE);
        slDyArrE[ip]->Draw();
        fDyArrE[ip]->SetLineColor(passArrE[ip]?kRed:kOrange+1);
        fDyArrE[ip]->SetLineWidth(2); fDyArrE[ip]->Draw("same");
        lyE->SetTextColor(passArrE[ip]?kRed:kOrange+1);
        lyE->DrawLatex(0.22,0.85,Form("#sigma_{y}=%.2f cm",sDyArrE[ip]));
      }
      cSlE->SaveAs("fcsTrkFcsResolutionSlices.png");
    }
  }

  // ── ECAL cluster vs MC truth impact (isolates sigma_pos, no track error) ──
  TH2F *hFcsMcDxVsE = (TH2F*)F->Get("FcsMcDxVsE");
  TH2F *hFcsMcDyVsE = (TH2F*)F->Get("FcsMcDyVsE");
  if(hFcsMcDxVsE && hFcsMcDxVsE->GetEntries()>0){
    int ne_mc = hFcsMcDxVsE->GetNbinsX();
    int ngrp_mc = 5;
    TGraphErrors *gMcDx = new TGraphErrors();
    TGraphErrors *gMcDy = new TGraphErrors();
    int nptMc=0;
    printf("\n=== ECAL Position Resolution: MC truth impact vs cluster (sigma_pos alone) ===\n");
    printf("  E_lo-hi  n_dx  sigma_dx  err_dx  sigma_dy  err_dy  [cm]\n");
    for(int igrp=0; igrp*ngrp_mc < ne_mc; igrp++){
      int ib1 = igrp*ngrp_mc + 1;
      int ib2 = ib1+ngrp_mc-1; if(ib2>ne_mc) ib2=ne_mc;
      double eCenter = 0.5*(hFcsMcDxVsE->GetXaxis()->GetBinLowEdge(ib1) +
                            hFcsMcDxVsE->GetXaxis()->GetBinUpEdge(ib2));
      TH1D *slDx = hFcsMcDxVsE->ProjectionY(Form("slMcDxG%d",igrp), ib1, ib2);
      TH1D *slDy = hFcsMcDyVsE->ProjectionY(Form("slMcDyG%d",igrp), ib1, ib2);
      if(slDx->GetEntries()<5) continue;
      const double fitWinMc = 4.0;
      TF1 *fDx = new TF1(Form("fMcDxG%d",igrp),"gaus",-20,20);
      double peakMcDx = slDx->GetBinCenter(slDx->GetMaximumBin());
      fDx->SetParameters(slDx->GetMaximum(), peakMcDx, 0.5);
      fDx->SetParLimits(2, 0.05, 4.0);
      fDx->SetRange(peakMcDx-fitWinMc, peakMcDx+fitWinMc);
      slDx->Fit(fDx,"QRBN");
      TF1 *fDy = new TF1(Form("fMcDyG%d",igrp),"gaus",-20,20);
      double peakMcDy = slDy->GetBinCenter(slDy->GetMaximumBin());
      fDy->SetParameters(slDy->GetMaximum(), peakMcDy, 0.5);
      fDy->SetParLimits(2, 0.05, 4.0);
      fDy->SetRange(peakMcDy-fitWinMc, peakMcDy+fitWinMc);
      slDy->Fit(fDy,"QRBN");
      double sDx = fDx->GetParameter(2), eDx = fDx->GetParError(2);
      double sDy = fDy->GetParameter(2), eDy = fDy->GetParError(2);
      if(sDx<0.05 || sDx>4 || sDy<0.05 || sDy>4) continue;
      printf("  %6.1f  %5.0f   %7.4f  %6.4f  %7.4f  %6.4f\n",
             eCenter,(double)slDx->GetEntries(),sDx,eDx,sDy,eDy);
      gMcDx->SetPoint(nptMc, eCenter, sDx); gMcDx->SetPointError(nptMc, 0, eDx);
      gMcDy->SetPoint(nptMc, eCenter, sDy); gMcDy->SetPointError(nptMc, 0, eDy);
      nptMc++;
    }

    TF1 *fMcResX = new TF1("fMcResX","[0]/sqrt(x)+[1]", 1.0, 50.0);
    fMcResX->SetParameters(0.5,0.2); fMcResX->SetLineColor(kBlue); fMcResX->SetLineWidth(2);
    TF1 *fMcResY = new TF1("fMcResY","[0]/sqrt(x)+[1]", 1.0, 50.0);
    fMcResY->SetParameters(0.5,0.2); fMcResY->SetLineColor(kRed);  fMcResY->SetLineWidth(2);
    gMcDx->SetMarkerStyle(20); gMcDx->SetMarkerColor(kBlue); gMcDx->SetLineColor(kBlue);
    gMcDy->SetMarkerStyle(24); gMcDy->SetMarkerColor(kRed);  gMcDy->SetLineColor(kRed);

    c1->Clear(); c1->Divide(2,2);
    c1->cd(1); hFcsMcDxVsE->Draw("colz");
    c1->cd(2); hFcsMcDyVsE->Draw("colz");
    c1->cd(3);
    gMcDx->SetTitle("#sigma(dx) vs E (MC truth);E_{cluster} [GeV];#sigma(dx) [cm]");
    gMcDx->Draw("AP");
    if(nptMc>=2){ gMcDx->Fit("fMcResX","R"); fMcResX->Draw("same"); }
    TLatex *lMcX = new TLatex(); lMcX->SetNDC(); lMcX->SetTextSize(0.048);
    lMcX->DrawLatex(0.18,0.82, Form("dx: A=%.3f B=%.3f cm",
        fMcResX->GetParameter(0), fMcResX->GetParameter(1)));
    TF1 *fThMcX = new TF1("fThMcX","0.5/sqrt(x)+0.2", 1.0, 50.0);
    fThMcX->SetLineColor(kGray+2); fThMcX->SetLineStyle(2); fThMcX->Draw("same");
    lMcX->DrawLatex(0.18,0.75,"-- theory: 0.50/#sqrt{E}+0.20");
    c1->cd(4);
    gMcDy->SetTitle("#sigma(dy) vs E (MC truth);E_{cluster} [GeV];#sigma(dy) [cm]");
    gMcDy->Draw("AP");
    if(nptMc>=2){ gMcDy->Fit("fMcResY","R"); fMcResY->Draw("same"); }
    TLatex *lMcY = new TLatex(); lMcY->SetNDC(); lMcY->SetTextSize(0.048);
    lMcY->DrawLatex(0.18,0.82, Form("dy: A=%.3f B=%.3f cm",
        fMcResY->GetParameter(0), fMcResY->GetParameter(1)));
    TF1 *fThMcY = (TF1*)fThMcX->Clone("fThMcY"); fThMcY->Draw("same");
    lMcY->DrawLatex(0.18,0.75,"-- theory: 0.50/#sqrt{E}+0.20");

    printf("  MC dx fit: A=%.3f+/-%.3f B=%.3f+/-%.3f cm\n",
        fMcResX->GetParameter(0),fMcResX->GetParError(0),
        fMcResX->GetParameter(1),fMcResX->GetParError(1));
    printf("  MC dy fit: A=%.3f+/-%.3f B=%.3f+/-%.3f cm\n",
        fMcResY->GetParameter(0),fMcResY->GetParError(0),
        fMcResY->GetParameter(1),fMcResY->GetParError(1));

    c1->SaveAs("fcsTrkMcResolution.png");
  }

  // ── ECAL cluster vs gamma MC truth (straight-line: exact sigma_pos) ──
  TH2F *hFcsGamDxVsE = (TH2F*)F->Get("FcsGamDxVsE");
  TH2F *hFcsGamDyVsE = (TH2F*)F->Get("FcsGamDyVsE");
  if(hFcsGamDxVsE && hFcsGamDxVsE->GetEntries()>0){
    int ne_gam = hFcsGamDxVsE->GetNbinsX();
    int ngrp_gam = 5;
    TGraphErrors *gGamDx = new TGraphErrors();
    TGraphErrors *gGamDy = new TGraphErrors();
    int nptGam=0;
    // storage for multi-panel slice plot
    TH1D *slDxArrG[12]; TH1D *slDyArrG[12];
    TF1  *fDxArrG[12];  TF1  *fDyArrG[12];
    double eCenArrG[12], sDxArrG[12], sDyArrG[12];
    int   passArrG[12], nSlG=0;
    for(int i=0;i<12;i++){slDxArrG[i]=0;slDyArrG[i]=0;fDxArrG[i]=0;fDyArrG[i]=0;passArrG[i]=0;}
    printf("\n=== ECAL Position Resolution: gamma truth (straight-line, exact sigma_pos) ===\n");
    printf("  E_lo-hi  n_dx  sigma_dx  err_dx  sigma_dy  err_dy  [cm]\n");
    for(int igrp=0; igrp*ngrp_gam < ne_gam; igrp++){
      int ib1 = igrp*ngrp_gam + 1;
      int ib2 = ib1+ngrp_gam-1; if(ib2>ne_gam) ib2=ne_gam;
      double eCenter = 0.5*(hFcsGamDxVsE->GetXaxis()->GetBinLowEdge(ib1) +
                            hFcsGamDxVsE->GetXaxis()->GetBinUpEdge(ib2));
      TH1D *slDx = hFcsGamDxVsE->ProjectionY(Form("slGamDxG%d",igrp), ib1, ib2);
      TH1D *slDy = hFcsGamDyVsE->ProjectionY(Form("slGamDyG%d",igrp), ib1, ib2);
      if(slDx->GetEntries()<5) continue;
      const double fitWinGam = 4.0;
      TF1 *fDx = new TF1(Form("fGamDxG%d",igrp),"gaus",-20,20);
      double peakGamDx = slDx->GetBinCenter(slDx->GetMaximumBin());
      fDx->SetParameters(slDx->GetMaximum(), peakGamDx, 0.5);
      fDx->SetParLimits(2, 0.05, 4.0);
      fDx->SetRange(peakGamDx-fitWinGam, peakGamDx+fitWinGam);
      slDx->Fit(fDx,"QRBN");
      TF1 *fDy = new TF1(Form("fGamDyG%d",igrp),"gaus",-20,20);
      double peakGamDy = slDy->GetBinCenter(slDy->GetMaximumBin());
      fDy->SetParameters(slDy->GetMaximum(), peakGamDy, 0.5);
      fDy->SetParLimits(2, 0.05, 4.0);
      fDy->SetRange(peakGamDy-fitWinGam, peakGamDy+fitWinGam);
      slDy->Fit(fDy,"QRBN");
      double sDx = fDx->GetParameter(2), eDx = fDx->GetParError(2);
      double sDy = fDy->GetParameter(2), eDy = fDy->GetParError(2);
      int pass = (sDx>=0.05 && sDx<4 && sDy>=0.05 && sDy<4);
      if(nSlG<12){ slDxArrG[nSlG]=slDx; slDyArrG[nSlG]=slDy;
                   fDxArrG[nSlG]=fDx;   fDyArrG[nSlG]=fDy;
                   eCenArrG[nSlG]=eCenter; sDxArrG[nSlG]=sDx; sDyArrG[nSlG]=sDy;
                   passArrG[nSlG]=pass; nSlG++; }
      if(!pass) continue;
      printf("  %6.1f  %5.0f   %7.4f  %6.4f  %7.4f  %6.4f\n",
             eCenter,(double)slDx->GetEntries(),sDx,eDx,sDy,eDy);
      gGamDx->SetPoint(nptGam, eCenter, sDx); gGamDx->SetPointError(nptGam, 0, eDx);
      gGamDy->SetPoint(nptGam, eCenter, sDy); gGamDy->SetPointError(nptGam, 0, eDy);
      nptGam++;
    }

    TF1 *fGamResX = new TF1("fGamResX","[0]/sqrt(x)+[1]", 1.0, 50.0);
    fGamResX->SetParameters(0.5,0.2); fGamResX->SetLineColor(kBlue); fGamResX->SetLineWidth(2);
    TF1 *fGamResY = new TF1("fGamResY","[0]/sqrt(x)+[1]", 1.0, 50.0);
    fGamResY->SetParameters(0.5,0.2); fGamResY->SetLineColor(kRed);  fGamResY->SetLineWidth(2);
    gGamDx->SetMarkerStyle(20); gGamDx->SetMarkerColor(kBlue); gGamDx->SetLineColor(kBlue);
    gGamDy->SetMarkerStyle(24); gGamDy->SetMarkerColor(kRed);  gGamDy->SetLineColor(kRed);

    c1->Clear(); c1->Divide(2,2);
    c1->cd(1); hFcsGamDxVsE->Draw("colz");
    c1->cd(2); hFcsGamDyVsE->Draw("colz");
    c1->cd(3);
    gGamDx->SetTitle("#sigma(dx) vs E (#gamma truth);E_{cluster} [GeV];#sigma(dx) [cm]");
    gGamDx->Draw("AP");
    if(nptGam>=2){ gGamDx->Fit("fGamResX","R"); fGamResX->Draw("same"); }
    TLatex *lGamX = new TLatex(); lGamX->SetNDC(); lGamX->SetTextSize(0.048);
    lGamX->DrawLatex(0.18,0.82, Form("dx: A=%.3f B=%.3f cm",
        fGamResX->GetParameter(0), fGamResX->GetParameter(1)));
    TF1 *fThGamX = new TF1("fThGamX","0.5/sqrt(x)+0.2", 1.0, 50.0);
    fThGamX->SetLineColor(kGray+2); fThGamX->SetLineStyle(2); fThGamX->Draw("same");
    lGamX->DrawLatex(0.18,0.75,"-- theory: 0.50/#sqrt{E}+0.20");
    c1->cd(4);
    gGamDy->SetTitle("#sigma(dy) vs E (#gamma truth);E_{cluster} [GeV];#sigma(dy) [cm]");
    gGamDy->Draw("AP");
    if(nptGam>=2){ gGamDy->Fit("fGamResY","R"); fGamResY->Draw("same"); }
    TLatex *lGamY = new TLatex(); lGamY->SetNDC(); lGamY->SetTextSize(0.048);
    lGamY->DrawLatex(0.18,0.82, Form("dy: A=%.3f B=%.3f cm",
        fGamResY->GetParameter(0), fGamResY->GetParameter(1)));
    TF1 *fThGamY = (TF1*)fThGamX->Clone("fThGamY"); fThGamY->Draw("same");
    lGamY->DrawLatex(0.18,0.75,"-- theory: 0.50/#sqrt{E}+0.20");

    printf("  gamma dx fit: A=%.3f+/-%.3f B=%.3f+/-%.3f cm\n",
        fGamResX->GetParameter(0),fGamResX->GetParError(0),
        fGamResX->GetParameter(1),fGamResX->GetParError(1));
    printf("  gamma dy fit: A=%.3f+/-%.3f B=%.3f+/-%.3f cm\n",
        fGamResY->GetParameter(0),fGamResY->GetParError(0),
        fGamResY->GetParameter(1),fGamResY->GetParError(1));

    c1->SaveAs("fcsTrkGamResolution.png");

    // Multi-panel: one column per energy bin, top row=dx, bottom row=dy
    if(nSlG>0){
      TCanvas *cSlG = new TCanvas("cGamSlice","",270*nSlG,540);
      cSlG->Divide(nSlG,2);
      TLatex *lxG=new TLatex(); lxG->SetNDC(); lxG->SetTextSize(0.08);
      TLatex *lyG=new TLatex(); lyG->SetNDC(); lyG->SetTextSize(0.08);
      for(int ip=0; ip<nSlG; ip++){
        double cxG, rxG, sxG;
        // dx
        cSlG->cd(ip+1);
        gPad->SetLeftMargin(0.18); gPad->SetBottomMargin(0.18);
        slDxArrG[ip]->SetTitle(Form("E=%.0f GeV;dx [cm];",eCenArrG[ip]));
        slDxArrG[ip]->SetLineColor(kBlue);
        cxG=fDxArrG[ip]->GetParameter(1); sxG=sDxArrG[ip];
        rxG=sxG*6>1.5?sxG*6:1.5; if(rxG>10)rxG=10;
        slDxArrG[ip]->GetXaxis()->SetRangeUser(cxG-rxG,cxG+rxG);
        slDxArrG[ip]->Draw();
        fDxArrG[ip]->SetLineColor(passArrG[ip]?kBlue:kRed);
        fDxArrG[ip]->SetLineWidth(2); fDxArrG[ip]->Draw("same");
        lxG->SetTextColor(passArrG[ip]?kBlue:kRed);
        lxG->DrawLatex(0.22,0.85,Form("#sigma_{x}=%.2f cm",sDxArrG[ip]));
        // dy
        cSlG->cd(ip+1+nSlG);
        gPad->SetLeftMargin(0.18); gPad->SetBottomMargin(0.18);
        slDyArrG[ip]->SetTitle(Form("E=%.0f GeV;dy [cm];",eCenArrG[ip]));
        slDyArrG[ip]->SetLineColor(kRed);
        cxG=fDyArrG[ip]->GetParameter(1); sxG=sDyArrG[ip];
        rxG=sxG*6>1.5?sxG*6:1.5; if(rxG>10)rxG=10;
        slDyArrG[ip]->GetXaxis()->SetRangeUser(cxG-rxG,cxG+rxG);
        slDyArrG[ip]->Draw();
        fDyArrG[ip]->SetLineColor(passArrG[ip]?kRed:kOrange+1);
        fDyArrG[ip]->SetLineWidth(2); fDyArrG[ip]->Draw("same");
        lyG->SetTextColor(passArrG[ip]?kRed:kOrange+1);
        lyG->DrawLatex(0.22,0.85,Form("#sigma_{y}=%.2f cm",sDyArrG[ip]));
      }
      cSlG->SaveAs("fcsTrkGamResolutionSlices.png");
    }
  }

  // ── FCSTRK efficiency page ──────────────────────────────────────────────
  // Denominator: MC electrons with ECAL cluster match (dr < 15 cm)
  // Numerator A: above + BLCVtx (type=4) match
  // Numerator B: above + FCSTRK (type=5) match
  TH1F *hFCSEffDen    = (TH1F*)F->Get("fcsEffDen");
  TH1F *hFCSEffNum    = (TH1F*)F->Get("fcsEffNum");
  TH1F *hBLCVtxEffNum = (TH1F*)F->Get("blcVtxEffNum");
  if(hFCSEffDen && hFCSEffNum && hBLCVtxEffNum && hFCSEffDen->GetEntries()>0){
    TH1F *hBLCEff  = (TH1F*)hBLCVtxEffNum->Clone("hBLCEff");
    TH1F *hFCSEff  = (TH1F*)hFCSEffNum->Clone("hFCSEff");
    hBLCEff->Divide(hFCSEffDen);
    hFCSEff->Divide(hFCSEffDen);

    // Print summary numbers (integrated)
    double den  = hFCSEffDen->GetEntries();
    double nBLC = hBLCVtxEffNum->GetEntries();
    double nFCS = hFCSEffNum->GetEntries();
    printf("FCSTRK efficiency (MC e with ECAL cluster as denominator):\n");
    printf("  Denominator (MC e + ECAL cluster): %.0f\n", den);
    printf("  BLCVtx match:  %.0f / %.0f = %.1f%%\n", nBLC, den, den>0?100.*nBLC/den:0);
    printf("  FCSTRK match:  %.0f / %.0f = %.1f%%\n", nFCS, den, den>0?100.*nFCS/den:0);

    TCanvas *cEff = new TCanvas("cFCSEff","FCSTRK Efficiency",800,600);
    hBLCEff->SetLineColor(kBlue); hBLCEff->SetMarkerColor(kBlue); hBLCEff->SetMarkerStyle(20);
    hFCSEff->SetLineColor(kRed);  hFCSEff->SetMarkerColor(kRed);  hFCSEff->SetMarkerStyle(21);
    hBLCEff->GetYaxis()->SetRangeUser(0,1.2);
    hBLCEff->SetTitle("FCSTRK Efficiency (denom = MC e^{-} + ECAL cluster); p_{MC} [GeV/c]; Efficiency");
    hBLCEff->Draw("E"); hFCSEff->Draw("E same");
    TLegend *leg = new TLegend(0.6,0.2,0.88,0.45);
    leg->AddEntry(hBLCEff,Form("BLCVtx (%.1f%%)",den>0?100.*nBLC/den:0),"lp");
    leg->AddEntry(hFCSEff,Form("FCSTRK (%.1f%%)",den>0?100.*nFCS/den:0),"lp");
    leg->Draw();
    cEff->SaveAs("fcsTrkEfficiency.png");
  }
}
