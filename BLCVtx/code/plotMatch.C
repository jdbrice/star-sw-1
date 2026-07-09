#include <iostream>

static const int mNType=5;
const char *TTYPE[mNType]={"Global","Beamline","Primary","FwdVtx","BLCVtx"};
// Colors per track type: black, blue, red, green, magenta
static const int TCOL[mNType]={kBlack, kBlue, kRed, kGreen+2, kMagenta+1};

static const int mNCut=3;
const char *EH[2]={"Ecal","Hcal"};
const char *NSTB[4][2]={{"North","South"},{"Top","Bottom"},{"North","South"},{"R<70","R>70"}};
const char *CUT[4][mNCut+1]={{"Same Event","Mixed Event","Same-Mixed North","Same-Mixed South"},
                             {"Same Event","Mixed Event","Same-Mixed Top","Same-Mixed Bottom"},
                             {"Same Event","Mixed Event","Same-Mixed North","Same-Mixed South"},
                             {"Same Event","Mixed Event","Same-Mixed R<70","Same-Mixed R>70"}};

void plotMatch(char* data="202604",int run=0, int set=0){
  //void plotMatch(char* data=".",int run=1,int set=0){
  char file[100];
  if(run==0){
    sprintf(file,"%s/hist_match/all.match.root",data);
  }else{
    sprintf(file,"%s/hist_match/%d.%d.match.root",data,run,set);
  }
  printf("Reading %s\n",file);
  TFile *F = new TFile(file,"old");
  
  TH1F *hTrkQPt[mNType], *hTrkEta[mNType];
  TH1F *hTrkPhi[mNType], *hTrkDcaZ[mNType];
  // MC truth-matching histograms (only present in simulation match files)
  TH1F *hTrkDqPt[mNType], *hTrkDpTrel[mNType];
  TH2F *hTrkQPtRecovsTru[mNType];

  TH2F *hXY[4];
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

  hXY[0] = (TH2F*)F->Get("xyEcal");
  hXY[1] = (TH2F*)F->Get("xyHcal");
  hXY[2] = (TH2F*)F->Get("xyETrk");
  hXY[3] = (TH2F*)F->Get("xyHTrk");
  hE[2]  = (TH1F*)F->Get("ETrk");
  hET[2] = (TH1F*)F->Get("PTTrk");
  for(int eh=0; eh<2; eh++){
    hE[eh]  = (TH1F*)F->Get(Form("E%s",EH[eh]));
    hET[eh] = (TH1F*)F->Get(Form("ET%s",EH[eh]));
    for(int cut=0; cut<mNCut-1; cut++){
      hxdx[eh][cut] = (TH2F*)F->Get(Form("xdx%s%s",EH[eh],CUT[0][cut]));
      hydx[eh][cut] = (TH2F*)F->Get(Form("ydx%s%s",EH[eh],CUT[0][cut]));
      hpdx[eh][cut] = (TH2F*)F->Get(Form("pdx%s%s",EH[eh],CUT[0][cut]));
      hrdx[eh][cut] = (TH2F*)F->Get(Form("rdx%s%s",EH[eh],CUT[0][cut]));
      hxdy[eh][cut] = (TH2F*)F->Get(Form("xdy%s%s",EH[eh],CUT[1][cut]));
      hydy[eh][cut] = (TH2F*)F->Get(Form("ydy%s%s",EH[eh],CUT[1][cut]));
      hpdy[eh][cut] = (TH2F*)F->Get(Form("pdy%s%s",EH[eh],CUT[1][cut]));
      hrdy[eh][cut] = (TH2F*)F->Get(Form("rdy%s%s",EH[eh],CUT[1][cut]));
      hxdp[eh][cut] = (TH2F*)F->Get(Form("xdp%s%s",EH[eh],CUT[2][cut]));
      hydp[eh][cut] = (TH2F*)F->Get(Form("ydp%s%s",EH[eh],CUT[2][cut]));
      hpdp[eh][cut] = (TH2F*)F->Get(Form("pdp%s%s",EH[eh],CUT[2][cut]));
      hrdp[eh][cut] = (TH2F*)F->Get(Form("rdp%s%s",EH[eh],CUT[2][cut]));
      hxdr[eh][cut] = (TH2F*)F->Get(Form("xdr%s%s",EH[eh],CUT[3][cut]));
      hydr[eh][cut] = (TH2F*)F->Get(Form("ydr%s%s",EH[eh],CUT[3][cut]));
      hpdr[eh][cut] = (TH2F*)F->Get(Form("pdr%s%s",EH[eh],CUT[3][cut]));
      hrdr[eh][cut] = (TH2F*)F->Get(Form("rdr%s%s",EH[eh],CUT[3][cut]));
      for(int nstb=0; nstb<2; nstb++){
        hdx[eh][nstb][cut] = (TH1F*)F->Get(Form("%s%sdx%s",EH[eh],NSTB[0][nstb],CUT[0][cut]));
        hdy[eh][nstb][cut] = (TH1F*)F->Get(Form("%s%sdy%s",EH[eh],NSTB[1][nstb],CUT[1][cut]));
        hdp[eh][nstb][cut] = (TH1F*)F->Get(Form("%s%sdp%s",EH[eh],NSTB[2][nstb],CUT[2][cut]));
        hdr[eh][nstb][cut] = (TH1F*)F->Get(Form("%s%sdr%s",EH[eh],NSTB[3][nstb],CUT[3][cut]));
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
  c1->SaveAs("fcsTrkMatchEcal.png");

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
  c1->SaveAs("fcsTrkMatchHcal.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hXY[0]->Draw("colz");
  c1->cd(2); hXY[1]->Draw("colz");
  c1->cd(3); hXY[2]->Draw("colz");
  c1->cd(4); hXY[3]->Draw("colz");
  c1->SaveAs("fcsTrkxy.png");

  c1->Clear();
  c1->Divide(3,2);
  c1->cd(1); hET[0]->Draw();
  c1->cd(2); hET[1]->Draw();
  c1->cd(3); hET[2]->Draw();
  c1->cd(4); hE[0]->Draw();
  c1->cd(5); hE[1]->Draw();
  c1->cd(6); hE[2]->Draw();
  c1->SaveAs("fcsTrkEt.png");

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
  c1->SaveAs("fcsTrkEcaldx.png");
  c1->Clear();

  c1->Divide(2,2);
  c1->cd(1); hxdy[0][2]->Draw("colz");
  c1->cd(2); hydy[0][2]->Draw("colz");
  c1->cd(3); hpdy[0][2]->Draw("colz");
  c1->cd(4); hrdy[0][2]->Draw("colz");
  c1->SaveAs("fcsTrkEcaldy.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdp[0][2]->Draw("colz");
  c1->cd(2); hydp[0][2]->Draw("colz");
  c1->cd(3); hpdp[0][2]->Draw("colz");
  c1->cd(4); hrdp[0][2]->Draw("colz");
  c1->SaveAs("fcsTrkEcaldp.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdr[0][2]->Draw("colz");
  c1->cd(2); hydr[0][2]->Draw("colz");
  c1->cd(3); hpdr[0][2]->Draw("colz");
  c1->cd(4); hrdr[0][2]->Draw("colz");
  c1->SaveAs("fcsTrkEcaldr.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdx[1][2]->Draw("colz");
  c1->cd(2); hydx[1][2]->Draw("colz");
  c1->cd(3); hpdx[1][2]->Draw("colz");
  c1->cd(4); hrdx[1][2]->Draw("colz");
  c1->SaveAs("fcsTrkHcaldx.png");
  c1->Clear();

  c1->Divide(2,2);
  c1->cd(1); hxdy[1][2]->Draw("colz");
  c1->cd(2); hydy[1][2]->Draw("colz");
  c1->cd(3); hpdy[1][2]->Draw("colz");
  c1->cd(4); hrdy[1][2]->Draw("colz");
  c1->SaveAs("fcsTrkHcaldy.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdp[1][2]->Draw("colz");
  c1->cd(2); hydp[1][2]->Draw("colz");
  c1->cd(3); hpdp[1][2]->Draw("colz");
  c1->cd(4); hrdp[1][2]->Draw("colz");
  c1->SaveAs("fcsTrkHcaldp.png");

  c1->Clear();
  c1->Divide(2,2);
  c1->cd(1); hxdr[1][2]->Draw("colz");
  c1->cd(2); hydr[1][2]->Draw("colz");
  c1->cd(3); hpdr[1][2]->Draw("colz");
  c1->cd(4); hrdr[1][2]->Draw("colz");
  c1->SaveAs("fcsTrkHcaldr.png");

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

    // Legend + resolution summary
    c1->cd(5);
    TLatex *rLtx = new TLatex();
    rLtx->SetNDC(); rLtx->SetTextSize(0.055);
    rLtx->DrawLatex(0.05,0.92,"MC truth resolution");
    float rYpos=0.78;
    for(int tt=0; tt<mNType; tt++){
      if(!hTrkDqPt[tt]||hTrkDqPt[tt]->GetEntries()==0) continue;
      float sigRes=hTrkDqPt[tt]->GetRMS();
      float biasRes=hTrkDqPt[tt]->GetMean();
      rLtx->SetTextColor(TCOL[tt]);
      rLtx->DrawLatex(0.05, rYpos, Form("%s: #sigma(1/p_{T})=%.3f bias=%+.3f",
                                        TTYPE[tt], sigRes, biasRes));
      rYpos -= 0.13;
    }
    #undef DRAWTYPES_TRUTH

    c1->SaveAs("fcsTrkResolution.png");
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

  bool hasBLCVtx = hBLCVtxTrkDcaZ && hBLCVtxTrkDcaZ->GetEntries()>0;
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
}
