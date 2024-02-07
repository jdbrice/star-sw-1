
// C++ headers
#include <iostream>

// ROOT headers
#include <iostream>
#include <vector>
#include <stdio.h>


// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TChain.h"
#include <TList.h>
#include <TObject.h>

#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAttCanvas.h"
#include "TLegend.h"
#include "TVector2.h"

void SetLeg(TLegend* l,float txtsize=0.03,int color =1, int borderstyle =2){
    l->SetBorderSize(1);
    l->SetFillColor(0);
    l->SetLineColor(color);
    l->SetTextSize(txtsize);
    l->SetLineStyle(borderstyle);
}

void fillLegend(TLegend* l, TH1F* h1, TH1F* h2, const char* entry1, const char* entry2){
    

  float entries1 = 0; float mean1 = 0;
  float entries2 = 0; float mean2 = 0;

  entries1 = h1->GetEntries(); mean1 = h1->GetMean();
  entries2 = h2->GetEntries(); mean2 = h2->GetMean();

  h1->SetStats(0);
  h2->SetStats(0);

  SetLeg(l);
  l->SetTextSize(.035);
  l->SetFillColor(0);
  l->SetBorderSize(3);
  l->SetLineColor(1);
  l->AddEntry(h1,"Legend","");
  l->AddEntry(h1,Form("%s - ",entry1),"L");
  l->AddEntry(h1,Form("Entries = %.4e", entries1),"");
  l->AddEntry(h1,Form("Mean = %.2e", mean1),"");
  l->AddEntry(h2,Form("%s - ",entry2),"L");
  l->AddEntry(h2,Form("Entries = %.4e", entries2),"");
  l->AddEntry(h2,Form("Mean = %.2e", mean2),"");
  
}



void drawHist(TH1F* histMC, TH1F* hist, float minX = 0, float maxY = 0, bool legOption = true ){
    
    TLegend* leg= new TLegend(.6,.6,.9,.9);
    fillLegend(leg, histMC, hist, "MCTracks", "RecoTracks");
    
    if(minX && minX) histMC->GetXaxis()->SetRangeUser( minX, maxY);
    histMC->SetTitle("");
    histMC->SetStats(0);
    histMC->SetLineColor(6);
    histMC->Draw();
    hist->SetLineColor(4);
    hist->Draw("same");
    if(legOption) leg->Draw();
    
    
}


void analyzeStFwdAnalysisMaker( const char* fwdFile = "StFwdAnalysisMaker.root"){
    

    cout <<"Reading in "<< fwdFile << " for analysis."<< endl;
    
    TFile* file = TFile::Open(fwdFile);

    TH1F* eta = (TH1F*)file->Get("eta");
    TH1F* etaMC = (TH1F*)file->Get("etaMC");
    TH1F* phi = (TH1F*)file->Get("phi");
    TH1F* phiMC = (TH1F*)file->Get("phiMC");
    TH1F* px = (TH1F*)file->Get("px");
    TH1F* pxMC = (TH1F*)file->Get("pxMC");
    TH1F* py = (TH1F*)file->Get("py");
    TH1F* pyMC = (TH1F*)file->Get("pyMC");
    TH1F* pz = (TH1F*)file->Get("pz");
    TH1F* pzMC = (TH1F*)file->Get("pzMC");
    TH1F* pt = (TH1F*)file->Get("pt");
    TH1F* ptMC = (TH1F*)file->Get("ptMC");
    
    TH1F* charge = (TH1F*)file->Get("charge");
    TH1F* chargeMC = (TH1F*)file->Get("chargeMC");
    
    TH1F* ecalMatchPerTrack = (TH1F*)file->Get("ecalMatchPerTrack");
    TH1F* hcalMatchPerTrack = (TH1F*)file->Get("hcalMatchPerTrack");

    TH1F* fwdMultAll = (TH1F*)file->Get("fwdMultAll");
    TH1F* fwdMultGood = (TH1F*)file->Get("fwdMultGood");
    TH1F* fwdMultFST = (TH1F*)file->Get("fwdMultFST");
    TH1F* nHitsFit = (TH1F*)file->Get("nHitsFit");
    
    
    
    
    

// ================================================================= Histogram Plotting ===================================================
    
    TCanvas* c1 = new TCanvas("c1","Jets",800,600);
    c1->Divide(1,2);
    c1->cd(1);
    c1_1->SetLogy();
    drawHist(etaMC, eta);
    c1->cd(2);
    c1_2->SetLogy();
    drawHist(phiMC, phi);
    c1->Print("out_fwdanalysis.pdf(");
    
    TCanvas* c1 = new TCanvas("c1","Jets",800,600);
    c1->Divide(2,2);
    c1->cd(1);
    c1_1->SetLogy();
    drawHist(pxMC, px, -2, 5);
    c1->cd(2);
    c1_2->SetLogy();
    drawHist(pyMC, py, -2, 2);
    c1->cd(3);
    c1_3->SetLogy();
    drawHist(pzMC, pz);
    c1->cd(4);
    c1_4->SetLogy();
    drawHist(ptMC, pt, 0, 100);
    c1->Print("out_fwdanalysis.pdf");
    
    
    TCanvas* c1 = new TCanvas("c1","Jets",800,600);
    c1->Divide(1,2);
    c1->cd(1);
    ecalMatchPerTrack->SetLineColor(4);
    ecalMatchPerTrack->Draw();
    c1->cd(2);
    hcalMatchPerTrack->SetLineColor(4);
    hcalMatchPerTrack->Draw();
    c1->Print("out_fwdanalysis.pdf");
    
    TCanvas* c1 = new TCanvas("c1","Jets",800,600);
    c1->Divide(2,2);
    c1->cd(1);
    fwdMultAll->SetLineColor(4);
    fwdMultAll->Draw();
    c1->cd(2);
    fwdMultGood->SetLineColor(4);
    fwdMultGood->Draw();
    c1->cd(3);
    fwdMultFST->SetLineColor(4);
    fwdMultFST->Draw();
    c1->cd(4);
    nHitsFit->SetLineColor(4);
    nHitsFit->Draw();
    c1->Print("out_fwdanalysis.pdf");
    
    TCanvas* c1 = new TCanvas("c1","Jets",800,600);
    c1->Divide(1,1);
    c1->cd(1);
    drawHist(chargeMC, charge, 0,0, false);
    c1->Print("out_fwdanalysis.pdf)");
    

    
}
