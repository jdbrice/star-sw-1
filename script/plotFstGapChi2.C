// plotFstGapChi2.C
//
// The headline figure for fstgap/index.html: chi2/ndf of the leave-one-out
// refit, off vs on, split by whether the REMOVED hit was on an outer or an
// inner FST sensor, versus nPointsUsed.
//
// Only outer hits move when the gap fix is applied, so a real improvement has
// to concentrate there; the inner rows are the control. Slicing on
// nPointsUsed fixes the track composition, so the comparison is not a
// selection effect.
//
// Numbers hard-coded from the final campaign (152 matched files per side,
// ~58M rows each) rather than re-read, so the figure cannot silently drift
// from the numbers quoted on the page. Source: gapfix_test/CONCLUSIONS.txt
// section 6, produced by script/compareGapfixCampaign.C.
//
// Usage: root4star -b -q 'script/plotFstGapChi2.C'

void plotFstGapChi2(const char* outdir = "fstgap"){

    gSystem->mkdir(outdir, kTRUE);
    gStyle->SetOptStat(0);

    const int N = 9;
    double np[9]   = {6,7,8,9,10,11,12,13,14};
    // outer sensors: the hits the fix moves
    double oOff[9] = {4.8910,5.3672,6.1670,6.7764,7.3607,7.6763,8.3562,9.2429,10.5867};
    double oOne[9] = {4.4773,4.9747,5.8695,6.4931,7.0600,7.6016,8.2601,9.2858,10.4755};
    double oEo[9]  = {0.0209,0.0173,0.0181,0.0191,0.0212,0.0231,0.0277,0.0346,0.0466};
    double oEn[9]  = {0.0190,0.0157,0.0166,0.0176,0.0192,0.0218,0.0259,0.0333,0.0442};
    // inner sensor: control, untouched by the fix
    double iOff[9] = {4.7310,5.8801,6.7925,7.0331,7.3043,7.0605,6.9577,7.5292,8.8704};
    double iOne[9] = {4.6716,5.8282,6.7428,6.9890,7.2673,7.0344,6.9527,7.5422,8.8994};
    double iEo[9]  = {0.0123,0.0099,0.0085,0.0082,0.0087,0.0089,0.0091,0.0103,0.0133};
    double iEn[9]  = {0.0120,0.0097,0.0084,0.0082,0.0086,0.0089,0.0090,0.0103,0.0134};

    // ---- left panel: chi2/ndf itself, off vs on ----
    TGraphErrors* gOo = new TGraphErrors(N, np, oOff, 0, oEo);
    TGraphErrors* gOn = new TGraphErrors(N, np, oOne, 0, oEn);
    TGraphErrors* gIo = new TGraphErrors(N, np, iOff, 0, iEo);
    TGraphErrors* gIn = new TGraphErrors(N, np, iOne, 0, iEn);

    // ---- right panel: the difference, which is the actual result ----
    double dO[9], dOe[9], dI[9], dIe[9];
    for (int i = 0; i < N; i++){
        dO[i]  = oOne[i] - oOff[i];
        dOe[i] = sqrt(oEo[i]*oEo[i] + oEn[i]*oEn[i]);
        dI[i]  = iOne[i] - iOff[i];
        dIe[i] = sqrt(iEo[i]*iEo[i] + iEn[i]*iEn[i]);
    }
    TGraphErrors* gDO = new TGraphErrors(N, np, dO, 0, dOe);
    TGraphErrors* gDI = new TGraphErrors(N, np, dI, 0, dIe);

    TCanvas* c = new TCanvas("c_gapchi2", "", 1100, 450);
    c->Divide(2,1);

    c->cd(1);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetGridy();
    TH1F* fr1 = gPad->DrawFrame(5.4, 4.0, 14.6, 11.5);
    fr1->GetXaxis()->SetTitle("nPointsUsed in the leave-one-out refit");
    fr1->GetYaxis()->SetTitle("#chi^{2}/ndf");
    fr1->GetYaxis()->SetTitleOffset(1.3);
    gOo->SetMarkerStyle(24); gOo->SetMarkerColor(kRed+1);  gOo->SetLineColor(kRed+1);  gOo->SetMarkerSize(1.2);
    gOn->SetMarkerStyle(20); gOn->SetMarkerColor(kRed+1);  gOn->SetLineColor(kRed+1);  gOn->SetMarkerSize(1.2);
    gIo->SetMarkerStyle(24); gIo->SetMarkerColor(kBlue+1); gIo->SetLineColor(kBlue+1); gIo->SetMarkerSize(1.2);
    gIn->SetMarkerStyle(20); gIn->SetMarkerColor(kBlue+1); gIn->SetLineColor(kBlue+1); gIn->SetMarkerSize(1.2);
    gOo->Draw("PL same"); gOn->Draw("PL same");
    gIo->Draw("PL same"); gIn->Draw("PL same");
    TLegend* l1 = new TLegend(0.16,0.66,0.55,0.88);
    l1->SetBorderSize(0); l1->SetFillStyle(0);
    l1->AddEntry(gOo, "outer sensors, fix OFF", "pl");
    l1->AddEntry(gOn, "outer sensors, fix ON",  "pl");
    l1->AddEntry(gIo, "inner sensor, OFF (control)", "pl");
    l1->AddEntry(gIn, "inner sensor, ON  (control)", "pl");
    l1->Draw();

    c->cd(2);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetGridy();
    TH1F* fr2 = gPad->DrawFrame(5.4, -0.50, 14.6, 0.12);
    fr2->GetXaxis()->SetTitle("nPointsUsed in the leave-one-out refit");
    fr2->GetYaxis()->SetTitle("#Delta#chi^{2}/ndf  (ON - OFF)");
    fr2->GetYaxis()->SetTitleOffset(1.3);
    gDO->SetMarkerStyle(20); gDO->SetMarkerColor(kRed+1);  gDO->SetLineColor(kRed+1);  gDO->SetMarkerSize(1.3); gDO->SetLineWidth(2);
    gDI->SetMarkerStyle(20); gDI->SetMarkerColor(kBlue+1); gDI->SetLineColor(kBlue+1); gDI->SetMarkerSize(1.3); gDI->SetLineWidth(2);
    gDO->Draw("PL same"); gDI->Draw("PL same");
    TLine* z = new TLine(5.4, 0, 14.6, 0); z->SetLineStyle(2); z->Draw();
    TLegend* l2 = new TLegend(0.42,0.20,0.88,0.36);
    l2->SetBorderSize(0); l2->SetFillStyle(0);
    l2->AddEntry(gDO, "outer sensors (hits that move)", "pl");
    l2->AddEntry(gDI, "inner sensor (control)", "pl");
    l2->Draw();

    c->SaveAs(Form("%s/fstGapChi2.png", outdir));
    printf("wrote %s/fstGapChi2.png\n", outdir);
}
