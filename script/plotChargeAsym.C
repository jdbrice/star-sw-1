// plotChargeAsym.C -- plots for the per-quadrant charge asymmetry.
//
//   root -l -b -q 'script/plotChargeAsym.C("/gpfs01/star/pwg_tasks/FwdCalib/akio/qc_20261007","FstFttFlipTest/chargeasym")'
//
// Reads merged_<type>.root (hadd of the afterburner's FwdDetResidual_<type>.root)
// and writes two plots:
//
//   chargeasym_vs_type.png    A per quadrant as a function of track type, plus
//                             the north-south swing. The headline: the asymmetry
//                             exists for Global (no vertex) and the vertex
//                             constraint multiplies the swing by about four.
//   chargeasym_vs_plane.png   A per detector plane, one pad per quadrant, one
//                             curve per track type. Flat across planes means the
//                             asymmetry is a property of the TRACK, not of a
//                             particular detector.
//
// A = (N+ - N-)/(N+ + N-) from PlaneUsageQuadCharge/hPlaneUsage_<quad>_<Pos|Neg>.
// Missing types are skipped, so this can be run before every merge has finished.
//
// CINT: no local std::vector/map, unique loop variable names, 0 not nullptr.

void plotChargeAsym( TString dir    = "/gpfs01/star/pwg_tasks/FwdCalib/akio/qc_20261007",
                     TString outdir = "FstFttFlipTest/chargeasym" )
{
    gStyle->SetOptStat(0);
    gStyle->SetPadGridY(1);

    const int NT = 5, NQ = 4;
    const char* types[NT] = {"Global","BLC","BLCVtx","Primary","FCSConstrained"};
    const char* quads[NQ] = {"Ntop","Nbot","Stop","Sbot"};
    int qcol[NQ]  = {kRed+1, kBlue+1, kGreen+2, kMagenta+1};
    int qmark[NQ] = {20, 21, 22, 23};
    int tcol[NT]  = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1};
    int tmark[NT] = {20, 21, 22, 23, 29};

    // ---- gather -------------------------------------------------------------
    // asym[type][quad] at FST1, and asymP[type][quad][plane] across planes
    double asym[NT][NQ];
    double swing[NT];
    int    have[NT];
    const int MAXP = 24;
    double asymP[NT][NQ][MAXP];
    int    nplane = 0;
    TString planeName[MAXP];

    for (int it = 0; it < NT; it++) {
        have[it] = 0; swing[it] = 0;
        for (int iq = 0; iq < NQ; iq++) {
            asym[it][iq] = 0;
            for (int ip = 0; ip < MAXP; ip++) asymP[it][iq][ip] = -999;
        }
        TFile* f = TFile::Open(Form("%s/merged_%s.root", dir.Data(), types[it]));
        if (!f || f->IsZombie()) { printf("  skip %s (no merged file)\n", types[it]); continue; }
        int okq = 0;
        for (int jq = 0; jq < NQ; jq++) {
            TH1F* hp = (TH1F*)f->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Pos",quads[jq]));
            TH1F* hn = (TH1F*)f->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Neg",quads[jq]));
            if (!hp || !hn) continue;
            okq++;
            int nb = hp->GetNbinsX(); if (nb > MAXP) nb = MAXP;
            if (nplane == 0) {
                nplane = nb;
                for (int kb = 1; kb <= nb; kb++) planeName[kb-1] = hp->GetXaxis()->GetBinLabel(kb);
            }
            int b1 = hp->GetXaxis()->FindBin("FST1");
            double np1 = hp->GetBinContent(b1), nn1 = hn->GetBinContent(b1);
            asym[it][jq] = (np1+nn1 > 0) ? (np1-nn1)/(np1+nn1) : 0;
            for (int kp = 1; kp <= nb; kp++) {
                double np = hp->GetBinContent(kp), nn = hn->GetBinContent(kp);
                // only fill where there is real occupancy; 0/0 would read as A=0
                if (np + nn > 1000) asymP[it][jq][kp-1] = (np-nn)/(np+nn);
            }
        }
        if (okq == NQ) {
            have[it] = 1;
            swing[it] = 0.5*(asym[it][2]+asym[it][3]) - 0.5*(asym[it][0]+asym[it][1]);
        }
        f->Close();
    }

    // Only join the points when every type is present; otherwise a line would
    // run straight from Global to FCSConstrained across the missing ones and
    // read as a trend through them.
    int nhave = 0;
    for (int ih = 0; ih < NT; ih++) nhave += have[ih];
    TString dopt = (nhave == NT) ? "LP same" : "P same";

    // ---- plot 1: asymmetry vs track type ------------------------------------
    TCanvas* c1 = new TCanvas("c1","charge asymmetry vs track type",1100,520);
    c1->Divide(2,1);

    c1->cd(1);
    gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14);
    TH2F* fr1 = new TH2F("fr1",";;A = (N^{+} - N^{-}) / (N^{+} + N^{-})", NT,0,NT, 100,-0.45,0.35);
    for (int ib = 0; ib < NT; ib++) fr1->GetXaxis()->SetBinLabel(ib+1, types[ib]);
    fr1->GetXaxis()->SetLabelSize(0.045);
    fr1->GetYaxis()->SetTitleOffset(1.5);
    fr1->Draw();
    TLine* z1 = new TLine(0,0,NT,0); z1->SetLineStyle(2); z1->Draw();
    TLegend* lg1 = new TLegend(0.60,0.70,0.89,0.89); lg1->SetBorderSize(0); lg1->SetFillStyle(0);
    for (int iq2 = 0; iq2 < NQ; iq2++) {
        TGraph* g = new TGraph();
        int n = 0;
        for (int it2 = 0; it2 < NT; it2++)
            if (have[it2]) { g->SetPoint(n, it2+0.5, asym[it2][iq2]); n++; }
        if (n == 0) continue;
        g->SetLineColor(qcol[iq2]); g->SetMarkerColor(qcol[iq2]);
        g->SetMarkerStyle(qmark[iq2]); g->SetMarkerSize(1.3); g->SetLineWidth(2);
        g->Draw(dopt);
        lg1->AddEntry(g, quads[iq2], (nhave == NT) ? "lp" : "p");
    }
    lg1->Draw();

    c1->cd(2);
    gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14);
    TH2F* fr2 = new TH2F("fr2",";;north-south swing  #LTA#GT_{S} - #LTA#GT_{N}", NT,0,NT, 100,0,0.5);
    for (int ic = 0; ic < NT; ic++) fr2->GetXaxis()->SetBinLabel(ic+1, types[ic]);
    fr2->GetXaxis()->SetLabelSize(0.045);
    fr2->GetYaxis()->SetTitleOffset(1.7);
    fr2->Draw();
    TGraph* gs = new TGraph();
    int ns = 0;
    for (int it3 = 0; it3 < NT; it3++)
        if (have[it3]) { gs->SetPoint(ns, it3+0.5, swing[it3]); ns++; }
    gs->SetLineColor(kBlack); gs->SetMarkerColor(kBlack);
    gs->SetMarkerStyle(20); gs->SetMarkerSize(1.4); gs->SetLineWidth(2);
    gs->Draw(dopt);
    TLatex* tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.035);
    if (nhave == NT) tx->DrawLatex(0.20,0.80,"vertex constraint multiplies the swing by ~4");
    else             tx->DrawLatex(0.20,0.80,Form("PARTIAL: only %d of %d types merged", nhave, NT));

    c1->SaveAs(Form("%s/chargeasym_vs_type.png", outdir.Data()));

    // ---- plot 2: asymmetry vs plane, one pad per quadrant --------------------
    TCanvas* c2 = new TCanvas("c2","charge asymmetry vs plane",1200,900);
    c2->Divide(2,2);
    for (int iq3 = 0; iq3 < NQ; iq3++) {
        c2->cd(iq3+1);
        gPad->SetBottomMargin(0.20); gPad->SetLeftMargin(0.12); gPad->SetTopMargin(0.08);
        TH2F* fr = new TH2F(Form("fr_%d",iq3), Form("%s;;A",quads[iq3]),
                            nplane,0,nplane, 100,-0.55,0.72);
        for (int ip2 = 0; ip2 < nplane; ip2++) fr->GetXaxis()->SetBinLabel(ip2+1, planeName[ip2]);
        fr->GetXaxis()->LabelsOption("v");
        fr->GetXaxis()->SetLabelSize(0.05);
        fr->Draw();
        TLine* zz = new TLine(0,0,nplane,0); zz->SetLineStyle(2); zz->Draw();
        TLegend* lg = new TLegend(0.12,0.73,0.90,0.89); lg->SetBorderSize(0); lg->SetFillStyle(0);
        lg->SetTextSize(0.034); lg->SetNColumns(3);
        for (int it4 = 0; it4 < NT; it4++) {
            if (!have[it4]) continue;
            TGraph* gp = new TGraph();
            int np2 = 0;
            for (int ip3 = 0; ip3 < nplane; ip3++)
                if (asymP[it4][iq3][ip3] > -900) { gp->SetPoint(np2, ip3+0.5, asymP[it4][iq3][ip3]); np2++; }
            if (np2 == 0) continue;
            gp->SetLineColor(tcol[it4]); gp->SetMarkerColor(tcol[it4]);
            gp->SetMarkerStyle(tmark[it4]); gp->SetMarkerSize(1.1); gp->SetLineWidth(2);
            gp->Draw("LP same");  // planes are always complete within a type
            lg->AddEntry(gp, types[it4], "lp");
        }
        lg->Draw();
    }
    c2->SaveAs(Form("%s/chargeasym_vs_plane.png", outdir.Data()));

    printf("\n  wrote %s/chargeasym_vs_type.png and chargeasym_vs_plane.png\n", outdir.Data());
    printf("  types included:");
    for (int it5 = 0; it5 < NT; it5++) if (have[it5]) printf(" %s", types[it5]);
    printf("\n\n");
}
