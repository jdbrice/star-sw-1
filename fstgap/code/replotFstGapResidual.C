// replotFstGapResidual.C
//
// Re-renders the residual figure from the cached histograms written by
// plotFstGapResidual.C -- no second pass over ~58M rows per side.
//
// Adds an ON/OFF ratio panel under each distribution. Without it the change is
// invisible: a 4.9% narrowing of the core is real and significant but the two
// curves sit on top of each other on a log plot, so the eye reads "no
// difference" where the numbers say otherwise.
//
// Usage: root4star -b -q 'script/replotFstGapResidual.C'
void replotFstGapResidual(const char* dir = "fstgap"){
    gStyle->SetOptStat(0);
    TFile* f = TFile::Open(Form("%s/gapResidualHistos.root", dir));
    if (!f){ printf("no cache\n"); return; }
    // hres_[pass*2+region], region 0=inner 1=outer
    TH1F* h[4];
    for (int i = 0; i < 4; i++){
        h[i] = (TH1F*) f->Get(Form("hres_%d", i));
        if (!h[i]){ printf("missing hres_%d\n", i); return; }
        h[i]->SetDirectory(0);
        if (h[i]->Integral() > 0) h[i]->Scale(1.0/h[i]->Integral());
    }
    TCanvas* c = new TCanvas("c_res2", "", 1000, 560);
    c->Divide(2,1);
    for (int ip = 0; ip < 2; ip++){
        int reg = (ip == 0) ? 1 : 0;              // outer first
        TH1F* a = h[0*2 + reg];                    // OFF
        TH1F* b = h[1*2 + reg];                    // ON
        c->cd(ip+1);
        TPad* pT = new TPad(Form("pT%d",ip), "", 0, 0.32, 1, 1);
        TPad* pB = new TPad(Form("pB%d",ip), "", 0, 0.0,  1, 0.34);
        pT->SetBottomMargin(0.02); pT->SetLeftMargin(0.14); pT->SetLogy(); pT->Draw();
        pB->SetTopMargin(0.02); pB->SetBottomMargin(0.30); pB->SetLeftMargin(0.14);
        pB->SetGridy(); pB->Draw();

        pT->cd();
        a->SetLineColor(kRed+1);  a->SetLineWidth(2); a->SetFillStyle(0);
        b->SetLineColor(kBlue+1); b->SetLineWidth(2); b->SetFillStyle(0);
        a->SetTitle(reg ? "OUTER sensors (r>16.5 cm)" : "INNER sensor (r<16.5 cm) -- control");
        a->GetYaxis()->SetTitle("fraction of rows / bin");
        a->GetYaxis()->SetTitleSize(0.05); a->GetYaxis()->SetLabelSize(0.045);
        a->GetXaxis()->SetLabelSize(0);
        a->Draw("hist"); b->Draw("hist same");
        TLegend* lg = new TLegend(0.17, 0.72, 0.50, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0);
        lg->AddEntry(a, "gap fix OFF", "l");
        lg->AddEntry(b, "gap fix ON",  "l");
        lg->Draw();

        pB->cd();
        TH1F* r = (TH1F*) b->Clone(Form("ratio%d", ip));
        r->Divide(a);
        r->SetLineColor(kBlack); r->SetLineWidth(2); r->SetFillStyle(0);
        r->SetTitle("");
        r->GetYaxis()->SetTitle("ON / OFF");
        r->GetYaxis()->SetRangeUser(0.80, 1.20);
        r->GetYaxis()->SetNdivisions(505);
        r->GetYaxis()->SetTitleSize(0.11); r->GetYaxis()->SetLabelSize(0.10);
        r->GetYaxis()->SetTitleOffset(0.55);
        r->GetXaxis()->SetTitle("r#Delta#phi [cm]");
        r->GetXaxis()->SetTitleSize(0.12); r->GetXaxis()->SetLabelSize(0.10);
        r->Draw("hist");
        TLine* one = new TLine(-1, 1, 1, 1); one->SetLineStyle(2); one->Draw();
    }
    c->SaveAs(Form("%s/fstGapResidual.png", dir));
    printf("wrote %s/fstGapResidual.png\n", dir);
}
