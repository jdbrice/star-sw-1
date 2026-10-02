// plotFstResidCmp.C -- FST r*dphi residual, Run 22 vs Run 24, one panel per FST disk.
//
// Both years are taken from the PICKUP arm (fttNoAdd=0) so the comparison is like for
// like: adding FTT hits to the fit widens the FST residual, so mixing arms would flatter
// whichever year came from the alignment arm. FST hits seed the track, so this is a fit
// residual, not an independent one -- its width reflects how well the whole forward fit
// is determined, which is exactly the point being made.
//
// Areas are normalised because the two samples differ in size (4.2M vs 1.5M tracks).
//
// CINT: unique loop variable names, fixed arrays (see CLAUDE.md).

void plotFstResidCmp(const char* f22 = "fttadd_Primary.root",
                     const char* f24 = "run24iter1_Primary.root",
                     const char* out = "FstFttFlipTest/plots/fst_resid_run22_run24.png") {
    gStyle->SetOptStat(0);
    TFile* fa = TFile::Open(f22);
    TFile* fb = TFile::Open(f24);
    if (!fa || fa->IsZombie() || !fb || fb->IsZombie()) { printf("cannot open inputs\n"); return; }

    TCanvas* cfr = new TCanvas("cfr", "cfr", 1200, 420);
    cfr->Divide(3, 1);
    for (int id = 0; id < 3; id++) {
        cfr->cd(id + 1);
        gPad->SetGridx(); gPad->SetLogy();
        TH1F* ha = (TH1F*)fa->Get(Form("FST/disk%d/h_fst_d%d_rdphi", id, id));
        TH1F* hb = (TH1F*)fb->Get(Form("FST/disk%d/h_fst_d%d_rdphi", id, id));
        if (!ha || !hb) { printf("  missing disk %d\n", id); continue; }
        TH1F* a = (TH1F*)ha->Clone(Form("a%d", id)); a->SetDirectory(0);
        TH1F* b = (TH1F*)hb->Clone(Form("b%d", id)); b->SetDirectory(0);
        if (a->Integral() > 0) a->Scale(1.0 / a->Integral());
        if (b->Integral() > 0) b->Scale(1.0 / b->Integral());
        a->SetLineColor(kGray + 2);  a->SetLineWidth(2);
        b->SetLineColor(kAzure + 2); b->SetLineWidth(2);
        a->SetTitle(Form("FST disk %d:  r#upoint#Delta#phi;r#upoint#Delta#phi [cm];fraction / bin", id));
        a->GetXaxis()->SetRangeUser(-1.5, 1.5);
        double ymax = TMath::Max(a->GetMaximum(), b->GetMaximum());
        a->GetYaxis()->SetRangeUser(ymax * 1e-4, ymax * 3);
        a->Draw("hist");
        b->Draw("hist same");
        TLegend* lg = new TLegend(0.13, 0.74, 0.62, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038);
        lg->AddEntry(a, Form("Run 22  rms %.3f cm", ha->GetRMS()), "l");
        lg->AddEntry(b, Form("Run 24  rms %.3f cm", hb->GetRMS()), "l");
        lg->Draw();
        TLatex t; t.SetNDC(); t.SetTextSize(0.036);
        t.DrawLatex(0.13, 0.69, Form("narrower by #times%.1f", (hb->GetRMS() > 0) ? ha->GetRMS() / hb->GetRMS() : 0));
    }
    cfr->SaveAs(out);
    printf("wrote %s\n", out);
}
