// plotBlindCrosscond.C
//
// Draws the cross-conditioned FTT-blind Dx/Dy/XY diagnostic histograms from
// fwd_blind_diag.root (see StRoot/StFwdTrackMaker/include/Tracker/FwdTracker.h,
// findFttStripsNearProjectedState()) as RAW COUNTS (no area normalization),
// with the off-axis 1D range widened to +/-35cm (2026-07-20).
//
// Usage: root4star -l -b -q 'plotBlindCrosscond.C("fwd_blind_diag.root","outdir")'

void plotBlindCrosscond(const char* infile = "fwd_blind_diag.root", const char* outdir = "20260720_blinddiag_crosscond") {
    gSystem->mkdir(outdir, true);
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", infile); return; }

    const int nDisk = 4;

    // ---- 1D Dx/Dy, both orientations, all 4 disks -----------------------
    const char* pngBase[4] = {"Dx_V", "Dy_V", "Dx_H", "Dy_H"};
    for (int d = 0; d < nDisk; d++) {
        TString allNames[4];
        TString matchedNames[4];
        allNames[0] = Form("hBlindDxAll_V_disk%d", d); matchedNames[0] = Form("hBlindDxMatched_V_disk%d", d);
        allNames[1] = Form("hBlindDyAll_V_disk%d", d); matchedNames[1] = Form("hBlindDyMatched_V_disk%d", d);
        allNames[2] = Form("hBlindDxAll_H_disk%d", d); matchedNames[2] = Form("hBlindDxMatched_H_disk%d", d);
        allNames[3] = Form("hBlindDyAll_H_disk%d", d); matchedNames[3] = Form("hBlindDyMatched_H_disk%d", d);

        for (int ip = 0; ip < 4; ip++) {
            TH1F *hAll = (TH1F*)f->Get(allNames[ip]);
            TH1F *hMatched = (TH1F*)f->Get(matchedNames[ip]);
            if (!hAll || !hMatched) { printf("Missing %s or %s\n", allNames[ip].Data(), matchedNames[ip].Data()); continue; }

            TCanvas *c = new TCanvas(Form("c_%d_%d", d, ip), "c", 700, 550);
            c->SetLogy(1);

            // hAll is drawn first and its frame sets the visible y-range --
            // ROOT does NOT auto-rescale that frame for a later "same" draw,
            // so if hMatched's peak bin is taller than hAll's, it silently
            // gets clipped by the frame's top edge (looked like "red poking
            // above the gray" -- it wasn't above, it was overflowing the
            // frame). Fix: set the frame's max from whichever histogram is
            // actually taller, with headroom. Fix y-min at 1e1 (requested) so
            // low-count "Matched" bins near the axis floor stay visible
            // instead of being squashed flat by an auto-min set too high.
            double ymax = TMath::Max(hAll->GetMaximum(), hMatched->GetMaximum());
            hAll->SetMinimum(10);
            hAll->SetMaximum(ymax * 2.0);

            hAll->SetFillColor(kGray);
            hAll->SetLineColor(kGray+2);
            hAll->Draw("hist");

            hMatched->SetLineColor(kRed);
            hMatched->SetLineWidth(2);
            hMatched->Draw("hist same");

            TLegend *leg = new TLegend(0.62, 0.75, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(hAll, "All candidates", "f");
            leg->AddEntry(hMatched, "Matched", "l");
            leg->Draw();

            c->SaveAs(Form("%s/blind_%s_disk%d.png", outdir, pngBase[ip], d));
            delete c;
        }
    }

    // ---- 2D XY, both orientations, all 4 disks (all vs matched side by side) ----
    for (int d = 0; d < nDisk; d++) {
        const char* oris[2] = {"V", "H"};
        for (int io = 0; io < 2; io++) {
            const char* ori = oris[io];
            TH2F *hAll2D = (TH2F*)f->Get(Form("hBlindXYAll_%s_disk%d", ori, d));
            TH2F *hMatched2D = (TH2F*)f->Get(Form("hBlindXYMatched_%s_disk%d", ori, d));
            if (!hAll2D || !hMatched2D) { printf("Missing XY %s disk%d\n", ori, d); continue; }

            TCanvas *c = new TCanvas(Form("c2_%d_%d", d, io), "c2", 1100, 500);
            c->Divide(2, 1);
            c->cd(1);
            gPad->SetLogz(1);
            hAll2D->Draw("colz");
            c->cd(2);
            gPad->SetLogz(1);
            hMatched2D->Draw("colz");

            c->SaveAs(Form("%s/blind_XY_%s_disk%d.png", outdir, ori, d));
            delete c;
        }
    }

    printf("Done -- wrote plots to %s/\n", outdir);
}
