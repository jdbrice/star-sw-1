// plotTrackMatch.C -- reproduce the 2026-07-06 trackmatch plots, same panels, same
// ranges, same 4 cm binning, for a chosen dataset and track type.
//
// Two products per detector, named as on the old page so they overlay one for one:
//   fcsTrkMatch<Det>.png   dX, dY, R*dPhi, dR, each with Same / Mixed / Same-Mixed
//   fcsTrk<Det>d{x,y,p,r}.png   the 2D dependence on fcsX, fcsY, fcsR*fcsPhi, fcsR
//
// picoMatch books the Same-Mixed slots but never fills them, so they are built here:
// mixed is normalised to same over |d| 80-160 cm, where same/mixed is flat, and
// subtracted. The macro also prints the match count and the signal-to-background, which
// a width alone does not show.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

double tmSig, tmBg, tmSB, tmRms;

TH1D* tmSub(TH1* hs, TH1* hm, int rebin, const char* tag) {
    if (!hs || !hm) return 0;
    TH1D* d = (TH1D*)hs->Clone(Form("tms_%s", tag)); d->SetDirectory(0);
    TH1D* m = (TH1D*)hm->Clone(Form("tmm_%s", tag)); m->SetDirectory(0);
    if (rebin > 1) { d->Rebin(rebin); m->Rebin(rebin); }
    double aS = 0, aM = 0;
    for (int ib = 1; ib <= d->GetNbinsX(); ib++) {
        double x = fabs(d->GetXaxis()->GetBinCenter(ib));
        if (x < 80 || x > 160) continue;
        aS += d->GetBinContent(ib); aM += m->GetBinContent(ib);
    }
    m->Scale(aM > 0 ? aS / aM : 1.0);
    d->Add(m, -1.0);
    // signal and background inside |d|<15 cm, where the match lives
    double sig = 0, bg = 0;
    for (int jb = 1; jb <= d->GetNbinsX(); jb++) {
        double x = fabs(d->GetXaxis()->GetBinCenter(jb));
        if (x > 15) continue;
        sig += d->GetBinContent(jb);
        bg  += m->GetBinContent(jb);
    }
    tmSig = sig; tmBg = bg; tmSB = (bg > 0) ? sig / bg : 0;
    TH1D* q = (TH1D*)d->Clone(Form("tmq_%s", tag)); q->SetDirectory(0);
    q->GetXaxis()->SetRangeUser(-30, 30);
    tmRms = q->GetRMS();
    return d;
}

void plotTrackMatch(const char* file, const char* type, const char* det,
                    const char* outdir, const char* label, int rebin = 8) {
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", file); return; }
    gStyle->SetOptStat(0);
    const char* vv[4] = {"dx", "dy", "dp", "dr"};
    const char* vt[4] = {"dX", "dY", "fcsR*dPhi", "dR"};
    const char* n1[4] = {"North", "Top", "North", "R<70"};
    const char* n2[4] = {"South", "Bottom", "South", "R>70"};

    printf("\n=== %s  [%s, %s] ===\n", label, type, det);
    printf("  %-16s %11s %11s %11s %8s %8s\n", "panel", "same", "mixed", "signal", "S/B", "rms");

    TCanvas* c = new TCanvas(Form("tmm_%s_%s", det, type), "match", 1000, 760);
    c->Divide(2, 2);
    for (int iv = 0; iv < 4; iv++) {
        c->cd(iv + 1); gPad->SetLeftMargin(0.13);
        TH1* hs1 = (TH1*)f->Get(Form("%s%s%sSame Event_%s", det, n1[iv], vv[iv], type));
        TH1* hm1 = (TH1*)f->Get(Form("%s%s%sMixed Event_%s", det, n1[iv], vv[iv], type));
        TH1* hs2 = (TH1*)f->Get(Form("%s%s%sSame Event_%s", det, n2[iv], vv[iv], type));
        TH1* hm2 = (TH1*)f->Get(Form("%s%s%sMixed Event_%s", det, n2[iv], vv[iv], type));
        if (!hs1 || !hm1) continue;
        TH1D* same = (TH1D*)hs1->Clone(Form("tmsa_%s%d", det, iv)); same->SetDirectory(0);
        TH1D* mix  = (TH1D*)hm1->Clone(Form("tmmi_%s%d", det, iv)); mix->SetDirectory(0);
        if (rebin > 1) { same->Rebin(rebin); mix->Rebin(rebin); }
        double aS = 0, aM = 0;
        for (int ib = 1; ib <= same->GetNbinsX(); ib++) {
            double x = fabs(same->GetXaxis()->GetBinCenter(ib));
            if (x < 80 || x > 160) continue;
            aS += same->GetBinContent(ib); aM += mix->GetBinContent(ib);
        }
        mix->Scale(aM > 0 ? aS / aM : 1.0);
        TH1D* d1 = tmSub(hs1, hm1, rebin, Form("%s%d_1", det, iv));
        double s1 = tmSig, b1 = tmBg, sb1 = tmSB, r1 = tmRms;
        TH1D* d2 = tmSub(hs2, hm2, rebin, Form("%s%d_2", det, iv));
        double s2 = tmSig, b2 = tmBg, sb2 = tmSB;
        printf("  %-16s %11.4g %11.4g %11.4g %8.3f %8.2f\n",
               Form("%s %s", n1[iv], vt[iv]), same->Integral(), mix->Integral(), s1, sb1, r1);
        printf("  %-16s %11s %11s %11.4g %8.3f %8.2f\n",
               Form("%s %s", n2[iv], vt[iv]), "", "", s2, sb2, tmRms);

        same->SetTitle(Form("%s%s-Trk %s;%s[cm];", det, n1[iv], vt[iv], vt[iv]));
        same->SetLineColor(kBlack); same->SetLineWidth(2);
        mix->SetLineColor(kBlue);   mix->SetLineWidth(2);
        same->GetXaxis()->SetRangeUser(-200, 200);
        same->Draw("hist"); mix->Draw("hist same");
        if (d1) { d1->SetLineColor(kRed);   d1->SetLineWidth(2); d1->Draw("hist same"); }
        if (d2) { d2->SetLineColor(kGreen+2); d2->SetLineWidth(2); d2->Draw("hist same"); }
        TLegend* lg = new TLegend(0.14, 0.66, 0.52, 0.88);
        lg->SetFillColor(0); lg->SetBorderSize(0); lg->SetTextSize(0.035);
        lg->AddEntry(same, "Same Event", "l");
        lg->AddEntry(mix,  "Mixed Event", "l");
        if (d1) lg->AddEntry(d1, Form("Same-Mixed %s", n1[iv]), "l");
        if (d2) lg->AddEntry(d2, Form("Same-Mixed %s", n2[iv]), "l");
        lg->Draw();
    }
    c->SaveAs(Form("%s/fcsTrkMatch%s.png", outdir, det));

    // the 2D dependence, one canvas per residual, four variables -- the old fcsTrk<Det>d*.png
    gStyle->SetOptStat(1111);
    const char* xv[4] = {"x", "y", "p", "r"};
    const char* xt[4] = {"fcsX", "fcsY", "fcsR*fcsPhi", "fcsR"};
    for (int jv = 0; jv < 4; jv++) {
        TCanvas* c2 = new TCanvas(Form("tm2_%s_%s_%d", det, type, jv), "dep", 1000, 760);
        c2->Divide(2, 2);
        for (int kx = 0; kx < 4; kx++) {
            c2->cd(kx + 1);
            TH2* h2 = (TH2*)f->Get(Form("%s%s%sSame Event_%s", xv[kx], vv[jv], det, type));
            if (!h2) continue;
            h2->SetTitle(Form("%s(%s) %s;%s[cm];%s[cm]", vt[jv], xt[kx], det, xt[kx], vt[jv]));
            h2->Draw("colz");
        }
        c2->SaveAs(Form("%s/fcsTrk%s%s.png", outdir, det, vv[jv]));
    }
    gStyle->SetOptStat(0);
}
