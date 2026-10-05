// fcsDyVsY.C -- dy as a function of fcsY, which separates a SCALE from an OFFSET.
//   flat        -> a pure translation between the FCS and track frames
//   sloped      -> a y scale; slope is the fractional mismatch
// Mixed event is normalised per slice over |d| 60-95 cm and subtracted.
// Lives in script/ rather than the scratchpad because the scratchpad does not
// survive a session restart.
// CINT: unique loop names, fixed arrays (see CLAUDE.md).

double fdSlope, fdInter;

void fcsDyVsY(const char* f, const char* type, const char* lab, const char* det = "Ecal") {
    fdSlope = 0; fdInter = 0;
    TFile* fp = TFile::Open(f);
    if (!fp || fp->IsZombie()) { printf("  %s: cannot open\n", lab); return; }
    TH2* hs = (TH2*)fp->Get(Form("ydy%sSame Event_%s", det, type));
    TH2* hm = (TH2*)fp->Get(Form("ydy%sMixed Event_%s", det, type));
    if (!hs || !hm) { printf("  %s %s %s: missing\n", lab, det, type); return; }
    const int kN = 4;
    double ylo[kN] = {-95, -40, 12, 40};
    double yhi[kN] = {-40, -12, 40, 95};
    double sx = 0, sy = 0, sxx = 0, sxy = 0; int n = 0;
    printf("\n  %s [%s] %s   dy vs fcsY\n", lab, type, det);
    for (int i = 0; i < kN; i++) {
        int b0 = hs->GetXaxis()->FindBin(ylo[i] + 1e-6);
        int b1 = hs->GetXaxis()->FindBin(yhi[i] - 1e-6);
        TH1D* ps = hs->ProjectionY(Form("fds_%s_%s_%d", det, type, i), b0, b1); ps->SetDirectory(0);
        TH1D* pm = hm->ProjectionY(Form("fdm_%s_%s_%d", det, type, i), b0, b1); pm->SetDirectory(0);
        double aS = 0, aM = 0;
        for (int jb = 1; jb <= ps->GetNbinsX(); jb++) {
            double x = fabs(ps->GetXaxis()->GetBinCenter(jb));
            if (x < 60 || x > 95) continue;
            aS += ps->GetBinContent(jb); aM += pm->GetBinContent(jb);
        }
        pm->Scale(aM > 0 ? aS / aM : 1.0);
        ps->Add(pm, -1.0);
        if (ps->Integral() < 400) { printf("   %5.0f..%-6.0f  %10s\n", ylo[i], yhi[i], "--"); continue; }
        int pb = ps->GetMaximumBin();
        double pk = ps->GetXaxis()->GetBinCenter(pb);
        TF1* g = new TF1(Form("fdg_%s_%s_%d", det, type, i), "gaus", pk - 8, pk + 8);
        ps->Fit(g, "QNR");
        double yc = 0.5 * (ylo[i] + yhi[i]), v = g->GetParameter(1);
        printf("   %5.0f..%-6.0f  %10.2f  %12.0f\n", ylo[i], yhi[i], v, ps->Integral());
        sx += yc; sy += v; sxx += yc * yc; sxy += yc * v; n++;
    }
    if (n >= 3) {
        fdSlope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
        fdInter = (sy - fdSlope * sx) / n;
        printf("   -> slope %+.5f per cm of y, intercept %+.2f cm\n", fdSlope, fdInter);
    } else {
        printf("   -> too few usable slices\n");
    }
    fp->Close();
}
