// fwoPairPlot.C -- the outer overlap width from the pair centroid.
//
//   root -l -b -q 'script/fwoPairPlot.C("fwo_pair.root")'
//
// Two panels. LEFT: the pair map itself, my distance from the wedge edge
// against the partner's. Pairs of strips at the same physical angle lie on
// d_me + d_nb = const, because the two grids run opposite ways through the
// overlap band, so the ridge position IS the width: W = <s> + 1.
// RIGHT: the ridge projected onto s, with the accidental floor subtracted,
// drawn twice -- all cells, and with the edge strip dropped from both sides.
//
// Why drop strip 0. Its occupancy is 12% ABOVE the bulk while its coincidence
// fraction is about half the plateau, which is what an edge-truncated cluster
// looks like: reconstructed on the last strip, true position outside the
// silicon, no partner possible. Those hits sit at d=0 and drag the centroid
// down, and a PAIR carries the contamination twice.
void fwoPairPlot(const char *fn = "fwo_pair.root")
{
    gStyle->SetOptStat(0); gStyle->SetPalette(55);
    TFile *f = TFile::Open(fn);
    TH2F *h = (TH2F *) f->Get("hPair");
    if (!h) { printf("no hPair in %s\n", fn); return; }

    const int N = 12;
    double eff[N] = {1.120,1.070,1.088,1.069,0.983,1.001,0.981,0.985,0.928,0.972,0.965,0.963};
    double cor[N][N];
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
        cor[i][j] = h->GetBinContent(i+1, j+1) / (eff[i]*eff[j]);
    double bkg = 0; int nb = 0;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
        if (i+j >= 14) { bkg += cor[i][j]; nb++; }
    bkg /= nb;

    TCanvas *c = new TCanvas("c", "pair centroid", 1200, 520);
    c->Divide(2,1);

    c->cd(1); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.15); gPad->SetBottomMargin(0.13);
    h->SetTitle("OUTER pair map: the ridge is the width;"
                "my strip number from the wedge edge;partner strip number");
    h->GetXaxis()->SetRangeUser(-0.5, 7.5); h->GetYaxis()->SetRangeUser(-0.5, 7.5);
    h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
    h->Draw("COLZ");
    TLine *rg = new TLine(-0.5, 4.5, 4.5, -0.5);
    rg->SetLineColor(kRed); rg->SetLineWidth(3); rg->SetLineStyle(2); rg->Draw();
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.036); tl.SetTextColor(kRed);
    tl.DrawLatex(0.45, 0.78, "d_{me} + d_{nb} = 4");
    tl.SetTextColor(kBlack);
    tl.DrawLatex(0.17, 0.20, "five flat cells on the ridge:");
    tl.DrawLatex(0.17, 0.15, "4091 4184 4437 4195 4169");

    c->cd(2); gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    const int S = 9;
    double xs[S], yA[S], yB[S], eA[S], eB[S];
    double nA = 0, mA = 0, nB = 0, mB = 0;
    for (int s = 0; s < S; s++) {
        double a = 0, b = 0, ra = 0, rb = 0;
        for (int i = 0; i < N; i++) {
            int j = s - i; if (j < 0 || j >= N) continue;
            a += cor[i][j] - bkg;  ra += h->GetBinContent(i+1,j+1);
            if (i >= 1 && j >= 1) { b += cor[i][j] - bkg; rb += h->GetBinContent(i+1,j+1); }
        }
        xs[s] = s; yA[s] = a; yB[s] = b;
        eA[s] = (ra > 0) ? a/sqrt(ra) : 0;  eB[s] = (rb > 0) ? b/sqrt(rb) : 0;
        nA += a; mA += s*a; nB += b; mB += s*b;
    }
    mA /= nA; mB /= nB;

    TH2F *fr = new TH2F("fr", "the ridge projected: s = d_{me} + d_{nb};"
                              "s = d_{me} + d_{nb};coincidences, accidentals subtracted",
                        10, -0.5, 8.5, 10, 0, 21000);
    fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05);
    fr->Draw();
    TGraphErrors *gA = new TGraphErrors(S, xs, yA, 0, eA);
    gA->SetMarkerStyle(24); gA->SetMarkerSize(1.4); gA->SetMarkerColor(kGray+2);
    gA->SetLineColor(kGray+2); gA->SetLineWidth(2); gA->Draw("PL same");
    TGraphErrors *gB = new TGraphErrors(S, xs, yB, 0, eB);
    gB->SetMarkerStyle(20); gB->SetMarkerSize(1.5); gB->SetLineWidth(2); gB->Draw("PL same");

    TLegend *lg = new TLegend(0.47, 0.70, 0.97, 0.88);
    lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.037);
    lg->AddEntry(gB, Form("strip 0 dropped:  W = %.2f", mB+1), "pl");
    lg->AddEntry(gA, Form("all strips:          W = %.2f", mA+1), "pl");
    lg->Draw();

    TLatex t2; t2.SetNDC(); t2.SetTextSize(0.036);
    t2.DrawLatex(0.17, 0.60, "W = #LTs#GT + 1");
    t2.DrawLatex(0.17, 0.52, "yield-profile method on the");
    t2.DrawLatex(0.17, 0.47, "same detector: 5.20 - 5.35");
    t2.DrawLatex(0.17, 0.39, "two estimators with nothing");
    t2.DrawLatex(0.17, 0.34, "in common, agreeing");
    t2.DrawLatex(0.17, 0.17, "200,956 events");

    c->SaveAs("FstOuterSensorOverlap/plots/fwo_pair.png");
    printf("\n  all strips      <s> = %.3f   W = %.3f\n", mA, mA+1);
    printf("  strip 0 dropped <s> = %.3f   W = %.3f\n", mB, mB+1);
    printf("  accidental floor %.0f per cell\n", bkg);
    printf("wrote FstOuterSensorOverlap/plots/fwo_pair.png\n");
}
