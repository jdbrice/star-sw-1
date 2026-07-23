// plotFwdResidualQuadCharge.C
// Plots the quadrant x charge plane-usage split and per-disk hit-position
// maps added to StFwdResidualMaker on 2026-07-14 (see proposal_mixed.txt,
// "charge+quadrant hit-map added to StFwdResidualMaker") -- the direct test
// of the FTT-disk2-south-half dead-zone hypothesis for the dilepton
// north-charge-sign mass-shape asymmetry.
//
// Input: a merged (hadd -k) StFwdResidualMaker output file for ONE track
// type (e.g. merged_Global.root, merged_BLCVtx.root -- each StFwdResidualMaker
// instance only fills its own mResidualTrackType's quadrant/charge histos).
//
// Usage:
//   root4star -b -q 'plotFwdResidualQuadCharge.C("merged_Global.root","residual_qc/Global")'

const char* qcQuadName[4]  = {"Ntop", "Nbot", "Stop", "Sbot"};
const char* qcChargeName[2] = {"Pos", "Neg"};
const char* qcFttOriName[3] = {"x", "y", "uv"};
// bins 1-15 of the 20-bin scheme (bin0=Vtx skipped, bin16=EPD/17=ECal/18=HCal/19=AllTrk skipped)
const char* qcPlaneName[15] = {
    "FST1","FST2","FST3",
    "FTT1x","FTT1y","FTT1uv","FTT2x","FTT2y","FTT2uv",
    "FTT3x","FTT3y","FTT3uv","FTT4x","FTT4y","FTT4uv"
};
const int qcPlaneBin[15] = {1,2,3, 4,5,6,7,8,9, 10,11,12,13,14,15}; // 0-indexed bin content, AllTrk = bin 19

void mkdirp(const char* dir) { gSystem->mkdir(dir, kTRUE); }

// Sequential blue->green->red palette for the Pos/Neg count maps (panels
// 1-2). No white midpoint -- white-in-the-middle reads as "diverging /
// centered on zero", which is confusing for a plain event-count map.
// Deliberately NOT factored into a helper function: CreateGradientColorTable
// apparently doesn't reliably survive its local stops/red/green/blue arrays
// going out of scope when called from a separate function that returns
// before the subsequent Draw("colz") -- observed the diverging (asymmetry)
// palette bleeding into the Pos/Neg panels even after refactoring this into
// setSequentialPalette()/setDivergingPalette() helpers. Inlined directly at
// each call site instead (matching the original, empirically-working
// pattern) so the arrays stay in scope for the whole of plotXYAsymmetry().
#define SET_SEQUENTIAL_PALETTE() { \
    const int nStops = 3; \
    double stops[nStops] = {0.0, 0.5, 1.0}; \
    double red[nStops]   = {0.0, 0.0, 1.0}; \
    double green[nStops] = {0.0, 1.0, 0.0}; \
    double blue[nStops]  = {1.0, 0.0, 0.0}; \
    TColor::CreateGradientColorTable(nStops, stops, red, green, blue, 255); \
}

// Diverging red(-1)/white(0)/blue(+1) palette for the asymmetry panel --
// kLightTemperature/kBird enums don't exist in ROOT5.
#define SET_DIVERGING_PALETTE() { \
    const int nStops = 3; \
    double stops[nStops] = {0.0, 0.5, 1.0}; \
    double red[nStops]   = {0.0, 1.0, 1.0}; \
    double green[nStops] = {0.0, 1.0, 0.0}; \
    double blue[nStops]  = {1.0, 1.0, 0.0}; \
    TColor::CreateGradientColorTable(nStops, stops, red, green, blue, 255); \
}

// One canvas: 2x2 pads (quadrants), each pad = usage fraction vs plane for
// Pos (blue) and Neg (red) overlaid. Directly shows which quadrant+plane has
// a charge-dependent usage deficit.
void plotUsageFraction(TFile* f, const char* outdir, const char* tag) {
    TCanvas* c = new TCanvas(Form("c_usage_%s", tag), "", 1200, 900);
    c->Divide(2, 2);

    printf("\n=== Usage-fraction asymmetry summary (%s) ===\n", tag);
    printf("%-8s %-8s %8s %8s %8s %8s %6s\n", "Quad", "Plane", "Pos", "Neg", "PosErr", "NegErr", "Zscore");

    double worstZ = 0; TString worstLabel;

    for (int q = 0; q < 4; q++) {
        c->cd(q+1);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.22); gPad->SetGridy();

        TH1F* hPos = (TH1F*)f->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Pos", qcQuadName[q]));
        TH1F* hNeg = (TH1F*)f->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Neg", qcQuadName[q]));
        if (!hPos || !hNeg) { printf("WARN: missing histos for quadrant %s\n", qcQuadName[q]); continue; }

        double nPos = hPos->GetBinContent(20); // AllTrk
        double nNeg = hNeg->GetBinContent(20);

        TH1F* hFracPos = new TH1F(Form("hFracPos_%s_%s", tag, qcQuadName[q]), Form("%s;;Usage fraction", qcQuadName[q]), 15, 0, 15);
        TH1F* hFracNeg = new TH1F(Form("hFracNeg_%s_%s", tag, qcQuadName[q]), "", 15, 0, 15);

        for (int p = 0; p < 15; p++) {
            double nP = hPos->GetBinContent(qcPlaneBin[p] + 1);
            double nN = hNeg->GetBinContent(qcPlaneBin[p] + 1);
            double fP = (nPos > 0) ? nP / nPos : 0;
            double fN = (nNeg > 0) ? nN / nNeg : 0;
            double eP = (nPos > 0) ? sqrt(fP * (1 - fP) / nPos) : 0;
            double eN = (nNeg > 0) ? sqrt(fN * (1 - fN) / nNeg) : 0;
            hFracPos->SetBinContent(p+1, fP); hFracPos->SetBinError(p+1, eP);
            hFracNeg->SetBinContent(p+1, fN); hFracNeg->SetBinError(p+1, eN);
            hFracPos->GetXaxis()->SetBinLabel(p+1, qcPlaneName[p]);

            double sigma = sqrt(eP*eP + eN*eN);
            double z = (sigma > 0) ? (fP - fN) / sigma : 0;
            printf("%-8s %-8s %8.4f %8.4f %8.4f %8.4f %6.2f\n", qcQuadName[q], qcPlaneName[p], fP, fN, eP, eN, z);
            if (fabs(z) > fabs(worstZ)) { worstZ = z; worstLabel = Form("%s/%s", qcQuadName[q], qcPlaneName[p]); }
        }

        hFracPos->SetLineColor(kBlue+1); hFracPos->SetMarkerColor(kBlue+1); hFracPos->SetMarkerStyle(20);
        hFracNeg->SetLineColor(kRed+1);  hFracNeg->SetMarkerColor(kRed+1);  hFracNeg->SetMarkerStyle(21);
        hFracPos->GetXaxis()->LabelsOption("v");
        hFracPos->SetMaximum(1.05);
        hFracPos->SetMinimum(0);
        hFracPos->SetStats(0);
        hFracPos->SetTitle(Form("%s (Ntrk: +%.0f / -%.0f)", qcQuadName[q], nPos, nNeg));
        hFracPos->Draw("E1");
        hFracNeg->Draw("E1 same");

        TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->AddEntry(hFracPos, "q>0", "lp");
        leg->AddEntry(hFracNeg, "q<0", "lp");
        leg->SetBorderSize(0);
        leg->Draw();
    }
    c->Print(Form("%s_usageFraction.png", outdir));
    printf("=== Largest |Z| asymmetry: %s  Z=%.2f ===\n\n", worstLabel.Data(), worstZ);
    delete c;
}

// For each of the 15 disk/orientation combos: 3-panel canvas (Pos, Neg,
// asymmetry = (Pos-Neg)/(Pos+Neg) with a symmetric red/blue palette). The
// asymmetry map is what would visually reveal a charge-dependent dead
// region (e.g. FTT disk2 south-bottom).
//
// dirName/histPrefix/pngTag select which set of 2D maps to use:
//   "HitXYByCharge"/"h2FstHitXY"/"h2FttHitXY" -- (x,y) of the hit actually
//     used, i.e. strip-quantized along the strip's long axis (banding).
//   "ProjXYByCharge"/"h2FstProjXY"/"h2FttProjXY" -- (x,y) of the track's
//     continuous projection at that disk, gated on a hit having been used
//     there. This is the one to use for spatial dead-zone hunting -- added
//     2026-07-15 per the user's suggestion, since the hit-position maps'
//     strip-axis quantization made real spatial structure hard to see.
void plotXYAsymmetry(TFile* f, const char* outdir, const char* tag,
                      const char* dirName, const char* fstPrefix, const char* fttPrefix,
                      const char* pngTag) {
    TString names[15];
    int idx = 0;
    for (int d = 0; d < 3; d++) names[idx++] = Form("%s_disk%d", fstPrefix, d);
    for (int d = 0; d < 4; d++)
        for (int o = 0; o < 3; o++)
            names[idx++] = Form("%s_disk%d_%s", fttPrefix, d, qcFttOriName[o]);

    for (int i = 0; i < 15; i++) {
        TH2F* hPos = (TH2F*)f->Get(Form("%s/%s_Pos", dirName, names[i].Data()));
        TH2F* hNeg = (TH2F*)f->Get(Form("%s/%s_Neg", dirName, names[i].Data()));
        if (!hPos || !hNeg) { printf("WARN: missing %s histos for %s\n", dirName, names[i].Data()); continue; }
        if (hPos->GetEntries() < 20 && hNeg->GetEntries() < 20) continue; // skip empty (e.g. unused seed source)

        // FST (i<3) is a much smaller detector than FTT (i>=3) -- zoom the
        // displayed axis range accordingly. This only changes the drawn
        // range, not the underlying binning/data (booked out to +-65cm in
        // StFwdResidualMaker.cxx for both).
        double range = (i < 3) ? 30.0 : 75.0;

        TH2F* hAsym = (TH2F*)hPos->Clone(Form("hAsym_%s_%s", tag, names[i].Data()));
        hAsym->Reset();
        int nbx = hPos->GetNbinsX(), nby = hPos->GetNbinsY();
        for (int bx = 1; bx <= nbx; bx++) {
            for (int by = 1; by <= nby; by++) {
                double p = hPos->GetBinContent(bx, by), n = hNeg->GetBinContent(bx, by);
                double sum = p + n;
                if (sum < 5) continue; // not enough stats in this bin to say anything
                hAsym->SetBinContent(bx, by, (p - n) / sum);
            }
        }

        gStyle->SetOptStat(0);

        // ROOT only has ONE active global colz palette at a time -- with a
        // single TCanvas::Divide() and multiple pads, the LAST
        // CreateGradientColorTable() call before Print() silently wins for
        // EVERY pad's colz rendering, regardless of which palette was
        // "current" when each pad's Draw() was called (confirmed with a
        // minimal reproducer: the palette applied at TCanvas::Print()'s
        // final paint pass, not at each pad's own Draw() time -- TExec
        // per-histogram doesn't help either, since TExec only fires on
        // interactive mouse events, never during batch Print()). So the
        // Pos/Neg panels were always silently getting the diverging
        // asymmetry-panel palette. Work around this by rendering each panel
        // to its own canvas/PNG (each gets an uncontested Print() with the
        // correct palette already set) and stitching the three side by side
        // with ImageMagick.
        int W = 500, H = 500;
        TString tmpPos = Form("%s_%s_%s_tmpPos.png", outdir, pngTag, names[i].Data());
        TString tmpNeg = Form("%s_%s_%s_tmpNeg.png", outdir, pngTag, names[i].Data());
        TString tmpAsym = Form("%s_%s_%s_tmpAsym.png", outdir, pngTag, names[i].Data());

        TCanvas* cPos = new TCanvas(Form("cPos_%s_%s", tag, names[i].Data()), "", W, H);
        gPad->SetRightMargin(0.15);
        gPad->SetLogz();
        SET_SEQUENTIAL_PALETTE();
        hPos->GetXaxis()->SetRangeUser(-range, range);
        hPos->GetYaxis()->SetRangeUser(-range, range);
        hPos->SetTitle(Form("%s q>0;x [cm];y [cm]", names[i].Data()));
        hPos->Draw("colz");
        cPos->Print(tmpPos);
        delete cPos;

        TCanvas* cNeg = new TCanvas(Form("cNeg_%s_%s", tag, names[i].Data()), "", W, H);
        gPad->SetRightMargin(0.15);
        gPad->SetLogz();
        SET_SEQUENTIAL_PALETTE();
        hNeg->GetXaxis()->SetRangeUser(-range, range);
        hNeg->GetYaxis()->SetRangeUser(-range, range);
        hNeg->SetTitle(Form("%s q<0;x [cm];y [cm]", names[i].Data()));
        hNeg->Draw("colz");
        cNeg->Print(tmpNeg);
        delete cNeg;

        // NOTE: no SetLogz() here -- asymmetry spans [-1,1] including 0/negative,
        // a log scale doesn't apply to it.
        TCanvas* cAsym = new TCanvas(Form("cAsym_%s_%s", tag, names[i].Data()), "", W, H);
        gPad->SetRightMargin(0.15);
        SET_DIVERGING_PALETTE();
        hAsym->GetXaxis()->SetRangeUser(-range, range);
        hAsym->GetYaxis()->SetRangeUser(-range, range);
        hAsym->SetMinimum(-1); hAsym->SetMaximum(1);
        hAsym->SetTitle(Form("%s asymmetry (P-N)/(P+N);x [cm];y [cm]", names[i].Data()));
        hAsym->Draw("colz");
        cAsym->Print(tmpAsym);
        delete cAsym;

        TString finalPng = Form("%s_%s_%s.png", outdir, pngTag, names[i].Data());
        gSystem->Exec(Form("convert +append %s %s %s %s", tmpPos.Data(), tmpNeg.Data(), tmpAsym.Data(), finalPng.Data()));
        gSystem->Exec(Form("rm -f %s %s %s", tmpPos.Data(), tmpNeg.Data(), tmpAsym.Data()));

        delete hAsym;
    }
}

void plotFwdResidualQuadCharge(const char* fname, const char* outdirbase) {
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) { printf("ERROR: cannot open %s\n", fname); return; }

    TString outdir(outdirbase);
    Ssiz_t slash = outdir.Last('/');
    if (slash != kNPOS) mkdirp(outdir(0, slash).Data());

    gErrorIgnoreLevel = kWarning; // silence routine TCanvas::Print info lines

    plotUsageFraction(f, outdirbase, outdirbase);
    plotXYAsymmetry(f, outdirbase, outdirbase, "HitXYByCharge", "h2FstHitXY", "h2FttHitXY", "hitXY");
    plotXYAsymmetry(f, outdirbase, outdirbase, "ProjXYByCharge", "h2FstProjXY", "h2FttProjXY", "projXY");

    f->Close();
}
