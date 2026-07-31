// checkFstPhi16.C
//
// Real data shows a sharp spike at exactly nRawHitsPhi = 16 (FstSlowSim/fstclus/index.html):
// 0.115%/0.164% of clusters, against a ~0.018% plateau at 9-15 and 0.010% at 17.
//
// nRawHitsPhi is the length of a contiguous run of occupied phi strips (the
// cluster maker's step-2 merge condition is unconditionally true, so adjacent
// phi strips always chain). From Calibrations/fst/fstMapping, one APV chip
// covers 32 phi strips x 4 r strips = 128 channels, so 16 phi strips is exactly
// HALF a chip = 64 channels.
//
// If the spike is a readout-block effect, those 16-strip runs must be ALIGNED to
// fixed 16-strip boundaries, so the cluster's meanPhiStrip (mod 16) piles up at
// one value. If instead they are random coincidences, meanPhiStrip mod 16 is
// flat. That is the test here.
//
// Usage:
//   root4star -b -q 'script/checkFstPhi16.C("<glob>", 0, "FstSlowSim/fstclus")'

void checkFstPhi16(const char* glob, int nEvents = 0, const char* outdir = "FstSlowSim/fstclus"){

    gSystem->mkdir(outdir, kTRUE);
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    gStyle->SetOptStat(0);

    TChain* ch = new TChain("MuDst");
    if (ch->Add(glob) <= 0){ printf("no files matched %s\n", glob); return; }

    TClonesArray* fstArr = 0;
    ch->SetBranchStatus("*", 0);
    ch->SetBranchStatus("FstHit*", 1);
    ch->SetBranchAddress("FstHit", &fstArr);

    TH1F* hMod16 = new TH1F("hMod16", ";meanPhiStrip mod 16;clusters with nRawHitsPhi = 16",
                            16, -0.5, 15.5);
    TH1F* hMod16c = new TH1F("hMod16c", ";meanPhiStrip mod 16;clusters with nRawHitsPhi = 1",
                             16, -0.5, 15.5);
    TH1F* hPhi16 = new TH1F("hPhi16", ";meanPhiStrip;clusters with nRawHitsPhi = 16",
                            128, -0.5, 127.5);
    TH1F* hApv16 = new TH1F("hApv16", ";apv;clusters with nRawHitsPhi = 16", 8, -0.5, 7.5);
    hMod16->SetDirectory(0); hMod16c->SetDirectory(0);
    hPhi16->SetDirectory(0); hApv16->SetDirectory(0);

    double n16 = 0, n1 = 0;

    Long64_t ne = ch->GetEntries();
    if (nEvents > 0 && nEvents < ne) ne = nEvents;
    printf("events to read: %lld\n", ne);

    for (Long64_t iev = 0; iev < ne; iev++){
        ch->GetEntry(iev);
        if (!fstArr) continue;
        int nh = fstArr->GetEntriesFast();
        for (int ih = 0; ih < nh; ih++){
            StMuFstHit* hit = (StMuFstHit*) fstArr->UncheckedAt(ih);
            if (!hit) continue;
            int nrp = (int) hit->getNRawHitsPhi();
            int mps = (int) hit->getMeanPhiStrip();
            if (mps < 0 || mps > 127) continue;

            if (nrp == 16){
                hMod16->Fill(mps % 16);
                hPhi16->Fill(mps);
                hApv16->Fill((int) hit->getApv() % 8);
                n16 += 1;
            } else if (nrp == 1){
                hMod16c->Fill(mps % 16);
                n1 += 1;
            }
        }
    }

    printf("\n===============================================================\n");
    printf(" Alignment test for the nRawHitsPhi = 16 spike\n");
    printf("===============================================================\n");
    printf("  clusters with nRawHitsPhi = 16 : %.0f\n", n16);
    printf("  clusters with nRawHitsPhi =  1 : %.0f   (flat-reference control)\n", n1);
    printf("\n  meanPhiStrip mod 16 |   nPhi=16      nPhi=1 (control)\n");
    for (int im = 0; im < 16; im++){
        double a = hMod16->GetBinContent(im+1);
        double b = hMod16c->GetBinContent(im+1);
        printf("        %2d           |  %7.3f%%      %7.3f%%\n", im,
               n16>0?100.0*a/n16:0, n1>0?100.0*b/n1:0);
    }
    printf("\n  A flat 6.25%% in both columns means the 16-strip runs are NOT\n");
    printf("  aligned to a hardware block; a single dominant row means they are.\n");

    printf("\n  apv distribution of the nPhi=16 clusters:\n");
    for (int ja = 0; ja < 8; ja++)
        printf("    apv %d : %6.2f%%\n", ja, n16>0?100.0*hApv16->GetBinContent(ja+1)/n16:0);

    TCanvas* cv = new TCanvas("c_phi16", "", 1000, 400);
    cv->Divide(2,1);
    cv->cd(1);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    if (n16 > 0) hMod16->Scale(1.0/n16);
    if (n1  > 0) hMod16c->Scale(1.0/n1);
    hMod16->SetLineColor(kRed+1);  hMod16->SetLineWidth(2); hMod16->SetFillStyle(0);
    hMod16c->SetLineColor(kBlack); hMod16c->SetLineWidth(2); hMod16c->SetFillStyle(0);
    hMod16->SetMinimum(0);
    hMod16->SetTitle("phi alignment of the 16-strip clusters");
    hMod16->GetYaxis()->SetTitle("fraction of clusters");
    hMod16->Draw("hist"); hMod16c->Draw("hist same");
    TLegend* lg = new TLegend(0.42,0.74,0.88,0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(hMod16,  "nRawHitsPhi = 16", "l");
    lg->AddEntry(hMod16c, "nRawHitsPhi = 1 (control)", "l");
    lg->Draw();

    cv->cd(2);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
    hPhi16->SetLineColor(kRed+1); hPhi16->SetLineWidth(1); hPhi16->SetFillStyle(0);
    hPhi16->SetTitle("meanPhiStrip of the 16-strip clusters");
    hPhi16->Draw("hist");
    cv->SaveAs(Form("%s/fstPhi16.png", outdir));
    printf("\nwrote %s/fstPhi16.png\n", outdir);
}
