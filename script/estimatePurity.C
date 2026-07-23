// estimatePurity.C
//
// Estimates the PURITY of FTT hit matching (of the hits assigned to a
// track, what fraction are the genuinely correct one) directly from the
// FST-blind excess check.
//
// Logic: for a given track-disk search, the real (correct) FTT hit -- if
// present and passing the loose off-axis-conditioned gate -- contributes
// AT MOST one entry to the "All candidates" histogram (it's a single
// physical hit object encountered once in the search loop), landing in the
// narrow peak at dx/dy=0 (sub-mm strip precision). So the background-
// subtracted integral of that peak directly counts "track-searches where
// the real hit was present and passed the loose gate" -- call this N_real.
// Since the real hit trivially also passes the tighter precise-coordinate
// gate (it's AT zero) and wins "closest" selection (nothing else is closer
// than the true hit), N_real is a good proxy for "how many of the entries
// in the Matched histogram are the genuinely correct hit", i.e.
//   purity = N_real / N_matched_total
// This is NOT full efficiency (that needs a denominator of "tracks that
// should have had a real hit available", which this diagnostic alone
// doesn't give -- see the writeup).
//
// Usage: root4star -l -b -q 'estimatePurity.C("fwd_blind_diag.root")'

void estimatePurity(const char* infile = "fwd_blind_diag.root") {
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", infile); return; }

    const char* oriName[2] = {"V", "H"};
    double sumReal = 0, sumMatched = 0;

    for (int d = 0; d < 4; d++) {
        for (int io = 0; io < 2; io++) {
            TString allName = (io == 0) ? Form("hBlindDxAll_V_disk%d", d) : Form("hBlindDyAll_H_disk%d", d);
            TString matName = (io == 0) ? Form("hBlindDxMatched_V_disk%d", d) : Form("hBlindDyMatched_H_disk%d", d);
            TH1F *hAll = (TH1F*)f->Get(allName);
            TH1F *hMat = (TH1F*)f->Get(matName);
            if (!hAll || !hMat) { printf("Missing %s or %s\n", allName.Data(), matName.Data()); continue; }

            // local background from sidebands just outside the +/-1.0cm peak window
            double loSum = 0; int loN = 0;
            double hiSum = 0; int hiN = 0;
            for (int ib = 1; ib <= hAll->GetNbinsX(); ib++) {
                double x = hAll->GetXaxis()->GetBinCenter(ib);
                if (x >= -3.0 && x < -1.0) { loSum += hAll->GetBinContent(ib); loN++; }
                if (x >= 1.0 && x < 3.0)   { hiSum += hAll->GetBinContent(ib); hiN++; }
            }
            double bgPerBin = (loSum + hiSum) / (loN + hiN);

            // integrate background-subtracted excess over +/-1.0cm (wide enough
            // to safely capture the full peak width seen in the plots)
            double nReal = 0;
            for (int ib = 1; ib <= hAll->GetNbinsX(); ib++) {
                double x = hAll->GetXaxis()->GetBinCenter(ib);
                if (fabs(x) < 1.0) nReal += (hAll->GetBinContent(ib) - bgPerBin);
            }
            double nMatched = hMat->Integral();
            double purity = (nMatched > 0) ? nReal / nMatched : 0;

            printf("disk%d %s: N_real=%.0f  N_matched=%.0f  purity=%.3f (%.1f%%)\n",
                   d, oriName[io], nReal, nMatched, purity, 100*purity);

            sumReal += nReal;
            sumMatched += nMatched;
        }
    }
    printf("\nOVERALL (sum over 4 disks x 2 orientations): N_real=%.0f  N_matched=%.0f  purity=%.3f (%.1f%%)\n",
           sumReal, sumMatched, sumReal/sumMatched, 100*sumReal/sumMatched);
}
