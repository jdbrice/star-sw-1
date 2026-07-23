// checkAllExcessQuad.C
//
// Quadrant-split version of checkAllExcess.C -- tests whether the
// getOrientation() H/V swap (StFttDb.cxx, 2026-07-21) helps specifically in
// quadrants A/C (odd rob, independently confirmed backward via online QA
// for disk1/quadA) versus B/D (even rob, swapped only by assumption, not
// independently confirmed).
//
// Usage: root4star -l -b -q 'checkAllExcessQuad.C("fwd_blind_diag.root")'

void checkAllExcessQuad(const char* infile = "fwd_blind_diag.root") {
    TFile *f = TFile::Open(infile);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", infile); return; }

    const char* qName[4] = {"A", "B", "C", "D"};
    const char* oriName[2] = {"V", "H"};

    for (int d = 0; d < 4; d++) {
        for (int q = 0; q < 4; q++) {
            for (int io = 0; io < 2; io++) {
                TString histName = (io == 0) ? Form("hBlindDxAll_V_disk%d_quad%s", d, qName[q])
                                              : Form("hBlindDyAll_H_disk%d_quad%s", d, qName[q]);
                TH1F *h = (TH1F*)f->Get(histName);
                if (!h) { printf("Missing %s\n", histName.Data()); continue; }

                double loSum = 0; int loN = 0;
                double hiSum = 0; int hiN = 0;
                double peakSum = 0; int peakN = 0;
                for (int ib = 1; ib <= h->GetNbinsX(); ib++) {
                    double x = h->GetXaxis()->GetBinCenter(ib);
                    if (x >= -3.0 && x < -0.6) { loSum += h->GetBinContent(ib); loN++; }
                    if (x >= 0.6 && x < 3.0)   { hiSum += h->GetBinContent(ib); hiN++; }
                    if (fabs(x) < 0.6)          { peakSum += h->GetBinContent(ib); peakN++; }
                }
                double lo = loSum / loN;
                double hi = hiSum / hiN;
                double bg = (loSum + hiSum) / (loN + hiN);
                double bgErr = sqrt(loSum + hiSum) / (loN + hiN);
                double peak = peakSum / peakN;
                double excess = (bg > 0) ? peak / bg : 0;
                double sigma = (bgErr > 0) ? (peak - bg) / bgErr : 0;

                printf("disk%d quad%s %s: bg=%.1f+/-%.1f (lo=%.1f hi=%.1f)  peak=%.1f  excess=%.2fx  sig=%.1fsig\n",
                       d, qName[q], oriName[io], bg, bgErr, lo, hi, peak, excess, sigma);
            }
        }
        printf("\n");
    }
}
