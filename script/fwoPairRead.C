// fwoPairRead.C -- the pair-centroid measurement of the FST overlap width.
//
//   root -l -b -q 'script/fwoPairRead.C("fwo_pair.root")'
//
// THE OBSERVABLE (akio's, 2026-11-10). Fold both ends of the wedge to a
// distance from its edge: d = c for the low end, d = 127 - c for the high end
// (equivalently 2*127 - i - j on the high side). The two strip grids run
// OPPOSITE ways through the overlap band, so every pair of strips at the same
// physical angle satisfies
//
//     s = d_me + d_nb = constant,     W = <s> + 1
//
// The sum is therefore a direct measure of the overlap width, and its centroid
// is much sharper than an edge read off a falling yield profile: a centroid
// does not care about the overall normalisation, only about the symmetry of
// the ridge.
//
// WHY A 12-SLOT WINDOW. With W ~ 5.3 strips, my strip 0 pairs with the
// neighbour's strip ~4. The original hFine kept only strips 0-3, so that cell
// was missing and the s = 4 ridge was sampled at three of its five cells,
// asymmetrically -- which biases the centroid low. Five or six a side is the
// minimum; 12 also buys a clean region to measure the accidental floor in.
//
// THE TWO CORRECTIONS, in the order they matter.
//
//  1. ACCIDENTALS. Cells far from the ridge are same-event combinatorics, not
//     overlap. Measured here from d_me + d_nb >= kBkgFrom and subtracted as a
//     flat level AFTER the efficiency correction (accidentals scale as
//     occ(i)*occ(j), so dividing by eff(i)*eff(j) flattens them).
//
//  2. PER-STRIP EFFICIENCY, APPLIED TWICE. Strips 0 and 127 are visibly less
//     efficient than their neighbours, and a PAIR carries the suppression
//     twice -- mine and the partner's -- so the weight is eff(i)*eff(j), not
//     a single factor. eff comes from hPhiStripPerSensor, the per-strip
//     occupancy, normalised to the flat bulk of the wedge; the assumption is
//     that illumination is uniform in phi across one wedge, which is what
//     makes a departure at the edge an efficiency and not a real gradient.
//
// THE INNER IS THE CASE THAT NEEDS THIS. Its overlap is 1-2 strips, so it
// lives almost entirely IN strips 0 and 1 -- the suppressed ones. The outer
// can afford to drop them and lose little; the inner cannot, because dropping
// them removes the signal itself. So for the inner the efficiency correction
// is not a refinement, it is the measurement. Both routes are printed below
// (drop the edge strips, or correct for them); if they disagree, the uniform
// illumination assumption is what is wrong.

const int kN        = 12;   // slots per side in the pair map
const int kBkgFrom  = 14;   // s at or above this is taken as pure accidental

void fwoPairOne(TFile *f, const char *hname, const char *occName,
                int sensorLo, int sensorHi, const char *label)
{
    TH2F *h = (TH2F *) f->Get(hname);
    if (!h) { printf("%s: absent\n", hname); return; }
    TH2F *hp = (TH2F *) f->Get(occName);
    if (!hp) { printf("%s: absent\n", occName); return; }

    // ---- per-strip efficiency, from the occupancy, normalised to the bulk --
    // hPhiStripPerSensor x axis is 132 bins over -2..130, so strip c is bin c+3.
    double occ[128];
    for (int c = 0; c < 128; c++) {
        double v = 0;
        for (int sv = sensorLo; sv <= sensorHi; sv++)
            v += hp->GetBinContent(c + 3, hp->GetYaxis()->FindBin(sv));
        occ[c] = v;
    }
    double bulk = 0; int nb = 0;
    for (int c = 20; c <= 107; c++) { bulk += occ[c]; nb++; }
    bulk /= nb;
    // fold: slot d sees strip d at one end and 127-d at the other
    double eff[kN];
    for (int d = 0; d < kN; d++) eff[d] = (occ[d] + occ[127 - d]) / (2.0 * bulk);

    printf("\n================ %s  (%s)\n", label, hname);
    printf("per-strip efficiency from the occupancy, bulk = 1:\n   d   :");
    for (int d = 0; d < kN; d++) printf(" %6d", d);
    printf("\n   eff :");
    for (int d = 0; d < kN; d++) printf(" %6.3f", eff[d]);
    printf("\n");

    // ---- the pair map, raw and corrected ---------------------------------
    double raw[kN][kN], cor[kN][kN];
    for (int i = 0; i < kN; i++)
        for (int j = 0; j < kN; j++) {
            raw[i][j] = h->GetBinContent(i + 1, j + 1);
            double w = eff[i] * eff[j];
            cor[i][j] = (w > 0) ? raw[i][j] / w : 0;
        }

    printf("\nraw pair map, my d (down) against partner d (across):\n      ");
    for (int j = 0; j < kN; j++) printf(" %7d", j);
    printf("\n");
    for (int i = 0; i < kN; i++) {
        printf("   %2d ", i);
        for (int j = 0; j < kN; j++) printf(" %7.0f", raw[i][j]);
        printf("\n");
    }

    // ---- accidental floor, measured where there is no ridge ---------------
    double bkg = 0; int nbk = 0;
    for (int i = 0; i < kN; i++)
        for (int j = 0; j < kN; j++)
            if (i + j >= kBkgFrom) { bkg += cor[i][j]; nbk++; }
    bkg = nbk ? bkg / nbk : 0;

    // ---- the sum distribution, three ways ---------------------------------
    // A: raw, everything            B: raw, both d >= 1 (edge strips dropped)
    // C: efficiency corrected, accidental subtracted
    double sA[2 * kN], sB[2 * kN], sC[2 * kN];
    for (int k = 0; k < 2 * kN; k++) { sA[k] = sB[k] = sC[k] = 0; }
    for (int i = 0; i < kN; i++)
        for (int j = 0; j < kN; j++) {
            sA[i + j] += raw[i][j];
            if (i >= 1 && j >= 1) sB[i + j] += raw[i][j];
            double v = cor[i][j] - bkg;
            if (v > 0) sC[i + j] += v;
        }

    printf("\naccidental floor after efficiency correction: %.0f per cell"
           "  (from s >= %d)\n", bkg, kBkgFrom);
    printf("\n s          :");
    for (int k = 0; k < 10; k++) printf(" %8d", k);
    printf("\n raw        :");
    for (int k = 0; k < 10; k++) printf(" %8.0f", sA[k]);
    printf("\n d>=1 only  :");
    for (int k = 0; k < 10; k++) printf(" %8.0f", sB[k]);
    printf("\n eff-corr   :");
    for (int k = 0; k < 10; k++) printf(" %8.0f", sC[k]);
    printf("\n");

    // The centroid must be taken over the PEAK, with the floor subtracted and
    // WITHOUT clipping. Averaging over the whole 0..22 window instead gives a
    // number dominated by the accidental sea -- it came out 6.09 for the outer
    // where the ridge is plainly at 4 -- and clipping negative cells at zero
    // biases it upward on top of that. Restricting to s <= kSMax keeps the
    // ridge and its shoulders and drops the sea.
    const int kSMax = 8;
    printf("\n   %-28s %8s %8s %7s\n", "", "<s>", "W=<s>+1", "rms");
    for (int t = 0; t < 2; t++) {
        int lo = t;                       // t=0 all d, t=1 drop the edge strip
        double num = 0, den = 0;
        for (int sv = 0; sv <= kSMax; sv++)
            for (int i = lo; i < kN; i++) {
                int j = sv - i;
                if (j < lo || j >= kN) continue;
                double v = cor[i][j] - bkg;
                num += sv * v; den += v;
            }
        if (den <= 0) continue;
        double m = num / den, var = 0;
        for (int sv = 0; sv <= kSMax; sv++)
            for (int i = lo; i < kN; i++) {
                int j = sv - i;
                if (j < lo || j >= kN) continue;
                var += (cor[i][j] - bkg) * (sv - m) * (sv - m);
            }
        printf("   %-28s %8.3f %8.3f %7.2f\n",
               t ? "eff-corrected, drop d=0" : "eff-corrected, all d",
               m, m + 1, sqrt(var / den));
    }
    printf("\n   A centroid is only unbiased if the ridge sits well inside the\n"
           "   window on BOTH sides. For the outer, ridge at s=4 in a 0-11\n"
           "   window, it does. For the inner the ridge is at s~0, hard against\n"
           "   the s>=0 wall, so <s> is pushed up by the tail and cannot be\n"
           "   turned into a width -- see the note in the page.\n");
}

void fwoPairRead(const char *fn = "fwo_pair.root")
{
    TFile *f = TFile::Open(fn);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", fn); return; }
    printf("W = <s> + 1, with s = d_me + d_nb and d folded from both wedge ends.\n");
    fwoPairOne(f, "hPair",  "hPhiStripPerSensor", 1, 2, "OUTER sensors");
    fwoPairOne(f, "hPairI", "hPhiStripPerSensor", 0, 0, "INNER sensor");
    f->Close();
}
