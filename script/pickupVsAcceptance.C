// pickupVsAcceptance.C -- separate "the sTGC is not there" from "the sTGC is there
// and we missed it".
//
// The sTGC has a design acceptance gap below the beampipe (and the four pentagons
// cluster above and to the sides of it), while FST covers nearly 2pi. So FST-seeded
// tracks projecting into the gap can never pick up an sTGC hit, and a raw pickup
// probability is diluted by geometry that is working exactly as designed.
//
//   numerator   h2FttProjXY_disk<d>_<ori>_<q>  projection of tracks that DID match
//   denominator h2FttProjAllXY_disk<d>         projection of EVERY track
//
// Both are binned at the PROJECTION position, so the per-cell ratio is a real
// efficiency. "Active" cells are those where matching demonstrably happens at all
// (numerator above a floor), i.e. where the detector is; restricting to those removes
// the design gap from the comparison. If top and bottom then agree, the deficit is
// acceptance. If bottom is still well below top, something is wrong in B/C.
//
// CINT: unique loop variable names, no std::map (see CLAUDE.md).

void pickupVsAcceptance(const char* fname = "acc_Primary.root", double ymid = 0.0,
                        int rebin = 4, double minDen = 20.0, double minNum = 3.0) {
    TFile* fa = TFile::Open(fname);
    if (!fa || fa->IsZombie()) { printf("cannot open %s\n", fname); return; }

    const char* on[2] = { "x", "y" };
    const char* cn[2] = { "Pos", "Neg" };
    const char* oname[2] = { "x (V strips)", "y (H strips)" };

    printf("\n=== pickup INSIDE the active area, top vs bottom ===\n");
    printf("  %s, rebin %d (cell %.1f cm), active cell: den>%.0f and num>%.0f\n",
           fname, rebin, 1.625*rebin, minDen, minNum);
    printf("  orientations kept separate: a track can match a V and an H hit on the\n");
    printf("  same plane, so summing them can exceed 1 per track.\n\n");
    printf("  %-4s %-14s %8s %7s   %8s %7s   %s\n",
           "disk","orientation","TOPpick","cells","BOTpick","cells","bot/top");

    for (int id = 0; id < 4; id++) {
        TH2* den0 = (TH2*)fa->Get(Form("ProjXYByCharge/h2FttProjAllXY_disk%d", id));
        if (!den0) { printf("  disk %d: no denominator -- rerun with the new library\n", id); continue; }
        TH2* den = (TH2*)den0->Clone(Form("dn%d", id)); den->SetDirectory(0); den->Rebin2D(rebin, rebin);

        for (int io = 0; io < 2; io++) {
            TH2* num = 0;
            for (int ic = 0; ic < 2; ic++) {
                TH2* t = (TH2*)fa->Get(Form("ProjXYByCharge/h2FttProjXY_disk%d_%s_%s", id, on[io], cn[ic]));
                if (!t) continue;
                if (!num) { num = (TH2*)t->Clone(Form("nm%d_%d", id, io)); num->SetDirectory(0); }
                else num->Add(t);
            }
            if (!num) continue;
            num->Rebin2D(rebin, rebin);

            double nT = 0, dT = 0, nB = 0, dB = 0; int cT = 0, cB = 0;
            for (int ix = 1; ix <= den->GetNbinsX(); ix++) {
                for (int iy = 1; iy <= den->GetNbinsY(); iy++) {
                    double dd = den->GetBinContent(ix, iy);
                    double nn = num->GetBinContent(ix, iy);
                    if (dd < minDen || nn < minNum) continue;   // active cells only
                    double yc = den->GetYaxis()->GetBinCenter(iy);
                    if (yc > ymid) { nT += nn; dT += dd; cT++; }
                    else            { nB += nn; dB += dd; cB++; }
                }
            }
            double pT = (dT > 0) ? nT / dT : 0, pB = (dB > 0) ? nB / dB : 0;
            printf("   %d   %-14s %8.4f %7d   %8.4f %7d   %s\n", id, oname[io], pT, cT, pB, cB,
                   (pT > 0 && pB > 0) ? Form("%.2f", pB / pT) : "  --");
        }
    }
    printf("\n  bot/top ~1  => the raw deficit is the DESIGN acceptance gap, B/C are fine\n");
    printf("  bot/top <<1 => the sTGC is there and being missed: a real B/C problem\n");
}
