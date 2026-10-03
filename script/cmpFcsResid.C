// cmpFcsResid.C -- FCS cluster-to-track residual, with and without aligned FTT hits
// on the track.  Two picoMatch outputs that differ ONLY in the afterburner's NOADD
// flag (same Run24 iteration-1 geometry, same 90 input files, same 647117 events).
//
// picoMatch books a "Same-Mixed" slot but never fills it, so the mixed-event
// subtraction is done here: mixed is normalised to same in the tails, where no real
// correlation lives, then subtracted.  The peak that survives is the real match.
//
// CINT: unique loop variable names, fixed arrays, no std::map (see CLAUDE.md).

double frMu, frSig, frErrMu, frErrSig, frSignif, frEnt, frRawRms;

// normalise mixed to same over |v| in [lo,hi] and subtract
TH1D* frSubtract(TH1* same, TH1* mixed, double lo, double hi) {
    if (!same || !mixed) return 0;
    TH1D* d = (TH1D*)same->Clone(Form("%s_sub", same->GetName()));
    d->SetDirectory(0);
    double sS = 0, sM = 0;
    for (int ib = 1; ib <= same->GetNbinsX(); ib++) {
        double x = fabs(same->GetXaxis()->GetBinCenter(ib));
        if (x < lo || x > hi) continue;
        sS += same->GetBinContent(ib);
        sM += mixed->GetBinContent(ib);
    }
    double k = (sM > 0) ? sS / sM : 1.0;
    TH1D* m = (TH1D*)mixed->Clone(Form("%s_scl", mixed->GetName()));
    m->SetDirectory(0);
    m->Scale(k);
    d->Add(m, -1.0);
    return d;
}

void frFit(TH1* h) {
    frMu = 0; frSig = 0; frErrMu = 0; frErrSig = 0; frSignif = 0; frEnt = 0;
    if (!h) return;
    frEnt = h->Integral();
    int pb = h->GetMaximumBin();
    double pk = h->GetXaxis()->GetBinCenter(pb);
    double w = 30.0;                              // the real peak is a few cm; fit wide
    TF1* fn = new TF1(Form("f_%s", h->GetName()), "gaus(0)+pol1(3)", pk - w, pk + w);
    fn->SetParameters(h->GetBinContent(pb), pk, 5.0, 0, 0);
    fn->SetParLimits(1, pk - w, pk + w);
    fn->SetParLimits(2, 0.5, 60.0);
    h->Fit(fn, "QNR");
    double A = fn->GetParameter(0), eA = fn->GetParError(0);
    frMu = fn->GetParameter(1);
    frErrMu = fn->GetParError(1);
    frSig = fabs(fn->GetParameter(2));
    frErrSig = fn->GetParError(2);
    frSignif = (eA > 0) ? A / eA : 0;
}

void frOne(TFile* f, const char* eh, const char* nstb, const char* xy, const char* type) {
    TH1* hs = (TH1*)f->Get(Form("%s%s%sSame Event_%s", eh, nstb, xy, type));
    TH1* hm = (TH1*)f->Get(Form("%s%s%sMixed Event_%s", eh, nstb, xy, type));
    if (!hs || !hm) { printf("   missing %s%s%s_%s\n", eh, nstb, xy, type); frSig = 0; return; }
    frRawRms = hs->GetRMS();
    TH1D* sub = frSubtract(hs, hm, 80.0, 200.0);
    frFit(sub);
}

void cmpFcsResid(const char* fOff, const char* fOn, const char* type = "Primary") {
    TFile* a = TFile::Open(fOff);
    TFile* b = TFile::Open(fOn);
    if (!a || a->IsZombie() || !b || b->IsZombie()) { printf("cannot open inputs\n"); return; }
    gStyle->SetOptStat(0);

    printf("\n=== FCS cluster - track residual, %s tracks, Run24 iter-1 geometry ===\n", type);
    printf("    arm A = FTT hits NOT on the track (NOADD=1)\n");
    printf("    arm B = aligned FTT hits ON the track (NOADD=0)\n");
    printf("    mixed-event subtracted, normalised over |d| in 80-200 cm\n\n");
    printf("  %-14s %16s %16s %8s %7s  %9s %9s\n",
           "panel", "A sigma", "B sigma", "B/A", "n_sig", "A ent", "B ent");

    const char* ehs[2] = {"Ecal", "Hcal"};
    const char* nsx[2] = {"North", "South"};
    const char* nsy[2] = {"Top", "Bottom"};
    for (int ie = 0; ie < 2; ie++) {
        for (int iq = 0; iq < 2; iq++) {
            frOne(a, ehs[ie], nsx[iq], "dx", type);
            double sa = frSig, ma = frMu, ea = frEnt, siga = frSignif, esa = frErrSig;
            frOne(b, ehs[ie], nsx[iq], "dx", type);
            double sb = frSig, mb = frMu;
            double eaS = frErrSig;
            if (sa > 0 && sb > 0 && siga > 3) {
                double de = sqrt(esa * esa + eaS * eaS);
                printf("  %-6s %-5s dx %8.3f +-%5.3f %8.3f +-%5.3f %8.3f %7.1f  %9.0f %9.0f\n",
                       ehs[ie], nsx[iq], sa, esa, sb, eaS, sb / sa,
                       de > 0 ? (sb - sa) / de : 0, ea, frEnt);
            }
            else
                printf("  %-6s %-5s dx    weak/absent peak (A %.1f sigma)\n", ehs[ie], nsx[iq], siga);
        }
        for (int jq = 0; jq < 2; jq++) {
            frOne(a, ehs[ie], nsy[jq], "dy", type);
            double sa2 = frSig, ma2 = frMu, ea2 = frEnt, siga2 = frSignif, esa2 = frErrSig;
            frOne(b, ehs[ie], nsy[jq], "dy", type);
            double sb2 = frSig, mb2 = frMu;
            double eaS2 = frErrSig;
            if (sa2 > 0 && sb2 > 0 && siga2 > 3) {
                double de2 = sqrt(esa2 * esa2 + eaS2 * eaS2);
                printf("  %-6s %-5s dy %8.3f +-%5.3f %8.3f +-%5.3f %8.3f %7.1f  %9.0f %9.0f\n",
                       ehs[ie], nsy[jq], sa2, esa2, sb2, eaS2, sb2 / sa2,
                       de2 > 0 ? (sb2 - sa2) / de2 : 0, ea2, frEnt);
            }
            else
                printf("  %-6s %-5s dy    weak/absent peak (A %.1f sigma)\n", ehs[ie], nsy[jq], siga2);
        }
    }
    printf("\n  sigma in cm. B/A < 1 means the aligned FTT hits sharpened the FCS match.\n");
}
