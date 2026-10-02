// cmpTrees.C -- same input file, same arguments, two source trees. Used to check that
// merging Daniel's dev into akio202607 did not change the physics.
//
// The expectation is AGREEMENT, not drift: dev's tracker fixes originated in fcstrk11,
// went out through fcstrk12 and came back when Daniel merged them, so akio202607 already
// contained them before the merge. A difference here is a red flag, not something to
// explain away.
//
// CINT: unique loop variable names, no std::map (see CLAUDE.md).

void ctLine(const char* tag, double a, double b, const char* unit) {
    double d = b - a;
    const char* flag = (fabs(d) < 1e-9) ? "  identical"
                     : (fabs(a) > 0 && fabs(d / a) < 0.01) ? "  <1%"
                     : "   <-- DIFFERS";
    printf("  %-34s %12.5f %12.5f %12.5f %-6s%s\n", tag, a, b, d, unit, flag);
}

void cmpTrees(const char* res11, const char* res15,
              const char* bd11,  const char* bd15) {
    TFile* r1 = TFile::Open(res11); TFile* r2 = TFile::Open(res15);
    if (!r1 || r1->IsZombie() || !r2 || r2->IsZombie()) { printf("cannot open residual files\n"); return; }

    printf("\n=== fcstrk11 vs fcstrk15 (merged), same file and arguments ===\n");
    printf("  %-34s %12s %12s %12s\n", "quantity", "fcstrk11", "fcstrk15", "diff");

    const char* tn[6] = {"Global","BLC","Primary","FwdVtx","BLCVtx","FCSConstrained"};
    for (int it = 0; it < 6; it++) {
        TH1F* h1 = (TH1F*)r1->Get(Form("PlaneUsage/hPlaneUsage_%s", tn[it]));
        TH1F* h2 = (TH1F*)r2->Get(Form("PlaneUsage/hPlaneUsage_%s", tn[it]));
        if (!h1 || !h2) continue;
        ctLine(Form("tracks %s", tn[it]), h1->GetBinContent(20), h2->GetBinContent(20), "");
    }
    for (int id = 0; id < 3; id++) {
        TH1F* h1 = (TH1F*)r1->Get(Form("FST/disk%d/h_fst_d%d_rdphi", id, id));
        TH1F* h2 = (TH1F*)r2->Get(Form("FST/disk%d/h_fst_d%d_rdphi", id, id));
        if (!h1 || !h2) continue;
        ctLine(Form("FST d%d rdphi mean", id), h1->GetMean(), h2->GetMean(), "cm");
        ctLine(Form("FST d%d rdphi rms",  id), h1->GetRMS(),  h2->GetRMS(),  "cm");
    }

    TFile* b1 = TFile::Open(bd11); TFile* b2 = TFile::Open(bd15);
    if (b1 && !b1->IsZombie() && b2 && !b2->IsZombie()) {
        printf("\n  --- FTT blind residual, per disk ---\n");
        for (int jd = 0; jd < 4; jd++) {
            TH1F* x1 = (TH1F*)b1->Get(Form("hBlindDxAll_V_disk%d", jd));
            TH1F* x2 = (TH1F*)b2->Get(Form("hBlindDxAll_V_disk%d", jd));
            TH1F* y1 = (TH1F*)b1->Get(Form("hBlindDyAll_H_disk%d", jd));
            TH1F* y2 = (TH1F*)b2->Get(Form("hBlindDyAll_H_disk%d", jd));
            if (x1 && x2) {
                ctLine(Form("disk%d dx entries", jd), x1->GetEntries(), x2->GetEntries(), "");
                ctLine(Form("disk%d dx mean",    jd), x1->GetMean(),    x2->GetMean(),    "cm");
            }
            if (y1 && y2) ctLine(Form("disk%d dy mean", jd), y1->GetMean(), y2->GetMean(), "cm");
        }
    }
    printf("\n  Any line marked DIFFERS needs a named dev commit to explain it.\n");
}
