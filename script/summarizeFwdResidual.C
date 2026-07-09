// summarizeFwdResidual.C
// Scans a set of StFwdResidualMaker output files (one per dataset x track-type)
// and writes an HTML table flagging planes with a statistically significant
// residual offset or with usage noticeably below its peers.
//
// Usage: root4star -b -q 'summarizeFwdResidual.C("residual/summary.html")'

bool significantOffset(double mean, double rms, int n, double floorCm) {
    if (n < 20) return false;
    if (fabs(mean) < floorCm) return false;
    if (rms <= 0) return false;
    double se = rms / sqrt((double)n);
    return fabs(mean) / se > 3.0;
}

// Builds a TGraphErrors from a TProfile-X of h, keeping only x-bins with
// >=minBinN raw entries AND a strictly positive error. A profile bin with
// exactly 1 entry has undefined spread and ROOT reports its error as 0; less
// obviously, a bin with >=2 entries can *also* come out with error exactly 0
// if those entries happen to be numerically identical (seen in practice with
// n=2-3 low-stats bins). Either way a zero-error point gets infinite weight
// in a chi2 fit and can return a wildly wrong, artificially tiny parameter
// error even with plenty of other points -- excluding any such bin avoids
// that regardless of which of the two ways it arose.
TGraphErrors* profileToGraph(TH2F* h, const char* gname, int minBinN = 2) {
    if (!h) return 0;
    TProfile* p = h->ProfileX(Form("%s_pfxtmp", gname));
    p->SetDirectory(0);
    int nb = p->GetNbinsX();
    double vx[200], vy[200], vey[200];
    int n = 0;
    for (int i = 1; i <= nb; i++) {
        if (p->GetBinEntries(i) < minBinN) continue;
        if (p->GetBinError(i) <= 0) continue;
        vx[n] = p->GetBinCenter(i);
        vy[n] = p->GetBinContent(i);
        vey[n] = p->GetBinError(i);
        n++;
    }
    delete p;
    if (n < 5) return 0;
    return new TGraphErrors(n, vx, vy, 0, vey);
}

// Linear-fits residual vs position (h's x-axis) and returns the slope/error.
// Returns false if there isn't enough data for a stable fit. If the line is a
// poor description of the points (chi2/ndf > 1), the reported error is
// inflated by sqrt(chi2/ndf) (the standard PDG scale-factor approach) so a
// bad fit doesn't masquerade as a tight, highly-significant slope -- profile
// bin errors here come from only a handful of raw hits each and routinely
// underestimate the true point-to-point scatter.
bool fitSlope(TH2F* h, const char* gname, double& slope, double& slopeErr) {
    if (!h || h->GetEntries() < 50) return false;
    TGraphErrors* g = profileToGraph(h, gname);
    if (!g) return false;
    TFitResultPtr r = g->Fit("pol1", "Q0S");
    bool ok = r.Get() && r->IsValid() && r->Ndf() > 0;
    if (ok) {
        slope = r->Parameter(1);
        double chi2ndf = r->Chi2() / r->Ndf();
        double scale = (chi2ndf > 1.0) ? sqrt(chi2ndf) : 1.0;
        slopeErr = r->ParError(1) * scale;
    }
    delete g;
    return ok;
}

bool significantSlope(double slope, double slopeErr, double floor) {
    if (slopeErr <= 0) return false;
    if (fabs(slope) < floor) return false;
    return fabs(slope) / slopeErr > 3.0;
}

// Translation test for FST: a rigid (x,y) shift of the disk projects onto the
// *azimuthal* residual r*dphi as A*sin(phi-phi0) (not a straight line in x,
// y, r, or r*phi -- those axes mix in the radial component, which is why the
// earlier vs-rphi/vs-x linear fits weren't the right probe). Using the
// identity A*sin(phi-phi0) = a*cos(phi) + b*sin(phi), amplitude/phase are
// linear in (a,b), so this is a robust two-parameter linear fit -- no
// nonlinear-phase convergence risk. amp = sqrt(a^2+b^2) = |translation|.
bool fitSine(TH2F* h, const char* gname, double& amp, double& ampErr) {
    if (!h || h->GetEntries() < 50) return false;
    TGraphErrors* g = profileToGraph(h, gname);
    if (!g) return false;
    TF1 fsin(Form("%s_fsin", gname), "[0]*cos(x)+[1]*sin(x)+[2]", -TMath::Pi(), TMath::Pi());
    TFitResultPtr r = g->Fit(&fsin, "Q0S");
    bool ok = r.Get() && r->IsValid() && r->Ndf() > 0;
    if (ok) {
        double a = r->Parameter(0), b = r->Parameter(1);
        double ea = r->ParError(0), eb = r->ParError(1);
        // Seen in practice (even with ~45 well-spread points and a sane
        // chi2/ndf): TF1::Fit() on this cos/sin/const model occasionally
        // returns a ParError() around 1e-9 for both parameters -- not exactly
        // zero, but numerically meaningless (three orders of magnitude below
        // any legitimate error in this table, ~1e-3 to 1e-2). ROOT auto-
        // detects the formula as linear-in-parameters and routes it through a
        // fast linear solver whose error propagation doesn't always come back
        // properly populated. Floor at 1e-6 (well above the bogus ~1e-9
        // values, well below real ones): report no result rather than let
        // amp/epsilon read as infinitely significant -- the same failure mode
        // the chi2/ndf safeguards below exist to catch, just via a different
        // mechanism that a plain "<=0" check doesn't.
        if (ea < 1e-6 || eb < 1e-6) { delete g; return false; }
        double chi2ndf = r->Chi2() / r->Ndf();
        double scale = (chi2ndf > 1.0) ? sqrt(chi2ndf) : 1.0;
        double A = sqrt(a*a + b*b);
        amp = A;
        ampErr = (A > 0) ? sqrt(a*a*ea*ea + b*b*eb*eb) / A * scale : 0;
    }
    delete g;
    return ok;
}

void writeAmpCell(FILE* fp, bool have, double amp, double ampErr, bool flagged) {
    const char* bg = flagged ? "#ffd6d6" : (have ? "#eaffea" : "#f0f0f0");
    fprintf(fp, "<td style=\"background:%s\">", bg);
    if (have) fprintf(fp, "A=%.4f%s<br>&plusmn;%.4f", amp, flagged ? " &#9888;" : "", ampErr);
    else fprintf(fp, "n/a");
    fprintf(fp, "</td>\n");
}

void writeSlopeCell(FILE* fp, bool have, double slope, double slopeErr, bool flagged) {
    const char* bg = flagged ? "#ffd6d6" : (have ? "#eaffea" : "#f0f0f0");
    fprintf(fp, "<td style=\"background:%s\">", bg);
    if (have) fprintf(fp, "slope=%.4f%s<br>&plusmn;%.4f", slope, flagged ? " &#9888;" : "", slopeErr);
    else fprintf(fp, "n/a");
    fprintf(fp, "</td>\n");
}

void writeCell(FILE* fp, double usage, bool flagUsage, double mean, double rms, int n, bool haveRes, bool flagOffset) {
    bool flagged = flagOffset || flagUsage;
    const char* bg = flagged ? "#ffd6d6" : "#eaffea";
    fprintf(fp, "<td style=\"background:%s\">", bg);
    if (usage >= 0) fprintf(fp, "usage=%.0f%%%s<br>", usage*100, flagUsage ? " &#9888;" : "");
    else fprintf(fp, "usage=n/a<br>");
    if (haveRes) {
        fprintf(fp, "&mu;=%.3f cm%s<br>", mean, flagOffset ? " &#9888;" : "");
        fprintf(fp, "&sigma;=%.3f cm", rms);
    } else {
        fprintf(fp, "&mu;=n/a<br>&sigma;=n/a");
    }
    fprintf(fp, "</td>\n");
}

void summarizeFwdResidual(const char* outHtml = "residual/summary.html") {

    const int NSET = 6;
    const char* dsFile[NSET] = {
        "st_fwd_23081015_raw_6000057.FwdDetResidual_Global.root",
        "st_fwd_23081015_raw_6000057.FwdDetResidual_BLC.root",
        "st_fwd_23081015_raw_6000057.FwdDetResidual_BLCVtx.root",
        "pythia.JPsi.vz0.run1.FwdDetResidual_Global.root",
        "pythia.JPsi.vz0.run1.FwdDetResidual_BLC.root",
        "pythia.JPsi.vz0.run1.FwdDetResidual_BLCVtx.root"
    };
    const char* dsLabel[NSET] = {
        "data / Global", "data / BLC", "data / BLCVtx",
        "pythia / Global", "pythia / BLC", "pythia / BLCVtx"
    };
    const char* dsType[NSET] = {"Global","BLC","BLCVtx","Global","BLC","BLCVtx"};

    FILE* fp = fopen(outHtml, "w");
    fprintf(fp, "\n<HTML>\n<HEAD>\n<TITLE> OGAWA Akio Home page </TITLE>\n</HEAD>\n");
    fprintf(fp, "<BODY BGCOLOR=\"#dfdfff\" TEXT=\"black\" LINK=\"blue\" VLINK=\"darkblue\">\n");
    fprintf(fp, "<p><a href=\"index.html\">&larr; back to residual index</a></p>\n");
    fprintf(fp, "<H1>FST/FTT Residual &amp; Plane-Usage Outlier Summary</H1><HR>\n");
    fprintf(fp, "<p>Auto-generated by <code>script/summarizeFwdResidual.C</code>. Two tables below: plane "
                "usage/offset, and residual-vs-position slope. A cell is flagged (red, &#9888;) if:</p>\n<ul>\n");
    fprintf(fp, "<li>usage &mdash; FST usage &lt;90%%, or an FTT plane's x/y usage differs from the mean of the other 3 FTT planes by &gt;15 percentage points</li>\n");
    fprintf(fp, "<li>residual offset &mdash; |mean| is &gt;3&sigma;/&radic;N from zero AND exceeds a physical floor (0.01 cm FST, 0.05 cm FTT), so large-N statistical flukes don't flag</li>\n");
    fprintf(fp, "<li>residual slope &mdash; see the slope table's own intro below for its criteria</li>\n</ul>\n");
    fprintf(fp, "<p><b>Why pythia's FST usage runs 75&ndash;84%% while data is always 100%%:</b> "
                "not a StFwdResidualMaker or tracking bug. Traced via targeted instrumentation to "
                "<code>ObjExporter::trackPosition()</code> (StFwdTrackMaker/include/Tracker/ObjExporter.h) "
                "returning its explicit failure sentinel <code>(-990,-990,-990)</code> whenever GenFit's "
                "<code>extrapolateToPlane()</code> throws while projecting to an individual FST z-plane -- "
                "the FST hit is still there (raw seed-point counts checked directly: always 3/3), it's the "
                "projection used for z-matching that's missing, so <code>findClosestZ()</code> can't pair "
                "the two and the plane-usage bin never gets filled. <code>script/sim.C</code> (pythia) uses "
                "much looser fit convergence than the real-data afterburner "
                "(<code>setFitMinIterations(10)/setFitMaxIterations(20)</code> vs. "
                "<code>40</code>/<code>100</code> in <code>fwd_afterburner_db.C</code>), so more MC tracks "
                "end up in fit states where this extrapolation throws. A config difference between the two "
                "macros, not a defect in either.</p>\n"
                "<p><b>Why data/BLCVtx-FST1 is flagged:</b> purely the &gt;3&sigma;/&radic;N significance "
                "test -- with &sigma;=0.188&nbsp;cm and N=4250, the standard error is only ~0.0029&nbsp;cm, "
                "so the measured &mu;=-0.011&nbsp;cm (about 110&nbsp;&mu;m) clears the 3&sigma; bar despite "
                "being tiny in absolute terms. Statistically real, physically negligible.</p>\n");

    fprintf(fp, "<table border=1 cellpadding=4>\n<tr><th>Set</th><th>FST1</th><th>FST2</th><th>FST3</th>"
                "<th>FTT1-x</th><th>FTT1-y</th><th>FTT2-x</th><th>FTT2-y</th>"
                "<th>FTT3-x</th><th>FTT3-y</th><th>FTT4-x</th><th>FTT4-y</th></tr>\n");

    int nFlagged = 0, nRows = 0;
    for (int is = 0; is < NSET; is++) {
        TFile* f = TFile::Open(dsFile[is]);
        if (!f || f->IsZombie()) {
            fprintf(fp, "<tr><td>%s</td><td colspan=10>missing: %s</td></tr>\n", dsLabel[is], dsFile[is]);
            printf("WARN: cannot open %s\n", dsFile[is]);
            continue;
        }

        TH1F* hu = (TH1F*)f->Get(Form("PlaneUsage/hPlaneUsage_%s", dsType[is]));
        double all = hu ? hu->GetBinContent(20) : 0;

        // FST: usage[3], mean/rms/n[3]
        double fstUsage[3], fstMean[3], fstRms[3]; int fstN[3]; bool fstHave[3], fstFlagOff[3], fstFlagUse[3];
        for (int d = 0; d < 3; d++) {
            fstUsage[d] = (hu && all > 0) ? hu->GetBinContent(2+d) / all : -1;
            fstHave[d] = false; fstMean[d] = 0; fstRms[d] = 0; fstN[d] = 0; fstFlagOff[d] = false;
            TH1F* h = (TH1F*)f->Get(Form("FST/disk%d/h_fst_d%d_rdphi", d, d));
            if (h && h->GetEntries() > 0) {
                fstHave[d] = true; fstN[d] = (int)h->GetEntries();
                fstMean[d] = h->GetMean(); fstRms[d] = h->GetRMS();
                fstFlagOff[d] = significantOffset(fstMean[d], fstRms[d], fstN[d], 0.01);
            }
            fstFlagUse[d] = (fstUsage[d] >= 0 && fstUsage[d] < 0.90);
        }

        // FTT: x and y usage/mean/rms/n[4]
        double fttxUsage[4], fttxMean[4], fttxRms[4]; int fttxN[4]; bool fttxHave[4], fttxFlagOff[4], fttxFlagUse[4];
        double fttyUsage[4], fttyMean[4], fttyRms[4]; int fttyN[4]; bool fttyHave[4], fttyFlagOff[4], fttyFlagUse[4];
        double sumX = 0, sumY = 0; int cntX = 0, cntY = 0;
        for (int p = 0; p < 4; p++) {
            fttxUsage[p] = (hu && all > 0) ? hu->GetBinContent(4 + 3*p + 1) / all : -1;
            fttyUsage[p] = (hu && all > 0) ? hu->GetBinContent(4 + 3*p + 2) / all : -1;
            fttxHave[p] = false; fttxMean[p] = 0; fttxRms[p] = 0; fttxN[p] = 0; fttxFlagOff[p] = false;
            fttyHave[p] = false; fttyMean[p] = 0; fttyRms[p] = 0; fttyN[p] = 0; fttyFlagOff[p] = false;

            TH1F* hx = (TH1F*)f->Get(Form("FTT/plane%d/h_ftt_p%d_dx", p, p));
            if (hx && hx->GetEntries() > 0) {
                fttxHave[p] = true; fttxN[p] = (int)hx->GetEntries();
                fttxMean[p] = hx->GetMean(); fttxRms[p] = hx->GetRMS();
                fttxFlagOff[p] = significantOffset(fttxMean[p], fttxRms[p], fttxN[p], 0.05);
            }
            TH1F* hy = (TH1F*)f->Get(Form("FTT/plane%d/h_ftt_p%d_dy", p, p));
            if (hy && hy->GetEntries() > 0) {
                fttyHave[p] = true; fttyN[p] = (int)hy->GetEntries();
                fttyMean[p] = hy->GetMean(); fttyRms[p] = hy->GetRMS();
                fttyFlagOff[p] = significantOffset(fttyMean[p], fttyRms[p], fttyN[p], 0.05);
            }
            if (fttxUsage[p] >= 0) { sumX += fttxUsage[p]; cntX++; }
            if (fttyUsage[p] >= 0) { sumY += fttyUsage[p]; cntY++; }
        }
        double avgX = (cntX > 0) ? sumX/cntX : 0;
        double avgY = (cntY > 0) ? sumY/cntY : 0;
        for (int p = 0; p < 4; p++) {
            fttxFlagUse[p] = (fttxUsage[p] >= 0 && fabs(fttxUsage[p]-avgX) > 0.15);
            fttyFlagUse[p] = (fttyUsage[p] >= 0 && fabs(fttyUsage[p]-avgY) > 0.15);
        }

        fprintf(fp, "<tr><td><b>%s</b></td>\n", dsLabel[is]);
        for (int d = 0; d < 3; d++) {
            writeCell(fp, fstUsage[d], fstFlagUse[d], fstMean[d], fstRms[d], fstN[d], fstHave[d], fstFlagOff[d]);
            if (fstFlagUse[d] || fstFlagOff[d]) nFlagged++;
        }
        for (int p = 0; p < 4; p++) {
            writeCell(fp, fttxUsage[p], fttxFlagUse[p], fttxMean[p], fttxRms[p], fttxN[p], fttxHave[p], fttxFlagOff[p]);
            if (fttxFlagUse[p] || fttxFlagOff[p]) nFlagged++;
            writeCell(fp, fttyUsage[p], fttyFlagUse[p], fttyMean[p], fttyRms[p], fttyN[p], fttyHave[p], fttyFlagOff[p]);
            if (fttyFlagUse[p] || fttyFlagOff[p]) nFlagged++;
        }
        fprintf(fp, "</tr>\n");
        nRows++;
        f->Close();
    }
    fprintf(fp, "</table>\n");
    fprintf(fp, "<p>%d flagged cell(s) across %d row(s).</p>\n", nFlagged, nRows);

    fprintf(fp, "<H3>Misalignment-Mode Tests</H3>\n");
    fprintf(fp, "<p>The three tables below replace an earlier single \"self-axis slope\" table that fit "
                "everything with one straight line -- misleading, since a rigid detector offset doesn't "
                "actually show up as a simple slope for every axis/geometry combination. Each table below "
                "targets the actual functional form a specific misalignment mode produces, derived from "
                "first-order geometry:</p>\n<ul>\n"
                "<li><b>Translation</b> (rigid x,y shift) -- for FTT (Cartesian), this is just a constant "
                "&Delta;x/&Delta;y offset, already covered by the mean/offset table above. For FST (polar), "
                "the same shift decomposes differently in the radial vs. azimuthal directions as &phi; "
                "changes, so it appears in r&middot;&Delta;&phi; as <b>A&middot;sin(&phi;-&phi;<sub>0</sub>)</b>, "
                "not a straight line in x, r, or r&middot;&phi;.</li>\n"
                "<li><b>Rotation</b> about the beam axis by &delta;&theta; -- FST: r&middot;&Delta;&phi; = "
                "r&middot;&delta;&theta;, linear in r, flat in &phi;. FTT: &Delta;x=-y&middot;&delta;&theta;, "
                "&Delta;y=+x&middot;&delta;&theta; -- linear, but in the <i>other</i> (cross) coordinate.</li>\n"
                "<li><b>Length-scale / gain error</b> (e.g. wrong strip pitch calibration) -- FTT only: "
                "&Delta;x&prop;x, &Delta;y&prop;y, linear in the <i>same</i> (self) coordinate. This is what "
                "the old single table was actually testing for FTT.</li>\n</ul>\n"
                "<p>All fits use the same artifact safeguards as before: profile x-bins need &ge;2 raw "
                "entries, and the reported error is inflated by &radic;(&chi;&sup2;/ndf) when &gt;1 (PDG "
                "scale-factor method) so a poor-fitting line doesn't masquerade as high significance.</p>\n");

    fprintf(fp, "<H4>Translation test (FST only -- sine amplitude vs &phi;)</H4>\n");
    fprintf(fp, "<table border=1 cellpadding=4>\n<tr><th>Set</th><th>FST1</th><th>FST2</th><th>FST3</th></tr>\n");
    int nAmpFlagged = 0, nAmpRows = 0;
    for (int is = 0; is < NSET; is++) {
        TFile* f = TFile::Open(dsFile[is]);
        if (!f || f->IsZombie()) {
            fprintf(fp, "<tr><td>%s</td><td colspan=3>missing: %s</td></tr>\n", dsLabel[is], dsFile[is]);
            continue;
        }
        fprintf(fp, "<tr><td><b>%s</b></td>\n", dsLabel[is]);
        for (int d = 0; d < 3; d++) {
            TH2F* h2 = (TH2F*)f->Get(Form("FST/disk%d/h2_fst_d%d_rdphi_vs_phi", d, d));
            double amp = 0, err = 0;
            bool have = fitSine(h2, Form("fstsin_%d_%d", is, d), amp, err);
            bool flagged = have && err > 0 && amp >= 0.01 && amp/err > 3.0;
            writeAmpCell(fp, have, amp, err, flagged);
            if (flagged) nAmpFlagged++;
        }
        fprintf(fp, "</tr>\n");
        nAmpRows++;
        f->Close();
    }
    fprintf(fp, "</table>\n<p>%d flagged cell(s) across %d row(s). (Needs the new h2_fst_..._vs_phi "
                "histogram -- n/a for files produced before this was added.)</p>\n", nAmpFlagged, nAmpRows);

    fprintf(fp, "<H4>Rotation test (FST vs r; FTT cross-axis)</H4>\n");
    fprintf(fp, "<table border=1 cellpadding=4>\n<tr><th>Set</th><th>FST1 (vs r)</th><th>FST2 (vs r)</th><th>FST3 (vs r)</th>"
                "<th>FTT1 &Delta;y(x)</th><th>FTT1 &Delta;x(y)</th><th>FTT2 &Delta;y(x)</th><th>FTT2 &Delta;x(y)</th>"
                "<th>FTT3 &Delta;y(x)</th><th>FTT3 &Delta;x(y)</th><th>FTT4 &Delta;y(x)</th><th>FTT4 &Delta;x(y)</th></tr>\n");
    int nRotFlagged = 0, nRotRows = 0;
    for (int is = 0; is < NSET; is++) {
        TFile* f = TFile::Open(dsFile[is]);
        if (!f || f->IsZombie()) {
            fprintf(fp, "<tr><td>%s</td><td colspan=10>missing: %s</td></tr>\n", dsLabel[is], dsFile[is]);
            continue;
        }
        fprintf(fp, "<tr><td><b>%s</b></td>\n", dsLabel[is]);
        for (int d = 0; d < 3; d++) {
            TH2F* h2 = (TH2F*)f->Get(Form("FST/disk%d/h2_fst_d%d_rdphi_vs_r", d, d));
            double slope = 0, err = 0;
            bool have = fitSlope(h2, Form("fstrot_%d_%d", is, d), slope, err);
            bool flagged = have && significantSlope(slope, err, 0.0005);
            writeSlopeCell(fp, have, slope, err, flagged);
            if (flagged) nRotFlagged++;
        }
        for (int p = 0; p < 4; p++) {
            TH2F* hyx = (TH2F*)f->Get(Form("FTT/plane%d/h2_ftt_p%d_dy_vs_x", p, p));
            double slopeYX = 0, errYX = 0;
            bool haveYX = fitSlope(hyx, Form("rotyx_%d_%d", is, p), slopeYX, errYX);
            bool flagYX = haveYX && significantSlope(slopeYX, errYX, 0.001);
            writeSlopeCell(fp, haveYX, slopeYX, errYX, flagYX);
            if (flagYX) nRotFlagged++;

            TH2F* hxy = (TH2F*)f->Get(Form("FTT/plane%d/h2_ftt_p%d_dx_vs_y", p, p));
            double slopeXY = 0, errXY = 0;
            bool haveXY = fitSlope(hxy, Form("rotxy_%d_%d", is, p), slopeXY, errXY);
            bool flagXY = haveXY && significantSlope(slopeXY, errXY, 0.001);
            writeSlopeCell(fp, haveXY, slopeXY, errXY, flagXY);
            if (flagXY) nRotFlagged++;
        }
        fprintf(fp, "</tr>\n");
        nRotRows++;
        f->Close();
    }
    fprintf(fp, "</table>\n<p>%d flagged cell(s) across %d row(s).</p>\n", nRotFlagged, nRotRows);

    fprintf(fp, "<H4>Length-scale / gain test (FTT self-axis; FST has no clean analog here)</H4>\n");
    fprintf(fp, "<table border=1 cellpadding=4>\n<tr><th>Set</th>"
                "<th>FTT1-x (vs x)</th><th>FTT1-y (vs y)</th><th>FTT2-x (vs x)</th><th>FTT2-y (vs y)</th>"
                "<th>FTT3-x (vs x)</th><th>FTT3-y (vs y)</th><th>FTT4-x (vs x)</th><th>FTT4-y (vs y)</th></tr>\n");
    int nSlopeFlagged = 0, nSlopeRows = 0;
    for (int is = 0; is < NSET; is++) {
        TFile* f = TFile::Open(dsFile[is]);
        if (!f || f->IsZombie()) {
            fprintf(fp, "<tr><td>%s</td><td colspan=8>missing: %s</td></tr>\n", dsLabel[is], dsFile[is]);
            continue;
        }

        fprintf(fp, "<tr><td><b>%s</b></td>\n", dsLabel[is]);
        for (int p = 0; p < 4; p++) {
            TH2F* h2ttx = (TH2F*)f->Get(Form("FTT/plane%d/h2_ftt_p%d_dx_vs_x", p, p));
            double slopeX = 0, errX = 0;
            bool haveX = fitSlope(h2ttx, Form("fttx_%d_%d", is, p), slopeX, errX);
            bool flagX = haveX && significantSlope(slopeX, errX, 0.001);
            writeSlopeCell(fp, haveX, slopeX, errX, flagX);
            if (flagX) nSlopeFlagged++;

            TH2F* h2tty = (TH2F*)f->Get(Form("FTT/plane%d/h2_ftt_p%d_dy_vs_y", p, p));
            double slopeY = 0, errY = 0;
            bool haveY = fitSlope(h2tty, Form("ftty_%d_%d", is, p), slopeY, errY);
            bool flagY = haveY && significantSlope(slopeY, errY, 0.001);
            writeSlopeCell(fp, haveY, slopeY, errY, flagY);
            if (flagY) nSlopeFlagged++;
        }
        fprintf(fp, "</tr>\n");
        nSlopeRows++;
        f->Close();
    }
    fprintf(fp, "</table>\n");
    fprintf(fp, "<p>%d flagged cell(s) across %d row(s). Pattern seen previously: a small "
                "(0.001&ndash;0.006 cm/cm) but significant self-axis slope in nearly every FTT plane, both "
                "data and pythia, all three track types -- consistent enough across independent samples to "
                "look like a real length-scale/gain effect rather than noise, worth a dedicated look once "
                "statistics improve.</p>\n", nSlopeFlagged, nSlopeRows);
    fprintf(fp, "<p><b>Checked and rejected:</b> an earlier ad-hoc scan flagged "
                "<code>pythia/Global FTT plane3 &Delta;y vs r&middot;&phi;</code> as a huge outlier "
                "(slope&asymp;-0.015, &chi;&sup2;/ndf&asymp;20) -- not a real effect. Pythia is an idealized "
                "MC with no real misalignment applied, so a true geometry-driven slope there would mean a "
                "reconstruction bug, but &chi;&sup2;/ndf&asymp;20 (vs. &asymp;1.2&ndash;8 for the fits in the "
                "tables above) means the line was a bad fit to begin with -- a low-statistics artifact "
                "(N=276 across 38 bins, several with just 1 raw hit and wild swings, e.g. one bin at "
                "+1.9&nbsp;cm), not a bug.</p>\n");

    fprintf(fp, "<HR>OGAWA, Akio<ADDRESS><A HREF=\"mailto://akio@bnl.gov\">akio@bnl.gov</A></ADDRESS>\n</BODY></HTML>\n");
    fclose(fp);
    printf("Wrote %s (offset: %d/%d rows; translation: %d/%d rows; rotation: %d/%d rows; length-scale: %d/%d rows)\n",
           outHtml, nFlagged, nRows, nAmpFlagged, nAmpRows, nRotFlagged, nRotRows, nSlopeFlagged, nSlopeRows);
}
