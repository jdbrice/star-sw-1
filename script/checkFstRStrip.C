// checkFstRStrip.C
//
// Follow-up to script/checkFstInnerOuter.C (see fst_wedge_testA_results.txt).
//
// That macro showed the 12-fold wedge structure in mean r*dphi lives almost
// entirely in the INNER FST sensors (r < 16.5 cm), by a factor 6-18 over the
// outer pair. That already rules out both rigid-body explanations:
//   - per-wedge ROTATION  -> d(phi) constant -> r*dphi grows AS r
//   - per-wedge TRANSLATION -> r*dphi constant, independent of r
// Neither gives inner >> outer.
//
// This macro splits finer: by the 8 FST r-strips, whose radii are fixed by
// kFstrStart + half pitch (kFstStripPitchR = 2.875 cm):
//   inner rStrip 0-3 -> r = 6.4375, 9.3125, 12.1875, 15.0625 cm
//   outer rStrip 4-7 -> r = 17.9375, 20.8125, 23.6875, 26.5625 cm
// and reports, per (disk, rStrip), how big the wedge-to-wedge structure is
// in BOTH r*dphi and dphi.
//
// Reading the result:
//   dphi structure flat vs r      -> rotation-like (r*dphi would grow as r)
//   r*dphi structure flat vs r    -> translation-like
//   both falling steeply with r   -> neither; something specific to the
//                                    inner sensor, which is what the
//                                    inner/outer split already suggested
//
// The per-sector pattern is irregular rather than a clean two-level
// alternation, so the size of the structure is quantified as the RMS of the
// 12 sector means about their own average (a global per-(disk,rStrip) offset
// is removed first -- that is a rotation of the whole disk and is absorbed by
// tracking anyway). The statistical noise floor on that RMS is reported next
// to it, so a small RMS can be told apart from a measurement limit.
//
// Usage:
//   root4star -b -q 'script/checkFstRStrip.C'
//   root4star -b -q 'script/checkFstRStrip.C("<glob>","<outdir>")'

void checkFstRStrip(
        const char* fname  = "/gpfs01/star/pwg_tasks/FwdCalib/akio/pico_20260704/alignment/*.FwdAlignment_BLCVtx.root",
        const char* outdir = "residual/data/alignment/rstrip" ){

    TChain* tr = new TChain("alignTree");
    int nAdded = tr->Add(fname);
    if (nAdded <= 0){ printf("No files matched %s\n", fname); return; }
    printf("Matched %d files\n", nAdded);
    gSystem->mkdir(outdir, kTRUE);
    gStyle->SetOptStat(0);

    Int_t   bDet = 0;
    Float_t bHx = 0, bHy = 0, bHz = 0, bPx = 0, bPy = 0;
    tr->SetBranchAddress("detType", &bDet);
    tr->SetBranchAddress("hitX",  &bHx);
    tr->SetBranchAddress("hitY",  &bHy);
    tr->SetBranchAddress("hitZ",  &bHz);
    tr->SetBranchAddress("projX", &bPx);
    tr->SetBranchAddress("projY", &bPy);

    double diskZ[3]; diskZ[0] = 151.750; diskZ[1] = 165.248; diskZ[2] = 178.781;

    // kFstrStart[i] + 0.5*kFstStripPitchR
    double stripR[8];
    stripR[0] =  6.4375; stripR[1] =  9.3125; stripR[2] = 12.1875; stripR[3] = 15.0625;
    stripR[4] = 17.9375; stripR[5] = 20.8125; stripR[6] = 23.6875; stripR[7] = 26.5625;

    // IMPORTANT: r*dphi must be truncated before averaging. At small r, phi is
    // ill-defined -- a projection off by ~1 cm at r=6.4 cm is already 160 mrad,
    // and a projection landing across the beam axis gives dphi up to pi. The
    // raw mean is then dominated by tails and is meaningless (untruncated, the
    // rStrip-0 "structure" comes out at several cm, i.e. hundreds of mrad).
    // checkFstInnerOuter.C truncated implicitly, via its TProfile y-range of
    // +-3 cm (TProfile silently drops out-of-range fills). Do the same here
    // explicitly, and carry a tighter +-0.5 cm window alongside so the result
    // can be checked for sensitivity to the choice.
    // A fixed r*dphi cut is a DIFFERENT dphi window at each radius (3 cm is
    // 470 mrad at r=6.4 cm but 170 mrad at r=17.9 cm), which biases exactly
    // the r-dependence being measured here. Window C is a constant-ANGLE cut
    // and is therefore the control: it treats every r-strip identically, and
    // it is the one to trust for the vs-r trend.
    const double kCutA = 3.0;    // matches checkFstInnerOuter.C
    const double kCutB = 0.5;    // matches the axis range of the published plots
    const double kCutC = 0.030;  // rad -- constant-angle control window

    // accumulators, flat: idx = (disk*8 + rStrip)*12 + sector  -> 288 cells
    double aN[288], aS[288], aQ[288];
    double bN[288], bS[288], bQ[288];
    double cN[288], cS[288];
    for (int ia = 0; ia < 288; ia++){
        aN[ia] = 0; aS[ia] = 0; aQ[ia] = 0;
        bN[ia] = 0; bS[ia] = 0; bQ[ia] = 0;
        cN[ia] = 0; cS[ia] = 0;
    }
    Long64_t nOff = 0;

    Long64_t nEnt = tr->GetEntries();
    printf("Total rows: %lld\n", nEnt);

    for (Long64_t ie = 0; ie < nEnt; ie++){
        tr->GetEntry(ie);
        if (bDet != 0) continue;

        int dsel = -1;
        for (int db = 0; db < 3; db++)
            if (fabs(bHz - diskZ[db]) < 5.0) dsel = db;
        if (dsel < 0) continue;

        double rr = sqrt(bHx*bHx + bHy*bHy);

        // nearest r-strip; radii are quantised so this is unambiguous
        int rs = -1; double rbest = 1.0;   // require within 1 cm of a strip radius
        for (int ir = 0; ir < 8; ir++){
            double dd = fabs(rr - stripR[ir]);
            if (dd < rbest){ rbest = dd; rs = ir; }
        }
        if (rs < 0){ nOff++; continue; }

        double phih = atan2(bHy, bHx);
        double phip = atan2(bPy, bPx);
        double dph  = phih - phip;
        while (dph >  TMath::Pi()) dph -= 2*TMath::Pi();
        while (dph < -TMath::Pi()) dph += 2*TMath::Pi();
        double rdp = rr * dph;

        double p360 = phih * 180.0 / TMath::Pi();
        if (p360 < 0) p360 += 360.0;
        int sec = (int)(p360 / 30.0); if (sec > 11) sec = 11;

        int cell = (dsel*8 + rs)*12 + sec;
        if (fabs(rdp) < kCutA){ aN[cell] += 1; aS[cell] += rdp; aQ[cell] += rdp*rdp; }
        if (fabs(rdp) < kCutB){ bN[cell] += 1; bS[cell] += rdp; bQ[cell] += rdp*rdp; }
        if (fabs(dph) < kCutC){ cN[cell] += 1; cS[cell] += dph*1000.0; }  // mrad

        if ((ie % 20000000) == 0 && ie > 0) printf("  ... %lld / %lld\n", ie, nEnt);
    }
    printf("Rows not matching any strip radius: %lld\n", nOff);

    // ---- report ----
    printf("\n=================================================================================\n");
    printf(" Wedge-to-wedge structure per r-strip (RMS of the 12 sector means, offset removed)\n");
    printf("=================================================================================\n");
    printf(" 'noise' is the statistical floor: sqrt(<err_i^2>) over the 12 sector means.\n");
    printf(" A structure RMS at or below 'noise' is not a measurement.\n\n");
    printf(" |r*dphi| < %.1f cm (matches checkFstInnerOuter.C) and < %.1f cm.\n\n", kCutA, kCutB);
    printf(" disk rStrip   r[cm]   rows/sec |  RMS rdphi  noise  RMSdphi | tight  | CONTROL |dphi|<%.0fmrad\n", kCutC*1000);
    printf("                                |    [cm]      [cm]   [mrad] | [mrad] | RMSdphi[mrad] RMSrdphi[cm]\n");

    // for the summary plot
    TGraph* grRdp[3];
    TGraph* grDp [3];
    for (int dg = 0; dg < 3; dg++){ grRdp[dg] = new TGraph(); grDp[dg] = new TGraph(); }

    for (int dc = 0; dc < 3; dc++){
        for (int rc = 0; rc < 8; rc++){
            double mu[12], er[12];
            double tot = 0, cnt = 0, minN = 1e18;
            for (int sc = 0; sc < 12; sc++){
                int cl = (dc*8 + rc)*12 + sc;
                if (aN[cl] < 50){ mu[sc] = 0; er[sc] = 0; continue; }
                mu[sc] = aS[cl]/aN[cl];
                double var = aQ[cl]/aN[cl] - mu[sc]*mu[sc];
                if (var < 0) var = 0;
                er[sc] = sqrt(var/aN[cl]);
                tot += mu[sc]; cnt += 1;
                if (aN[cl] < minN) minN = aN[cl];
            }
            if (cnt < 6) continue;
            double avg = tot/cnt;

            double ss = 0, ee = 0;
            for (int sd = 0; sd < 12; sd++){
                int cl2 = (dc*8 + rc)*12 + sd;
                if (aN[cl2] < 50) continue;
                ss += (mu[sd]-avg)*(mu[sd]-avg);
                ee += er[sd]*er[sd];
            }
            double rmsRdp = sqrt(ss/cnt);
            double noise  = sqrt(ee/cnt);
            double rmsDp  = rmsRdp / stripR[rc] * 1000.0;   // mrad

            // same, in the tight window
            double tt = 0, tc = 0, tavg = 0;
            for (int sg = 0; sg < 12; sg++){
                int cl4 = (dc*8 + rc)*12 + sg;
                if (bN[cl4] < 50) continue;
                tavg += bS[cl4]/bN[cl4]; tc += 1;
            }
            double rmsB = 0;
            if (tc >= 6){
                tavg /= tc;
                for (int sh = 0; sh < 12; sh++){
                    int cl5 = (dc*8 + rc)*12 + sh;
                    if (bN[cl5] < 50) continue;
                    double mb = bS[cl5]/bN[cl5];
                    tt += (mb-tavg)*(mb-tavg);
                }
                rmsB = sqrt(tt/tc);
            }

            // control window: constant angle, so directly comparable across r
            double uavg = 0, uc = 0, uu = 0;
            for (int sk = 0; sk < 12; sk++){
                int cl6 = (dc*8 + rc)*12 + sk;
                if (cN[cl6] < 50) continue;
                uavg += cS[cl6]/cN[cl6]; uc += 1;
            }
            double rmsC = 0;
            if (uc >= 6){
                uavg /= uc;
                for (int sm = 0; sm < 12; sm++){
                    int cl7 = (dc*8 + rc)*12 + sm;
                    if (cN[cl7] < 50) continue;
                    double mc = cS[cl7]/cN[cl7];
                    uu += (mc-uavg)*(mc-uavg);
                }
                rmsC = sqrt(uu/uc);
            }

            printf("   %d    %d    %7.4f %9.0f |  %8.4f %7.4f %8.3f | %6.3f | %10.4f %12.4f\n",
                   dc, rc, stripR[rc], minN, rmsRdp, noise, rmsDp,
                   rmsB/stripR[rc]*1000.0,
                   rmsC, rmsC/1000.0*stripR[rc]);

            // plot the CONTROL window -- the only r-unbiased one
            grRdp[dc]->SetPoint(grRdp[dc]->GetN(), stripR[rc], rmsC/1000.0*stripR[rc]);
            grDp [dc]->SetPoint(grDp [dc]->GetN(), stripR[rc], rmsC);
        }
        printf("\n");
    }

    // ---- front/back (even/odd sector) alternation, in the CONTROL window ----
    // Even sector index = the wedges GEANT places at one z, odd = the other
    // (see fst_wedge_testA_results.txt section 5). The parity is FLIPPED on
    // disk 2, so a front/back-driven effect must change sign on FST1.
    printf("=================================================================================\n");
    printf(" Front/back alternation, CONTROL window (constant angle), r*dphi [cm]\n");
    printf("=================================================================================\n");
    printf(" amp = 0.5*(<even sectors> - <odd sectors>); resid = RMS about that alternation.\n");
    printf(" A front/back effect gives |amp| >> resid, and amp flipping sign on FST1.\n\n");
    printf(" disk rStrip    r[cm]   amp [cm]   resid [cm]   |amp|/resid\n");
    for (int dp = 0; dp < 3; dp++){
        for (int rp = 0; rp < 8; rp++){
            double ev = 0, od = 0, ne = 0, no = 0;
            for (int sp = 0; sp < 12; sp++){
                int cl8 = (dp*8 + rp)*12 + sp;
                if (cN[cl8] < 50) continue;
                double vc = cS[cl8]/cN[cl8]/1000.0*stripR[rp];  // mrad -> cm
                if (sp % 2 == 0){ ev += vc; ne += 1; } else { od += vc; no += 1; }
            }
            if (ne < 3 || no < 3) continue;
            double amp = 0.5*(ev/ne - od/no);
            double ctr = 0.5*(ev/ne + od/no);
            double rr2 = 0, nr2 = 0;
            for (int sq = 0; sq < 12; sq++){
                int cl9 = (dp*8 + rp)*12 + sq;
                if (cN[cl9] < 50) continue;
                double vd = cS[cl9]/cN[cl9]/1000.0*stripR[rp];
                double pred = ctr + ((sq % 2 == 0) ? amp : -amp);
                rr2 += (vd-pred)*(vd-pred); nr2 += 1;
            }
            double resid = sqrt(rr2/nr2);
            printf("   %d    %d     %7.4f  %+8.4f    %8.4f     %8.2f\n",
                   dp, rp, stripR[rp], amp, resid, resid > 0 ? fabs(amp)/resid : 0);
        }
        printf("\n");
    }

    // ---- per-sector detail for the inner strips, where the effect lives ----
    printf("=================================================================================\n");
    printf(" Per-sector mean r*dphi [cm], inner r-strips only (window A, +-3 cm)\n");
    printf("=================================================================================\n");
    for (int dj = 0; dj < 3; dj++){
        printf("\n FST%d  sector:", dj);
        for (int se = 0; se < 12; se++) printf(" %6d", se*30);
        printf("\n");
        for (int re = 0; re < 4; re++){
            printf("   rStrip %d (r=%5.2f):", re, stripR[re]);
            for (int sf = 0; sf < 12; sf++){
                int cl3 = (dj*8 + re)*12 + sf;
                double vv = aN[cl3] > 50 ? aS[cl3]/aN[cl3] : 0;
                printf(" %+6.3f", vv);
            }
            printf("\n");
        }
    }

    // ---- plot ----
    TCanvas* cv = new TCanvas("c_rstrip", "", 900, 400);
    cv->Divide(2,1);
    int col[3]; col[0] = kBlue+1; col[1] = kRed+1; col[2] = kGreen+2;
    cv->cd(1);
    TH1F* fr1 = gPad->DrawFrame(0, 0, 30, 0.10);
    fr1->GetXaxis()->SetTitle("r [cm]");
    fr1->GetYaxis()->SetTitle("wedge structure RMS of r#Delta#phi [cm]");
    for (int dh = 0; dh < 3; dh++){
        grRdp[dh]->SetMarkerStyle(20); grRdp[dh]->SetMarkerColor(col[dh]);
        grRdp[dh]->SetLineColor(col[dh]); grRdp[dh]->SetLineWidth(2);
        grRdp[dh]->Draw("LP same");
    }
    cv->cd(2);
    TH1F* fr2 = gPad->DrawFrame(0, 0, 30, 10.0);
    fr2->GetXaxis()->SetTitle("r [cm]");
    fr2->GetYaxis()->SetTitle("wedge structure RMS of #Delta#phi [mrad]");
    for (int di = 0; di < 3; di++){
        grDp[di]->SetMarkerStyle(20); grDp[di]->SetMarkerColor(col[di]);
        grDp[di]->SetLineColor(col[di]); grDp[di]->SetLineWidth(2);
        grDp[di]->Draw("LP same");
    }
    cv->SaveAs(Form("%s/fstWedgeStructure_vsR.png", outdir));
    printf("\nWrote %s/fstWedgeStructure_vsR.png\n", outdir);
}
