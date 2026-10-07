// chargeAsym.C -- per-quadrant charge asymmetry per track type.
//
//   root -l -b -q 'script/chargeAsym.C("/gpfs01/star/pwg_tasks/FwdCalib/akio/qc_20261007")'
//
// Input: merged_<type>.root written by hadd from the afterburner's
// FwdDetResidual_<type>.root files, e.g.
//   hadd -f $O/merged_FCSConstrained.root $W/*FwdDetResidual_FCSConstrained.root
// (do NOT pass -j N to hadd -- the N is parsed as an input filename.)
//
// The asymmetry is read from the FST1 bin of PlaneUsageQuadCharge. NOT bin 1:
// bin 1 is "Vtx", which is empty for Global tracks because they have no vertex,
// and reading it silently gives an all-zero Global row.

void chargeAsym( TString dir = "/gpfs01/star/pwg_tasks/FwdCalib/akio/qc_20261007",
                 TString plane = "FST1" ) {
    const char* types[5] = {"Global","BLC","BLCVtx","Primary","FCSConstrained"};
    const char* quads[4] = {"Ntop","Nbot","Stop","Sbot"};

    printf("\n  Charge asymmetry  A = (N+ - N-)/(N+ + N-)   from the %s bin\n", plane.Data());
    printf("  %-16s %9s %9s %9s %9s %11s\n","type",quads[0],quads[1],quads[2],quads[3],"N-S swing");
    for (int t = 0; t < 5; t++) {
        TFile* f = TFile::Open(Form("%s/merged_%s.root", dir.Data(), types[t]));
        if (!f || f->IsZombie()) { printf("  %-16s  (no merged file)\n", types[t]); continue; }
        double a[4] = {0,0,0,0};
        bool ok = true;
        printf("  %-16s", types[t]);
        for (int q = 0; q < 4; q++) {
            TH1F* hp = (TH1F*)f->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Pos",quads[q]));
            TH1F* hn = (TH1F*)f->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Neg",quads[q]));
            if (!hp || !hn) { printf(" %9s","--"); ok = false; continue; }
            int b = hp->GetXaxis()->FindBin(plane.Data());
            double np = hp->GetBinContent(b), nn = hn->GetBinContent(b);
            a[q] = (np + nn > 0) ? (np - nn)/(np + nn) : 0;
            printf(" %+9.4f", a[q]);
        }
        // south mean minus north mean
        if (ok) printf(" %+11.4f\n", 0.5*(a[2]+a[3]) - 0.5*(a[0]+a[1]));
        else    printf("\n");
        f->Close();
    }

    printf("\n  raw %s counts, millions (q+ / q-)\n", plane.Data());
    for (int t = 0; t < 5; t++) {
        TFile* f = TFile::Open(Form("%s/merged_%s.root", dir.Data(), types[t]));
        if (!f || f->IsZombie()) continue;
        printf("  %-16s", types[t]);
        for (int q = 0; q < 4; q++) {
            TH1F* hp = (TH1F*)f->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Pos",quads[q]));
            TH1F* hn = (TH1F*)f->Get(Form("PlaneUsageQuadCharge/hPlaneUsage_%s_Neg",quads[q]));
            if (!hp || !hn) continue;
            int b = hp->GetXaxis()->FindBin(plane.Data());
            printf("  %s %.1f/%.1f", quads[q], hp->GetBinContent(b)/1e6, hn->GetBinContent(b)/1e6);
        }
        printf("\n");
        f->Close();
    }
    printf("\n");
}
