// planarAB.C -- compare space-point and planar arms of the forward tracker.
//
//   root4star -l -b -q 'script/planarAB.C("A_sp.picoDst.root","A_pl.picoDst.root")'
//
// Reads StPicoFwdTrack from two picoDst files written by the SAME build with
// only TrackFitter:useSpacePoints differing, and prints chi2 and the number of
// seed points per track type. That is the number the planar-mode PR needs: the
// 2026-10-01 campaign measured Primary chi2 mean 4.332 -> 16.737 (+286%), but
// that predates the 17.55 mm FTT strip-offset fix and the Geometry_ftt DB
// tables, so it has to be re-made before it can be quoted.
//
// Track types are StFwdTrack::StFwdTrackType, packed in the low 3 bits of
// StPicoFwdTrack's vtxIndex byte:
//   0 Global  1 Beamline  2 Primary  3 FwdVtx  4 BLCVtx  5 FCSConstrained

void planarAB( const char *fSp, const char *fPl )
{
    gROOT->Macro("$STAR/StRoot/macros/loadMuDst.C");
    gSystem->Load("libStPicoEvent");
    gSystem->Load("libStPicoDstMaker");

    const char *tn[6] = { "Global", "Beamline", "Primary",
                          "FwdVtx", "BLCVtx", "FCSConstr" };
    // [arm][type]; arm 0 = space points, arm 1 = planar
    double sumChi2[2][6], sumSeed[2][6];
    long   n[2][6];
    for ( int ia = 0; ia < 2; ia++ )
        for ( int it = 0; it < 6; it++ ) { sumChi2[ia][it]=0; sumSeed[ia][it]=0; n[ia][it]=0; }

    const char *fn[2] = { fSp, fPl };
    for ( int iarm = 0; iarm < 2; iarm++ ) {
        TFile *f = TFile::Open( fn[iarm] );
        if ( !f || f->IsZombie() ) { printf("planarAB: cannot open %s\n", fn[iarm]); return; }
        TTree *t = (TTree*) f->Get("PicoDst");
        if ( !t ) { printf("planarAB: no PicoDst tree in %s\n", fn[iarm]); return; }

        TClonesArray *arr = 0;
        t->SetBranchAddress( "FwdTracks", &arr );   // "FwdTracks", plural
        Long64_t nev = t->GetEntries();
        printf("%s : %lld events\n", fn[iarm], nev);
        for ( Long64_t iev = 0; iev < nev; iev++ ) {
            t->GetEntry( iev );
            if ( !arr ) continue;
            for ( int itr = 0; itr < arr->GetEntriesFast(); itr++ ) {
                StPicoFwdTrack *tr = (StPicoFwdTrack*) arr->At( itr );
                if ( !tr ) continue;
                int ty = (int) tr->trackType();   // vtxIndexRaw() is not CINT-callable
                if ( ty < 0 || ty > 5 ) continue;
                sumChi2[iarm][ty] += tr->chi2();
                sumSeed[iarm][ty] += tr->numberOfSeedPoints();
                n[iarm][ty]++;
            }
        }
        f->Close();
    }

    printf("\n%-11s | %22s | %22s | %s\n", "track type",
           "space points", "planar", "ratio");
    printf("%-11s | %7s %7s %6s | %7s %7s %6s | %6s %6s\n", "",
           "N", "chi2", "nSeed", "N", "chi2", "nSeed", "chi2", "nSeed");
    printf("------------+------------------------+------------------------+--------------\n");
    for ( int it = 0; it < 6; it++ ) {
        if ( n[0][it] == 0 && n[1][it] == 0 ) continue;
        double c0 = n[0][it] ? sumChi2[0][it]/n[0][it] : 0;
        double c1 = n[1][it] ? sumChi2[1][it]/n[1][it] : 0;
        double s0 = n[0][it] ? sumSeed[0][it]/n[0][it] : 0;
        double s1 = n[1][it] ? sumSeed[1][it]/n[1][it] : 0;
        printf("%-11s | %7ld %7.3f %6.3f | %7ld %7.3f %6.3f | %6.2f %6.2f\n",
               tn[it], n[0][it], c0, s0, n[1][it], c1, s1,
               c0 > 0 ? c1/c0 : 0, s0 > 0 ? s1/s0 : 0);
    }
    printf("\nThe 2026-10-01 result to beat: Primary chi2 ratio was 3.86\n");
    printf("(4.332 space points -> 16.737 planar) on pre-17.55mm geometry.\n");
}
