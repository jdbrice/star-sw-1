// testFttSlowSim.C -- one-command regression test for StFttSlowSimMaker.
//
//   root4star -b -q 'StRoot/StFttSimMaker/macro/testFttSlowSim.C("some.fzd")'
//   root4star -b -q 'StRoot/StFttSimMaker/macro/testFttSlowSim.C("some.fzd",200)'
//        ... 200 events instead of 1000, about a minute.
//   root4star -b -q 'StRoot/StFttSimMaker/macro/testFttSlowSim.C("some.fzd",300,900)'
//        ... third argument overrides mAdcPerKeV, to confirm the charge check
//            really does fail when the gain is wrong (900 should FAIL, while
//            the geometry checks still pass -- the two are independent).
//
// Needs any GEANT .fzd with sTGC hits; no other input. The real-data numbers
// it compares against are built in, so it runs anywhere.
//
// Runs the REAL chain  StFttDbMaker -> StFttSlowSimMaker -> StFttClusterMaker
// on a GEANT .fzd and checks two independent things, printing PASS/FAIL:
//
//   GEOMETRY  every truth hit should come back as a cluster at the same place.
//             Decodes volume_id the same way the maker does and compares the
//             cluster's measured coordinate against truth.
//   CHARGE    the cluster peak-strip ADC should match real FTT data. That is
//             what mAdcPerKeV is anchored on (the peak strip is the center
//             strip, which this +-1-strip model represents faithfully).
//
// sumAdc and the 4-7 strip tail are reported but NOT failed on: a +-1-strip
// model cannot make the wide clusters a shallow track leaves in real data.
// See ftt_sim_hit_maker/index.html#decharge .

// Real-data reference, measured on run 23 p+p FTT clusters (5.6M clusters).
// Hardcoded so this macro is self-contained; pass realRef to re-measure them
// from a fttChargeSharing.root if you have one (see
// StRoot/StFttDbMaker/macro and ftt_sim_hit_maker/charge_sharing.html).
const double kRealPeakMean   = 226.3;   // cluster peak-strip ADC
const double kRealPeakRms    = 195.8;
const double kRealSumAdcN2   = 508.3;   // cluster sumAdc, nStrips >= 2
const double kRealNStrips4Up = 20.8;    // % of clusters with 4 or more strips

void testFttSlowSim( const char* inFile = "ele.pt1.5.vz0.run1.fzd", int nevents = 1000,
                     float adcPerKeV = -1, const char* realRef = "" ) {
    TH1::AddDirectory( kFALSE );
    gRandom->SetSeed( 20261017 );   // reproducible

    TString _chain = "fzin sdt20211016 MakeEvent StEvent McEvent bigbig";
    gSystem->Load( "libStarRoot.so" );
    gROOT->LoadMacro( "bfc.C" );
    bfc( -1, _chain, inFile );
    gSystem->Load( "libStFttDbMaker.so" );
    gSystem->Load( "libStFttSimMaker.so" );
    gSystem->Load( "libStFttClusterMaker.so" );
    chain->AddMaker( new StFttDbMaker() );
    StFttSlowSimMaker* sim = new StFttSlowSimMaker();
    if ( adcPerKeV > 0 ) { sim->setAdcPerKeV( adcPerKeV );
                           printf("\n  (mAdcPerKeV overridden to %.0f)\n", adcPerKeV); }
    chain->AddMaker( sim );
    StFttClusterMaker* clu = new StFttClusterMaker();
    clu->SetTimeCut( 1 /*kTimeCutModeAcceptAll*/, -9999, 9999 );
    chain->AddMaker( clu );
    if ( chain->Init() ) { printf("\nFATAL: chain Init failed\n"); return; }

    const int kAgmlToDbQuadrant[4] = { 0, 3, 2, 1 };
    TH1D* hres = new TH1D( "hres", "", 400, -20, 20 );
    TH1D* hPk  = new TH1D( "hPk",  "", 200, 0, 1024 );   // SAME binning as the real hPeakAdc
    TH1D* hS2  = new TH1D( "hS2",  "", 200, 0, 4000 );
    TH1D* hNs  = new TH1D( "hNs",  "", 12, -0.5, 11.5 );
    long nTruth = 0, nMatched = 0, nRaw = 0;
    int nev = 0;

    while ( nev < nevents ) {
        if ( chain->Make() ) break;
        nev++;
        StEvent* ev  = (StEvent*)chain->GetInputDS( "StEvent" );
        StFttDb* db  = (StFttDb*)chain->GetDataSet( "fttDb" );
        St_g2t_fts_hit* t = (St_g2t_fts_hit*)chain->GetDataSet( "geant/g2t_stg_hit" );
        if ( !ev || !ev->fttCollection() || !db || !t ) { chain->Clear(); continue; }
        StFttCollection* fc = ev->fttCollection();
        nRaw += fc->rawHits().size();

        for ( size_t i = 0; i < fc->clusters().size(); i++ ) {
            StFttCluster* c = fc->clusters()[i];
            hNs->Fill( c->nStrips() );
            if ( c->nStrips() >= 2 ) hS2->Fill( c->sumAdc() );
            double pk = 0;
            for ( size_t ih = 0; ih < c->rawHits().size(); ih++ )
                if ( c->rawHits()[ih]->adc() > pk ) pk = c->rawHits()[ih]->adc();
            hPk->Fill( pk );
        }

        // Truth position is the dE-weighted centroid of the crossing's GEANT
        // steps, the same quantity the maker injects at -- the charge-weighted
        // centre is what the strips see. Collapse the step list first.
        const int kM = 4000; int ktr[kM], kvol[kM]; double kde[kM], kx[kM], ky[kM]; int nk = 0;
        for ( int ig = 0; ig < t->GetNRows(); ig++ ) {
            g2t_fts_hit_st* g = (g2t_fts_hit_st*)t->At( ig ); if ( !g ) continue;
            int f = -1;
            for ( int kq = 0; kq < nk; kq++ ) if ( ktr[kq]==g->track_p && kvol[kq]==g->volume_id ) { f = kq; break; }
            if ( f < 0 && nk < kM ) { f = nk++; ktr[f]=g->track_p; kvol[f]=g->volume_id; kde[f]=0; kx[f]=0; ky[f]=0; }
            if ( f >= 0 ) { kde[f] += g->de; kx[f] += g->de*g->x[0]; ky[f] += g->de*g->x[1]; }
        }
        for ( int ic2 = 0; ic2 < nk; ic2++ ) {
            if ( kde[ic2] <= 0 ) continue;
            double tx = kx[ic2]/kde[ic2], ty = ky[ic2]/kde[ic2];
            int st = kvol[ic2] / 100, aq = (kvol[ic2] / 10) % 10;
            if ( st < 0 || st > 3 || aq < 0 || aq > 3 ) continue;
            int dq = kAgmlToDbQuadrant[aq];
            int orient = ( kvol[ic2] % 2 == 0 ) ? 1 /*kFttVertical*/ : 0 /*kFttHorizontal*/;
            float dx, sx, dy, sy, dz, sz;
            db->getGloablOffset_ClusterPoint( (UChar_t)st, (UChar_t)dq, dx, sx, dy, sy, dz, sz );
            double lx = (tx*10.0 - dx)*sx, ly = (ty*10.0 - dy)*sy;
            if ( lx < 0 || ly < 0 ) continue;             // no strips under it, nothing to match
            double meas = ( orient == 1 ) ? lx : ly;      // V measures local x, H measures local y
            nTruth++;
            double best = 1e9;
            for ( size_t ic = 0; ic < fc->clusters().size(); ic++ ) {
                StFttCluster* c = fc->clusters()[ic];
                if ( c->plane() != st || c->quadrant() != dq ) continue;
                if ( (int)c->orientation() != orient ) continue;
                double d = c->x() - meas;
                if ( fabs(d) < fabs(best) ) best = d;
            }
            if ( best < 1e8 ) { hres->Fill( best ); if ( fabs(best) < 5.0 ) nMatched++; }
        }
        chain->Clear();
    }

    // real-data reference: the built-in numbers, or re-measured if a file was given
    double realPeakMean = kRealPeakMean, realPeakRms = kRealPeakRms;
    double realSumAdcN2 = kRealSumAdcN2, realNStrips4Up = kRealNStrips4Up;
    const char* refSource = "built-in (run 23 p+p)";
    if ( realRef && realRef[0] ) {
        TFile* fr = TFile::Open( realRef );
        TH1F* rPk = ( fr && !fr->IsZombie() ) ? (TH1F*)fr->Get("hPeakAdc")   : 0;
        TH1F* rS2 = ( fr && !fr->IsZombie() ) ? (TH1F*)fr->Get("hSumAdcN2p") : 0;
        TH1F* rNs = ( fr && !fr->IsZombie() ) ? (TH1F*)fr->Get("hNStrips")   : 0;
        if ( rPk && rS2 && rNs ) {
            realPeakMean = rPk->GetMean(); realPeakRms = rPk->GetRMS();
            realSumAdcN2 = rS2->GetMean();
            realNStrips4Up = 100.*rNs->Integral( rNs->FindBin(4), rNs->GetNbinsX() ) / rNs->GetEntries();
            refSource = realRef;
        } else {
            printf("\n  WARNING: could not read hPeakAdc/hSumAdcN2p/hNStrips from %s,"
                   " using the built-in numbers\n", realRef);
        }
    }

    double matchFrac = nTruth ? 100.*nMatched/nTruth : 0;
    double resMean = hres->GetMean(), resRms = hres->GetRMS();
    double pkMean = hPk->GetMean(), pkRms = hPk->GetRMS();

    printf("\n==================== StFttSlowSimMaker test ====================\n");
    printf("  input %s, %d events, %ld raw hits, %.0f clusters\n", inFile, nev, nRaw, hNs->GetEntries());
    printf("\n  GEOMETRY\n");
    printf("    truth hits on strips      %ld\n", nTruth);
    printf("    matched within 5 mm       %.2f%%            expect > 90%%      %s\n",
           matchFrac, matchFrac > 90.0 ? "PASS" : "FAIL");
    printf("    residual mean             %+.3f mm          expect |mean|<0.2 %s\n",
           resMean, fabs(resMean) < 0.2 ? "PASS" : "FAIL");
    printf("    residual rms              %.3f mm           (strip pitch 3.2 mm)\n", resRms);

    printf("\n  CHARGE  (peak-strip ADC is what mAdcPerKeV is anchored on)\n");
    double dm = 100.*(pkMean - realPeakMean)/realPeakMean;
    double dr = 100.*(pkRms  - realPeakRms )/realPeakRms;
    printf("    reference                 %s\n", refSource);
    printf("    peak ADC mean   MC %6.1f   real %6.1f   %+5.1f%%   expect |d|<5%% %s\n",
           pkMean, realPeakMean, dm, fabs(dm) < 5.0 ? "PASS" : "FAIL");
    printf("    peak ADC rms    MC %6.1f   real %6.1f   %+5.1f%%   expect |d|<5%% %s\n",
           pkRms, realPeakRms, dr, fabs(dr) < 5.0 ? "PASS" : "FAIL");

    printf("\n  FOR INFORMATION ONLY -- known limits of the +-1-strip model, not failures\n");
    printf("    sumAdc(nStrips>=2) MC %6.1f   real %6.1f   %+5.1f%%\n",
           hS2->GetMean(), realSumAdcN2, 100.*(hS2->GetMean()-realSumAdcN2)/realSumAdcN2);
    printf("    nStrips>=4 fraction  MC %5.1f%%   real %5.1f%%   (model stops at 3 by construction)\n",
           100.*hNs->Integral(hNs->FindBin(4),hNs->GetNbinsX())/hNs->GetEntries(), realNStrips4Up);
    printf("================================================================\n");
}
