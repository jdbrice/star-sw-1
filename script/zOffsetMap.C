// zOffsetMap.C -- measure the signed z offset of each sensitive layer from the
// nominal plane z that reconstruction currently assigns, for BOTH detectors.
//
// This is the constant the fix needs, and it must be MEASURED, not read off the
// XML placement order: nothing guarantees that the g2t volume_id's layer digit
// follows the order the <Placement> lines appear in, nor that a quadrant's
// rotation preserves the sign of local z.
//
// FTT: per (station, quadrant, layer) from geant/g2t_stg_hit, against
//      StFttDb::idealPlaneZLocations_Quad[A-D].
// FST: per (disk, wedge, sensor) from geant/g2t_fsi_hit, against the one z per
//      disk that StFstHitMaker.cxx:167-169 hardcodes. Also prints the z that the
//      StFstDb sensor matrix already carries, which decides whether the FST fix
//      is "add an offset" or "stop inventing z".
//
// Usage: root4star -b -q 'script/zOffsetMap.C("file.fzd", 400)'
void zOffsetMap(const char* fzd = "ele2pi_np1.dev2022sm.vz0.run3.fzd",
                int nev = 400,
                const char* tag = "dev2022sm",
                const char* sdt = "20211112.000004")
{
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");
    TString ch = Form("fzin sdt%s %s agml usexgeom misalign geant gstar "
                      "fstSlowSim fstRawHit fstCluster fstHit "
                      "MakeEvent StEvent McEvent", sdt, tag);
    bfc(-1, ch, fzd);
    StChain* chain = (StChain*)StMaker::GetChain();
    gSystem->Load("libStFttDbMaker.so");
    chain->AddMaker(new StFttDbMaker());
    if (chain->Init() != kStOK) { printf("zOffsetMap: Init FAILED\n"); return; }

    StFttDb* fttDb = (StFttDb*)chain->GetDataSet("fttDb");
    printf("zOffsetMap: fttDb = %p\n", (void*)fttDb);

    // --- accumulators. Fixed-size arrays, no local std::map (CINT). ---
    // FTT: [station 0-3][quad 0-3][layer 0-1]  (layer index = layerDigit-1)
    double fSum[4][4][2], fSum2[4][4][2];
    int    fN[4][4][2];
    for (int ia = 0; ia < 4; ia++) for (int ib = 0; ib < 4; ib++) for (int ic = 0; ic < 2; ic++)
        { fSum[ia][ib][ic] = 0; fSum2[ia][ib][ic] = 0; fN[ia][ib][ic] = 0; }

    // FST: [disk 0-2][wedge 0-11 within disk][sensor 0-2]
    double sSum[3][12][3], sSum2[3][12][3];
    int    sN[3][12][3];
    for (int id = 0; id < 3; id++) for (int iw = 0; iw < 12; iw++) for (int is = 0; is < 3; is++)
        { sSum[id][iw][is] = 0; sSum2[id][iw][is] = 0; sN[id][iw][is] = 0; }

    const double kFttZ[4][4] = {  // [quad A,B,C,D][station] from StFttDb
        {312.342, 329.953, 347.637, 365.422},
        {312.165, 330.094, 347.543, 365.409},
        {312.317, 329.723, 347.346, 365.311},
        {312.098, 329.731, 347.444, 365.455} };
    const double kFstZ[3] = { 151.750, 165.248, 178.781 };  // StFstHitMaker.cxx:167-169

    for (int iev = 0; iev < nev; iev++) {
        if (chain->Make() != kStOK) break;

        St_g2t_fts_hit* stg = (St_g2t_fts_hit*)chain->GetDataSet("geant/g2t_stg_hit");
        if (stg && fttDb) {
            for (int ks = 0; ks < stg->GetNRows(); ks++) {
                g2t_fts_hit_st* h = (g2t_fts_hit_st*)stg->At(ks);
                int vid = h->volume_id;
                int st  = fttDb->getStationFromG2TVolumeId(vid);
                int qd  = fttDb->getQuadFromG2TVolumeId(vid);
                int ly  = fttDb->getLayerFromG2TVolumeId(vid);
                if (st < 0 || qd < 0 || ly < 1 || ly > 2) continue;
                double dz = TMath::Abs(h->x[2]) - kFttZ[qd][st];
                fSum[st][qd][ly-1] += dz; fSum2[st][qd][ly-1] += dz*dz; fN[st][qd][ly-1]++;
            }
        }

        St_g2t_fts_hit* fsi = (St_g2t_fts_hit*)chain->GetDataSet("geant/g2t_fsi_hit");
        if (fsi) {
            for (int kf = 0; kf < fsi->GetNRows(); kf++) {
                g2t_fts_hit_st* h2 = (g2t_fts_hit_st*)fsi->At(kf);
                int vid2 = h2->volume_id;
                // FST volume_id digits: the disk digit is 4-6, see StFstDb
                int dk = StFstDb::getDiskFromG2TVolumeId ? 0 : 0;  // placeholder, set below
                dk = 0;
                // use the StFstDb decode through the dataset
                StFstDb* fstDb = (StFstDb*)chain->GetDataSet("fstDb");
                if (!fstDb) continue;
                dk = fstDb->getDiskFromG2TVolumeId(vid2);
                int wg = fstDb->getWedgeFromG2TVolumeId(vid2);
                int sn = fstDb->getSensorFromG2TVolumeId(vid2);
                if (dk < 1 || dk > 3 || wg < 1 || wg > 36 || sn < 0 || sn > 2) continue;
                int wInDisk = (wg - 1) % 12;
                double dz2 = TMath::Abs(h2->x[2]) - kFstZ[dk-1];
                sSum[dk-1][wInDisk][sn] += dz2; sSum2[dk-1][wInDisk][sn] += dz2*dz2;
                sN[dk-1][wInDisk][sn]++;
            }
        }
        chain->Clear();
    }

    printf("\n==================== FTT: g2t layer z - StFttDb nominal ====================\n");
    printf("volume_id = 100*station + 10*pentagon + layer; orientation = (volume_id%%2==0) ? V : H\n");
    printf(" st qd  layer   N      mean dz [cm]   rms\n");
    for (int js = 0; js < 4; js++) {
        for (int jq = 0; jq < 4; jq++) {
            for (int jl = 0; jl < 2; jl++) {
                int n = fN[js][jq][jl];
                if (!n) { printf("  %d  %d    %d   %5d        --\n", js, jq, jl+1, n); continue; }
                double m = fSum[js][jq][jl]/n;
                double r = TMath::Sqrt(TMath::Max(0.0, fSum2[js][jq][jl]/n - m*m));
                printf("  %d  %d    %d   %5d     %+8.4f     %6.4f\n", js, jq, jl+1, n, m, r);
            }
        }
    }

    printf("\n==================== FST: g2t sensor z - StFstHitMaker nominal ====================\n");
    printf(" disk wedgeInDisk sensor    N      mean dz [cm]   rms\n");
    for (int jd = 0; jd < 3; jd++) {
        for (int jw = 0; jw < 12; jw++) {
            for (int jn = 0; jn < 3; jn++) {
                int n2 = sN[jd][jw][jn];
                if (!n2) continue;
                double m2 = sSum[jd][jw][jn]/n2;
                double r2 = TMath::Sqrt(TMath::Max(0.0, sSum2[jd][jw][jn]/n2 - m2*m2));
                printf("   %d       %2d        %d     %5d     %+8.4f     %6.4f\n",
                       jd+1, jw, jn, n2, m2, r2);
            }
        }
    }

    // --- does the StFstDb sensor matrix already know the front/back z? ---
    printf("\n==================== FST: z carried by the StFstDb sensor matrix ====================\n");
    printf("If these already alternate front/back, the fix is to STOP inventing z,\n");
    printf("not to add an offset.  wedge: global z of sensor-local origin [cm]\n");
    for (int mw = 1; mw <= 36; mw++) {
        double l0[3] = {0,0,0}, g0[3] = {0,0,0};
        bool ok = StFstDb::sensorLocalToGlobal(mw, 0, l0, g0);
        double z0 = ok ? g0[2] : -999;
        bool ok1 = StFstDb::sensorLocalToGlobal(mw, 1, l0, g0);
        double z1 = ok1 ? g0[2] : -999;
        printf("  wedge %2d   sensor0 z %9.4f   sensor1 z %9.4f\n", mw, z0, z1);
    }
}
