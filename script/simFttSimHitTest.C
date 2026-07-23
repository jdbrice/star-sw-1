// simFttSimHitTest.C
//
// Test driver for StFttSimHitMaker (see proposal_ftt_sim_maker.txt, option 6,
// and status_FttSimHitMaker.txt). Reads a particle-gun .fzd (produced by
// runSimFlat.C), reconstructs FTT through the REAL StFttClusterMaker /
// StFttClusterPointMaker chain (mUseGeantData=false) fed by StFttSimHitMaker
// instead of the GEANT-truth blur (MakeGeantPoints), then compares the
// resulting StFttPoint global positions directly against GEANT truth
// (St_g2t_fts_hit) event by event.
//
// Usage:
//   root4star -b -q 'script/simFttSimHitTest.C("ele.pt1.5.vz0.run1.fzd", 5)'

#include <map>
#include <utility>

class StFttSimHitMaker;
class StFttClusterMaker;
class StFttClusterPointMaker;
class StFttDbMaker;

void simFttSimHitTest(const char* inFile = "ele.pt1.5.vz0.run1.fzd", int nevents = 5, int debug = 0) {

    TString _geom = ""; // use the fGeom.root cache
    TString _chain = Form("fzin %s sdt20211016 MakeEvent StEvent McEvent bigbig", _geom.Data());
    printf("Chain: \n%s\n", _chain.Data());

    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile);

    cout << "Using the Geometry cache: fGeom.root (empty _geom tag + cwd fGeom.root)" << endl;

    // ------------------------------------------------------------------
    // Wire in the FTT chain manually (same makers/order as
    // fwd_afterburner_db.C's real-data chain), with StFttSimHitMaker in
    // place of StFttRawHitMaker/StFttHitCalibMaker:
    //   StFttDbMaker -> StFttSimHitMaker -> StFttClusterMaker -> StFttClusterPointMaker
    // ------------------------------------------------------------------
    gSystem->Load("libStFttDbMaker.so");
    gSystem->Load("libStFttSimHitMaker.so");
    gSystem->Load("libStFttClusterMaker.so");
    gSystem->Load("libStFttClusterPointMaker.so");

    StFttDbMaker* fttDbMk = new StFttDbMaker();
    chain->AddMaker(fttDbMk);

    StFttSimHitMaker* fttSimHit = new StFttSimHitMaker();
    fttSimHit->SetDebug(debug);
    chain->AddMaker(fttSimHit);

    StFttClusterMaker* fttClu = new StFttClusterMaker();
    fttClu->SetDebug(debug);
    fttClu->SetTimeCut(1 /*kTimeCutModeAcceptAll*/, -9999, 9999);
    chain->AddMaker(fttClu);

    StFttClusterPointMaker* fttCP = new StFttClusterPointMaker();
    fttCP->SetDebug(debug);
    // mUseGeantData left at its constructor default (false) -- take the REAL
    // MakeLocalPoints/MakeGlobalPoints path, the whole point of this test.
    chain->AddMaker(fttCP);

    Int_t iInit = chain->Init();
    cout << "CHAIN INIT DONE? (good==0): " << iInit << endl;
    if (iInit) chain->Fatal(iInit, "on init");
    chain->PrintInfo();

    int nClustersDetailPrinted = 0; // cap the per-strip breakdown at 10 clusters total

    for (int i = 0; i < nevents; i++) {
        cout << "--------->START EVENT: " << i << endl;
        chain->Clear();
        if (kStOK != chain->Make()) break;

        StEvent* event = (StEvent*)chain->GetInputDS("StEvent");
        St_g2t_fts_hit* g2t = (St_g2t_fts_hit*)chain->GetDataSet("geant/g2t_stg_hit");

        if (!event || !event->fttCollection() || !g2t) {
            cout << "  missing StEvent/fttCollection/g2t -- skipping analysis this event" << endl;
            continue;
        }

        StFttCollection* fttColl = event->fttCollection();
        cout << "  nRawHits=" << fttColl->rawHits().size()
             << " nClusters=" << fttColl->clusters().size()
             << " nPoints=" << fttColl->points().size() << endl;

        // Build truth map: (plane_id, quadrant_id, orientation) -> (x,y,z) cm,
        // same volume_id decoding as StFttSimHitMaker / MakeGeantPoints
        std::map<std::pair<int,int>, int> track_vol_count;
        for (int k = 0; k < g2t->GetNRows(); k++) {
            g2t_fts_hit_st* git = (g2t_fts_hit_st*)g2t->At(k);
            if (!git) continue;
            int track_id = git->track_p;
            int volume_id = git->volume_id;
            if (++track_vol_count[std::make_pair(track_id, volume_id)] > 1) continue;
            int plane_id = (volume_id - 1) / 100;
            int quadrant_id = ((volume_id - (100 * plane_id)) / 10) - 1;
            int orientation = (volume_id % 2 == 0) ? 1 /*kFttVertical*/ : 0 /*kFttHorizontal*/;
            printf("  TRUTH track=%d plane=%d quad=%d orient=%d  global=(%7.3f,%7.3f,%7.3f) cm\n",
                   track_id, plane_id, quadrant_id, orientation, git->x[0], git->x[1], git->x[2]);
        }

        for (size_t ic = 0; ic < fttColl->clusters().size(); ic++) {
            StFttCluster* c = fttColl->clusters()[ic];
            printf("  CLUSTER plane=%d quad=%d row=%d orient=%d nStrips=%d sumAdc=%.0f x=%.2fmm sigma=%.3fmm maxStripLength=%.2f\n",
                   c->plane(), c->quadrant(), c->row(), c->orientation(), c->nStrips(), c->sumAdc(), c->x(), c->sigma(), c->maxStripLength());
            if (nClustersDetailPrinted < 10) {
                printf("    strips:");
                for (size_t ih = 0; ih < c->rawHits().size(); ih++) {
                    StFttRawHit* rh = c->rawHits()[ih];
                    printf(" [strip=%d adc=%d]", (int)rh->strip(), (int)rh->adc());
                }
                printf("\n");
                nClustersDetailPrinted++;
            }
        }

        for (size_t ip = 0; ip < fttColl->points().size(); ip++) {
            StFttPoint* p = fttColl->points()[ip];
            printf("  POINT   plane=%d quad=%d  global=(%7.3f,%7.3f,%7.3f) cm  nClusters=%d\n",
                   p->plane(), p->quadrant(), p->xyz().x(), p->xyz().y(), p->xyz().z(),
                   p->nClusters());
        }

        cout << "<---------- END EVENT" << endl;
    }
}
