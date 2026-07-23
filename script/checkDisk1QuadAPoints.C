// checkDisk1QuadAPoints.C
//
// Follow-up to checkDisk1QuadA.C: feb 2,4,6 (Vertical/X strips) are
// confirmed hardware-dead for disk1(plane1)/quadA (zero raw hits). This
// checks what StFttClusterPointMaker actually builds there -- specifically
// whether StFttPoints still get created from the alive H(Y) clusters
// (feb 1,3,5), each carrying a real Y but only a maxStripCenter() X
// placeholder (StFttClusterPointMaker.cxx:216), and dumps the actual X
// values so we can see the real placeholder spread/location, not just
// reason about it from the code.
//
// Usage: root4star -l -b -q 'checkDisk1QuadAPoints.C("file.MuDst.root", 2000)'

void loadLibsCD1QAP();

void checkDisk1QuadAPoints(const Char_t *fileList = "st_fwd_23081004_raw_6000003.MuDst.root",
                            size_t nEvents = 2000) {
    loadLibsCD1QAP();

    StChain *chain = new StChain("StChain");
    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", 1);
    TChain &muDstChain = *muDstMaker.chain();
    printf("MuDst file has %d events available in tree\n", muDstChain.GetEntries());

    St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
    StMuDst2StEventMaker *mu2ev = new StMuDst2StEventMaker();

    StFttDbMaker *fttDbMk = new StFttDbMaker();
    chain->AddMaker(fttDbMk);
    StFttHitCalibMaker *ftthcm = new StFttHitCalibMaker();
    StFttClusterMaker *fttclu = new StFttClusterMaker();
    fttclu->SetTimeCut(2, -40, 100); // matches current production config
    StFttClusterPointMaker *fttCP = new StFttClusterPointMaker();

    Int_t iInit = chain->Init();
    if (iInit) chain->Fatal(iInit, "on init");

    size_t nEntries = muDstChain.GetEntries();
    if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

    const int targetPlane = 1; // disk1, 0-based
    const int targetQuad  = 0; // quadA, 0-based

    int nEventsSeen = 0;
    long nPointsTotal = 0;
    long nPointsHOnly = 0; // nClusters==1, that cluster is Horizontal (Y-loop-built point)
    long nPointsVOnly = 0; // nClusters==1, that cluster is Vertical (X-loop-built point)
    int nPrinted = 0;

    TH1F *hXlocal  = new TH1F("hXlocal",  "placeholder local x (cell units), disk1/quadA;local x [cell];points", 200, 0, 400);
    TH1F *hXglobal = new TH1F("hXglobal", "placeholder global X, disk1/quadA;global X [cm];points", 200, -10, 50);

    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) { printf("event %d: Make() non-OK, stopping\n", (int)iev); break; }
        nEventsSeen++;

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttPoint &pts = event->fttCollection()->points();
        for (size_t ip = 0; ip < pts.size(); ip++) {
            StFttPoint *pt = pts[ip];
            if (pt->plane() != targetPlane || pt->quadrant() != targetQuad) continue;
            nPointsTotal++;

            if (pt->nClusters() == 1) {
                // StFttPoint::cluster(i) is indexed by ORIENTATION SLOT
                // (0=H,1=V,2=DiagH,3=DiagV), not insertion order -- must
                // check each slot for null, not assume cluster(0) is "the"
                // single cluster.
                StFttCluster *cH = pt->cluster(0);
                StFttCluster *cV = pt->cluster(1);
                if (cH != 0) {
                    nPointsHOnly++;
                    hXlocal->Fill(pt->x());
                    hXglobal->Fill(pt->xyz().x());
                    if (nPrinted < 30) {
                        printf("evt=%d local x=%.3f y=%.3f sigX=%.3f sigY=%.3f global=(%.3f,%.3f,%.3f)\n",
                               (int)iev, pt->x(), pt->y(), pt->sigmaX(), pt->sigmaY(),
                               pt->xyz().x(), pt->xyz().y(), pt->xyz().z());
                        nPrinted++;
                    }
                } else if (cV != 0) {
                    nPointsVOnly++;
                }
            }
        }
    }
    printf("\nprocessed %d events\n", nEventsSeen);
    printf("disk1/quadA points: total=%ld  H-only(X=placeholder)=%ld  V-only=%ld\n",
           nPointsTotal, nPointsHOnly, nPointsVOnly);
    printf("placeholder local x:  mean=%.3f rms=%.3f range=[%.1f,%.1f]\n",
           hXlocal->GetMean(), hXlocal->GetRMS(), hXlocal->GetXaxis()->GetXmin(), hXlocal->GetXaxis()->GetXmax());
    printf("placeholder global X: mean=%.3f rms=%.3f  underflow=%.0f overflow=%.0f\n",
           hXglobal->GetMean(), hXglobal->GetRMS(), hXglobal->GetBinContent(0), hXglobal->GetBinContent(hXglobal->GetNbinsX()+1));
    printf("global X histogram (nonzero 1cm-ish bins):\n");
    for (int b = 1; b <= hXglobal->GetNbinsX(); b++) {
        if (hXglobal->GetBinContent(b) > 0)
            printf("  X=%.2f : %.0f\n", hXglobal->GetBinCenter(b), hXglobal->GetBinContent(b));
    }
}

void loadLibsCD1QAP() {
    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();

    gSystem->Load("StarMagField");
    gSystem->Load("StMagF");
    gSystem->Load("StDetectorDbMaker");
    gSystem->Load("StTpcDb");
    gSystem->Load("StDaqLib");
    gSystem->Load("StDbBroker");
    gSystem->Load("StDbUtilities");
    gSystem->Load("St_db_Maker");

    gSystem->Load("StEvent");
    gSystem->Load("StEventMaker");

    gSystem->Load("St_base.so");
    gSystem->Load("StUtilities.so");
    gSystem->Load("libPhysics.so");

    gSystem->Load("StarClassLibrary");
    gSystem->Load("StStrangeMuDstMaker");
    gSystem->Load("StMuDSTMaker");

    gSystem->Load("StFttDbMaker");
    gSystem->Load("StFttHitCalibMaker");
    gSystem->Load("StFttClusterMaker");
    gSystem->Load("StFttClusterPointMaker");

    gSystem->Load("StStarLogger.so");
}
