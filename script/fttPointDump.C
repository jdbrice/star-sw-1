// fttPointDump.C -- one line per StFttPoint, for diffing two builds of the
// FTT reconstruction against each other.
//
// Runs only the FTT chain off a MuDst (StFttDbMaker -> StFttHitCalibMaker ->
// StFttClusterMaker -> StFttClusterPointMaker), no tracking, so it is cheap and
// the output is deterministic: clustering depends on strip INDICES only, so two
// builds that differ in the strip->position map produce the same points in the
// same order, and a line-by-line diff isolates the position change exactly.
//
// Columns: event, plane, quadrant, nClusters, lowest occupied orientation slot,
// local (x,y) in mm and global (x,y,z) in cm.
//
//   root4star -l -b -q 'script/fttPointDump.C("file.MuDst.root", 200, "pts.txt")'
void loadLibsPD();

void fttPointDump(const Char_t *fileList = "st_physics.MuDst.root",
                  size_t nEvents = 200, const char* outTxt = "ftt_points.txt")
{
    loadLibsPD();

    StChain *chain = new StChain("StChain");
    StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", fileList, "MuDst.root", 1);
    TChain &muDstChain = *muDstMaker.chain();
    printf("fttPointDump: %d events in tree\n", muDstChain.GetEntries());

    St_db_Maker *dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
    StMuDst2StEventMaker *mu2ev = new StMuDst2StEventMaker();

    StFttDbMaker *fttDbMk = new StFttDbMaker();
    chain->AddMaker(fttDbMk);
    StFttHitCalibMaker *ftthcm = new StFttHitCalibMaker();
    StFttClusterMaker *fttclu = new StFttClusterMaker();
    fttclu->SetTimeCut(2, -40, 100);   // production config
    StFttClusterPointMaker *fttCP = new StFttClusterPointMaker();

    Int_t iInit = chain->Init();
    if (iInit) chain->Fatal(iInit, "on init");

    FILE *fp = fopen(outTxt, "w");
    if (!fp) { printf("fttPointDump: cannot write %s\n", outTxt); return; }

    size_t nEntries = muDstChain.GetEntries();
    if (nEntries > nEvents && nEvents > 0) nEntries = nEvents;

    long nPts = 0;
    int  nEvtSeen = 0;
    for (size_t iev = 0; iev < nEntries; iev++) {
        chain->Clear();
        if (kStOK != chain->Make()) { printf("event %d: Make() non-OK, stopping\n", (int)iev); break; }
        nEvtSeen++;

        StEvent *event = (StEvent*)chain->GetInputDS("StEvent");
        if (!event || !event->fttCollection()) continue;

        StSPtrVecFttPoint &pts = event->fttCollection()->points();
        for (size_t ip = 0; ip < pts.size(); ip++) {
            StFttPoint *pt = pts[ip];
            // cluster(i) is indexed by ORIENTATION SLOT (0=H,1=V,2=DiagH,3=DiagV),
            // not insertion order, so every slot has to be probed for null.
            int slot = -1;
            for (int is = 0; is < 4; is++) {
                if (0 == pt->cluster(is)) continue;
                if (slot < 0) slot = is;
            }
            fprintf(fp, "%6d %d %d %d %2d %10.4f %10.4f   %9.4f %9.4f %9.4f\n",
                    (int)iev, (int)pt->plane(), (int)pt->quadrant(), pt->nClusters(), slot,
                    pt->x(), pt->y(), pt->xyz().x(), pt->xyz().y(), pt->xyz().z());
            nPts++;
        }
    }
    fclose(fp);
    printf("fttPointDump: %d events, %ld points -> %s\n", nEvtSeen, nPts, outTxt);
}

void loadLibsPD() {
    gSystem->Load("libStarClassLibrary.so");
    gSystem->Load("libStarRoot.so");
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    gSystem->Load("StarMagField");       gSystem->Load("StMagF");
    gSystem->Load("StDetectorDbMaker");  gSystem->Load("StTpcDb");
    gSystem->Load("StDaqLib");           gSystem->Load("StDbBroker");
    gSystem->Load("StDbUtilities");      gSystem->Load("St_db_Maker");
    gSystem->Load("StEvent");            gSystem->Load("StEventMaker");
    gSystem->Load("St_base.so");         gSystem->Load("StUtilities.so");
    gSystem->Load("libPhysics.so");      gSystem->Load("StarClassLibrary");
    gSystem->Load("StStrangeMuDstMaker");gSystem->Load("StMuDSTMaker");
    gSystem->Load("StFttDbMaker");       gSystem->Load("StFttHitCalibMaker");
    gSystem->Load("StFttClusterMaker");  gSystem->Load("StFttClusterPointMaker");
    gSystem->Load("StStarLogger.so");
}
