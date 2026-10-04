// Probe: what FCS detector offsets does the DB give for a given run/timestamp,
// versus the compiled-in values the analysis has been using (setDbAccess(0))?
void fcsoff(int run, int date, int time) {
    gSystem->Load("libPhysics"); gSystem->Load("St_base"); gSystem->Load("StChain");
    gSystem->Load("St_Tables"); gSystem->Load("StUtilities"); gSystem->Load("StEvent");
    gSystem->Load("StDbLib"); gSystem->Load("StDbBroker"); gSystem->Load("libStDb_Tables");
    gSystem->Load("St_db_Maker"); gSystem->Load("StFcsDbMaker");

    StChain* chain = new StChain("fcsprobe");
    St_db_Maker* dbMk = new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb", "StarDb");
    dbMk->SetDateTime(date, time);
    StFcsDbMaker* mk = new StFcsDbMaker();
    mk->setDbAccess(1);
    chain->Init();
    dbMk->InitRun(run);
    int st = mk->InitRun(run);
    printf("  StFcsDbMaker::InitRun returned %d\n", st);
    StFcsDb* db = (StFcsDb*)mk->GetDataSet("fcsDb");
    if (!db) { printf("  no fcsDb\n"); return; }
    printf("\n  RUN %d  (db timestamp %d %06d)  mDbAccess=1\n", run, date, time);
    for (int i = 0; i < 4; i++) {
        StThreeVectorD o = db->getDetectorOffset(i);
        printf("    det=%d  %9.4f %9.4f %9.4f\n", i, o.x(), o.y(), o.z());
    }
}
