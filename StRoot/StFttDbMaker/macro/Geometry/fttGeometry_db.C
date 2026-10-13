// fttGeometry_db.C -- push the Geometry/ftt survey tables to the offline MySQL DB.
//
//   root4star -b -q 'fttGeometry_db.C("run22iter3")'            dry run: read + print
//   root4star -b -q 'fttGeometry_db.C("run22iter3","writedb")'  actually store
//
// Run the dry form first. It prints every row it would store and touches nothing.
//
// The numbers are NOT duplicated here. They are read from ftt/*.C, the same files
// St_db_Maker reads when they are placed in ./StarDb/Geometry/ftt, so there is one
// copy of every number. See README for what each tag means.
//
// NAMESPACE: this macro writes Geometry_ftt. The tables used to live under
// Geometry/stgc; the rename to Geometry/ftt is PR #23 (code half) plus #20 (the
// StgmGeo1.xml <Misalign> paths). Do not push to Geometry_ftt until those have
// landed, or production will keep reading Geometry/stgc and never see these. Table
// NAMES are unchanged -- only the namespace moved -- so stgcOnTpc and
// stationOnStgc keep their stgc-flavoured names on purpose.
//
// beginTime follows Gene Van Buren's rule (database_timestamp.txt), per RUN YEAR:
//
//     run year   ideal        survey/estimate   real data
//     2021       2020-12-10   2020-12-15        any time after
//     2022       2021-10-10   2021-10-15        any time after
//     2023       2023-04-10   2023-04-15        any time after
//     2024       2023-12-10   2023-12-15        any time after
//     2025       2024-12-10   2024-12-15        any time after
//
// The Run 22 iterations sit at 2021-11-1x, which is the "real data" range and so
// correct for measured constants. The reserved Run 22 slots hold the identity
// entry; FTT has nothing for the survey slot, having no survey-only table.
// Still empty on purpose, and listed when the macro is given a bad tag: the Run
// 23 and Run 25 survey slots and the Run 25 ideal slot. Filling those decides
// what those runs reconstruct with, and an unpaired ideal entry hijacks its
// year -- see README.
//
// Unlike the FST tables, the three FTT tables do NOT all have an entry at every
// timestamp: pentOnStation was revised on its own several times. The map below
// therefore gives a stamp PER TABLE, picking the latest entry at or before the
// tag's time -- exactly what St_db_Maker would resolve to. Getting this wrong is
// the easy mistake here: storing pentOnStation alone leaves the other two at
// whatever was there before.

void fttGeometry_db( const char* tag = "", const char* opt = "" )
{
    gROOT->Macro("LoadLogger.C");
    gSystem->Load("St_base.so");
    gSystem->Load("libStDb_Tables.so");
    gSystem->Load("StDbLib.so");

    TString option(opt);
    const int writedb = option.Contains("writedb");

    const char *tables[3] = { "stgcOnTpc", "stationOnStgc", "pentOnStation" };
    TString stamp[3];     // per-table file stamp
    TString storeTime, what;

    TString t(tag);
    if ( t == "run22ideal" ) {
        stamp[0]="20211010.000000"; stamp[1]="20211010.000000"; stamp[2]="20211010.000000";
        storeTime = "2021-10-10 00:00:00";
        what = "Run 22 ideal: identity, no shift, no rotation";
    } else if ( t == "run22iter1" ) {
        stamp[0]="20211112.000001"; stamp[1]="20211112.000001"; stamp[2]="20211112.000001";
        storeTime = "2021-11-12 00:00:00";
        what = "Run 22 iteration 1: per-quadrant x/y from FST->sTGC residuals";
    } else if ( t == "run22iter2" ) {
        stamp[0]="20211112.000001"; stamp[1]="20211112.000001"; stamp[2]="20211112.000002";
        storeTime = "2021-11-12 00:00:02";
        what = "Run 22 iteration 2: adds per-quadrant z (t2); x/y unchanged from iter1";
    } else if ( t == "run22iter3" ) {
        stamp[0]="20211112.000001"; stamp[1]="20211112.000001"; stamp[2]="20211112.000003";
        storeTime = "2021-11-12 00:00:03";
        what = "Run 22 iteration 3: iter2 + measured deltas (FINAL Run 22)";
    } else if ( t == "run23ideal" ) {
        stamp[0]="20230410.000000"; stamp[1]="20230410.000000"; stamp[2]="20230410.000000";
        storeTime = "2023-04-10 00:00:00";
        what = "Run 23 ideal (link to the Run 22 ideal; was Daniel's 20221220 placeholder)";
    } else if ( t == "run24ideal" ) {
        stamp[0]="20231210.000000"; stamp[1]="20231210.000000"; stamp[2]="20231210.000000";
        storeTime = "2023-12-10 00:00:00";
        what = "Run 24 ideal (link to the Run 22 ideal)";
    } else if ( t == "run24seed" ) {
        stamp[0]="20231215.000000"; stamp[1]="20231215.000000"; stamp[2]="20231215.000000";
        storeTime = "2023-12-15 00:00:00";
        what = "Run 24 estimate: Run 22 finals copied verbatim, not a Run 24 measurement";
    } else if ( t == "run24iter1" ) {
        stamp[0]="20231215.000000"; stamp[1]="20231215.000000"; stamp[2]="20231215.000001";
        storeTime = "2023-12-15 00:00:01";
        what = "Run 24 iteration 1: first actual Run 24 measurement (90-job campaign)";
    } else {
        std::cout << "fttGeometry_db: unknown tag \"" << tag << "\"\n"
                  << "  known tags:\n"
                  << "    run22ideal   2021-10-10 00:00:00   identity\n"
                  << "    run22iter1   2021-11-12 00:00:00   per-quadrant x/y\n"
                  << "    run22iter2   2021-11-12 00:00:02   + per-quadrant z\n"
                  << "    run22iter3   2021-11-12 00:00:03   + deltas  (FINAL Run 22)\n"
                  << "    run23ideal   2023-04-10 00:00:00   identity (link)\n"
                  << "    run24ideal   2023-12-10 00:00:00   identity (link)\n"
                  << "    run24seed    2023-12-15 00:00:00   Run 22 finals carried over\n"
                  << "    run24iter1   2023-12-15 00:00:01   first Run 24 measurement\n"
                  << "\n  NOT SET, deliberately: the Run 23 and Run 25 SURVEY slots\n"
                  << "  (2023-04-15, 2024-12-15) and the Run 25 IDEAL slot (2024-12-10).\n"
                  << "  An ideal entry in a year with nothing later HIJACKS that year to\n"
                  << "  identity -- measured: a Run 25 ideal link took Run 25 from the Run\n"
                  << "  24 alignment to zeros. Add them only in ideal+estimate PAIRS.\n"
                  << "\n  Pass \"all\" to push every slot in order.\n";
        if ( t != "all" ) return;
    }

    if ( t == "all" ) {
        const char *every[8] = { "run22ideal","run22iter1","run22iter2","run22iter3",
                                 "run23ideal","run24ideal","run24seed","run24iter1" };
        for ( int i = 0; i < 8; i++ ) fttGeometry_db( every[i], opt );
        return;
    }

    std::cout << "fttGeometry_db: tag=" << tag << "\n"
              << "                beginTime=\"" << storeTime << "\"\n"
              << "                " << what << "\n"
              << "                writedb=" << writedb << "\n";

    TString here = gSystem->DirName( __FILE__ );

    for ( int it = 0; it < 3; it++ ) {

        TString fn = Form( "%s/ftt/%s.%s.C", here.Data(), tables[it], stamp[it].Data() );
        if ( gSystem->AccessPathName( fn.Data() ) ) {
            std::cout << "ERROR: cannot read " << fn << " -- nothing stored\n";
            return;
        }

        gROOT->ProcessLine( Form( ".L %s", fn.Data() ) );
        TDataSet *ds = (TDataSet *) gROOT->ProcessLineFast( "CreateTable()" );
        if ( !ds ) {
            std::cout << "ERROR: CreateTable() returned 0 for " << fn << "\n";
            gROOT->ProcessLine( Form( ".U %s", fn.Data() ) );
            return;
        }

        St_Survey *tableSet = (St_Survey *) ds;
        Survey_st *rows = tableSet->GetTable();
        const int nrows  = tableSet->GetNRows();

        std::cout << "\n  " << tables[it] << ": " << nrows << " row(s) from "
                  << gSystem->BaseName( fn.Data() ) << "\n";
        for ( int i = 0; i < nrows; i++ ) {
            printf("    Id=%3d  t=(%9.5f,%9.5f,%9.5f)\n",
                   rows[i].Id, rows[i].t0, rows[i].t1, rows[i].t2 );
        }

        if ( writedb ) {
            gSystem->Setenv("DB_ACCESS_MODE","write");
            std::cout << "    DB_ACCESS_MODE=" << gSystem->Getenv("DB_ACCESS_MODE") << "\n";
            StDbManager    *mgr  = StDbManager::Instance();
            StDbConfigNode *node = mgr->initConfig("Geometry_ftt");
            StDbTable      *table = node->addDbTable( tables[it] );
            mgr->setStoreTime( storeTime.Data() );
            table->SetTable( (char*)rows, nrows );
            mgr->storeDbTable( table );
            gSystem->Unsetenv("DB_ACCESS_MODE");
            std::cout << "    stored " << tables[it] << " at " << storeTime << "\n";
        }

        gROOT->ProcessLine( Form( ".U %s", fn.Data() ) );
    }

    if ( writedb ) std::cout << "\nDone with database upload.\n";
    else           std::cout << "\nDry run only -- pass \"writedb\" to store.\n";
}
