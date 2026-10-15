// fttGeometry_db.C -- push the Geometry/ftt survey tables to the offline MySQL DB.
//
//   root4star -b -q 'fttGeometry_db.C("run22iter3")'            dry run: read + print
//   root4star -b -q 'fttGeometry_db.C("run22iter3","writedb")'  actually store
//   root4star -b -q 'fttGeometry_db.C("all","writedb")'         store every slot
//   root4star -b -q 'fttGeometry_db.C("run22ideal","writedb","stgcOnTpc")'
//                                                     one table only, for a first push
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
// Every entry for a run year sits on that year's ideal or estimate date, with the
// time of day carrying the rest: minute 00 is the baseline kept in this directory
// (00:00:00 reserved against a real sTGC survey, 00:00:01+ our best guesses),
// minute 01 is the alignment program's lane. See README.
//
// Unlike the FST tables, the three FTT tables do NOT all have an entry at every
// timestamp: pentOnStation was revised on its own several times. The map below
// therefore gives a stamp PER TABLE, picking the latest entry at or before the
// tag's time -- exactly what St_db_Maker would resolve to. Getting this wrong is
// the easy mistake here: storing pentOnStation alone leaves the other two at
// whatever was there before.

// Returns false if anything failed, so "all" stops instead of reporting success
// it did not get.
bool fttGeometry_db( const char* tag = "", const char* opt = "", const char* onlyTable = "" )
{
    gROOT->Macro("LoadLogger.C");
    gSystem->Load("St_base.so");
    gSystem->Load("libStDb_Tables.so");
    gSystem->Load("StDbLib.so");

    TString option(opt);
    const int writedb = option.Contains("writedb");

    // onlyTable: restrict to ONE of the three tables. Only for a first push, to
    // test write permission on the smallest possible change. All three are
    // chained by StFttDb::loadGeometry(), so a partial update leaves the
    // geometry inconsistent at that timestamp -- finish the others at once.
    TString only(onlyTable);
    if ( only.Length() ) std::cout << "\n*** onlyTable=\"" << only
                                   << "\" -- partial, finish the other two ***\n";

    const char *tables[3] = { "stgcOnTpc", "stationOnStgc", "pentOnStation" };
    TString stamp[3];     // per-table file stamp
    TString storeTime, what;

    TString t(tag);
    if ( t == "run22ideal" ) {
        stamp[0]="20211010.000000"; stamp[1]="20211010.000000"; stamp[2]="20211010.000000";
        storeTime = "2021-10-10 00:00:00";
        what = "Run 22 ideal: identity, no shift, no rotation";
    } else if ( t == "run22iter1" ) {
        stamp[0]="20211015.000001"; stamp[1]="20211015.000001"; stamp[2]="20211015.000001";
        storeTime = "2021-10-15 00:00:01";
        what = "Run 22 iteration 1: per-quadrant x/y from FST->sTGC residuals";
    } else if ( t == "run22iter2" ) {
        stamp[0]="20211015.000001"; stamp[1]="20211015.000001"; stamp[2]="20211015.000002";
        storeTime = "2021-10-15 00:00:02";
        what = "Run 22 iteration 2: adds per-quadrant z (t2); x/y unchanged from iter1";
    } else if ( t == "run22iter3" ) {
        stamp[0]="20211015.000001"; stamp[1]="20211015.000001"; stamp[2]="20211015.000003";
        storeTime = "2021-10-15 00:00:03";
        what = "Run 22 iteration 3: iter2 + measured deltas (FINAL Run 22)";
    } else if ( t == "run23ideal" ) {
        stamp[0]="20230410.000000"; stamp[1]="20230410.000000"; stamp[2]="20230410.000000";
        storeTime = "2023-04-10 00:00:00";
        what = "Run 23 ideal (link to the Run 22 ideal; was Daniel's 20221220 placeholder)";
    } else if ( t == "run23est" ) {
        stamp[0]="20230415.000001"; stamp[1]="20230415.000001"; stamp[2]="20230415.000001";
        storeTime = "2023-04-15 00:00:01";
        what = "Run 23 estimate: Run 22 iter 3 carried over (best guess, no Run 23 measurement)";
    } else if ( t == "run24ideal" ) {
        stamp[0]="20231210.000000"; stamp[1]="20231210.000000"; stamp[2]="20231210.000000";
        storeTime = "2023-12-10 00:00:00";
        what = "Run 24 ideal (link to the Run 22 ideal)";
    } else if ( t == "run24est" ) {
        stamp[0]="20231215.000001"; stamp[1]="20231215.000001"; stamp[2]="20231215.000001";
        storeTime = "2023-12-15 00:00:01";
        what = "Run 24 estimate: Run 22 iter 3 carried over, not a Run 24 measurement";
    } else if ( t == "run24iter1" ) {
        stamp[0]="20231215.000001"; stamp[1]="20231215.000001"; stamp[2]="20231215.000002";
        storeTime = "2023-12-15 00:00:02";
        what = "Run 24 iteration 1: first actual Run 24 measurement (90-job campaign)";
    } else if ( t == "run25ideal" ) {
        stamp[0]="20241210.000000"; stamp[1]="20241210.000000"; stamp[2]="20241210.000000";
        storeTime = "2024-12-10 00:00:00";
        what = "Run 25 ideal (link to the Run 22 ideal)";
    } else if ( t == "run25est" ) {
        stamp[0]="20241215.000001"; stamp[1]="20241215.000001"; stamp[2]="20241215.000001";
        storeTime = "2024-12-15 00:00:01";
        what = "Run 25 estimate: Run 24 iter 1 carried over (best guess)";
    } else {
        std::cout << "fttGeometry_db: unknown tag \"" << tag << "\"\n"
                  << "  known tags:\n"
                  << "    run22ideal   2021-10-10 00:00:00   identity\n"
                  << "    run22iter1   2021-10-15 00:00:01   per-quadrant x/y\n"
                  << "    run22iter2   2021-10-15 00:00:02   + per-quadrant z\n"
                  << "    run22iter3   2021-10-15 00:00:03   + deltas  (BEST for Run 22)\n"
                  << "    run23ideal   2023-04-10 00:00:00   identity (link)\n"
                  << "    run23est     2023-04-15 00:00:01   Run 22 iter 3 carried over\n"
                  << "    run24ideal   2023-12-10 00:00:00   identity (link)\n"
                  << "    run24est     2023-12-15 00:00:01   Run 22 iter 3 carried over\n"
                  << "    run24iter1   2023-12-15 00:00:02   first Run 24 measurement\n"
                  << "    run25ideal   2024-12-10 00:00:00   identity (link)\n"
                  << "    run25est     2024-12-15 00:00:01   Run 24 iter 1 carried over\n"
                  << "\n  The 00:00:0x lane is the baseline kept in this directory; the\n"
                  << "  00:01:0x lane on the same dates is the alignment program's.\n"
                  << "  Reserved and empty: 00:00:00 on each estimate date, against the\n"
                  << "  day there is a real sTGC survey.\n"
                  << "  Add ideal and estimate slots only in PAIRS -- an ideal entry in a\n"
                  << "  year with nothing later HIJACKS that year to identity.\n"
                  << "\n  Pass \"all\" to push every slot in order.\n";
        if ( t != "all" ) return false;
    }

    if ( t == "all" ) {
        const char *every[11] = { "run22ideal","run22iter1","run22iter2","run22iter3",
                                  "run23ideal","run23est",
                                  "run24ideal","run24est","run24iter1",
                                  "run25ideal","run25est" };
        for ( int i = 0; i < 11; i++ ) {
            if ( !fttGeometry_db( every[i], opt, onlyTable ) ) {
                std::cout << "\n*** STOPPING at tag \"" << every[i]
                          << "\" -- nothing further was attempted ***\n";
                return false;
            }
        }
        return true;
    }

    std::cout << "fttGeometry_db: tag=" << tag << "\n"
              << "                beginTime=\"" << storeTime << "\"\n"
              << "                " << what << "\n"
              << "                writedb=" << writedb << "\n";

    TString here = gSystem->DirName( __FILE__ );

    for ( int it = 0; it < 3; it++ ) {

        if ( only.Length() && only != tables[it] ) continue;

        TString fn = Form( "%s/ftt/%s.%s.C", here.Data(), tables[it], stamp[it].Data() );
        if ( gSystem->AccessPathName( fn.Data() ) ) {
            std::cout << "ERROR: cannot read " << fn << " -- nothing stored\n";
            return false;
        }

        gROOT->ProcessLine( Form( ".L %s", fn.Data() ) );
        TDataSet *ds = (TDataSet *) gROOT->ProcessLineFast( "CreateTable()" );
        if ( !ds ) {
            std::cout << "ERROR: CreateTable() returned 0 for " << fn << "\n";
            gROOT->ProcessLine( Form( ".U %s", fn.Data() ) );
            return false;
        }

        St_Survey *tableSet = (St_Survey *) ds;
        Survey_st *rows = tableSet->GetTable();
        const int nrows  = tableSet->GetNRows();

        std::cout << "\n  " << tables[it] << ": " << nrows << " row(s) from "
                  << gSystem->BaseName( fn.Data() ) << "\n";
        for ( int i = 0; i < nrows; i++ ) {
            printf("    Id=%3d  t=(%9.5f,%9.5f,%9.5f)  r00,r11,r22=(%6.3f,%6.3f,%6.3f)\n",
                   rows[i].Id, rows[i].t0, rows[i].t1, rows[i].t2,
                   rows[i].r00, rows[i].r11, rows[i].r22 );
        }

        if ( writedb ) {
            gSystem->Setenv("DB_ACCESS_MODE","write");
            std::cout << "    DB_ACCESS_MODE=" << gSystem->Getenv("DB_ACCESS_MODE") << "\n";
            StDbManager    *mgr  = StDbManager::Instance();
            StDbConfigNode *node = mgr->initConfig("Geometry_ftt");
            StDbTable      *table = node->addDbTable( tables[it] );
            mgr->setStoreTime( storeTime.Data() );
            table->SetTable( (char*)rows, nrows );
            // storeDbTable returns false when the write is refused -- e.g.
            // "INSERT command denied to user ...". Without this check the macro
            // printed "stored ..." and "Done with database upload" after a write
            // that never happened.
            const bool ok = mgr->storeDbTable( table );
            gSystem->Unsetenv("DB_ACCESS_MODE");
            if ( !ok ) {
                std::cout << "\n*** WRITE FAILED for " << tables[it] << " at "
                          << storeTime << " -- NOTHING WAS STORED ***\n"
                          << "    see the MysqlDb/StDbSql errors above\n";
                gROOT->ProcessLine( Form( ".U %s", fn.Data() ) );
                return false;
            }
            std::cout << "    stored " << tables[it] << " at " << storeTime << "\n";
        }

        gROOT->ProcessLine( Form( ".U %s", fn.Data() ) );
    }

    if ( writedb ) std::cout << "\nDone with database upload.\n";
    else           std::cout << "\nDry run only -- pass \"writedb\" to store.\n";
    return true;
}
