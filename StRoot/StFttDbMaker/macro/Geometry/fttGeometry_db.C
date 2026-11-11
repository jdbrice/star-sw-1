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

// Is this slot already in the DB, and does it match the file?
//
//   0  nothing stored at this beginTime     -> safe to store
//   1  complete and identical to the file   -> skip, nothing to do
//   2  complete but the content DIFFERS     -> refuse; we cannot overwrite
//   3  PARTIAL, n rows where N are expected -> refuse loudly, needs a DBA
//  -1  could not ask                        -> refuse, never store blind
//
// Asked with a direct SQL query, deliberately NOT through St_db_Maker. St_db_Maker
// caches per timestamp and served the PRE-write state once already (the data-window
// push, 2026-10-16); a cached "nothing there" would turn this guard into exactly the
// blind insert it exists to prevent.
//
// State 3 is the one that matters. We hold INSERT but not DELETE on Geometry_ftt,
// while StDbSql rolls a partially failed multi-row write back with a DELETE -- so a
// 16-row write that dies at row 7 leaves orphans only the DBAs can remove, and a
// retry would then be a partial insert on top of a partial insert.
//
// Host is the write master from conf/dbLoadBalancerLocalConfig_BNL.xml.
// NOTE the plain double array rather than Survey_st*: CINT parses this signature
// when the macro is loaded, which is before the function body has done
// gSystem->Load("libStDb_Tables.so"), so a Survey_st argument here is "Unknown type
// ... handled as int" and the call then fails to resolve at all.
int fttGeoSlotState( const char *table, const char *beginTime,
                     double *want, int nrows )
{
    TSQLServer *sdb = TSQLServer::Connect( "mysql://robinson.star.bnl.gov:3306/Geometry_ftt",
                                           gSystem->Getenv("USER"), "" );
    if ( !sdb ) { std::cout << "    cannot reach the DB to check\n"; return -1; }

    TString q = Form( "select elementID,Id,r00,r01,r02,r10,r11,r12,r20,r21,r22,t0,t1,t2"
                      " from %s where beginTime='%s' and flavor='ofl' and deactive=0"
                      " order by elementID", table, beginTime );
    TSQLResult *res = sdb->Query( q.Data() );
    if ( !res ) { delete sdb; std::cout << "    query failed\n"; return -1; }

    const int ngot = res->GetRowCount();
    int state = 0;

    if ( ngot == 0 ) {
        state = 0;
    } else if ( ngot != nrows ) {
        std::cout << "    *** " << ngot << " row(s) present where " << nrows
                  << " are expected ***\n";
        state = 3;
    } else {
        int ndiff = 0;
        for ( int jr = 0; jr < nrows; jr++ ) {
            TSQLRow *rw = res->Next();
            if ( !rw ) { ndiff++; break; }
            for ( int jf = 0; jf < 14; jf++ ) {
                const double got = atof( rw->GetField(jf) );
                const double exp = want[ jr * 14 + jf ];
                if ( fabs( got - exp ) > 1e-9 ) {
                    if ( ndiff < 5 )
                        printf("    DIFF row %d field %d: db=%.10f file=%.10f\n",
                               jr, jf, got, exp);
                    ndiff++;
                }
            }
            delete rw;
        }
        state = ndiff ? 2 : 1;
    }

    delete res;
    delete sdb;
    return state;
}

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

    // Does the table carry an element index? stgcOnTpc does not (one row,
    // elementID 0); stationOnStgc and pentOnStation do, and their element lists
    // live in Geometry_ftt.StationIDs (1-4) and PentIDs (1-16).
    //
    // This MATTERS and is easy to get wrong. In Geometry_fst the Nodes row holds
    // an elementID RANGE ("1-36"), so StDbSql::readTableInfo fills the table's
    // element list for you. Geometry_ftt uses the other supported style --
    // elementID='None' plus indexName -- which populates only the element NAME.
    // The list is then left null, and StDbTable::StreamAccessor quietly fills it
    // with 0..n-1. Those do not match the *IDs tables, so the rows would be
    // stored one off: elementID 0 invisible to every read, 1..3 served as
    // Station 1..3, Station 4 missing. We hold INSERT but not DELETE, so that
    // could not be cleaned up afterwards. Pass the list explicitly instead.
    const int indexed[3] = { 0, 1, 1 };
    TString stamp[3];     // per-table file stamp
    TString storeTime, what;

    TString t(tag);
    if ( t == "run22ideal" ) {
        stamp[0]="20211010.000000"; stamp[1]="20211010.000000"; stamp[2]="20211010.000000";
        storeTime = "2021-10-10 00:00:00";
        what = "Run 22 ideal: identity, no shift, no rotation";
    } else if ( t == "run22survey" ) {
        stamp[0]="20211011.000000"; stamp[1]="20211011.000000"; stamp[2]="20211011.000000";
        storeTime = "2021-10-11 00:00:00";
        what = "Run 22 SURVEY: the built detector, no alignment folded in";
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
    } else if ( t == "run22iter4" ) {
        stamp[0]="20211015.000001"; stamp[1]="20211015.000001"; stamp[2]="20211015.000004";
        storeTime = "2021-10-15 00:00:04";
        what = "Run 22 iter 3 MINUS the 17.55 mm strip-frame constant (FINAL Run 22)";
    } else if ( t == "run23ideal" ) {
        stamp[0]="20230410.000000"; stamp[1]="20230410.000000"; stamp[2]="20230410.000000";
        storeTime = "2023-04-10 00:00:00";
        what = "Run 23 ideal (link to the Run 22 ideal; was Daniel's 20221220 placeholder)";
    } else if ( t == "run23est" ) {
        stamp[0]="20230415.000001"; stamp[1]="20230415.000001"; stamp[2]="20230415.000001";
        storeTime = "2023-04-15 00:00:01";
        what = "Run 23 estimate: Run 22 iter 3 carried over (best guess, no Run 23 measurement)";
    } else if ( t == "run23est2" ) {
        stamp[0]="20230415.000001"; stamp[1]="20230415.000001"; stamp[2]="20230415.000002";
        storeTime = "2023-04-15 00:00:02";
        what = "Run 23 estimate, 17.55 mm corrected (Run 22 iter4 carried over)";
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
    } else if ( t == "run24iter2" ) {
        stamp[0]="20231215.000001"; stamp[1]="20231215.000001"; stamp[2]="20231215.000003";
        storeTime = "2023-12-15 00:00:03";
        what = "Run 24 iter 1 MINUS the 17.55 mm strip-frame constant (FINAL Run 24)";
    } else if ( t == "run25ideal" ) {
        stamp[0]="20241210.000000"; stamp[1]="20241210.000000"; stamp[2]="20241210.000000";
        storeTime = "2024-12-10 00:00:00";
        what = "Run 25 ideal (link to the Run 22 ideal)";
    } else if ( t == "run25est" ) {
        stamp[0]="20241215.000001"; stamp[1]="20241215.000001"; stamp[2]="20241215.000001";
        storeTime = "2024-12-15 00:00:01";
        what = "Run 25 estimate: Run 24 iter 1 carried over (best guess)";
    } else if ( t == "run25est2" ) {
        stamp[0]="20241215.000001"; stamp[1]="20241215.000001"; stamp[2]="20241215.000002";
        storeTime = "2024-12-15 00:00:02";
        what = "Run 25 estimate, 17.55 mm corrected (Run 24 iter2 carried over)";
    } else {
        std::cout << "fttGeometry_db: unknown tag \"" << tag << "\"\n"
                  << "  known tags:\n"
                  << "    run22ideal   2021-10-10 00:00:00   identity\n"
                  << "    run22survey  2021-10-11 00:00:00   survey, no alignment\n"
                  << "    run22iter1   2021-10-15 00:00:01   per-quadrant x/y\n"
                  << "    run22iter2   2021-10-15 00:00:02   + per-quadrant z\n"
                  << "    run22iter3   2021-10-15 00:00:03   + deltas\n"
                  << "    run22iter4   2021-10-15 00:00:04   -17.55 mm  (BEST for Run 22)\n"
                  << "    run23ideal   2023-04-10 00:00:00   identity (link)\n"
                  << "    run23est     2023-04-15 00:00:01   Run 22 iter 3 carried over\n"
                  << "    run23est2    2023-04-15 00:00:02   -17.55 mm  (BEST for Run 23)\n"
                  << "    run24ideal   2023-12-10 00:00:00   identity (link)\n"
                  << "    run24est     2023-12-15 00:00:01   Run 22 iter 3 carried over\n"
                  << "    run24iter1   2023-12-15 00:00:02   first Run 24 measurement\n"
                  << "    run24iter2   2023-12-15 00:00:03   -17.55 mm  (BEST for Run 24)\n"
                  << "    run25ideal   2024-12-10 00:00:00   identity (link)\n"
                  << "    run25est     2024-12-15 00:00:01   Run 24 iter 1 carried over\n"
                  << "    run25est2    2024-12-15 00:00:02   -17.55 mm  (BEST for Run 25)\n"
                  << "\n  The 00:00:0x lane is the baseline kept in this directory; the\n"
                  << "  The next day, <year>-<mm>-16 00:00:0x, is the alignment program's lane.\n"
                  << "  Reserved and empty: 00:00:00 on each estimate date, against the\n"
                  << "  day there is a real sTGC survey.\n"
                  << "  Add ideal and estimate slots only in PAIRS -- an ideal entry in a\n"
                  << "  year with nothing later HIJACKS that year to identity.\n"
                  << "\n  Pass \"all\" to push every slot in order.\n";
        if ( t != "all" ) return false;
    }

    if ( t == "all" ) {
        const char *every[16] = { "run22ideal","run22survey",
                                  "run22iter1","run22iter2","run22iter3","run22iter4",
                                  "run23ideal","run23est","run23est2",
                                  "run24ideal","run24est","run24iter1","run24iter2",
                                  "run25ideal","run25est","run25est2" };
        for ( int i = 0; i < 16; i++ ) {
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

        // The element id of row i is that row's own Id. Every file in ftt/ numbers
        // its rows 1..N in the order AGML's <Misalign row="N"/> expects, which is
        // exactly how StationIDs and PentIDs are numbered, so Id == elementID and
        // the stored rows can be checked against the file by eye.
        int elems[128];
        for ( int ie = 0; ie < nrows && ie < 128; ie++ ) elems[ie] = (int) rows[ie].Id;

        std::cout << "\n  " << tables[it] << ": " << nrows << " row(s) from "
                  << gSystem->BaseName( fn.Data() ) << "\n";
        if ( indexed[it] )
            std::cout << "    indexed: elementID taken from each row's Id, "
                      << elems[0] << ".." << elems[nrows-1] << "\n";
        else
            std::cout << "    not indexed: elementID 0\n";
        for ( int i = 0; i < nrows; i++ ) {
            printf("    Id=%3d  t=(%9.5f,%9.5f,%9.5f)  r00,r11,r22=(%6.3f,%6.3f,%6.3f)\n",
                   rows[i].Id, rows[i].t0, rows[i].t1, rows[i].t2,
                   rows[i].r00, rows[i].r11, rows[i].r22 );
        }

        // Flatten the file's rows into the 14 columns the check reads back, in the
        // order the query selects them: elementID, Id, r00..r22, t0..t2.
        double want[16 * 14];
        for ( int jw = 0; jw < nrows && jw < 16; jw++ ) {
            double *w = &want[jw * 14];
            w[0]  = indexed[it] ? (double) rows[jw].Id : 0.0;
            w[1]  = (double) rows[jw].Id;
            w[2]  = rows[jw].r00; w[3]  = rows[jw].r01; w[4]  = rows[jw].r02;
            w[5]  = rows[jw].r10; w[6]  = rows[jw].r11; w[7]  = rows[jw].r12;
            w[8]  = rows[jw].r20; w[9]  = rows[jw].r21; w[10] = rows[jw].r22;
            w[11] = rows[jw].t0;  w[12] = rows[jw].t1;  w[13] = rows[jw].t2;
        }
        const int state = fttGeoSlotState( tables[it], storeTime.Data(), want, nrows );
        if ( state == 1 ) {
            std::cout << "    ALREADY IN THE DB and identical -- nothing to do\n";
            gROOT->ProcessLine( Form( ".U %s", fn.Data() ) );
            continue;
        }
        if ( state != 0 ) {
            std::cout << "\n*** REFUSING to store " << tables[it] << " at " << storeTime
                      << " ***\n";
            if ( state == 2 )
                std::cout << "    an entry is already there and its content DIFFERS.\n"
                          << "    We hold INSERT but not DELETE, so it cannot be replaced:\n"
                          << "    give the new numbers a later beginTime, or ask the DBAs.\n";
            else if ( state == 3 )
                std::cout << "    the entry is PARTIAL. A previous write stopped half way\n"
                          << "    and could not roll itself back. Do not retry -- ask the\n"
                          << "    DBAs to remove the orphan rows first.\n";
            else
                std::cout << "    could not read the current state, so storing would be blind.\n";
            gROOT->ProcessLine( Form( ".U %s", fn.Data() ) );
            return false;
        }

        if ( writedb ) {
            gSystem->Setenv("DB_ACCESS_MODE","write");
            std::cout << "    DB_ACCESS_MODE=" << gSystem->Getenv("DB_ACCESS_MODE") << "\n";
            StDbManager    *mgr  = StDbManager::Instance();
            StDbConfigNode *node = mgr->initConfig("Geometry_ftt");
            StDbTable      *table = node->addDbTable( tables[it] );
            mgr->setStoreTime( storeTime.Data() );
            // stgcOnTpc keeps the no-idList form: it stored elementID 0 on
            // 2026-11-06 and must stay consistent with that row.
            if ( indexed[it] ) table->SetTable( (char*)rows, nrows, elems );
            else               table->SetTable( (char*)rows, nrows );
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
    else           std::cout << "\nDry run only -- pass \"writedb\" to store.\n"
                             << "(the \"already there?\" check above ran for real)\n";
    return true;
}
