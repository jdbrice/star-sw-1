// fstGeometry_db.C -- push the Geometry/fst survey tables to the offline MySQL DB.
//
//   root4star -b -q 'fstGeometry_db.C("survey22")'            dry run: read + print
//   root4star -b -q 'fstGeometry_db.C("survey22","writedb")'  store one slot
//   root4star -b -q 'fstGeometry_db.C("all","writedb")'       store every slot
//   root4star -b -q 'fstGeometry_db.C("bogus")'              list the tags
//
// Run the dry form first. It prints every row it would store, and touches
// nothing.
//
// The numbers are NOT duplicated here. They are read from fst/*.C, the same
// files St_db_Maker itself reads when they are placed in ./StarDb/Geometry/fst,
// so there is exactly one copy of every number and no way for the DB and the
// local-override path to drift apart. See README for what each file means.
//
// WHAT THESE TABLES CONTAIN: not alignment iterations. Every non-identity number
// in fst/*.C is THE SURVEY -- $STAR/StarDb/Geometry/fst ships it as hssOnFst at
// 20211110.000001. Our 20211015.000000 is byte-identical to it, and our
// 20211015.000001 is the same survey re-decomposed. FST alignment has never been
// iterated past the survey, which is why these sit in Gene's SURVEY slot. The
// later run years are symlinks to the same two files. See README.
//
// beginTime follows Gene Van Buren's rule (database_timestamp.txt). The date is
// per RUN YEAR, not per calendar date of the work, and each year reserves two
// slots before the run starts:
//
//     run year   ideal        survey/estimate   real data
//     2021       2020-12-10   2020-12-15        any time after
//     2022       2021-10-10   2021-10-15        any time after
//     2023       2023-04-10   2023-04-15        any time after
//     2024       2023-12-10   2023-12-15 (*)    any time after
//     2025       2024-12-10   2024-12-15        any time after
//
//   (*) database_timestamp.txt as received reads "2024-12-15" for the 2024
//       survey slot, which cannot be right: the ideal slot for the same year is
//       2023-12-10 and Run 24 began in early 2024, so a 2024-12-15 estimate
//       would land after the run it describes. Read as 2023-12-15. CONFIRM WITH
//       GENE before pushing any Run 24 table.
//
// The first entry at a given date uses time 00:00:00; later revisions of the
// same numbers step the time, 00:00:01, 00:01:00, and so on.
//
// All four FST tables are St_Survey. They are stored together, because
// StFstDb::loadGeometry() chains all four and a partial update leaves the
// geometry inconsistent.

void fstGeometry_db( const char* tag = "", const char* opt = "" )
{
    gROOT->Macro("LoadLogger.C");
    gSystem->Load("St_base.so");
    gSystem->Load("libStDb_Tables.so");
    gSystem->Load("StDbLib.so");

    TString option(opt);
    const int writedb = option.Contains("writedb");

    // ---- tag -> (which files, what beginTime) --------------------------------
    TString t(tag);
    TString stamp;      // the timestamp part of the fst/*.C filenames to read
    TString storeTime;  // the beginTime to store under, per the rule above
    TString what;

    // tag -> (file stamp, beginTime, description). One row per DB slot.
    // The later years reuse the Run 22 files: fst/ has symlinks at those stamps,
    // so the same survey is served without a second copy of the numbers.
    const int NTAG = 9;
    const char *tg[NTAG] = { "ideal22", "survey22flat", "survey22",
                             "ideal23", "survey23",
                             "ideal24", "survey24",
                             "ideal25", "survey25" };
    const char *st[NTAG] = { "20211010.000000", "20211015.000000", "20211015.000001",
                             "20230410.000000", "20230415.000000",
                             "20231210.000000", "20231215.000000",
                             "20241210.000000", "20241215.000000" };
    const char *bt[NTAG] = { "2021-10-10 00:00:00", "2021-10-15 00:00:00", "2021-10-15 00:00:01",
                             "2023-04-10 00:00:00", "2023-04-15 00:00:00",
                             "2023-12-10 00:00:00", "2023-12-15 00:00:00",
                             "2024-12-10 00:00:00", "2024-12-15 00:00:00" };
    const char *ds[NTAG] = {
        "Run 22 ideal: identity, misalignment off",
        "Run 22 survey, exactly as the release ships it (all of it in hssOnFst)",
        "Run 22 survey, re-decomposed (common shift in fstOnTpc) -- what we run",
        "Run 23 ideal (link to the Run 22 ideal)",
        "Run 23 survey (link to the Run 22 re-decomposed survey)",
        "Run 24 ideal (link to the Run 22 ideal)",
        "Run 24 survey (link to the Run 22 re-decomposed survey)",
        "Run 25 ideal (link to the Run 22 ideal)",
        "Run 25 survey (link to the Run 22 re-decomposed survey)" };

    int itag = -1;
    for ( int i = 0; i < NTAG; i++ ) if ( t == tg[i] ) itag = i;
    if ( itag < 0 ) {
        std::cout << "fstGeometry_db: unknown tag \"" << tag << "\"\n  known tags:\n";
        for ( int i = 0; i < NTAG; i++ )
            printf("    %-13s %-22s %s\n", tg[i], bt[i], ds[i]);
        std::cout << "\n  survey22 and survey22flat are the SAME alignment decomposed\n"
                  << "  differently -- push ONE of them. survey22 is what we run.\n"
                  << "  Pass \"all\" to push every slot in order.\n";
        if ( t != "all" ) return;
    }

    // "all" loops over every slot; a single tag does just that one.
    int ifirst = ( t == "all" ) ? 0 : itag;
    int ilast  = ( t == "all" ) ? NTAG - 1 : itag;
    // survey22flat is redundant with survey22; skip it in "all".
    for ( int islot = ifirst; islot <= ilast; islot++ ) {
        if ( t == "all" && TString(tg[islot]) == "survey22flat" ) continue;
        stamp     = st[islot];
        storeTime = bt[islot];
        what      = ds[islot];
        fstGeometry_db_one( stamp, storeTime, what, tg[islot], writedb );
    }
    return;
}

// Store (or just print) the four tables at one stamp/beginTime.
void fstGeometry_db_one( TString stamp, TString storeTime, TString what,
                         const char *tagname, int writedb )
{
    std::cout << "\n==== " << tagname << "  files=fst/*." << stamp << ".C\n"
              << "     beginTime=\"" << storeTime << "\"\n"
              << "     " << what << "\n"
              << "     writedb=" << writedb << "\n";

    const char *tables[4] = { "fstOnTpc", "hssOnFst", "fstWedgeOnHss", "fstSensorOnWedge" };

    // Where this macro lives, so it works from any working directory.
    TString here = gSystem->DirName( __FILE__ );

    for ( int it = 0; it < 4; it++ ) {

        TString fn = Form( "%s/fst/%s.%s.C", here.Data(), tables[it], stamp.Data() );
        if ( gSystem->AccessPathName( fn.Data() ) ) {
            std::cout << "ERROR: cannot read " << fn << " -- nothing stored\n";
            return;
        }

        // Each file defines CreateTable(); load one at a time and unload after,
        // or the second file's definition collides with the first.
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
            printf("    Id=%3d  t=(%9.5f,%9.5f,%9.5f)  r00,r11,r22=(%6.3f,%6.3f,%6.3f)\n",
                   rows[i].Id, rows[i].t0, rows[i].t1, rows[i].t2,
                   rows[i].r00, rows[i].r11, rows[i].r22 );
        }

        if ( writedb ) {
            gSystem->Setenv("DB_ACCESS_MODE","write");
            std::cout << "    DB_ACCESS_MODE=" << gSystem->Getenv("DB_ACCESS_MODE") << "\n";
            StDbManager   *mgr  = StDbManager::Instance();
            StDbConfigNode *node = mgr->initConfig("Geometry_fst");
            StDbTable     *table = node->addDbTable( tables[it] );
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
