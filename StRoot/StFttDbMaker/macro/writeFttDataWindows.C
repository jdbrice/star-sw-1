// timestamps:
// before sTGC installed: 2021-10-01 00:00:00

void writeFttDataWindows( bool debug = true, TString storeTime = "2021-10-01 00:00:00" ) {

    // if you want to use root.exe instead of root4star, uncomment block below:

    if ( debug )
        printf( "DEBUG MODE > NOT SAVING TO DB\n" );

   printf( "START-TIME: %s\n\n", storeTime.Data() );;
   gSystem->AddDynamicPath("/usr/lib/mysql");
   gSystem->AddDynamicPath("/usr/lib64/mysql");
   gSystem->AddDynamicPath("$OPTSTAR/lib/mysql/");
   gSystem->Load("libmysqlclient");

    // Load all required libraries
    gROOT->Macro("LoadLogger.C");
    gSystem->Load("St_base.so");
    gSystem->Load("libStDb_Tables.so");
    gSystem->Load("StDbLib.so");

    if ( !debug ){
        printf( "DB_ACCESS_MODE : write !!!!\n\n" );
        gSystem->Setenv("DB_ACCESS_MODE","write");
    }

    // Initialize db manager
    StDbManager* mgr = StDbManager::Instance();
    StDbConfigNode* node = mgr->initConfig("Calibrations_ftt");
    StDbTable* dbtable = node->addDbTable("fttDataWindowsB");
    
    // TString storeTime = "2022-03-06 18:28:40"; // beginTime timestamp in MySQL format: "YYYY-MM-DD HH:mm:ss"
    mgr->setStoreTime(storeTime.Data());

    // Create your c-struct
    fttDataWindowsB_st table;

    /* fob(1-96) x vmm(1-4) = index 1 - 384 */
    for ( int i = 1; i < 385; i++ ){
        table.uuid[i] = i; /* fob(1-96) x vmm(1-4) = index 1 - 384 */
        table.mode[i] = 0; /* 0 = timebin, 1 = bcid */
        table.anchor[i] = 0; /* calibrated time anchor for BCID */
        table.min[i] = -65; /* time window min > -32768 */
        table.max[i] = 100; /* time window max < 32768 */
    }

    
    printf( "TABLE Content:\n\n" );
    for ( int i = 0; i < 385; i++ ){
        printf( "(uuid=%d, mode=%d, anchor=%d, min=%d, max=%d)\n", i, table.mode[i], table.anchor[i], table.min[i], table.max[i] );
    }
    if ( debug ){
        printf( "NOT WRITING, set debug to false to save\n\n" );
        return;
    }


    
    // Fill structure with data 
    // table.feb[0] = 0; // sample setup for a single channel, please add more channels!

    // Store data to the StDbTable
    dbtable->SetTable((char*)&table, 1);

    // uncomment next line to set "sim" flavor. "ofl" flavor is set by default, no need to set it.
    // dbtable->setFlavor("sim");

    // Store table to database
    mgr->storeDbTable(dbtable);
}