//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable

void daq_track(    int n = 2000,
                    const char *inFile = "/star-sw/FwdWork/star-sw-fwd/st_fwd_25127021_raw_0100003.daq",
                    // const char *inFile = "/gpfs01/star/pwg_tasks/FwdCalib/wangzhen/daq/run23/MultCheck/st_physics_24157073_raw_0500028.daq",
                    std::string configFile = "daq/daq_track.xml",
                    // std::string configFile = "daq/ftt_track.xml",
                    const char *geom = "y2023") {
    TString _chain;
    gSystem->Load( "libStarRoot.so" );

    // Simplest chain with fst, fcs, ftt and fwdTracker
    //bfc( 100, "in, y2023, useXgeom, AgML, db, StEvent, MuDST, evout, cmudst, tree, epdHit, trgd", "/gpfs01/star/pwg_tasks/FwdCalib/DAQ/zfa/st_physics_23072017_raw_2500008.daq" )
    // _chain = Form("in, %s, useXgeom, AgML, db, StEvent, MakeEvent, MuDST, trgd, btof, ftt, fttQA, fstMuRawHit, fwdTrack, evout, cmudst, tree", geom);
    _chain = Form("in, %s, useXgeom, AgML, db, StEvent, MakeEvent, MuDST, trgd, btof, fcs, fst, ftt, fttQA, fstMuRawHit, fwdTrack, evout, cmudst, tree", geom);
    // _chain = Form("in, db, StEvent, MakeEvent, MuDST, trgd, fcs, fst, ftt, fttQA, fstMuRawHit, fwdTrack, evout, cmudst, tree");
    
    // needed in this wonky spack environment 
    gROOT->SetMacroPath(".:/star-sw/StRoot/macros:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");

    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile );

    //  Extra configuration  for the Forward Tracking
    StFwdTrackMaker *fwdTrack = (StFwdTrackMaker*) chain->GetMaker("fwdTrack");
    if ( fwdTrack ){ //if it is in the chain
        fwdTrack->SetConfigFile( configFile );
        // write debug histograms and ttree?
        fwdTrack->SetGenerateTree( true );
        fwdTrack->SetGenerateHistograms( true );
        // write out wavefront OBJ files
        // fwdTrack->SetVisualize( true );
        fwdTrack->SetVisualize( false );
        fwdTrack->SetDebug();
    }

    // StFttDb *fttDb = (StFttDb*)chain->GetMaker("fttDb");
    // fttDb->SetDebug();
    StFttRawHitMaker *fttRawHit = (StFttRawHitMaker*)chain->GetMaker("stgc");
    fttRawHit->SetDebug();
    StFttClusterMaker *fttClu = (StFttClusterMaker*) chain->GetMaker("stgcCluster");
    // fttClu->SetDebug();
    StFttPointMaker *fttPoint = (StFttPointMaker*) chain->GetMaker("stgcPoint");
    // fttPoint->SetDebug();


    StFttHitCalibMaker *fttCalib = (StFttHitCalibMaker*) chain->GetMaker("fttHitCalib");
    fttCalib->SetMode( 0 );

    gSystem->Load("StFcsTrackMatchMaker");
    StFcsTrackMatchMaker *match = new StFcsTrackMatchMaker();
    match->setMaxDistance(6,10);
    match->setFileName("fcstrk.root");
    match->SetDebug();
    chain->AddAfter( "fwdTrack", match);

    StMuDstMaker * muDstMaker = (StMuDstMaker*)chain->GetMaker( "MuDst" );
    chain->AddAfter( "FcsTrkMatch", muDstMaker );


    gSystem->Load("StFwdUtils.so");
    StFwdAnalysisMaker * fwdAna = new StFwdAnalysisMaker();
    chain->AddAfter("MuDst", fwdAna);


    St_db_Maker *db_maker = (St_db_Maker *)chain->GetMaker("db");
    if(db_maker){
        db_maker->SetAttr("blacklist", "tpc");
        db_maker->SetAttr("blacklist", "svt");
        db_maker->SetAttr("blacklist", "ssd");
        db_maker->SetAttr("blacklist", "ist");
        db_maker->SetAttr("blacklist", "pxl");
        db_maker->SetAttr("blacklist", "pp2pp");
        db_maker->SetAttr("blacklist", "ftpc");
        db_maker->SetAttr("blacklist", "emc");
        db_maker->SetAttr("blacklist", "eemc");
        db_maker->SetAttr("blacklist", "mtd");
        db_maker->SetAttr("blacklist", "pmd");
        // db_maker->SetAttr("blacklist", "tof");
    }
    

    // Initialize the chain
    chain->Init();

    // if (fwdTrack){
    //     fwdTrack->test();
    //     return;
    // }
    //_____________________________________________________________________________
    //
    // MAIN EVENT LOOP
    //_____________________________________________________________________________
    for (int i = 0; i < n; i++) {
        chain->Clear();

        if (kStOK != chain->Make())
            break;
    }
}
