//usr/bin/env root4star -l -b -q $0'("'${1:-/gpfs01/star/pwg/mrosales/jetFinderTest2024/star-sw/Jet_Data_NoFilter_500/pythia_jet_vz0_run100.fzd}'",'${2:-100}')'; exit $?
// that is a valid shebang to run script as executable, but with only two arg

// generate some input data using genfzd

TFile *output = 0;

bool RunFttChain = false; // we use GEANT directly
bool RunFstChain = false; // we use GEANT directly
bool RunFcsChain = true;
bool RunFwdChain = true;
bool RunMuDstMaker = true;

// When true, run FTT through the real StFttDb reconstruction chain
// (StFttDbMaker -> StFttSlowSimMaker -> StFttClusterMaker ->
// StFttClusterPointMaker, mUseGeantData=false) in place of reading GEANT
// truth directly (MakeGeantPoints, mUseGeantData=true) -- exercises the
// same geometry/clustering code real data uses.
// Set false for the old/original GEANT-truth-direct behavior.
bool UseFttSimHitMaker = true;

// When true, run FST through the real raw-hit reconstruction chain
//   StFstSlowSimMaker -> StFstRawHitMaker -> StFstClusterMaker -> StFstHitMaker
// instead of StFstFastSimMaker, which is a parallel reimplementation of the
// clustering and hit positioning. The BFC options already exist and
// BigFullChain.h states the recipe: use fstSlowSim WITH fstRawHit,fstCluster,
// fstHit and WITHOUT fstFastSim. StFstSlowSimMaker publishes "fstRawAdcSimu",
// which StFstRawHitMaker already knew how to consume.
//
// Implies RunFstChain (the BFC fst tags have to be in the chain string) and
// switches the tracker's FST input from GEANT truth to StEvent, since the
// point is to read back what the real clustering produced. Read the SCOPE
// section of StFstSlowSimMaker.h first: v0 is one raw hit per GEANT hit, so
// every cluster is single-strip and nothing about cluster size or resolution
// is meaningful yet -- it is a geometry/chain closure tool.
bool UseFstSlowSim = false;

bool UseCachedGeom = true;
bool UseConstBz = false;
bool UseZeroB = false;
bool doFwdResidual = true;  // FST/FTT alignment residuals (StFwdResidualMaker); no MuDst file needed
// Unbiased (hit-removed) FST/FTT residuals (StFwdAlignmentMaker) -- see
// proposal_alignment_path.txt. Off by default: costs ~(num FST+FTT planes)
// extra refits per Global track, not meant for routine running. Completely
// independent of doFwdResidual/StFwdResidualMaker above -- separate maker,
// separate output file, doesn't touch StEvent/MuDst/PicoDst.
bool doFwdAlignment = false;
// StFwdTrack::StFwdTrackType: 0=Global 1=BLC 2=Primary 3=FwdVtx 4=BLCVtx 5=FCSConstrained
// One StFwdResidualMaker instance per entry, all in the same chain/event loop -- the
// tracker already fits every type each event, so there's no need to re-run sim() once
// per type (that used to cost 3x the wall time for no physics reason).
const int kNResidualTypes = 3;
int residualTypesToRun[kNResidualTypes] = {0, 1, 4}; // Global, BLC, BLCVtx
bool fieldOff = false;  // BFC MagF tag: reconstruct with StarMagField factor=0 instead of the hardcoded ReverseField default

TString _fttChain = "fttSim";
TString _fcsChain = "fcsSim fcsWFF fcsCluster";
TString _fstChain = "fstFastSim";
TString _fwdTrackChain = "fwdTrack";
TString _geom = "y2022 agml usexgeom";

void DisableTrackFitting() {
    // Disable track fitting
    StFwdTrackMaker * fwdTrack = (StFwdTrackMaker*) chain->GetMaker( "fwdTrack" );
    assert( fwdTrack );
    fwdTrack->setTrackFittingOff();
}

void DoOnlyGlobalTrackFitting() {
    // Disable track fitting
    StFwdTrackMaker * fwdTrack = (StFwdTrackMaker*) chain->GetMaker( "fwdTrack" );
    if ( fwdTrack ){
        fwdTrack->setConfigKeyValue("TrackFitter:refit", false);
        fwdTrack->setConfigKeyValue("TrackFitter:doGlobalTrackFitting", true);
        fwdTrack->setConfigKeyValue("TrackFitter:doBeamlineTrackFitting", true);
        fwdTrack->setConfigKeyValue("TrackFitter:doPrimaryTrackFitting", true);
        fwdTrack->setConfigKeyValue("TrackFitter:doSecondaryTrackFitting", true);
        // skip finding fwd vertices
        fwdTrack->setConfigKeyValue("TrackFitter:findFwdVertices", true);
    }
}

void sim(Int_t n=1000, Int_t run=1, const char* pid="JPsi", float vz=0.0) {
  //void sim(Int_t n=1000, Int_t run=1, const char* pid="ele", float vz=0.0) {
    // report all of the parameters passed in
    char inFile[200];
    TString spid(pid);
    if (spid.Contains(".")) {
        // particle gun: pid = "ele.pt1.5" → ele.pt1.5.vz0.run1.fzd
        sprintf(inFile, "%s.vz%d.run%i.fzd", pid, (int)vz, run);
    } else {
        // Pythia: pid = "JPsi" → pythia.JPsi.vz0.run1.fzd
        sprintf(inFile, "pythia.%s.vz%d.run%i.fzd", pid, (int)vz, run);
    }
    cout << "inFile = " << inFile << endl;
    cout << "# of Events = " << n << endl;

    // to use the geom cache (skip agml build which is faster)
    // set the _geom string to "" and make sure the cache file ("fGeom.root") is present
    if (UseCachedGeom)
        _geom = "";

    // Setup the chain for reading an FZD
    TString _chain = "";

    // FST slow sim replaces the fast sim with the real raw-hit chain, and needs
    // the fst tags in the chain string -- so it forces RunFstChain on rather
    // than silently doing nothing when the two toggles disagree.
    if (UseFstSlowSim) {
        RunFstChain = true;
        _fstChain = "fstSlowSim fstRawHit fstCluster fstHit";
        cout << "UseFstSlowSim: FST through the real raw-hit chain, _fstChain = "
             << _fstChain << endl;
    }

    // Now turn off parts of the chain that we don't need
    if (!RunFttChain)
        _fttChain = "";
    if (!RunFcsChain)
        _fcsChain = "";
    if (!RunFstChain)
        _fstChain = "";
    if (!RunFwdChain)
        _fwdTrackChain = "";

    if (RunFcsChain && RunFwdChain){
        _fwdTrackChain = "fwdTrack fcsTrackMatch";
    }
    
    // Form the complete chain
    const char* _fieldTag = fieldOff ? "FieldOff" : "ReverseField";
    _chain = Form("fzin %s sdt20211016 %s %s %s %s MakeEvent StEvent McEvent %s bigbig CMuDST tree", _geom.Data(), _fttChain.Data(), _fcsChain.Data(), _fstChain.Data(), _fwdTrackChain.Data(), _fieldTag); 
    // Note, I dont include the PicoWrite and PicoVtxless in chain because they load a bunch of things I dont want (and somehow cannot remove with -options)
    printf("Chain: \n%s\n", _chain.Data());
    
    gSystem->Load( "libStarRoot.so" );
    // Preload StFttDbMaker BEFORE bfc builds the chain. bfc loads
    // libStFwdTrackMaker, which has an unresolved reference to StFttDb
    // statics (e.g. StFttDb::X_shift_QuadC); with RunFttChain=false the
    // _fttChain tag is blank, so bfc never loads StFttDbMaker itself and the
    // dlopen fails with "undefined symbol: _ZN7StFttDb13X_shift_QuadCE",
    // tripping StBFChain's libraryload assertion. The symbol exists in the
    // current library -- this is purely a load-order problem.
    // (St_base/StChain first -- StFttDbMaker inherits StMaker, so loading it
    // bare here fails with "undefined symbol: _ZNK7StMaker7GetNameEv".)
    gSystem->Load( "St_base" );
    gSystem->Load( "StChain" );
    gSystem->Load( "libStFttDbMaker.so" );
    gROOT->LoadMacro("bfc.C");
    bfc(-1, _chain, inFile);

    // The "FieldOff" chain tag alone doesn't reliably zero the field in this particular
    // tag combination (empirically still reports Scale factor=1, real Bz~5kG) -- force it
    // directly instead of relying on tag parsing.
    if (fieldOff && StarMagField::Instance()) {
        StarMagField::Instance()->SetFactor(0.0);
        cout << "fieldOff: forced StarMagField::SetFactor(0.0), factor now = "
             << StarMagField::Instance()->GetFactor() << endl;
    }

    TString outfile(inFile);
    outfile.ReplaceAll(".fzd",".root");
    cout << "output file=" <<outfile<<endl;
    chain->SetOutputFile(outfile);
    StMuDstMaker * muDstMaker = (StMuDstMaker*)chain->GetMaker( "MuDst" );    
    cout << "mudst file=" << muDstMaker->GetFile() << endl;

    if ( UseConstBz )
        StarMagField::setConstBz(true);

    // FCS setup, if included
    if (RunFcsChain) {

        StFcsDbMaker* fcsdbmkr = (StFcsDbMaker*) chain->GetMaker("fcsDbMkr");
        cout << "fcsdbmkr="<<fcsdbmkr<<endl;
        StFcsDb* fcsdb = (StFcsDb*) chain->GetDataSet("fcsDb");
        cout << "fcsdb="<<fcsdb<<endl;
        fcsdb->forceFixGain();
        fcsdb->forceFixGainCorrection();
        fcsdb->setDbAccess(0);
	fcsdb->InitRun(0);
	for(int det=0; det<4; det++){
	  StThreeVectorD off = fcsdb->getDetectorOffset(det);
	  printf("AAA FCSDB Det=%1d OFF = %6.2f %6.2f %6.2f\n",det,off.x(),off.y(),off.z());
	}

        // Configure FCS simulator
        StFcsFastSimulatorMaker *fcssim = (StFcsFastSimulatorMaker*) chain->GetMaker("fcsSim");
        fcssim->setDebug(0); // was 1 -- per-hit ECAL/HCAL/Pres LOG_INFO dump, ~500k lines/run
        //fcssim->setLeakyHcal(0);

        StFcsWaveformFitMaker *fcsWFF= (StFcsWaveformFitMaker*) chain->GetMaker("StFcsWaveformFitMaker");
        fcsWFF->setEnergySelect(0);

        StFcsClusterMaker *fcsclu = (StFcsClusterMaker*) chain->GetMaker("StFcsClusterMaker");
        fcsclu->setDebug(0); // was 1 -- triggers full hit/cluster/point dump (mFcsCollection->print(3)) every event
    }

    gSystem->Load("StFwdUtils.so");
    gSystem->Load("StFwdResidualMaker.so");
    gSystem->Load("StFwdAlignmentMaker.so");

    // Configure FST FastSim -- not present at all when the slow sim is driving
    if (RunFstChain && !UseFstSlowSim){ // otherwise it is not loaded
        StFstFastSimMaker *fstFastSim = (StFstFastSimMaker*) chain->GetMaker( "fstFastSim" );;
        if (fstFastSim) {
            printf("fstFastSim = %p\n", fstFastSim);
            TString qaoutname(gSystem->BaseName(inFile));
            qaoutname.ReplaceAll(".fzd", ".FastSimu.QA.root");
            
            // if (SiIneff)
            //     fstFastSim->SetInEfficiency(0.1); // inefficiency of Si

            fstFastSim->SetQAFileName(qaoutname);
        }
    }

    if ( UseFttSimHitMaker ) {
        // Real chain: StFttDbMaker -> StFttSlowSimMaker -> StFttClusterMaker ->
        // StFttClusterPointMaker (mUseGeantData=false), same makers/order as
        // fwd_afterburner_db.C's real-data chain, with StFttSlowSimMaker in
        // place of StFttRawHitMaker/StFttHitCalibMaker.
        gSystem->Load( "libStFttDbMaker.so" );
        gSystem->Load( "libStFttSimMaker.so" );
        gSystem->Load( "libStFttClusterMaker.so" );
        gSystem->Load( "libStFttClusterPointMaker.so" );

        // STOPGAP. The real fix for the ordering is in the slow-sim PR branches in
        // fcstrk13, which put each slow sim into its own bfc chain so that bfc
        // orders it: fcstrk13 commits 2f3eb26a77 ("put FttSimHit in fttChain,
        // before FttClu") and 1e6c787aab ("put fstSlowSim in fstChain"), review
        // #26. fcstrk11 is behind both -- it has no fttSlowSim entry at all and
        // fstSlowSim with an empty chain field -- so the hand-assembly below stands
        // in until those land. Once they do, delete it and add "fttSlowSim" to
        // _fttChain; only StFttClusterPointMaker will still need adding by hand,
        // since bfc's FttPoint is StFttPointMaker.
        //
        // ALL FOUR go in with AddBefore("fwdTrack"), in order. AddMaker()
        // APPENDS to the end of the chain, and bfc has already put fwdTrack at
        // position 9 of 12 -- so the three makers below used to run AFTER the
        // tracker (indeed after StTreeMaker and StPicoDstMaker; the logs show
        // StFttClusterMaker::Init firing dead last) while only the point maker
        // was inserted before it. The point maker therefore ran on whatever
        // clusters the PREVIOUS event had left behind, and nothing the FTT
        // produced could reach the tracker in the same event.
        //
        // AddBefore inserts immediately before the named maker, so four calls
        // in this order give  db -> slowSim -> cluster -> point -> fwdTrack.
        // There is no bfc option for either StFttSlowSimMaker or
        // StFttClusterPointMaker (fttSim is the FAST simulator and FttPoint is
        // StFttPointMaker), which is why this is assembled by hand.
        // With RunFwdChain=false there is no fwdTrack to insert before, and
        // AddBefore on a name that is not in the chain does not place the maker.
        // Appending is the right answer in that case -- nothing downstream
        // consumes the points -- so decide once and use it for all four.
        bool fttHaveAnchor = ( chain->GetMaker("fwdTrack") != 0 );
        if ( !fttHaveAnchor )
            cout << "UseFttSimHitMaker: no fwdTrack in the chain, appending the FTT makers" << endl;

        StFttDbMaker * fttDbMk = new StFttDbMaker();
        if ( fttHaveAnchor ) chain->AddBefore("fwdTrack", fttDbMk); else chain->AddMaker(fttDbMk);

        StFttSlowSimMaker * fttSlowSim = new StFttSlowSimMaker();
        if ( fttHaveAnchor ) chain->AddBefore("fwdTrack", fttSlowSim); else chain->AddMaker(fttSlowSim);

        StFttClusterMaker * fttClu = new StFttClusterMaker();
        fttClu->SetTimeCut( 1 /*kTimeCutModeAcceptAll*/, -9999, 9999 ); // no real timing in MC
        if ( fttHaveAnchor ) chain->AddBefore("fwdTrack", fttClu); else chain->AddMaker(fttClu);

        StFttClusterPointMaker * fttClusterPointMaker = new StFttClusterPointMaker("fttClusterPointMaker");
        fttClusterPointMaker->SetDebug(0); // was 1 -- extra per-hit debug prints
        // mUseGeantData left at its constructor default (false) -- take the
        // real MakeLocalPoints/MakeGlobalPoints path through StFttDb.
        if ( fttHaveAnchor ) chain->AddBefore("fwdTrack", fttClusterPointMaker); else chain->AddMaker(fttClusterPointMaker);
    } else {
        // Old/original behavior: read GEANT truth directly in global
        // coordinates (StFttClusterPointMaker::MakeGeantPoints()) -- never
        // calls into StFttDb's real geometry transform. Kept for backward
        // compatibility (UseFttSimHitMaker = false).
        gSystem->Load( "StFttDbMaker" );
        gSystem->Load( "libStFttSimMaker" );
        gSystem->Load( "libStFttClusterPointMaker" );
        // make an StFttClusterPointMaker
        StFttClusterPointMaker * fttClusterPointMaker = new StFttClusterPointMaker("fttClusterPointMaker");
        fttClusterPointMaker->SetDebug(0); // was 1 -- extra per-hit debug prints
        fttClusterPointMaker->setUseGeantData( true );
        chain->AddBefore("fwdTrack", fttClusterPointMaker);
    }
        
    // Configure the Forward Tracker
        StFwdTrackMaker * fwdTrack = (StFwdTrackMaker*) chain->GetMaker( "fwdTrack" );

        if ( fwdTrack ){
            if ( _geom == "" ){
                cout << "Using the Geometry cache: fGeom.root" << endl;
                fwdTrack->setGeoCache( "fGeom.root" );
            }

            fwdTrack->setOutputFilename( TString::Format( "%s.output.root", inFile ).Data() );

            // Fitter
            fwdTrack->setFitDebugLvl( 0 );
            fwdTrack->setFitMinIterations( 10 );
            fwdTrack->setFitMaxIterations( 20 );
            
            fwdTrack->setDeltaPval( 1e-1 );
            fwdTrack->setRelChi2Change( 1e-6 );
            
            // StFwdHitLoader::DataSource is { GEANT=0, STEVENT=1, MUDST=2, IGNORE=3 },
            // so the 1 here is STEVENT. It used to be commented "IGNORE", which is
            // the opposite of what it does -- the value was right, the label was not.
            // STEVENT is correct either way round: with UseFttSimHitMaker the point
            // maker fills StEvent from the real clustering, without it the same maker
            // fills StEvent from GEANT truth. Only the contents differ.
            fwdTrack->setFttHitSource( 1 /*StFwdHitLoader::STEVENT*/ );
            // FST: GEANT truth, unless the slow sim is running the real chain, in
            // which case read back the StFstHits that StFstHitMaker produced.
            fwdTrack->setFstHitSource( UseFstSlowSim ? 1 /*STEVENT*/ : 0 /*GEANT*/ );

            // DisableTrackFitting();
            // DoOnlyGlobalTrackFitting();
            // fwdTrack->setTrackFittingOff();
            fwdTrack->setConfigKeyValue( "TrackFitter:refit", true );
            
            if ( UseZeroB ){
                cout << "Setting B = 0" << endl;
                fwdTrack->setZeroB( true );
            }
            if ( UseConstBz ){
                cout << "Setting Bz = const everywhere" << endl;
                fwdTrack->setConstBz( true );
            }

            
            cout << "fwd tracker setup" << endl;
        }
    
    bool doFitQA = false;
    if ( doFitQA ){
        StFwdFitQAMaker *fwdFitQA = new StFwdFitQAMaker();
        fwdFitQA->SetDebug();
        TString fitqaoutname(gSystem->BaseName(inFile));
        fitqaoutname.ReplaceAll(".fzd", ".FwdFitQA.root");
        fwdFitQA->setOutputFilename( fitqaoutname );
        chain->AddAfter("fwdTrack", fwdFitQA);
    }

    bool doFwdAna = false;
    if (!RunFcsChain && doFwdAna ){
        StFwdAnalysisMaker *fwdAna = new StFwdAnalysisMaker();
        fwdAna->SetDebug();
        chain->AddAfter("fwdTrack", fwdAna);
    }

    StFwdResidualMaker *fwdResiduals[kNResidualTypes] = {NULL, NULL, NULL};
    if (doFwdResidual && RunFwdChain){
        const char* residualTypeName[6] = {"Global","BLC","Primary","FwdVtx","BLCVtx","FCSConstrained"};
        for (int i = 0; i < kNResidualTypes; i++) {
            int rt = residualTypesToRun[i];
            TString residualName(gSystem->BaseName(inFile));
            residualName.ReplaceAll(".fzd", Form(".FwdDetResidual_%s.root", residualTypeName[rt]));
            fwdResiduals[i] = new StFwdResidualMaker(residualName, Form("fwdResidual_%s", residualTypeName[rt]));
            fwdResiduals[i]->setTrackType((UChar_t)rt);
            chain->AddAfter("MuDst", fwdResiduals[i]);  // must run after MuDst is fully populated, not just after fwdTrack
        }
    }

    const int kNAlignTypes = 2;
    int alignTypesToRun[kNAlignTypes] = {0, 4}; // Global, BLCVtx
    StFwdAlignmentMaker *fwdAlignments[kNAlignTypes] = {NULL, NULL};
    if (doFwdAlignment && RunFwdChain){
        // Unbiased (hit-removed) residuals -- see proposal_alignment_path.txt.
        // Completely separate output/maker from doFwdResidual above; does not
        // read or write anything it touches. Global (clean baseline) + BLCVtx
        // (to test the FST1/FST2 translation signal found there).
        StFwdTrackMaker *fwdTrackForAlign = (StFwdTrackMaker*) chain->GetMaker("fwdTrack");
        const char* alignTypeName[6] = {"Global","BLC","Primary","FwdVtx","BLCVtx","FCSConstrained"};
        for (int i = 0; i < kNAlignTypes; i++) {
            int rt = alignTypesToRun[i];
            TString alignName(gSystem->BaseName(inFile));
            alignName.ReplaceAll(".fzd", Form(".FwdAlignment_%s.root", alignTypeName[rt]));
            fwdAlignments[i] = new StFwdAlignmentMaker(alignName, Form("fwdAlignment_%s", alignTypeName[rt]));
            fwdAlignments[i]->setTrackMaker(fwdTrackForAlign);
            fwdAlignments[i]->setTrackType((UChar_t)rt);
            chain->AddAfter("fwdTrack", fwdAlignments[i]);  // only needs fwdTrack's results, not MuDst
        }
    }

    //StMuDstMaker * muDstMaker = (StMuDstMaker*)chain->GetMaker( "MuDst" );    
    // if (RunFcsChain) {
    //     // FwdTrack and FcsCluster assciation
    //     gSystem->Load("StFcsTrackMatchMaker");
    //     StFcsTrackMatchMaker *match = new StFcsTrackMatchMaker();
    //     match->setMaxDistance(6,10);
    //     match->setFileName("fcstrk.root");
    //     match->SetDebug();
    //     chain->AddMaker(match);

    //     if ( doFwdAna ){
    //         StFwdAnalysisMaker *fwdAna = new StFwdAnalysisMaker();
    //         fwdAna->SetDebug();
    //         chain->AddAfter("FcsTrkMatch", fwdAna);
    //     }

    //     // Produce MuDst output
    //     if ( muDstMaker )
    //         chain->AddAfter( "FcsTrkMatch", muDstMaker );
    // } else {
    //     if ( muDstMaker && doFwdAna )
    //         chain->AddAfter( "fwdAna", muDstMaker );
    // }

    // The PicoDst
    gSystem->Load("libStPicoEvent");
    gSystem->Load("libStPicoDstMaker");
    StPicoDstMaker *picoMk = new StPicoDstMaker(StPicoDstMaker::IoWrite);
    cout << "picoMk = " << picoMk << endl;
    picoMk->setVtxMode(StPicoDstMaker::Vtxless);

    StMemStat stmem;
    stmem.PrintMem("MEM before Chain::Init");
chain_loop:
    chain->Init();
    stmem.PrintMem("MEM after Chain::Init");
    cout << "mudst file=" << muDstMaker->GetFile() << endl;

    //_____________________________________________________________________________
    //
    // MAIN EVENT LOOP
    //_____________________________________________________________________________
    for (int i = 0; i < n; i++) {

        cout << "--------->START EVENT: " << i << endl;

        if (i > 1)
            stmem.PrintMem("MEM before Chain::Clear + Make");
        chain->Clear();
        if (kStOK != chain->Make())
            break;

        if (i > 1)
            stmem.PrintMem("MEM after Chain::Clear + Make");

        // StMuDst * mds = muDstMaker->muDst();
        // StMuFwdTrackCollection * ftc = mds->muFwdTrackCollection();
        // cout << "Number of StMuFwdTracks: " << ftc->numberOfFwdTracks() << endl;
        // for ( size_t iTrack = 0; iTrack < ftc->numberOfFwdTracks(); iTrack++ ){
        //     StMuFwdTrack * muFwdTrack = ftc->getFwdTrack( iTrack );
        //     cout << "muFwdTrack->mPt = " << muFwdTrack->momentum().Pt() << endl;

        // }
        cout << "<---------- END EVENT" << endl;
    } // event loop

    for (int i = 0; i < kNResidualTypes; i++) {
        if (fwdResiduals[i]) fwdResiduals[i]->Finish();
    }
    for (int i = 0; i < kNAlignTypes; i++) {
        if (fwdAlignments[i]) fwdAlignments[i]->Finish();
    }

    stmem.PrintMem("MEM after event loop");
    // delete chain;
    stmem.PrintMem("MEM after delete chain");
}
