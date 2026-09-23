//usr/bin/env root4star -l root -l -q  $0; exit $?
//usr/bin/env root4star -l -b -q $0'("'${1:-st_physics_23055058_raw_1500001.MuDst.root}'",'${2:-100}')'; exit $?
// that is a valid shebang to run script as executable, but with only one argd

#include <typeinfo.h>

// Fast fwd tracking without DB
// bool runDb = false;
// bool runFttChain = true;
// bool runFcsChain = false; 
// bool runFwdChain = true;
// bool refillMuDst = false;
// bool runFwdQa = false;
// bool runFitQa = false;
// bool runPico = true;

// For EPD QA only
// bool runDb = false;
// bool runFttChain = false;
// bool runFcsChain = true;
// bool runFwdChain = false;
// bool refillMuDst = false;
// bool runFwdQa = false;
// bool runFitQa = true;

// Tracking without FCS (but with DB)
// bool runDb = true;
// bool runFttChain = true;
// bool runFcsChain = true; 
// bool runFwdChain = true;
// bool refillMuDst = false;
// bool runFwdQa = false;
// bool runFitQa = false;
// bool runPico = true;


// Memory Baseline
bool runDb       = true;
bool runFttChain = true;
bool runFcsChain = true;
bool runFwdChain = true;
bool refillMuDst = false;
bool runFwdQa    = true;
bool runFitQa    = true;
bool runPico     = true;
bool runFcsQa    = true;
#include "StMemStat.h"


void loadLibs();
void fwd_afterburner_mod( 	const Char_t * fileList = "st_physics_23037002_raw_1000064.MuDst.root",
						size_t nEvents = 1400,
						size_t nSkip   = 0 ){
  cout << "FileList: " << fileList << endl;
  cout << "nEvents: "  << nEvents  << endl;
  cout << "nSkip: "    << nSkip    << endl;

  // First load some shared libraries we need
  loadLibs();

  //Uncomment to reduce message printing
  //gMessMgr->SetLimit("I", 0);
  //gMessMgr->SetLimit("Q", 0);
  //gMessMgr->SetLimit("W", 0);

  // create the chain
  StChain *chain  = new StChain("StChain");

  const char* inMuDstFile = fileList;
  // create the StMuDstMaker
  StMuDstMaker *muDstMaker = new StMuDstMaker(  	0,
							0,
							"",
							inMuDstFile,
							"MuDst.root",
							1
							);
  TChain& muDstChain = *muDstMaker->chain();
  printf( "MuDst file has %d events available in tree\n", muDstChain.GetEntries());
  muDstChain.SetCacheSize(0);   // disable input TTreeCache: avoids ~128 MB-class per-event VSize jumps
  printf( "Input TTreeCache size set to 0 (disabled)\n");
	
  /*******************************************************************************************/
  // Initialize the database
  muDstChain.GetEntry(0);
  int run = muDstMaker->muDst()->event()->runNumber();
  if (runDb){
    cout << endl << "============  Data Base =========" << endl;
    St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
    time_t tt = (time_t)muDstMaker->muDst()->event()->eventInfo().time();
    printf("Run number from 1st event: %d time: %d\n", run, tt);
    struct tm* gmt = gmtime(&tt);
    int date = (gmt->tm_year + 1900) * 10000 + (gmt->tm_mon + 1) * 100 + gmt->tm_mday;
    int itime = gmt->tm_hour * 10000 + gmt->tm_min * 100 + gmt->tm_sec;
    printf("GMT date: %d time: %d\nhour:%d|min:%d|sec:%d\n", date, itime,gmt->tm_hour,gmt->tm_min,gmt->tm_sec);
    dbMk->SetDateTime(date, itime);
    // things will run fine without a timestamp set, but FCS DB will give bad values ...
  }
  /*******************************************************************************************/
	

  /*******************************************************************************************/
  // Create the StMuDst2StEventMaker
  StMuDst2StEventMaker * mu2ev = new StMuDst2StEventMaker();
  mu2ev->SetActive(true);
  /*******************************************************************************************/

  /*******************************************************************************************/
  // Setup Fcs Database if needed
  if ( (runFcsChain && runDb) || runFitQa){
    StFcsDbMaker * fcsDbMkr = new StFcsDbMaker();
    //chain->AddMaker(fcsDb);
    // fcsDb->SetDebug();
  }
  /*******************************************************************************************/
	

  /*******************************************************************************************/
  // FTT chain
  if (runFttChain){
    StFttDbMaker * fttDbMk = new StFttDbMaker();
    //chain->AddMaker(fttDbMk);
    StFttRawHitMaker* fttrawhitmkr = new StFttRawHitMaker();
    fttrawhitmkr->setReadMuDst(1);
    StFttHitCalibMaker * ftthcm = new StFttHitCalibMaker();
    StFttClusterMaker * fttclu = new StFttClusterMaker();
    fttclu->SetTimeCut(1, -40, 40);
    StFttClusterPointMaker *fttCP = new StFttClusterPointMaker();
    //StFttPointMaker * fttpoint = new StFttPointMaker();
  }
  /*******************************************************************************************/

  /*******************************************************************************************/
  // FCS Chain
  if (runFcsChain){
    //StFcsWaveformFitMaker *fcsWFF = new StFcsWaveformFitMaker();  //Turn on if need to reanalyze waveforms or reapply gains
    //fcsWFF->setEnergySelect(0);  // This should only be used for simulated data, for real data this done in the database
    //fcsWFF->setAnaWaveform(false); // This skips waveform analysis, and only apply new gain from DB
    StFcsRawHitMaker* hit = new StFcsRawHitMaker();
    hit->setReadMuDst(1); //assuming reading from MuDst, so always true
    //StFcsClusterMaker *fcsclu = new StFcsClusterMaker();  //Turn on re-cluster which may be needed after a new gain is applied
  }
  /*******************************************************************************************/

  /*******************************************************************************************/
  // FwdTrackMaker Chain
  StFwdTrackMaker *fwdTrack = NULL;
  if (runFwdChain){
    // FwdTrackMaker
    fwdTrack = new StFwdTrackMaker();
    fwdTrack->SetDebug(1);
    fwdTrack->setGeoCache( "fGeom.root" );
    fwdTrack->setSeedFindingWithFst();
    fwdTrack->setTrackRefit( true );
    fwdTrack->setFillAlignment( true );
    fwdTrack->setAlignmentOutputFilename( "align_test.root" );

    // Fitter Options
    fwdTrack->setFitDebugLvl( 0 );
    fwdTrack->setFitMinIterations( 40 );
    fwdTrack->setFitMaxIterations( 100 );
		
    // fwdTrack->setDeltaPval( 1e-9 );
    // fwdTrack->setRelChi2Change( 1e-9 );

    // fwdTrack->setSeedFindingOff();
    // fwdTrack->setTrackFittingOff();
    fwdTrack->setFstHitSource( 2 /* = MUDST */);
    fwdTrack->setFttHitSource( 1 /* = STEVENT */);


    // fwdTrack->setConfigKeyValue("TrackFitter:doBeamlineTrackFitting", false);
    // fwdTrack->setConfigKeyValue("TrackFitter:doPrimaryTrackFitting", false);
    // fwdTrack->setConfigKeyValue("TrackFitter:doSecondaryTrackFitting", false);
    // skip finding fwd vertices
  }



  if (runFwdChain && runFcsChain){
    // FwdTrack and FcsCluster assciation
    gSystem->Load("StFcsTrackMatchMaker");
    StFcsTrackMatchMaker *match = new StFcsTrackMatchMaker();
    match->setMaxDistance(6,10);
    match->setFileName("fcstrk.root");
  }

  
  if (runFwdQa){
    StFwdQAMaker *fwdQA = new StFwdQAMaker();
    fwdQA->SetDebug(2);
    TString fwdqaname( gSystem->BaseName(inMuDstFile) );
    fwdqaname.ReplaceAll(".MuDst.root", "_FwdHists.root");
    //cout << fwdqaname.Data() << endl;
    //fwdQA->setTreeFilename(fwdqaname);
    fwdQA->setLocalOutputFile(fwdqaname.Data());
    gSystem->Load("StFwdUtils.so");
    StFwdAnalysisMaker * fwdAna = new StFwdAnalysisMaker();
    fwdAna->setMuDstInput();
  }


  // The PicoDst
  if (runPico){
    gSystem->Load("libStPicoEvent");
    gSystem->Load("libStPicoDstMaker");
    StPicoDstMaker *picoMk = (StMaker*) (new StPicoDstMaker(StPicoDstMaker::IoWrite, inMuDstFile, "picoDst"));
    cout << "picoMk = " << picoMk << endl;
    picoMk->setVtxMode(StPicoDstMaker::Vtxless);
  }

  if ( runFitQa && runFwdChain){
    StFwdFitQAMaker *fwdFitQA = new StFwdFitQAMaker();
    fwdFitQA->SetDebug();
    TString fitqaoutname(gSystem->BaseName(inMuDstFile));
    fitqaoutname.ReplaceAll(".MuDst.root", ".FwdFitQA.root");
    fwdFitQA->setOutputFilename( fitqaoutname );
  }



  if(runFcsChain && runFcsQa){
    StEpdDbMaker* epddb = new StEpdDbMaker();
    StEpdHitMaker* epdhitmkr = new StEpdHitMaker();
    epdhitmkr->setReadMuDst();
    gSystem->Load("StFwdData");
    TString foriternum(inMuDstFile);
    Ssiz_t last_ = foriternum.Last('_');
    TString filename = "StFcsRun22Qa_";
    TString iternum = foriternum(last_+1,7);
    iternum.ReplaceAll(".list","");
    TString runnum;
    if( foriternum.Contains("_raw_") ){ runnum = foriternum(last_-12,12); }
    else{ runnum = foriternum(last_-8,8); }
    //TString filenametree = "FcsRun22Pi0Ana_";
    filename += runnum + "_" + iternum + ".root";
    HistManager* treehists = new HistManager();
    treehists->InitFile(filename.Data(),"RECREATE");

    gSystem->Load("StFwdAna");
    StFwdAnaData* fwdanadata = new StFwdAnaData();
    fwdanadata->setTreeOnBit(0);
    fwdanadata->setRandomSeed(time(0));
    fwdanadata->setEpdNmipCut(0.7);
    //fwdanadata->setIgnoreTrig();
    //Below is list of all triggers for FCS Run 22
    fwdanadata->AddTrig("fcsJPsi");
    fwdanadata->AddTrig("fcsJPDE1");
    fwdanadata->AddTrig("fcsJPDE0");
    fwdanadata->AddTrig("fcsJPBC1");
    fwdanadata->AddTrig("fcsJPBC0");
    fwdanadata->AddTrig("fcsJPA1");
    fwdanadata->AddTrig("fcsJPA0");
    fwdanadata->AddTrig("fcsJP2");
    fwdanadata->AddTrig("fcsEM0");
    fwdanadata->AddTrig("fcsEM1");
    fwdanadata->AddTrig("fcsEM2");
    fwdanadata->AddTrig("fcsEM3");
    fwdanadata->AddTrig("fcsEM0_tpc");
    fwdanadata->AddTrig("fcsEM1_tpc");
    fwdanadata->AddTrig("fcsEM2_tpc");
    fwdanadata->AddTrig("fcsEM3_tpc");
    fwdanadata->AddTrig("fcsEHT-N/S");
    fwdanadata->AddTrig("fcsDYAsy");
    fwdanadata->AddTrig("fcsDY");
    fwdanadata->AddTrig("fcsDiJPAsy");
    fwdanadata->AddTrig("fcsDiJP");
    //Don't look at hadron triggers for pi0 analysis
    //fwdanadata->AddTrig("fcsHad0");
    //fwdanadata->AddTrig("fcsHad1");
    //fwdanadata->AddTrig("fcsHad2");
    
    StFwdAnaDataMaker* fwddatamkr = new StFwdAnaDataMaker();
    fwddatamkr->setAnaData(fwdanadata);
    fwddatamkr->setPolDataFilename("Run22PolForJobs.txt");
    fwddatamkr->setFcsTrigFilename("FcsSortedTrig.txt");
    //fwddatamkr->setOutFilename(filenametree.Data();//Only do this if an external hist manager was not declared and the file was not initialized like above 
    fwddatamkr->setHistManager(treehists);
  
    fwddatamkr->addAna(new StFwdAnaPolarization());
    //fwddatamkr->addAna(new StFwdAnaSpin());

    fwddatamkr->addAna(new StFwdAnaFstRun22Qa());
    fwddatamkr->addAna(new StFwdAnaFttRun22Qa());
  
    StFwdAnaFcsRun22Qa* fcsqa = new StFwdAnaFcsRun22Qa();
    fcsqa->setFcsAdcTbOn(false);
    fcsqa->setEpdAdcQaOn(false);
    fcsqa->setEpdTacQaOn(false);
    fcsqa->setBestMassOn(false);
    fwddatamkr->addAna(fcsqa);
    //fwddatamkr->addAna(new StMuFcsAnaCheckFillClusPoint());  //This is for checking mudst clusters/points against clustermaker/pointmaker
  
    StFwdAnaEpdQaAndVert* epdqa = new StFwdAnaEpdQaAndVert();
    epdqa->setEpdTacAdcOn(false);
    fwddatamkr->addAna(epdqa);
    
    fwddatamkr->addAna(new StFwdAnaVertex());
  }
  /*******************************************************************************************/

  // gMessMgr->MemoryOff();

  /*******************************************************************************************/
  // Initialize chain
  chain->SetDebug(1);
  Int_t iInit = chain->Init();
  chain->SetDebug(1);
  cout << "CHAIN INIT DONE? (good==0): " << iInit << endl;
  // ensure that the chain initializes
  
  if ( iInit )
    chain->Fatal(iInit,"on init");
  
  // print the chain status
  chain->PrintInfo();


  // Read 1st event from MuDst to get run number and event time 
  // makes sure that the DB is set with the correct timestamps
  //fcsDb->InitRun(run); //not sure why I need to call this separately...
  
  //chain->InitRun(run);//This is a noop since StChain has no InitRun and StMaker::InitRun() does nothing, InitRun() of submakers are called in StMaker::Make()
  //StFcsDb* fcsDb = (StFcsDb*) chain->GetDataSet("fcsDb");
  //fcsDb->InitRun(run); //not sure why I need to call this separately...
  //StFcsDb* fcsdb = (StFcsDb*) chain->GetDataSet("fcsDb");
  //make sure we get good values
  /*
  for(int d=0; d<4; d++){
    StThreeVectorD off=fcsDb->getDetectorOffset(d);
    printf("FCS Offset d=%1d %8.2f  %8.2f  %8.2f\n",d,off.x(),off.y(),off.z());
  }
  printf("FCS Gain 352 R17c0 = %8.4f\n",fcsDb->getGainCorrection(0,352));
  printf("FCS Gain 374 R18c0 = %8.4f\n",fcsDb->getGainCorrection(0,374));

  */
  StMemStat stmem;
  stmem.PrintMem("BEFORE Event Loop");

  chain->EventLoop(nSkip,nEvents);
  
  /*******************************************************************************************/
  // OPTIONAL: skip first nSkip events (advance the input MuDst without running any
  // downstream maker work, so per-event accumulation from those events does NOT happen).
  // Use this to test whether a crash is event-content driven (always at the same input
  // event #) or accumulation driven (always after N processed events).
  /*******************************************************************************************/
  /*
  size_t totalEntries = muDstChain.GetEntries();
  if (nSkip >= totalEntries) {
    cout << "ERROR: nSkip (" << nSkip << ") >= total entries (" << totalEntries << "). Nothing to process." << endl;
    return;
  }
  if (nSkip > 0) {
    cout << "Skipping first " << nSkip << " events (advancing input only, no chain processing)..." << endl;
    for (size_t s = 0; s < nSkip; s++) {
      chain->Clear();
      if (kStOK != muDstMaker->Make()) {
	cout << "ERROR: muDstMaker->Make() failed during skip at s=" << s << endl;
	break;
      }
    }
    stmem.PrintMem(TString::Format("After Skip of %zu events:", nSkip).Data());
    cout << "Skip complete. Main loop will now process input events starting at #" << nSkip << endl;
    }*/
  
  /*******************************************************************************************/
  // MAIN EVENT LOOP
  /*******************************************************************************************/
  
  /*size_t availableAfterSkip = totalEntries - nSkip;
  size_t nEntries = availableAfterSkip;
  if (nEvents > 0 && nEvents < availableAfterSkip) {
    nEntries = nEvents;
    cout << "Limiting to " << nEntries << " events (input #" << nSkip
	 << " .. #" << (nSkip + nEntries - 1) << ")." << endl;
  }
  double initialUsedHeap = stmem.Used();
  size_t numProcessed = 0;
  for (int i = 0; i < nEntries; i++) {
    size_t inputEv = nSkip + i;   // absolute input event number
    printf("Processing event %d of %d (input #%zu)\n", i, nEntries, inputEv);
    if (i > 0) // skip first event to make it consistent
      stmem.Start();
    chain->Clear();
    
    if (fwdTrack)
    fwdTrack->SetDebug(1);

    if (kStOK != chain->Make())
      break;

    if (refillMuDst){
      StEvent *mStEvent = static_cast<StEvent *>(muDstMaker->GetInputDS("StEvent"));
      // muDstMaker->fillFwdTrack( mStEvent);
      fwdQA->Make();
    }
    stmem.PrintMem(TString::Format("After Event %zu:", inputEv).Data());
    if (i > 0)
      stmem.Stop();
    // MipMaker->Make();
    // picoMk->Make();
    cout << "EVENT #" << i << " (input #" << inputEv << ") COMPLETED" << endl;
    double currentUsedHeap = stmem.Used();
    double deltaUsed = currentUsedHeap - initialUsedHeap;
    cout << "Memory used after event #" << inputEv << ": " << currentUsedHeap << " MB (delta: " << deltaUsed << " MB)" << endl;

    }*/
  stmem.PrintMem("After Event Loop");
  stmem.Summary();
  /*******************************************************************************************/

  // Chain Finish
  cout << "FINISH up" << endl;
  chain->Finish();

  // delete chain;
}



void loadLibs(){	
  //gSystem->Load("libStarClassLibrary.so");
  //gSystem->Load("libStarRoot.so");
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
	gSystem->Load("StarMagField");
 
	gSystem->Load("libGeom");
	gSystem->Load("St_g2t");
	
	// Added for Run16 And beyond
	gSystem->Load("libGeom.so");
	
	gSystem->Load("St_base.so");
	gSystem->Load("StUtilities.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("StarAgmlUtil.so");
	gSystem->Load("StarAgmlLib.so");
	gSystem->Load("libStarGeometry.so");
	gSystem->Load("libGeometry.so");
	
	gSystem->Load("xgeometry");
 
	gSystem->Load("St_geant_Maker");


	// needed since I use the StMuTrack
	gSystem->Load("StarClassLibrary");
	gSystem->Load("StStrangeMuDstMaker");
	gSystem->Load("StMuDSTMaker");
	gSystem->Load("StBTofCalibMaker");
	gSystem->Load("StVpdCalibMaker");
	gSystem->Load("StBTofMatchMaker");
	gSystem->Load("StFcsDbMaker");	

	/*******************************************************************************************/
	// loading libraries
	gSystem->Load("StFstUtil");
	gSystem->Load("StEpdDbMaker");
	gSystem->Load("StEpdHitMaker");
	gSystem->Load("StFcsDbMaker");
	gSystem->Load("StFcsRawHitMaker");
	gSystem->Load("StFcsWaveformFitMaker");
	gSystem->Load("StFcsClusterMaker");
	gSystem->Load("StFcsPointMaker");
	gSystem->Load( "StFttDbMaker" );
	gSystem->Load("StFttRawHitMaker");
	gSystem->Load( "StFttHitCalibMaker" );
	gSystem->Load( "StFttClusterMaker" );
	gSystem->Load( "StFttClusterPointMaker" );
	gSystem->Load( "StFttPointMaker" );
	gSystem->Load("libStarGeneratorUtil.so");
	gSystem->Load("libgenfit2");
	gSystem->Load("libKiTrack");
	gSystem->Load("libXMLIO.so");
	gSystem->Load( "StFwdTrackMaker.so" );
	gSystem->Load( "StFwdUtils.so" );
	gSystem->Load("libStEpdUtil.so");

	gSystem->Load("StStarLogger.so");

	/*******************************************************************************************/


}
