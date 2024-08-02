//usr/bin/env root4star -l -b -q  $0; exit $?
// that is a valid shebang to run script as executable

#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;

StChain *chain_;
void analyze_picos(   int n = 10,//92,
			   const char *inFile = "picos_small.list",//"tracking_clustering.picoDst.root", //"test4.picoDst.root",
                    bool realisticSim = false, // enables data-like mode, real track finding and fitting without MC seed
                    //std::string configFile = "StRoot/StFwdTrackMaker/macro/event/event_track.xml",
			const char *geom = "y2023", TString outfile = "out/Jpsi_test_allcomponents.root") {
  //test
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  
  gSystem->Load("StPicoEvent");//
  gSystem->Load("StPicoDstMaker");//
  gSystem->Load("libStarRoot.so" );
  gSystem->Load("libStarClassLibrary.so");
  gSystem->Load("libStarRoot.so");
  // gSystem->Load("St_db_Maker");//testing out   
  // gSystem->Load( "libStFcsDbMaker.so" ) ; //testing out
  
  chain_ = new StChain("in, y2023, useXgeom, AgML, db, picoDst, picoEvt, fwdTrack");
  
  StPicoDstMaker *picoMaker =0X0; 
  StPicoDstMaker::PicoIoMode IoMode = 2;
  cout<<"Input file is "<<inFile<<endl;
  picoMaker = new StPicoDstMaker(IoMode,inFile,"picoDst");
  
  // needed in this wonky spack environment 
  gROOT->SetMacroPath(".:/star-sw/StRoot/macros:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
  /*
  // gSystem->Load("St_db_Maker");//testing out
  St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb"); //testing out    
  dbMk->AddMaker(dbMk);//testing out
  */
  //gSystem->Load( "libStFcsDbMaker.so" ) ; //testing out
  //StFcsDbMaker * fcsDb = new StFcsDbMaker();
  //fcsDb->setDbAccess(1);//testing out
  // chain_->AddMaker(fcsDb);
  // fcsDb->SetDebug(1);
  
    /* from other file
    StFcsDbMaker* fcsdbmkr = (StFcsDbMaker*) chain->GetMaker("fcsDbMkr");
    cout << "fcsdbmkr="<<fcsdbmkr<<endl;
    fcsdbmkr->setDbAccess(1);

    StFcsDb* fcsdb = (StFcsDb*) chain_->GetDataSet("fcsDb");
    cout << "fcsdb="<<fcsdb<<endl;
    //fcsdb->readGainFromText();                                                                                         
    //fcsdb->readGainCorrFromText();                                                                                     
    fcsdb->forceUniformGain(0.0053);
    fcsdb->forceUniformGainCorrection(1.0);
    */
    
/*
    gSystem->Load("StFcsTrackMatchMaker");
    StFcsTrackMatchMaker *match = new StFcsTrackMatchMaker();
    match->setMaxDistance(6,10);
	//TString fcstrkFile(outfile); fcstrkFile.ReplaceAll(".root",".trackMatch.root");
    match->setFileName("fcstrk.root");//fcstrkFile.Data());//"fcstrk.root");
    match->SetDebug(1);
    chain->AddMaker(match);
    //chain->AddAfter("fwdTrack", match);
*/
    gSystem->Load("StPicoAnalysisTest");
    StPicoAnalysisTest *dilep = new StPicoAnalysisTest("ana",picoMaker);
    TString dilepfile(outfile); dilepfile.ReplaceAll(".root",".dilep.root");
    dilep->setFileName(dilepfile.Data());
    chain_->AddMaker(dilep);//maybe?

    chain_->Init();

    int numevents = picoMaker->chain()->GetEntries();

    cout << numevents << " entries in the chain!" << endl;

    for (Int_t i=0; i<n; i++){
      if(i%1==0){
	cout << "Working on eventNumber " << i << endl;
      }
      chain_->Clear();
      int iret = chain_->Make(i);
      if (iret) { cout << "Bad return code!" << iret << endl; break;}
      //      total++;                                                                                                                                                                                            
    }

    //      cout<<"Number of passed events is "<<nminbias<<endl;                                                                                                                                                
    cout << endl;
    cout << "****************************************** " << endl;
    cout << "Work done... now its time to close up shop!"<< endl;
    cout << "****************************************** " << endl;
    chain_->Finish();
    cout << "****************************************** " << endl;

    delete chain_;

    
}
