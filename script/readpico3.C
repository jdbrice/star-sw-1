#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TText.h"

#include "StPicoEvent/StPicoDstReader.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoFcsHit.h"
#include "StPicoEvent/StPicoFcsCluster.h"
#include "StPicoEvent/StPicoFwdTrack.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

void InitMatch(int, int set=-1);
void RunMatch(StPicoDst*, StFcsDb*,int);
void EndMatch();
void InitDilepton(int, int set=-1);
void RunDilepton(StPicoDst*, StFcsDb*);
void EndDilepton();

static const int mDebug=0;

void readpico3(int n=0, int opt=11, int run=1, int set=0) {
  StFcsDbMaker *fcsDbMk = new StFcsDbMaker();
  fcsDbMk->Init();
  StFcsDb* fcsDb = dynamic_cast<StFcsDb*>(fcsDbMk->GetDataSet("fcsDb"));
  fcsDb->setDbAccess(0);
  fcsDb->InitRun(run);
  std::cout << "FcsDb initializaed for run " << run << std::endl;
  for(int det=0; det<4; det++){
    StThreeVectorD off=fcsDb->getDetectorOffset(det);
    printf("det=%d offset=%6.3f %6.3f %6.3f\n",det,off.x(),off.y(),off.z());
  }

  char inFile[200];
  if (run==0) {
    sprintf(inFile,"picolist_flat.lis");          // electron MC (default)
  } else if(run==2) {
    sprintf(inFile,"picolist_gamma.lis");         // gamma particle gun MC
  } else if(run < 10000000) {
    sprintf(inFile,"picolist_pythia/%d.list",set);
  } else {
    sprintf(inFile,"picolist/%d.%d.lis",run,set);
  }
  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("FcsHits",1);
  picoReader->SetStatus("FcsClusters",1);
  picoReader->SetStatus("FwdTracks",1);
  picoReader->SetStatus("McTrack",1);  // MC truth for momentum resolution studies
  std::cout << "Status has been set" << std::endl;

  if( !picoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "Events in Tree:  "  << eventsInTree << std::endl;
  Long64_t events2read = picoReader->chain()->GetEntries();
  std::cout << "Events in chain: " << events2read << std::endl;
  if(n>0 && events2read>n){
    events2read=n;
    std::cout << "Limit # of event to read to " << n << std::endl;
  }

  if(opt%10    >0)  InitMatch(run,set);
  if(opt%100/10>0)  InitDilepton(run,set);

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {
    if(iEvent%10==0) std::cout << "Working on event #[" << iEvent<< "/" <<
			  events2read << "] = " <<
			  Form("%3.3f %%",(double)iEvent/(double)events2read*100) << std::endl;

    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..."<< std::endl;
      break;
    }

    // Retrieve picoDst
    StPicoDst *dst = picoReader->picoDst();

    // Retrieve event information
    StPicoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..."<< std::endl;
      break;
    }

    if(opt%10    >0)  RunMatch(dst, fcsDb, iEvent);
    if(opt%100/10>0)  RunDilepton(dst,fcsDb);

  } //event loop

  picoReader->Finish();
  if(opt%10    >0)  EndMatch();
  if(opt%100/10>0)  EndDilepton();
}
