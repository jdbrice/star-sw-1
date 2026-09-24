//These classes have to be loaded in the source file since dictionary generation of StSPtrVecFtt* can't happen in the header
#include "StEvent/StFttCollection.h"
#include "StEvent/StFttRawHit.h"
#include "StEvent/StFttCluster.h"
#include "StEvent/StFttPoint.h"

#include "StFwdAnaFttRun22Qa.h"

ClassImp(StFwdAnaFttRun22Qa)

StFwdAnaFttRun22Qa::StFwdAnaFttRun22Qa()
{
  memset(mH2S_FttRawHit_adcVch,   0, sizeof(mH2S_FttRawHit_adcVch));
  memset(mH2S_FttRawHit_bcidVch,  0, sizeof(mH2S_FttRawHit_bcidVch));
  memset(mH2S_FttRawHit_dbcidVch, 0, sizeof(mH2S_FttRawHit_dbcidVch));
  memset(mH2S_FttRawHit_tbVch,    0, sizeof(mH2S_FttRawHit_tbVch));
  memset(mH2S_FttRawHit_timeVch,  0, sizeof(mH2S_FttRawHit_timeVch));
}

StFwdAnaFttRun22Qa::~StFwdAnaFttRun22Qa()
{}

UInt_t StFwdAnaFttRun22Qa::LoadHists(TFile* file, HistManager* histman, StFwdAnaData* anadata)
{
  UInt_t nloaded = 0;
  std::stringstream ss_histname;
  std::stringstream ss_histtitle;
  for( unsigned int ivmm=0; ivmm<StFttDb::nVMM; ++ivmm ){
    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_FttRawHit_adcVch_"<< ivmm;
    ss_histtitle << "Raw Hit ADC vs. channel for VMM ID "<< ivmm << ";ch;adc";
    nloaded += histman->AddH2S(file,mH2S_FttRawHit_adcVch[ivmm],ss_histname.str().c_str(),ss_histtitle.str().c_str(), StFttDb::nChPerVMM,0,StFttDb::nChPerVMM, StFttDb::maxADC-1,0,StFttDb::maxADC-1);  //@[August 10, 2026] > For some reason in StFttDb maxADC is set to 1025 even though it is 10 bits which means max should be 1024

    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_FttRawHit_bcidVch_"<< ivmm;
    ss_histtitle << "Raw Hit BCID vs. channel for VMM ID "<< ivmm << ";ch;bcid";
    nloaded += histman->AddH2S(file,mH2S_FttRawHit_bcidVch[ivmm],ss_histname.str().c_str(),ss_histtitle.str().c_str(), StFttDb::nChPerVMM,0,StFttDb::nChPerVMM, StFttDb::maxBCID-1,0,StFttDb::maxBCID-1);  //@[August 10, 2026] > For some reason in StFttDb maxBCID is set to 4097 even though it is 12 bits which means max should be 4096

    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_FttRawHit_dbcidVch_"<< ivmm;
    ss_histtitle << "Raw Hit Delta BCID vs. channel for VMM ID "<< ivmm << ";ch;dbcid";
    nloaded += histman->AddH2S(file,mH2S_FttRawHit_dbcidVch[ivmm],ss_histname.str().c_str(),ss_histtitle.str().c_str(), StFttDb::nChPerVMM,0,StFttDb::nChPerVMM, 2048,-4096,4096); //@[August 10, 2026] > Delta BCID is presumably the difference in BCID which has max of 4096 so largest difference is 4096. Use 2048 bins to reduce data size

    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_FttRawHit_tbVch_"<< ivmm;
    ss_histtitle << "Raw Hit timebin(tb) vs. channel for VMM ID "<< ivmm << ";ch;tb";
    nloaded += histman->AddH2S(file,mH2S_FttRawHit_tbVch[ivmm],ss_histname.str().c_str(),ss_histtitle.str().c_str(), StFttDb::nChPerVMM,0,StFttDb::nChPerVMM, 938,StFttDb::minTb,StFttDb::maxTb);  //@[August 10, 2026] > minTB=-32768-1000 and maxTB=32768+1000, don't want 67532 bins so use 938 bins instead (67532/938==72)

    ss_histname.str("");
    ss_histtitle.str("");
    ss_histname << "H2S_FttRawHit_timeVch_"<< ivmm;
    ss_histtitle << "Raw Hit time vs. channel for VMM ID "<< ivmm << ";ch;time";
    nloaded += histman->AddH2S(file,mH2S_FttRawHit_timeVch[ivmm],ss_histname.str().c_str(),ss_histtitle.str().c_str(), StFttDb::nChPerVMM,0,StFttDb::nChPerVMM, 2048,-2048,2048);  //@[August 10, 2026] > From StFttHitCalibMaker/HitCalibHelper.h which wraps around for diferences larger than 2048?
  }

  nloaded += histman->AddH2S(file,mH2S_FttRawHit_nchsVvmm,"H2S_FttRawHit_nchsVvmm", "Raw Hit NChs vs VMM ID;VMM ID;NChs",StFttDb::nVMM,0,StFttDb::nVMM, StFttDb::nChPerVMM,0,StFttDb::nChPerVMM);

  return nloaded;
}


Int_t StFwdAnaFttRun22Qa::DoMake(StFwdAnaData* anadata)
{
  StFttCollection* fttcoll = anadata->fttColl();
  if( fttcoll==0 ){
    LOG_WARN << "StFwdAnaFttRun22Qa::No FTT Hit Collection" << endm;
    return kStWarn;
  }

  //loop over raw hits
  StSPtrVecFttRawHit& hits = fttcoll->rawHits();
  //if( hits==0 ){ LOG_WARN << "StFwdAnaFttRun22Qa::DoMake() No Ftt Hits found" << endm; return kStWarn; }
  //std::cout << "|fttrawhits:"<<hits.size() <<"|vmmsize:"<<StFttDb::nVMM<< std::endl;

  Short_t nchhitcounter[StFttDb::nVMM] = {0};   //For counting number of channels hit for a given vmm id
  for( unsigned int ihit=0; ihit<hits.size(); ++ihit ){
    StFttRawHit* fttraw = hits.at(ihit);
    if( fttraw==0 ){ continue; }
    size_t vmmid = StFttDb::vmmId(fttraw);
    UChar_t ch = fttraw->channel();

    //std::cout << " +|vmmid:"<<vmmid << std::endl;
    ((TH2*)mH2S_FttRawHit_adcVch[vmmid])->Fill(ch,fttraw->adc());
    ((TH2*)mH2S_FttRawHit_bcidVch[vmmid])->Fill(ch,fttraw->bcid());
    ((TH2*)mH2S_FttRawHit_dbcidVch[vmmid])->Fill(ch,fttraw->dbcid());
    ((TH2*)mH2S_FttRawHit_tbVch[vmmid])->Fill(ch,fttraw->tb());
    ((TH2*)mH2S_FttRawHit_timeVch[vmmid])->Fill(ch,fttraw->time());

    ++nchhitcounter[vmmid];  //Since vmm id is unique up to channel number this counts number or channels or number of hits in a given vmm id
  }
  //std::cout << "end for loop" << std::endl;
  for( unsigned int i=0; i<StFttDb::nVMM; ++i ){
    mH2S_FttRawHit_nchsVvmm->Fill(i,nchhitcounter[i]);
  }
  
  return kStOk;
}

void StFwdAnaFttRun22Qa::Paint_RawHitAdc(TCanvas* canvas, int plane, int quad, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(5,5);

  int planequadkey = StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*quad + StFttDb::nQuadPerPlane*StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*plane; //Unpacking vmm_id() to match input "plane", "quad" to correct array index, There are 24 VMMs per quad so plot one vmm per quad per plane on each pad (there are 6 fobs per quad and 4 vmms per fob so a total of 24 vmms per quad). Counting starts from 0
  for( int ipad=0; ipad<24; ++ipad){
    canvas->cd(ipad+1);
    if( mH2S_FttRawHit_adcVch[planequadkey+ipad] ){ mH2S_FttRawHit_adcVch[planequadkey+ipad]->Draw("colz"); }
  }

  canvas->Print(savename);    
}

void StFwdAnaFttRun22Qa::Paint_RawHitBcid(TCanvas* canvas, int plane, int quad, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(5,5);

  int planequadkey = StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*quad + StFttDb::nQuadPerPlane*StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*plane; //Unpacking vmm_id() to match input "plane", "quad" to correct array index, There are 24 VMMs per quad so plot one vmm per quad per plane on each pad (there are 6 fobs per quad and 4 vmms per fob so a total of 24 vmms per quad). Counting starts from 0
  for( int ipad=0; ipad<24; ++ipad){
    canvas->cd(ipad+1);
    if( mH2S_FttRawHit_bcidVch[planequadkey+ipad] ){ mH2S_FttRawHit_bcidVch[planequadkey+ipad]->Draw("colz"); }
  }

  canvas->Print(savename);    
}

void StFwdAnaFttRun22Qa::Paint_RawHitDBcid(TCanvas* canvas, int plane, int quad, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(5,5);

  int planequadkey = StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*quad + StFttDb::nQuadPerPlane*StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*plane; //Unpacking vmm_id() to match input "plane", "quad" to correct array index, There are 24 VMMs per quad so plot one vmm per quad per plane on each pad (there are 6 fobs per quad and 4 vmms per fob so a total of 24 vmms per quad). Counting starts from 0 so [0-24]
  for( int ipad=0; ipad<24; ++ipad){
    canvas->cd(ipad+1);
    if( mH2S_FttRawHit_dbcidVch[planequadkey+ipad] ){ mH2S_FttRawHit_dbcidVch[planequadkey+ipad]->Draw("colz"); }
  }

  canvas->Print(savename);    
}

void StFwdAnaFttRun22Qa::Paint_RawHitTb(TCanvas* canvas, int plane, int quad, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(5,5);

  int planequadkey = StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*quad + StFttDb::nQuadPerPlane*StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*plane; //Unpacking vmm_id() to match input "plane", "quad" to correct array index, There are 24 VMMs per quad so plot one vmm per quad per plane on each pad (there are 6 fobs per quad and 4 vmms per fob so a total of 24 vmms per quad). Counting starts from 0
  for( int ipad=0; ipad<24; ++ipad){
    canvas->cd(ipad+1);
    if( mH2S_FttRawHit_tbVch[planequadkey+ipad] ){ mH2S_FttRawHit_tbVch[planequadkey+ipad]->Draw("colz"); }
  }

  canvas->Print(savename);    
}

void StFwdAnaFttRun22Qa::Paint_RawHitTime(TCanvas* canvas, int plane, int quad, const char* savename) const
{
  canvas->Clear();
  canvas->Divide(5,5);

  int planequadkey = StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*quad + StFttDb::nQuadPerPlane*StFttDb::nFobPerQuad*StFttDb::nVMMPerFob*plane; //Unpacking vmm_id() to match input "plane", "quad" to correct array index, There are 24 VMMs per quad so plot one vmm per quad per plane on each pad (there are 6 fobs per quad and 4 vmms per fob so a total of 24 vmms per quad). Counting starts from 0
  for( int ipad=0; ipad<24; ++ipad){
    canvas->cd(ipad+1);
    if( mH2S_FttRawHit_timeVch[planequadkey+ipad] ){ mH2S_FttRawHit_timeVch[planequadkey+ipad]->Draw("colz"); }
  }

  canvas->Print(savename);    
}

void StFwdAnaFttRun22Qa::Paint_RawHitNchs(TCanvas* canvas, const char* savename) const
{
  canvas->Clear();
  //canvas->cd();
  if( mH2S_FttRawHit_nchsVvmm ){ mH2S_FttRawHit_nchsVvmm->Draw("colz"); }
  canvas->Print(savename);
}

