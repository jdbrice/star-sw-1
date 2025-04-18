/*
  FPE_OFF
  root.exe -q -b -x kfpAnalysis.C
*/
class StGoodTrigger;
//void kfpAnalysis(Int_t N = 1000000, const Char_t *input = "/net/l404/data/fisyak/Pico/2016/125/17125034/st_physics_17125034_raw_5500079.picoDst.root", const Char_t *output = "picoAna.root", const Char_t *triggerSet = "y2016") {
//void kfpAnalysis(Int_t N = 1000000, const Char_t *input = "st_physics_adc_17125034_raw_1000007.femtoDst.root", const Char_t *output = "picoAna.root", const Char_t *triggerSet = "y2016") {
//void kfpAnalysis(Int_t N = 1000000, const Char_t *input = "/star/data01/pwg_tasks/picoDs/*picoDst.root", const Char_t *output = "picoAna.root", const Char_t *triggerSet = "y2011") {
//void kfpAnalysis(Int_t N = 1000000, const Char_t *input = "./st_mtd_19110005_raw_3500041.picoDst.root", const Char_t *output = "picoAna.root", const Char_t *triggerSet = "y2018") {
//void kfpAnalysis(Int_t N = 100000, const Char_t *input = "/gpfs01/star/pwg_tasks/tfg02/2010/11GeV/*picoDst.root", const Char_t *output = "picoAna2011AuAu11.root", const Char_t *triggerSet = "y2011") {
//void kfpAnalysis(Int_t N = 10000000, const Char_t *input = "/net/l401/data/scratch1/reco/2020/TFG19m/RF/11p5GeV/*picoDst.root", const Char_t *output = "picoAna2020AuAu11p5GeV.root", const Char_t *triggerSet = "y2020") {
//void kfpAnalysis(Int_t N = 1000, const Char_t *input = "/gpfs01/star/pwg_tasks/tfg02/2010/11GeV/st_physics_11148001_raw_1010001.picoDst.root", const Char_t *output = "picoAna2011AuAu11.root", const Char_t *triggerSet = "y2011") {
//void kfpAnalysis(Int_t N = 1000, const Char_t *input = "/gpfs01/star/pwg/fisyak/Pico/2010AuAu11/11148001.picoDst.root", const Char_t *output = "picoAna2011AuAu11.root", const Char_t *triggerSet = "y2011") {
//void kfpAnalysis(Int_t N = 10000000, const Char_t *input = "/net/l401/data/scratch1/reco/2020/TFG19m/RF/11p5GeV.B/347/20347034/hlt_20347034_13_02_000.picoDst.root", const Char_t *output = "Ana2020AuAu11p5GeV.root", const Char_t *triggerSet = "y2020", Bool_t idNdx = kFALSE) {
//void kfpAnalysis(Int_t N = 10000000, const Char_t *input = "./*.picoDst.root", const Char_t *output = "Ana.root", const Char_t *triggerSet = "y2022", Bool_t idNdx = kFALSE) {
void kfpAnalysis(Int_t N = 10000000, 
		 //		 const Char_t *input = "/star/data102/reco/production_OO_200GeV_2021/ReversedFullField/P23ib/2021/136/22136010/st_physics_22136010_raw_1500002.picoDst.root", 
		 const Char_t *input = "*.picoDst.root", 
		 const Char_t *output = "Ana.root", const Char_t *triggerSet = "y2019", Bool_t idNdx = kFALSE) {
#if !defined(__CINT__)
  std::cout << "This code cannot be compiled" << std::endl;
#else
  //  gSystem->SetFPEMask(kInvalid | kDivByZero | kOverflow );
  Bool_t isPico = kTRUE;
  if (TString(input).Contains("MuDst")) isPico = kFALSE;
  gROOT->LoadMacro("lMuDst.C");
  TString Chain("r");
  Chain += triggerSet;
  //  if (! isPico) Chain += ",RMuDst";
  if (! isPico) {Chain += ",RMuDst";} //,PicoWrite"; isPico = kTRUE;}
  else          {Chain += ",RpicoDst";}
  Chain += ",kfpAna,mysql,detDb,nodefault,LdEdxY2,quiet,MuDST";
  //  lMuDst(0,input,"ry2016,RpicoDst,mysql,PicoAnalysis,quiet,nodefault",output);
  lMuDst(-1,input,Chain,output);
//________________________________________________________________________________
// from  /gpfs01/star/pwg/mzyzak/2019_2_2/Template/femtoAnalysis.C
//   StKFParticleInterface::instance()->SetChiPrimaryCut(10);
//   StKFParticleInterface::instance()->SetPtCutCharm(0.7);
//   StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(2);
//   StKFParticleInterface::instance()->SetSoftKaonPIDMode();
//   StKFParticleInterface::instance()->SetSoftTofPidMode();
  std::cout << "KFParticleAnalysis: running analysis for triggerSet " << triggerSet << "." << std::endl; 
  StKFParticleAnalysisMaker* kfpAnalysis = (StKFParticleAnalysisMaker*) StMaker::GetTopChain()->Maker("KFParticleAnalysis");
  if (! isPico) kfpAnalysis->ProcessSignal();  // requires MC info
  kfpAnalysis->AnalyseDsPhiPi();
  kfpAnalysis->CollectPIDHistograms();
  kfpAnalysis->CollectPVHistograms();
  kfpAnalysis->CollectTrackHistograms();
  
  kfpAnalysis->AddDecayToReconstructionList( 310);    // K0
  //  kfpAnalysis->AddDecayToReconstructionList( 100321); // K -> 3pi
  //  kfpAnalysis->AddDecayToReconstructionList(-100321); 
  //  kfpAnalysis->AddDecayToReconstructionList( 200321); // K+3piK
  //  kfpAnalysis->AddDecayToReconstructionList(-200321);
  kfpAnalysis->AddDecayToReconstructionList( 3122); // Lambda
  kfpAnalysis->AddDecayToReconstructionList(-3122); // Lambda_bar
  kfpAnalysis->AddDecayToReconstructionList( 3312); // Xi-
  kfpAnalysis->AddDecayToReconstructionList(-3312); // Xi+
  kfpAnalysis->AddDecayToReconstructionList( 3334); // Omega-
  kfpAnalysis->AddDecayToReconstructionList(-3334); // Omega+

  kfpAnalysis->AddDecayToReconstructionList( 22);   // gamma
  kfpAnalysis->AddDecayToReconstructionList( 111);  // pi0
  kfpAnalysis->AddDecayToReconstructionList( 333);  // phi
#if 0
  //   
  kfpAnalysis->AddDecayToReconstructionList( 113);  // rho
  kfpAnalysis->AddDecayToReconstructionList( 313);  // K*0
  kfpAnalysis->AddDecayToReconstructionList(-313);  // K*0_bar
  kfpAnalysis->AddDecayToReconstructionList( 323);  // K*+
  kfpAnalysis->AddDecayToReconstructionList(-323);  // K*-
  
  kfpAnalysis->AddDecayToReconstructionList( 3212);  // Sigma0
  kfpAnalysis->AddDecayToReconstructionList( 3124); // Lambda*
  kfpAnalysis->AddDecayToReconstructionList( 3124); // Lambda*
  kfpAnalysis->AddDecayToReconstructionList( 3224); // Sigma*+
  kfpAnalysis->AddDecayToReconstructionList( 3114); // Sigma*-
  kfpAnalysis->AddDecayToReconstructionList( 3324); // Xi*0
  kfpAnalysis->AddDecayToReconstructionList(-3324); // Xi*0
#endif  
//   Hyoer Nuclears 
  kfpAnalysis->AddDecayToReconstructionList( 3004); // H3L
  kfpAnalysis->AddDecayToReconstructionList( 3005); // H4L
  kfpAnalysis->AddDecayToReconstructionList( 3006); // He4L
  kfpAnalysis->AddDecayToReconstructionList( 3007); // He5L
#if 0 
  kfpAnalysis->AddDecayToReconstructionList( 3000); // LL
  kfpAnalysis->AddDecayToReconstructionList( 3001); // Lppi
  kfpAnalysis->AddDecayToReconstructionList( 3003); // Ln
  //  kfpAnalysis->AddDecayToReconstructionList(-3003); // Ln_bar
  kfpAnalysis->AddDecayToReconstructionList( 3103); // Lnn
  kfpAnalysis->AddDecayToReconstructionList( 3203); // LLn
  kfpAnalysis->AddDecayToReconstructionList( 3008); // H4LL
  kfpAnalysis->AddDecayToReconstructionList( 3009); // H4LL
  
  kfpAnalysis->AddDecayToReconstructionList( 3012);  // H3L_{dppi}
  kfpAnalysis->AddDecayToReconstructionList( 3013);  // H4L_{tppi}
#endif  
  
#if 0  
  kfpAnalysis->AddDecayToReconstructionList( 9001); // pi+pi+
  kfpAnalysis->AddDecayToReconstructionList(-9001); // pi-pi-
  kfpAnalysis->AddDecayToReconstructionList( 9002); // pi+K+
  kfpAnalysis->AddDecayToReconstructionList(-9002); // pi-K-
  kfpAnalysis->AddDecayToReconstructionList( 9003); // K+K+
  kfpAnalysis->AddDecayToReconstructionList(-9003); // K-K-
  kfpAnalysis->AddDecayToReconstructionList( 9004); // pK+
  kfpAnalysis->AddDecayToReconstructionList(-9004); // p-K-
  kfpAnalysis->AddDecayToReconstructionList( 1003004); // H3L*
  kfpAnalysis->AddDecayToReconstructionList( 1003005); // H4L*
  kfpAnalysis->AddDecayToReconstructionList( 1003006); // He4L*
  kfpAnalysis->AddDecayToReconstructionList( 1003007); // He5L*
#endif  

  StMaker *dbMk = chain->GetMaker("db");
  if (dbMk) {
    dbMk->SetDebug(1);
  }
  chain->Init();
#if 0
  //   StKFParticleInterface::instance()->SetUsedx2(kTRUE); // old dE/dx calibration before SL
  //  if(isPico)
  //  {
  if (triggerSet != "y2022") {
    StKFParticleInterface::instance()->CleanLowPVTrackEvents();
    StMuDst::SetMaxTrackDca(10.0); // dEdxW24 for pp510GeV_2022
  }
  //     StKFParticleInterface::instance()->UseHFTTracksOnly();
  //}
  
  StKFParticleInterface::instance()->SetSoftKaonPIDMode();
  StKFParticleInterface::instance()->SetSoftTofPidMode();

//   StKFParticleInterface::instance()->SetChiPrimaryMaxCut(5.e3);
//  StKFParticleInterface::instance()->SetChiPrimaryCut(12);
  StKFParticleInterface::instance()->SetChiPrimaryCut(18.f);
  StKFParticleInterface::instance()->SetChiPrimaryCutFragments(8.f);
  
  StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(1);
  StKFParticleInterface::instance()->SetLCut(0.f);
  
  //  StKFParticleInterface::instance()->SetChiPrimaryCut2D(8);
  StKFParticleInterface::instance()->SetChiPrimaryCut2D(0);
  StKFParticleInterface::instance()->SetChi2Cut2D(3);
  StKFParticleInterface::instance()->SetLdLCut2D(5);
  
  StKFParticleInterface::instance()->SetChi2CutXiOmega(3);
  StKFParticleInterface::instance()->SetChi2TopoCutXiOmega(3);
  StKFParticleInterface::instance()->SetLdLCutXiOmega(5);  

  StKFParticleInterface::instance()->SetChi2CutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetLdLCutCharmManybodyDecays(5);
#else /* Maksym   /gpfs01/star/pwg/mzyzak/kaons/Template/femtoAnalysis.C Oct 27 2024 */
  if(isPico)
  {
    StKFParticleInterface::instance()->CleanLowPVTrackEvents();
//     StKFParticleInterface::instance()->UseHFTTracksOnly();
  }
  
  StKFParticleInterface::instance()->SetSoftKaonPIDMode();
  StKFParticleInterface::instance()->SetSoftTofPidMode();


//   StKFParticleInterface::instance()->SetChiPrimaryMaxCut(5.e3);
  StKFParticleInterface::instance()->SetChiPrimaryCut(18.f);
  // StKFParticleInterface::instance()->SetChiPrimaryCutFragments(0.f); // TODO
  StKFParticleInterface::instance()->SetChiPrimaryCutFragments(18.f);   // TODO
  // StKFParticleInterface::instance()->SetSecondaryCuts(3.f, 5.f, 0.f);   // TODO

  StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(1);
  StKFParticleInterface::instance()->SetLCut(0.f);
  
  StKFParticleInterface::instance()->SetChiPrimaryCut2D(0);
  StKFParticleInterface::instance()->SetChi2Cut2D(3);
  StKFParticleInterface::instance()->SetLdLCut2D(5.f); // TODO
  // StKFParticleInterface::instance()->SetLdLCut2D(0); // TODO
  
  StKFParticleInterface::instance()->SetChi2CutXiOmega(3);
  StKFParticleInterface::instance()->SetChi2TopoCutXiOmega(3);
  StKFParticleInterface::instance()->SetLdLCutXiOmega(5);  
  
  StKFParticleInterface::instance()->SetChi2CutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetLdLCutCharmManybodyDecays(5);

  StKFParticleInterface::instance()->SetBeamSpot("2020_5AGeV");
  //  StKFParticleInterface::instance()->SetBeamSpot();
#endif
  
//   StKFParticleInterface::instance()->SetSecondaryCuts(3, 3, 5);
//________________________________________________________________________________

  if (idNdx) {
    cout << "StKFParticleInterface::instance()->SetdEdXType(2); // dNdx" << endl;
    StKFParticleInterface::instance()->SetdEdXType(2); // dNdx
  }
  TTree *tree = 0;
  if (! isPico) {
    StKFParticleAnalysisMaker *ana = ( StKFParticleAnalysisMaker *) chain->Maker("KFParticleAnalysis");
    if (! ana) return;
    ana->AnalyseMuDst();
    StMuDstMaker * MuMk = (StMuDstMaker *) StMaker::GetTopChain()->Maker("MuDst");
    if (! MuMk) return;
    MuMk->SetStatus("*",1);
    tree = MuMk->chain();
    // enum PicoVtxMode {NotSet=0, Default=1, Vpd=2, VpdOrDefault=3, Mtd=4, FXT=5};
    StMuDst::setVtxMode(3);    // 2019 AuAu 19GeV
    StMuDst::SetVxXYrange(-0.4,0.4, -0.6, 0.0);
    StMuDst::SetVxZrange(-150.,150.);
  } else {
    StPicoDstMaker * picoMk = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
    if (! picoMk) return;
    picoMk->SetStatus("*",1);
    tree = picoMk->chain();
  }
  if (! tree ) {
    cout << "No MuDst/PicoDst tree. Exit." << endl;
    return;
  }
  Long64_t nentries = tree->GetEntries();
  cout << "no. events in tree. " <<nentries << endl;
#if 1
  if (nentries <= 0) {
    cout << "Tree is empty. Exit" << endl;
    return;
  }
  Long64_t nevent = N;
  nevent = TMath::Min(nevent,nentries);
  cout << nentries << " events in chain " << nevent << " will be read." << endl;
#else
  Long64_t nevent= N;
#endif
  //  new StGoodTrigger(triggerSet);
  //  chain->SetAttr(".Privilege",1,"StPicoDstMaker::*");
  chain->EventLoop(nevent);
#endif
  
}
