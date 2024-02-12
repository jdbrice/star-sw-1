#include "StFcsDiLeptonMaker.h"
#include "TDataSetIter.h"
#include "StDAQMaker/StDAQReader.h"

#include "StRoot/StEvent/StEvent.h"
#include "StRoot/St_base/StMessMgr.h"
#include "StRoot/StEvent/StFcsCollection.h"
#include "StRoot/StEvent/StFcsHit.h"
#include "StRoot/StEvent/StFcsCluster.h"
#include "StRoot/StFcsDbMaker/StFcsDb.h"
#include "StRoot/StEvent/StFwdTrack.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"

#include <string.h>
#include <time.h>

ClassImp(StFcsDiLeptonMaker)

StFcsDiLeptonMaker::StFcsDiLeptonMaker(const char *name):StMaker(name){}

StFcsDiLeptonMaker::~StFcsDiLeptonMaker(){};
 
Int_t StFcsDiLeptonMaker::Init(){  
  mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
  if(!mFcsDb){
    LOG_WARN << "StFcsDiLeptonMaker::Init cannot get fcsDb" << endm;
    return kStErr;
  }

  if(mFilename){
    LOG_INFO << "StFcsDiLeptonMaker::Init - Opening "<<mFilename<<endm;
    mFile=new TFile(mFilename,"RECREATE");
  }

  mE_allclusters = new TH1F("mE_allclusters","",100,0.0,15.0);

  const char* nameCut[mNCut] = {"All","ETOT","HTOT","Cone","SigmaMax","TrackMatch","ChargeSign"};
  
  for(int cut=0; cut<mNCut; cut++){
    mETot[cut]  = new TH1F(Form("RETot_%s", nameCut[cut]),Form("Epair/ETOT_%s",nameCut[cut]),50,0.0,1.1);
    mHTot[cut]  = new TH1F(Form("RHTot_%s", nameCut[cut]),Form("Epair/HTOT_%s",nameCut[cut]),50,0.0,3.0);
    mCone[cut]  = new TH1F(Form("RCone_%s", nameCut[cut]),Form("Epair/Cone_%s (R=%3.1f)",nameCut[cut],mConeR),50,0.0,1.1);
    mSigmax[cut]= new TH1F(Form("Sigmax_%s",nameCut[cut]),Form("SigmaMax_%s",  nameCut[cut]),50,0.0,2.0);
    mPToverET[cut] = new TH1F(Form("EToverPT_%s",nameCut[cut]),Form("EcalET/TrackPT_%s",nameCut[cut]),50,0.0,3.0);
    mChargeSum[cut]= new TH1F(Form("ChargeSum_%s",nameCut[cut]),Form("ChargeSum_%s",nameCut[cut]),5,-2.5,2.5); 

    mET[cut]   = new TH1F(Form("ET_%s"  ,nameCut[cut]),Form("ET_%s"   ,nameCut[cut]),50,0.0,5.0);
    mEZ[cut]   = new TH1F(Form("EZ_%s"  ,nameCut[cut]),Form("EZ_%s"   ,nameCut[cut]),50,0.0,120.0);
    mM[cut]    = new TH1F(Form("M_%s"   ,nameCut[cut]),Form("Mass_%s" ,nameCut[cut]),50,0.0,10.0);
    mZ[cut]    = new TH1F(Form("Z_%s"   ,nameCut[cut]),Form("Zll_%s"  ,nameCut[cut]),50,0.0,1.0);
    mCosT[cut] = new TH1F(Form("CosT_%s",nameCut[cut]),Form("CosT_%s" ,nameCut[cut]),50,-1.0,1.0);
    mPhi[cut]  = new TH1F(Form("Phi_%s" ,nameCut[cut]),Form("Phi_%s"  ,nameCut[cut]),50,-M_PI,M_PI);

    mXFPT[cut] = new TH2F(Form("XFPT_%s",nameCut[cut]),Form("XFPT12_%s;xF;ET",nameCut[cut]),50,0.0,0.5,50,0.0,8.0);
    mET12[cut] = new TH2F(Form("ET12_%s",nameCut[cut]),Form("ET12_%s;ET1;ET2",nameCut[cut]),50,0.0,8.0,50,0.0,8.0);
    mXY[cut]   = new TH2F(Form("XY_%s"  ,nameCut[cut]),Form("XY12_%s;X;Y"    ,nameCut[cut]),50,-130,130,50,-110,110);
    mPTET[cut] = new TH2F(Form("PTET_%s",nameCut[cut]),Form("ETvsPT_%s; ET(Ecal); TrkPT",nameCut[cut]),50,0,8,50,0,8);
  }
  return kStOk;
}

Int_t StFcsDiLeptonMaker::Make(){
  mFcsCollection=0;
  
  //Isaac's test:                                                                                                                    
  //StEvent *event = (StEvent *)GetDataSet("StEvent");//redundant/unnecessary since unused, but leaving for now.
  
  //Isaac's test:
  StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
  if (!stEvent) {
    LOG_INFO << "No StEvent found" << endm;
    return kStErr;
  }
  /*
  StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
  if (!ftc)
    return;
  */
  LOG_INFO << "Checking FcsCollection" << endm;
  mFcsCollection = stEvent->fcsCollection();
  if (!mFcsCollection) {
    LOG_INFO << "No StFcsCollection found" << endm;
    return kStErr;
  }

  //Isaac test:
  if (mFcsCollection->numberOfClusters(0)!=0 || mFcsCollection->numberOfClusters(1)!=0) {
    LOG_INFO << "have a cluster: " <<  mFcsCollection->numberOfClusters(0) << " " << mFcsCollection->numberOfClusters(1) << endm;
  }

  //No clusters found
  if(mFcsCollection->numberOfClusters(0)==0 || mFcsCollection->numberOfClusters(1)==0) return kStOK;


  //Find highest ET clusters for north and south
  StFcsCluster* highest[2]={0,0};
  for(int ns=0; ns<2; ns++){
    //ISAAC TEST! //uncomment this if you ever comment out the "no clusters found" check above
    //if (mFcsCollection->numberOfClusters(ns)==0) {continue;} 
    StSPtrVecFcsCluster& ecal= mFcsCollection->clusters(ns);        
    //sort by ET
    std::sort(ecal.begin(), ecal.end(), [](StFcsCluster* a, StFcsCluster* b) {
        return b->fourMomentum().perp() < a->fourMomentum().perp();
      });    
    
    //keep highest
    //comment out "if...mETCut)", if you want no Et cut for tests
    if(ecal[0]->fourMomentum().perp() > mETCut) highest[ns] = ecal[0];
          
    //ISAAC TEST: leaving commented out unless needed for tests to avoid an expensive loop
    /*
    for (int i = 0; i < ecal.size(); ++ i) {
      mE_allclusters->Fill(ecal[i]->energy());//fills energy of all clusters on both sides
    }
    */
  }  

  
  //  LOG_INFO << "highest[0] = " << highest[0] << " highest[1] " << highest[1] << endm;
  
  //No lepton pair candidates found above mETCut
  if(highest[0]==0 || highest[1]==0) return kStOK; 

  //3 vectors for lepton candidates
  StThreeVectorD VN = mFcsDb->getStarXYZ(highest[0]);
  StThreeVectorD VS = mFcsDb->getStarXYZ(highest[1]);

  //Getting TOT & Cone for isolation cut
  float tot[2] = {0,0}; //eh
  float cone[2]= {0,0}; //ns
  double eta[2],phi[2];
  eta[0]=VN.pseudoRapidity();
  phi[0]=VN.phi();
  eta[1]=VS.pseudoRapidity();
  phi[1]=VS.phi();
  for(int eh=0; eh<2; eh++){    
    for(int ns=0; ns<2; ns++){
      int det=eh*2 + ns;
      StSPtrVecFcsHit& hits= mFcsCollection->hits(det);    
      int n=mFcsCollection->numberOfHits(det);
      for(int i=0; i<n; i++) {
	tot[eh] += hits[i]->energy();
	StThreeVectorD v = mFcsDb->getStarXYZ(hits[i]);
	double e=v.pseudoRapidity();
	double p=v.phi();
	double deta = e-eta[ns];
	double dphi = p-phi[ns];
	while(dphi> M_PI) {dphi -= 2*M_PI;}
	while(dphi<-M_PI) {dphi += 2*M_PI;}
	double dr = sqrt(deta*deta + dphi*dphi);
	if(dr < mConeR) cone[ns] += hits[i]->energy();
      }
    }
  }

  //2 body decay kinematics
  StLorentzVectorD ln = highest[0]->fourMomentum();
  StLorentzVectorD ls = highest[1]->fourMomentum();
  StLorentzVectorD di = ln + ls;
  StLorentzVectorD bln = ln.boost(-di);
  StLorentzVectorD bls = ls.boost(-di);
  double EN  = ln.e();
  double ES  = ls.e();
  double E   = di.e();
  double ETN = ln.perp();
  double ETS = ls.perp();
  double ET  = di.perp();
  double EZ  = di.pz();
  double M   = di.m();
  double Z   = abs(EN-ES)/(EN+ES);
  double CosTN = bln.cosTheta();
  double CosTS = bls.cosTheta();
  double CosT  = CosTN; //take north one for now... When we have tracking, take positive charged
  double Phi   = di.phi(); 
  LOG_DEBUG << Form("AAA CosTheta N=%7.4f S=%7.4f",CosTN,CosTS) << endm;

  //Ecal cluster SigmaMax
  double SigmaMaxN = highest[0]->sigmaMax();
  double SigmaMaxS = highest[1]->sigmaMax();
  
  //Ratio of DiLepton candidate to TOT  
  double ratioETOT = E/tot[0];
  double ratioHTOT = 9.99;
  if(tot[1]>0) ratioHTOT=E/tot[1];

  //Ratio of DiLepton candidates to cone
  double ratioConeN = EN/cone[0];
  double ratioConeS = ES/cone[1];
  
  //top pT associated track
  StFwdTrack *trk1=0, *trk2=0;
  float pt1=0,pt2=0,r1=0,r2=0;
  int cg1=0, cg2=0;
  if(highest[0]->tracks().size()>0) {
    trk1=highest[0]->tracks()[0]; 
    pt1=trk1->momentum().perp(); 
    r1=ETN/pt1;
    cg1=trk1->charge();
  }
  if(highest[1]->tracks().size()>0) {
    trk2=highest[1]->tracks()[0]; 
    pt2=trk2->momentum().perp(); 
    r2=ETS/pt2;
    cg2=trk2->charge();
  } 
  LOG_INFO << Form("trk1 pt=%6.2f et=%6.2f R=%6.4f cg=%2d",pt1,ETN,r1,cg1)<<endm;
  LOG_INFO << Form("trk2 pt=%6.2f et=%6.2f R=%6.4f cg=%2d",pt2,ETS,r2,cg2)<<endm;

  for(int cut=0; cut<mNCut; cut++){
    if(cut==1 && ratioETOT<mETotCut) break;
    if(cut==2 && ratioHTOT<mHTotCut) break;
    if(cut==3 && (ratioConeN<mConeCut || ratioConeS<mConeCut)) break;    
    if(cut==4 && (SigmaMaxN > mSigmaMaxCut || SigmaMaxS > mSigmaMaxCut) ) break;
    if(cut==5 && (r1 < mETPTCutLow || r2 < mETPTCutLow || r1 > mETPTCutHigh || r2 > mETPTCutHigh )) break;    
    if(cut==6 && cg1 + cg2 != 0) break;

    mETot[cut]->Fill(ratioETOT);
    mHTot[cut]->Fill(ratioHTOT);
    mCone[cut]->Fill(ratioConeN);
    mCone[cut]->Fill(ratioConeS);
    mSigmax[cut]->Fill(SigmaMaxN); 
    mSigmax[cut]->Fill(SigmaMaxS); 
    if(trk1) mPToverET[cut]->Fill(r1);
    if(trk2) mPToverET[cut]->Fill(r2);
    if(trk1 && trk2) mChargeSum[cut]->Fill(cg1 + cg2);

    mET  [cut]->Fill(ET);
    mEZ  [cut]->Fill(EZ);
    mM   [cut]->Fill(M);
    mZ   [cut]->Fill(Z);
    mCosT[cut]->Fill(CosT);
    mPhi [cut]->Fill(Phi);

    mET12[cut]->Fill(ETN,ETS);
    mXFPT[cut]->Fill(EN/255.0,ETN);
    mXFPT[cut]->Fill(ES/255.0,ETS);
    mXY[cut]->Fill(VN.x(),VN.y());
    mXY[cut]->Fill(VS.x(),VS.y());
    if(trk1) mPTET[cut]->Fill(ETN,pt1);
    if(trk2) mPTET[cut]->Fill(ETS,pt2);
  }

  return kStOK; 
}
  
Int_t StFcsDiLeptonMaker::Finish(){
  mFile->Write();
  mFile->Close();
  printf("StFcsDiLeptonMaker::Finish - Closing %s\n",mFilename);
  return kStOK;
};

ClassImp(StFcsDiLeptonMaker);

