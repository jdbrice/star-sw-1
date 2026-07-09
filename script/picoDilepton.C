#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include "StPicoEvent/StPicoDstReader.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoFcsHit.h"
#include "StPicoEvent/StPicoFcsCluster.h"
#include "StPicoEvent/StPicoFwdTrack.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

enum {mNCut=7};
const char* nameCut[mNCut] = {"All","ETOT","HTOT","Cone","SigmaMax","TrackMatch","ChargeSign"};
static const int mNType=5;
const char* TTYPE[mNType]={"Global","Beamline","Primary","FwdVtx","BLCVtx"};
float mETCut=0.4;          //GeV for single electron
float mETotCut=0.3;        //E_lepton/ETOT ratio cut
float mHTotCut=0.5;        //E_lepton/HTOT ratio cut
float mConeR=0.5;          //Isolation Cone Radius
float mConeCut=0.6;        //E_lepton/Cone ratio cut
float mSigmaMaxCut=0.7;    //Cluster Sigma Max cut
//float mETPTCutLow=0.3;   //FcsET/TrackPT > ThresholdLow
float mETPTCutLow=0.00001;  //FcsET/TrackPT > ThresholdLow
//float mETPTCutHigh=1.7;    //FcsET/TrackPT < ThresholdHigh
float mETPTCutHigh=10000;    //FcsET/TrackPT < ThresholdHigh
//float mDChi2Cut=65;
float mDChi2Cut=1000000;
float mDRCut=20.0;      // radial dR = R(fcs)-R(trk) cut [cm]
float mRdphiCut=8.0;  // azimuthal arc Rdphi = R(fcs)*(phi(fcs)-phi(trk)) cut [cm]

//float mDNHitCut=4;
int mDNHitCut=4;
int mDebug=1;
static long nEvt=0, nEvtBothClus=0;
static long nEvtCut[mNType][7];
static long nTrkAll=0, nTrkFailStatus=0, nTrkFailNHit=0, nTrkFailChi2=0;
static long nTrkPassDR[mNType][2];
static long nTrkPassPT[mNType][2];

//Mixed-event background: cut level at which an event's north/south FCS cluster
//4-vectors are buffered and mixed with the previous qualifying event, per track
//type. Adjustable; default matches the tightest existing cut (ChargeSign).
int mMixCut=6;
static bool mPrevValid[mNType];
static TLorentzVector mPrevLN[mNType];
static TLorentzVector mPrevLS[mNType];
static long nMixFilled[mNType];

char filenameD[200];
TFile *mFileD;

//define histograms (all indexed [trackType][cut])
TH1F *mZVTX[mNType][mNCut];
TH1F *mZVTXA[mNType][mNCut];
TH1F *mZVTXD[mNType][mNCut];
TH1F *mETot[mNType][mNCut];
TH1F *mHTot[mNType][mNCut];
TH1F *mCone[mNType][mNCut];
TH1F *mSigmax[mNType][mNCut];
TH1F *mPToverET[mNType][mNCut];
TH1F *mChargeSum[mNType][mNCut];
TH1F *mET[mNType][mNCut];
TH1F *mEZ[mNType][mNCut];
TH1F *mM[mNType][mNCut];
TH1F *mMmix[mNType];
TH1F *mZ[mNType][mNCut];
TH1F *mCosT[mNType][mNCut];
TH1F *mPhi[mNType][mNCut];
TH2F *mET12[mNType][mNCut];
TH2F *mXFPT[mNType][mNCut];
TH2F *mXY[mNType][mNCut];
TH2F *mPTET[mNType][mNCut];

void InitDilepton(int run, int set=-1){
  if(set<0) sprintf(filenameD,"hist_dilep/%d.dilep.root",run);
  else       sprintf(filenameD,"hist_dilep/%d.%d.dilep.root",run,set);
  printf("Opening %s\n",filenameD);
  mFileD=new TFile(filenameD,"RECREATE");

  for(int tt=0; tt<mNType; tt++){
    mPrevValid[tt] = false;
    nMixFilled[tt] = 0;
    mMmix[tt] = new TH1F(Form("Mmix_%s",TTYPE[tt]),Form("Mixed-Event Mass %s (cut=%s)",TTYPE[tt],nameCut[mMixCut]),50,0.0,10.0);
    for(int cut=0; cut<mNCut; cut++){
      mETot[tt][cut]     = new TH1F(Form("RETot_%s_%s",    TTYPE[tt],nameCut[cut]),Form("Epair/ETOT %s %s",               TTYPE[tt],nameCut[cut]),50,0.0,1.1);
      mHTot[tt][cut]     = new TH1F(Form("RHTot_%s_%s",    TTYPE[tt],nameCut[cut]),Form("Epair/HTOT %s %s",               TTYPE[tt],nameCut[cut]),50,0.0,3.0);
      mCone[tt][cut]     = new TH1F(Form("RCone_%s_%s",    TTYPE[tt],nameCut[cut]),Form("Epair/Cone %s %s (R=%3.1f)",     TTYPE[tt],nameCut[cut],mConeR),50,0.0,1.1);
      mSigmax[tt][cut]   = new TH1F(Form("Sigmax_%s_%s",   TTYPE[tt],nameCut[cut]),Form("SigmaMax %s %s",                 TTYPE[tt],nameCut[cut]),50,0.0,2.0);
      mPToverET[tt][cut] = new TH1F(Form("EToverPT_%s_%s", TTYPE[tt],nameCut[cut]),Form("EcalET/TrackPT %s %s",           TTYPE[tt],nameCut[cut]),50,0.0,6.0);
      mChargeSum[tt][cut]= new TH1F(Form("ChargeSum_%s_%s",TTYPE[tt],nameCut[cut]),Form("ChargeSum %s %s",                TTYPE[tt],nameCut[cut]),5,-2.5,2.5);
      mZVTX[tt][cut]     = new TH1F(Form("ZVTX_%s_%s",     TTYPE[tt],nameCut[cut]),Form("ZVTX %s %s",                    TTYPE[tt],nameCut[cut]),50,-200,200);
      mZVTXA[tt][cut]    = new TH1F(Form("ZVTXA_%s_%s",    TTYPE[tt],nameCut[cut]),Form("ZVTXAvg %s %s",                 TTYPE[tt],nameCut[cut]),50,-200,200);
      mZVTXD[tt][cut]    = new TH1F(Form("ZVTXD_%s_%s",    TTYPE[tt],nameCut[cut]),Form("ZVTXDif %s %s",                 TTYPE[tt],nameCut[cut]),50,-200,200);
      mET[tt][cut]       = new TH1F(Form("ET_%s_%s",        TTYPE[tt],nameCut[cut]),Form("ET %s %s",                     TTYPE[tt],nameCut[cut]),50,0.0,5.0);
      mEZ[tt][cut]       = new TH1F(Form("EZ_%s_%s",        TTYPE[tt],nameCut[cut]),Form("EZ %s %s",                     TTYPE[tt],nameCut[cut]),50,0.0,120.0);
      mM[tt][cut]        = new TH1F(Form("M_%s_%s",          TTYPE[tt],nameCut[cut]),Form("Mass %s %s",                   TTYPE[tt],nameCut[cut]),50,0.0,10.0);
      mZ[tt][cut]        = new TH1F(Form("Z_%s_%s",          TTYPE[tt],nameCut[cut]),Form("Zll %s %s",                    TTYPE[tt],nameCut[cut]),50,0.0,1.0);
      mCosT[tt][cut]     = new TH1F(Form("CosT_%s_%s",      TTYPE[tt],nameCut[cut]),Form("CosT %s %s",                   TTYPE[tt],nameCut[cut]),50,-1.0,1.0);
      mPhi[tt][cut]      = new TH1F(Form("Phi_%s_%s",        TTYPE[tt],nameCut[cut]),Form("Phi %s %s",                    TTYPE[tt],nameCut[cut]),50,-M_PI,M_PI);
      mXFPT[tt][cut]     = new TH2F(Form("XFPT_%s_%s",      TTYPE[tt],nameCut[cut]),Form("XFPT12 %s %s;xF;ET",           TTYPE[tt],nameCut[cut]),50,0.0,0.5,50,0.0,8.0);
      mET12[tt][cut]     = new TH2F(Form("ET12_%s_%s",       TTYPE[tt],nameCut[cut]),Form("ET12 %s %s;ET1;ET2",            TTYPE[tt],nameCut[cut]),50,0.0,8.0,50,0.0,8.0);
      mXY[tt][cut]       = new TH2F(Form("XY_%s_%s",         TTYPE[tt],nameCut[cut]),Form("XY12 %s %s;X;Y",               TTYPE[tt],nameCut[cut]),50,-130,130,50,-110,110);
      mPTET[tt][cut]     = new TH2F(Form("PTET_%s_%s",       TTYPE[tt],nameCut[cut]),Form("ETvsPT %s %s; ET(Ecal); TrkPT",TTYPE[tt],nameCut[cut]),50,0,8,50,0,8);
    }
  }
}

void RunDilepton(StPicoDst *dst, StFcsDb* fcsDb){
  StPicoEvent *event = dst->event();
  nEvt++;

  //loop over FCS clusters and find higest ET Ecal clusters
  float maxet[2]={mETCut,mETCut};
  StPicoFcsCluster* highest[2]={0,0};
  int nClusters = dst->numberOfFcsClusters();
  for(int ic=0; ic<nClusters; ic++) {
    StPicoFcsCluster* c = dst->fcsCluster(ic);
    int det=c->detectorId();
    int ehp=fcsDb->ecalHcalPres(det);
    if(ehp!=0) continue; //only Ecal
    int ns=fcsDb->northSouth(det);
    //bug in pico dst fcs cluster fourMomentum() calc.... redo calc based on cluster xy
    //double et=c->fourMomentum().Et();
    StThreeVectorD xyz=fcsDb->getStarXYZfromColumnRow(det,c->x(),c->y());
    StLorentzVectorD lv=fcsDb->getLorentzVector(xyz,c->energy(),event->primaryVertex().Z());
    double et=lv.perp();
    if(et>maxet[ns]) {maxet[ns]=et; highest[ns]=c;}
  }

  //No lepton pair candidates found above mETCut
  if(highest[0]==0 || highest[1]==0) return;
  nEvtBothClus++;

  //3 vectors for lepton candidates
  StThreeVectorD FcsXyz[2];
  FcsXyz[0] = fcsDb->getStarXYZfromColumnRow(highest[0]->detectorId(),highest[0]->x(),highest[0]->y());
  FcsXyz[1] = fcsDb->getStarXYZfromColumnRow(highest[1]->detectorId(),highest[1]->x(),highest[1]->y());

  //Getting TOT & Cone for isolation cut
  float tot[2] = {0,0}; //eh
  float cone[2]= {0,0}; //ns
  double eta[2],phi[2];
  eta[0]=FcsXyz[0].pseudoRapidity();
  phi[0]=FcsXyz[0].phi();
  eta[1]=FcsXyz[1].pseudoRapidity();
  phi[1]=FcsXyz[1].phi();
  int nHits = dst->numberOfFcsHits();
  for(int ih=0; ih<nHits; ih++) {
    StPicoFcsHit* hit = dst->fcsHit(ih);
    int det=hit->detectorId();
    int ehp=fcsDb->ecalHcalPres(det);
    if(ehp>=2) continue; //only Ecal and Hcal
    int ns=fcsDb->northSouth(det);
    tot[ehp] += hit->energy();
    StThreeVectorD v = fcsDb->getStarXYZ(hit->detectorId(),hit->id());
    double e=v.pseudoRapidity();
    double p=v.phi();
    double deta = e-eta[ns];
    double dphi = p-phi[ns];
    while(dphi> M_PI) {dphi -= 2*M_PI;}
    while(dphi<-M_PI) {dphi += 2*M_PI;}
    double dr = sqrt(deta*deta + dphi*dphi);
    if(dr < mConeR) cone[ns] += hit->energy();
  }

  //2 body decay kinematics
  //bug in pico dst fcs cluster fourMomentum() calc.... redo calc based on cluster xy
  //TLorentzVector ln = highest[0]->fourMomentum();
  //TLorentzVector ls = highest[1]->fourMomentum();
  StLorentzVectorD stln=fcsDb->getLorentzVector(FcsXyz[0],highest[0]->energy(),event->primaryVertex().Z());
  StLorentzVectorD stls=fcsDb->getLorentzVector(FcsXyz[1],highest[1]->energy(),event->primaryVertex().Z());
  TLorentzVector ln(stln.px(),stln.py(),stln.pz(),stln.e());
  TLorentzVector ls(stls.px(),stls.py(),stls.pz(),stls.e());
  TLorentzVector di = ln + ls;
  //lab-frame copies for mixed-event background (ln/ls get boosted to the pair CM frame below)
  TLorentzVector lnLab = ln;
  TLorentzVector lsLab = ls;
  double EN  = ln.E();
  double ES  = ls.E();
  double E   = di.E();
  double ETN = ln.Perp();
  double ETS = ls.Perp();
  double FcsET[2]={ETN,ETS};
  double ET  = di.Perp();
  double EZ  = di.Pz();
  double M   = di.M();
  double Z   = abs(EN-ES)/(EN+ES);
  double Phi = di.Phi();
  TVector3 boost=-di.BoostVector();
  ln.Boost(boost);
  ls.Boost(boost);
  double CosTN = ln.CosTheta();
  double CosTS = ls.CosTheta();
  double CosT  = CosTN; //take north one for now... When we have tracking, take positive charged
  if(mDebug>0){
    printf("FCS VTX= %8.3f  %8.3f  %8.3f\n",event->primaryVertex().X(),event->primaryVertex().Y(),event->primaryVertex().Z());
    printf("FCS EN %8.3f %8.3f %8.3f E=%8.3f M=%8.3f ET=%8.3f Phi=%8.3f\n",
	   FcsXyz[0].x(),FcsXyz[0].y(),FcsXyz[0].z(),highest[0]->energy(),
	   0.0,FcsXyz[0].perp(),FcsXyz[0].phi());
    printf("FCS ES %8.3f %8.3f %8.3f E=%8.3f M=%8.3f ET=%8.3f Phi=%8.3f\n",
	   FcsXyz[1].x(),FcsXyz[1].y(),FcsXyz[1].z(),highest[1]->energy(),0.0,
	   FcsXyz[1].perp(),FcsXyz[1].phi());
    //printf("BST EN %8.3f %8.3f %8.3f E=%8.3f M=%8.3f ET=%8.3f Phi=%8.3f\n",ln.Px(),ln.Py(),ln.Pz(),ln.E(),ln.M(),ln.Perp(),ln.Phi());
    //printf("BST ES %8.3f %8.3f %8.3f E=%8.3f M=%8.3f ET=%8.3f Phi=%8.3f\n",ls.Px(),ls.Py(),ls.Pz(),ls.E(),ls.M(),ls.Perp(),ls.Phi());
    //printf("Dilep  %8.3f %8.3f %8.3f E=%8.3f M=%8.3f ET=%8.3f Phi=%8.3f\n",di.Px(),di.Py(),di.Pz(),di.E(),di.M(),di.Perp(),di.Phi());
    //printf("CosTheta N=%8.3f S=%8.3f\n",CosTN,CosTS);
  }

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

  //best-pT matched track per (trackType, north/south)
  StPicoFwdTrack *trk[mNType][2];
  double trkpt[mNType][2];
  double etpt[mNType][2];
  double dr[mNType][2];
  double rdphi[mNType][2];
  int cg[mNType][2];
  for(int tt=0; tt<mNType; tt++){
    for(int ns=0; ns<2; ns++){
      trk[tt][ns]=0; trkpt[tt][ns]=0; etpt[tt][ns]=0; dr[tt][ns]=0; rdphi[tt][ns]=0; cg[tt][ns]=0;
    }
  }

  int nTracks = dst->numberOfFwdTracks();
  if(mDebug>1){
    for(int it=0; it<nTracks; it++) {
      StPicoFwdTrack* t=dst->fwdTrack(it);
      if(t->status() < 2)  continue;
      if(t->chi2()==0.0) continue;
      TVector3 proj = t->ecalProjection();
      double tx=proj.X();
      double ty=proj.Y();
      printf("TRK1 id=%3d Typ=%1d chi2=%10.2f st=%1hhu nSd=%2d nFit=%2d Vtx=%3hhu Chg=%2d P=%8.2f %8.2f %8.2f E=%10.2f XY=%8.2f %8.2f DCA=%7.3f %7.1f Ecal=%1d Hcal=%1d\n",
             t->id(),t->trackType(),
             t->chi2(),t->status(),t->numberOfSeedPoints(),abs(t->numberOfFitPoints()),
             t->vertexIndex(), t->charge(), t->momentum().X(), t->momentum().Y(), t->momentum().Z(), t->momentum().Mag(),
	     tx,ty,
             t->dcaXY(),t->dcaZ(),
             t->numberOfEcalMatchIndices(),t->numberOfHcalMatchIndices()
             );
    }
  }
  for(int it=0; it<nTracks; it++) {
    StPicoFwdTrack* t=dst->fwdTrack(it);
    nTrkAll++;
    if(t->status() < 2)                        { nTrkFailStatus++; continue; }
    int tt = t->trackType();
    if(tt < 0 || tt >= mNType)                 { continue; }
    if(abs(t->numberOfFitPoints()) < mDNHitCut) { nTrkFailNHit++;   continue; }
    TVector3 proj = t->ecalProjection();
    double tx=proj.X();
    double ty=proj.Y();
    double pt=t->momentum().Perp();
    int ns=(tx>0) ? 1 : 0;
    double fcsx=FcsXyz[ns].x();
    double fcsy=FcsXyz[ns].y();
    double fcsR  = sqrt(fcsx*fcsx + fcsy*fcsy);
    double dR    = fcsR - sqrt(tx*tx + ty*ty);
    double dphi  = atan2(fcsy,fcsx) - atan2(ty,tx);
    while(dphi >  TMath::Pi()) dphi -= TMath::TwoPi();
    while(dphi < -TMath::Pi()) dphi += TMath::TwoPi();
    double Rdphi = fcsR * dphi;
    if(mDebug>0){
      printf("TRK2 ns=%d tt=%d tx=%7.2f ty=%7.2f fcsX=%7.2f fcsY=%7.2f dR=%7.2f Rdphi=%7.2f pt=%7.2f chi2=%8.2f\n",
             ns,tt,tx,ty,fcsx,fcsy,dR,Rdphi,pt,t->chi2());
    }
    if(t->chi2() > mDChi2Cut) { nTrkFailChi2++;   continue; }
    if(fabs(dR) < mDRCut && fabs(Rdphi) < mRdphiCut) {
      nTrkPassDR[tt][ns]++;
      if(pt > trkpt[tt][ns]) {
        nTrkPassPT[tt][ns]++;
        trk[tt][ns]=t; trkpt[tt][ns]=pt; etpt[tt][ns]=FcsET[ns]/pt; cg[tt][ns]=t->charge();
        dr[tt][ns]=dR; rdphi[tt][ns]=Rdphi;
      }
    }
  }
  if(mDebug>0){
    for(int tt=0; tt<mNType; tt++){
      for(int ns=0; ns<2; ns++){
        if(trk[tt][ns]){
	  printf("trk NS=%1d TT=%1d(%s) et=%6.2f M=%6.2f ET/pT=%6.4f trkid=%2d pt=%12.2f cg=%2d DCA=%7.3f %7.1f dR=%7.2f Rdphi=%7.2f\n",
	         ns,tt,TTYPE[tt],FcsET[ns],M,etpt[tt][ns],trk[tt][ns]->id(),trkpt[tt][ns],cg[tt][ns],
                 trk[tt][ns]->dcaXY(),trk[tt][ns]->dcaZ(),dr[tt][ns],rdphi[tt][ns]);
        }
      }
    }
  }

  //Apply cuts and fill histograms for each track type
  for(int tt=0; tt<mNType; tt++){
    for(int cut=0; cut<mNCut; cut++){
      if(cut==1 && ratioETOT<mETotCut) break;
      if(cut==2 && ratioHTOT<mHTotCut) break;
      if(cut==3 && (ratioConeN<mConeCut || ratioConeS<mConeCut)) break;
      if(cut==4 && (SigmaMaxN > mSigmaMaxCut || SigmaMaxS > mSigmaMaxCut) ) break;
      if(cut==5 && (etpt[tt][0] < mETPTCutLow || etpt[tt][1] < mETPTCutLow || etpt[tt][0] > mETPTCutHigh || etpt[tt][1] > mETPTCutHigh )) break;
      if(cut==6 && cg[tt][0] + cg[tt][1] != 0) break;
      nEvtCut[tt][cut]++;

      mETot[tt][cut]->Fill(ratioETOT);
      mHTot[tt][cut]->Fill(ratioHTOT);
      mCone[tt][cut]->Fill(ratioConeN);
      mCone[tt][cut]->Fill(ratioConeS);
      mSigmax[tt][cut]->Fill(SigmaMaxN);
      mSigmax[tt][cut]->Fill(SigmaMaxS);
      if(trk[tt][0]) mPToverET[tt][cut]->Fill(etpt[tt][0]);
      if(trk[tt][1]) mPToverET[tt][cut]->Fill(etpt[tt][1]);
      if(trk[tt][0] && trk[tt][1]) mChargeSum[tt][cut]->Fill(cg[tt][0] + cg[tt][1]);

      mET  [tt][cut]->Fill(ET);
      mEZ  [tt][cut]->Fill(EZ);
      mM   [tt][cut]->Fill(M);
      mZ   [tt][cut]->Fill(Z);
      mCosT[tt][cut]->Fill(CosT);
      mPhi [tt][cut]->Fill(Phi);

      //Mixed-event background: at the chosen cut level, pair this event's north/south
      //FCS cluster candidates with the previous event's (opposite side), then buffer
      //this event's candidates for the next one.
      if(cut==mMixCut){
        if(mPrevValid[tt]){
          mMmix[tt]->Fill( (lnLab + mPrevLS[tt]).M() ); //this-North + previous-South
          mMmix[tt]->Fill( (mPrevLN[tt] + lsLab).M() ); //previous-North + this-South
          nMixFilled[tt] += 2;
        }
        mPrevLN[tt] = lnLab;
        mPrevLS[tt] = lsLab;
        mPrevValid[tt] = true;
      }

      mET12[tt][cut]->Fill(ETN,ETS);
      mXFPT[tt][cut]->Fill(EN/255.0,ETN);
      mXFPT[tt][cut]->Fill(ES/255.0,ETS);
      mXY[tt][cut]->Fill(FcsXyz[0].x(),FcsXyz[0].y());
      mXY[tt][cut]->Fill(FcsXyz[1].x(),FcsXyz[1].y());
      if(trk[tt][0]) mPTET[tt][cut]->Fill(ETN,trkpt[tt][0]);
      if(trk[tt][1]) mPTET[tt][cut]->Fill(ETS,trkpt[tt][1]);

      if(trk[tt][0]) mZVTX[tt][cut]->Fill(trk[tt][0]->dcaZ());
      if(trk[tt][1]) mZVTX[tt][cut]->Fill(trk[tt][1]->dcaZ());
      if(trk[tt][0] && trk[tt][1]){
        mZVTXA[tt][cut]->Fill((trk[tt][0]->dcaZ()+trk[tt][1]->dcaZ())/2.0);
        mZVTXD[tt][cut]->Fill(trk[tt][0]->dcaZ()-trk[tt][1]->dcaZ());
      }
    }
  }
}

void EndDilepton(){
  printf("=== Dilepton Event Statistics ===\n");
  printf("Events processed              : %ld\n", nEvt);
  printf("Events with both FCS clusters : %ld\n", nEvtBothClus);
  for(int tt=0; tt<mNType; tt++){
    printf("--- %s ---\n",TTYPE[tt]);
    for(int cut=0; cut<mNCut; cut++){
      printf("  Cut[%d] %-12s : %ld\n", cut, nameCut[cut], nEvtCut[tt][cut]);
    }
    printf("  PassDR  ns=0:%ld  ns=1:%ld\n", nTrkPassDR[tt][0], nTrkPassDR[tt][1]);
    printf("  Mixed-event entries filled (cut=%s) : %ld\n", nameCut[mMixCut], nMixFilled[tt]);
    printf("  BestPT  ns=0:%ld  ns=1:%ld\n", nTrkPassPT[tt][0], nTrkPassPT[tt][1]);
  }
  printf("=== Track Selection Statistics ===\n");
  printf("Total tracks seen             : %ld\n", nTrkAll);
  printf("Survived status>=2            : %ld\n", nTrkAll-nTrkFailStatus);
  printf("Survived nHit>=%d             : %ld\n", mDNHitCut, nTrkAll-nTrkFailStatus-nTrkFailNHit);
  printf("Survived chi2<=%.0f           : %ld\n", mDChi2Cut, nTrkAll-nTrkFailStatus-nTrkFailNHit-nTrkFailChi2);
  printf("Writing and closing %s\n",filenameD);
  mFileD->Write();
  mFileD->Close();
}
