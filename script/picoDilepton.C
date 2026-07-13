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

// Pi0 (gamma-gamma) mass reconstruction -- self-contained addition, see
// RunPi0() at the end of this file. Forward-declared here so RunDilepton()
// (defined earlier in the file) can call it.
void RunPi0(StPicoDst *dst, StFcsDb* fcsDb, double zVertex);

enum {mNCut=7};
const char* nameCut[mNCut] = {"All","ETOT","HTOT","Cone","SigmaMax","TrackMatch","ChargeSign"};
static const int mNType=6;
const char* TTYPE[mNType]={"Global","Beamline","Primary","FwdVtx","BLCVtx","FCSTRK"};
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
int mDebug=0; // per-track/per-event printf is extremely verbose (multi-GB logs on
              // real data with condor's Log=Output=Error all pointed at the same
              // file); set to 1 only for interactive small-sample debugging.
static long nEvt=0, nEvtBothClus=0;
static long nEvtCut[mNType][7];
static long nTrkAll=0, nTrkFailStatus=0, nTrkFailNHit=0, nTrkFailChi2=0;
static long nTrkPassDR[mNType][2];
static long nTrkPassPT[mNType][2];

//Mixed-event background: cut level at which an event's north/south FCS cluster
//4-vectors are buffered and mixed with a pool of the last mMixPoolSize qualifying
//events, per track type. Adjustable; default matches the tightest existing cut
//(ChargeSign). Pool gives up to 2*mMixPoolSize mixed pairs per event (vs. 2 with
//a single previous-event buffer), reducing the mixed-bg statistical error by
//roughly sqrt(mMixPoolSize) relative to the old single-previous-event scheme.
//
//The pool is further split into mNZBin(z-vertex) x mNChBin(north-charge-sign)
//bins, and mixing only draws from the bin matching the current event's own
//zVertex and north-lepton charge sign. This fixes two real gaps in the
//original flat pool: (1) mixing pairs used to combine 4-vectors reconstructed
//with each event's own (possibly very different) zVertex, which is not a
//physically consistent random pairing given the wide (~30-40cm sigma) z-vertex
//spread in real data; (2) with only same-event pairs required to be
//opposite-sign (cg[0]+cg[1]==0), the pool had no charge bookkeeping, so a
//"north+" from this event could get mixed with a "south+" pulled from a past
//event that itself was (north-,south+) -- silently contaminating the
//OS-labeled mixed background with SS-equivalent pairs. Binning by north
//charge sign and only mixing within the same bin ("current north charge
//matches past north charge" -- south then matches automatically since both
//sides passed the same-event OS cut) closes that gap.
//nMixFilledBin/nSameFilledBin let plotDilep.C reweight+sum the 12 bins
//correctly (mixed-pair fill count vs. real same-event pair count actually
//landing in each bin) instead of assuming a uniform flat normalization.
int mMixCut=6;
static const int mMixPoolSize=10; // number of past qualifying events kept per (type,zbin,chargebin)
static const int mNZBin=6;                 // z-vertex bins, 25cm steps
static const float mZBinLo=-50.0;          // bin edges: [-50,-25),...,[75,100)
static const float mZBinWidth=25.0;
static const int mNChBin=2;                // 0: north charge>=0, 1: north charge<0
static int mPoolNext[mNType][mNZBin][mNChBin];                        // next slot to overwrite (circular)
static bool mPoolValid[mNType][mNZBin][mNChBin][mMixPoolSize];
static TLorentzVector mPoolLN[mNType][mNZBin][mNChBin][mMixPoolSize];
static TLorentzVector mPoolLS[mNType][mNZBin][mNChBin][mMixPoolSize];
static long nMixFilled[mNType];                          // total mixed-pair fills (all bins), for EndDilepton summary
static long nMixFilledBin[mNType][mNZBin][mNChBin];       // mixed-pair fills per bin
static long nSameFilledBin[mNType][mNZBin][mNChBin];      // same-event pairs landing in each bin

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
TH1F *mMmixBin[mNType][mNZBin][mNChBin];  // mixed-mass, per (type,zbin,chargebin)
TH1F *mMsameBin[mNType][mNZBin][mNChBin]; // same-event (OS) mass, per (type,zbin,chargebin)
// Same-event LIKE-sign mass, per (type,zbin,chargebin) -- chargebin here means
// north-lepton charge sign of the LS pair (++ -> bin0, -- -> bin1), same key
// as the OS pool. Real dileptons are essentially all opposite-sign, so LS
// pairs are a signal-free measure of the combinatorial rate and make a better
// mixed-event normalization anchor than the OS same-event count (which
// includes true signal, biasing a full-range-integral normalization -- see
// plotDilep.C).
TH1F *mMLSameBin[mNType][mNZBin][mNChBin];
TH1F *mBLCVtxZ[mNCut]; // event-level BLC vertex Z, filled per cut (BLCVtx track type only)
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

  for(int cut=0; cut<mNCut; cut++){
    mBLCVtxZ[cut] = new TH1F(Form("BLCVtxZ_%s",nameCut[cut]),Form("BLC vertex Z (event-level, cut=%s)",nameCut[cut]),100,-200,200);
  }

  for(int tt=0; tt<mNType; tt++){
    nMixFilled[tt] = 0;
    for(int zb=0; zb<mNZBin; zb++){
      for(int cb=0; cb<mNChBin; cb++){
        mPoolNext[tt][zb][cb] = 0;
        for(int k=0; k<mMixPoolSize; k++) mPoolValid[tt][zb][cb][k] = false;
        nMixFilledBin[tt][zb][cb] = 0;
        nSameFilledBin[tt][zb][cb] = 0;
        float zlo=mZBinLo+zb*mZBinWidth, zhi=zlo+mZBinWidth;
        const char* cname = (cb==0) ? "Np" : "Nm";
        mMmixBin[tt][zb][cb] = new TH1F(Form("MmixBin_%s_z%d_%s",TTYPE[tt],zb,cname),
          Form("Mixed-Event Mass %s z=[%.0f,%.0f) north-chg=%s (cut=%s)",TTYPE[tt],zlo,zhi,cname,nameCut[mMixCut]),50,0.0,10.0);
        mMsameBin[tt][zb][cb] = new TH1F(Form("MsameBin_%s_z%d_%s",TTYPE[tt],zb,cname),
          Form("Same-Event Mass %s z=[%.0f,%.0f) north-chg=%s (cut=%s)",TTYPE[tt],zlo,zhi,cname,nameCut[mMixCut]),50,0.0,10.0);
        mMLSameBin[tt][zb][cb] = new TH1F(Form("MLSameBin_%s_z%d_%s",TTYPE[tt],zb,cname),
          Form("Same-Event Like-Sign Mass %s z=[%.0f,%.0f) north-chg=%s (cuts 0-5)",TTYPE[tt],zlo,zhi,cname),50,0.0,10.0);
      }
    }
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

  // This picoDst was written with setVtxMode(StPicoDstMaker::Vtxless), so
  // event->primaryVertex() is never filled and always reads back as the
  // (-999,-999,-999) sentinel, not a real vertex. Using -999 as zVertex in
  // getLorentzVector() below nearly doubles the assumed cluster-to-vertex z
  // distance (real ~710cm -> ~1709cm), which collimates every reconstructed
  // direction toward pure z: ET (~sin theta) collapses while EZ (~E for a
  // forward-peaked direction) barely changes -- exactly the symptom seen.
  // MC (Pythia) events are generated at the nominal z=0 vertex, so substitute
  // 0 for the sentinel here instead of a real measured vertex.
  double zVertexPrimary = event->primaryVertex().Z();
  if (zVertexPrimary < -900) zVertexPrimary = 0;

  // BLCVtx (tt==4) and FCSTRK (tt==5, BLCVtx tracks refit with an FCS ECAL
  // cluster constraint) are both downstream of the BLC (beamline-constrained)
  // vertex fit, so their FCS cluster kinematics should use that fitted vertex
  // instead of the general primary vertex -- it's a real per-event fit
  // result (unaffected by the Vtxless sentinel above) and is the physically
  // consistent vertex for tracks anchored to it. When the BLC fit didn't
  // converge this event (haveBLCVtx=false), those two types are skipped
  // entirely below (kinValid[1]=false).
  bool haveBLCVtx = (event->blcVertexNTracks() > 0);
  double zVertexBLC = haveBLCVtx ? event->blcVertex().Z() : 0;

  //loop over FCS clusters and find higest ET Ecal clusters (selection uses
  //the primary-vertex zVertex uniformly across types -- this only decides
  //*which* 2 clusters are candidates, a choice only weakly sensitive to a
  //vertex-position assumption at the cm level)
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
    StLorentzVectorD lv=fcsDb->getLorentzVector(xyz,c->energy(),zVertexPrimary);
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

  //Getting TOT & Cone for isolation cut (vertex-independent: raw hit energies
  //in an eta/phi cone around the cluster's own position)
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

  //2-body decay kinematics, computed once per vertex source: v=0 uses
  //zVertexPrimary (track types Global/Beamline/Primary/FwdVtx), v=1 uses
  //zVertexBLC (BLCVtx/FCSTRK). kinValid[1]=false (BLC fit unavailable) means
  //v=1's arrays are never read -- every downstream use is guarded by the
  //per-track-type skip in the main cut loop below.
  TLorentzVector lnLab[2], lsLab[2];
  double EN[2],ES[2],E[2],ETN[2],ETS[2],ET[2],EZ[2],M[2],Z[2],Phi[2],CosT[2];
  bool kinValid[2] = {true, haveBLCVtx};
  double zv[2] = {zVertexPrimary, zVertexBLC};
  for(int v=0; v<2; v++){
    if(!kinValid[v]) continue;
    StLorentzVectorD stln=fcsDb->getLorentzVector(FcsXyz[0],highest[0]->energy(),zv[v]);
    StLorentzVectorD stls=fcsDb->getLorentzVector(FcsXyz[1],highest[1]->energy(),zv[v]);
    TLorentzVector ln(stln.px(),stln.py(),stln.pz(),stln.e());
    TLorentzVector ls(stls.px(),stls.py(),stls.pz(),stls.e());
    TLorentzVector di = ln + ls;
    lnLab[v]=ln; lsLab[v]=ls; //lab-frame copies for mixed-event background (ln/ls get boosted to the pair CM frame below)
    EN[v]=ln.E(); ES[v]=ls.E(); E[v]=di.E();
    ETN[v]=ln.Perp(); ETS[v]=ls.Perp();
    ET[v]=di.Perp(); EZ[v]=di.Pz(); M[v]=di.M();
    Z[v]=abs(EN[v]-ES[v])/(EN[v]+ES[v]);
    Phi[v]=di.Phi();
    TVector3 boost=-di.BoostVector();
    ln.Boost(boost);
    ls.Boost(boost);
    CosT[v]=ln.CosTheta(); //take north one for now... When we have tracking, take positive charged
  }
  if(mDebug>0){
    printf("FCS VTX= %8.3f  %8.3f  %8.3f  BLCVTX=%8.3f (nTrk=%d)\n",
           event->primaryVertex().X(),event->primaryVertex().Y(),event->primaryVertex().Z(),
           event->blcVertex().Z(),event->blcVertexNTracks());
    printf("FCS EN %8.3f %8.3f %8.3f E=%8.3f ET=%8.3f Phi=%8.3f\n",
	   FcsXyz[0].x(),FcsXyz[0].y(),FcsXyz[0].z(),highest[0]->energy(),FcsXyz[0].perp(),FcsXyz[0].phi());
    printf("FCS ES %8.3f %8.3f %8.3f E=%8.3f ET=%8.3f Phi=%8.3f\n",
	   FcsXyz[1].x(),FcsXyz[1].y(),FcsXyz[1].z(),highest[1]->energy(),FcsXyz[1].perp(),FcsXyz[1].phi());
  }

  //Ecal cluster SigmaMax (vertex-independent, cluster shape only)
  double SigmaMaxN = highest[0]->sigmaMax();
  double SigmaMaxS = highest[1]->sigmaMax();

  //Ratio of DiLepton candidate to TOT/Cone, per vertex source
  double ratioETOT[2]={0,0}, ratioHTOT[2]={9.99,9.99}, ratioConeN[2]={0,0}, ratioConeS[2]={0,0};
  for(int v=0; v<2; v++){
    if(!kinValid[v]) continue;
    ratioETOT[v] = E[v]/tot[0];
    if(tot[1]>0) ratioHTOT[v]=E[v]/tot[1];
    ratioConeN[v] = EN[v]/cone[0];
    ratioConeS[v] = ES[v]/cone[1];
  }

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
        trk[tt][ns]=t; trkpt[tt][ns]=pt; cg[tt][ns]=t->charge();
        dr[tt][ns]=dR; rdphi[tt][ns]=Rdphi;
      }
    }
  }
  //etpt depends on which vertex source this track type uses (ETN/ETS[v]) --
  //fill after the track loop once trkpt[][] is finalized.
  for(int tt=0; tt<mNType; tt++){
    int v = (tt>=4) ? 1 : 0;
    if(!kinValid[v]) continue;
    if(trk[tt][0]) etpt[tt][0] = ETN[v]/trkpt[tt][0];
    if(trk[tt][1]) etpt[tt][1] = ETS[v]/trkpt[tt][1];
  }
  if(mDebug>0){
    for(int tt=0; tt<mNType; tt++){
      for(int ns=0; ns<2; ns++){
        if(trk[tt][ns]){
	  printf("trk NS=%1d TT=%1d(%s) ET/pT=%6.4f trkid=%2d pt=%12.2f cg=%2d DCA=%7.3f %7.1f dR=%7.2f Rdphi=%7.2f\n",
	         ns,tt,TTYPE[tt],etpt[tt][ns],trk[tt][ns]->id(),trkpt[tt][ns],cg[tt][ns],
                 trk[tt][ns]->dcaXY(),trk[tt][ns]->dcaZ(),dr[tt][ns],rdphi[tt][ns]);
        }
      }
    }
  }

  //Apply cuts and fill histograms for each track type
  for(int tt=0; tt<mNType; tt++){
    int v = (tt>=4) ? 1 : 0;
    if(!kinValid[v]) continue; //no BLC vertex fit this event -- BLCVtx/FCSTRK skipped entirely
    for(int cut=0; cut<mNCut; cut++){
      if(cut==1 && ratioETOT[v]<mETotCut) break;
      if(cut==2 && ratioHTOT[v]<mHTotCut) break;
      if(cut==3 && (ratioConeN[v]<mConeCut || ratioConeS[v]<mConeCut)) break;
      if(cut==4 && (SigmaMaxN > mSigmaMaxCut || SigmaMaxS > mSigmaMaxCut) ) break;
      if(cut==5 && (etpt[tt][0] < mETPTCutLow || etpt[tt][1] < mETPTCutLow || etpt[tt][0] > mETPTCutHigh || etpt[tt][1] > mETPTCutHigh )) break;
      if(cut==6 && cg[tt][0] + cg[tt][1] != 0){
        //Like-sign pair (both cuts 0-5 passed): record as a signal-free
        //combinatorial-rate reference for plotDilep.C's mixed-event
        //normalization, binned the same way as the OS mixing pool.
        double zVtxForBin = (v==1) ? zVertexBLC : zVertexPrimary;
        int zbin = (int)floor((zVtxForBin - mZBinLo)/mZBinWidth);
        int chbin = (cg[tt][0] >= 0) ? 0 : 1;
        if(zbin>=0 && zbin<mNZBin) mMLSameBin[tt][zbin][chbin]->Fill(M[v]);
        break;
      }
      nEvtCut[tt][cut]++;

      mETot[tt][cut]->Fill(ratioETOT[v]);
      mHTot[tt][cut]->Fill(ratioHTOT[v]);
      mCone[tt][cut]->Fill(ratioConeN[v]);
      mCone[tt][cut]->Fill(ratioConeS[v]);
      mSigmax[tt][cut]->Fill(SigmaMaxN);
      mSigmax[tt][cut]->Fill(SigmaMaxS);
      if(trk[tt][0]) mPToverET[tt][cut]->Fill(etpt[tt][0]);
      if(trk[tt][1]) mPToverET[tt][cut]->Fill(etpt[tt][1]);
      if(trk[tt][0] && trk[tt][1]) mChargeSum[tt][cut]->Fill(cg[tt][0] + cg[tt][1]);

      mET  [tt][cut]->Fill(ET[v]);
      mEZ  [tt][cut]->Fill(EZ[v]);
      mM   [tt][cut]->Fill(M[v]);
      mZ   [tt][cut]->Fill(Z[v]);
      mCosT[tt][cut]->Fill(CosT[v]);
      mPhi [tt][cut]->Fill(Phi[v]);

      //Mixed-event background: at the chosen cut level, pair this event's north/south
      //FCS cluster candidates with every qualifying event currently in the pool
      //*within the same (zVertex, north-charge-sign) bin* (opposite side), then
      //buffer this event's candidates into that bin's pool for future events --
      //always done AFTER mixing, so an event is never mixed with itself. Events
      //with |zVertex|>=150 (outside the binned range) are skipped entirely for
      //mixing purposes.
      if(cut==mMixCut){
        double zVtxForBin = (v==1) ? zVertexBLC : zVertexPrimary;
        int zbin = (int)floor((zVtxForBin - mZBinLo)/mZBinWidth);
        int chbin = (cg[tt][0] >= 0) ? 0 : 1;
        if(zbin>=0 && zbin<mNZBin){
          for(int k=0; k<mMixPoolSize; k++){
            if(!mPoolValid[tt][zbin][chbin][k]) continue;
            double m1 = (lnLab[v] + mPoolLS[tt][zbin][chbin][k]).M(); //this-North + pool-South
            double m2 = (mPoolLN[tt][zbin][chbin][k] + lsLab[v]).M(); //pool-North + this-South
            mMmix[tt]->Fill(m1);
            mMmix[tt]->Fill(m2);
            mMmixBin[tt][zbin][chbin]->Fill(m1);
            mMmixBin[tt][zbin][chbin]->Fill(m2);
            nMixFilled[tt] += 2;
            nMixFilledBin[tt][zbin][chbin] += 2;
          }
          mMsameBin[tt][zbin][chbin]->Fill(M[v]);
          nSameFilledBin[tt][zbin][chbin]++;
          int slot = mPoolNext[tt][zbin][chbin];
          mPoolLN[tt][zbin][chbin][slot] = lnLab[v];
          mPoolLS[tt][zbin][chbin][slot] = lsLab[v];
          mPoolValid[tt][zbin][chbin][slot] = true;
          mPoolNext[tt][zbin][chbin] = (slot + 1) % mMixPoolSize;
        }
      }

      mET12[tt][cut]->Fill(ETN[v],ETS[v]);
      mXFPT[tt][cut]->Fill(EN[v]/255.0,ETN[v]);
      mXFPT[tt][cut]->Fill(ES[v]/255.0,ETS[v]);
      mXY[tt][cut]->Fill(FcsXyz[0].x(),FcsXyz[0].y());
      mXY[tt][cut]->Fill(FcsXyz[1].x(),FcsXyz[1].y());
      if(trk[tt][0]) mPTET[tt][cut]->Fill(ETN[v],trkpt[tt][0]);
      if(trk[tt][1]) mPTET[tt][cut]->Fill(ETS[v],trkpt[tt][1]);

      if(trk[tt][0]) mZVTX[tt][cut]->Fill(trk[tt][0]->dcaZ());
      if(trk[tt][1]) mZVTX[tt][cut]->Fill(trk[tt][1]->dcaZ());
      if(trk[tt][0] && trk[tt][1]){
        mZVTXA[tt][cut]->Fill((trk[tt][0]->dcaZ()+trk[tt][1]->dcaZ())/2.0);
        mZVTXD[tt][cut]->Fill(trk[tt][0]->dcaZ()-trk[tt][1]->dcaZ());
      }

      // Event-level BLC vertex Z (not track-type dependent, but only makes
      // sense to fill once per event -- piggyback on the BLCVtx (tt==4) pass
      // through this cut chain, same cut boundaries as everything else here).
      if(tt==4 && haveBLCVtx) mBLCVtxZ[cut]->Fill(event->blcVertex().Z());
    }
  }

  // Independent pi0 (gamma-gamma) mass reconstruction -- not tied to any
  // track type or cut chain above, so just called once per event here.
  RunPi0(dst, fcsDb, zVertexPrimary);
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
    for(int zb=0; zb<mNZBin; zb++){
      for(int cb=0; cb<mNChBin; cb++){
        printf("    zbin=[%4.0f,%4.0f) north-chg=%s : same=%-6ld mixed=%-6ld\n",
               mZBinLo+zb*mZBinWidth, mZBinLo+(zb+1)*mZBinWidth, (cb==0)?"+":"-",
               nSameFilledBin[tt][zb][cb], nMixFilledBin[tt][zb][cb]);
      }
    }
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

//=====================================================================
// Pi0 (gamma-gamma) mass reconstruction from FCS ECAL clusters.
// Simple inclusive combinatorial pairing: every ECAL cluster pair in
// the event (both N/S sides, no track match needed -- photons don't
// leave tracks), each cluster required only to pass a minimum energy
// cut. No isolation/shower-shape cuts, unlike the electron-candidate
// selection above -- deliberately loose to keep low-energy pi0 decays
// in the sample. Real pi0 signal shows up as a peak on top of
// combinatorial background; no background subtraction done here.
//=====================================================================
float mPi0MinClusterE = 0.5; //GeV, minimum FCS ECAL cluster energy for pi0 pairing
TH1F *mPi0Mass = 0; //lazily created on first call below

void RunPi0(StPicoDst *dst, StFcsDb* fcsDb, double zVertex){
  if(!mPi0Mass){
    mPi0Mass = new TH1F("Pi0Mass",
      "FCS ECAL cluster-pair mass (all pairs, E>0.5 GeV each, no N/S or track requirement);M_{#gamma#gamma} [GeV];Counts",
      100,0.0,1.0);
  }

  int nClusters = dst->numberOfFcsClusters();
  for(int i=0; i<nClusters; i++){
    StPicoFcsCluster* ci = dst->fcsCluster(i);
    if(fcsDb->ecalHcalPres(ci->detectorId())!=0) continue; //only Ecal
    if(ci->energy() < mPi0MinClusterE) continue;
    StThreeVectorD xyzi = fcsDb->getStarXYZfromColumnRow(ci->detectorId(),ci->x(),ci->y());
    StLorentzVectorD stlvi = fcsDb->getLorentzVector(xyzi,ci->energy(),zVertex);
    TLorentzVector lvi(stlvi.px(),stlvi.py(),stlvi.pz(),stlvi.e());
    for(int j=i+1; j<nClusters; j++){
      StPicoFcsCluster* cj = dst->fcsCluster(j);
      if(fcsDb->ecalHcalPres(cj->detectorId())!=0) continue; //only Ecal
      if(cj->energy() < mPi0MinClusterE) continue;
      StThreeVectorD xyzj = fcsDb->getStarXYZfromColumnRow(cj->detectorId(),cj->x(),cj->y());
      StLorentzVectorD stlvj = fcsDb->getLorentzVector(xyzj,cj->energy(),zVertex);
      TLorentzVector lvj(stlvj.px(),stlvj.py(),stlvj.pz(),stlvj.e());
      mPi0Mass->Fill((lvi+lvj).M());
    }
  }
}
