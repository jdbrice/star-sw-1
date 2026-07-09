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
#include "StPicoEvent/StPicoMcTrack.h"
#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

static const int mDebug=0;
static const int mNCut=3;
static const int mMaxTrack=200;
//int mTrackType=2;
int mTrackType=1;
double mChi2Cut=100000;
double mNHitCut=4;
//double mEcalXminCut=35;
double mEcalXminCut=15;
double mEcalXmaxCut=150;
//double mHcalXminCut=40;
double mHcalXminCut=15;
double mHcalXmaxCut=150;
double mYmaxCut=100;
double mDCut=25;        //dy cut for dx, dx cut for dy
int mMaxEvent=20000000;

double PI=TMath::Pi();
double TWOPI=2.0*PI;
double PIO2=PI/2.0;
double THREEPIO2=3.0*PI/2.0;
double ANG(double phi){
  while(phi < -PIO2) phi+=TWOPI;
  while(phi >= THREEPIO2) phi-=TWOPI;
  return phi;
}
double DANG(double phi){
  while(phi < -PI) phi+=TWOPI;
  while(phi >= PI) phi-=TWOPI;
  return phi;
}
double R(double x, double y){return sqrt(x*x + y*y);}
double PHI(double x, double y){return ANG(atan2(y,x));}
double DPHI(double phi1, double phi2){return DANG(phi1-phi2);}

char filenameM[200];
TFile *mFileM;

const char *EH[2]={"Ecal","Hcal"};
const char *NSTB[4][2]={{"North","South"},{"Top","Bottom"},{"North","South"},{"R<70","R>70"}};
const char *CUT[4][mNCut+1]={{"Same Event","Mixed Event","Same-Mixed North","Same-Mixed South"},
			     {"Same Event","Mixed Event","Same-Mixed Top","Same-Mixed Bottom"},
			     {"Same Event","Mixed Event","Same-Mixed North","Same-Mixed South"},
			     {"Same Event","Mixed Event","Same-Mixed R<70","Same-Mixed R>70"}};
TH1F *hBunchId; 
TH2F *hXY[4];
TH1F *hE[3];
TH1F *hET[3];
TH1F *hdx[2][2][3], *hdy[2][2][3], *hdp[2][2][3], *hdr[2][2][3];
TH2F *hxdx[2][3], *hydx[2][3], *hpdx[2][3], *hrdx[2][3];
TH2F *hxdy[2][3], *hydy[2][3], *hpdy[2][3], *hrdy[2][3];
TH2F *hxdp[2][3], *hydp[2][3], *hpdp[2][3], *hrdp[2][3];
TH2F *hxdr[2][3], *hydr[2][3], *hpdr[2][3], *hrdr[2][3];

static int ntrk[2]={0,0};
static double tx[2][2][mMaxTrack],ty[2][2][mMaxTrack];

// Track-quality histograms by track type (no tracktype cut; basic quality only)
// Types: 0=Global 1=Beamline 2=Primary 3=FwdVtx 4=BLCVtx
// BBB: BLCVtx = BLC tracks refitted to the beam-line-constrained forward vertex
static const int mNType=5;
const char *TTYPE[mNType]={"Global","Beamline","Primary","FwdVtx","BLCVtx"};
TH1F *hTrkQPt[mNType];       // q/pT = signed curvature (charge sign + momentum)
TH1F *hTrkEta[mNType];      // pseudorapidity (forward acceptance)
TH1F *hTrkPhi[mNType];      // azimuthal angle (acceptance gaps)
TH1F *hTrkDcaZ[mNType];     // DCA z-position (vertex quality)
// MC truth-matching histograms (filled only when McTracks available)
TH1F *hTrkDqPt[mNType];     // q/pT_reco - q/pT_true  (momentum bias + resolution)
TH1F *hTrkDpTrel[mNType];   // (pT_reco - pT_true)/pT_true  (relative pT residual)
TH2F *hTrkQPtRecovsTru[mNType]; // q/pT_reco vs q/pT_true (2D correlation)

// BLCVtx vertex-quality histograms (per-event, using BLC and BLCVtx tracks)
TH1F *hBLCVtxN;        // N good BLC tracks per event (used as vertex inputs)
TH1F *hBLCVtxZ;        // per-event avg BLC DCA-z = estimated z_vtx
TH1F *hBLCVtxZres;     // estimated z_vtx - event primary z (= MC truth vz in simulation)
TH1F *hBLCVtxZdiff;    // BLC dcaZ_i - dcaZ_j for each pair (vertex resolution proxy)
TH1F *hBLCVtxTrkDcaZ;  // DCA-z of BLCVtx tracks after refit to vertex (should be ~0)

void InitMatch(int run, int set=-1) {
  if(set<0) sprintf(filenameM,"hist_match/%d.match.root",run);
  else       sprintf(filenameM,"hist_match/%d.%d.match.root",run,set);
  printf("Opening %s\n",filenameM);
  mFileM=new TFile(filenameM,"RECREATE");

  hBunchId = new TH1F("bunchId","bunchId",120,0.0,120.0);
  hXY[0] = new TH2F("xyEcal","Ecal; x; y",100,-150,150,100,-100,100);
  hXY[1] = new TH2F("xyHcal","Hcal; x; y",100,-150,150,100,-100,100);
  hXY[2] = new TH2F("xyETrk","Trk Ecal Projection; x; y",100,-150,150,100,-100,100);
  hXY[3] = new TH2F("xyHTrk","Trk Hcal Projection; x; y",100,-150,150,100,-100,100);
  for(int eh=0; eh<2; eh++){
    hE[eh]=new TH1F(Form("E%s",EH[eh]),Form("E %s; E[GeV]",EH[eh]),100,0.0,40.0);
    hET[eh]=new TH1F(Form("ET%s",EH[eh]),Form("ET %s; ET[GeV]",EH[eh]),100,0.0,5.0);
    for(int cut=0; cut<mNCut; cut++){
      hxdx[eh][cut]=new TH2F(Form("xdx%s%s",EH[eh],CUT[0][cut]),Form("dX(fcsX) %s; fcsX[cm]; dx[cm]",      EH[eh]),100,-150,150,100,-100,100);
      hydx[eh][cut]=new TH2F(Form("ydx%s%s",EH[eh],CUT[0][cut]),Form("dX(fcsY) %s; fcsY[cm]; dx[cm]",      EH[eh]),100,-100,100,100,-100,100);
      hpdx[eh][cut]=new TH2F(Form("pdx%s%s",EH[eh],CUT[0][cut]),Form("dX(fcsR*fcsPhi) %s; fcsR*Phi[cd]; dx[cm]",EH[eh]),100,-150,150,100,-100,100);
      hrdx[eh][cut]=new TH2F(Form("rdx%s%s",EH[eh],CUT[0][cut]),Form("dX(fcsR) %s; fcsR[cm]; dx[cm]",      EH[eh]),100,   0,180,100,-100,100);
      hxdy[eh][cut]=new TH2F(Form("xdy%s%s",EH[eh],CUT[1][cut]),Form("dY(fcsX) %s; fcsX[cm]; dy[cm]",      EH[eh]),100,-150,150,100,-100,100);
      hydy[eh][cut]=new TH2F(Form("ydy%s%s",EH[eh],CUT[1][cut]),Form("dY(fcsY) %s; fcsY[cm]; dy[cm]",      EH[eh]),100,-100,100,100,-100,100);
      hpdy[eh][cut]=new TH2F(Form("pdy%s%s",EH[eh],CUT[1][cut]),Form("dY(fcsR*fcsPhi) %s; fcsR*Phi[cm]; dy[cm]",EH[eh]),100,-150,150,100,-100,100);
      hrdy[eh][cut]=new TH2F(Form("rdy%s%s",EH[eh],CUT[1][cut]),Form("dY(fcsR) %s; fcsR[cm]; dy[cm]",      EH[eh]),100,   0,180,100,-100,100);
      hxdp[eh][cut]=new TH2F(Form("xdp%s%s",EH[eh],CUT[2][cut]),Form("fcsR*dPhi(fcsX) %s; fcsX[cm]; fcsr*rdphi[cm]",     EH[eh]),100,-150,150,100,-100,100);
      hydp[eh][cut]=new TH2F(Form("ydp%s%s",EH[eh],CUT[2][cut]),Form("fcsR*dPhi(fcsY) %s; fcsY[cm]; fcsr*dphi[cm]",      EH[eh]),100,-100,100,100,-100,100);
      hpdp[eh][cut]=new TH2F(Form("pdp%s%s",EH[eh],CUT[2][cut]),Form("fcsR*dPhi(fcsR*fcsPhi) %s; fcsR*Phi[cm]; fcsr*dphi[cm]",EH[eh]),100,-150,150,100,-100,100);
      hrdp[eh][cut]=new TH2F(Form("rdp%s%s",EH[eh],CUT[2][cut]),Form("fcsR*dPhi(fcsR) %s; fcsR[cm]; fcsr*dphi[cm]",      EH[eh]),100,   0,180,100,-100,100);
      hxdr[eh][cut]=new TH2F(Form("xdr%s%s",EH[eh],CUT[3][cut]),Form("dR(fcsX) %s; fcsX[cm]; dr[cm]",      EH[eh]),100,-150,150,100,-100,100);
      hydr[eh][cut]=new TH2F(Form("ydr%s%s",EH[eh],CUT[3][cut]),Form("dR(fcsY) %s; fcsY[cm]; dr[cm]",      EH[eh]),100,-100,100,100,-100,100);
      hpdr[eh][cut]=new TH2F(Form("pdr%s%s",EH[eh],CUT[3][cut]),Form("dR(fcsR*fcsPhi) %s; fcsR*Phi[cm]; dr[cm]",EH[eh]),100,-150,150,100,-100,100);
      hrdr[eh][cut]=new TH2F(Form("rdr%s%s",EH[eh],CUT[3][cut]),Form("dR(fcdR) %s; fcsR[cm]; dr[cm]",      EH[eh]),100,   0,180,100,-100,100);
      for(int nstb=0; nstb<2; nstb++){
	printf("A %s %s %s\n",EH[eh],NSTB[0][nstb],CUT[0][cut]);
	hdx[eh][nstb][cut] = new TH1F(Form("%s%sdx%s",EH[eh],NSTB[0][nstb],CUT[0][cut]),
				      Form("%s%s-Trk dX %s; dX[cm]",EH[eh],NSTB[0][nstb],CUT[0][cut]),100,-200.0,200.0);
	printf("B %s %s %s\n",EH[eh],NSTB[1][nstb],CUT[1][cut]);
	hdy[eh][nstb][cut] = new TH1F(Form("%s%sdy%s",EH[eh],NSTB[1][nstb],CUT[1][cut]),
				      Form("%s%s-Trk dY %s; dY[cm]",EH[eh],NSTB[1][nstb],CUT[1][cut]),100,-200.0,200.0);      

	printf("C %s %s %s\n",EH[eh],NSTB[2][nstb],CUT[2][cut]);
	hdp[eh][nstb][cut] = new TH1F(Form("%s%sdp%s",EH[eh],NSTB[2][nstb],CUT[2][cut]),
				      Form("%s%s-Trk R*dPhi %s; fcsR*dPhi[cm]",EH[eh],NSTB[2][nstb],CUT[2][cut]),100,-200.0,200.0);

	printf("D %s %s %s\n",EH[eh],NSTB[3][nstb],CUT[3][cut]);
	hdr[eh][nstb][cut] = new TH1F(Form("%s%sdr%s",EH[eh],NSTB[3][nstb],CUT[3][cut]),
				      Form("%s%s-Trk dR %s; dR[cm]",EH[eh],NSTB[3][nstb],CUT[3][cut]),100,-200.0,200.0);      	
      }
    }
  }
  hE[2]=new TH1F(Form("ETrk"),Form("E Trk; p[GeV]"),100,0.0,40.0);
  hET[2]=new TH1F(Form("PTTrk"),Form("PT Trk; pT[GeV]"),100,0.0,5.0);

  // Track-quality histograms for all track types (basic quality cuts, no tracktype filter)
  for(int tt=0; tt<mNType; tt++){
    hTrkQPt[tt]    = new TH1F(Form("trkQPt%s",   TTYPE[tt]),
        Form("q/pT %s; q/pT [GeV/c]^{-1}",TTYPE[tt]),100,-5.0,5.0);
    hTrkEta[tt]    = new TH1F(Form("trkEta%s",   TTYPE[tt]),
        Form("#eta %s; #eta",            TTYPE[tt]),100, 1.5, 6.0);
    hTrkPhi[tt]    = new TH1F(Form("trkPhi%s",   TTYPE[tt]),
        Form("#phi %s; #phi [rad]",      TTYPE[tt]),100,-3.2, 3.2);
    hTrkDcaZ[tt]   = new TH1F(Form("trkDcaZ%s",  TTYPE[tt]),
        Form("DCA_Z %s; DCA_{z} [cm]",  TTYPE[tt]),100,-50.0,50.0);
    hTrkDqPt[tt]  = new TH1F(Form("trkDqPt%s",  TTYPE[tt]),
        Form("q/pT residual %s; q/pT_{reco}-q/pT_{true} [GeV/c]^{-1}",TTYPE[tt]),200,-2.0,2.0);
    hTrkDpTrel[tt]= new TH1F(Form("trkDpTrel%s",TTYPE[tt]),
        Form("#Deltap_{T}/p_{T} %s; (p_{T,reco}-p_{T,true})/p_{T,true}",TTYPE[tt]),200,-2.0,2.0);
    hTrkQPtRecovsTru[tt]=new TH2F(Form("trkQPtRecovsTru%s",TTYPE[tt]),
        Form("q/pT reco vs true %s; q/pT_{true} [GeV/c]^{-1}; q/pT_{reco} [GeV/c]^{-1}",TTYPE[tt]),
        100,-5,5,100,-5,5);
  }

  // BLCVtx vertex-quality histograms
  hBLCVtxN       = new TH1F("blcVtxN",      "N BLC tracks/event; N_{BLC}",                    21,-0.5,20.5);
  hBLCVtxZ       = new TH1F("blcVtxZ",      "BLC vertex z (avg DCA-z); z_{vtx} [cm]",        100,-30.0,30.0);
  hBLCVtxZres    = new TH1F("blcVtxZres",   "BLC z_{vtx} - event z_{vtx}; #Deltaz [cm]",     100,-15.0,15.0);
  hBLCVtxZdiff   = new TH1F("blcVtxZdiff",  "BLC dcaZ_{i} - dcaZ_{j} (pairs); #Deltaz [cm]",100,-30.0,30.0);
  hBLCVtxTrkDcaZ = new TH1F("blcVtxTrkDcaZ","BLCVtx track DCA-z (after refit); DCA_{z} [cm]",100,-10.0,10.0);
}

void RunMatch(StPicoDst *dst, StFcsDb* fcsDb, int iEvent){
  StPicoEvent *event = dst->event();
  
  //Select good tracks and keep it for this and previous event
  int same=iEvent%2;
  int nTracks = dst->numberOfFwdTracks();
  ntrk[same]=0;

  // AAA diagnostic: per-event track-type count and per-track q/pT for BLC
  if(mDebug>=1 && iEvent==0){
    // print primary vertex info from first event
    printf("AAA EVENT vtxX=%.4f vtxY=%.4f vtxZ=%.4f nPriTrk=%d\n",
           event->primaryVertex().X(), event->primaryVertex().Y(), event->primaryVertex().Z(),
           (int)event->numberOfPrimaryTracks());
  }
  if(mDebug>=1){
    int nType[mNType]={0};
    int nGood[mNType]={0};
    for(int it=0; it<nTracks; it++){
      StPicoFwdTrack* t=dst->fwdTrack(it);
      int tt=t->trackType(); if(tt<0||tt>=mNType) continue;
      nType[tt]++;
      bool good = t->chi2()>0 && t->chi2()<mChi2Cut
               && t->didFitConvergeFully()
               && abs(t->numberOfFitPoints())>=(int)mNHitCut;
      if(good) nGood[tt]++;
      // print every BLC (type=1) track so we can compare real data vs MC
      if(tt==1){
        TVector3 p=t->momentum();
        TVector3 pe=t->ecalProjection();
        float qpt = (p.Perp()>0) ? t->charge()/p.Perp() : 0;
        float eta = (p.Mag()>0)  ? p.PseudoRapidity() : 0;
        printf("AAA BLC ev=%4d it=%3d Cvg=%d NHit=%3d Chi2=%8.2f"
               " pT=%6.3f pZ=%7.3f eta=%5.2f q=%+d qpT=%7.3f dcaZ=%7.2f"
               " EcalXY=%6.1f %6.1f\n",
               iEvent,it,(int)t->didFitConvergeFully(),t->numberOfFitPoints(),t->chi2(),
               p.Perp(),p.Z(),eta,t->charge(),qpt,t->dcaZ(),pe.X(),pe.Y());
      }
      // print GLOBAL (0) and PRIMARY (2) tracks for comparison with BLC
      if((tt==0||tt==2) && good){
        TVector3 p=t->momentum();
        float qpt = (p.Perp()>0) ? t->charge()/p.Perp() : 0;
        printf("AAA %s ev=%4d it=%3d NHit=%3d Chi2=%8.2f"
               " pT=%6.3f pZ=%7.3f q=%+d qpT=%7.3f dcaZ=%7.2f\n",
               (tt==0?"GLO":"PRI"),
               iEvent,it,t->numberOfFitPoints(),t->chi2(),
               p.Perp(),p.Z(),t->charge(),qpt,t->dcaZ());
      }
    }
    printf("AAA EVSUM ev=%4d nTrk=%3d  G=%3d/%3d Blc=%3d/%3d Pri=%3d/%3d Vtx=%3d/%3d\n",
           iEvent,nTracks,
           nType[0],nGood[0], nType[1],nGood[1],
           nType[2],nGood[2], nType[3],nGood[3]);
  }

  // Per-event BLC DCA-z collection for vertex quality histograms
  std::vector<double> perEvtBLCdcaZ;

  for(int it=0; it<nTracks; it++) {
    StPicoFwdTrack* t=dst->fwdTrack(it);
    TVector3 pe = t->ecalProjection();
    TVector3 ph = t->hcalProjection();
    double x=ph.X();
    double y=ph.Y();
    double r=sqrt(x*x + y*y);
    if(mDebug>1) printf("ev=%4d trk=%3d Type=%1d NHit=%3d Cvg=%1d Chi2=%8.1f pxyz=%7.2f %7.2f %7.2f E=%7.2f ProjE=%7.2f %7.2f %7.2f (%7.2f %7.2f)\n",
			iEvent,it,t->trackType(),t->numberOfFitPoints(),(int)t->didFitConvergeFully(),t->chi2(),
			t->momentum().X(),t->momentum().Y(),t->momentum().Z(),t->momentum().Mag(),
			pe.X(),pe.Y(),pe.Z(),
			t->momentum().X()/t->momentum().Z()*750.0,
			t->momentum().Y()/t->momentum().Z()*750.0
			);
    // Fill track-quality histograms for all track types with basic quality cuts
    int tt = t->trackType();
    if( tt>=0 && tt<mNType &&
        t->chi2() > 0.0 &&
        t->chi2() < mChi2Cut &&
        t->didFitConvergeFully() &&
        abs(t->numberOfFitPoints()) >= mNHitCut ){
      TVector3 p = t->momentum();
      float qpt_reco = (p.Perp()>0) ? t->charge() * 1.0f / p.Perp() : 0;
      if(p.Perp()>0) hTrkQPt[tt]->Fill(qpt_reco);
      if(p.Mag()>0)  hTrkEta[tt]->Fill(p.PseudoRapidity());
      hTrkPhi[tt] ->Fill(p.Phi());
      hTrkDcaZ[tt]->Fill(t->dcaZ());
      // MC truth matching — only when McTracks are available (MC simulation)
      // idTruth is 1-based index into McTrack array
      int truthId = t->idTruth();
      if( truthId > 0 && truthId <= (int)dst->numberOfMcTracks() ){
        StPicoMcTrack *mc = dst->mcTrack(truthId - 1);
        if( mc && mc->pt() > 0 ){
          // Use MC charge (not reco charge) for truth q/pT
          float qpt_true = mc->charge() * 1.0f / mc->pt();
          hTrkDqPt[tt]        ->Fill(qpt_reco - qpt_true);
          hTrkDpTrel[tt]      ->Fill((p.Perp() - mc->pt()) / mc->pt());
          hTrkQPtRecovsTru[tt]->Fill(qpt_true, qpt_reco);
        }
      }
      // Collect BLC and BLCVtx DCA-z for per-event vertex quality
      if(tt==1) perEvtBLCdcaZ.push_back(t->dcaZ());
      if(tt==4) hBLCVtxTrkDcaZ->Fill(t->dcaZ());
    }

    if(t->trackType()==mTrackType &&
       t->chi2() > 0.0 &&
       t->chi2() < mChi2Cut &&
       t->didFitConvergeFully() &&
       abs(t->numberOfFitPoints())>=mNHitCut &&
       fabs(x)>mHcalXminCut &&
       fabs(x)<mHcalXmaxCut &&
       fabs(y)<mYmaxCut ){
      TVector3 p=t->momentum();
      hXY[2]->Fill(pe.X(),pe.Y());
      hXY[3]->Fill(ph.X(),ph.Y());
      hE[2]->Fill(p.Mag());
      hET[2]->Fill(p.Perp());
      tx[same][0][ntrk[same]]=pe.X();
      ty[same][0][ntrk[same]]=pe.Y();
      tx[same][1][ntrk[same]]=ph.X();
      ty[same][1][ntrk[same]]=ph.Y();
      ntrk[same]++;
      if(mDebug>0) printf("ev=%6d n=%3d trk=%3d Exy=%7.2f %7.2f Hxy=%7.2f %7.2f Chi2=%6.3f Type=%1d NHit=%2d Cvg=%1d\n",
			  iEvent,nTracks,it,pe.X(),pe.Y(),ph.X(),ph.Y(),
			  t->chi2(),t->trackType(),t->numberOfFitPoints(),(int)t->didFitConvergeFully());
    }
  }
  
  // Per-event BLCVtx vertex quality
  {
    int nBLC = (int)perEvtBLCdcaZ.size();
    hBLCVtxN->Fill(nBLC);
    if(nBLC > 0){
      double zsum = 0; for(double z : perEvtBLCdcaZ) zsum += z;
      double zavg = zsum / nBLC;
      hBLCVtxZ->Fill(zavg);
      hBLCVtxZres->Fill(zavg - event->primaryVertex().Z());
    }
    for(int ii=0; ii<nBLC; ii++)
      for(int jj=ii+1; jj<nBLC; jj++)
        hBLCVtxZdiff->Fill(perEvtBLCdcaZ[ii] - perEvtBLCdcaZ[jj]);
  }

  //loop over FCS clusters
  int nClusters = dst->numberOfFcsClusters();
  for(int ic=0; ic<nClusters; ic++) {
    StPicoFcsCluster* c = dst->fcsCluster(ic);
    int det=c->detectorId();
    int ehp=fcsDb->ecalHcalPres(det);
    StThreeVectorD xyz = fcsDb->getStarXYZfromColumnRow(det,c->x(),c->y());
    double fcsx=xyz.x();
    double fcsy=xyz.y();
    double fcse=c->energy();
    double fcset=c->fourMomentum().Et();
    if(ehp==0 && (fabs(fcsx)<mEcalXminCut || fabs(fcsx)>mEcalXmaxCut)) continue;
    if(ehp==1 && (fabs(fcsx)<mHcalXminCut || fabs(fcsx)>mHcalXmaxCut)) continue;
    double fcsp=PHI(fcsx,fcsy);
    double fcsr=R(fcsx,fcsy);
    hXY[ehp]->Fill(fcsx,fcsy);
    hE[ehp]->Fill(fcse);
    hET[ehp]->Fill(fcset);
    if(mDebug>0) printf("ev=%4d clu=%3d det=%1d ehp=%1d xyz=%7.2f %7.2f E=%6.3f SigMax/Min=%6.3f %6.3f Ctg=%1d\n",
			iEvent,ic,det,ehp,xyz.x(),xyz.y(),c->energy(),
			c->sigmaMax(),c->sigmaMin(),c->category());
    for(int i=0; i<2; i++){ //same event or mixed event
      int evt=(same+i)%2;
      if(iEvent==0 && i==1) continue; //no mixed event for 1st event
      for(int it=0; it<ntrk[evt]; it++) {
	if(fcsx * tx[evt][ehp][it] > 0){ //same side
	  //if(ehp==0 && fcse>0.5) continue;  //select MIP for ecal-trk
	  //if(ehp==1 && fcse<1.0) continue; //cut off low evenry for hcal-trk
	  double trkx=tx[evt][ehp][it];
	  double trky=ty[evt][ehp][it];
	  double trkp=PHI(trkx,trky);
	  double trkr=R(trkx,trky);
	  double fcsrp=fcsr*fcsp;
	  double dx=fcsx - trkx;
	  double dy=fcsy - trky;
	  double rdp=fcsr*DPHI(fcsp,trkp);
	  double dr=fcsr - trkr;
	  int nstb=0;
	  if(fcsx>0) {nstb=0;} else {nstb=1;}
	  if(fabs(dy)<mDCut) {
	    hdx[ehp][nstb][i]->Fill(dx);  
	    hxdx[ehp][i]->Fill(fcsx,dx); 
	    hydx[ehp][i]->Fill(fcsy,dx);
	    hpdx[ehp][i]->Fill(fcsrp,dx);
	    hrdx[ehp][i]->Fill(fcsr,dx);
	  }
	  if(fcsy>0) {nstb=0;} else {nstb=1;}
	  if(fabs(dx)<mDCut){
	    hdy[ehp][nstb][i]->Fill(dy);  
	    hxdy[ehp][i]->Fill(fcsx,dy); 
	    hydy[ehp][i]->Fill(fcsy,dy);
	    hpdy[ehp][i]->Fill(fcsrp,dy); 
	    hrdy[ehp][i]->Fill(fcsr,dy);
	  }
	  if(fcsx>0) {nstb=0;} else {nstb=1;}
	  if(fabs(dr)<mDCut) {
	    hdp[ehp][nstb][i]->Fill(rdp);  
	    hxdp[ehp][i]->Fill(fcsx,rdp); 
	    hydp[ehp][i]->Fill(fcsy,rdp);
	    hpdp[ehp][i]->Fill(fcsrp,rdp); 
	    hrdp[ehp][i]->Fill(fcsr,rdp);
	  }
	  if(fcsr<70) {nstb=0;} else {nstb=1;}
	  if(fabs(rdp)<mDCut) {
	    hdr[ehp][nstb][i]->Fill(dr);  
	    hxdr[ehp][i]->Fill(fcsx,dr); 
	    hydr[ehp][i]->Fill(fcsy,dr);
	    hpdr[ehp][i]->Fill(fcsrp,dr); 
	    hrdr[ehp][i]->Fill(fcsr,dr);
	  }
	}//same side
      }//fwd track loop
    }//same event or mixed
  }//fcs cluster loop    
}

void EndMatch(){						
  printf("Writing and closing %s\n",filenameM);
  mFileM->Write();
  mFileM->Close();
}
