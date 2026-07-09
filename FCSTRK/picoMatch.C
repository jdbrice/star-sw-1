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

// ECAL position resolution vs cluster energy (electron MC, Primary track projection)
// dx = fcsX - track_ecalProjection.X(); use mTrackType=2 (Primary) so sigma_track << sigma_pos
TH2F *hFcsEdxVsE;
TH2F *hFcsEdyVsE;

// ECAL cluster vs MC truth impact point
// MC truth: propagate mc->p() through Bz=0.5T helix (z<200 cm) then straight to z=710 cm
// dx_mc = fcsX - x_mc_truth(z=710); isolates sigma_pos alone (no track reco error)
TH2F *hFcsMcDxVsE;
TH2F *hFcsMcDyVsE;

// ECAL cluster vs gamma MC truth (straight-line, exact for neutral particles)
// dx_gam = fcsX - (px/pz)*zFCS; cleanest sigma_pos measurement
TH2F *hFcsGamDxVsE;
TH2F *hFcsGamDyVsE;

static int ntrk[2]={0,0};
static double tx[2][2][mMaxTrack],ty[2][2][mMaxTrack];
// MC truth ECAL impact point per event (stored alongside track projections)
static double mcx[2][mMaxTrack],mcy[2][mMaxTrack]; // x,y at z=710 cm from MC truth helix
static double mcsx[2][mMaxTrack],mcsy[2][mMaxTrack]; // dx/dz, dy/dz slope after B-field region
static double mce[2][mMaxTrack];                   // MC truth energy (= cluster energy to match)
static int    mcOrigIdx[2][mMaxTrack];             // original 0-based MC track index
static int    nmctrk[2]={0,0};
// Gamma MC truth ECAL impact point (straight-line, no B-field)
static double gamx[2][mMaxTrack],gamy[2][mMaxTrack],game[2][mMaxTrack];
static int    ngam[2]={0,0};
// Best-match cluster to each gamma (minimum dr within mDCut)
static double gamBestDx[2][mMaxTrack],gamBestDy[2][mMaxTrack];
static double gamBestE[2][mMaxTrack],gamBestDr[2][mMaxTrack];

// Track-quality histograms by track type (no tracktype cut; basic quality only)
// Types: 0=Global 1=Beamline 2=Primary 3=FwdVtx 4=BLCVtx 5=FCSTRK
// BBB: BLCVtx = BLC tracks refitted to the beam-line-constrained forward vertex
// FCSTRK = BLCVtx tracks refitted with FCS ECAL cluster as additional measurement
static const int mNType=6;
const char *TTYPE[mNType]={"Global","Beamline","Primary","FwdVtx","BLCVtx","FCSTRK"};
TH1F *hTrkQPt[mNType];       // q/pT = signed curvature (charge sign + momentum)
TH1F *hTrkEta[mNType];      // pseudorapidity (forward acceptance)
TH1F *hTrkPhi[mNType];      // azimuthal angle (acceptance gaps)
TH1F *hTrkP[mNType];        // total momentum |p| [GeV/c]
TH1F *hTrkPt[mNType];       // transverse momentum pT [GeV/c]
TH1F *hTrkP_noNSgap[mNType];   // |p|, phi away from NS gap (|phi|<0.9 or |phi|>2.5)
TH1F *hTrkPt_noNSgap[mNType];  // pT,  phi away from NS gap
TH2F *hTrkEtaPhi[mNType];   // 2D: eta vs phi (reveals FCS acceptance gaps)
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

// FCSTRK efficiency histograms
// Denominator: MC electrons (geantId=3) with a matched ECAL cluster (dr < mEcalDrMatch)
// Numerator:   above + matched FCSTRK (type=5) track via idTruth
// Fill vs MC electron total momentum p (= energy for massless electron)
static const double mEcalDrMatch = 15.0;  // cm — MC-cluster match radius
TH1F *hFCSEffDen;   // denominator: MC e with ECAL cluster
TH1F *hFCSEffNum;   // numerator:   MC e with ECAL cluster AND FCSTRK (type=5) match
TH1F *hBLCVtxEffNum;// for comparison: MC e with ECAL cluster AND BLCVtx (type=4) match
// Per-event bookkeeping (indexed by 0-based MC track index)
static bool  mcIsEcalMatched[mMaxTrack];  // has nearby ECAL cluster this event
static bool  mcHasFCSTRK[mMaxTrack];      // has matched FCSTRK track this event
static bool  mcHasBLCVtxE[mMaxTrack];    // has matched BLCVtx track this event
static double mcTotalP[mMaxTrack];        // total MC momentum (for histogram axis)

// Per-event highest-|pT| track selection histograms (0=BLCVtx type4, 1=FCSTRK type5)
// Filling ALL tracks inflates tails with secondaries; selecting 1/event gives true resolution.
// Selection: highest |pT| converged track in the event.
TH1F *hNTrkPerEvt[4];   // N tracks per event: Global/BLC/BLCVtx/FCSTRK
TH1F *hQPtSel[4];       // q/pT of selected track (charge misID: q/pT>0 for e-)
TH1F *hDqPtSel[4];      // q/pT residual of selected track (needs MC match)

void InitMatch(int run, int set=-1, int trkType=1) {
  mTrackType = trkType;
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

  hFcsEdxVsE = new TH2F("FcsEdxVsE",
      "ECAL #sigma_{x} vs E (Primary trk); E_{cluster} [GeV]; dx = fcsX-trkX [cm]",
      50, 0, 50,  80, -20, 20);
  hFcsEdyVsE = new TH2F("FcsEdyVsE",
      "ECAL #sigma_{y} vs E (Primary trk); E_{cluster} [GeV]; dy = fcsY-trkY [cm]",
      50, 0, 50,  80, -20, 20);
  hFcsMcDxVsE = new TH2F("FcsMcDxVsE",
      "ECAL #sigma_{x} vs E (MC truth); E_{cluster} [GeV]; dx = fcsX-mcX [cm]",
      50, 0, 50,  40, -20, 20);
  hFcsMcDyVsE = new TH2F("FcsMcDyVsE",
      "ECAL #sigma_{y} vs E (MC truth); E_{cluster} [GeV]; dy = fcsY-mcY [cm]",
      50, 0, 50,  40, -20, 20);
  hFcsGamDxVsE = new TH2F("FcsGamDxVsE",
      "ECAL #sigma_{x} vs E (gamma truth); E_{cluster} [GeV]; dx = fcsX-#gamma_{X} [cm]",
      50, 0, 50,  40, -20, 20);
  hFcsGamDyVsE = new TH2F("FcsGamDyVsE",
      "ECAL #sigma_{y} vs E (gamma truth); E_{cluster} [GeV]; dy = fcsY-#gamma_{Y} [cm]",
      50, 0, 50,  40, -20, 20);

  // Track-quality histograms for all track types (basic quality cuts, no tracktype filter)
  for(int tt=0; tt<mNType; tt++){
    hTrkQPt[tt]    = new TH1F(Form("trkQPt%s",   TTYPE[tt]),
        Form("q/pT %s; q/pT [GeV/c]^{-1}",TTYPE[tt]), 50,-5.0,5.0);
    hTrkEta[tt]    = new TH1F(Form("trkEta%s",   TTYPE[tt]),
        Form("#eta %s; #eta",            TTYPE[tt]), 50, 1.5, 6.0);
    hTrkPhi[tt]    = new TH1F(Form("trkPhi%s",   TTYPE[tt]),
        Form("#phi %s; #phi [rad]",      TTYPE[tt]), 50,-3.2, 3.2);
    hTrkP[tt]      = new TH1F(Form("trkP%s",     TTYPE[tt]),
        Form("|p| %s; |p| [GeV/c]",     TTYPE[tt]), 40, 0.0,40.0);
    hTrkPt[tt]     = new TH1F(Form("trkPt%s",    TTYPE[tt]),
        Form("p_{T} %s; p_{T} [GeV/c]", TTYPE[tt]), 40, 0.0, 5.0);
    hTrkP_noNSgap[tt]  = new TH1F(Form("trkP_noNSgap%s",  TTYPE[tt]),
        Form("|p| (away NS gap) %s; |p| [GeV/c]",  TTYPE[tt]), 40, 0.0,40.0);
    hTrkPt_noNSgap[tt] = new TH1F(Form("trkPt_noNSgap%s", TTYPE[tt]),
        Form("p_{T} (away NS gap) %s; p_{T} [GeV/c]",TTYPE[tt]), 40, 0.0, 5.0);
    hTrkEtaPhi[tt] = new TH2F(Form("trkEtaPhi%s",TTYPE[tt]),
        Form("#eta vs #phi %s; #eta; #phi [rad]",TTYPE[tt]),
        30,1.5,6.0, 32,-3.2,3.2);
    hTrkDcaZ[tt]   = new TH1F(Form("trkDcaZ%s",  TTYPE[tt]),
        Form("DCA_Z %s; DCA_{z} [cm]",  TTYPE[tt]), 50,-50.0,50.0);
    hTrkDqPt[tt]  = new TH1F(Form("trkDqPt%s",  TTYPE[tt]),
        Form("q/pT residual %s; q/pT_{reco}-q/pT_{true} [GeV/c]^{-1}",TTYPE[tt]),100,-2.0,2.0);
    hTrkDpTrel[tt]= new TH1F(Form("trkDpTrel%s",TTYPE[tt]),
        Form("#Deltap_{T}/p_{T} %s; (p_{T,reco}-p_{T,true})/p_{T,true}",TTYPE[tt]),100,-2.0,2.0);
    hTrkQPtRecovsTru[tt]=new TH2F(Form("trkQPtRecovsTru%s",TTYPE[tt]),
        Form("q/pT reco vs true %s; q/pT_{true} [GeV/c]^{-1}; q/pT_{reco} [GeV/c]^{-1}",TTYPE[tt]),
        50,-5,5,50,-5,5);
  }

  // BLCVtx vertex-quality histograms
  hBLCVtxN       = new TH1F("blcVtxN",      "N BLC tracks/event; N_{BLC}",                    21,-0.5,20.5);
  hBLCVtxZ       = new TH1F("blcVtxZ",      "BLC vertex z (avg DCA-z); z_{vtx} [cm]",        100,-50.0,50.0);
  hBLCVtxZres    = new TH1F("blcVtxZres",   "BLC z_{vtx} - event z_{vtx}; #Deltaz [cm]",     100,-50.0,50.0);
  hBLCVtxZdiff   = new TH1F("blcVtxZdiff",  "BLC dcaZ_{i} - dcaZ_{j} (pairs); #Deltaz [cm]",100,-50.0,50.0);
  hBLCVtxTrkDcaZ = new TH1F("blcVtxTrkDcaZ","BLCVtx track DCA-z (after refit); DCA_{z} [cm]",100,-50.0,50.0);

  // FCSTRK efficiency histograms (MC electrons with ECAL cluster as denominator)
  hFCSEffDen    = new TH1F("fcsEffDen",    "MC e^{-} with ECAL cluster (denominator); p_{MC} [GeV/c]", 50,0,50);
  hFCSEffNum    = new TH1F("fcsEffNum",    "MC e^{-} with ECAL cluster + FCSTRK match; p_{MC} [GeV/c]",50,0,50);
  hBLCVtxEffNum = new TH1F("blcVtxEffNum","MC e^{-} with ECAL cluster + BLCVtx match; p_{MC} [GeV/c]",50,0,50);

  // Per-event selected-track histograms (1 per event, highest |pT|)
  const char *SEL[4]={"Global","BLC","BLCVtx","FCSTRK"};
  for(int is=0; is<4; is++){
    hNTrkPerEvt[is] = new TH1F(Form("nTrkPerEvt%s",SEL[is]),
                                Form("N %s tracks/event; N tracks; Events",SEL[is]), 11,-0.5,10.5);
    hQPtSel[is]     = new TH1F(Form("qPtSel%s",SEL[is]),
                                Form("%s q/p_{T} (highest |p_{T}| per event); q/p_{T} [(GeV/c)^{-1}]",SEL[is]),
                                400,-10,10);
    hDqPtSel[is]    = new TH1F(Form("dqPtSel%s",SEL[is]),
                                Form("%s #Deltaq/p_{T} selected (highest |p_{T}|); #Delta(q/p_{T}) [(GeV/c)^{-1}]",SEL[is]),
                                400,-10,10);
  }
}

void RunMatch(StPicoDst *dst, StFcsDb* fcsDb, int iEvent){
  StPicoEvent *event = dst->event();

  //Select good tracks and keep it for this and previous event
  int same=iEvent%2;
  int nTracks = dst->numberOfFwdTracks();
  ntrk[same]=0;

  // Reset per-event FCSTRK efficiency bookkeeping
  int nMCTotal = dst->numberOfMcTracks();
  for(int im=0; im<nMCTotal && im<mMaxTrack; im++){
    mcIsEcalMatched[im] = false;
    mcHasFCSTRK[im]     = false;
    mcHasBLCVtxE[im]    = false;
    StPicoMcTrack *mc = dst->mcTrack(im);
    mcTotalP[im] = (mc) ? mc->p().Mag() : 0.0;
  }

  // Per-event track selection: collect candidates for highest-|pT| selection
  // sel[0]=BLCVtx(type4), sel[1]=FCSTRK(type5)
  struct SelTrack { float pt, qpt_reco, qpt_true; bool hasTrue; };
  std::vector<SelTrack> selCand[4];

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
      // print FCSTRK (type=5) tracks for FCSTRK momentum comparison
      if(tt==5){
        TVector3 p=t->momentum();
        float qpt = (p.Perp()>0) ? t->charge()/p.Perp() : 0;
        float eta = (p.Mag()>0)  ? p.PseudoRapidity() : 0;
        printf("AAA FCS ev=%4d it=%3d Cvg=%d NHit=%3d Chi2=%8.2f"
               " pT=%6.3f pZ=%7.3f eta=%5.2f q=%+d qpT=%7.3f\n",
               iEvent,it,(int)t->didFitConvergeFully(),t->numberOfFitPoints(),t->chi2(),
               p.Perp(),p.Z(),eta,t->charge(),qpt);
      }
    }
    printf("AAA EVSUM ev=%4d nTrk=%3d  G=%3d/%3d Blc=%3d/%3d Pri=%3d/%3d Vtx=%3d/%3d BLCVtx=%3d/%3d FCS=%3d/%3d\n",
           iEvent,nTracks,
           nType[0],nGood[0], nType[1],nGood[1],
           nType[2],nGood[2], nType[3],nGood[3],
           nType[4],nGood[4], nType[5],nGood[5]);
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
      hTrkPhi[tt]   ->Fill(p.Phi());
      hTrkP[tt]     ->Fill(p.Mag());
      hTrkPt[tt]    ->Fill(p.Perp());
      double phi = p.Phi();
      if(fabs(phi) < 0.9 || fabs(phi) > 2.5){
        hTrkP_noNSgap[tt] ->Fill(p.Mag());
        hTrkPt_noNSgap[tt]->Fill(p.Perp());
      }
      if(p.Mag()>0)  hTrkEtaPhi[tt]->Fill(p.PseudoRapidity(), p.Phi());
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
      if(tt==1 && t->dcaXY() < 10.0) perEvtBLCdcaZ.push_back(t->dcaZ());
      if(tt==4) hBLCVtxTrkDcaZ->Fill(t->dcaZ());

      // FCSTRK efficiency bookkeeping: record which MC electrons have BLCVtx or FCSTRK tracks
      int truthId2 = t->idTruth();
      if(truthId2 > 0 && truthId2 <= (int)dst->numberOfMcTracks() && truthId2 <= mMaxTrack){
        StPicoMcTrack *mc2 = dst->mcTrack(truthId2 - 1);
        if(mc2 && mc2->geantId() == 3){  // geantId=3 is electron
          if(tt == 4) mcHasBLCVtxE[truthId2-1] = true;
          if(tt == 5) mcHasFCSTRK[truthId2-1]  = true;
        }
      }

      // Collect per-event candidates for highest-|pT| selection (Global/BLC/BLCVtx/FCSTRK)
      if(tt==0 || tt==1 || tt==4 || tt==5){
        int is = (tt==0) ? 0 : (tt==1) ? 1 : (tt==4) ? 2 : 3;
        SelTrack st;
        st.pt       = p.Perp();
        st.qpt_reco = qpt_reco;
        st.hasTrue  = false;
        st.qpt_true = 0;
        int tid = t->idTruth();
        if(tid > 0 && tid <= (int)dst->numberOfMcTracks()){
          StPicoMcTrack *mcS = dst->mcTrack(tid - 1);
          if(mcS && mcS->pt() > 0){
            st.qpt_true = mcS->charge() * 1.0f / mcS->pt();
            st.hasTrue  = true;
          }
        }
        selCand[is].push_back(st);
      }
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

  // Per-event highest-|pT| selection: fill N-per-event and selected-track histograms
  for(int is=0; is<4; is++){
    hNTrkPerEvt[is]->Fill((int)selCand[is].size());
    if(selCand[is].empty()) continue;
    // find highest |pT|
    int bestIdx = 0;
    for(int ic=1; ic<(int)selCand[is].size(); ic++)
      if(selCand[is][ic].pt > selCand[is][bestIdx].pt) bestIdx = ic;
    const SelTrack &best = selCand[is][bestIdx];
    hQPtSel[is]->Fill(best.qpt_reco);
    if(best.hasTrue) hDqPtSel[is]->Fill(best.qpt_reco - best.qpt_true);
  }

  // Compute MC truth ECAL impact point for all MC tracks in this event.
  // Method: helix through constant Bz=0.5 T for z in [0, 200 cm], then straight to z=710 cm.
  // Conversion: R[cm] = pT[GeV/c] / (0.3 * B[T]) * 100
  // Only store if x at z=710 is in ECAL acceptance (|x|>15, |x|<150, |y|<100).
  nmctrk[same] = 0;
  {
    const double BzSigned = -0.5; // Tesla: simulation uses field=-5kGauss=-0.5T
    const double zB   = 200.0;   // cm — B-field exit point (approximate)
    const double zFCS = 710.0;   // cm — ECAL z
    int nMC = dst->numberOfMcTracks();
    for(int im=0; im<nMC && nmctrk[same]<mMaxTrack; im++){
      StPicoMcTrack *mc = dst->mcTrack(im);
      if(!mc) continue;
      if(mc->p().Mag() < 1.0) continue;  // skip shower secondaries: gun electrons have |p|≈E≥1 GeV/c; shower secondaries have |p|≪1 GeV/c
      TVector3 p = mc->p();
      double pT = p.Perp();
      double pz = p.Z();
      if(pT<=0 || pz<=0) continue;  // require forward, non-zero pT
      double q  = mc->charge();      // +1 or -1
      // radius of curvature in x-y plane [cm] — always positive
      double R  = pT / (0.3 * TMath::Abs(BzSigned)) * 100.0;
      // azimuthal angle of initial transverse momentum
      double phi0 = TMath::ATan2(p.Y(), p.X());
      // total turning angle in x-y by z=zB:
      // transverse arc length = z * (pT/pz), radius R_xy = pT/(0.3*B) [m]
      // phi_turn = arc/R = 0.3*B*z[m]/pz  (independent of pT!)
      double phiTurn = 0.3 * BzSigned * (zB/100.0) / pz * (double)q;  // sign: q*BzSigned
      // position at z=zB after helix
      // x(zB) = x0 + R*(sin(phi0+phiTurn) - sin(phi0))
      // y(zB) = y0 + R*(-cos(phi0+phiTurn) + cos(phi0))
      double xB = R * (TMath::Sin(phi0 + phiTurn) - TMath::Sin(phi0));
      double yB = R * (-TMath::Cos(phi0 + phiTurn) + TMath::Cos(phi0));
      // slope at z=zB: (dx/dz, dy/dz) from rotated momentum direction
      double sxB = TMath::Sin(phi0 + phiTurn) * (pT/pz);
      double syB = TMath::Cos(phi0 + phiTurn) * (pT/pz); // note: rotated direction
      // Actually: slope = (px_exit/pz); px_exit = pT*cos(phi0+phiTurn), py_exit = pT*sin(phi0+phiTurn)
      double sxB2 = TMath::Cos(phi0 + phiTurn) * (pT/pz);
      double syB2 = TMath::Sin(phi0 + phiTurn) * (pT/pz);
      // propagate straight from zB to zFCS
      double xFCS = xB + sxB2 * (zFCS - zB);
      double yFCS = yB + syB2 * (zFCS - zB);
      // acceptance cut
      if(fabs(xFCS)<mEcalXminCut || fabs(xFCS)>mEcalXmaxCut || fabs(yFCS)>mYmaxCut) continue;
      // store total kinetic energy as a proxy for cluster energy (massless electron approx)
      double eTotal = p.Mag();
      mcx[same][nmctrk[same]] = xFCS;
      mcy[same][nmctrk[same]] = yFCS;
      mcsx[same][nmctrk[same]] = sxB2;
      mcsy[same][nmctrk[same]] = syB2;
      mce[same][nmctrk[same]] = eTotal;
      mcOrigIdx[same][nmctrk[same]] = im;  // preserve original MC track index
      nmctrk[same]++;
    }
  }

  // Gamma MC truth: straight-line from vertex to z=710 cm.
  // Exact for neutral particles (no B-field bending). pdgId==22.
  ngam[same] = 0;
  {
    const double zFCS = 710.0;
    int nMC = dst->numberOfMcTracks();
    for(int im=0; im<nMC && ngam[same]<mMaxTrack; im++){
      StPicoMcTrack *mc = dst->mcTrack(im);
      if(!mc) continue;
      if(mc->geantId() != 1) continue;   // geantId=1 is photon (pdgId() not implemented)
      if(mc->idVtxStart() != 1) continue; // primary particles only (skip shower secondaries)
      TVector3 p = mc->p();
      if(p.Z() <= 0) continue;
      double xFCS = p.X()/p.Z() * zFCS;
      double yFCS = p.Y()/p.Z() * zFCS;
      if(fabs(xFCS)<mEcalXminCut || fabs(xFCS)>mEcalXmaxCut || fabs(yFCS)>mYmaxCut) continue;
      gamx[same][ngam[same]] = xFCS;
      gamy[same][ngam[same]] = yFCS;
      game[same][ngam[same]] = p.Mag();
      ngam[same]++;
    }
  }

  // Reset gamma best-match before cluster loop
  for(int ig=0; ig<ngam[same]; ig++) gamBestDr[same][ig] = 999.0;

  //loop over FCS clusters
  int nClusters = dst->numberOfFcsClusters();
  for(int ic=0; ic<nClusters; ic++) {
    StPicoFcsCluster* c = dst->fcsCluster(ic);
    int det=c->detectorId();
    int ehp=fcsDb->ecalHcalPres(det);
    StThreeVectorD xyz = fcsDb->getStarXYZfromColumnRow(det,c->x(),c->y());
    double fcsx=xyz.x();
    double fcsy=xyz.y();
    double fcsz=xyz.z(); // STAR z of cluster (shower-max depth, ~725 cm for ECAL)
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
	    if(ehp==0 && i==0) hFcsEdxVsE->Fill(fcse, dx);
	  }
	  if(fcsy>0) {nstb=0;} else {nstb=1;}
	  if(fabs(dx)<mDCut){
	    hdy[ehp][nstb][i]->Fill(dy);
	    hxdy[ehp][i]->Fill(fcsx,dy);
	    hydy[ehp][i]->Fill(fcsy,dy);
	    hpdy[ehp][i]->Fill(fcsrp,dy);
	    hrdy[ehp][i]->Fill(fcsr,dy);
	    if(ehp==0 && i==0) hFcsEdyVsE->Fill(fcse, dy);
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

    // MC truth (helix): match ECAL cluster to MC track impact in same event only
    // Correct truth from z=710 to actual shower-max z (fcsz) using stored slope
    if(ehp==0){
      for(int im=0; im<nmctrk[same]; im++){
        if(fcsx * mcx[same][im] <= 0) continue;
        double mcx_cor = mcx[same][im] + mcsx[same][im] * (fcsz - 710.0);
        double mcy_cor = mcy[same][im] + mcsy[same][im] * (fcsz - 710.0);
        double dx_mc = fcsx - mcx_cor;
        double dy_mc = fcsy - mcy_cor;
        if(fabs(dy_mc)<mDCut) hFcsMcDxVsE->Fill(fcse, dx_mc);
        if(fabs(dx_mc)<mDCut) hFcsMcDyVsE->Fill(fcse, dy_mc);
        // FCSTRK efficiency: mark this MC electron as having a matching ECAL cluster
        double dr_mc = sqrt(dx_mc*dx_mc + dy_mc*dy_mc);
        if(dr_mc < mEcalDrMatch){
          int origIdx = mcOrigIdx[same][im];
          if(origIdx >= 0 && origIdx < mMaxTrack) mcIsEcalMatched[origIdx] = true;
        }
      }
    }
    // Gamma truth: update best-match (closest cluster) per gamma
    // Scale truth from z=710 to actual shower-max z (fcsz): gamx was computed at z=710
    if(ehp==0){
      for(int ig=0; ig<ngam[same]; ig++){
        if(fcsx * gamx[same][ig] <= 0) continue;
        double gamx_cor = gamx[same][ig] * fcsz / 710.0;
        double gamy_cor = gamy[same][ig] * fcsz / 710.0;
        double dx_gam = fcsx - gamx_cor;
        double dy_gam = fcsy - gamy_cor;
        double dr = sqrt(dx_gam*dx_gam + dy_gam*dy_gam);
        if(dr < gamBestDr[same][ig]){
          gamBestDr[same][ig] = dr;
          gamBestDx[same][ig] = dx_gam;
          gamBestDy[same][ig] = dy_gam;
          gamBestE[same][ig]  = fcse;
        }
      }
    }
  }//fcs cluster loop

  // Fill gamma histograms with best-match cluster (one fill per gamma per event)
  for(int ig=0; ig<ngam[same]; ig++){
    if(gamBestDr[same][ig] >= mDCut) continue;
    hFcsGamDxVsE->Fill(gamBestE[same][ig], gamBestDx[same][ig]);
    hFcsGamDyVsE->Fill(gamBestE[same][ig], gamBestDy[same][ig]);
  }

  // FCSTRK efficiency fill: for each MC electron with ECAL cluster, count BLCVtx and FCSTRK
  for(int im=0; im<nmctrk[same]; im++){
    int origIdx = mcOrigIdx[same][im];
    if(origIdx < 0 || origIdx >= mMaxTrack) continue;
    if(!mcIsEcalMatched[origIdx]) continue;  // require ECAL cluster match
    // Only electrons (geantId=3) — nmctrk already requires forward pT>0 but any species
    StPicoMcTrack *mc = dst->mcTrack(origIdx);
    if(!mc || mc->geantId() != 3) continue;
    double p_mc = mcTotalP[origIdx];
    hFCSEffDen->Fill(p_mc);                                    // denominator
    if(mcHasBLCVtxE[origIdx]) hBLCVtxEffNum->Fill(p_mc);     // BLCVtx match
    if(mcHasFCSTRK[origIdx])  hFCSEffNum->Fill(p_mc);         // FCSTRK match
  }
}

void EndMatch(){						
  printf("Writing and closing %s\n",filenameM);
  mFileM->Write();
  mFileM->Close();
}
