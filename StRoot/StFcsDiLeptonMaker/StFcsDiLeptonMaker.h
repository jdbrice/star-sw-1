#ifndef STAR_StFcsDiLeptonMaker
#define STAR_StFcsDiLeptonMaker                              

#include "StMaker.h"

class StFcsDb;
class StFcsCollection;
class TH1F;
class TH2F;

class StFcsDiLeptonMaker : public StMaker {
public: 
  StFcsDiLeptonMaker(const char *name="FcsDiLepton");
  virtual ~StFcsDiLeptonMaker();
  virtual Int_t Init();
  virtual Int_t  Make();
  virtual Int_t Finish();

  void setFileName(char* file){mFilename=file;}

protected:
  
private:
  StFcsDb *mFcsDb=0; 
  StFcsCollection *mFcsCollection=0; 

  TFile* mFile=0;
  char* mFilename=0;
  
  enum {mNCut=7};
  float mETCut=0.8;          //GeV for single electron
  float mETotCut=0.3;        //E_lepton/ETOT ratio cut
  float mHTotCut=0.5;        //E_lepton/HTOT ratio cut
  float mConeR=0.5;          //Isolation Cone Radius
  float mConeCut=0.0;//0.6;TEMP!        //E_lepton/Cone ratio cut
  float mSigmaMaxCut=0.7;    //Cluster Sigma Max cut
  float mETPTCutLow=0.3;     //FcsET/TrackPT > ThresholdLow
  float mETPTCutHigh=1.7;    //FcsET/TrackPT < ThresholdHigh

  TH1F *mE_allclusters;//ISAAC
  
  TH1F *mETot[mNCut];
  TH1F *mHTot[mNCut];
  TH1F *mCone[mNCut];
  TH1F *mSigmax[mNCut];
  TH1F *mPToverET[mNCut];
  TH1F *mChargeSum[mNCut];

  TH1F *mET[mNCut]; 
  TH1F *mEZ[mNCut]; 
  TH1F *mM[mNCut]; 
  TH1F *mZ[mNCut]; 
  TH1F *mCosT[mNCut]; 
  TH1F *mPhi[mNCut]; 

  TH2F *mET12[mNCut]; 
  TH2F *mXFPT[mNCut]; 
  TH2F *mXY[mNCut]; 
  TH2F *mPTET[mNCut]; 

  ClassDef(StFcsDiLeptonMaker,1)
};

#endif

 
