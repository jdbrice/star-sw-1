/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  To generate a ROOT file of histograms related to the FST in the produced MuDsts from RHIC Run 22 that can be used to do quality assurance (QA) on the data.

  DESCRIPTION
  This class inherits from StFwdAnaVirtual and contains many histograms to be used for Quality Assurance (QA) of MuDst files that contain the Forward Silicon Tracker (FST) data. It uses #LoadHists() to ease histogram creation and management. #DoMake() is used to fill QA histograms related to event information.

  LOG
  @[July 27, 2026] > First instance that mostly copied from #StFstQAMaker
  @[August 14, 2026] > Added two boolean flags #mRawHitOn, and #mHitOn which can be used to toggle QA for FST Raw Hits and/or FST Hits. The main reason is because MuDst trees don't contain FST Raw Hit information by default so it helps to save on number of histograms generated. Also changed name of histograms to include "Fst" since looking at the ROOT file of histograms, both FST and FTT use the same "RawHit" and "Hit" style of classes so it was hard to distinguish at a first glance which histograms are coming from FST QA and which from FTT QA.
  @[August 15, 2026] > Added various paint functions for the QA histograms
*/

#ifndef STFWDANA_STFWDANAFSTRUN22QA_HH
#define STFWDANA_STFWDANAFSTRUN22QA_HH

//C/C++ Headers
#include <iostream>
#include <sstream>

//STAR Headers
#include "StEvent/StFstConsts.h"

//Custom headers in this folder
#include "StFwdAnaVirtual.h"

class StFwdAnaFstRun22Qa : public StFwdAnaVirtual
{
 public:
  StFwdAnaFstRun22Qa();
  ~StFwdAnaFstRun22Qa();
  
  virtual UInt_t LoadHists(TFile* file, HistManager* histman, StFwdAnaData* anadata);
  virtual Int_t DoMake(StFwdAnaData* anadata);
  
  void setRawHitQa(bool val){ mRawHitOn = val; }
  void setHitQa(bool val){ mHitOn = val; }

  void Paint_FstRawHitStrip_rVphi  (TCanvas* canvas, int disk, const char* savename) const;
  void Paint_FstHitStripMean_rVphi (TCanvas* canvas, int disk, const char* savename) const;
  void Paint_FstHit_rVphi          (TCanvas* canvas, int disk, const char* savename) const;
  void Paint_FstHit_apvVgeoid      (TCanvas* canvas, int disk, const char* savename) const;
  void Paint_FstHitGlobal_yVx      (TCanvas* canvas, int disk, const char* savename) const;
  void Paint_FstHitGlobal_rVphi    (TCanvas* canvas, int disk, const char* savename) const;
  void Paint_FstRawHit_adcVgeoid   (TCanvas* canvas, const char* savename) const;
  void Paint_FstRawHit_adcerrVgeoid(TCanvas* canvas, const char* savename) const;
  void Paint_FstRawHit_maxtbVapv   (TCanvas* canvas, const char* savename) const;
  void Paint_FstHit_adcVid         (TCanvas* canvas, const char* savename) const;
  void Paint_FstHit_adcerrVid      (TCanvas* canvas, const char* savename) const;
  void Paint_FstHit_maxtbVid       (TCanvas* canvas, const char* savename) const;
  void Paint_Fst_nrawhitsVid       (TCanvas* canvas, const char* savename) const;
  void Paint_Fst_nhitsVid          (TCanvas* canvas, const char* savename) const;
  void Paint_FstHit_nrawhitsVid    (TCanvas* canvas, const char* savename) const;
  void Paint_FstHit_nrawhitsrVid   (TCanvas* canvas, const char* savename) const;
  void Paint_FstHit_nrawhitsphiVid (TCanvas* canvas, const char* savename) const;

  void PaintFstRawHitQa(TCanvas* canvas, const char* savename="FstRawHitQa.png") const;
  void PaintFstHitQa(TCanvas* canvas, const char* savename="FstHitQa.png") const;
  
protected:
  bool mRawHitOn = false;                        ///< Flag to turn on/off QA histograms for FST raw hits. Since these aren't normally stored in MuDsts it is off by default. Idea is to keep number of histograms smaller if no data present. Since there is no way to check before histogram creation if data contains FST Raw Hits; using a boolean flag
  bool mHitOn = true;                            ///< Flag to turn on/off QA histograms for FST hits. Since these are "highest" level of FST reconstruction. Keep on by default
  
  //position
  TH1* mH2S_FstRawHitStrip_rVphi[kFstNumSensors];   ///< raw hit phistrip vs. rstrip per sensor
  TH1* mH2S_FstHitStripMean_rVphi[kFstNumSensors];  ///< hit mean phistrip vs. mean rstrip per sensor
  TH1* mH2S_FstHit_rVphi[kFstNumDisk];              ///< hit map in r vs. phi per disk
  TH1* mH2S_FstHit_apvVgeoid[kFstNumDisk];          ///< hit map in APV geometry Id vs. module geometry Id per disk
  TH1* mH2S_FstHitGlobal_yVx[kFstNumDisk];          ///< hit global x vs. y per disk
  TH1* mH2S_FstHitGlobal_rVphi[kFstNumDisk];        ///< hit global r vs. phi per disk
  //Charge
  TH1* mH2S_FstRawHit_adcVgeoid[kFstNumTimeBins];   ///< Charge (ADC) vs channel ID over all time bins
  TH1* mH2S_FstRawHit_adcerrVgeoid = 0;             ///< RMS noise vs channel ID
  TH1* mH2S_FstRawHit_maxtbVapv = 0;                ///< Raw hit max ADC time bin vs APV electronics ID [48*(ARC-1)+16*ARM+APV]
  TH1* mH2S_FstHit_adcVid = 0;                      ///< Charge vs sensorID
  TH1* mH2S_FstHit_adcerrVid = 0;                   ///< Charge uncertainty vs sensorID
  TH1* mH2S_FstHit_maxtbVid = 0;	                 ///< hit max ADC time bin vs sensorID
  //hit or raw hit number
  TH1* mH2S_Fst_nrawhitsVid = 0;                     ///< number of raw hits vs sensor Id
  TH1* mH2S_Fst_nhitsVid = 0;	                 ///< number of hits vs sensor Id
  //TProfile* numOfRawHits_EventId[kFstNumSensors];
  //cluster size
  TH1* mH2S_FstHit_nrawhitsVid = 0;                 ///< hit cluster size
  TH1* mH2S_FstHit_nrawhitsrVid = 0;                ///< hit cluster size in R direction
  TH1* mH2S_FstHit_nrawhitsphiVid = 0;              ///< hit cluster size in Phi direction  

  ClassDef(StFwdAnaFstRun22Qa,1);
};

#endif
  
