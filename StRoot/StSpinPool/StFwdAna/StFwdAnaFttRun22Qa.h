/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  To generate a ROOT file of histograms related to the FTT(sTGC) in the produced MuDsts from RHIC Run 22 that can be used to do quality assurance (QA) on the data.

  DESCRIPTION
  This class inherits from StFwdAnaVirtual and contains many histograms to be used for Quality Assurance (QA) of MuDst files that contain the small-strip Thin Gap Chamber (sTGC/FTT) data. It uses #LoadHists() to ease histogram creation and management. #DoMake() is used to fill QA histograms related to event information.

  LOG
  @[July 27, 2026] > First instance
  @[August 17, 2026] > Added various paint functions for the QA histograms
*/

#ifndef STFWDANA_STFWDANAFTTRUN22QA_HH
#define STFWDANA_STFWDANAFTTRUN22QA_HH

//C/C++ Headers
#include <iostream>
#include <sstream>

//STAR Headers
#include "StFttDbMaker/StFttDb.h"

//Custom headers in this folder
#include "StFwdAnaVirtual.h"

class StFwdAnaFttRun22Qa : public StFwdAnaVirtual
{
 public:
  StFwdAnaFttRun22Qa();
  ~StFwdAnaFttRun22Qa();
  
  virtual UInt_t LoadHists(TFile* file, HistManager* histman, StFwdAnaData* anadata);
  virtual Int_t DoMake(StFwdAnaData* anadata);

  void Paint_RawHitAdc(TCanvas* canvas, int quad, int plane, const char* savename) const;
  void Paint_RawHitBcid(TCanvas* canvas, int quad, int plane, const char* savename) const;
  void Paint_RawHitDBcid(TCanvas* canvas, int quad, int plane, const char* savename) const;
  void Paint_RawHitTb(TCanvas* canvas, int quad, int plane, const char* savename) const;
  void Paint_RawHitTime(TCanvas* canvas, int quad, int plane, const char* savename) const;
  void Paint_RawHitNchs(TCanvas* canvas, const char* savename) const;

protected:
  TH1* mH2S_FttRawHit_adcVch[StFttDb::nVMM];        ///< raw hit adc vs. channel for every hardware "vmm_id" from StFttDb
  TH1* mH2S_FttRawHit_bcidVch[StFttDb::nVMM];       ///< raw hit BCID vs channel for every hardware "vmm_id" from StFttDb
  TH1* mH2S_FttRawHit_dbcidVch[StFttDb::nVMM];      ///< raw hit Delta BCID vs channel for every hardware "vmm_id" from StFttDb
  TH1* mH2S_FttRawHit_tbVch[StFttDb::nVMM];         ///< raw hit timebin (tb) vs channel for every hardware "vmm_id" from StFttDb
  TH1* mH2S_FttRawHit_timeVch[StFttDb::nVMM];       ///< raw hit time vs channel for every hardware "vmm_id" from StFttDb

  TH1* mH2S_FttRawHit_nchsVvmm = 0;                     ///< raw hit number of channels hit vs every hardware "vmm_id" from StFttDb

  ClassDef(StFwdAnaFttRun22Qa,1)
};

#endif
  
