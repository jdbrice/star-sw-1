/*
  AUTHOR
  David Kapukchyan

  PURPOSE
  The purpose of this class is to act as a bridge between the STAR makers and the "FcsAna" strategy framework. It will also fill the basic data structure of #StFwdAnaData and handle the calls of different sub-analysis modules (#StMuFcsVirtualAna). It will also manage the histogram objects created by the various sub-analysis modules by utilizing #HistManager.
  
  DESCRIPTION
  In the "FcsAna" strategy framework this class contains a #HistManager, and a #StFwdAnaData that sub-analysis modules will utilize. Sub-analysis modules will inherit from #StMuFcsVirtualAna and should be added to this class. Contains a #StFwdAnaData, a #HistManager and a vector of #StMuFcsVirtualAna objects. It will fill #StFwdAnaData with the relevant information by utilizing the filenames #mPolDataFilename (For polarization information taken from https://wiki.bnl.gov/rhicspin/Results_(Polarimetry)), #mFcsTrigFilename (Generated from GetTrigList.csh which is a text file of all FCS triggers names and their offline id as well as the starting and ending run number)  , #mFilename (name of file to save histograms to, passed to #mHistManager). Inside of #Init(), it will call all the #StMuFcsVirtualAna::LoadHists(). Inside of #Make() it will call all the #StMuFcsVirtualAna::DoMake() methods. The rest of the code handles everything else so all you need to do is insert your algorithm for what you want "#Make()" to do using the #addAna() call. #InitRun() will grab the FCS database object.

    Creates a histogram to keep track of make calls. Grabs the MuDstMaker pointer and the related pointers and stores them in #StFwdAnaData. Fills the #FcsEventInfo object in #StFwdAnaData

  LOG
  @[January 9, 2026] > Copied from *StMuFcsTreeMaker* and modified to do the relevant things. Since I am keeping *StMuFcsTreeMaker* please see that log for all relevant developments of the analysis code.
  @[June 8, 2026] > Implemented #StFwdAnaData::mEvent
  @[June 30, 2026] > Changed name from StMuFcsAnaDataMaker to StFwdAnaDataMaker to be consistent with new naming that is more general for STAR forward analysis. Changed StMuFcsAnaData::mEvtInfo to StFwdAnaData::mEvtData. Changed StMuFcsVirtualAna to StFwdAnaVirtual
  @[July 24, 2026] > Implemented #LoadGraphs() and #FillGraphs() which is for doing run by run qa on the vertex histograms
  @[August 3, 2026] > Implemented loading of MuFst collection and the raw hits and hits in that collection, MuFtt collection, and MuFwdTrack collection in #Make()
  @[August 13, 2026] > Got rid of extraneous headers. #Make() now fills #StFwdAnaData with pointers related to FTT. Also, contains commented out code for getting FST collection from StEvent
*/


#ifndef STFWDANA_STFWDANADATAMAKER_HH
#define STFWDANA_STFWDANADATAMAKER_HH

//C/C++ Headers
#include <iostream>

//ROOT Headers
#include "TString.h"
#include "TPolyLine.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGeoPolygon.h"

//STAR Headers
#include "St_base/StMessMgr.h"       //STAR Logger
#include "St_base/Stypes.h"          //StMaker return values
#include "StChain/StMaker.h"         //STAR base class for analysis
#include "StEvent/StTriggerData.h"   //STAR trigger information
#include "StEvent/StTriggerId.h"

//Custom Headers
#include "StFwdAnaVirtual.h"
#include "StFwdAnaData.h"

class StFwdAnaDataMaker : public StMaker {
public:
  
  StFwdAnaDataMaker(const Char_t* name = "MuFcsAnaData");
  ~StFwdAnaDataMaker();
  
  UInt_t addAna(StFwdAnaVirtual* ana);
  
  StFwdAnaData* anaData()const{ return mAnaData; }
  HistManager* getHists()const{ return mHists; }
  void setPolDataFilename(const char* name) { mPolDataFilename = name; }
  void setFcsTrigFilename(const char* name) { mFcsTrigFilename = name; }
  void setOutFilename(const char* name) { mFilename = name; }
  void setAnaData( StFwdAnaData* anadata);
  void setHistManager( HistManager* hm );
  
  UInt_t LoadDataFromFile(TFile* file);       ///< Load tree and histograms from file, calls #StFwdAnaVirtual::LoadHists()
  UInt_t LoadGraphsFromFile(TFile* file);     ///< Load graphs from file, calls #StFwdAnaVirtual::LoadGraphs() and creates graph array
  void FillGraphs(int irun);                  ///< Call #StFwdAnaVirtualAna::FillGraphs() to fill graphs
  void WriteGraphs(Int_t option=0, Int_t bufsize=0);              ///< Write TGraphs in #mGraphs to current gDirectory

  virtual Int_t Init();
  virtual Int_t InitRun(int runnumber);
  virtual void Clear(Option_t* option="");    ///< Gets called before #Make() in #StChain::EventLoop()
  virtual Int_t Make();                       ///< Call #StFwdAnaVirtual::DoMake() on internal analysis modules
  virtual Int_t Finish();
  
protected:
  TString mPolDataFilename = "";                ///< File to read Polarization information
  TString mFcsTrigFilename = "";                ///< File to read Fcs trigger information
  TString mFilename = "";                       ///< File name to save to. Passed to #mHists
  StFwdAnaData* mAnaData = 0;                   ///< Data structure that contains the pointers to the needed makers
  HistManager* mHists = 0;                      ///< Manage loading and saving histograms
  TObjArray* mGraphs = 0;                       ///< Manage any graphs created

  TH1* mH1D_Entries = 0;                        ///< Number of events processed no cuts (i.e. "Make" calls)
  TH1* mH1F_RndmSpin = 0;                       ///< Single bin histogram to know if random spins are being used or grabbing from database

  TGraph* mG_Entries = 0;                       ///< Graph for number of entries vs. Run Index
  
private:
  bool mInternalHists = false;                  ///< Boolean to keep track if mHists was added externally or an internal one was created
  std::vector<StFwdAnaVirtual*> mAnaList;       ///< List of analysis modules whose calls are executed sequentially

  ClassDef(StFwdAnaDataMaker, 1)
};

#endif

