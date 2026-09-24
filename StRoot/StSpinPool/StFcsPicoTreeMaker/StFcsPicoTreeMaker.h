/*
  The purpose of this StMaker is to process a file and generate a StFcsPicoTree with StFcsPicoHit/Cluster/Point information. Those classes can be found in StRoot/StSpinPool/StFcsTreeManager. The reason for keeping them seperate is so that STAR libraries don't need to be loaded in order to read the tree
  
  @[April 3, 2023](David Kapukchyan) > First instance
 */

#ifndef StFcsPicoTreeMaker_H
#define StFcsPicoTreeMaker_H


//ROOT headers
#include "TFile.h"
#include "TString.h"

//STAR headers
#include "StMaker.h"

//Custom header
#include "StSpinPool/StFcsTreeManager/StFcsPicoTree.h"

//These classes need to be forward declared since including the headers doesn't compile right
class StFcsDb;
class StFcsCollection;

class StFcsPicoTreeMaker : public StMaker
{
  
 public:
  StFcsPicoTreeMaker(const char* name = "StFcsPicoTreeMaker");
  virtual ~StFcsPicoTreeMaker();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();

  const char* FileName(){ return mFileName.Data(); }
  void setFileName(const char* filename) { mFileName = filename; }
  StFcsPicoTree* setTree(const char* treename, const char* title="",Int_t splitlevel=99);
  
  void turnOffHits();
  void turnOffClusters();
  void turnOffPoints();
  
protected:
  TString mFileName = "";       //!< Used for input and output
  
  StFcsDb* mFcsDb = 0;
  StFcsCollection* mFcsColl = 0;
  
  StFcsPicoTree* mDataTree = 0;
  bool mDoHits = true;        //! Used to turn off filling tree with FCS hits
  bool mDoClus = true;        //! Used to turn off filling tree with FCS clusters
  bool mDoPoints = true;      //! Used to turn off filling tree with FCS points
  
private:
  TFile* mFile = 0;
  
  ClassDef(StFcsPicoTreeMaker,1);
};

#endif
