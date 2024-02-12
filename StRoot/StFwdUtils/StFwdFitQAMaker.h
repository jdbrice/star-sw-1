#ifndef ST_FWD_FIT_QA_MAKER_H
#define ST_FWD_FIT_QA_MAKER_H

#include <map>

#include "StChain/StMaker.h"
#include "TVector3.h"


class g2t_track_st;
class McFwdTrack {
  public:
  TVector3 p;
  int pid = -1;
  int id = -1;
  int q = 0;
  int numFST = 0;
  int numFTT = 0;
  g2t_track_st *g2track = nullptr;
};


class StFwdFitQAMaker : public StMaker
{
  public:
    StFwdFitQAMaker();
    ~StFwdFitQAMaker(){/* nada */};

    int Init();
    int Finish();
    int Make();
    void Clear(const Option_t *opts = "");
    void ProcessFwdTracks();
    void ProcessFwdMuTracks();

    // StEvent analyzed by default
    // call this to analyze the MuDst instead
    void setMuDstInput() { mAnalyzeMuDst = true; }

  protected:

    /**
     * @brief Map of <name (TString), histogram>
     * 
     */
    std::map<TString, TH1*> mHists;

    /**
     * @brief Add a histogram to the map
     * 
     * Convenience method to avoid duplicate typing name / possible mismatch
     * @param h1 TH1* histogram object
     * @return TH1* return the same object for ease of use
     */
    TH1 * addHist( TH1 * h1){
      mHists[h1->GetName()] = h1;
      return h1;
    }

    /**
     * @brief Get the Hist object from the map
     *  - Additional check and safety for missing histograms
     * @param n Histogram name
     * @return TH1* histogram if found, otherwise a 'nil' histogram with one bin
     */
    TH1* getHist( TString n ){
      if (mHists.count(n))
        return mHists[n];
      LOG_ERROR << "Attempting to access non-existing histogram" << endm;
      return new TH1F( "NULL", "NULL", 1, 0, 1 ); // returning nullptr can lead to seg fault, this fails softly
    }

    /**
     * @brief Control whether the analysis uses StEvent (default) or MuDst as input
     * 
     */
    bool mAnalyzeMuDst = false;
    vector<McFwdTrack> mcTracks;

  ClassDef(StFwdFitQAMaker, 0);
};

#endif
