
#ifndef ST_FTT_SLOW_SIM_MAKER_H
#define ST_FTT_SLOW_SIM_MAKER_H

class g2t_emc_hit_st;
class StFtsHit;
class StEvent;

#include "StChain/StMaker.h"
#include <vector>
#include <map>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

class StFttDb;
class StFttCollection;
class StRnDHit;


const Int_t FTT_MAX_HITS = 50000;
const Int_t FTT_MAX_CLUSTERS = 10000;
const Int_t FTT_MAX_POINTS = 4000;
struct FttMCData
{
    // event information
    Int_t    EVT;
    UShort_t N;

    //raw hit information
    UChar_t    sec[FTT_MAX_HITS];
    UChar_t    rdo[FTT_MAX_HITS];
    UChar_t    plane[FTT_MAX_HITS];
    UChar_t    quad[FTT_MAX_HITS];
    UChar_t    feb[FTT_MAX_HITS];
    UChar_t    febvmm[FTT_MAX_HITS];
    UChar_t    vmm[FTT_MAX_HITS];
    UChar_t    ch[FTT_MAX_HITS];
    UShort_t   bcid[FTT_MAX_HITS];
    Short_t    dbcid[FTT_MAX_HITS];
    Short_t    time[FTT_MAX_HITS];
    UShort_t   adc[FTT_MAX_HITS];
    Short_t    tb[FTT_MAX_HITS];
    UChar_t    row[FTT_MAX_HITS];
    UChar_t    strip[FTT_MAX_HITS];
    UChar_t    dir[FTT_MAX_HITS];
    
};

struct StFttSlowSimTreeData {
  // GEANT info
  vector<double> px, py, pz, ds, de, tof, pt, eta, phi;
  vector<int> vid, trackp;
  vector<short> plane, quad;
  vector<Int_t> rdo, feb, vmm, ch, bcid, dbcid, tb, ADC, row, strip, dir;
  vector<double> x, y, z;

  void clear() {
    x.clear(); y.clear(); z.clear();
    px.clear(); py.clear(); pz.clear();
    ds.clear(); de.clear(); tof.clear();
    vid.clear();
    plane.clear(); quad.clear(); 
  }
};

class StFttSlowSimMaker : public StMaker {
  public:
    explicit StFttSlowSimMaker(const Char_t *name = "fttSlowSim");
    virtual ~StFttSlowSimMaker() {}
    Int_t Make();
    int Init();
    
    int Finish();


  private:
    int GetQuad(double x_global, double y_global);
    void Global2Local_2D(double &x_local, double &y_local, double x_global, double y_global, int i_plane, int i_quad);// transform the global mc information to the local "signal like" information
    bool GetStripandRow_XY(float x_local, float y_local, int &strip_x, int &strip_y, int &row_x, int &row_y);
    void GetStripandRow_Diag(float x_local, float y_local, int &strip_dv, int &strip_dh, int &row_dv, int &row_dh);
    Int_t SampleCluster(int center_strip, double xhit_local, int row_x, int sec, int rdo, int is_diag); // 0 for XY and 1 for diag
    Int_t SampleCluster(int center_strip, double xhit_local, int row_x, int sec, int rdo, int is_diag, int i_evt); // 0 for XY and 1 for diag, for QA
    void FillThinGapChambers(StEvent *event);

    int iEvent;

    const double WIRE_WIDTH = 0.32;    // cm
    const double CLUSTER_SIGMA_STRIP = 0.01;       // 100 microns

    int sTGCNRealPoints = 0;
    int sTGCNGhostPoints = 0;

    // for sample MC points
    TF1* fClusterProfile;
    TF1* fClusterWidth;

    TH2 * hXY;

    TFile *mTreeFile = nullptr;
    TTree *mTree     = nullptr;
    StFttSlowSimTreeData mTreeData;

    StEvent*             mEvent;
    StFttCollection*     mFttCollection;
    StFttDb*             mFttDb;

    Bool_t  mDebug = kTRUE;//for debug
    // Bool_t  mDebug = kFALSE;

    //for QA
    TFile* mFile;
    std::map< string, TH1* >    mH1d;
    std::map< string, TH2* >    mH2d;
    void BookHistograms();
    void WriteHistograms(); // how to create a root files?

    //function to get the row of MC hit
    bool is_Group1(int &row_x, int &row_y, double x, double y);
    bool is_Group2(int &row_x, int &row_y, double x, double y);
    bool is_Group3(int &row_x, int &row_y, double x, double y);
    bool is_Group4(int &row_x, int &row_y, double x, double y);
    bool is_Group5(int &row_x, int &row_y, double x, double y);
    bool is_Group6(int &row_x, int &row_y, double x, double y);
    bool is_Group7(int &row_x, int &row_y, double x, double y);
    bool is_Group8(int &row_x, int &row_y, double x, double y);

    ClassDef(StFttSlowSimMaker, 0)
};


#endif
