/***************************************************************************
 *
 * StFttPointMaker.h
 *
 * Author: jdb 2021
 ***************************************************************************
 *
 * Description: StFttPointMaker - class to fill the points in StEvent
 *
 ***************************************************************************/
#ifndef STFTTPOINTMAKER_H
#define STFTTPOINTMAKER_H
#include "StMaker.h"
#include <vector>
#include <map>
#include "StFttDbMaker/StFttDb.h"


// class StFttDb;
class StEvent;
class StFttCollection;
class StFttCluster;
class StFttPoint;

class StFttPointMaker: public StMaker {

public:
    StFttPointMaker( const char* name = "stgcPoint" );

    ~StFttPointMaker();


    Int_t  Init();
    Int_t  InitRun( Int_t );
    Int_t  FinishRun( Int_t );
    Int_t  Finish();
    Int_t  Make();

private:
    void InjectTestData();
    void MakeLocalPoints(UChar_t Rob);
    void MakeGlobalPoints();
    bool GhostHitRejection_StripGroup(int row_x, int row_y, double x, double y);
    //return kTRUE when a diagoanl cluster was find to match (x,y)
    bool GhostHitRejection_DiagH(double x, double y, int Rob, int &i_cluster);
    bool GhostHitRejection_DiagV(double x, double y, int Rob, int &i_cluster);
    
    StEvent*             mEvent;
    StFttCollection*     mFttCollection;
    Bool_t               Debug_flag;
    Bool_t               mUseTestData;
    StFttDb*             mFttDb;
    std::vector<StFttPoint*> mFttPoint;

    // why can not using some thing like this:
    std::vector<StFttCluster *> clustersPerRob[StFttDb::nRob][StFttDb::nStripOrientations];//save the cluster for per quadrant
    
    inline bool is_Group1(int row_x, int row_y, double x, double y) const { return ( (StFttDb::X_StripGroupEdge[0] <= x && x <= StFttDb::X_StripGroupEdge[1]) && (StFttDb::Y_StripGroupEdge[0] <= y && y <= StFttDb::Y_StripGroupEdge[1]) && (row_x == 0) && (row_y == 0) ); }
    inline bool is_Group2(int row_x, int row_y, double x, double y) const { return ( (StFttDb::X_StripGroupEdge[1] <= x && x <= StFttDb::X_StripGroupEdge[4]) && (StFttDb::Y_StripGroupEdge[0] <= y && y <= StFttDb::Y_StripGroupEdge[1]) && (row_x == 0) && (row_y == 1)); }
    inline bool is_Group3(int row_x, int row_y, double x, double y) const { 
        return ( ( ( (StFttDb::X_StripGroupEdge[4] <= x && x <= StFttDb::X_StripGroupEdge[6]) && (StFttDb::Y_StripGroupEdge[0] <= y && y <= StFttDb::Y_StripGroupEdge[1]) ) || ( (StFttDb::X_StripGroupEdge[6]<= x && x <= StFttDb::X_StripGroupEdge[7]) && (StFttDb::Y_StripGroupEdge[0] <= y && y <= StFttDb::Y_StripGroupEdge[2]) ) ) && (row_x == 0) && (row_y == 2) );
    }
    inline bool is_Group4(int row_x, int row_y, double x, double y) const { return ( (StFttDb::X_StripGroupEdge[0] <= x && x <= StFttDb::X_StripGroupEdge[1]) && (StFttDb::Y_StripGroupEdge[1] <= y && y <= StFttDb::Y_StripGroupEdge[4]) && (row_x == 1) && (row_y == 0) ); }
    inline bool is_Group5(int row_x, int row_y, double x, double y) const {
        return ( ( ((StFttDb::X_StripGroupEdge[1] <= x && x <= StFttDb::X_StripGroupEdge[3]) && (StFttDb::Y_StripGroupEdge[1] <= y && y <= StFttDb::Y_StripGroupEdge[4])) || ((StFttDb::X_StripGroupEdge[3] <= x && x <= StFttDb::X_StripGroupEdge[4]) && (StFttDb::Y_StripGroupEdge[1] <= y && y <= StFttDb::Y_StripGroupEdge[5])) || ((StFttDb::Y_StripGroupEdge[4] <= x && x <= StFttDb::X_StripGroupEdge[5]) && (StFttDb::Y_StripGroupEdge[3] <= y && y <= StFttDb::Y_StripGroupEdge[5])) ) && (row_x == 1) && (row_y == 1) );
    }
    inline bool is_Group6(int row_x, int row_y, double x, double y) const { return ((StFttDb::X_StripGroupEdge[4] <= x && x <= StFttDb::X_StripGroupEdge[6]) && (StFttDb::Y_StripGroupEdge[1] <= y && y <= StFttDb::Y_StripGroupEdge[3])) && (row_x == 1) && (row_y == 2); }
    inline bool is_Group7(int row_x, int row_y, double x, double y) const { 
        return ( ( ( (StFttDb::Y_StripGroupEdge[4] <= y && y <= StFttDb::Y_StripGroupEdge[6]) && (StFttDb::X_StripGroupEdge[0] <= x && x <= StFttDb::X_StripGroupEdge[1]) ) || ( (StFttDb::Y_StripGroupEdge[6]<= y && y <= StFttDb::Y_StripGroupEdge[7]) && (StFttDb::X_StripGroupEdge[0] <= x && x <= StFttDb::X_StripGroupEdge[2]) ) ) && (row_x == 2) && (row_y == 0));
    }
    inline bool is_Group8(int row_x, int row_y, double x, double y) const { return ( ((StFttDb::Y_StripGroupEdge[4] <= y && y <= StFttDb::Y_StripGroupEdge[6]) && (StFttDb::X_StripGroupEdge[1] <= x && x <= StFttDb::X_StripGroupEdge[3])) && (row_x == 2) && (row_y == 1)); }

    ClassDef( StFttPointMaker, 0 )
};

#endif // STFTTPOINTMAKER_H
