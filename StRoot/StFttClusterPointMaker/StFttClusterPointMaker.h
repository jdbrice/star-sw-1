#ifndef STFTTCLUSTERPOINTMAKER_H
#define STFTTCLUSTERPOINTMAKER_H
#include "StMaker.h"
#include <vector>
#include <map>
#include "StFttDbMaker/StFttDb.h"


// class StFttDb;
class StEvent;
class StFttCollection;
class StFttCluster;
class StFttPoint;

class StFttClusterPointMaker: public StMaker {

public:
    StFttClusterPointMaker( const char* name = "stgcClusterPoint" );

    ~StFttClusterPointMaker();


    Int_t  Init();
    Int_t  InitRun( Int_t );
    Int_t  FinishRun( Int_t );
    Int_t  Finish();
    Int_t  Make();

    void setUseGeantData( bool useGeantData ) { mUseGeantData = useGeantData; }

    /** @brief Exchange which strip orientation measures X and which measures Y.
     *
     * The offline convention is that VERTICAL strips measure X and HORIZONTAL
     * strips measure Y -- the name refers to the direction the strips RUN, not
     * the coordinate they measure. The online/hardware map implies the
     * opposite assignment; as of 2026-09 the online side is believed correct.
     *
     * TEMPORARY DIAGNOSTIC. A real fix has to start from the geometry
     * description and be made consistently everywhere; this switch only
     * exchanges the assignment at point-building time so the hypothesis can be
     * tested on real data without touching StFttDb.
     *
     * NOTE, and it should be checked rather than assumed: an equivalent swap
     * was tried in 2026-07 by relabelling in StFttDb::getOrientation(), and at
     * 101-file statistics it LOWERED FTT matched-hit purity (14.9% -> 7.0%).
     * That test predates the VMM time-calibration and maxStripLength fixes, so
     * it is worth repeating, but it is evidence against this hypothesis and
     * should not be quietly forgotten. See proposal_next_step_20260725.txt.
     *  @param apply : true to exchange X and Y
    */
    void setApplyXYMirror( bool apply ) { mApplyXYMirror = apply; }

    /** @brief Enable this maker's own verbose point dump (mDebug).
     * Distinct from StMaker::SetDebug(); mDebug was previously fixed to false
     * in the constructor with no way to turn it on, which made the point
     * coordinates impossible to inspect without editing the source.
    */
    void setPointDebug( bool d ) { mDebug = d; }

private:
    void InjectTestData();
    void MakeLocalPoints(UChar_t Rob);
    void MakeGlobalPoints();
    void MakeGeantPoints(); // act as a slow-sim
    //in the ClusterPointMaker there no ghost hit rejection since all the clusters will be saved

    StEvent*             mEvent;
    StFttCollection*     mFttCollection;
    Bool_t               mDebug;
    Bool_t               mUseTestData;
    Bool_t               mUseGeantData; // if true, use the geant hits to make points
    Bool_t               mApplyXYMirror = kFALSE; // see setApplyXYMirror()
    StFttDb*             mFttDb;
    std::vector<StFttCluster *> clustersPerRob[StFttDb::nRob][StFttDb::nStripOrientations]; //save the cluster for per quadrant

    inline bool is_Group1(int row_x, int row_y, double x, double y) const { return ( (14.60 <= x && x <= 172.29) && (14.60 <= y && y <= 172.29) && (row_x == 0) && (row_y == 0) ); }
    inline bool is_Group2(int row_x, int row_y, double x, double y) const { return ( (172.29 <= x && x <= 360.09) && (14.60 <= y && y <= 172.29) && (row_x == 0) && (row_y == 1)); }
    inline bool is_Group3(int row_x, int row_y, double x, double y) const
    {
        return ( ( ( (360.09 <= x && x <= 504.2) && (14.60 <= y && y <= 172.29) ) || ( (504.2<= x && x <= 548.3) && (14.60 <= y && y <= 216.89) ) ) && (row_x == 0) && (row_y == 2) );
    }
    inline bool is_Group4(int row_x, int row_y, double x, double y) const { return ( (14.60 <= x && x <= 172.29) && (172.29 <= y && y <= 360.09) && (row_x == 1) && (row_y == 0) ); }
    inline bool is_Group5(int row_x, int row_y, double x, double y) const
    {
        return ( ( ((172.29 <= x && x <= 315.4) && (172.29 <= y && y <= 360.09)) || ((315.4 <= x && x <= 360.09) && (172.29 <= y && y <= 410.9)) || ((360.09 <= x && x <= 410.9) && (315.4 <= y && y <= 410.9)) ) && (row_x == 1) && (row_y == 1) );
    }
    inline bool is_Group6(int row_x, int row_y, double x, double y) const { return ((360.09 <= x && x <= 504.2) && (172.29 <= y && y <= 315.4)) && (row_x == 1) && (row_y == 2); }
    inline bool is_Group7(int row_x, int row_y, double x, double y) const
    {
        return ( ( ( (360.09 <= y && y <= 504.2) && (14.60 <= x && x <= 172.29) ) || ( (504.2<= y && y <= 548.3) && (14.60 <= x && x <= 216.89) ) ) && (row_x == 2) && (row_y == 0));
    }
    inline bool is_Group8(int row_x, int row_y, double x, double y) const { return ( ((360.09 <= y && y <= 504.2) && (172.29 <= x && x <= 315.4)) && (row_x == 2) && (row_y == 1)); }

    ClassDef( StFttClusterPointMaker, 1 )
};

#endif