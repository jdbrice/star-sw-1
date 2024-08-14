/***************************************************************************
 *
 * StFttPointMaker.cxx
 *
 * Author: jdb & Zhen 2024
 ***************************************************************************
 *
 * Description: StFttPointMaker - class to fill the StFttPoint in StEvent
 * Note: in the cluster and point local position, the unit is mm, for the 
 *       global position, unit is cm
 ***************************************************************************/
#include <vector>
#include <map>
#include <array>
#include <algorithm>


#include "StEvent.h"
#include "StEnumerations.h"

#include "StFttPointMaker.h"


#include "StEvent/StFttRawHit.h"
#include "StEvent/StFttCluster.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"

#include "StFttDbMaker/StFttDb.h"


//_____________________________________________________________
StFttPointMaker::StFttPointMaker( const char* name )
: StMaker( name ),
  mEvent( 0 ),          /// pointer to StEvent
  Debug_flag( false ),       /// print out of all full messages for debugging
//   Debug_flag( true ),       /// print out of all full messages for debugging
  mUseTestData( false ),
  mFttDb( nullptr )
{
    LOG_DEBUG << "StFttPointMaker::ctor"  << endm;
    LOG_INFO << "******** StFttPointMaker::StFttPointMaker = "<<name<<endm;
}

//_____________________________________________________________
StFttPointMaker::~StFttPointMaker()
{  /* no op */ }

//_____________________________________________________________
Int_t
StFttPointMaker::Init()
{
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttPointMaker::InitRun( Int_t runnumber )
{ 
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttPointMaker::FinishRun( Int_t runnumber )
{ 
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttPointMaker::Finish()
{ 
    return kStOk;
}


//_____________________________________________________________
Int_t
StFttPointMaker::Make()
{ 
    LOG_DEBUG << "StFttPointMaker::Make()" << endm;

    mEvent = (StEvent*)GetInputDS("StEvent");
    if(mEvent) {
        LOG_DEBUG<<"Found StEvent"<<endm;
    } else {
        return kStOk;
    }
    mFttCollection=mEvent->fttCollection();
    if(!mFttCollection) {
        return kStOk;
    }

    mFttDb = static_cast<StFttDb*>(GetDataSet("fttDb"));

    assert( mFttDb );

    if ( mUseTestData )
        InjectTestData();

    // next we will need them in even more detail
    // per strip group, but start here
    // Zhen: remove the defination to the .h files
    for (int i = 0; i<16; i++)
    {
        for (int j = 0; j<4; j++)
        {
            clustersPerRob[i][j].clear();
        }
    }

    for ( StFttCluster* clu : mFttCollection->clusters() ) 
    {
        // group clusters by quadrant, hor, vert, hdiag, vdiag
       

        UChar_t rob = mFttDb->rob( clu );
        if ( clu->nStrips() < 2 ) continue;// add cluster width limit
        if(Debug_flag)
        {
            LOG_INFO << "rob = " << (int)rob << endm;
            LOG_INFO << "direction = " << (int)clu->orientation() << endm;
            LOG_INFO << "cluster x = " << clu->x() << endm;
        }
        clustersPerRob[ (int)rob ][ clu->orientation() ].push_back( clu );// clustersPerRob[ rob ] [ orientation ]
        // clustersPerRob[ clu->orientation() ] [ rob ].push_back( clu ); 
    } // loop on hit

    for (int i = 0; i < 16; i++) 
    {
        if (Debug_flag) 
        {
            LOG_INFO << "nCluster kFttVertical = " << clustersPerRob[ i ][ kFttVertical ].size() << endm;
            LOG_INFO << "nCluster kFttHorizontal = " << clustersPerRob[ i ][ kFttHorizontal ].size() << endm;
            LOG_INFO << "nCluster kFttDiagonalV = " << clustersPerRob[ i ][ kFttDiagonalV ].size() << endm;
            LOG_INFO << "nCluster kFttDiagonalH = " << clustersPerRob[ i ][ kFttDiagonalH ].size() << endm;
        }

        MakeLocalPoints((UChar_t)i); // make local points for each Quad
    }
    MakeGlobalPoints();
    if(Debug_flag) LOG_INFO << "StFttPointMaker made " << mFttCollection->numberOfPoints() << " points this event" << endm;

    return kStOk;
}

void StFttPointMaker::InjectTestData(){
    mFttCollection->rawHits().clear();

    // TODO: inject clean strip hits to test cluster finder
    // StFttRawHit *hit = new StFttRawHit( sec, rdo, feb, vm, vmm[0].ch, vmm[0].adc, vmm[0].bcid, vmm[0].tb );
    // hit->setMapping( plane, quadrant, row, strip ) 
}

void StFttPointMaker::MakeGlobalPoints() {
    for ( StFttPoint * p : mFttCollection->points() ){

        float x = p->x();
        float y = p->y();
        float z = 0.0;
        StThreeVectorD global;

        // dx is a local shift
        float dx = 0, dy = 0, dz = 0;
        // sx is only {1,-1} -> reflected or normal
        float sx = 0, sy = 0, sz = 0;
        mFttDb->getGloablOffset( p->plane(), p->quadrant(), dx, sx, dy, sy, dz, sz );
        global.set( (x) * sx +dx, (y) * sy +dy, (z + dz) * sz );
        //convert from mm to the cm, just divided by 10, 24.08.07
        // global.set((global.x()),(global.y()),(global.z())); 
        global.set((global.x()/10.),(global.y()/10.),(global.z())); 
        p->setXYZ( global );
    }
}

// using strip group method to reject ghost hit
bool StFttPointMaker::GhostHitRejection_StripGroup( int row_x, int row_y, double x, double y) 
{
    // check all the strip groups to reject the ghost
    bool is_ghosthit = kFALSE;
    if ( is_Group1(row_x,row_y,x,y) || is_Group2(row_x,row_y,x,y) || is_Group3(row_x,row_y,x,y) || is_Group4(row_x,row_y,x,y) || is_Group5(row_x,row_y,x,y) || is_Group6(row_x,row_y,x,y) || is_Group6(row_x,row_y,x,y) || is_Group7(row_x,row_y,x,y) || is_Group8(row_x,row_y,x,y) ) 
    is_ghosthit = kTRUE;
    // is_ghosthit = is_Group1(x,y);
    // is_ghosthit = is_Group2(x,y);
    // is_ghosthit = is_Group3(x,y);
    // is_ghosthit = is_Group4(x,y);
    // is_ghosthit = is_Group5(x,y);
    // is_ghosthit = is_Group6(x,y);
    // is_ghosthit = is_Group7(x,y);
    // is_ghosthit = is_Group8(x,y);

    return is_ghosthit;
}

//using diagnoal horizontal strip to reject method, 
bool StFttPointMaker::GhostHitRejection_DiagH(double x, double y, int Rob, int &i_cluster)
{
    // TODO: how to confirm i_cluster
    bool is_pair = kFALSE;
    i_cluster = -999;
    //loop the diagonal cluster find a cluster can include this cluster
    // for (StFttCluster* clu_dx : *clustersPerRob[(UChar_t)Rob][kFttDiagonalH])
    size_t nclusters = clustersPerRob[(UChar_t)Rob][kFttDiagonalH].size();
    if(Debug_flag) LOG_INFO << "nDV clusters = " << nclusters << endm;
    double distance = -99.;
    double distance_prev = 999.;
    for (size_t iClu_DH = 0; iClu_DH<nclusters; iClu_DH++)
    {
        double intercept = 0;
        intercept = x+y;
        auto clu_dx=clustersPerRob[(UChar_t)Rob][kFttDiagonalH][iClu_DH];

        double LEdge = clu_dx->maxStripLeftEdge()*sqrt(2);
        double REdge = clu_dx->maxStripRightEdge()*sqrt(2);

        //for loose cluster cut and check performance, this number now select by hand
        // LEdge = LEdge-1.6*sqrt(2);
        // REdge = REdge+1.6*sqrt(2);

        if(Debug_flag)
        {
            LOG_INFO << "intercept = " << intercept << " LEdge = " << LEdge << " REdge = " << REdge << endm;
            LOG_INFO << "cluster x = " << clu_dx->x()*sqrt(2) << endm;
        }

        if(intercept >= LEdge && intercept <= REdge) 
        {
            is_pair = kTRUE;
            distance = abs(intercept-clu_dx->x()*sqrt(2));
            if (distance < distance_prev)
            {
                distance_prev = distance;
                i_cluster = iClu_DH;
            }
        }
        // LOG_INFO << "i_cluster = " << i_cluster << endm;
        // if(clu_dx->x()*sqrt(2)+1.60*3 > intercept && clu_dx->x()*sqrt(2)-1.60*3 < intercept)
        // {
        //     is_pair = kTRUE;
        //     distance = abs(intercept-clu_dx->x()*sqrt(2));
        //     if (distance < distance_prev)
        //     {
        //         distance_prev = distance;
        //         i_cluster = iClu_DH;
        //     }
        // }
        // LOG_INFO << "i_cluster = " << i_cluster << endm;
    }
    // LOG_INFO << "is pair = " << (int)is_pair << endm;
    return is_pair;
}

//using diagnoal vertical strip to reject method, 
bool StFttPointMaker::GhostHitRejection_DiagV(double x, double y, int Rob, int &i_cluster)
{
    if(Debug_flag) LOG_INFO << "starting do GhostHitRejection_DiagV " << endm;
    // TODO: how to confirm i_cluster
    bool is_pair = kFALSE;
    i_cluster = -999;
    //loop the diagonal cluster find a cluster can include this cluster
    size_t nclusters = clustersPerRob[(UChar_t)Rob][kFttDiagonalV].size();
    double distance = -99.;
    double distance_prev = 999.;
    for (size_t iClu_DV = 0; iClu_DV<nclusters; iClu_DV++)
    {
        double intercept = 0;
        intercept = x+y;
        auto clu_dx=clustersPerRob[(UChar_t)Rob][kFttDiagonalV][iClu_DV];

        double LEdge = clu_dx->maxStripLeftEdge()*sqrt(2);
        double REdge = clu_dx->maxStripRightEdge()*sqrt(2);
 
        //for loose cluster cut and check performance, this number now select by hand
                //  LEdge = LEdge-1.6*sqrt(2);
                //  REdge = REdge+1.6*sqrt(2);
          

        if(Debug_flag)
        {
            LOG_INFO << "intercept = " << intercept << " LEdge = " << LEdge << " REdge = " << REdge << endm;
            LOG_INFO << "cluster x = " << clu_dx->x()*sqrt(2) << endm;
        }

        //Method 1
        if(intercept >= LEdge && intercept <= REdge) 
        {
            is_pair = kTRUE;
            distance = abs(intercept-clu_dx->x()*sqrt(2));
            if (distance < distance_prev)
            {
                distance_prev = distance;
                i_cluster = iClu_DV;
            }
        }
        // LOG_INFO << "i_cluster = " << i_cluster << endm;
        //Method 2
        // if(clu_dx->x()*sqrt(2) < intercept+1.60*3 && clu_dx->x()*sqrt(2) > intercept-1.60*3)
        // {
        //     is_pair = kTRUE;
        //     distance = abs(intercept-clu_dx->x()*sqrt(2));
        //     if (distance < distance_prev)
        //     {
        //         distance_prev = distance;
        //         i_cluster = iClu_DV;
        //     }
        // }
        // LOG_INFO << "i_cluster = " << i_cluster << endm;
    }
    // LOG_INFO << "is pair = " << (int)is_pair << endm;
    return is_pair;
}

//--------------------------------------------------------------
//for  the loacl coordinate, if using the the center of pin hole as (0,0)
//center of first strip of V&H strips is 15.95mm
//center of first strip of dia strips is 19.42mm
void StFttPointMaker::MakeLocalPoints(UChar_t Rob)
{
    StFttPoint* point;
    double x = -999.;
    double y = -999.;
    size_t nClusters_X = 0;size_t nClusters_Y = 0;size_t nClusters_DX = 0;size_t nClusters_DY = 0;
    nClusters_X = clustersPerRob[(UChar_t)Rob][kFttVertical].size();
    nClusters_Y = clustersPerRob[(UChar_t)Rob][kFttHorizontal].size();

    if(Debug_flag)
    {
        LOG_INFO << "rob = " << (int)Rob << endm;
        LOG_INFO << "nClusterX = " << nClusters_X << " nClusterY = " << nClusters_Y << endm;
        LOG_INFO << "nCluster dV = " << clustersPerRob[(UChar_t)Rob][kFttDiagonalV].size()<< endm;
        LOG_INFO << "nCluster dH = " << clustersPerRob[(UChar_t)Rob][kFttDiagonalH].size()<< endm;
    }

    
    for ( size_t iClu_X = 0; iClu_X < nClusters_X; iClu_X++ )
    {
        // point = new StFttPoint();
        auto clu_x = clustersPerRob[(UChar_t)Rob][kFttVertical][iClu_X];
        x = clu_x->x();
        int Row_x = clu_x->row();

        if(Debug_flag)
        {
            LOG_INFO << "x cluster plane = " << (int)clu_x->plane() << " quad = " << (int)clu_x->quadrant() << endm;
            LOG_INFO << "start x loop, x = " << x << endm;
        }

        for ( size_t iClu_Y = 0; iClu_Y < nClusters_Y; iClu_Y++ )
        {
            point = new StFttPoint();
            auto clu_y = clustersPerRob[(UChar_t)Rob][kFttHorizontal][iClu_Y];
            y = clu_y->x();
            int Row_y = clu_y->row();
            if(Debug_flag)
            {
                if (iClu_Y == 0) LOG_INFO << "y cluster plane = " << (int)clu_y->plane() << " quad = " << (int)clu_y->quadrant() << endm;
                LOG_INFO << "start y loop, y = " << y << endm;
            }

            // get the x-y pair and check the region
            //using strip group to reject ghost hit
            if ( !GhostHitRejection_StripGroup(Row_x,Row_y,x,y))
                continue;
            
            if ( x < 1.e-5 || y < 1.e-5) continue;//remove the cluster with wrong number
            if ( x > 1.e5 || y > 1.e5) continue;//remove the cluster with wrong number

            if (x>y)//for the diagonal_V is fully cover the lower half 
            {
                bool is_pair = kFALSE;
                int i_cluster = -1;//the the index of d_V cluster
                if(Debug_flag)
                {
                    LOG_INFO << "x>y" << endm;
                }
                if( GhostHitRejection_DiagV(x,y,Rob,i_cluster) ) 
                {
                    // LOG_INFO << "debug Diag V1" << endm;
                    // point->setX(x/10.);// cm unit. 24.08.07
                    point->setX(x);
                    point->setSigmaX(clu_x->sigma());
                    // point->setY(y/10.);// cm unit. 24.08.07
                    point->setY(y);
                    point->setSigmaY(clu_y->sigma());
                    if(Debug_flag) LOG_INFO << "Point (X,Y) = " << x << ", " << y << endm;
                    point->setPlane(clu_x->plane());
                    point->setQuadrant(clu_x->quadrant());
                    point->addCluster(clu_x,kFttVertical);
                    point->addCluster(clu_y,kFttHorizontal);
                    auto clu_dv =  clustersPerRob[(UChar_t)Rob][kFttDiagonalV][i_cluster];
                    // if( ((x+y) - clu_dv->x()*sqrt(2)) < 10 && ((x+y) - clu_dv->x()*sqrt(2)) > 6 )
                    // {
                        // // LOG_INFO << "x = " << x << " y = " << y <<"cluster center = " << clu_dv->x() <<" cluster x * sqrt(2) = " << clu_dv->x()*sqrt(2) << " MaxADC strip center = " << clu_dv->SC_MaxStrip() << endm;
                    // }
                    // LOG_INFO << "match diagonal cluster is " << i_cluster << " DV cluster" << endm;
                    clu_dv->print();
                    point->addCluster(clu_dv,kFttDiagonalV);
                    // point->setD1(clu_dv->x()/10.); // cm unit. 24.08.07
                    point->setD1(clu_dv->x());
                    is_pair = kTRUE;
                    // LOG_INFO << "debug Diag V2" << endm;
                    if ( GhostHitRejection_DiagH(x,y,Rob,i_cluster) )
                    {
                        auto clu_dh =  clustersPerRob[(UChar_t)Rob][kFttDiagonalH][i_cluster];
                        // point->setD2(clu_dh->x()/10.); // cm unit. 24.08.07
                        point->setD2(clu_dh->x());
                        // LOG_INFO << "debug Diag H1" << endm;
                        point->addCluster(clu_dh,kFttDiagonalH);
                    }
                    
                } else if ( GhostHitRejection_DiagH(x,y,Rob,i_cluster) )
                {

                    // LOG_INFO << "debug Diag H2" << endm;
                    is_pair = kTRUE;
                    if(Debug_flag) LOG_INFO << "Point (X,Y) = " << x << ", " << y << endm;
                    // point->setX(x/10.);// cm unit. 24.08.07
                    // point->setY(y/10.);
                    point->setX(x);
                    point->setY(y);
                    point->setSigmaX(clu_x->sigma());
                    point->setSigmaY(clu_y->sigma());
                    point->setPlane(clu_x->plane());
                    point->setQuadrant(clu_x->quadrant());
                    // LOG_INFO << "debug Diag H3" << endm;
                    point->addCluster(clu_x,kFttVertical);
                    // LOG_INFO << "debug Diag H4" << endm;
                    point->addCluster(clu_y,kFttHorizontal);
                    if(Debug_flag) LOG_INFO << "match diagonal cluster is " << i_cluster << " DH cluster" << endm;
                    auto clu_dh =  clustersPerRob[(UChar_t)Rob][kFttDiagonalH][i_cluster];
                    clu_dh->print();
                    // point->setD1(clu_dh->x()/10.);// cm unit. 24.08.07
                    point->setD1(clu_dh->x());
                    point->addCluster(clu_dh,kFttDiagonalH);
                }
                else 
                    continue;

                // LOG_INFO << "is pair = " << is_pair << endm;
                if(is_pair)
                {
                    point->print();
                    LOG_INFO << "PAIR Point (X,Y) = " << point->x() << ", " << point->y() << endm;
                    mFttPoint.push_back(point);
                    mFttCollection->addPoint(point);
                }
            }
            if (x<y)//for the diagonal_H is fully cover the higher half 
            {
                if(Debug_flag)
                {
                    LOG_INFO << "x<y" << endm;
                }
                bool is_pair = kFALSE;
                int i_cluster = -1;//the the index of d_V cluster
                if( GhostHitRejection_DiagH(x,y,Rob,i_cluster) ) 
                {
                    is_pair = kTRUE;
                    // point->setX(x/10.);
                    // point->setY(y/10.);// cm unit. 24.08.07
                    point->setX(x);
                    point->setY(y);
                    point->setSigmaX(clu_x->sigma());
                    point->setSigmaY(clu_y->sigma());
                    point->setPlane(clu_x->plane());
                    point->setQuadrant(clu_x->quadrant());
                    point->addCluster(clu_x,kFttVertical);
                    point->addCluster(clu_y,kFttHorizontal);
                    auto clu_dh =  clustersPerRob[(UChar_t)Rob][kFttDiagonalH][i_cluster];
                    // if( ((x+y) - clu_dh->x()*sqrt(2)) < 10 && ((x+y) - clu_dh->x()*sqrt(2)) > 6 )
                    // {
                        // // LOG_INFO << "x = " << x << " y = " << y << "cluster center = " << clu_dh->x() << " cluster x * sqrt(2) = " << clu_dh->x()*sqrt(2) << " MaxADC strip center = " << clu_dh->SC_MaxStrip() << endm;
                    // }
                    // LOG_INFO << "match diagonal cluster is " << i_cluster << " DH cluster" << endm;
                    clu_dh->print();
                    // point->setD1(clu_dh->x()/10.);// cm unit. 24.08.07
                    point->setD1(clu_dh->x());
                    point->addCluster(clu_dh,kFttDiagonalH);
                    if (GhostHitRejection_DiagV(x,y,Rob,i_cluster))
                    {
                        auto clu_dv =  clustersPerRob[(UChar_t)Rob][kFttDiagonalV][i_cluster];
                        // point->setD2(clu_dv->x()/10.);// cm unit. 24.08.07
                        point->setD2(clu_dv->x());
                        point->addCluster(clu_dv,kFttDiagonalV);
                    }
                } else if ( GhostHitRejection_DiagV(x,y,Rob,i_cluster) )
                {
                    is_pair = kTRUE;
                    // point->setX(x/10.);// cm unit. 24.08.07
                    // point->setY(y/10.);
                    point->setX(x);
                    point->setY(y);
                    point->setSigmaX(clu_x->sigma());
                    point->setSigmaY(clu_y->sigma());
                    point->setPlane(clu_x->plane());
                    point->setQuadrant(clu_x->quadrant());
                    point->addCluster(clu_x,kFttVertical);
                    point->addCluster(clu_y,kFttHorizontal);
                    if(Debug_flag) LOG_INFO << "match diagonal cluster is " << i_cluster << " DV cluster" << endm;
                    auto clu_dv =  clustersPerRob[(UChar_t)Rob][kFttDiagonalV][i_cluster];
                    clu_dv->print();
                    point->setD1(clu_dv->x()/10.);
                    // point->setD1(clu_dv->x());
                    point->addCluster(clu_dv,kFttDiagonalV);
                } else continue;

                // LOG_INFO << "is pair = " << is_pair << endm;nn
                if(is_pair)
                {
                    if(Debug_flag) LOG_INFO << "PAIR Point (X,Y) = " << point->x() << ", " << point->y() << endm;
                    // point->print();
                    mFttPoint.push_back(point);
                    mFttCollection->addPoint(point);
                }
            }
            if (!point)
            {
                if(Debug_flag) LOG_INFO << "empty point !!!!!!" << endm;
                continue;
            }
            // point->print();
            // delete point;
        }
    }
    
}


