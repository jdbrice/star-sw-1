#include <vector>
#include <map>
#include <array>
#include <algorithm>

#include "StEvent.h"
#include "StEnumerations.h"

#include "StFttClusterPointMaker.h"

#include "StEvent/StFttRawHit.h"
#include "StEvent/StFttCluster.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFttCollection.h"

#include "StFttDbMaker/StFttDb.h"

StFttClusterPointMaker::StFttClusterPointMaker( const char* name )
: StMaker( name ),
  mEvent( 0 ),          /// pointer to StEvent
  mDebug( true ),       /// print out of all full messages for debugging
  mUseTestData( false ),
  mFttDb( nullptr )
{
    LOG_DEBUG << "StFttClusterPointMaker::ctor"  << endm;
}

//_____________________________________________________________
StFttClusterPointMaker::~StFttClusterPointMaker()
{  /* no op */ }

//_____________________________________________________________
Int_t
StFttClusterPointMaker::Init()
{
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttClusterPointMaker::InitRun( Int_t runnumber )
{
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttClusterPointMaker::FinishRun( Int_t runnumber )
{
    return kStOk;
}

//_____________________________________________________________
Int_t
StFttClusterPointMaker::Finish()
{
    return kStOk;
}

Int_t StFttClusterPointMaker::Make() {
    LOG_DEBUG << "StFttClusterPointMaker::Make()" << endm;

    mEvent = (StEvent*)GetInputDS("StEvent"); //get the event from stevent and make sure it exists
    if(mEvent) {
        if(mDebug){ LOG_DEBUG<<"Found StEvent"<<endm; }
    } else {
        return kStOk;
    }

    mFttCollection=mEvent->fttCollection(); //get the ftt collection from the event and make sure it exists
    if(!mFttCollection) { return kStOk; }

    mFttDb = static_cast<StFttDb*>(GetDataSet("fttDb")); //get the ftt database
    assert(mFttDb);

    //clear cluster vector
    for (int i = 0; i<16; i++) {
        for (int j = 0; j<4; j++) {
            clustersPerRob[i][j].clear();
        }
    }

    //organize clusters by rob and orientation
    for ( StFttCluster* clu : mFttCollection->clusters() ) {
        // group clusters by rob (1-16, 4 quadrants of 4 sTGC planes) and strip orientation
        UChar_t rob = mFttDb->rob( clu );
        if (mDebug){
            LOG_INFO << "rob = " << (int)rob << endm;
            LOG_INFO << "direction = " << (int)clu->orientation() << endm;
            LOG_INFO << "cluster x = " << clu->x() << endm;
        }
        clustersPerRob[ (int)rob ][ clu->orientation() ].push_back( clu );
    } // loop on cluster

    //loop over rob and make local points
    for (int i = 0; i < 16; i++) {
        if (mDebug){
            LOG_INFO << "Now at ROB " << i << endm;
            LOG_INFO << "nCluster kFttVertical = " << clustersPerRob[ i ][ kFttVertical ].size() << endm; //vertical is vertical strips, clusters_x
            LOG_INFO << "nCluster kFttHorizontal = " << clustersPerRob[ i ][ kFttHorizontal ].size() << endm;
            LOG_INFO << "nCluster kFttDiagonalV = " << clustersPerRob[ i ][ kFttDiagonalV ].size() << endm;
            LOG_INFO << "nCluster kFttDiagonalH = " << clustersPerRob[ i ][ kFttDiagonalH ].size() << endm;
        }
        MakeLocalPoints((UChar_t)i);
    }

    //make the global points
    MakeGlobalPoints();
    LOG_INFO << "StFttClusterPointMaker made " << mFttCollection->numberOfPoints() << " points this event" << endm;

    return kStOk;
}

void StFttClusterPointMaker::InjectTestData(){
    mFttCollection->rawHits().clear();

    // TODO: inject clean strip hits to test cluster finder
    // StFttRawHit *hit = new StFttRawHit( sec, rdo, feb, vm, vmm[0].ch, vmm[0].adc, vmm[0].bcid, vmm[0].tb );
    // hit->setMapping( plane, quadrant, row, strip )
}

void StFttClusterPointMaker::MakeLocalPoints(UChar_t Rob) {
    LOG_INFO << "Making local points" << endm;

    //initialize
    StFttPoint* point;
    size_t nClusters_X = 0;size_t nClusters_Y = 0;size_t nClusters_DX = 0;size_t nClusters_DY = 0;

    //set nClusters
    nClusters_X = clustersPerRob[(UChar_t)Rob][kFttVertical].size();
    nClusters_Y = clustersPerRob[(UChar_t)Rob][kFttHorizontal].size();
    nClusters_DX = clustersPerRob[(UChar_t)Rob][kFttDiagonalV].size();
    nClusters_DY = clustersPerRob[(UChar_t)Rob][kFttDiagonalH].size();

    if(mDebug)
    {
        LOG_INFO << "rob = " << (int)Rob << endm;
        LOG_INFO << "nClusterX = " << nClusters_X << " nClusterY = " << nClusters_Y << endm;
    }

    //loop over x clusters; clustersPerRob[Rob][kFttVertical][iClu_X]
    for (int iClu_X=0; iClu_X < nClusters_X; iClu_X++) {
        StFttCluster* clu_x = clustersPerRob[(UChar_t)Rob][kFttVertical][iClu_X];

        point = new StFttPoint();

        std::vector<std::vector<float>> covMatrix(2,std::vector<float>(2,0));
        covMatrix[0][0] = clu_x->sigma()*clu_x->sigma();
        covMatrix[1][1] = clu_x->maxStripLength()*clu_x->maxStripLength()/12.;
        covMatrix[0][1] = 0;
        covMatrix[1][0] = 0;

        point->setX( clu_x->x() );
        point->setY( mFttDb->YX_StripGroupEdge[clu_x->row()]+clu_x->maxStripLength()/2. );
        point->setSigmaX(clu_x->sigma());
        point->setSigmaY(clu_x->maxStripLength()/sqrt(12));
        point->setSigmaXY(0);
        point->setCov(covMatrix);

        point->setQuadrant( clu_x->quadrant() );
        point->setPlane( clu_x->plane() );
        point->addCluster(clu_x,kFttVertical);

        mFttCollection->addPoint(point);
    }

    //loop over y clusters; clustersPerRob[Rob][kFttHorizontal][iClu_Y]
    for (int iClu_Y=0; iClu_Y < nClusters_Y; iClu_Y++) {
        StFttCluster* clu_y = clustersPerRob[(UChar_t)Rob][kFttHorizontal][iClu_Y];

        point = new StFttPoint();

        std::vector<std::vector<float>> covMatrix(2,std::vector<float>(2,0));
        covMatrix[0][0] = clu_y->maxStripLength()*clu_y->maxStripLength()/12.;
        covMatrix[1][1] = clu_y->sigma()*clu_y->sigma();
        covMatrix[0][1] = 0;
        covMatrix[1][0] = 0;

        point->setX(  mFttDb->YX_StripGroupEdge[clu_y->row()]+clu_y->maxStripLength()/2. );
        point->setY( clu_y->x() );
        point->setSigmaX(clu_y->maxStripLength()/sqrt(12.));
        point->setSigmaY(clu_y->sigma());
        point->setSigmaXY(0);
        point->setCov(covMatrix);

        point->setQuadrant( clu_y->quadrant() );
        point->setPlane( clu_y->plane() );
        point->addCluster(clu_y,kFttHorizontal);

        mFttCollection->addPoint(point);
    }

    //loop over dx clusters; clustersPerRob[Rob][kFttDiagonalV][iClu_DX]
    for (int iClu_DX=0; iClu_DX<nClusters_DX; iClu_DX++){
        StFttCluster* clu_dx = clustersPerRob[(UChar_t)Rob][kFttDiagonalV][iClu_DX];

        point = new StFttPoint();

        double x_prime = clu_dx->x();
        double y_prime = clu_dx->maxStripLength()/2.; //needs to be + or - depending on row 3 or row 4
        if (clu_dx->row() == 3) {y_prime = -y_prime;}
        if (clu_dx->row() == 4) {y_prime = y_prime;}

        double xvar_prime = (clu_dx->sigma())*(clu_dx->sigma());
        double yvar_prime = clu_dx->maxStripLength()*clu_dx->maxStripLength()/12.;

        double x = mFttDb->D_StripGroupEdge[0]+(sqrt(2)/2)*(x_prime-y_prime);
        double y = mFttDb->D_StripGroupEdge[0]+(sqrt(2)/2)*(x_prime+y_prime);

        double C_diag = (xvar_prime + yvar_prime)/2;
        double C_off = (xvar_prime - yvar_prime)/2;

        std::vector<std::vector<float>> covMatrix(2,std::vector<float>(2,0));
        covMatrix[0][0] = C_diag;
        covMatrix[1][1] = C_diag;
        covMatrix[0][1] = C_off;
        covMatrix[1][0] = C_off;

        point->setX( x );
        point->setY( y );
        point->setSigmaX(sqrt(C_diag));
        point->setSigmaY(sqrt(C_diag));
        point->setSigmaXY(C_off);
        point->setCov(covMatrix);

        point->setQuadrant( clu_dx->quadrant() );
        point->setPlane( clu_dx->plane() );
        point->addCluster(clu_dx,kFttDiagonalV);

        mFttCollection->addPoint(point);
    }

    //loop over dy clusters; clustersPerRob[Rob][kFttDiagonalH][iClu_DY]
    for (int iClu_DY=0; iClu_DY<nClusters_DY; iClu_DY++){
        StFttCluster* clu_dy = clustersPerRob[(UChar_t)Rob][kFttDiagonalH][iClu_DY];

        point = new StFttPoint();

        double x_prime = clu_dy->x();
        double y_prime = clu_dy->maxStripLength()/2.;
        if (clu_dy->row() == 3) {y_prime = y_prime;}
        if (clu_dy->row() == 4) {y_prime = -y_prime;}

        double xvar_prime = (clu_dy->sigma())*(clu_dy->sigma());
        double yvar_prime = clu_dy->maxStripLength()*clu_dy->maxStripLength()/12.;

        double x = mFttDb->D_StripGroupEdge[0]+(sqrt(2)/2)*(x_prime-y_prime);
        double y = mFttDb->D_StripGroupEdge[0]+(sqrt(2)/2)*(x_prime+y_prime);

        double C_diag = (xvar_prime + yvar_prime)/2;
        double C_off = (xvar_prime - yvar_prime)/2;

        std::vector<std::vector<float>> covMatrix(2,std::vector<float>(2,0));
        covMatrix[0][0] = C_diag;
        covMatrix[1][1] = C_diag;
        covMatrix[0][1] = C_off;
        covMatrix[1][0] = C_off;

        point->setX( x );
        point->setY( y );
        point->setSigmaX(sqrt(C_diag));
        point->setSigmaY(sqrt(C_diag));
        point->setSigmaXY(C_off);
        point->setCov(covMatrix);

        point->setQuadrant( clu_dy->quadrant() );
        point->setPlane( clu_dy->plane() );
        point->addCluster(clu_dy,kFttDiagonalV);

        mFttCollection->addPoint(point);
    }
}

void StFttClusterPointMaker::MakeGlobalPoints() {
    for ( StFttPoint * p : mFttCollection->points() ){
        if (mDebug && !p) {
            LOG_INFO << "Point is NULL" << endm;
        }
        float x=p->x(); float y=p->y(); float z=0; //local coordinates
        float dx = 0, dy = 0, dz = 0; //offset to global coordiantes
        float sx = 1, sy = 1, sz = 1; //reflection to global coordinates

        mFttDb->getGloablOffset_ClusterPoint( p->plane(), p->quadrant(), dx, sx, dy, sy, dz, sz );

        StThreeVectorD global;
        global.set( (x*sx)+dx, (y*sy)+dy, z+dz );

        if (p->quadrant() == 1 || p->quadrant() == 3) {
            p->setSigmaXY(-p->sigmaXY());
            std::vector<std::vector<float>> new_cov(2,std::vector<float>(2,0));
            new_cov[0][0] = p->cov()[0][0];
            new_cov[1][1] = p->cov()[1][1];
            new_cov[0][1] = -p->cov()[0][1];
            new_cov[1][0] = -p->cov()[1][0];
            p->setCov(new_cov);
        }

        if (mDebug) {
            LOG_INFO << "Global x: " << global.x() << " y: " << global.y() << " z: " << global.z() << endm;
        }

        p->setXYZ( global );

        if (mDebug) {
            LOG_INFO << "Point x: " << p->xyz().x() << " y: " << p->xyz().y() << " z: " << p->xyz().z() << endm;
            }
    }
}
