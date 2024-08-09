
#include "StFttSlowSimMaker.h"

#include "StEvent/StEvent.h"
#include "St_base/StMessMgr.h"

#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StFttCollection.h"
#include "StThreeVectorF.hh"

#include "StEvent/StFttRawHit.h"
#include "StFttDbMaker/StFttDb.h"

#include "TCanvas.h"
#include "TCernLib.h"
#include "TF1.h"
#include "TH2F.h"
#include "TLine.h"
#include "TString.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include <array>

#include "StarGenerator/UTIL/StarRandom.h"

namespace FttGlobal {
    const bool verbose = true;
}

StFttSlowSimMaker::StFttSlowSimMaker(const Char_t *name)
    : StMaker{name}
    // : StMaker{name}, mDebug(kFALSE), mFttDb( nullptr )
    {
        LOG_DEBUG << "StFttSlowSimMaker start"  << endm;
        LOG_INFO << "******** StFttSlowSimMaker::StFttSlowSimMaker = "<<name<<endm;
    }

int StFttSlowSimMaker::Init() {
    iEvent = 0;

    hXY = new TH2F( "hXY", ";X;Y;", 400, -100, 100, 400, -100, 100 );
    //initialize the cluster shape function, that need to fit to the real case
    fClusterProfile = new TF1("fClusterProfile", "gausn");// signal shape
    fClusterProfile->SetParameters(1.0, 0, 1.4 * 3.2);
    fClusterProfile->SetRange(-100, 100);
    fClusterProfile->SetNpx(1.e6);

    fClusterWidth = new TF1("fClusterWidth", "gausn");//width 
    fClusterWidth->SetParameters(1.0, 1.6, 1.4);
    fClusterWidth->SetRange(1, 10);
    fClusterWidth->SetNpx(1.e6);


    auto cwd = gDirectory;
    
    
    mTreeData.x.push_back(0.5f);
    mTreeData.y.push_back(0.5f);
    mTreeData.z.push_back(0.5f);

    mTreeFile = new TFile("StFttSlowSim.root", "RECREATE");
    mTreeFile->cd();
    mTree = new TTree("fttss", "ftt slow sim tree");
    // mTree->Branch("N",         &mTreeData. N, "N/I");
    mTree->Branch("xx",         &mTreeData. x  );
    mTree->Branch("xy",         &mTreeData. y  );
    mTree->Branch("xz",         &mTreeData. z  );

    mTree->Branch("px",         &mTreeData. px  );
    mTree->Branch("py",         &mTreeData. py  );
    mTree->Branch("pz",         &mTreeData. pz  );

    mTree->Branch("pt",         &mTreeData. pt  );
    mTree->Branch("eta",        &mTreeData. eta  );
    mTree->Branch("phi",        &mTreeData. phi  );

    mTree->Branch("de",         &mTreeData. de  );
    mTree->Branch("ds",         &mTreeData. ds  );
    mTree->Branch("tof",        &mTreeData. tof  );
    mTree->Branch("vid",        &mTreeData. vid  );
    mTree->Branch("trackp",     &mTreeData. trackp  );

    mTree->Branch("plane",        &mTreeData. plane  );
    mTree->Branch("quad",         &mTreeData. quad  );

    mTree->Branch("rdo",        &mTreeData. rdo);
    mTree->Branch("feb",        &mTreeData. feb);
    mTree->Branch("vmm",        &mTreeData. vmm);
    mTree->Branch("ch",         &mTreeData. ch);
    mTree->Branch("bcid",       &mTreeData. bcid);
    mTree->Branch("dbcud",      &mTreeData. dbcid);
    mTree->Branch("tb",         &mTreeData. tb);
    mTree->Branch("ADC",        &mTreeData. ADC);
    mTree->Branch("row",        &mTreeData. row);
    mTree->Branch("strip",      &mTreeData. strip);
    mTree->Branch("dir",        &mTreeData. dir);
    

    // mTree->SetAutoFlush(0);

    cwd->cd();

    mFile = new TFile( "fttQA.root", "RECREATE" );
    mFile->cd();
    mFile->Close();

    BookHistograms();

    return StMaker::Init();
}

Int_t StFttSlowSimMaker::Make() {
    if (FttGlobal::verbose) {
        LOG_INFO << "StFttSlowSimMaker::Make" << endm;
    }

    //initialize
    // Init();

    // Get the existing StEvent, or add one if it doesn't exist.
    StEvent *mEvent = static_cast<StEvent *>(GetDataSet("StEvent"));
    if (!mEvent) {
        mEvent = new StEvent;
        AddData(mEvent);
        LOG_DEBUG << "Creating StEvent" << endm;
    }

    mFttCollection=mEvent->fttCollection();
    if(!mFttCollection) {
        mFttCollection = new StFttCollection();
        mEvent->setFttCollection(mFttCollection);
        LOG_DEBUG <<"Added NEW StFttCollection"<<endm;
    } else {
        mFttCollection = mEvent->fttCollection();
        LOG_DEBUG <<"Found EXISTING StFttCollection"<<endm;
    }

    mFttDb = static_cast<StFttDb*>(GetDataSet("fttDb"));
    assert( mFttDb );

    // Read the g2t table
    St_g2t_fts_hit *hitTable = static_cast<St_g2t_fts_hit *>(GetDataSet("g2t_stg_hit"));
    if (!hitTable) {
        LOG_INFO << "g2t_stg_hit table is empty" << endm;
        return kStOk;
    } // if !hitTable

    mTreeData.clear();

    const int nhits = hitTable->GetNRows();
    const g2t_fts_hit_st *hit = hitTable->GetTable();
    // mTreeData.N = nhits;



    TString histoname;
    LOG_INFO << "g2t_stg_hit table has " << nhits << endm;
    for (int i = 0; i < nhits; i++) {
        LOG_INFO << "***************** now at hit ***************** " << i << endm;

        hit = (g2t_fts_hit_st *)hitTable->At(i);
        if (0 == hit)
            continue;
        //what's the unit of this? cm?  * 10 to scale from cm to mm
        float xhit = hit->x[0]*10;
        float yhit = hit->x[1]*10;
        float zhit = hit->x[2]*10;
        int volume_id = hit->volume_id;// what is the volume mean?
        hXY->Fill( xhit, yhit );
        LOG_INFO  << " xhit = " << xhit 
                  << " yhit = " << yhit 
                  << endm;


        mTreeData.x.push_back( xhit );
        mTreeData.y.push_back( yhit );
        mTreeData.z.push_back( zhit );

        TVector3 p( hit->p[0], hit->p[1], hit->p[2] );
        mTreeData.px.push_back( hit->p[0] );
        mTreeData.py.push_back( hit->p[1] );
        mTreeData.pz.push_back( hit->p[2] );

        mTreeData.pt.push_back( p.Pt() );
        mTreeData.eta.push_back( p.Eta() );
        mTreeData.phi.push_back( p.Phi() );

        //What's the ds and de mean? does them not relate to the fwd information
        mTreeData.ds.push_back( hit->ds );
        mTreeData.de.push_back( hit->de );
        mTreeData.tof.push_back( hit->tof );
        mTreeData.vid.push_back( hit->volume_id );

        mTreeData.trackp.push_back( hit->track_p );
        LOG_INFO << "after load the MC hits " << endm;


        // lets calculate the detector level info
        // it start from 0, is this correct? need to ask Daniel
        int plane = (volume_id - 1) / 4;
        int quad =  (volume_id - 1) % 4;
        int quad_global = GetQuad(xhit,yhit); // zhen add it using the suvery data to select the quad.
        LOG_INFO << "quad Dan = " << quad << " quad Zhen = " << quad_global << endm;
        LOG_INFO << "plane = " << plane << endm;
        if (quad_global < -0.0001) continue;

        mTreeData.plane.push_back( plane );
        mTreeData.quad.push_back( quad_global );

        int sec = plane + 1;
        int rdo = quad + 1;
        int rob = quad + ( plane * StFttDb::nQuadPerPlane ) + 1;
        LOG_INFO  << " sec = " << sec 
                  << " rdo = " << rdo
                  << " rob = " << rob 
                  << endm;

        // ----------- important check needed: --------------
        // 1. simulation data gives a local position or global position?
        // 2. now assume it give a global position
        // 3. where is the z? z is the center of the sTGC or something else?
        // ----------- important check needed: --------------

        // ----------- important note: --------------
        //for FttPointMaker, current version is using the mm to get all the things, the output below is mm
        // ----------- important note  --------------

        // ----------- get the center strip --------------
        // get the local position
        double xhit_local = -999;
        double yhit_local = -999;
        double zhit_local = -999;
        Global2Local_2D(xhit_local,yhit_local,xhit,yhit,plane,quad_global); // convert the global position to sTGC module local position

        LOG_INFO << "after get the center strip " << endm;
        // get the strip information for center strip of MC point
        int strip_x = -999; // vertical strip
        int strip_y = -999; // horizontal strip
        int strip_dv = -999; // diagnoal vertical
        int strip_dh = -999; // diagnoal vertical
        int row_x = -999; // vertical strip
        int row_y = -999; // horizontal strip
        int row_dv = -999; // diagnoal vertical
        int row_dh = -999; // diagnoal vertical
        LOG_INFO  << " xhit_local = " << xhit_local
                  << " yhit_local = " << yhit_local 
                  << endm;
        
        //make sure there is correct row and strip information
        bool is_vaild = GetStripandRow_XY( xhit_local, yhit_local, strip_x, strip_y, row_x, row_y);
        if ( !is_vaild ) continue;
        GetStripandRow_Diag( xhit_local, yhit_local, strip_dv, strip_dh, row_dv, row_dh);
        LOG_INFO << "after GetStripandRow_XY " << endm;
        LOG_INFO  << " row_x = " << row_x
                  << " strip_x = " << strip_x
                  << " row_y = " << row_y 
                  << " strip_y = " << strip_y 
                  << " row_dv = " << row_dv
                  << " strip_dv = " << strip_dv
                  << " row_dh = " << row_dh 
                  << " strip_dh = " << strip_dh 
                  << endm;


        //QA for the local X and local Y
        histoname = Form("LocalXY_Position_Plane_%d_Quad_%d",plane,quad_global);
        mH2d[ histoname.Data() ]->Fill(xhit_local,yhit_local);
        histoname = Form("LocalXY_Strip_Plane_%d_Quad_%d",plane,quad_global);
        mH2d[ histoname.Data() ]->Fill(strip_x,strip_y);
        // histoname = Form("hitmap_Plane%d_Quad%d_hit%d",sec-1,rdo-1,i);
        // TString histotitle = Form("hitmap_Plane%d_Quad%d_evt%d;Strip;Row",sec-1,rdo-1,i);
        // mH2d[ histoname.Data() ] = new TH2D(histoname.Data(),histotitle.Data(),200,0,200,6,0,6);
        // histoname = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d",plane,quad_global);
        // mH2d[ histoname.Data() ]->Fill(xhit,strip_x);
        // mH2d[ histoname.Data() ]->Fill(yhit,strip_y);
        // mH2d[ histoname.Data() ]->Fill(sqrt(xhit*xhit+yhit*yhit),strip_x);
        // mH2d[ histoname.Data() ]->Fill(sqrt(xhit*xhit+yhit*yhit),strip_y);
        histoname = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d_X",plane,quad_global);
        mH2d[ histoname.Data() ]->Fill(xhit,strip_x);
        histoname = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d_Y",plane,quad_global);
        mH2d[ histoname.Data() ]->Fill(yhit,strip_y);
        histoname = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d_Diag",plane,quad_global);
        mH2d[ histoname.Data() ]->Fill(((xhit+yhit)/sqrt(2)-23.6)/3.2,strip_dv);
        mH2d[ histoname.Data() ]->Fill(((xhit+yhit)/sqrt(2)-23.6)/3/2,strip_dh);


        //get the electronic information
        //get X, Y, diagnoalV, diagnoal H information separately 
        int feb = -999;
        int vmm = -999;
        int ch = -999;
        // get the X information and loop for the hits
        mFttDb->reverseHardwareMap(feb, vmm, ch, row_x, strip_x);//strip and row
        LOG_INFO << "after reverseHardwareMap X " << endm;
        SampleCluster(strip_x, xhit_local, row_x, sec, rdo, 0);
        // SampleCluster(strip_x, xhit_local, row_x, sec, rdo, 0, i);
        LOG_INFO << "after SampleCluster X " << endm;
        // get the Y information
        mFttDb->reverseHardwareMap(feb, vmm, ch, row_x, strip_y);
        SampleCluster(strip_y, yhit_local, row_y, sec, rdo, 0);
        // SampleCluster(strip_y, yhit_local, row_y, sec, rdo, 0, i);
        LOG_INFO << "after SampleCluster Y " << endm;
        // get the Dv information
        double local_diag = sqrt(xhit_local*xhit_local+yhit_local*yhit_local);
        if (row_dv >= 0 )
        {
            mFttDb->reverseHardwareMap(feb,vmm, ch, row_dv, strip_dv);
            SampleCluster(strip_dv, local_diag , row_dv, sec, rdo, 1);
            // SampleCluster(strip_dv, local_diag , row_dv, sec, rdo, 1, i);
            LOG_INFO << "after SampleCluster Dv " << endm;
        }
        // get the DH information
        if( row_dh >= 0 )
        {
            mFttDb->reverseHardwareMap(feb,vmm, ch, row_dv, strip_dh);
            SampleCluster(strip_dh, local_diag, row_dh, sec, rdo, 1);
            // SampleCluster(strip_dh, local_diag, row_dh, sec, rdo, 1, i);
            LOG_INFO << "after SampleCluster Dh " << endm;
        }
        LOG_INFO << "get the strip information for center strip of MC point " << endm;

        // below is the part now is done
        // figure out strip center X, Y, D1, D2
        // generate cluster in X, Y, D1, D2
        // function to generate cluster -> determines the # of strips
        // for each one, loop over hit strips 
        //  -> convert to electronic IDs
        //  -> write a raw hit


        // if ( feb % 2 != 0 ) { // odd
        // below is the part write the tree and save some histograms for the QA

    }

    LOG_INFO << "mTreeData.x.size() == " << mTreeData.x.size() << endm;

    mTree->Fill();
    iEvent++;

    return kStOk;
}

//     D(3) |  A(0)  
//   ------------------
//     C(2) |  B(1)  
int StFttSlowSimMaker::GetQuad(double x_global, double y_global)
{
    if(x_global > 8 && y_global > 94) return 0; //quadA
    if(x_global > 112 && y_global < 83) return 1; //quadA
    if(x_global < -107 && y_global < 82) return 2; //quadA
    if(x_global < -3 && y_global > 94) return 3; //quadA
    
    return -999;
}

void StFttSlowSimMaker::Global2Local_2D(double &x_local, double &y_local, double x_global, double y_global, int i_plane, int i_quad)
{
    // the input position is global position, do not need to add the quad information

    float dx = -999.;
    float dy = -999.;
    float dz = -999.;
    float sx = 0.;
    float sy = 0.;
    float sz = 0.;

    mFttDb->getGloablOffset(i_plane,i_quad,dx,sx,dy,sy,dz,sz);
    LOG_INFO << " sx = " << sx 
             << " dx = " << dx 
             << " sy = " << sy 
             << " dy = " << dy 
             << endm;


    //for different quad there already have the shift
    x_local = (x_global-dx)*sx;
    y_local = (y_global-dy)*sy;
    
}

//--------------- using the strip group information to get the row infomation
bool StFttSlowSimMaker::is_Group1(int &row_x, int &row_y, double x, double y)
{
    if( (14.60 <= x && x <= 172.29) && (14.60 <= y && y <= 172.29) ) { row_x = 0; row_y = 0; return kTRUE;}
    else return kFALSE;
}
bool StFttSlowSimMaker::is_Group2(int &row_x, int &row_y, double x, double y)
{
    if( (172.29 <= x && x <= 360.09) && (14.60 <= y && y <= 172.29) ) {row_x = 0; row_y = 1; return kTRUE;}
    else return kFALSE;
}
bool StFttSlowSimMaker::is_Group3(int &row_x, int &row_y, double x, double y)
{
    if( ( ( (360.09 <= x && x <= 504.2) && (14.60 <= y && y <= 172.29) ) || ( (504.2<= x && x <= 548.3) && (14.60 <= y && y <= 216.89) ) ) ) {row_x = 0; row_y = 2; return kTRUE;}
    else return kFALSE;
}
bool StFttSlowSimMaker::is_Group4(int &row_x, int &row_y, double x, double y)
{
    if( (14.60 <= x && x <= 172.29) && (172.29 <= y && y <= 360.09) ){row_x = 1; row_y = 0; return kTRUE;}
    else return kFALSE;
}
bool StFttSlowSimMaker::is_Group5(int &row_x, int &row_y, double x, double y)
{
    if( ( ((172.29 <= x && x <= 315.4) && (172.29 <= y && y <= 360.09)) || ((315.4 <= x && x <= 360.09) && (172.29 <= y && y <= 410.9)) || ((360.09 <= x && x <= 410.9) && (315.4 <= y && y <= 410.9)) ) ) {row_x = 1; row_y = 1; return kTRUE;}
    else return kFALSE;
}
bool StFttSlowSimMaker::is_Group6(int &row_x, int &row_y, double x, double y)
{
    if( (360.09 <= x && x <= 504.2) && (172.29 <= y && y <= 315.4) ) {row_x = 1; row_y = 2; return kTRUE;}
    else return kFALSE;
}
bool StFttSlowSimMaker::is_Group7(int &row_x, int &row_y, double x, double y)
{
    if( ( ( (360.09 <= y && y <= 504.2) && (14.60 <= x && x <= 172.29) ) || ( (504.2<= y && y <= 548.3) && (14.60 <= x && x <= 216.89) ) )  ) {row_x = 2; row_y = 0; return kTRUE;}
    else return kFALSE;
}
bool StFttSlowSimMaker::is_Group8(int &row_x, int &row_y, double x, double y)
{
    if( ((360.09 <= y && y <= 504.2) && (172.29 <= x && x <= 315.4)) ) { row_x = 2; row_y = 1; return kTRUE;}
    return kFALSE;
}

bool StFttSlowSimMaker::GetStripandRow_XY(float x_local, float y_local, int &strip_x, int &strip_y, int &row_x, int &row_y)
{
    row_x = -999;
    row_y = -999;
    strip_x = -999;
    strip_y = -999;

    strip_x = (x_local-StFttDb::X_StripGroupEdge[0])/3.2; //frist strip index is 0;
    strip_y = (y_local-StFttDb::X_StripGroupEdge[0])/3.2; //frist strip index is 0;

    if ( is_Group1( row_x, row_y, x_local, y_local) ) {return kTRUE; }
    if ( is_Group2( row_x, row_y, x_local, y_local) ) {return kTRUE; }
    if ( is_Group3( row_x, row_y, x_local, y_local) ) {return kTRUE; }
    if ( is_Group4( row_x, row_y, x_local, y_local) ) {return kTRUE; }
    if ( is_Group5( row_x, row_y, x_local, y_local) ) {return kTRUE; }
    if ( is_Group6( row_x, row_y, x_local, y_local) ) {return kTRUE; }
    if ( is_Group7( row_x, row_y, x_local, y_local) ) {return kTRUE; }
    if ( is_Group8( row_x, row_y, x_local, y_local) ) {return kTRUE; }

    if( mDebug )
    {
        LOG_INFO << "******** no X & Y strip & row information for current hit!!!!!!!!! ********"<<endm;
    }
    
    return kFALSE;
}
// if this hit do not have the dv or dh hit, the row and strip will be -999
void StFttSlowSimMaker::GetStripandRow_Diag(float x_local, float y_local, int &strip_dv, int &strip_dh, int &row_dv, int &row_dh)
{
    row_dv = -999;
    row_dh = -999;
    strip_dv = -999;
    strip_dh = -999;

    double intecpt = (x_local+y_local)/sqrt(2);
    if ( x_local > y_local )
    {
        strip_dv = (intecpt-StFttDb::DiagStripShift)/StFttDb::stripPitch; //frist strip index is 0;
        strip_dh = (intecpt-StFttDb::DiagStripShift)/StFttDb::stripPitch; //frist strip index is 0;

        row_dv = 3;
        if ( strip_dh < 58 )
        {
            row_dh = 4;
        }
    }
    else if ( x_local < y_local )
    {
        strip_dv = (intecpt-StFttDb::DiagStripShift)/StFttDb::stripPitch; //frist strip index is 0;
        strip_dh = (intecpt-StFttDb::DiagStripShift)/StFttDb::stripPitch; //frist strip index is 0;
        row_dh = 3;
        if ( strip_dv < 58 )
        {
            row_dv = 4;
        }
        
    }
    if( mDebug )
    {
        LOG_INFO << " row_dv = " << row_dv 
                 << " row_dg = " << row_dh
                 <<endm;
    }
    
    return;
}

// function used to sample the cluster
Int_t StFttSlowSimMaker::SampleCluster(int center_strip, double xhit_local, int row_x, int sec, int rdo, int is_diag) // 0 for XY and 1 for diag
{

    if(center_strip < 0) return 0;

    int nStrips = fClusterWidth->GetRandom();
    int half_width = nStrips/2;
    int is_odd = nStrips%2;
    double factor = -999; // used to scale the ADC

    double left_right_selection = gRandom->Uniform(-1,1); //first go to left side when it <0;
    int i_left = -999;
    int i_right = -999;
    if( is_odd == 1)//odd
    {
        i_left = center_strip-half_width;
        i_right = center_strip+half_width;
    } else if(is_odd == 0)
    {
        if (left_right_selection > 0) // left < right
        {
            i_left = center_strip-half_width+1;
            i_right = center_strip+half_width;
        } else if (left_right_selection < 0) //left > right
        {
            i_left = center_strip-half_width;
            i_right = center_strip+half_width-1;
        }
    }
    
    int ADC = gRandom->Uniform(100,1000);// sample the maximum ADC for this cluster
    int BCID = 1; // Zhen add this for test, now is constant 

    if (mDebug)
    {
        LOG_INFO << " the ininatial information for sample ADC of the MC FttRawHit is ::" 
             << " center strip = " << center_strip
             << " i_left = " << i_left
             << " i_right = " << i_right
             << " nStrips = " << nStrips
             << " half_width = " << half_width
             << " StFttDb::FirstStripEdge[is_diag] = " << StFttDb::FirstStripEdge[is_diag] 
             << " center_strip*StFttDb::stripPitch = " << center_strip*StFttDb::stripPitch
             << " xhit_local = " << xhit_local
             << " xhit_local-StFttDb::FirstStripEdge[is_diag] = " << xhit_local-StFttDb::FirstStripEdge[is_diag]
             << " center_strip*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]) = " << center_strip*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]) 
             << endm;
    }

    double integral_lowEdge = center_strip*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]);
    double integral_higEdge = (center_strip+1)*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]);
    factor = ADC/fClusterProfile->Integral(integral_lowEdge,integral_higEdge);

    for (int i = i_left; i <= i_right; i++)
    {
        int feb = -999;
        int vmm = -999;
        int ch = -999;

        integral_lowEdge = i*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]);
        integral_higEdge = (i+1)*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]);
        ADC = fClusterProfile->Integral(integral_lowEdge,integral_higEdge)*factor;
        mFttDb->reverseHardwareMap(feb, vmm, ch, row_x, i);//strip and row

        if (mDebug)
        {
            LOG_INFO << " the information for sample ADC of the MC FttRawHit is ::" 
                 << " integral_lowEdge = " << integral_lowEdge 
                 << " integral_higEdge = " << integral_higEdge 
                 << " ADC =  " << ADC
                 << " factor = " << factor
                 << " Integral = " << fClusterProfile->Integral(integral_lowEdge,integral_higEdge)
                 << endm;

            LOG_INFO << " the information for the MC FttRawHit is ::" 
                 << " feb = " << feb 
                 << " vmm = " << vmm 
                 << " ch =  " << ch
                 << " strip = " << i
                 << " row = " << row_x
                 << endm;
        }
        //for QA things
        TString histoname = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d",sec-1,rdo-1);
        mH2d[ histoname.Data() ]->Fill(xhit_local,i);

        //create a new raw hit
        StFttRawHit *hit = new StFttRawHit( sec, rdo, feb, vmm, ch, ADC, BCID, 1, 1 );// Zhen add it : now for the BCID, Tb, bcid_delta are set as a constant
        mFttCollection->addRawHit( hit ); //add it to the FttCollection
    }

    return 1;
}

// function used to sample the cluster. QA version
Int_t StFttSlowSimMaker::SampleCluster(int center_strip, double xhit_local, int row_x, int sec, int rdo, int is_diag, int i_evt) // 0 for XY and 1 for diag
{

    TString histoname = Form("hitmap_Plane%d_Quad%d_hit%d",sec-1,rdo-1,i_evt);

    if(center_strip < 0) return 0;

    int nStrips = fClusterWidth->GetRandom();
    int half_width = nStrips/2;
    int is_odd = nStrips%2;
    double factor = -999; // used to scale the ADC

    double left_right_selection = gRandom->Uniform(-1,1); //first go to left side when it <0;
    int i_left = -999;
    int i_right = -999;
    if( is_odd == 1)//odd
    {
        i_left = center_strip-half_width;
        i_right = center_strip+half_width;
    } else if(is_odd == 0)
    {
        if (left_right_selection > 0) // left < right
        {
            i_left = center_strip-half_width+1;
            i_right = center_strip+half_width;
        } else if (left_right_selection < 0) //left > right
        {
            i_left = center_strip-half_width;
            i_right = center_strip+half_width-1;
        }
    }
    
    int ADC = gRandom->Uniform(100,1000);// sample the maximum ADC for this cluster
    int BCID = 1; // Zhen add this for test, now is constant 

    LOG_INFO << " the ininatial information for sample ADC of the MC FttRawHit is ::" 
                 << " center strip = " << center_strip
                 << " i_left = " << i_left
                 << " i_right = " << i_right
                 << " nStrips = " << nStrips
                 << " half_width = " << half_width
                 << " StFttDb::FirstStripEdge[is_diag] = " << StFttDb::FirstStripEdge[is_diag] 
                 << " center_strip*StFttDb::stripPitch = " << center_strip*StFttDb::stripPitch
                 << " xhit_local = " << xhit_local
                 << " xhit_local-StFttDb::FirstStripEdge[is_diag] = " << xhit_local-StFttDb::FirstStripEdge[is_diag]
                 << " center_strip*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]) = " << center_strip*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]) 
                 << endm;
                 

    double integral_lowEdge = center_strip*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]);
    double integral_higEdge = (center_strip+1)*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]);
    factor = ADC/fClusterProfile->Integral(integral_lowEdge,integral_higEdge);

    for (int i = i_left; i <= i_right; i++)
    {
        int feb = -999;
        int vmm = -999;
        int ch = -999;

        integral_lowEdge = i*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]);
        integral_higEdge = (i+1)*StFttDb::stripPitch-(xhit_local-StFttDb::FirstStripEdge[is_diag]);
        ADC = fClusterProfile->Integral(integral_lowEdge,integral_higEdge)*factor;
        mFttDb->reverseHardwareMap(feb, vmm, ch, row_x, i);//strip and row

        if (mDebug)
        {
            LOG_INFO << " the information for sample ADC of the MC FttRawHit is ::" 
                 << " integral_lowEdge = " << integral_lowEdge 
                 << " integral_higEdge = " << integral_higEdge 
                 << " ADC =  " << ADC
                 << " factor = " << factor
                 << " Integral = " << fClusterProfile->Integral(integral_lowEdge,integral_higEdge)
                 << endm;

            LOG_INFO << " the information for the MC FttRawHit is ::" 
                 << " feb = " << feb 
                 << " vmm = " << vmm 
                 << " ch =  " << ch
                 << " strip = " << i
                 << " row = " << row_x
                 << endm;
        }
        mH2d[ histoname.Data() ]->Fill(i,row_x,ADC);

        //create a new raw hit
        StFttRawHit *hit = new StFttRawHit( sec, rdo, feb, vmm, ch, ADC, BCID, 1, 1 );// Zhen add it : now for the BCID, Tb, bcid_delta are set as a constant
        mFttCollection->addRawHit( hit ); //add it to the FttCollection
    }

    return 1;
}

//for QA 
void StFttSlowSimMaker::BookHistograms()
{

    TString name;
    TString title;
    //2D histograms for the QA 
    // for local position
    
    
    // TString label = StFttDb::orientationLabels[i];
    for (size_t i = 0; i < StFttDb::nPlane; i++)
    {
        for (size_t j = 0; j < StFttDb::nQuadPerPlane ; j++)
        {
            //TH2
            name = Form("LocalXY_Position_Plane_%d_Quad_%d",i,j);
            title =  Form("LocalXY_Position_Plane_%d_Quad_%d; x (mm), y (mm)",i,j);
            mH2d[ name.Data() ] = new TH2D(name.Data(), title.Data(), 3000, 0, 600, 3000, 0, 600);
            name = Form("LocalXY_Strip_Plane_%d_Quad_%d",i,j);
            title = Form("LocalXY_Strip_Plane_%d_Quad_%d; iStrip_{x}; iStrip_{Y}",i,j);
            mH2d[ name.Data() ] = new TH2D( name.Data(), title.Data(), 200, 0, 200, 200, 0, 200);
            name = Form("vmm_Plane%d_Quad%d",i,j);
            title = Form("vmm_Plane%d_Quad%d; vmm; Counts",i,j);
            mH2d[ "vmm" ] = new TH2D(name.Data(), title.Data(), 6,0,6,10,0,10);
            name = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d",i,j);
            title = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d; x (mm); strip ID (mm)",i,j);
            mH2d[ name.Data() ] = new TH2D(name.Data(), title.Data(), 1400,0,700,200,0,200);
            name = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d_X",i,j);
            title = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d_X; x (mm); strip ID (mm)",i,j);
            mH2d[ name.Data() ] = new TH2D(name.Data(), title.Data(), 1400,0,700,200,0,200);
            name = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d_Y",i,j);
            title = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d_Y; x (mm); strip ID (mm)",i,j);
            mH2d[ name.Data() ] = new TH2D(name.Data(), title.Data(), 1400,0,700,200,0,200);
            name = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d_Diag",i,j);
            title = Form("MCPosition_vs_StirpID_Plane_%d_Quad_%d_Diag; x (mm); strip ID (mm)",i,j);
            mH2d[ name.Data() ] = new TH2D(name.Data(), title.Data(), 1400,0,700,200,0,200);

            //TH1
            name = Form("feb_Plane%d_Quad%d",i,j);
            title = Form("feb_Plane%d_Quad%d; feb; Counts",i,j);
            mH1d[ "feb" ] = new TH1D(name.Data(), title.Data(), 10,0,10);
        }
    }
    
    return;
}
void StFttSlowSimMaker::WriteHistograms ()
{
    LOG_DEBUG << "StFttQAMaker::WriteHistograms()" << endm;
    for (const auto& kv : mH1d)
    {
        // if( kv.second->GetEntries() > 0 ) kv.second->Write();
        kv.second->Write();
    }
    for( const auto& kv : mH2d ) {
        // if( kv.second->GetEntries() > 0 ) kv.second->Write();
        kv.second->Write();
    }
}

Int_t StFttSlowSimMaker::Finish(){
    
    auto cwd = gDirectory;

    mTreeFile->cd();
    mTree->Write();
    mTreeFile->Write();
    mTreeFile->Close();

    LOG_DEBUG << "StFttQAMaker::WriteHistograms()" << endm;
    TFile * fOut = new TFile( "FttSlowSim.root", "RECREATE" );
    fOut->cd();

    hXY->Write();
    WriteHistograms();

    fOut->Write();
    fOut->Close();

    cwd->cd();

    return kStOk;

    
}
