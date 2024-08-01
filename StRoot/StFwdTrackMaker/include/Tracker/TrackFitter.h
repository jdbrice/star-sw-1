#ifndef TRACK_FITTER_H
#define TRACK_FITTER_H

#include "GenFit/ConstField.h"
#include "GenFit/EventDisplay.h"
#include "GenFit/Exception.h"
#include "GenFit/FieldManager.h"
#include "GenFit/KalmanFitStatus.h"
#include "GenFit/GblFitter.h"
#include "StFwdTrackMaker/StFwdGbl.h"

#include "GenFit/KalmanFitter.h"
#include "GenFit/KalmanFitterInfo.h"
#include "GenFit/KalmanFitterRefTrack.h"
#include "GenFit/MaterialEffects.h"
#include "GenFit/PlanarMeasurement.h"
#include "GenFit/RKTrackRep.h"
#include "GenFit/SpacepointMeasurement.h"
#include "GenFit/StateOnPlane.h"
#include "GenFit/TGeoMaterialInterface.h"
#include "GenFit/Track.h"
#include "GenFit/TrackPoint.h"

#include "Criteria/SimpleCircle.h"

#include "TDatabasePDG.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TNtuple.h"

#include <vector>
#include <memory>

#include "StFwdTrackMaker/Common.h"

#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StFwdTrackMaker/include/Tracker/TrackFitter.h"
#include "StFwdTrackMaker/include/Tracker/STARField.h"
#include "StFwdTrackMaker/include/Tracker/FwdGeomUtils.h"

#include "StarGenerator/UTIL/StarRandom.h"

#include "St_base/StMessMgr.h"
#include "St_db_Maker/St_db_Maker.h"

#include "tables/St_Survey_Table.h"
#include "StMaker.h"

/* Cass for fitting tracks(space points) with GenFit
 *
 */
class TrackFitter {

// Accessors and options
  public:
    std::shared_ptr<genfit::Track> getTrack() { return mFitTrack; }

  public:

    /**
     * @brief Construct a new Track Fitter object
     *
     * @param _mConfig : Config object
     * @param geoCache : Geometry cache filename
     */
    TrackFitter(FwdTrackerConfig _mConfig, TString geoCache) : mConfig(_mConfig), mGeoCache(geoCache), mFitTrack(nullptr) {}

    /**
     * @brief Setup the tracker object
     * Load geometry
     * Setup Material Effects
     * Setup the magnetic field
     * Setup the fitter
     * Setup the fit planes
     */
    void setup() {

        // Determine if we are using misaligned geometry to 
        mMisaligned = mConfig.get<bool>("Geometry:misaligned",false);       
        mMcTracking = mConfig.get<bool>("TrackFitter:McTracking",false);       
        
        mAlignmentOutput = new TFile("alignment.root", "RECREATE");
        mAlignmentOutput->cd();

        int BufSize = (int)pow(2., 16.);

        if(mMcTracking)
        {
          const char* varlist = "vid:mcx:mcy:mcz:rcx:rcy:rcz:trackx:tracky:mcxlab:mcylab:rcxlab:rcylab:trackxlab:trackylab"; // 15 variables
          mAlignmentInfo = new TNtuple("AlignmentInfo", "AlignmentInfo", varlist, BufSize);
        }
        else
        {
          const char* varlist = "vid:rcx:rcy:trackx:tracky:rcxlab:rcylab:trackxlab:trackylab"; // 9 variables, z positions will always be the same
          mAlignmentInfo = new TNtuple("AlignmentInfo", "AlignmentInfo", varlist, BufSize);
        }

        // the geometry manager that GenFit will use
        TGeoManager * gMan = nullptr;

        // Setup the Geometry used by GENFIT
        LOG_INFO << "StFwdTrackMaker is loading the geometry cache: " << mConfig.get<string>("Geometry", mGeoCache.Data()).c_str() << endm;
        TGeoManager::Import(mConfig.get<string>("Geometry", mGeoCache.Data()).c_str());
        gMan = gGeoManager;
        // Set up the material interface and set material effects on/off from the config
        genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

        // Set Material Stepper debug level
        genfit::MaterialEffects::getInstance()->setDebugLvl( mConfig.get<int>("TrackFitter.MaterialEffects:DebugLvl", 0) );

        genfit::MaterialEffects::getInstance()->setEnergyLossBetheBloch( mConfig.get<int>("TrackFitter.MaterialEffects.EnergyLossBetheBloch", true) );
        genfit::MaterialEffects::getInstance()->setNoiseBetheBloch( mConfig.get<int>("TrackFitter.MaterialEffects.NoiseBetheBloch", true) );
        genfit::MaterialEffects::getInstance()->setNoiseCoulomb( mConfig.get<int>("TrackFitter.MaterialEffects.NoiseCoulomb", true) );
        genfit::MaterialEffects::getInstance()->setEnergyLossBrems( mConfig.get<int>("TrackFitter.MaterialEffects.EnergyLossBrems", true) );
        genfit::MaterialEffects::getInstance()->setNoiseBrems( mConfig.get<int>("TrackFitter.MaterialEffects.NoiseBrems", true) );
        genfit::MaterialEffects::getInstance()->ignoreBoundariesBetweenEqualMaterials( mConfig.get<int>("TrackFitter.MaterialEffects.ignoreBoundariesBetweenEqualMaterials", true) );

        // do this last to override
        genfit::MaterialEffects::getInstance()->setNoEffects( !mConfig.get<bool>("TrackFitter:MaterialEffects", false)); // negated, true means defaul is all effects on (noEffects off)
        if (!mConfig.get<bool>("TrackFitter:MaterialEffects", false)){
            LOG_INFO << "Turning OFF GenFit Material Effects in stepper" << endm;
        }

        // Determine which Magnetic field to use
        // Either constant field or real field from StarFieldAdaptor
        if (mConfig.get<bool>("TrackFitter:constB", false)) {
            mBField = std::unique_ptr<genfit::AbsBField>(new genfit::ConstField(0., 0., 5.)); // 0.5 T Bz
            LOG_INFO << "StFwdTrackMaker: Tracking with constant magnetic field" << endm;
        } else if (mConfig.get<bool>("TrackFitter:zeroB", false)) {
            mBField = std::unique_ptr<genfit::AbsBField>(new genfit::ConstField(0., 0., 0.)); // ZERO FIELD
            LOG_INFO << "StFwdTrackMaker: Tracking with ZERO magnetic field" << endm;
        } else {
            mBField = std::unique_ptr<genfit::AbsBField>(new StarFieldAdaptor());
            LOG_INFO << "StFwdTrackMaker: Tracking with StarFieldAdapter" << endm;
        }
        // we must have one of the two available fields at this point
        // note, the pointer is still bound to the lifetime of the TackFitter
        genfit::FieldManager::getInstance()->init(mBField.get());

        // initialize the main mFitter using a KalmanFitter with reference tracks
        mFitter = std::unique_ptr<genfit::AbsKalmanFitter>(new genfit::KalmanFitterRefTrack());

        // Here we load several options from the config,
        // to customize the mFitter behavior
        mFitter->setMaxFailedHits(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MaxFailedHits", -1)); // default -1, no limit
        mFitter->setDebugLvl(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:DebugLvl", 0)); // default 0, no output
        mFitter->setMaxIterations(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MaxIterations", 40)); // default 4 iterations
        mFitter->setMinIterations(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MinIterations", 3)); // default 0 iterations

        if(mConfig.get<bool>("TrackFitter:refitGBL",false)) {
            mRefitGBL = true;
            mGblFitter = std::unique_ptr<StFwdGbl>(new StFwdGbl());
            mGblFitter->setGBLOptions("",true,false);
            mGblFitter->setMP2Options(0.0,1,"millefile.dat",0.);
            mGblFitter->setDebugLvl(1);
            mGblFitter->setAllFSTFlag(mConfig.get<bool>("TrackFitter:allFST", false));
            mGblFitter->setIncludeVertexFlag(mIncludeVertexInFit);
            mGblFitter->beginRun();
        }

        // FwdGeomUtils looks into the loaded geometry and gets detector z locations if present
        FwdGeomUtils fwdGeoUtils( gMan );

        // these default values are the default if the detector is
        // a) not found in the geometry
        // b) not provided in config

        // NOTE: these defaults are needed since the geometry file might not include FST (bug being worked on separately)
        mFSTZLocations = fwdGeoUtils.fstZ(
            mConfig.getVector<double>("TrackFitter.Geometry:fst", 
                {151.750, 165.248, 178.781 }
                // 144.633,158.204,171.271
            )
        );

        if ( fwdGeoUtils.fstZ( 0 ) < 1.0 ) { // returns 0.0 on failure
            LOG_WARN << "Using FST z-locations from config or defautl, may not match hits" << endm;
        }

        //-----  LOAD ALIGNMENT MATRICES  -----//
        St_db_Maker *dbMk=new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");

        dbMk->SetDebug();
        dbMk->SetDateTime(20211026,7); // event or run start time, set to your liking
        dbMk->SetFlavor("ofl");

        dbMk->Init();
        dbMk->Make();
        
        // Load in FST alignment Tables

        // FST Survey Tables
        Survey_st *fstOnTpc;
        Survey_st *hssOnFst;
        Survey_st *fstWedgeOnHss;
        Survey_st *fstSensorOnWedge;

        // FST Matrices for Transformations
        TMatrixD MfstOnTpc(4,4);
        TMatrixD MhssOnFst(4,4);
        TMatrixD MfstWedgeOnHss(4,4);
        TMatrixD MfstSensorOnWedge(4,4);
        MfstOnTpc.Zero();
        MhssOnFst.Zero();
        MfstWedgeOnHss.Zero();
        MfstSensorOnWedge.Zero();


        TDataSet *DB = 0;
        DB = dbMk->GetDataBase("Geometry/fst/fstOnTpc");
        if (!DB) {
          std::cout << "ERROR: no fstOnTpc table found in db, or malformed local db config" << std::endl;
        }

        St_Survey *dataset = 0;
        dataset = (St_Survey*) DB->Find("fstOnTpc");

        if (dataset) {

            Int_t rows = dataset->GetNRows();
            if (rows > 1) {
              std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
            }

            TDatime val[2];
            dbMk->GetValidity((TTable*)dataset,val);
            std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - "
              << val[1].GetDate() << "." << val[1].GetTime() << " ] "
              << std::endl;

            fstOnTpc = dataset->GetTable();
            
            MfstOnTpc(0,0) = fstOnTpc[0].r00;
            MfstOnTpc(0,1) = fstOnTpc[0].r01;
            MfstOnTpc(0,2) = fstOnTpc[0].r02;
            MfstOnTpc(1,0) = fstOnTpc[0].r10;
            MfstOnTpc(1,1) = fstOnTpc[0].r11;
            MfstOnTpc(1,2) = fstOnTpc[0].r12;
            MfstOnTpc(2,0) = fstOnTpc[0].r20;
            MfstOnTpc(2,1) = fstOnTpc[0].r21;
            MfstOnTpc(2,2) = fstOnTpc[0].r22;
            MfstOnTpc(0,3) = fstOnTpc[0].t0 ;
            MfstOnTpc(1,3) = fstOnTpc[0].t1 ;
            MfstOnTpc(2,3) = fstOnTpc[0].t2 ;
            MfstOnTpc(3,3) = 1.0            ;

        } else {
            std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        }
      
        DB = dbMk->GetDataBase("Geometry/fst/hssOnFst");
        if (!DB) {
          std::cout << "ERROR: no hssOnFst table found in db, or malformed local db config" << std::endl;
        }

        dataset = (St_Survey*) DB->Find("hssOnFst");

        if (dataset) {

            Int_t rows = dataset->GetNRows();
            if (rows > 1) {
              std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
            }

            TDatime val[2];
            dbMk->GetValidity((TTable*)dataset,val);
            std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - "
              << val[1].GetDate() << "." << val[1].GetTime() << " ] "
              << std::endl;

            hssOnFst = dataset->GetTable();

        } else {
            std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        }

        DB = dbMk->GetDataBase("Geometry/fst/fstWedgeOnHss");
        if (!DB) {
          std::cout << "ERROR: no fstWedgeOnHss table found in db, or malformed local db config" << std::endl;
        }

        dataset = (St_Survey*) DB->Find("fstWedgeOnHss");

        if (dataset) {

            Int_t rows = dataset->GetNRows();
            if (rows > 1) {
              std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
            }

            TDatime val[2];
            dbMk->GetValidity((TTable*)dataset,val);
            std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - "
              << val[1].GetDate() << "." << val[1].GetTime() << " ] "
              << std::endl;

            fstWedgeOnHss = dataset->GetTable();

        } else {
            std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        }

        DB = dbMk->GetDataBase("Geometry/fst/fstSensorOnWedge");
        if (!DB) {
          std::cout << "ERROR: no fstSensorOnWedge table found in db, or malformed local db config" << std::endl;
        }

        dataset = (St_Survey*) DB->Find("fstSensorOnWedge");

        if (dataset) {

            Int_t rows = dataset->GetNRows();
            if (rows > 1) {
              std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
            }

            TDatime val[2];
            dbMk->GetValidity((TTable*)dataset,val);
            std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - "
              << val[1].GetDate() << "." << val[1].GetTime() << " ] "
              << std::endl;

            fstSensorOnWedge = dataset->GetTable();

        } else {
            std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        }


        // Now add the Si detector planes at the desired location
        std::stringstream sstr;
        sstr << "Adding FST Planes at: ";
        string delim = "";
        for (auto z : mFSTZLocations) {
            mFSTPlanes.push_back(
                genfit::SharedPlanePtr(
                    // these normals make the planes face along z-axis
                    new genfit::DetPlane(TVector3(0, 0, z), TVector3(1, 0, 0), TVector3(0, 1, 0) )
                )
            );
        }
        LOG_DEBUG  << sstr.str() << endm;

        // Create GenFit planes for FST 
 
        // default FST sensor z locations, found externally from ideal geometry
        for (int is = 0; is < 108; is++) {

            double z;

            // FST half
            int h = (is / 18) % 2; // 0 (left +x half), 1 (right -x half) 

            // FST disk (integer division rounds down)
            int d = is / 36; // 0-2

            // FST wedge
            int w = is / 3; // 0-35

            // FST sensor
            int s = is % 3; // 0 (inner), 1 (outer), 2 (outer)
            int ds = (s == 0)? 0 : 1; // +0 for inner, +1 for outer

            int defaultZidx = d * 4 + 2 * (w % 2) + ds;
            z = fstDefaultZ[defaultZidx];


            double angle = TMath::Pi()*5./12. - double(w)*TMath::Pi()*1./6.;

            double phiAxisShift = 0.0;
            if(d == 0 || d == 2)
            {
              if(w%2 == 0)
              {
                if(s == 1)      phiAxisShift = +8.0 * TMath::Pi() / 180.0;
                else if(s == 2) phiAxisShift = -8.0 * TMath::Pi() / 180.0;
              }
              else if(w%2 == 1)
              {
                if(s == 1)      phiAxisShift = -8.0 * TMath::Pi() / 180.0;
                else if(s == 2) phiAxisShift = +8.0 * TMath::Pi() / 180.0;
              }
            }
            else if(d == 1)
            {
              if(w%2 == 0)
              {
                if(s == 1)      phiAxisShift = -8.0 * TMath::Pi() / 180.0;
                else if(s == 2) phiAxisShift = +8.0 * TMath::Pi() / 180.0;
              }
              else if(w%2 == 1)
              {
                if(s == 1)      phiAxisShift = +8.0 * TMath::Pi() / 180.0;
                else if(s == 2) phiAxisShift = -8.0 * TMath::Pi() / 180.0;
              }
            }

            angle += phiAxisShift;

            TMatrixD RotZ(4,4);
            RotZ(0,0) =  TMath::Cos(angle);
            RotZ(0,1) = -TMath::Sin(angle);
            RotZ(1,0) =  TMath::Sin(angle);
            RotZ(1,1) =  TMath::Cos(angle);
            RotZ(2,2) =  1.0              ;
            RotZ(3,3) =  1.0              ;

            double ua[4] = {1.,0.,0.,1.};
            double va[4] = {0.,1.,0.,1.};
            double oa[4] = {0.,0.,0.,1.};
            if(s == 0) oa[0] = 10.75;
            if(s == 1 || s == 2) oa[0] = 22.25;
            TMatrixD u4(4,1,ua); // default u (corresponds to x on the sensor, (x,0,0) = (r,0,0) along x-axis) 
            TMatrixD v4(4,1,va); // default v (corresponds to y on the sensor)
            TMatrixD o4(4,1,oa); // default origin

            // create matrices from alignment tables
            MhssOnFst(0,0) = 1;//hssOnFst[h].r00;
            MhssOnFst(0,1) = 0;//hssOnFst[h].r01;
            MhssOnFst(0,2) = 0;//hssOnFst[h].r02;
            MhssOnFst(1,0) = 0;//hssOnFst[h].r10;
            MhssOnFst(1,1) = 1;//hssOnFst[h].r11;
            MhssOnFst(1,2) = 0;//hssOnFst[h].r12;
            MhssOnFst(2,0) = 0;//hssOnFst[h].r20;
            MhssOnFst(2,1) = 0;//hssOnFst[h].r21;
            MhssOnFst(2,2) = 1;//hssOnFst[h].r22;
            MhssOnFst(0,3) = 0;//hssOnFst[h].t0 ;
            MhssOnFst(1,3) = 0;//hssOnFst[h].t1 ;
            MhssOnFst(2,3) = 0;//hssOnFst[h].t2 ;
            MhssOnFst(3,3) = 1.0            ;

            MfstWedgeOnHss(0,0) = fstWedgeOnHss[w].r00;
            MfstWedgeOnHss(0,1) = fstWedgeOnHss[w].r01;
            MfstWedgeOnHss(0,2) = fstWedgeOnHss[w].r02;
            MfstWedgeOnHss(1,0) = fstWedgeOnHss[w].r10;
            MfstWedgeOnHss(1,1) = fstWedgeOnHss[w].r11;
            MfstWedgeOnHss(1,2) = fstWedgeOnHss[w].r12;
            MfstWedgeOnHss(2,0) = fstWedgeOnHss[w].r20;
            MfstWedgeOnHss(2,1) = fstWedgeOnHss[w].r21;
            MfstWedgeOnHss(2,2) = fstWedgeOnHss[w].r22;
            MfstWedgeOnHss(0,3) = fstWedgeOnHss[w].t0 ;
            MfstWedgeOnHss(1,3) = fstWedgeOnHss[w].t1 ;
            MfstWedgeOnHss(2,3) = fstWedgeOnHss[w].t2 ;
            MfstWedgeOnHss(3,3) = 1.0                 ;

            MfstSensorOnWedge(0,0) = fstSensorOnWedge[is].r00;
            MfstSensorOnWedge(0,1) = fstSensorOnWedge[is].r01;
            MfstSensorOnWedge(0,2) = fstSensorOnWedge[is].r02;
            MfstSensorOnWedge(1,0) = fstSensorOnWedge[is].r10;
            MfstSensorOnWedge(1,1) = fstSensorOnWedge[is].r11;
            MfstSensorOnWedge(1,2) = fstSensorOnWedge[is].r12;
            MfstSensorOnWedge(2,0) = fstSensorOnWedge[is].r20;
            MfstSensorOnWedge(2,1) = fstSensorOnWedge[is].r21;
            MfstSensorOnWedge(2,2) = fstSensorOnWedge[is].r22;
            MfstSensorOnWedge(0,3) = fstSensorOnWedge[is].t0 ;
            MfstSensorOnWedge(1,3) = fstSensorOnWedge[is].t1 ;
            MfstSensorOnWedge(2,3) = fstSensorOnWedge[is].t2 ;
            MfstSensorOnWedge(3,3) = 1.0                     ;

            // Rotate and Translate plane normal vectors and origin
            TMatrixD M = MfstOnTpc * MhssOnFst * MfstWedgeOnHss * MfstSensorOnWedge * RotZ;
            u4 = M * u4;
            v4 = M * v4;
            o4 = M * o4;

            cout << "Sensor Matrix " << endl;
            M.Print();

            // save this inverse matrix for use with misaligned simulated data
            //TMatrixD MnoRot = MfstOnTpc * MhssOnFst * MfstWedgeOnHss * MfstSensorOnWedge;
            mInverseMFst[is].ResizeTo(4,4);
            mInverseMFst[is] = M.Invert();

            TVector3 u(u4(0,0),u4(1,0),u4(2,0));
            TVector3 v(v4(0,0),v4(1,0),v4(2,0));
            TVector3 o(o4(0,0),o4(1,0),o4(2,0));

            mFSTSensorPlanes.push_back(
                genfit::SharedPlanePtr(
                    // these normals make the planes face along z-axis
                    new genfit::DetPlane(TVector3(0,0,z)+o, u, v)
                )
            );
            mFSTSensorPlanes[is]->Print();

            MhssOnFst.Zero();
            MfstWedgeOnHss.Zero();
            MfstSensorOnWedge.Zero();
        }

        // Load in STGC alignment Tables
       
        // STGC Survey Tables
        Survey_st *stgcOnTpc;
        Survey_st *stationOnStgc;
        Survey_st *pentOnStation;

        // STGC Matrices for Transformations
        TMatrixD MstgcOnTpc(4,4);
        TMatrixD MstationOnStgc(4,4);
        TMatrixD MpentOnStation(4,4);
        MstgcOnTpc.Zero();
        MstationOnStgc.Zero();
        MpentOnStation.Zero();


        //DB = dbMk->GetDataBase("Geometry/stgc/stgcOnTpc");
        //if (!DB) {
        //  std::cout << "ERROR: no stgcOnTpc table found in db, or malformed local db config" << std::endl;
        //}

        //dataset = (St_Survey*) DB->Find("stgcOnTpc");

        //if (dataset) {

        //    Int_t rows = dataset->GetNRows();
        //    if (rows > 1) {
        //      std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
        //    }

        //    TDatime val[2];
        //    dbMk->GetValidity((TTable*)dataset,val);
        //    std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - "
        //      << val[1].GetDate() << "." << val[1].GetTime() << " ] "
        //      << std::endl;

        //    stgcOnTpc = dataset->GetTable();

        //    MstgcOnTpc(0,0) = stgcOnTpc[0].r00;
        //    MstgcOnTpc(0,1) = stgcOnTpc[0].r01;
        //    MstgcOnTpc(0,2) = stgcOnTpc[0].r02;
        //    MstgcOnTpc(1,0) = stgcOnTpc[0].r10;
        //    MstgcOnTpc(1,1) = stgcOnTpc[0].r11;
        //    MstgcOnTpc(1,2) = stgcOnTpc[0].r12;
        //    MstgcOnTpc(2,0) = stgcOnTpc[0].r20;
        //    MstgcOnTpc(2,1) = stgcOnTpc[0].r21;
        //    MstgcOnTpc(2,2) = stgcOnTpc[0].r22;
        //    MstgcOnTpc(0,3) = stgcOnTpc[0].t0 ;
        //    MstgcOnTpc(1,3) = stgcOnTpc[0].t1 ;
        //    MstgcOnTpc(2,3) = stgcOnTpc[0].t2 ;
        //    MstgcOnTpc(3,3) = 1.0            ;

        //} else {
        //    std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        //}

        //DB = dbMk->GetDataBase("Geometry/stgc/stationOnStgc");
        //if (!DB) {
        //  std::cout << "ERROR: no stationOnStgc table found in db, or malformed local db config" << std::endl;
        //}

        //dataset = (St_Survey*) DB->Find("stationOnStgc");

        //if (dataset) {

        //    Int_t rows = dataset->GetNRows();
        //    if (rows > 1) {
        //      std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
        //    }

        //    TDatime val[2];
        //    dbMk->GetValidity((TTable*)dataset,val);
        //    std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - "
        //      << val[1].GetDate() << "." << val[1].GetTime() << " ] "
        //      << std::endl;

        //    stationOnStgc = dataset->GetTable();

        //} else {
        //    std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        //}

        //DB = dbMk->GetDataBase("Geometry/stgc/pentOnStation");
        //if (!DB) {
        //  std::cout << "ERROR: no pentOnStation table found in db, or malformed local db config" << std::endl;
        //}

        //dataset = (St_Survey*) DB->Find("pentOnStation");

        //if (dataset) {

        //    Int_t rows = dataset->GetNRows();
        //    if (rows > 1) {
        //      std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
        //    }

        //    TDatime val[2];
        //    dbMk->GetValidity((TTable*)dataset,val);
        //    std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - "
        //      << val[1].GetDate() << "." << val[1].GetTime() << " ] "
        //      << std::endl;

        //    pentOnStation = dataset->GetTable();

        //} else {
        //    std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        //}

        // Create GenFit planes for STGC  
      
        // default STGC z locations (Optimal for MC where hit is registered on the first gas plane hit)

        for (int ip = 0; ip < 16; ip++) { // 16 pentagons

            // STGC plane
            int p = ip/4;

            // STGC quadrant
            int q = ip%4; 

            double z = stgcDefaultZ[p];

            // Currently, the rotation is not needed for STGC coordinate system
            //TMatrixD RotZ(4,4);
            //RotZ(0,0) =  TMath::Cos(angle);
            //RotZ(0,1) = -TMath::Sin(angle);
            //RotZ(1,0) =  TMath::Sin(angle);
            //RotZ(1,1) =  TMath::Cos(angle);
            //RotZ(2,2) =  1.0              ;
            //RotZ(3,3) =  1.0              ;

            double ua[4] = {1.,0.,0.,1.};
            double va[4] = {0.,1.,0.,1.};
            double oa[4] = {0.,0.,0.,1.};
 
            // Shift STGC origins according to quadrant 
            switch(q) {
              case 0 : // +x,+y
                oa[0] = mCenterPent;
                oa[1] = 5.9+mCenterPent;
                break;
              case 1 : // -x,+y
                oa[0] = -mCenterPent;
                oa[1] = 5.9+mCenterPent;
                break;
              case 2 : // -x,-y
                oa[0] = -6.5-mCenterPent;
                oa[1] = 5.9-mCenterPent;
                break;
              case 3 : // +x,-y
                oa[0] = 6.5+mCenterPent;
                oa[1] = 5.9-mCenterPent;
                break;
            }

            TMatrixD u4(4,1,ua); // default u (corresponds to x on the pentagon) 
            TMatrixD v4(4,1,va); // default v (corresponds to y on the pentagon)
            TMatrixD o4(4,1,oa); // default origin

            // create matrices from alignment tables
            //MstationOnStgc(0,0) = stationOnStgc[p].r00;
            //MstationOnStgc(0,1) = stationOnStgc[p].r01;
            //MstationOnStgc(0,2) = stationOnStgc[p].r02;
            //MstationOnStgc(1,0) = stationOnStgc[p].r10;
            //MstationOnStgc(1,1) = stationOnStgc[p].r11;
            //MstationOnStgc(1,2) = stationOnStgc[p].r12;
            //MstationOnStgc(2,0) = stationOnStgc[p].r20;
            //MstationOnStgc(2,1) = stationOnStgc[p].r21;
            //MstationOnStgc(2,2) = stationOnStgc[p].r22;
            //MstationOnStgc(0,3) = stationOnStgc[p].t0 ;
            //MstationOnStgc(1,3) = stationOnStgc[p].t1 ;
            //MstationOnStgc(2,3) = stationOnStgc[p].t2 ;
            //MstationOnStgc(3,3) = 1.0                 ;

            //MpentOnStation(0,0) = pentOnStation[ip].r00;
            //MpentOnStation(0,1) = pentOnStation[ip].r01;
            //MpentOnStation(0,2) = pentOnStation[ip].r02;
            //MpentOnStation(1,0) = pentOnStation[ip].r10;
            //MpentOnStation(1,1) = pentOnStation[ip].r11;
            //MpentOnStation(1,2) = pentOnStation[ip].r12;
            //MpentOnStation(2,0) = pentOnStation[ip].r20;
            //MpentOnStation(2,1) = pentOnStation[ip].r21;
            //MpentOnStation(2,2) = pentOnStation[ip].r22;
            //MpentOnStation(0,3) = pentOnStation[ip].t0 ;
            //MpentOnStation(1,3) = pentOnStation[ip].t1 ;
            //MpentOnStation(2,3) = pentOnStation[ip].t2 ;
            //MpentOnStation(3,3) = 1.0                  ;

            //// Rotate and Translate plane normal vectors and origin
            //TMatrixD M = MstgcOnTpc * MstationOnStgc * MpentOnStation;// * RotZ;
            TMatrixD M(4,4);
            M.UnitMatrix();

            u4 = M * u4;
            v4 = M * v4;
            o4 = M * o4;

            cout << "STGC Pentagon Matrix " << ip << endl;
            M.Print();

            // save this inverse matrix for use with misaligned simulated data
            mInverseMStgc[ip].ResizeTo(4,4);
            mInverseMStgc[ip] = M.Invert();

            TVector3 u(u4(0,0),u4(1,0),u4(2,0));
            TVector3 v(v4(0,0),v4(1,0),v4(2,0));
            TVector3 o(o4(0,0),o4(1,0),o4(2,0));

            mFTTPlanes.push_back(
                genfit::SharedPlanePtr(
                    // these normals make the planes face along z-axis
                    new genfit::DetPlane(TVector3(0,0,z)+o, u, v)
                )
            );
            mFTTPlanes[ip]->Print();

            MstationOnStgc.Zero();
            MpentOnStation.Zero();
        }

        //-----  LOAD ALIGNMENT MATRICES  -----//

        // get default vertex values used in simulation from the config
        mVertexSigmaXY = mConfig.get<double>("TrackFitter.Vertex:sigmaXY", 1.0);
        mVertexSigmaZ = mConfig.get<double>("TrackFitter.Vertex:sigmaZ", 30.0);
        mIncludeVertexInFit = mConfig.get<bool>("TrackFitter.Vertex:includeInFit", false);
        mSmearMcVertex = mConfig.get<bool>("TrackFitter.Vertex:smearMcVertex", false);

        mRSize = mConfig.get<float>("SiRasterizer:r",3.0);
        mPhiSize = mConfig.get<float>("SiRasterizer:phi",0.004);

        if ( mGenHistograms )
            makeHistograms();
    }

    TMatrixDSym makeSiCovMat(TVector3 hit) {
        // we can calculate the CovMat since we know the det info, but in future we should probably keep this info in the hit itself
    
        float rSize = mRSize; //xfg.get<float>("SiRasterizer:r", 3.0);
        float phiSize = mPhiSize; //xfg.get<float>("SiRasterizer:phi", 0.004);
    
        // measurements on a plane only need 2x2
        // for Si geom we need to convert from cylindrical to cartesian coords
        TMatrixDSym cm(2);
        TMatrixD T(2, 2);
        TMatrixD J(2, 2);
        const float x = hit.X();
        const float y = hit.Y();
        const float R = sqrt(x * x + y * y);
        const float cosphi = x / R;
        const float sinphi = y / R;
        const float sqrt12 = sqrt(12.);
    
        const float dr = rSize / sqrt12;
        const float dphi = (phiSize) / sqrt12;
    
        // Setup the Transposed and normal Jacobian transform matrix;
        // note, the si fast sim did this wrong
        // row col
        T(0, 0) = cosphi;
        T(0, 1) = -R * sinphi;
        T(1, 0) = sinphi;
        T(1, 1) = R * cosphi;
    
        J(0, 0) = cosphi;
        J(0, 1) = sinphi;
        J(1, 0) = -R * sinphi;
        J(1, 1) = R * cosphi;
    
        TMatrixD cmcyl(2, 2);
        cmcyl(0, 0) = dr * dr;
        cmcyl(1, 1) = dphi * dphi;
    
        TMatrixD r = T * cmcyl * J;
    
        // note: float sigmaX = sqrt(r(0, 0));
        // note: float sigmaY = sqrt(r(1, 1));
    
        cm(0, 0) = r(0, 0);
        cm(1, 1) = r(1, 1);
        cm(0, 1) = r(0, 1);
        cm(1, 0) = r(1, 0);
    
        return cm;
    }
   
    /**
     * @brief end the Gbl Fitter
     * 
     */
    void finish() {
        mGblFitter->endRun();
    }

    /**
     * @brief Writing alignment histograms
     * 
     */
    void writeAlignment() {
        mAlignmentOutput->cd();
        mAlignmentInfo->Write();
    }


    /**
     * @brief Prepare QA histograms
     * 
     */
    void makeHistograms() {
        std::string n = "";
        mHist["ECalProjPosXY"] = new TH2F("ECalProjPosXY", ";X;Y", 1000, -500, 500, 1000, -500, 500);
        mHist["ECalProjSigmaXY"] = new TH2F("ECalProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 50, 0, 0.5, 50, 0, 0.5);
        mHist["ECalProjSigmaR"] = new TH1F("ECalProjSigmaR", ";#sigma_{XY} (cm) at ECAL", 50, 0, 0.5);

        mHist["SiProjPosXY"] = new TH2F("SiProjPosXY", ";X;Y", 1000, -500, 500, 1000, -500, 500);
        mHist["SiProjSigmaXY"] = new TH2F("SiProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 150, 0, 15, 150, 0, 15);

        mHist["VertexProjPosXY"] = new TH2F("VertexProjPosXY", ";X;Y", 100, -5, 5, 100, -5, 5);
        mHist["VertexProjSigmaXY"] = new TH2F("VertexProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 150, 0, 20, 150, 0, 20);

        mHist["VertexProjPosZ"] = new TH1F("VertexProjPosZ", ";Z;", 100, -50, 50);
        mHist["VertexProjSigmaZ"] = new TH1F("VertexProjSigmaZ", ";#sigma_{Z};", 100, 0, 10);

        mHist["SiWrongProjPosXY"] = new TH2F("SiWrongProjPosXY", ";X;Y", 1000, -500, 500, 1000, -500, 500);
        mHist["SiWrongProjSigmaXY"] = new TH2F("SiWrongProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 50, 0, 0.5, 50, 0, 0.5);

        mHist["SiDeltaProjPosXY"] = new TH2F("SiDeltaProjPosXY", ";X;Y", 1000, 0, 20, 1000, 0, 20);

        mHist["FstDiffZVsR"] = new TH2F( "FstDiffZVsR", ";R;dz", 400, 0, 40, 500, -5, 5 );
        mHist["FstDiffZVsPhiSliceInner"] = new TH2F( "FstDiffZVsPhiSliceInner", ";slice;dz", 15, 0, 15, 500, -5, 5 );
        mHist["FstDiffZVsPhiSliceOuter"] = new TH2F( "FstDiffZVsPhiSliceOuter", ";slice;dz", 15, 0, 15, 500, -5, 5 );

        mHist["FstDiffZVsPhiOuter"] = new TH2F( "FstDiffZVsPhiOuter", ";slice;dz", 628, 0, TMath::Pi()*2, 500, -5, 5 );

        mHist["CorrFstDiffZVsPhiSliceInner"] = new TH2F( "CorrFstDiffZVsPhiSliceInner", ";slice;dz", 15, 0, 15, 500, -5, 5 );
        mHist["CorrFstDiffZVsPhiSliceOuter"] = new TH2F( "CorrFstDiffZVsPhiSliceOuter", ";slice;dz", 15, 0, 15, 500, -5, 5 );


        n = "seed_curv";
        mHist[n] = new TH1F(n.c_str(), ";curv", 1000, 0, 10000);
        n = "seed_pT";
        mHist[n] = new TH1F(n.c_str(), ";pT (GeV/c)", 500, 0, 10);
        n = "seed_eta";
        mHist[n] = new TH1F(n.c_str(), ";eta", 500, 0, 5);

        n = "delta_fit_seed_pT";
        mHist[n] = new TH1F(n.c_str(), ";#Delta( fit, seed ) pT (GeV/c)", 500, -5, 5);
        n = "delta_fit_seed_eta";
        mHist[n] = new TH1F(n.c_str(), ";#Delta( fit, seed ) eta", 500, 0, 5);
        n = "delta_fit_seed_phi";
        mHist[n] = new TH1F(n.c_str(), ";#Delta( fit, seed ) phi", 500, -5, 5);

        n = "FitStatus";
        mHist[n] = new TH1F(n.c_str(), ";", 5, 0, 5);
        FwdTrackerUtils::labelAxis(mHist[n]->GetXaxis(), {"Total", "Pass", "Fail", "GoodCardinal", "Exception"});

        n = "FitDuration";
        mHist[n] = new TH1F(n.c_str(), "; Duraton (ms)", 5000, 0, 50000);

        n = "FailedFitDuration";
        mHist[n] = new TH1F(n.c_str(), "; Duraton (ms)", 500, 0, 50000);
    }

    /**
     * @brief writes histograms stored in map only if mGenHistograms is true
     * 
     */
    void writeHistograms() {
        if ( !mGenHistograms )
            return;
        for (auto nh : mHist) {
            nh.second->SetDirectory(gDirectory);
            nh.second->Write();
        }
    }

    /**
     * @brief Convert the 3x3 covmat to 2x2 by dropping z
     *
     * @param h : hit with cov matrix
     * @return TMatrixDSym : cov matrix 2x2
     */
    TMatrixDSym CovMatPlane(KiTrack::IHit *h){
        TMatrixDSym cm(2);
        cm(0, 0) = static_cast<FwdHit*>(h)->_covmat(0, 0);
        cm(1, 1) = static_cast<FwdHit*>(h)->_covmat(1, 1);
        cm(0, 1) = static_cast<FwdHit*>(h)->_covmat(0, 1);
        //cm.Print();
        return cm;
    }
    TMatrixDSym CovMatPlaneFst(TMatrixDSym input){
        TMatrixDSym cm(2);
        cm(0, 0) = input(0, 0);
        cm(1, 1) = input(1, 1);
        cm(0, 1) = input(0, 1);
        return cm;
    }


    /**
     * @brief Fit points to a simple circle
     *
     * @param trackSeed : seed points to fit
     * @param i0 : index of the first hit
     * @param i1 : index of the second hit
     * @param i2 : index of the third hit
     * @return float : curvature
     */
    float fitSimpleCircle(Seed_t trackSeed, size_t i0, size_t i1, size_t i2) {
        float curv = 0;

        // ensure that no index is outside of range for FST or FTT volumes
        if (i0 > 12 || i1 > 12 || i2 > 12)
            return 0;

        try {
            KiTrack::SimpleCircle sc(trackSeed[i0]->getX(), trackSeed[i0]->getY(), trackSeed[i1]->getX(), trackSeed[i1]->getY(), trackSeed[i2]->getX(), trackSeed[i2]->getY());
            curv = sc.getRadius();
        } catch (KiTrack::InvalidParameter &e) {
            // if we got here we failed to get  a valid seed. We will still try to move forward but the fit will probably fail
            LOG_WARN << e.what() << endl;
            LOG_WARN << "Circle fit failed, FWD track fit will likely fail" << endm;
        }

        //  make sure the curv is valid
        if (isinf(curv)){
            curv = 0.0;
        }

        return curv;
    } // fitSimpleCircle

    /**
     * @brief Determines the seed state to start the fit
     *
     * @param trackSeed : seed points
     * @param seedPos : position seed (return)
     * @param seedMom : momenum seed (return)
     * @return float : curvature
     */
    float seedState(Seed_t trackSeed, TVector3 &seedPos, TVector3 &seedMom) {
        // we require at least 3 hits, so this should be gauranteed
        LOG_DEBUG << "Seed state with " << trackSeed.size() << " seed points " << endm;
        if(trackSeed.size() < 3){
            // failure
            return 0.0;
        }

        // we want to use the LAST 3 hits, since silicon doesnt have R information
        TVector3 p0, p1, p2;
        // use the closest hit to the interaction point for the seed pos
        FwdHit *hit_closest_to_IP = static_cast<FwdHit *>(trackSeed[0]);

        // maps from <key=vol_id> to <value=index in trackSeed>
        std::map<size_t, int> vol_map;

        // init the map
        for (size_t i = 0; i < 13; i++)
            vol_map[i] = -1;

        vector<size_t> idx;
        for (size_t i = 0; i < trackSeed.size(); i++) {
            auto fwdHit = static_cast<FwdHit *>(trackSeed[i]);
            if ( !fwdHit ) continue;
            //cout << "fwdHit->_vid  = " << fwdHit->_vid << endl;
            if (vol_map[ abs(fwdHit->_vid) ] == -1)
                idx.push_back(fwdHit->_vid);
            vol_map[abs(fwdHit->_vid)] = (int)i;

            // find the hit closest to IP for the initial position seed
            if (hit_closest_to_IP == nullptr || hit_closest_to_IP->getZ() > fwdHit->getZ())
                hit_closest_to_IP = fwdHit;
        }

        // now get an estimate of the pT from several overlapping simple circle fits
        // enumerate the available partitions (example for FTT)
        // 12 11 10
        // 12 11 9
        // 12 10 9
        // 11 10 9
        vector<float> curvs;

        if (idx.size() < 3){
            return 0.0;
        }

        if ( idx.size() == 3 ){
            curvs.push_back(fitSimpleCircle(trackSeed, vol_map[idx[0]], vol_map[idx[1]], vol_map[idx[2]]));
        } else if ( idx.size() >= 4 ){
            curvs.push_back(fitSimpleCircle(trackSeed, vol_map[idx[0]], vol_map[idx[1]], vol_map[idx[2]]));
            curvs.push_back(fitSimpleCircle(trackSeed, vol_map[idx[0]], vol_map[idx[1]], vol_map[idx[3]]));
            curvs.push_back(fitSimpleCircle(trackSeed, vol_map[idx[0]], vol_map[idx[2]], vol_map[idx[3]]));
            curvs.push_back(fitSimpleCircle(trackSeed, vol_map[idx[1]], vol_map[idx[2]], vol_map[idx[3]]));
        }

        // average them and exclude failed fits
        float mcurv = 0;
        float nmeas = 0;

        for (size_t i = 0; i < curvs.size(); i++) {
            if (curvs[i] > 10) {
                mcurv += curvs[i];
                nmeas += 1.0;
            }
        }

        if (nmeas >= 1)
            mcurv = mcurv / nmeas;
        else
            mcurv = 10000;

        // Now lets get eta information
        p0.SetXYZ(trackSeed[vol_map[idx[0]]]->getX(), trackSeed[vol_map[idx[0]]]->getY(), trackSeed[vol_map[idx[0]]]->getZ());
        p1.SetXYZ(trackSeed[vol_map[idx[1]]]->getX(), trackSeed[vol_map[idx[1]]]->getY(), trackSeed[vol_map[idx[1]]]->getZ());
        if ( abs(p0.X() - p1.X()) < 1e-6 ){
            p1.SetXYZ(trackSeed[vol_map[idx[2]]]->getX(), trackSeed[vol_map[idx[2]]]->getY(), trackSeed[vol_map[idx[2]]]->getZ());
        }

        LOG_DEBUG << TString::Format( "Fwd SeedState: p0 (%f, %f, %f), p1 (%f, %f, %f)", p0.X(), p0.Y(), p0.Z(), p1.X(), p1.Y(), p1.Z() ) << endm;

        const double K = 0.00029979; //K depends on the units used for Bfield
        double pt = mcurv * K * 5; // pT from average measured curv
        double dx = (p1.X() - p0.X());
        double dy = (p1.Y() - p0.Y());
        double dz = (p1.Z() - p0.Z());
        double phi = TMath::ATan2(dy, dx);
        double Rxy = sqrt(dx * dx + dy * dy);
        double theta = TMath::ATan2(Rxy, dz);
        if (abs(dx) < 1e-6 || abs(dy) < 1e-6){
            phi = TMath::ATan2( p1.Y(), p1.X() );
        }
        Rxy = sqrt( p0.X()*p0.X() + p0.Y()*p0.Y() );
        theta = TMath::ATan2(Rxy, p0.Z());

        LOG_DEBUG << TString::Format( "pt=%f, dx=%f, dy=%f, dz=%f, phi=%f, theta=%f", pt, dx, dy, dz, phi, theta ) << endm;
        // double eta = -log( tantheta / 2.0 );
        // these starting conditions can probably be improvd, good study for students

        seedMom.SetPtThetaPhi(pt, theta, phi);
        seedPos.SetXYZ(hit_closest_to_IP->getX(), hit_closest_to_IP->getY(), hit_closest_to_IP->getZ());

        return mcurv;
    }//seedState

    /**
     * @brief Get track projection to given FST plane
     * 
     * @param fstPlane : fst sensor index
     * @param fitTrack : track to project
     * @return TVector2 
     */
    TVector2 projectToFst(size_t fstPlane, std::shared_ptr<genfit::Track> fitTrack) {
        if (fstPlane > 107) {
            TVector2 nil;
            return nil;
        }

        auto detSi = mFSTSensorPlanes[fstPlane];
        genfit::MeasuredStateOnPlane tst = fitTrack->getFittedState(1);
        auto TCM = fitTrack->getCardinalRep()->get6DCov(tst);

        //  can get the track length if needed
        fitTrack->getCardinalRep()->extrapolateToPlane(tst, detSi);

        TCM = fitTrack->getCardinalRep()->get6DCov(tst);

        // can get the projected positions if needed
        double x = tst.getPos().X();
        double y = tst.getPos().Y();
        double z = tst.getPos().Z();
        // and the uncertainties
        LOG_DEBUG << "Track Uncertainty at FST (plane=" << fstPlane << ") @ x= " << x << ", y= " << y << ", z= " << z << " : " << sqrt(TCM(0, 0)) << ", " << sqrt(TCM(1, 1)) << endm;

        return TVector2(x,y);
    }

    /**
     * @brief convert FST hit to a planar 2D hit
     * 
     * @param h : fst hit
     * @return TVector2 
     */
    TVector2 convertFstHitPlaneToLab(FwdHit *h){
        int fstPlane = h->getSensor();
        auto detSi = mFSTSensorPlanes[fstPlane];
        
        TVector3 hitLab = detSi->toLab(TVector2(h->getX(),h->getY()));

        return TVector2(hitLab.X(),hitLab.Y());    
    }


    /**
     * @brief Get track projection to given FTT plane
     * 
     * @param ftt_plane : plane index
     * @param fitTrack : track to project
     * @return TVector2 
     */
    TVector2 projectToFtt(size_t ftt_plane, std::shared_ptr<genfit::Track> fitTrack) {
        if (ftt_plane > 15) {
            TVector2 nil;
            return nil;
        }

        auto detFtt = mFTTPlanes[ftt_plane];
        genfit::MeasuredStateOnPlane tst = fitTrack->getFittedState(1);
        auto TCM = fitTrack->getCardinalRep()->get6DCov(tst);

        //  can get the track length if needed
        // double len = fitTrack->getCardinalRep()->extrapolateToPlane(tst, detSi, false, true);
        fitTrack->getCardinalRep()->extrapolateToPlane(tst, detFtt);

        TCM = fitTrack->getCardinalRep()->get6DCov(tst);

        // can get the projected positions if needed
        double x = tst.getPos().X();
        double y = tst.getPos().Y();
        double z = tst.getPos().Z();
        // and the uncertainties
        LOG_DEBUG << "Track Uncertainty at FTT (plane=" << ftt_plane << ") @ x= " << x << ", y= " << y << ", z= " << z << " : " << sqrt(TCM(0, 0)) << ", " << sqrt(TCM(1, 1)) << endm;

        return TVector2(x,y);
    }

    /**
     * @brief Get the Fst Plane object for a given hit
     *
     * @param h : hit
     * @return genfit::SharedPlanePtr
     */
    genfit::SharedPlanePtr getFstPlane( FwdHit * h ){

        size_t planeId  = h->getSector();
        size_t sensorId = h->getSensor();

        int sensorIdx = sensorId;

        if(!mMisaligned){
          size_t moduleId = mModuleMap[planeId][(int(sensorId)/3)%12];
          sensorId = mSensorMap[sensorId%3];

          sensorIdx = planeId*36 + moduleId*3 + sensorId;
        }

        auto planeCorr = mFSTSensorPlanes[sensorIdx];

        return planeCorr;

    } // GetFST PLANE

    /**
     * @brief Refit a track with additional FST hits
     *
     * Takes a previously fit track re-fits it with the newly added silicon hits
     * @param originalTrack : original fit track
     * @param fstHits : new FST hits to add
     * @return TVector3 : momentum
     */
    TVector3 refitTrackWithFstHits(genfit::Track *originalTrack, Seed_t fstHits) {
        static const TVector3 pOrig = originalTrack->getCardinalRep()->getMom(originalTrack->getFittedState(1, originalTrack->getCardinalRep()));

        if (originalTrack->getFitStatus(originalTrack->getCardinalRep())->isFitConverged() == false) {
            // in this case the original track did not converge so we should not refit.
            // probably never get here due to previous checks
            return pOrig;
        }

        // Setup the Track Reps
        auto trackRepPos = new genfit::RKTrackRep(mPdgPositron);
        auto trackRepNeg = new genfit::RKTrackRep(mPdgElectron);

        // get the space points on the original track
        auto trackPoints = originalTrack->getPointsWithMeasurement();

        if ((trackPoints.size() < (mFTTZLocations.size() +1) && mIncludeVertexInFit) || trackPoints.size() < mFTTZLocations.size() ) {
            // we didnt get enough points for a refit
            return pOrig;
        }

        // points are in reverse order, so n-1 is closest to IP (double check this)
        //TVectorD rawCoords = trackPoints[trackPoints.size()-1]->getRawMeasurement()->getRawHitCoords();
        //double z = mFTTZLocations[0]; //first FTT plane, used if we dont have PV in fit
        //if (mIncludeVertexInFit)
        //    z = rawCoords(2);

        //TVector3 seedPos(rawCoords(0), rawCoords(1), z);
        //TVector3 seedMom = pOrig;

        TVectorD rawCoords(3);
        if (!mIncludeVertexInFit)
        {
            int seedPlaneId = trackPoints[3]->getRawMeasurement()->getDetId();
            cout << "Seed Plane ID for silicon refit = " << seedPlaneId << endl;
            //double syshift = 5.9;
            //double sxshift = 0.0;
            //int squadId = -1;
            //squadId = seedPlaneId % 4;
            //if(squadId == 0) { sxshift =      + mCenterPent; syshift += mCenterPent; }
            //if(squadId == 1) { sxshift =      - mCenterPent; syshift += mCenterPent; }
            //if(squadId == 2) { sxshift = -6.5 - mCenterPent; syshift -= mCenterPent; }
            //if(squadId == 3) { sxshift = +6.5 + mCenterPent; syshift -= mCenterPent; }  

            TVectorD rawHitCoords = trackPoints[3]->getRawMeasurement()->getRawHitCoords();       
            double seedArr[4] = {rawHitCoords(0),rawHitCoords(1),rawHitCoords(2),1.0};
            TMatrixD seed4D(4,1,seedArr);
            seed4D = mInverseMStgc[seedPlaneId]*seed4D;

            double z = mFTTZLocations[0]; //first FTT plane, used if we dont have PV in fit
            rawCoords(0) = seed4D(0,0);// + sxshift;
            rawCoords(1) = seed4D(1,0);// + syshift;
            rawCoords(2) = seed4D(2,0);
        }
        if (mIncludeVertexInFit)
        {
            rawCoords = trackPoints[0]->getRawMeasurement()->getRawHitCoords();
        }

        TVector3 seedPos(rawCoords(0), rawCoords(1), rawCoords(2));
        TVector3 seedMom = pOrig;

        // Create the ref track using the seed state
        auto mFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
        mFitTrack->addTrackRep(trackRepNeg);

        size_t firstFTTIndex = 0;
        if (mIncludeVertexInFit) {
            // clone the PRIMARY VERTEX into this track
            mFitTrack->insertPoint(new genfit::TrackPoint(trackPoints[0]->getRawMeasurement(), mFitTrack));
            firstFTTIndex = 1; // start on hit index 1 below
        }

        // initialize the hit coords on plane
        TVectorD hitCoords(2);
        hitCoords[0] = 0;
        hitCoords[1] = 0;

        size_t planeId(0);
        int hitId(5);

        std::vector<double> siZ;
        for ( int infst = 0; infst < fstHits.size(); infst++ ) {
            auto h = fstHits[infst];
            if ( nullptr == h ) continue; // if no Si hit in this plane, skip

            int is = static_cast<FwdHit*>(h)->getSensor();

            double z;

            // FST disk (integer division rounds down)
            int d = is / 36; // 0-2

            // FST wedge
            int w = is / 3; // 0-35

            // FST sensor
            int s = is % 3; // 0 (inner), 1 (outer), 2 (outer)
            int ds = (s == 0)? 0 : 1; // +0 for inner, +1 for outer

            int defaultZidx = d * 4 + 2 * (w % 2) + ds;
            z = fstDefaultZ[defaultZidx];

            siZ.push_back(z);
        }

        std::vector<size_t> idx(siZ.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&siZ](size_t i1, size_t i2) {return siZ[i1] < siZ[i2];});

        // add the hits to the track
        for ( int infst = 0; infst < fstHits.size(); infst++ ) {
            auto h = fstHits[idx[infst]];
            if ( nullptr == h ) continue; // if no Si hit in this plane, skip

            LOG_DEBUG << "Hit index = " << infst << ",   Order index = " << idx[infst] << "   Z = " << siZ[idx[infst]] << endm;
                      
            planeId = h->getSector();
            auto plane = getFstPlane( static_cast<FwdHit*>(h) );
            int sensorId = static_cast<FwdHit*>(h)->getSensor();

            double hArr[4] = {h->getX(),h->getY(),h->getZ(),1.0};
            TMatrixD h4D(4,1,hArr);
            h4D = mInverseMFst[sensorId]*h4D;

            if(sensorId % 3 == 0) h4D(0,0) -= 10.75;
            if(sensorId % 3 != 0) h4D(0,0) -= 22.25;
            hitCoords[0] = h4D(0,0);
            hitCoords[1] = h4D(1,0);

            genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitCoords, CovMatPlane(h), sensorId+1000, ++hitId, nullptr);

            if (mFSTPlanes.size() <= planeId) {
                LOG_WARN << "invalid VolumId -> out of bounds DetPlane, vid = " << planeId << endm;
                return pOrig;
            }
            measurement->setPlane(plane, sensorId+1000);
            mFitTrack->insertPoint(new genfit::TrackPoint(measurement, mFitTrack));

            TVector3 hitXYZ( h->getX(), h->getY(), h->getZ() );
            float phi = hitXYZ.Phi();
            if ( phi < 0 ) phi = TMath::Pi() * 2 + phi;
            double phi_slice = phi / (TMath::Pi() / 6.0); // 2pi/12
            int phi_index = ((int)phi_slice);
            double dz = (h->getZ() - plane->getO().Z());

            double r  = sqrt( pow(hitXYZ.x(), 2) + pow(hitXYZ.y(), 2) );

            auto planeCorr = mFSTSensorPlanes[sensorId];
     
            double cdz = (h->getZ() - planeCorr->getO().Z());

            if (mGenHistograms){
                
                ((TH2*)mHist[ "FstDiffZVsR" ])->Fill( r, dz );

                if ( sensorId % 3 == 0 ) {// inner
                    mHist["FstDiffZVsPhiSliceInner"]->Fill( phi_slice, dz );
                    mHist["CorrFstDiffZVsPhiSliceInner"]->Fill( phi_slice, cdz );
                } else {
                    mHist["FstDiffZVsPhiSliceOuter"]->Fill( phi_slice, dz );
                    mHist["CorrFstDiffZVsPhiSliceOuter"]->Fill( phi_slice, cdz );
                    mHist["FstDiffZVsPhiOuter"]->Fill( phi, dz );
                }
            }

         }

        for (int iorder = 0; iorder < 4; iorder++) {
            for (size_t i = firstFTTIndex; i < trackPoints.size(); i++) {
                // clone the track points into this track
                planeId = trackPoints[i]->getRawMeasurement()->getDetId();
                if(int(planeId/4) != iorder) continue;
                mFitTrack->insertPoint(new genfit::TrackPoint(trackPoints[i]->getRawMeasurement(), mFitTrack));
            }
        }

        try {
            //Track RE-Fit with GENFIT2
            // check consistency of all points
            mFitTrack->checkConsistency();

            // do the actual track fit
            mFitter->processTrack(mFitTrack);

            mFitTrack->checkConsistency();

            // this chooses the lowest chi2 fit result as cardinal
            mFitTrack->determineCardinalRep();

        } catch (genfit::Exception &e) {
            // will be caught below by converge check
            LOG_WARN << "Track fit exception : " << e.what() << endm;
        }

        if (mFitTrack->getFitStatus(mFitTrack->getCardinalRep())->isFitConverged() == false) {
            // Did not converge
            cout << "Track fit did not converge" << endl;
            return pOrig;
        } else { // we did converge, return new momentum

            try {
                // causes seg fault
                // auto cardinalRep = mFitTrack->getCardinalRep();
                // auto cardinalStatus = mFitTrack->getFitStatus(cardinalRep);
                // mFitStatus = *cardinalStatus; // save the status of last fit
            } catch (genfit::Exception &e) {
            }

            static const TVector3 p = mFitTrack->getCardinalRep()->getMom(mFitTrack->getFittedState(1, mFitTrack->getCardinalRep()));
            
            if(mRefitGBL)
            {
              cout << "refit with GBL" << endl;
              mGblFitter->setSuccessfulFitFlag(false);
              refitTrackWithGBL(mFitTrack);
            }

            if(mGblFitter->getSuccessfulFitFlag())
            {
              LOG_INFO << "Successful GBL refit" << endl;
              int nhits = 0;
              for (auto h : fstHits)
              {
                if ( nullptr == h ) continue; // if no Si hit in this plane, skip
                nhits++;
              }
              for (auto h : fstHits)
              {
                if ( nullptr == h ) continue; // if no Si hit in this plane, skip

                int disk = h->getSector();
                int sid = static_cast<FwdHit*>(h)->getSensor();
                cout << "sid = " << sid << endl;

                auto plane = getFstPlane( static_cast<FwdHit*>(h) );
                TVector3 origin = plane->getO();
                double mcZ;
                double rcZ = origin.Z();

                TVectorD rcCoords(2);
                TVectorD mcCoords(2);
                TVectorD trackCoords(2);

                TVector3 rcCoordsLab;
                TVector3 mcCoordsLab;
                TVector3 trackCoordsLab;

                rcCoords[0] = h->getX();
                rcCoords[1] = h->getY();

                for (auto mch : static_cast<FwdHit*>(h)->_mcTrack->mFttHits)
                {

                  cout << "mch->getZ() = " << mch->getZ() << endl;
                  if(mch->getZ() > 200.0) continue;
                  if(static_cast<FwdHit*>(mch)->getSensor() != sid) continue;
                  double arr[4] = {mch->getX(),mch->getY(),mch->getZ(),1.0};
                  TMatrixD mc4D(4,1,arr);
                  mc4D = mInverseMFst[sid] * mc4D;

                  mcCoordsLab = TVector3(mch->getX(),mch->getY(),mch->getZ()); 

                  mcCoords[0] = mc4D(0,0);
                  mcCoords[1] = mc4D(1,0);
                  mcZ = mc4D(2,0);
                }


                cout << "mcZ = " << mcZ << "    rcZ = " << rcZ << "       mcZ-rcZ = " << mcZ-rcZ << endl;

                //double trackZ = 0.0;
                int stateid = -1;
                for( int tp = 0; tp < mFitTrack->getNumPointsWithMeasurement(); tp++)
                {
                  genfit::TrackPoint* point_meas_temp = mFitTrack->getPointWithMeasurement(tp);
                  genfit::PlanarMeasurement* measPlanar = dynamic_cast<genfit::PlanarMeasurement*>(point_meas_temp->getRawMeasurement(0));
                  if (measPlanar) stateid = measPlanar->getPlaneId();
                  if (stateid == sid + 1000)
                  {
                    //trackZ = mFSTSensorPlanes[sid]->getO()->Z();
                    genfit::KalmanFitterInfo* fi = dynamic_cast<genfit::KalmanFitterInfo*>(point_meas_temp->getFitterInfo(mFitTrack->getCardinalRep()));
                    genfit::ReferenceStateOnPlane* reference = new genfit::ReferenceStateOnPlane(*fi->getReferenceState());
                    TVectorD state = reference->getState();
                    genfit::AbsMeasurement* raw_meas = point_meas_temp->getRawMeasurement(0);
                    std::unique_ptr<const genfit::AbsHMatrix> HitHMatrix(raw_meas->constructHMatrix(mFitTrack->getCardinalRep()));
                    TVectorD planeState(HitHMatrix->Hv(state));

                    trackCoords[0] = planeState(0);
                    trackCoords[1] = planeState(1);

                    trackCoordsLab = mFSTSensorPlanes[sid]->toLab(TVector2(trackCoords[0],trackCoords[1])); 

                    break;
                  }
                }
                 
                double mcR = TMath::Sqrt(mcCoords[0]*mcCoords[0]+mcCoords[1]*mcCoords[1]);
                double mcP = TMath::ATan2(mcCoords[1],mcCoords[0]);
                double rcR = TMath::Sqrt(rcCoords[0]*rcCoords[0]+rcCoords[1]*rcCoords[1]);
                double rcP = TMath::ATan2(rcCoords[1],rcCoords[0]);
                double trackR = TMath::Sqrt(trackCoords[0]*trackCoords[0]+trackCoords[1]*trackCoords[1]);
                double trackP = TMath::ATan2(trackCoords[1],trackCoords[0]);

                //cout << "trackZ = " << trackZ << "    rcZ = " << rcZ << "       trackZ-rcZ = " << trackZ-rcZ << endl;
                if(TMath::Abs(mcR-rcR) > 2.875/2) cout << "MC R DOES NOT MATCH THE RC STRIP" << endl;

                double mcPhiStrip;
                int    mcPhiStripInt;
                double mcPhiStripPos;

                double rcPhiStrip;
                int    rcPhiStripInt;
                double rcPhiStripPos;

                double trackPhiStrip;
                int    trackPhiStripInt;
                double trackPhiStripPos;

                double phiGap = TMath::Pi() / 180. ;
                double stripWidth = TMath::Pi() / 6.0 / 128. ;

                if(sid % 3 == 0) // inner sensor
                {
                  //if(mcP >= 0.0) 
                  //{  
                  //  mcPhiStrip = mcP / stripWidth ; 
                  //  rcPhiStrip = rcP / stripWidth ; 
                  //  trackPhiStrip = trackP / stripWidth ; 
                  //}
                  //if(mcP < 0.0) 
                  //{
                  mcPhiStrip = (mcP+2*TMath::Pi()) / stripWidth ;
                  rcPhiStrip = (rcP+2*TMath::Pi()) / stripWidth ;
                  trackPhiStrip = (trackP+2*TMath::Pi()) / stripWidth ;
                  //}
                  mcPhiStripInt = int(mcPhiStrip) ;
                  mcPhiStripPos = (mcPhiStrip - double(mcPhiStripInt) - 0.5) * stripWidth;
                  cout << "mcP = " << mcP << "  stripWidth = " << stripWidth << "   mcPhiStrip = " << mcPhiStrip << "   mcPhiStripInt " << mcPhiStripInt << "   mcPhiStripPos = " << mcPhiStripPos << endl;

                  rcPhiStripInt = int(rcPhiStrip) ;
                  rcPhiStripPos = (rcPhiStrip - double(rcPhiStripInt) - 0.5) * stripWidth;
                  cout << "rcP = " << rcP << "  stripWidth = " << stripWidth << "   rcPhiStrip = " << rcPhiStrip << "   rcPhiStripInt " << rcPhiStripInt << "   rcPhiStripPos = " << rcPhiStripPos << endl;

                  trackPhiStripInt = int(trackPhiStrip) ;
                  trackPhiStripPos = (trackPhiStrip - double(rcPhiStripInt) - 0.5) * stripWidth;
                  cout << "trackP = " << trackP << "  stripWidth = " << stripWidth << "   trackPhiStrip = " << trackPhiStrip << "   trackPhiStripInt " << trackPhiStripInt << "   trackPhiStripPos = " << trackPhiStripPos << endl;

                }
                if(sid % 3 != 0 && mcP < 0.0)
                {
                  mcPhiStrip = (mcP+2*TMath::Pi()+phiGap/2.) / stripWidth ;
                  rcPhiStrip = (rcP+2*TMath::Pi()+phiGap/2.) / stripWidth ;
                  trackPhiStrip = (trackP+2*TMath::Pi()+phiGap/2.) / stripWidth ;

                  mcPhiStripInt = int(mcPhiStrip) ;
                  mcPhiStripPos = (mcPhiStrip - double(mcPhiStripInt) - 0.5) * stripWidth;
                  cout << "mcP = " << mcP << "  stripWidth = " << stripWidth << "   mcPhiStrip = " << mcPhiStrip << "   mcPhiStripInt " << mcPhiStripInt << "   mcPhiStripPos = " << mcPhiStripPos << endl;

                  rcPhiStripInt = int(rcPhiStrip) ;
                  rcPhiStripPos = (rcPhiStrip - double(rcPhiStripInt) - 0.5) * stripWidth;
                  cout << "rcP = " << rcP << "  stripWidth = " << stripWidth << "   rcPhiStrip = " << rcPhiStrip << "   rcPhiStripInt " << rcPhiStripInt << "   rcPhiStripPos = " << rcPhiStripPos << endl;

                  trackPhiStripInt = int(trackPhiStrip) ;
                  trackPhiStripPos = (trackPhiStrip - double(rcPhiStripInt) - 0.5) * stripWidth;
                  cout << "trackP = " << trackP << "  stripWidth = " << stripWidth << "   trackPhiStrip = " << trackPhiStrip << "   trackPhiStripInt " << trackPhiStripInt << "   trackPhiStripPos = " << trackPhiStripPos << endl;
                }
                if(sid % 3 != 0 && mcP > 0.0)
                {
                  mcPhiStrip = (mcP+2*TMath::Pi()-phiGap/2.) / stripWidth ;
                  rcPhiStrip = (rcP+2*TMath::Pi()-phiGap/2.) / stripWidth ;
                  trackPhiStrip = (trackP+2*TMath::Pi()-phiGap/2.) / stripWidth ;

                  mcPhiStripInt = int(mcPhiStrip) ;
                  mcPhiStripPos = (mcPhiStrip - double(mcPhiStripInt) - 0.5) * stripWidth;
                  cout << "mcP = " << mcP << "  stripWidth = " << stripWidth << "   mcPhiStrip = " << mcPhiStrip << "   mcPhiStripInt " << mcPhiStripInt << "   mcPhiStripPos = " << mcPhiStripPos << endl;

                  rcPhiStripInt = int(rcPhiStrip) ;
                  rcPhiStripPos = (rcPhiStrip - double(rcPhiStripInt) - 0.5) * stripWidth;
                  cout << "rcP = " << rcP << "  stripWidth = " << stripWidth << "   rcPhiStrip = " << rcPhiStrip << "   rcPhiStripInt " << rcPhiStripInt << "   rcPhiStripPos = " << rcPhiStripPos << endl;

                  trackPhiStripInt = int(trackPhiStrip) ;
                  trackPhiStripPos = (trackPhiStrip - double(rcPhiStripInt) - 0.5) * stripWidth;
                  cout << "trackP = " << trackP << "  stripWidth = " << stripWidth << "   trackPhiStrip = " << trackPhiStrip << "   trackPhiStripInt " << trackPhiStripInt << "   trackPhiStripPos = " << trackPhiStripPos << endl;
                }

                int is = -1;
                for(int i = 0; i < 4; i++)
                {
                  if(rcR > 5.0 + i * 2.875 && rcR <= 5.0 + (i + 1) * 2.875 && sid % 3 == 0 ) is = i;
                  if(rcR > 5.0 + (i + 4) * 2.875 && rcR <= 5.0 + (i + 5) * 2.875 && sid % 3 != 0 ) is = i + 4;
                }
                if(is < 0) continue;
                  
                float arr[15];
                int idx = 0;
                arr[idx++] = float(stateid);
                cout << "stateid = " << stateid << endl;
                arr[idx++] = float(mcCoords[0]);
                cout << "mcx = " << mcCoords[0] << endl;
                arr[idx++] = float(mcCoords[1]);
                cout << "mcy = " << mcCoords[1] << endl;
                arr[idx++] = float(mcZ);
                cout << "mcz = " << mcZ << endl;
                arr[idx++] = float(rcCoords[0]);
                cout << "rcx = " << rcCoords[0] << endl;
                arr[idx++] = float(rcCoords[1]);
                cout << "rcy = " << rcCoords[1] << endl;
                arr[idx++] = float(rcZ);
                cout << "rcz = " << rcZ << endl;
                arr[idx++] = float(trackCoords[0]);
                cout << "trackx = " << trackCoords[0] << endl;
                arr[idx++] = float(trackCoords[1]);
                cout << "tracky = " << trackCoords[1] << endl;
                arr[idx++] = float(mcCoordsLab[0]);
                cout << "mcxlab = " << mcCoordsLab[0] << endl;
                arr[idx++] = float(mcCoordsLab[1]);
                cout << "mcylab = " << mcCoordsLab[1] << endl;
                arr[idx++] = float(rcCoordsLab[0]);
                cout << "rcxlab = " << rcCoordsLab[0] << endl;
                arr[idx++] = float(rcCoordsLab[1]);
                cout << "rcylab = " << rcCoordsLab[1] << endl;
                arr[idx++] = float(trackCoordsLab[0]);
                cout << "trackxlab = " << trackCoordsLab[0] << endl;
                arr[idx++] = float(trackCoordsLab[1]);
                cout << "trackylab = " << trackCoordsLab[1] << endl;

                mAlignmentInfo->Fill(arr);
              }

              auto finalTrackPoints = mFitTrack->getPointsWithMeasurement();

              double mcZ;
              double rcZ;
              TVectorD rcCoords(2);
              TVectorD mcCoords(2);
              TVectorD trackCoords(2);

              for (auto mch : static_cast<FwdHit*>(fstHits[0])->_mcTrack->mFttHits)
              {
                if(mch->getZ() < 200.0) continue;
                int sector = mch->getSector();
                cout << "mch sector = " << sector << endl;

                int quadId = static_cast<FwdHit*>(mch)->getSensor();
                cout << "diskId = " << sector << ",    quadId = " << quadId << endl;
                int stgId = sector*4 + quadId;

                double xshift = 0.0;
                if(quadId == 2) xshift =  6.5;
                if(quadId == 3) xshift = -6.5;
                mcCoords[0] = mch->getX()+xshift;
                mcCoords[1] = mch->getY()-5.9;
                mcZ = mch->getZ();

                int stateid = -1;
                for( int tp = 0; tp < finalTrackPoints.size(); tp++)
                {
                  stateid = finalTrackPoints[tp]->getRawMeasurement()->getDetId();
                  if (stgId == stateid)
                  {
                    TVector3 rcCoordsLab;
                    TVector3 trackCoordsLab;

                    genfit::TrackPoint* point_meas_temp = finalTrackPoints[tp];
                    genfit::KalmanFitterInfo* fi = dynamic_cast<genfit::KalmanFitterInfo*>(point_meas_temp->getFitterInfo(mFitTrack->getCardinalRep()));
                    genfit::ReferenceStateOnPlane* reference = new genfit::ReferenceStateOnPlane(*fi->getReferenceState());
                    TVectorD state = reference->getState();
                    genfit::AbsMeasurement* raw_meas = point_meas_temp->getRawMeasurement(0);
                    rcCoords[0] = raw_meas->getRawHitCoords()[0];
                    rcCoords[1] = raw_meas->getRawHitCoords()[1];

                    rcCoordsLab = mFTTPlanes[stateid]->toLab(TVector2(rcCoords[0],rcCoords[1])); 

                    std::unique_ptr<const genfit::AbsHMatrix> HitHMatrix(raw_meas->constructHMatrix(mFitTrack->getCardinalRep()));
                    TVectorD planeState(HitHMatrix->Hv(state));

                    trackCoords[0] = planeState(0);
                    trackCoords[1] = planeState(1);
  
                    trackCoordsLab = mFTTPlanes[stateid]->toLab(TVector2(trackCoords[0],trackCoords[1])); 

                    //trackZ = mFSTSensorPlanes[sid]->getO()->Z();
                    rcZ = mFTTPlanes[stgId]->getO().Z();

                    float arr[15];
                    int idx = 0;
                    arr[idx++] = float(stateid);
                    cout << "stateid = " << stateid << endl;
                    arr[idx++] = float(mcCoords[0]);
                    cout << "mcx = " << mcCoords[0] << endl;
                    arr[idx++] = float(mcCoords[1]);
                    cout << "mcy = " << mcCoords[1] << endl;
                    arr[idx++] = float(mcZ);
                    cout << "mcz = " << mcZ << endl;
                    arr[idx++] = float(rcCoords[0]);
                    cout << "rcx = " << rcCoords[0] << endl;
                    arr[idx++] = float(rcCoords[1]);
                    cout << "rcy = " << rcCoords[1] << endl;
                    arr[idx++] = float(rcZ);
                    cout << "rcz = " << rcZ << endl;
                    arr[idx++] = float(trackCoords[0]);
                    cout << "trackx = " << trackCoords[0] << endl;
                    arr[idx++] = float(trackCoords[1]);
                    cout << "tracky = " << trackCoords[1] << endl;
                    arr[idx++] = float(mch->getX());
                    cout << "mcxlab = " << mch->getX() << endl;
                    arr[idx++] = float(mch->getY());
                    cout << "mcylab = " << mch->getY() << endl;
                    arr[idx++] = float(rcCoordsLab[0]);
                    cout << "rcxlab = " << rcCoordsLab[0] << endl;
                    arr[idx++] = float(rcCoordsLab[1]);
                    cout << "rcylab = " << rcCoordsLab[1] << endl;
                    arr[idx++] = float(trackCoordsLab[0]);
                    cout << "trackxlab = " << trackCoordsLab[0] << endl;
                    arr[idx++] = float(trackCoordsLab[1]);
                    cout << "trackylab = " << trackCoordsLab[1] << endl;

                    mAlignmentInfo->Fill(arr);

                    break;
                  }
                }
              }
            }
            return p;
        }
        return pOrig;
    } // refit with Si hits

    //TVector3 refitTrackWithFstHits(genfit::Track *originalTrack, Seed_t fstHits) {
    //    TVector3 pOrig = originalTrack->getCardinalRep()->getMom(originalTrack->getFittedState(1, originalTrack->getCardinalRep()));

    //    if (originalTrack->getFitStatus(originalTrack->getCardinalRep())->isFitConverged() == false) {
    //        // in this case the original track did not converge so we should not refit. 
    //        // probably never get here due to previous checks
    //        return pOrig;
    //    }

    //    // Setup the Track Reps
    //    auto trackRepPos = new genfit::RKTrackRep(mPdgPositron);
    //    auto trackRepNeg = new genfit::RKTrackRep(mPdgElectron);

    //    // get the space points on the original track
    //    auto trackPoints = originalTrack->getPointsWithMeasurement();
    //    
    //    if ((trackPoints.size() < (mFTTZLocations.size() +1) && mIncludeVertexInFit) || trackPoints.size() < mFTTZLocations.size() ) {
    //        // we didnt get enough points for a refit
    //        return pOrig;
    //    }

    //    TVectorD rawCoords = trackPoints[0]->getRawMeasurement()->getRawHitCoords();
    //    double z = mFSTZLocations[0]; //first FTT plane, used if we dont have PV in fit
    //    if (mIncludeVertexInFit)
    //        z = rawCoords(2);

    //    TVector3 seedPos(rawCoords(0), rawCoords(1), z);
    //    TVector3 seedMom = pOrig;

    //    // Create the ref track using the seed state
    //    auto mFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
    //    mFitTrack->addTrackRep(trackRepNeg);

    //    genfit::Track &fitTrack = *mFitTrack;

    //    size_t firstFTTIndex = 0;
    //    if (mIncludeVertexInFit) {
    //        // clone the PRIMARY VERTEX into this track
    //        fitTrack.insertPoint(new genfit::TrackPoint(trackPoints[0]->getRawMeasurement(), &fitTrack));
    //        firstFTTIndex = 1; // start on hit index 1 below
    //    }

    //    // initialize the hit coords on plane
    //    TVectorD hitCoords(2);
    //    hitCoords[0] = 0;
    //    hitCoords[1] = 0;

    //    size_t planeId(0);
    //    int hitId(5);

    //    // add the hits to the track
    //    for (auto h : fstHits) {
    //        if ( nullptr == h ) continue; // if no Si hit in this plane, skip

    //        hitCoords[0] = h->getX();
    //        hitCoords[1] = h->getY();
    //        genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitCoords, CovMatPlane(h), h->getSector(), ++hitId, nullptr);

    //        planeId = h->getSector();

    //        if (mFSTPlanes.size() <= planeId) {
    //            LOG_WARN << "invalid VolumId -> out of bounds DetPlane, vid = " << planeId << endm;
    //            return pOrig;
    //        }

    //        auto plane = getFstPlane( static_cast<FwdHit*>(h) );

    //        measurement->setPlane(plane, planeId);
    //        fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));


    //        TVector3 hitXYZ( h->getX(), h->getY(), h->getZ() );
    //        float phi = hitXYZ.Phi();
    //        if ( phi < 0 ) phi = TMath::Pi() * 2 + phi;
    //        double phi_slice = phi / (TMath::Pi() / 6.0); // 2pi/12
    //        int phi_index = ((int)phi_slice);
    //        double dz = (h->getZ() - plane->getO().Z());

    //        double r  =sqrt( pow(hitXYZ.x(), 2) + pow(hitXYZ.y(), 2) );

    //        size_t idx = phi_index % 2;
    //        auto planeCorr = mFSTPlanesInner[planeId + idx];
    //        if ( r > 16 ){
    //            planeCorr = mFSTPlanesOuter[planeId + idx];
    //        }
    //        double cdz = (h->getZ() - planeCorr->getO().Z());

    //        if (mGenHistograms){
    //            ((TH2*)mHist[ "FstDiffZVsR" ])->Fill( r, dz );

    //            if ( r < 16 ) {// inner
    //                mHist["FstDiffZVsPhiSliceInner"]->Fill( phi_slice, dz );
    //                mHist["CorrFstDiffZVsPhiSliceInner"]->Fill( phi_slice, cdz );
    //            } else {
    //                mHist["FstDiffZVsPhiSliceOuter"]->Fill( phi_slice, dz );
    //                mHist["CorrFstDiffZVsPhiSliceOuter"]->Fill( phi_slice, cdz );
    //                mHist["FstDiffZVsPhiOuter"]->Fill( phi, dz );
    //            }
    //        } // gen histograms 
    //    } // for fstHits
    //    // start at 0 if PV not included, 1 otherwise 
    //    for (size_t i = firstFTTIndex; i < trackPoints.size(); i++) {
    //        // clone the track points into this track
    //        fitTrack.insertPoint(new genfit::TrackPoint(trackPoints[i]->getRawMeasurement(), &fitTrack));
    //    }

    //    try {
    //        //Track RE-Fit with GENFIT2
    //        // check consistency of all points
    //        fitTrack.checkConsistency();

    //        // do the actual track fit
    //        mFitter->processTrack(&fitTrack);

    //        fitTrack.checkConsistency();

    //        // this chooses the lowest chi2 fit result as cardinal
    //        fitTrack.determineCardinalRep(); 

    //    } catch (genfit::Exception &e) {
    //        // will be caught below by converge check
    //        LOG_WARN << "Track fit exception : " << e.what() << endm;
    //    }

    //    if (fitTrack.getFitStatus(fitTrack.getCardinalRep())->isFitConverged() == false) {
    //        // Did not converge
    //        return pOrig;
    //    } else { // we did converge, return new momentum
    //        
    //        try {
    //            // causes seg fault
    //            auto cardinalRep = fitTrack.getCardinalRep();
    //            auto cardinalStatus = fitTrack.getFitStatus(cardinalRep);
    //            mFitStatus = *cardinalStatus; // save the status of last fit
    //        } catch (genfit::Exception &e) {
    //        }

    //        TVector3 p = fitTrack.getCardinalRep()->getMom(fitTrack.getFittedState(1, fitTrack.getCardinalRep()));
    //        return p;
    //    }
    //    return pOrig;
    //} // refit with Si hits

    /**
     * @brief Refit a track with GBL method
     * 
     * Takes a previously fit track and fits it with GBL method and provides output to alignment data files 
     * @param originalTrack : original fit track
     * @return TVector3 : momentum
     */
    void refitTrackWithGBL( genfit::Track *originalTrack ) {
        // mem leak, global track is overwritten without delete.
        static const TVector3 pOrig = originalTrack->getCardinalRep()->getMom(originalTrack->getFittedState(1, originalTrack->getCardinalRep()));
        
        // auto cardinalStatus = originalTrack->getFitStatus(originalTrack->getCardinalRep());

        if (originalTrack->getFitStatus(originalTrack->getCardinalRep())->isFitConverged() == false) {
            // in this case the original track did not converge so we should not refit. 
            // probably never get here due to previous checks
            return;// pOrig;
        }

        try {
            // check consistency of all points
            originalTrack->checkConsistency();

            // do the actual track fit
            mGblFitter->processTrackWithRep(originalTrack, originalTrack->getCardinalRep());

            originalTrack->checkConsistency();

            // this chooses the lowest chi2 fit result as cardinal
            originalTrack->determineCardinalRep();
 
          } catch (genfit::Exception &e) {
             // will be caught below by converge check
             LOG_WARN << "Track fit exception : " << e.what() << endm;
         }
 
         if (originalTrack->getFitStatus(originalTrack->getCardinalRep())->isFitConverged() == false) {
             LOG_WARN << "GBL fit did not converge" << endm;
             return;//pOrig;
         } else { // we did converge, return new momentum

            try {
                // causes seg fault
                auto cardinalRep = originalTrack->getCardinalRep();
                auto cardinalStatus = originalTrack->getFitStatus(cardinalRep);
                mFitStatus = *cardinalStatus; // save the status of last fit
            } catch (genfit::Exception &e) {
                LOG_WARN << "Failed to get cardinal status from converged fit" << endm;
            }

            return; //originalTrack->getCardinalRep()->getMom(originalTrack->getFittedState(1, originalTrack->getCardinalRep()));
        }
        return; //pOrig;
    } //refitwith GBL

    void fitSpacePoints( Seed_t trackSeed, double *Vertex = 0, TVector3 *seedMomentum = 0 ){

        auto seedPos = TVector3(0, 0, 0);
        auto seedMom = TVector3(0, 0, 10);
        auto theTrackRep = new genfit::RKTrackRep(mPdgMuon);
        // setup track for fit with positive and negative reps
        // mFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
        mFitTrack = std::make_shared<genfit::Track>(theTrackRep, seedPos, seedMom);

        // Create Spaces points from each seed
        vector<genfit::SpacepointMeasurement*> spoints;
        for (size_t i = 0; i < trackSeed.size(); i++) {
            auto seed = trackSeed[i];
            TMatrixDSym cm(3);
            cm(0, 0) = 0.01;
            cm(1, 1) = 0.01;
            cm(2, 2) = 0.01;
            auto rhc = TVectorD( 3 );
            rhc[0] = seed->getX();
            rhc[1] = seed->getY();
            rhc[2] = seed->getZ();
            auto spoint = new genfit::SpacepointMeasurement(rhc, cm, 0, i, nullptr);
            spoints.push_back(spoint);
        }

        LOG_INFO << "Fitting track with " << spoints.size() << " space points" << endm;
        // try adding the points to track and fitting
        try {
            for ( size_t i = 0; i < spoints.size(); i++ ){
                mFitTrack->insertPoint(new genfit::TrackPoint(spoints[i], mFitTrack.get()));
            }

            mFitter->processTrack(mFitTrack.get());
            mFitTrack->determineCardinalRep();

        } catch (genfit::Exception &e) {
            LOG_ERROR << "GenFit failed to fit track with: " << e.what() << endm;
        }

        try {
            mFitTrack->checkConsistency();

            mFitTrack->determineCardinalRep();
            auto cardinalRep = mFitTrack->getCardinalRep();
            // static const TVector3 p = cardinalRep->getMom(mFitTrack->getFittedState(1, cardinalRep));
            // sucess, return momentum
            LOG_INFO << "Fit status: " << mFitTrack->getFitStatus(cardinalRep)->isFitConverged() << endm;
            return;
        } catch (genfit::Exception &e) {
            LOG_ERROR << "GenFit failed to fit track consistency with: " << e.what() << endm;
        }
        return ;
    }

    /**
     * @brief Generic method for fitting space points with GenFit
     *
     * @param spoints : spacepoints
     * @param seedPos : seed position
     * @param seedMom : seed momentum
     * @return TVector3 : momentum from fit
     */
    TVector3 fitSpacePoints( vector<genfit::SpacepointMeasurement*> spoints, TVector3 &seedPos, TVector3 &seedMom ){

        // setup track reps
        // auto trackRepPos = new genfit::RKTrackRep(mPdgPositron);
        // auto trackRepNeg = new genfit::RKTrackRep(mPdgElectron);


        // return TVector3(0, 0, 0);
        auto theTrackRep = new genfit::RKTrackRep(mPdgMuon);
        // setup track for fit with positive and negative reps
        // mFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
        mFitTrack = std::make_shared<genfit::Track>(theTrackRep, seedPos, seedMom);
        // mFitTrack->addTrackRep(trackRepNeg);
        //

        LOG_INFO << "Fitting track with " << spoints.size() << " space points" << endm;
        // try adding the points to track and fitting
        try {
            for ( size_t i = 0; i < spoints.size(); i++ ){
                mFitTrack->insertPoint(new genfit::TrackPoint(spoints[i], mFitTrack.get()));
            }
            // do the fit against the two possible fits
            // mFitter->processTrackWithRep(mFitTrack, trackRepPos);
            // mFitter->processTrackWithRep(mFitTrack, trackRepNeg);
            mFitter->processTrack(mFitTrack.get());
            mFitTrack->determineCardinalRep();

        } catch (genfit::Exception &e) {
            LOG_ERROR << "GenFit failed to fit track with: " << e.what() << endm;
        }

        try {
            mFitTrack->checkConsistency();

            mFitTrack->determineCardinalRep();
            auto cardinalRep = mFitTrack->getCardinalRep();
            static const TVector3 p = cardinalRep->getMom(mFitTrack->getFittedState(1, cardinalRep));
            // sucess, return momentum
            return p;
        } catch (genfit::Exception &e) {
            LOG_ERROR << "GenFit failed to fit track with: " << e.what() << endm;
        }
        return TVector3(0, 0, 0);
    }


    /**
     * @brief setup the track from the given seed and optional primary vertex
     * @param trackSeed : seed points
     * @param seedMom : seed momentum
     * @param seedPos : seed position
     * @param Vertex : primary vertex
     */
    void setupTrack(Seed_t trackSeed, TVector3 seedMom, TVector3 seedPos, double *Vertex = 0) {
        // If we use the PV, use that as the start pos for the track
        TVectorD pv(3);
        if (Vertex != 0) {
            seedPos.SetXYZ(Vertex[0], Vertex[1], Vertex[2]);
            pv[0] = Vertex[0];
            pv[1] = Vertex[1];
            pv[2] = Vertex[2];
        }

        LOG_INFO << "Vertex position: " << endm;
        pv.Print();
        
        // create the track representations
        // Note that multiple track reps differing only by charge results in a silent failure of GenFit
        auto theTrackRep = new genfit::RKTrackRep(mPdgMuon);

        // Create the track
        mFitTrack = std::make_shared<genfit::Track>(theTrackRep, seedPos, seedMom);
        // now add the points to the track

        size_t planeId(0);     // detector plane ID
        size_t quadId(0);      // detector quad ID
        int hitId(0);       // hit ID
        /******************************************************************************************************************
        * Include the Primary vertex if desired
        ******************************************************************************************************************/
        if (mIncludeVertexInFit) {
            LOG_DEBUG << "Including vertex in fit" << endm;
            TMatrixDSym hitCov3(3);
            hitCov3(0, 0) = mVertexSigmaXY * mVertexSigmaXY;
            hitCov3(1, 1) = mVertexSigmaXY * mVertexSigmaXY;
            hitCov3(2, 2) = mVertexSigmaZ * mVertexSigmaZ;

            genfit::SpacepointMeasurement *measurement = new genfit::SpacepointMeasurement(pv, hitCov3, 9999, ++hitId, nullptr);
            mFitTrack->insertPoint(new genfit::TrackPoint(measurement, mFitTrack.get()));
            //LOG_INFO << "Added vertex to track" << endm;
        }


        // initialize the hit coords on plane
        TVectorD hitCoords(2);
        hitCoords[0] = 0;
        hitCoords[1] = 0;

        /******************************************************************************************************************
	* sort hits to add by their z-location
 	******************************************************************************************************************/
        //LOG_INFO << "About to sort hits by their z-location" << endm;
        std::vector<double> hitZ;
        for (auto h : trackSeed) {
            if ( nullptr == h ) continue; // if no Si hit in this plane, skip
            //cout << "static_cast<FwdHit*>(h)->_covmat(0, 0) = " <<  static_cast<FwdHit*>(h)->_covmat(0, 0) << endl;
            //cout << "static_cast<FwdHit*>(h)->_covmat(0, 1) = " <<  static_cast<FwdHit*>(h)->_covmat(0, 1) << endl;
            //cout << "static_cast<FwdHit*>(h)->_covmat(1, 1) = " <<  static_cast<FwdHit*>(h)->_covmat(1, 1) << endl;

            double z;

            if(h->getZ() < 200.0)
            {
                int is = static_cast<FwdHit*>(h)->getSensor();
                //LOG_INFO << "is = " << is << endm;
                // FST disk (integer division rounds down)
                int d = is / 36; // 0-2
                //LOG_INFO << "d = " << d << endm;

                // FST wedge
                int w = is / 3; // 0-35
                //LOG_INFO << "w = " << w << endm;

                // FST sensor
                int s = is % 3; // 0 (inner), 1 (outer), 2 (outer)
                //LOG_INFO << "s = " << s << endm;
                int ds = (s == 0)? 0 : 1; // +0 for inner, +1 for outer

                int defaultZidx = d * 4 + 2 * (w % 2) + ds;
                //LOG_INFO << "defaultZidx = " << defaultZidx << endm;
                z = fstDefaultZ[defaultZidx];
            }
            else
            {
                size_t diskId = h->getSector();
                z = stgcDefaultZ[diskId];
            }
            hitZ.push_back(z);
        }

        std::vector<size_t> idx(hitZ.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&hitZ](size_t i1, size_t i2) {return hitZ[i1] < hitZ[i2];});
        
        //LOG_INFO << "Finished sorting" << endm;
      
        /******************************************************************************************************************
	* loop over the hits, add them to the track
 	******************************************************************************************************************/
        //LOG_INFO << "Loop over the hits and add them to the track" << endm;
        for ( int ih = 0; ih < trackSeed.size(); ih++ ) {
            auto h = trackSeed[idx[ih]];
            size_t diskId = h->getSector();

            double z = h->getZ();
            //cout << "Z = " << z << endl;
            bool isFst = z < 200;

            TMatrixDSym covMatFst(2); 
            double hArr[4] = {h->getX(),h->getY(),h->getZ(),1.0};
            TMatrixD h4D(4,1,hArr);
            if(!isFst) {
                quadId = static_cast<FwdHit*>(h)->getSensor();
                planeId = diskId*4 + quadId;
                h4D = mInverseMStgc[planeId]*h4D;
                switch(quadId) {
                  case 0 : // +x,+y
                    h4D(0,0) -= mCenterPent;
                    h4D(1,0) -= (5.9+mCenterPent);
                    break;
                  case 1 : // -x,+y
                    h4D(0,0) -= -mCenterPent;
                    h4D(1,0) -= (5.9+mCenterPent);
                    break;
                  case 2 : // -x,-y
                    h4D(0,0) -= (-6.5-mCenterPent);
                    h4D(1,0) -= (5.9-mCenterPent);
                    break;
                  case 3 : // +x,-y
                    h4D(0,0) -= (6.5+mCenterPent);
                    h4D(1,0) -= (5.9-mCenterPent);
                    break;
                }
            }
            if(isFst) {                
                planeId = static_cast<FwdHit*>(h)->getSensor();
                h4D = mInverseMFst[planeId]*h4D;
                //cout << "x = " << h4D(0,0) << endl;
                //cout << "y = " << h4D(1,0) << endl;
                //cout << "z = " << h4D(2,0) << endl;
                //covMatFst = makeSiCovMat(TVector3(h4D(0,0),h4D(1,0),h4D(2,0)));
                //covMatFst.Print();
                int s = planeId%3; 
                if(s == 0) h4D(0,0) -= 10.75;
                if(s != 0) h4D(0,0) -= 22.25;
            }

            hitCoords[0] = h4D(0,0);
            hitCoords[1] = h4D(1,0);

            //cout << "static_cast<FwdHit*>(h)->_covmat(0, 0) = " <<  static_cast<FwdHit*>(h)->_covmat(0, 0) << endl;
            //cout << "static_cast<FwdHit*>(h)->_covmat(0, 1) = " <<  static_cast<FwdHit*>(h)->_covmat(0, 1) << endl;
            //cout << "static_cast<FwdHit*>(h)->_covmat(1, 1) = " <<  static_cast<FwdHit*>(h)->_covmat(1, 1) << endl;


            genfit::PlanarMeasurement *measurement;
            if(!isFst)
            {
              //cout << "About to Create a planar measurement for FTT hit" << endl;
              //CovMatPlane(h);
              //cout << "Called CovMatPlane(h);" << endl;
              measurement = new genfit::PlanarMeasurement(hitCoords, CovMatPlane(h), planeId, ++hitId, nullptr);
              //cout << "Created a planar measurement for FTT hit!" << endl;
            }
            if(isFst ) measurement = new genfit::PlanarMeasurement(hitCoords, makeSiCovMat(TVector3(h4D(0,0),h4D(1,0),h4D(2,0)))/*covMatFst*//*CovMatPlaneFst(covMatFst)*/, planeId+1000, ++hitId, nullptr);
             
            //cout << "planeId = " << planeId << endl;
            //cout << "diskId = " << diskId << endl;

            //cout << "Before selecting plane" << endl;
            genfit::SharedPlanePtr plane;
            if(!isFst) plane = mFTTPlanes[planeId];
            if(isFst ) {
                plane = mFSTSensorPlanes[planeId];
                planeId += 1000; // offset for FST plane ID for alignment parameter ID
            }
            measurement->setPlane(plane, planeId);
            mFitTrack->insertPoint(new genfit::TrackPoint(measurement, mFitTrack.get()));

            if (abs(h->getZ() - plane->getO().Z()) > 0.05) {
                LOG_WARN << "Z Mismatch h->z = " << h->getZ() << ", plane->z = "<< plane->getO().Z() <<", diff = " << h->getZ() - plane->getO().Z() << endm;
            }
        } // loop on trackSeed
    } // setupTrack

    /** @brief performs the fit on a track
     *  @param t : track to fit
    */
    void performFit( std::shared_ptr<genfit::Track> t ){
        /******************************************************************************************************************
		 * Do the fit
		 ******************************************************************************************************************/
        try {
            // do the fit
            mFitter->processTrack(t.get());
            // find track rep with smallest chi2
            t->determineCardinalRep();
            // update the seed
            t->udpateSeed();
        } catch (genfit::Exception &e) {
            LOG_ERROR << "Exception on fit update" << e.what() << endm;
        }
        LOG_INFO << "Track fit update complete!" << endm;
    }
    /**
     * @brief Primary track fitting routine
     *
     * @param trackSeed :
     * @param Vertex : Primary Vertex
     * @param seedMomentum : seed momentum (can be from MC)
     * @return void : the results can be accessed via the getTrack() method
     */
    long long fitTrack(Seed_t trackSeed, double *Vertex = 0, TVector3 *seedMomentum = 0, bool gblRefit = false) {
        long long itStart = FwdTrackerUtils::nowNanoSecond();
        LOG_DEBUG << "Fitting track with " << trackSeed.size() << " seed points" << endm;
        
        // get the seed info from our hits
        static TVector3 seedMom, seedPos;
        LOG_DEBUG << "Getting seed state" << endm;
        // returns track curvature if needed
        seedState(trackSeed, seedPos, seedMom);

        if (seedMomentum != nullptr) {
            seedMom = *seedMomentum;
            LOG_DEBUG << "Using provided seedMomentum: " << TString::Format( "(pt=%f, eta=%f, phi=%f)", seedMom.Pt(), seedMom.Eta(), seedMom.Phi() ) << endm;
        }

        setupTrack(trackSeed, seedMom, seedPos, Vertex);
        LOG_INFO << "Ready to fit with " << mFitTrack->getNumPoints() << " track points" << endm;
        //LOG_DEBUG << "Ready to fit with " << mFitTrack->getNumPoints() << " track points" << endm;

        /******************************************************************************************************************
		 * Do the fit
		 ******************************************************************************************************************/
        performFit( mFitTrack );

        genfit::Track gblTrack(*mFitTrack); 

        LOG_INFO << "About to check if we should refit with GBL " << endm;
        //cout << "mRefitGBL = " << mRefitGBL << endl;
        //cout << "gblRefit = " << gblRefit << endl;
        if(mRefitGBL && gblRefit)
        {
          LOG_INFO << "We are going to attempt GBL refit" << endm;
          mGblFitter->setSuccessfulFitFlag(false);
          //refitTrackWithGBL(&gblTrack);
          refitTrackWithGBL(&*mFitTrack);
        
          LOG_INFO << "refitTrackWithGBL() was processed " << endm;
          //cout << "Tried to refit track with GBL " << endl;
          if(mGblFitter->getSuccessfulFitFlag())
          {
            LOG_INFO << "Successful GBL refit" << endm;
            int stateid = -1;
            for( int tp = 0; tp < gblTrack.getNumPointsWithMeasurement(); tp++)
            {
              genfit::TrackPoint* point_meas_temp = gblTrack.getPointWithMeasurement(tp);
              genfit::PlanarMeasurement* measPlanar = dynamic_cast<genfit::PlanarMeasurement*>(point_meas_temp->getRawMeasurement(0));
              if (measPlanar) stateid = measPlanar->getPlaneId();
              //trackZ = mFSTSensorPlanes[sid]->getO()->Z();
              genfit::KalmanFitterInfo* fi = dynamic_cast<genfit::KalmanFitterInfo*>(point_meas_temp->getFitterInfo(gblTrack.getCardinalRep()));
              genfit::ReferenceStateOnPlane* reference = new genfit::ReferenceStateOnPlane(*fi->getReferenceState());
              TVectorD state = reference->getState();
              genfit::AbsMeasurement* raw_meas = point_meas_temp->getRawMeasurement(0);
              std::unique_ptr<const genfit::AbsHMatrix> HitHMatrix(raw_meas->constructHMatrix(gblTrack.getCardinalRep()));
              TVectorD planeState(HitHMatrix->Hv(state));

              TVectorD rcCoords(2);
              TVectorD trackCoords(2);
              rcCoords[0] = raw_meas->getRawHitCoords()[0];
              rcCoords[1] = raw_meas->getRawHitCoords()[1];
              
              float x = rcCoords[0];
              float y = rcCoords[1];
              if(stateid >= 1000 && stateid < 9000) {
                if((stateid-1000)%3 == 0) x+=10.75;
                if((stateid-1000)%3 != 0) x+=22.25;
                float r = TMath::Sqrt(x*x+y*y);
                if(stateid >= 1000 && stateid < 1036 && r >= 22.250) continue; 
                if(stateid >= 1036 && stateid < 1072 && r <=  7.875) continue; 
                if(stateid >= 1036 && stateid < 1072 && r >= 25.125) continue; 
                if(stateid >= 1072 && stateid < 1108 && r <=  7.875) continue; 
              }

              trackCoords[0] = planeState(0);               
              trackCoords[1] = planeState(1);
 
 
              TVector2 rcCoordsLocal(rcCoords[0],rcCoords[1]);
              TVector2 trackCoordsLocal(trackCoords[0],trackCoords[1]);

              TVector3 rcCoordsLab;
              TVector3 trackCoordsLab;

              if(stateid < 0 && stateid > 9000) continue;

              if(stateid >= 0 && stateid < 20) 
              {
                 rcCoordsLab = mFTTPlanes[stateid]->toLab(rcCoordsLocal);
                 trackCoordsLab = mFTTPlanes[stateid]->toLab(trackCoordsLocal);
              }
              if(stateid > 20 && stateid < 9000)
              {
                 rcCoordsLab = mFSTSensorPlanes[stateid-1000]->toLab(rcCoordsLocal);
                 trackCoordsLab = mFSTSensorPlanes[stateid-1000]->toLab(trackCoordsLocal);
              }
            
              TVector3 mcCoordsLab;
              TVectorD mcCoords(2); 
              double mcZ = 0.0; 

              //if(mMcTracking)
              //{
              //  auto h = trackSeed[idx[0]];
              //  
              //  for (auto mch : static_cast<FwdHit*>(h)->_mcTrack->mFstHits)
              //  {
              //    if(mch->getZ() < 200.0)
              //    { 
              //      if(static_cast<FwdHit*>(mch)->getSensor() != stateid-1000) continue;
              //      double arr[4] = {mch->getX(),mch->getY(),mch->getZ(),1.0};
              //      //cout << "x = " << mch->getX() << ", y = " << mch->getY() << ", z = " << mch->getZ() << endl;
              //      TMatrixD mc4D(4,1,arr);
              //      mc4D = mInverseMFst[stateid-1000] * mc4D; // This is the MC position in the local frame where origins coincide between local and lab

              //      mcCoordsLab = TVector3(mch->getX(),mch->getY(),mch->getZ()); 

              //      if((stateid-1000)%3 == 0) mcCoords[0] = mc4D(0,0)-10.75; // we must shift the origin accoding to sensor
              //      if((stateid-1000)%3 != 0) mcCoords[0] = mc4D(0,0)-22.25;
              //      mcCoords[1] = mc4D(1,0);
              //      mcZ = mc4D(2,0);
              //    }
              //  }
              //  for (auto mch : static_cast<FwdHit*>(h)->_mcTrack->mFttHits)
              //  {
              //    if(mch->getZ() > 200.0)
              //    {
              //      int sector = mch->getSector();
              //      //cout << "mch sector = " << sector << endl;

              //      int quadId = static_cast<FwdHit*>(mch)->getSensor();
              //      //cout << "diskId = " << sector << ",    quadId = " << quadId << endl;
              //      int stgId = sector*4 + quadId;
              //      
              //      if(stgId != stateid) continue;
              //       
              //      mcCoordsLab = TVector3(mch->getX(),mch->getY(),mch->getZ());

              //      double xshift = 0.0;
              //      double yshift = 0.0;
              //      switch(quadId) {
              //        case 0 : // +x,+y
              //          xshift -= mCenterPent;
              //          yshift -= (5.9+mCenterPent);
              //          break;
              //        case 1 : // -x,+y
              //          xshift -= -mCenterPent;
              //          yshift -= (5.9+mCenterPent);
              //          break;
              //        case 2 : // -x,-y
              //          xshift -= (-6.5-mCenterPent);
              //          yshift -= (5.9-mCenterPent);
              //          break;
              //        case 3 : // +x,-y
              //          xshift -= (6.5+mCenterPent);
              //          yshift -= (5.9-mCenterPent);
              //          break;
              //      }
              //      mcCoords[0] = mch->getX()+xshift;
              //      mcCoords[1] = mch->getY()+yshift;
              //      mcZ = mch->getZ();
              //    }
              //  }
              //}
              //
              //if(mMcTracking)
              //{
              //  float arr[15];
              //  int idx = 0;
              //  arr[idx++] = float(stateid);
              //  //cout << "stateid = " << stateid << endl;
              //  arr[idx++] = float(mcCoords[0]);
              //  //cout << "mcx = " << mcCoords[0] << endl;
              //  arr[idx++] = float(mcCoords[1]);
              //  //cout << "mcy = " << mcCoords[1] << endl;
              //  arr[idx++] = float(mcZ);
              //  //cout << "mcz = " << mcZ << endl;
              //  arr[idx++] = float(rcCoords[0]);
              //  //cout << "rcx = " << rcCoords[0] << endl;
              //  arr[idx++] = float(rcCoords[1]);
              //  //cout << "rcy = " << rcCoords[1] << endl;
              //  arr[idx++] = float(rcCoordsLab[2]);
              //  //cout << "rcz = " << rcCoordsLab[2] << endl;
              //  arr[idx++] = float(trackCoords[0]);
              //  //cout << "trackx = " << trackCoords[0] << endl;
              //  arr[idx++] = float(trackCoords[1]);
              //  //cout << "tracky = " << trackCoords[1] << endl;
              //  arr[idx++] = float(mcCoordsLab[0]);
              //  //cout << "mcxlab = " << mcCoordsLab[0] << endl;
              //  arr[idx++] = float(mcCoordsLab[1]);
              //  //cout << "mcylab = " << mcCoordsLab[1] << endl;
              //  arr[idx++] = float(rcCoordsLab[0]);
              //  //cout << "rcxlab = " << rcCoordsLab[0] << endl;
              //  arr[idx++] = float(rcCoordsLab[1]);
              //  //cout << "rcylab = " << rcCoordsLab[1] << endl;
              //  arr[idx++] = float(trackCoordsLab[0]);
              //  //cout << "trackxlab = " << trackCoordsLab[0] << endl;
              //  arr[idx++] = float(trackCoordsLab[1]);
              //  //cout << "trackylab = " << trackCoordsLab[1] << endl;


              //  //cout << endl << "LOCAL RESIDUALS" << endl;
              //  //cout << "MC-RC X = " << mcCoords[0]-rcCoords[0] << endl;
              //  //cout << "MC-RC Y = " << mcCoords[1]-rcCoords[1] << endl;
              //  //cout << "Track-RC X = " << trackCoords[0]-rcCoords[0] << endl;
              //  //cout << "Track-RC Y = " << trackCoords[1]-rcCoords[1] << endl;
              //  //cout << "Track-MC X = " << trackCoords[0]-mcCoords[0] << endl;
              //  //cout << "Track-MC Y = " << trackCoords[1]-mcCoords[1] << endl << endl;
              //  //cout << "LAB RESIDUALS" << endl;
              //  //cout << "MC-RC X = " << mcCoordsLab[0]-rcCoordsLab[0] << endl;
              //  //cout << "MC-RC Y = " << mcCoordsLab[1]-rcCoordsLab[1] << endl;
              //  //cout << "Track-RC X = " << trackCoordsLab[0]-rcCoordsLab[0] << endl;
              //  //cout << "Track-RC Y = " << trackCoordsLab[1]-rcCoordsLab[1] << endl;
              //  //cout << "Track-MC X = " << trackCoordsLab[0]-mcCoordsLab[0] << endl;
              //  //cout << "Track-MC Y = " << trackCoordsLab[1]-mcCoordsLab[1] << endl << endl;

              //  mAlignmentInfo->Fill(arr);
              //}
              if(!mMcTracking)
              {
                float arr[9];
                int idx = 0;
                arr[idx++] = float(stateid);
                cout << "stateid = " << stateid << endl;
                arr[idx++] = float(rcCoords[0]);
                cout << "rcx = " << rcCoords[0] << endl;
                arr[idx++] = float(rcCoords[1]);
                cout << "rcy = " << rcCoords[1] << endl;
                arr[idx++] = float(trackCoords[0]);
                cout << "trackx = " << trackCoords[0] << endl;
                arr[idx++] = float(trackCoords[1]);
                cout << "tracky = " << trackCoords[1] << endl;
                arr[idx++] = float(rcCoordsLab[0]);
                cout << "rcxlab = " << rcCoordsLab[0] << endl;
                arr[idx++] = float(rcCoordsLab[1]);
                cout << "rcylab = " << rcCoordsLab[1] << endl;
                arr[idx++] = float(trackCoordsLab[0]);
                cout << "trackxlab = " << trackCoordsLab[0] << endl;
                arr[idx++] = float(trackCoordsLab[1]);
                cout << "trackylab = " << trackCoordsLab[1] << endl;
                mAlignmentInfo->Fill(arr);
              }
            }
          }
        }
        long long duration = (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6; // milliseconds
        return duration;
    } // fitTrack

    //TVector3 fitTrack(Seed_t trackSeed, double *Vertex = 0, TVector3 *seedMomentum = 0) {
    //    long long itStart = FwdTrackerUtils::nowNanoSecond();
    //    if (mGenHistograms) this->mHist["FitStatus"]->Fill("Total", 1);

    //    // The PV information, if we want to use it
    //    TVectorD pv(3);

    //    StarRandom rand = StarRandom::Instance();
    //    LOG_DEBUG << "Setting up the vertex info" << endm;
    //    if (0 == Vertex) { // randomized from simulation
    //        pv[0] = mVertexPos[0] + rand.gauss(mVertexSigmaXY);
    //        pv[1] = mVertexPos[1] + rand.gauss(mVertexSigmaXY);
    //        pv[2] = mVertexPos[2] + rand.gauss(mVertexSigmaZ);
    //    } else {
    //        pv[0] = Vertex[0];
    //        pv[1] = Vertex[1];
    //        pv[2] = Vertex[2];
    //    }

    //    // get the seed info from our hits
    //    TVector3 seedMom, seedPos;
    //    LOG_DEBUG << "Getting seed state" << endm;
    //    // returns track curvature if needed
    //    seedState(trackSeed, seedPos, seedMom);

    //    if (seedMomentum != nullptr) {
    //        seedMom = *seedMomentum;
    //    }

    //    // If we use the PV, use that as the start pos for the track
    //    if (mIncludeVertexInFit) {
    //        LOG_DEBUG << "Primary Vertex in fit (seed pos) @ " << TString::Format( "(%f, %f, %f)", pv[0], pv[1], pv[2] ).Data()  << endm;
    //        seedPos.SetXYZ(pv[0], pv[1], pv[2]);
    //    }

    //    if (mFitTrack){
    //        delete mFitTrack;
    //    }

    //    // create the track representations
    //    auto trackRepPos = new genfit::RKTrackRep(mPdgMuon);
    //    auto trackRepNeg = new genfit::RKTrackRep(mPdgAntiMuon);
    //    
    //    // Create the track
    //    mFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
    //    mFitTrack->addTrackRep(trackRepNeg);

    //    // TODO: TVector3 can fault on Eta() if Pt=0... Find a better fallback in this case for the seed 
    //    if ( fabs(seedMom.Z() / seedMom.Y()) > 1e10 ){
    //        seedMom.SetXYZ( 0.1, 0.1, -1 );
    //    }
    //    LOG_DEBUG << "seedPos : (" << seedPos.X() << ", " << seedPos.Y() << ", " << seedPos.Z() << " )" << endm;
    //    LOG_DEBUG << ", seedMom : (" << seedMom.X() << ", " << seedMom.Y() << ", " << seedMom.Z() << " )" << endm;
    //    LOG_DEBUG << ", seedMom : (" << seedMom.Pt() << ", " << seedMom.Eta() << ", " << seedMom.Phi() << " )" << endm;

    //    genfit::Track &fitTrack = *mFitTrack;

    //    size_t planeId(0);     // detector plane ID
    //    int hitId(0);       // hit ID

    //    // initialize the hit coords on plane
    //    TVectorD hitCoords(2);
    //    hitCoords[0] = 0;
    //    hitCoords[1] = 0;

    //    /******************************************************************************************************************
    //    * Include the Primary vertex if desired
    //    ******************************************************************************************************************/
    //    if (mIncludeVertexInFit) {

    //        TMatrixDSym hitCov3(3);
    //        hitCov3(0, 0) = mVertexSigmaXY * mVertexSigmaXY;
    //        hitCov3(1, 1) = mVertexSigmaXY * mVertexSigmaXY;
    //        hitCov3(2, 2) = mVertexSigmaZ * mVertexSigmaZ;

    //        genfit::SpacepointMeasurement *measurement = new genfit::SpacepointMeasurement(pv, hitCov3, 0, ++hitId, nullptr);
    //        fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
    //    }

    //    /******************************************************************************************************************
    //    	 * loop over the hits, add them to the track
    //    	 ******************************************************************************************************************/
    //    for (auto h : trackSeed) {
    //        
    //        const bool isFTT = h->getZ() > 200;
    //        hitCoords[0] = h->getX();
    //        hitCoords[1] = h->getY();
    //        
    //        genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitCoords, CovMatPlane(h), h->getSector(), ++hitId, nullptr);

    //        planeId = h->getSector();

    //        genfit::SharedPlanePtr plane;
    //        if (isFTT && mFTTPlanes.size() <= planeId) {
    //            LOG_ERROR << "invalid VolumId -> out of bounds DetPlane, vid = " << planeId << endm;
    //            return TVector3(0, 0, 0);
    //        }

    //        if (isFTT)
    //            plane = mFTTPlanes[planeId];
    //        else 
    //            plane = getFstPlane( static_cast<FwdHit*>(h) );

    //        measurement->setPlane(plane, planeId);
    //        fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

    //        if (abs(h->getZ() - plane->getO().Z()) > 0.05) {
    //            LOG_WARN << "Z Mismatch h->z = " << h->getZ() << ", plane->z = "<< plane->getO().Z() <<", diff = " << abs(h->getZ() - plane->getO().Z()) << endm;
    //        }
    //    } // loop on trackSeed


    //    /******************************************************************************************************************
    //    	 * Do the fit
    //    	 ******************************************************************************************************************/
    //    try {
    //        // do the fit
    //        mFitter->processTrackWithRep(&fitTrack, trackRepPos);
    //        mFitter->processTrackWithRep(&fitTrack, trackRepNeg);

    //    } catch (genfit::Exception &e) {
    //        if (mGenHistograms) mHist["FitStatus"]->Fill("Exception", 1);
    //    }

    //    TVector3 p(0, 0, 0);

    //    /******************************************************************************************************************
    //    	 * Now check the fit
    //    	 ******************************************************************************************************************/
    //    try {
    //        //check
    //        fitTrack.checkConsistency();

    //        // find track rep with smallest chi2
    //        fitTrack.determineCardinalRep();
    //        auto cardinalRep = fitTrack.getCardinalRep();
    //        auto cardinalStatus = fitTrack.getFitStatus(cardinalRep);
    //        mFitStatus = *cardinalStatus; // save the status of last fit

    //        // Delete any previous track rep
    //        if (mTrackRep)
    //            delete mTrackRep;

    //        // Clone the cardinal rep for persistency
    //        mTrackRep = cardinalRep->clone(); // save the result of the fit
    //        if (fitTrack.getFitStatus(cardinalRep)->isFitConverged() && mGenHistograms ) {
    //            this->mHist["FitStatus"]->Fill("GoodCardinal", 1);
    //        }

    //        if (fitTrack.getFitStatus(trackRepPos)->isFitConverged() == false &&
    //            fitTrack.getFitStatus(trackRepNeg)->isFitConverged() == false) {
    //        
    //            LOG_WARN << "FWD Track GenFit Failed" << endm;

    //            p.SetXYZ(0, 0, 0);
    //            long long duration = (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6; // milliseconds
    //            if (mGenHistograms) {
    //                this->mHist["FitStatus"]->Fill("Fail", 1);
    //                this->mHist["FailedFitDuration"]->Fill(duration);
    //            }
    //            return p;
    //        } // neither track rep converged

    //        p = cardinalRep->getMom(fitTrack.getFittedState(1, cardinalRep));
    //        mQ = cardinalRep->getCharge(fitTrack.getFittedState(1, cardinalRep));
    //        mP = p;

    //        LOG_DEBUG << "track fit p = " << TString::Format( "(%f, %f, %f), q=%f", p.X(), p.Y(), p.Z(), mQ ).Data() << endm;

    //    } catch (genfit::Exception &e) {
    //        LOG_WARN << "Exception on track fit: " << e.what() << endm;
    //        p.SetXYZ(0, 0, 0);

    //        long long duration = (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6; // milliseconds
    //        if (mGenHistograms) {
    //            this->mHist["FitStatus"]->Fill("Exception", 1);
    //            this->mHist["FailedFitDuration"]->Fill(duration);
    //        }

    //        return p;
    //    } // try/catch 

    //    long long duration = (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6; // milliseconds
    //    if (mGenHistograms) {
    //        this->mHist["FitStatus"]->Fill("Pass", 1);
    //        this->mHist["delta_fit_seed_pT"]->Fill(p.Pt() - seedMom.Pt());
    //        this->mHist["delta_fit_seed_eta"]->Fill(p.Eta() - seedMom.Eta());
    //        this->mHist["delta_fit_seed_phi"]->Fill(p.Phi() - seedMom.Phi());
    //        this->mHist["FitDuration"]->Fill(duration);
    //    }
    //    return p;
    //}

    int getCharge() {
        return (int)mQ;
    }

    // Store the planes for FTT and FST
    vector<genfit::SharedPlanePtr> mFTTPlanes;
    vector<genfit::SharedPlanePtr> mFSTPlanes;
    vector<genfit::SharedPlanePtr> mFSTSensorPlanes;
    vector<genfit::SharedPlanePtr> mFSTPlanesInner;
    vector<genfit::SharedPlanePtr> mFSTPlanesOuter;

    void SetIncludeVertex( bool vert ) { mIncludeVertexInFit = vert; }

  protected:
    std::unique_ptr<genfit::AbsBField> mBField;

    FwdTrackerConfig mConfig; // main config object
    TString mGeoCache;

    // optional histograms, off by default
    std::map<std::string, TH1 *> mHist;
    bool mGenHistograms = false;

    TNtuple *mAlignmentInfo;
    bool mMisaligned = false;

    int mModuleMap[3][12] = {{1,6,0,11,5,10,4,9,3,8,2,7},
                             {6,0,11,5,10,4,9,3,8,2,7,1},
                             {1,6,0,11,5,10,4,9,3,8,2,7}};

    int mSensorMap[3] = {2,0,1};

    TMatrixD mInverseMFst[108];
    TMatrixD mInverseMStgc[16];
    
    // Main GenFit fitter instance
    std::unique_ptr<genfit::AbsKalmanFitter> mFitter = nullptr;
    std::unique_ptr<StFwdGbl> mGblFitter = nullptr;

    // PDG codes for the default plc type for fits
    const int mPdgPiPlus = 211;
    const int mPdgPiMinus = -211;
    const int mPdgPositron = 11;
    const int mPdgElectron = -11;
    const int mPdgMuon = 13;
    const int mPdgAntiMuon = -13;


    // det z locations loaded from geom or config
    vector<double> mFSTZLocations, mFTTZLocations;

    // parameter ALIASED from mConfig wrt PV vertex
    double mVertexSigmaXY = 1;
    double mVertexSigmaZ = 30;
    vector<double> mVertexPos;
    bool mIncludeVertexInFit = false;
    bool mSmearMcVertex = false;
    bool mRefitGBL = false;
    double mCenterPent = 60.2361/2.0;
    double fstDefaultZ[12] = {150.008101+0.0009308,151.403100+0.000946,153.491899+0.0010220,152.096900+0.0010770,
                              166.989901+0.0010680,165.594901+0.001053,163.506101+0.0009766,164.901101+0.0009918,
                              177.039106+0.0009308,178.434106+0.000946,180.522905+0.0010220,179.127906+0.0010070};
    double stgcDefaultZ[4] = {281.082-0.01508, 304.062+0.004874, 325.028+0.03889, 348.068-0.001104};

    // GenFit state
    genfit::FitStatus mFitStatus;
    genfit::AbsTrackRep *mTrackRep;

    // Fit results
    TVector3 mP;
    double mQ;

    TFile *mAlignmentOutput = nullptr;

    bool mMcTracking = false;

    float mRSize = 3.0;
    float mPhiSize = 0.004;
    // GenFit state - resused 
    //
    std::shared_ptr<genfit::Track> mFitTrack;
};

#endif
