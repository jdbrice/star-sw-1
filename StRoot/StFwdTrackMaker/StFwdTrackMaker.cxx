#include "StFwdTrackMaker/StFwdTrackMaker.h"
#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StFwdTrackMaker/include/Tracker/FwdTracker.h"
#include "StFwdTrackMaker/include/Tracker/TrackFitter.h"
#include "StFwdTrackMaker/include/Tracker/FwdGeomUtils.h"
#include "StFwdTrackMaker/include/Tracker/ObjExporter.h"
#include "StFwdTrackMaker/StFwdGbl.h"

#include "KiTrack/IHit.h"
#include "GenFit/Track.h"
#include "GenFit/GFRaveVertexFactory.h"

#include "TMath.h"

#include <limits>
#include <map>
#include <string>
#include <string>
#include <vector>

#include "StBFChain/StBFChain.h"

#include "StEvent/StEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StHelixModel.h"
#include "StEvent/StPrimaryTrack.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StTrack.h"
#include "StEvent/StTrackGeometry.h"
#include "StEvent/StTrackNode.h"
#include "StEvent/StPrimaryVertex.h"
#include "StEvent/StEnumerations.h"
#include "StEvent/StTrackDetectorInfo.h"
#include "StEvent/StFttPoint.h"
#include "StEvent/StFcsHit.h"
#include "StEvent/StFcsCluster.h"
#include "StEvent/StFttCollection.h"
#include "StEvent/StFcsCollection.h"
#include "StEvent/StTriggerData.h"
#include "StEvent/StFstHitCollection.h"
#include "StEvent/StFstHit.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StChain/StChainOpt.h"

#include "StEventUtilities/StEventHelper.h"

#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_g2t_event_Table.h"

#include "StarMagField/StarMagField.h"

#include "St_base/StMessMgr.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"

#include <SystemOfUnits.h>
#include <exception>

#include "TROOT.h"
#include "TLorentzVector.h"

#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "StFstUtil/StFstCollection.h"

#include "StEvent/StFwdTrack.h"
#include "GenFit/AbsMeasurement.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFstCollection.h"
#include "StMuDSTMaker/COMMON/StMuFstHit.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"

#include "St_db_Maker/St_db_Maker.h"
#include "tables/St_Survey_Table.h"

#include "sys/types.h"
#include "sys/sysinfo.h"

FwdSystem* FwdSystem::sInstance = nullptr;
TMVA::Reader * BDTCrit2::reader = nullptr;
float BDTCrit2::Crit2_RZRatio = -999;
float BDTCrit2::Crit2_DeltaRho = -999;
float BDTCrit2::Crit2_DeltaPhi = -999;
float BDTCrit2::Crit2_StraightTrackRatio = -999;



#include "stdlib.h"
#include "stdio.h"
#include "string.h"

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

double getValueRAM(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    const double megabyte = 1024 * 1024;
    double rMb = ((double)result) / megabyte;
    return rMb;
}

//void reportMem(){
//    struct sysinfo memInfo;
//
//    sysinfo (&memInfo);
//    long long totalVirtualMem = memInfo.totalram;
//    //Add other values in next statement to avoid int overflow on right hand side...
//    totalVirtualMem += memInfo.totalswap;
//    totalVirtualMem *= memInfo.mem_unit;
//
//    long long virtualMemUsed = memInfo.totalram - memInfo.freeram;
//    //Add other values in next statement to avoid int overflow on right hand side...
//    virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
//    virtualMemUsed *= memInfo.mem_unit;
//
//    long long totalPhysMem = memInfo.totalram;
//    //Multiply in next statement to avoid int overflow on right hand side...
//    totalPhysMem *= memInfo.mem_unit;
//    long long physMemUsed = memInfo.totalram - memInfo.freeram;
//    //Multiply in next statement to avoid int overflow on right hand side...
//    physMemUsed *= memInfo.mem_unit;
//    const double megabyte = 1024 * 1024;
//
//    LOG_INFO << "VIRTUAL MEM USED % = " << ( (double)virtualMemUsed / (double)totalVirtualMem ) << endm;
//    LOG_INFO << "RAM USED % = " << ((double)physMemUsed / (double)totalPhysMem) << endm;
//    // memf << ( (double)virtualMemUsed / (double)totalVirtualMem ) << endl;
//    // memf << ( (double)virtualMemUsed ) << endl;
//    memf << (getValueRAM()) << endl;
//
//}

void FwdTreeData::clear(){
    // FTT
    fttN = 0;
    fttX.clear(); 
    fttY.clear(); 
    fttZ.clear();
    fttVolumeId.clear();
    fttPt.clear();
    fttTrackId.clear();
    fttVertexId.clear();

    // FST Info
    fstN = 0;
    fstX.clear(); 
    fstY.clear(); 
    fstZ.clear();
    fstTrackId.clear();
    
    // FCS info
    fcsN = 0;
    fcsX.clear(); 
    fcsY.clear(); 
    fcsZ.clear();
    fcsE.clear();
    fcsDet.clear();

    // EPD Info
    epdX.clear(); 
    epdY.clear(); 
    epdZ.clear();
    epdE.clear();

    // Reco Track info
    rcN = 0;
    rcPt.clear();
    rcEta.clear();
    rcPhi.clear();
    rcQuality.clear();
    rcTrackId.clear();
    rcNumFST.clear();
    rcCharge.clear();
    rcNumFTT.clear();
    rcNumPV.clear();
    rcEpdX.clear();
    rcEpdY.clear();

    mcN = 0;
    mcPt.clear();
    mcEta.clear();
    mcPhi.clear();
    mcVertexId.clear();
    mcCharge.clear();
    mcNumFtt.clear();
    mcNumFst.clear();

    // MonteCarlo Vertex info
    vmcN = 0;
    vmcX.clear();
    vmcY.clear();
    vmcZ.clear();

    // Reco vertex info
    vrcN = 0;
    vrcX.clear();
    vrcY.clear();
    vrcZ.clear();

    thdN = 0;
    thdX.clear();
    thdY.clear();
    thdP.clear();
    thdR.clear();
    thaX.clear();
    thaY.clear();
    thaZ.clear();


    Crits.clear();
    CritTrackIds.clear();
}

//_______________________________________________________________________________________
class GenfitUtils{
    public:

    // For now, accept anything we are passed, no matter what it is or how bad it is
    template<typename T> static bool accept( T ) { return true; }

    
}; // GenfitUtils

// Basic sanity cuts on genfit tracks
template<> bool GenfitUtils::accept( genfit::Track *track )
{
    // This also gets rid of failed fits (but may need to explicitly
    // for fit failure...)
    if (track->getNumPoints() <= 0 ) return false; // fit may have failed

    auto cardinal = track->getCardinalRep();

    // Check that the track fit converged
    auto status = track->getFitStatus( cardinal );
    if ( !status->isFitConverged() ) {
    return false;
    }


    // Next, check that all points on the track have fitter info
    // (may be another indication of a failed fit?)
    for ( auto point : track->getPoints() ) {
    if ( !point->hasFitterInfo(cardinal) ) {
    return false;
    }
    }

    // Following line fails with an exception, because some tracks lack 
    //   forward update, or prediction in fitter info at the first point
    //
    // genfit::KalmanFitterInfo::getFittedState(bool) const of 
    //                         GenFit/fitters/src/KalmanFitterInfo.cc:250

    // Fitted state at the first point
    // const auto &atFirstPoint = track->getFittedState();

    // Getting the fitted state from a track occasionally fails, because
    // the first point on the fit doesn't have forward/backward fit
    // information.  So we want the first point with fit info...

    genfit::TrackPoint* first = nullptr;
    unsigned int ipoint = 0;
    for ( ipoint = 0; ipoint < track->getNumPoints(); ipoint++ ) {
    first = track->getPointWithFitterInfo( ipoint );
    if ( first ) break;
    } 

    // No points on the track have fit information
    if ( !first ) {
        LOG_WARN << "No fit information on fwd genfit track" << endm;
        return false;
    }

    auto& fittedState= track->getFittedState(ipoint);

    TVector3 momentum = fittedState.getMom();
    double   pt       = momentum.Perp();

    if (pt < 0.10 ) return false; // below this

    return true;

};

//______________________________________________________________________________________
class SiRasterizer {
  public:
    SiRasterizer() {}
    SiRasterizer(FwdTrackerConfig &_cfg) { setup(_cfg); }
    ~SiRasterizer() {}
    void setup(FwdTrackerConfig &_cfg) {
        cfg = _cfg;
        mRasterR = cfg.get<double>("SiRasterizer:r", 3.0);
        mRasterPhi = cfg.get<double>("SiRasterizer:phi", 0.004);
    }

    bool active() {
        return cfg.get<bool>("SiRasterizer:active", false);
    }

    TVector3 raster(TVector3 p0) {
        TVector3 p = p0;
        double r = p.Perp();
        double phi = p.Phi();
        const double minR = 5.0;
        // 5.0 is the r minimum of the Si
        p.SetPerp(minR + (std::floor((r - minR) / mRasterR) * mRasterR + mRasterR / 2.0));
        p.SetPhi(-TMath::Pi() + (std::floor((phi + TMath::Pi()) / mRasterPhi) * mRasterPhi + mRasterPhi / 2.0));
        return p;
    }

    FwdTrackerConfig cfg;
    double mRasterR, mRasterPhi;
};

//  Wrapper class around the forward tracker
class ForwardTracker : public ForwardTrackMaker {
  public:
    // Replaces original initialization.  Config file and hitloader
    // will be provided by the maker.
    void initialize( TString geoCache, bool genHistograms ) {
        nEvents = 1; // only process single event

        // Create the forward system...
        FwdSystem::sInstance = new FwdSystem();

        // make our quality plotter
        mQualityPlotter = new QualityPlotter(mConfig);
        if ( genHistograms )
            mQualityPlotter->makeHistograms(mConfig.get<size_t>("TrackFinder:nIterations", 1));

        // initialize the track fitter
        mTrackFitter = new TrackFitter(mConfig, geoCache);
        mTrackFitter->setGenerateHistograms(genHistograms);
        mTrackFitter->setup();

        ForwardTrackMaker::initialize( geoCache, genHistograms );
    }

    void finish() {

        if ( mGenHistograms ){
            mQualityPlotter->finish();
            writeEventHistograms();
        }

        if (FwdSystem::sInstance){
            delete FwdSystem::sInstance;
            FwdSystem::sInstance = 0;
        }
        if (mQualityPlotter){
            delete mQualityPlotter;
            mQualityPlotter = 0;
        }
        if (mTrackFitter){
            mTrackFitter->finish();
            mTrackFitter->writeAlignment();
            delete mTrackFitter;
            mTrackFitter= 0;
        }
    }
};

//________________________________________________________________________
StFwdTrackMaker::StFwdTrackMaker() : StMaker("fwdTrack"), mGenHistograms(false), mGenTree(false), mForwardTracker(nullptr), mForwardData(nullptr), mGeoCache(""){
    SetAttr("useFtt",1);                 // Default Ftt on 
    SetAttr("useFst",1);                 // Default Fst on
    SetAttr("useFcs",1);                 // Default Fcs on
    SetAttr("config", "");     // Default configuration file (user may override before Init())
    SetAttr("fillEvent",1); // fill StEvent
};

int StFwdTrackMaker::Finish() {
    auto prevDir = gDirectory;
    if ( mGenHistograms ) {
        
        // output file name
        string name = mFwdConfig.get<string>("Output:url", "fwdTrackerOutput.root");
        LOG_INFO << "Saving StFwdTrackMaker Histograms to ROOT file: " << name << endm;
        TFile *fOutput = new TFile(name.c_str(), "RECREATE");
        fOutput->cd();

        fOutput->mkdir("StFwdTrackMaker");
        fOutput->cd("StFwdTrackMaker");
        for (auto nh : mHistograms) {
            nh.second->SetDirectory(gDirectory);
            nh.second->Write();
        }
        fOutput->cd("");
    }

    mForwardTracker->finish();

    prevDir->cd();

    if (mGenTree) {
        mTreeFile->cd();
        mTree->Write();
        mTreeFile->Write();
    }
    return kStOk;
}

void StFwdTrackMaker::LoadConfiguration() {
    if (mConfigFile.length() < 5){    
        LOG_INFO << "Forward Tracker is using default config for ";
        if ( defaultConfig == defaultConfigData ){
            LOG_INFO << " DATA" << endm;
        } else {
            LOG_INFO << " Simulation" << endm;
        }
        mFwdConfig.load( defaultConfig, true );
    } else {
        LOG_INFO << "Forward Tracker is using config from file : " <<  mConfigFile << endm;
        mFwdConfig.load( mConfigFile );
    }
    configLoaded = true;
}

//________________________________________________________________________
int StFwdTrackMaker::Init() {
    if ( !configLoaded ){
        LoadConfiguration();
    }

    if (mGenTree) {
        mTreeFile = new TFile("fwdtree.root", "RECREATE");
        mTree = new TTree("fwd", "fwd tracking tree");
        mTree->Branch("fttN",         &mTreeData. fttN, "fttN/I");
        mTree->Branch("fttX",         &mTreeData. fttX  );
        mTree->Branch("fttY",         &mTreeData. fttY  );
        mTree->Branch("fttZ",         &mTreeData. fttZ  );
        
        mTree->Branch("fttTrackId",   &mTreeData. fttTrackId   );
        mTree->Branch("fttVolumeId",  &mTreeData. fttVolumeId  );
        mTree->Branch("fttPt",        &mTreeData. fttPt        );
        mTree->Branch("fttVertexId",  &mTreeData. fttVertexId  );

        mTree->Branch("fstN",         &mTreeData. fstN, "fstN/I");
        mTree->Branch("fstX",         &mTreeData. fstX  );
        mTree->Branch("fstY",         &mTreeData. fstY  );
        mTree->Branch("fstZ",         &mTreeData. fstZ  );
        mTree->Branch("fstTrackId",   &mTreeData. fstTrackId  );

        mTree->Branch("fcsN",         &mTreeData. fcsN, "fcsN/I");
        mTree->Branch("fcsX",         &mTreeData. fcsX  );
        mTree->Branch("fcsY",         &mTreeData. fcsY  );
        mTree->Branch("fcsZ",         &mTreeData. fcsZ  );
        mTree->Branch("fcsE",         &mTreeData. fcsE  );
        mTree->Branch("fcsDet",       &mTreeData. fcsDet  );

        mTree->Branch("epdX",         &mTreeData. epdX  );
        mTree->Branch("epdY",         &mTreeData. epdY  );
        mTree->Branch("epdZ",         &mTreeData. epdZ  );
        mTree->Branch("epdE",         &mTreeData. epdE  );

        // mc tracks
        mTree->Branch("mcN",        &mTreeData. mcN, "mcN/I");
        mTree->Branch("mcPt",       &mTreeData. mcPt );
        mTree->Branch("mcEta",      &mTreeData. mcEta );
        mTree->Branch("mcPhi",      &mTreeData. mcPhi );
        mTree->Branch("mcCharge",   &mTreeData. mcCharge );
        mTree->Branch("mcVertexId", &mTreeData. mcVertexId );
        mTree->Branch("mcNumFtt",   &mTreeData. mcNumFtt );
        mTree->Branch("mcNumFst",   &mTreeData. mcNumFst );

        // mcverts
        mTree->Branch("vmcN",       &mTreeData. vmcN, "vmcN/I");
        mTree->Branch("vmcX",       &mTreeData. vmcX );
        mTree->Branch("vmcY",       &mTreeData. vmcY );
        mTree->Branch("vmcZ",       &mTreeData. vmcZ );

        // rcverts
        mTree->Branch("vrcN",       &mTreeData. vrcN, "vrcN/I");
        mTree->Branch("vrcX",       &mTreeData. vrcX );
        mTree->Branch("vrcY",       &mTreeData. vrcY );
        mTree->Branch("vrcZ",       &mTreeData. vrcZ );

        // rc tracks
        mTree->Branch("rcN",        &mTreeData. rcN, "rcN/I");
        mTree->Branch("rcPt",       &mTreeData. rcPt );
        mTree->Branch("rcEta",      &mTreeData. rcEta );
        mTree->Branch("rcPhi",      &mTreeData. rcPhi );
        mTree->Branch("rcCharge",   &mTreeData. rcCharge );
        mTree->Branch("rcTrackId",  &mTreeData. rcTrackId );
        mTree->Branch("rcNumFST",   &mTreeData. rcNumFST );
        mTree->Branch("rcNumFTT",   &mTreeData. rcNumFTT );
        mTree->Branch("rcNumPV",    &mTreeData. rcNumPV );
        mTree->Branch("rcQuality",  &mTreeData. rcQuality );

        mTree->Branch("rcEpdX",     &mTreeData. rcEpdX );
        mTree->Branch("rcEpdY",     &mTreeData. rcEpdY );

        mTree->Branch("thdN",         &mTreeData. thdN, "thdN/I");
        mTree->Branch("thdX",         &mTreeData. thdX  );
        mTree->Branch("thdY",         &mTreeData. thdY  );
        mTree->Branch("thaX",         &mTreeData. thaX  );
        mTree->Branch("thaY",         &mTreeData. thaY  );
        mTree->Branch("thaZ",         &mTreeData. thaZ  );
        mTree->Branch("thdP",         &mTreeData. thdP  );
        mTree->Branch("thdR",         &mTreeData. thdR  );

        // track projections
        mTree->Branch("tprojN",       &mTreeData. tprojN,  "tprojN/I");
        mTree->Branch("tprojIdD",     &mTreeData. tprojIdD);
        mTree->Branch("tprojIdT",     &mTreeData. tprojIdT);
        mTree->Branch("tprojX",       &mTreeData. tprojX);
        mTree->Branch("tprojY",       &mTreeData. tprojY);
        mTree->Branch("tprojZ",       &mTreeData. tprojZ);
        mTree->Branch("tprojPx",      &mTreeData. tprojPx);
        mTree->Branch("tprojPy",      &mTreeData. tprojPy);
        mTree->Branch("tprojPz",      &mTreeData. tprojPz);

        std::string path = "TrackFinder.Iteration[0].SegmentBuilder";
        std::vector<string> paths = mFwdConfig.childrenOf(path);

        if (mTreeData.saveCrit){
            for (string p : paths) {
                string name = mFwdConfig.get<string>(p + ":name", "");
                mTreeData.Crits[name]; // create the entry
                mTree->Branch(name.c_str(), &mTreeData.Crits[name]);
                mTree->Branch((name + "_trackIds").c_str(), &mTreeData.CritTrackIds[name]);

                if ( name == "Crit2_RZRatio" ){
                    string n = name + "_x1";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);

                    n = name + "_y1";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);

                    n = name + "_z1";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);

                    n = name + "_x2";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);

                    n = name + "_y2";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);

                    n = name + "_z2";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);

                    n = name + "_h1";
                    mTreeData.CritTrackIds[(n)]; mTree->Branch(n.c_str(), &mTreeData.CritTrackIds[n]);
                    n = name + "_h2";
                    mTreeData.CritTrackIds[(n)]; mTree->Branch(n.c_str(), &mTreeData.CritTrackIds[n]);
                    n = name + "_h3";
                    mTreeData.CritTrackIds[(n)]; mTree->Branch(n.c_str(), &mTreeData.CritTrackIds[n]);
                }

                if ( name == "Crit2_BDT" ){
                    string n = name + "_DeltaPhi";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);
                    n = name + "_DeltaRho";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);
                    n = name + "_RZRatio";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);
                    n = name + "_StraightTrackRatio";
                    mTreeData.Crits[(n)]; mTree->Branch(n.c_str(), &mTreeData.Crits[n]);
                }
            }

            // Three hit criteria
            path = "TrackFinder.Iteration[0].ThreeHitSegments";
            paths = mFwdConfig.childrenOf(path);

            for (string p : paths) {
                string name = mFwdConfig.get<string>(p + ":name", "");
                mTreeData.Crits[name]; // create the entry
                mTree->Branch(name.c_str(), &mTreeData.Crits[name]);
                mTree->Branch((name + "_trackIds").c_str(), &mTreeData.CritTrackIds[name]);
            }
        }

        mTree->SetAutoFlush(0);
    } // gen tree

    if ( mGeoCache == "" ){
        /// Instantiate and cache the geometry 
        GetDataBase("VmcGeometry");

        mGeoCache = GetChainOpt()->GetFileOut();  
        if ( mGeoCache=="" ) 
            mGeoCache = GetChainOpt()->GetFileIn();

        // Strip out @ symbol
        mGeoCache = mGeoCache.ReplaceAll("@",""); 
        // Strip off the last extention in the mGeoCache
        mGeoCache = mGeoCache( 0, mGeoCache.Last('.') );
        // Append geom.root to the extentionless mGeoCache
        mGeoCache+=".geom.root";
    }
    // create an SiRasterizer in case we need it 
    mSiRasterizer = std::shared_ptr<SiRasterizer>( new SiRasterizer(mFwdConfig));
    mForwardTracker = std::shared_ptr<ForwardTracker>(new ForwardTracker( ));
    mForwardTracker->setConfig(mFwdConfig);

    // only save criteria values if we are generating a tree.
    mForwardTracker->setSaveCriteriaValues(mGenTree);

    mForwardData = std::shared_ptr<FwdDataSource>(new FwdDataSource());
    mForwardTracker->setData(mForwardData);
    mForwardTracker->initialize( mGeoCache, mGenHistograms );

    // geometry should be available from here (mForwardTracker will initialize cache if needed)
    if (gGeoManager) {
        FwdGeomUtils fwdGeoUtils( gGeoManager );
        // get the z-locations from geometry model and fallback to the defaults
        auto fstZ = fwdGeoUtils.fstZ( {151.750000, 165.248001, 178.781006} );
        mFstZFromGeom.assign( fstZ.begin(), fstZ.end() );
        auto fttZ = fwdGeoUtils.fttZ( {280.904999, 303.704987, 326.605011, 349.404999} );
        mFttZFromGeom.assign( fttZ.begin(), fttZ.end() );
    }

    if ( mGenHistograms ){
        mHistograms["fwdVertexZ"] = new TH1D("fwdVertexZ", "FWD Vertex (RAVE);z", 1000, -50, 50);
        mHistograms["fwdVertexXY"] = new TH2D("fwdVertexXY", "FWD Vertex (RAVE);x;y", 100, -1, 1, 100, -1, 1);
        mHistograms["fwdVertexDeltaZ"] = new TH2D("fwdVertexDeltaZ", "FWD Vertex - MC Vertex;#Delta z", 100, -1, 1, 100, -1, 1);
        
        mHistograms["McEventEta"] = new TH1D("McEventEta", ";MC Track Eta", 1000, -5, 5);
        mHistograms["McEventPt"] = new TH1D("McEventPt", ";MC Track Pt (GeV/c)", 1000, 0, 10);
        mHistograms["McEventPhi"] = new TH1D("McEventPhi", ";MC Track Phi", 1000, 0, 6.2831852);

        // these are tracks within 2.5 < eta < 4.0
        mHistograms["McEventFwdEta"] = new TH1D("McEventFwdEta", ";MC Track Eta", 1000, -5, 5);
        mHistograms["McEventFwdPt"] = new TH1D("McEventFwdPt", ";MC Track Pt (GeV/c)", 1000, 0, 10);
        mHistograms["McEventFwdPhi"] = new TH1D("McEventFwdPhi", ";MC Track Phi", 1000, 0, 6.2831852);

        // create mHistograms
        mHistograms["nMcTracks"] = new TH1I("nMcTracks", ";# MC Tracks/Event", 1000, 0, 1000);
        mHistograms["nMcTracksFwd"] = new TH1I("nMcTracksFwd", ";# MC Tracks/Event", 1000, 0, 1000);
        mHistograms["nMcTracksFwdNoThreshold"] = new TH1I("nMcTracksFwdNoThreshold", ";# MC Tracks/Event", 1000, 0, 1000);

        mHistograms["nHitsSTGC"] = new TH1I("nHitsSTGC", ";# STGC Hits/Event", 1000, 0, 1000);
        mHistograms["nHitsFSI"] = new TH1I("nHitsFSI", ";# FSIT Hits/Event", 1000, 0, 1000);

        mHistograms["stgc_volume_id"] = new TH1I("stgc_volume_id", ";stgc_volume_id", 50, 0, 50);
        mHistograms["fsi_volume_id"] = new TH1I("fsi_volume_id", ";fsi_volume_id", 50, 0, 50);

        mHistograms["fsiHitDeltaR"] = new TH1F("fsiHitDeltaR", "FSI; delta r (cm); ", 500, -5, 5);
        mHistograms["fsiHitDeltaPhi"] = new TH1F("fsiHitDeltaPhi", "FSI; delta phi; ", 500, -5, 5);

        // there are 4 stgc stations
        for (int i = 0; i < 4; i++) {
            mHistograms[TString::Format("stgc%dHitMap", i).Data()] = new TH2F(TString::Format("stgc%dHitMap", i), TString::Format("STGC Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);

            mHistograms[TString::Format("stgc%dHitMapPrim", i).Data()] = new TH2F(TString::Format("stgc%dHitMapPrim", i), TString::Format("STGC Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);
            mHistograms[TString::Format("stgc%dHitMapSec", i).Data()] = new TH2F(TString::Format("stgc%dHitMapSec", i), TString::Format("STGC Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);
        }

        // There are 3 silicon stations
        for (int i = 0; i < 3; i++) {
            mHistograms[TString::Format("fsi%dHitMap", i).Data()] = new TH2F(TString::Format("fsi%dHitMap", i), TString::Format("FSI Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);
            mHistograms[TString::Format("fsi%dHitMapZ", i).Data()] = new TH2F(TString::Format("fsi%dHitMapZ", i), TString::Format("FSI Layer %d; x (cm); y(cm)", i), 200, -100, 100, 200, -100, 100);

            mHistograms[TString::Format("fsi%dHitMapR", i).Data()] = new TH1F(TString::Format("fsi%dHitMapR", i), TString::Format("FSI Layer %d; r (cm); ", i), 500, 0, 50);
            mHistograms[TString::Format("fsi%dHitMapPhi", i).Data()] = new TH1F(TString::Format("fsi%dHitMapPhi", i), TString::Format("FSI Layer %d; phi; ", i), 320, 0, TMath::Pi() * 2 + 0.1);
        }

    } // mGenHistograms

    //-----  LOAD ALIGNMENT MATRICES  -----//
    St_db_Maker *dbMk=new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
    dbMk->SetDebug();
    dbMk->SetDateTime(20211026,7); // event or run start time, set to your liking
    dbMk->SetFlavor("ofl");

    dbMk->Init();
    dbMk->Make();
    
    Survey_st *fstOnTpc;
    Survey_st *hssOnFst;
    Survey_st *fstWedgeOnHss;
    Survey_st *fstSensorOnWedge;

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


    // default FST sensor z locations, found externally from ideal geometry
    for (int is = 0; is < 108; is++) {
        
        // FST half
        int h = (is / 18) % 2; // 0 (left +x half), 1 (right -x half) 
    
        // FST wedge
        int w = is / 3; // 0-35

        if(h == 0)
        {
          // create matrices from alignment tables
          MhssOnFst(0,0) = 1.0; 
          MhssOnFst(0,1) = 0.0; 
          MhssOnFst(0,2) = 0.0;
          MhssOnFst(1,0) = 0.0; 
          MhssOnFst(1,1) = 1.0; 
          MhssOnFst(1,2) = 0.0;
          MhssOnFst(2,0) = 0.0; 
          MhssOnFst(2,1) = 0.0; 
          MhssOnFst(2,2) = 1.0;
          MhssOnFst(0,3) = 0.0;
          MhssOnFst(1,3) = 0.0;
          MhssOnFst(2,3) = 0.0;//hssOnFst[h].t2 ;
          MhssOnFst(3,3) = 1.0            ;
        }
        if(h == 1)
        {
          MhssOnFst(0,0) = 1.0;//hssOnFst[h].r00; 
          MhssOnFst(0,1) = 0.0;//hssOnFst[h].r01; 
          MhssOnFst(0,2) = 0.0;//hssOnFst[h].r02;
          MhssOnFst(1,0) = 0.0;//hssOnFst[h].r10; 
          MhssOnFst(1,1) = 1.0;//hssOnFst[h].r11; 
          MhssOnFst(1,2) = 0.0;//hssOnFst[h].r12;
          MhssOnFst(2,0) = 0.0;//hssOnFst[h].r20; 
          MhssOnFst(2,1) = 0.0;//hssOnFst[h].r21; 
          MhssOnFst(2,2) = 1.0;//hssOnFst[h].r22;
          MhssOnFst(0,3) = 0.0;
          MhssOnFst(1,3) = 0.0;
          MhssOnFst(2,3) = 0.0;
          MhssOnFst(3,3) = 1.0            ;
        }
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

        //if(is != mMisSensor)
        //{
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
        //}
        //else
        //{       
        //  double cost = TMath::Cos(mDeltaGamma);//0.002);    
        //  double sint = TMath::Sin(mDeltaGamma);//0.002);    
        //  double cos75 = TMath::Cos(1.309);    
        //  double sin75 = TMath::Sin(1.309);  
        //  double dxu =  cos75*mDeltaU; //0.005; // dv = 100 um   
        //  double dyu =  sin75*mDeltaU; //0.005; // dv = 100 um   
        //  double dxv = -sin75*mDeltaV; //0.005; // dv = 100 um   
        //  double dyv =  cos75*mDeltaV; //0.005; // dv = 100 um   
        //  double dx = dxu + dxv; 
        //  double dy = dyu + dyv; 
        //  MfstSensorOnWedge(0,0) = cost;//fstSensorOnWedge[is].r00; 
        //  MfstSensorOnWedge(0,1) = -sint;//fstSensorOnWedge[is].r01; 
        //  MfstSensorOnWedge(0,2) = fstSensorOnWedge[is].r02;
        //  MfstSensorOnWedge(1,0) = sint;//fstSensorOnWedge[is].r10; 
        //  MfstSensorOnWedge(1,1) = cost;//fstSensorOnWedge[is].r11; 
        //  MfstSensorOnWedge(1,2) = fstSensorOnWedge[is].r12;
        //  MfstSensorOnWedge(2,0) = fstSensorOnWedge[is].r20; 
        //  MfstSensorOnWedge(2,1) = fstSensorOnWedge[is].r21; 
        //  MfstSensorOnWedge(2,2) = fstSensorOnWedge[is].r22;
        //  MfstSensorOnWedge(0,3) = dx;//fstSensorOnWedge[is].t0 ;
        //  MfstSensorOnWedge(1,3) = dy;//fstSensorOnWedge[is].t1 ;
        //  MfstSensorOnWedge(2,3) = fstSensorOnWedge[is].t2 ;
        //  MfstSensorOnWedge(3,3) = 1.0                     ;
        //}

        // Rotate and Translate plane normal vectors and origin
        TMatrixD M = MfstOnTpc * MhssOnFst * MfstWedgeOnHss * MfstSensorOnWedge;
        M.Print();

        // save this inverse matrix for use with misaligned simulated data
        mInverseM[is].ResizeTo(4,4);
        mInverseM[is] = M.Invert();
        //mInverseM[is].Print();

        MhssOnFst.Zero();
        MfstWedgeOnHss.Zero();
        MfstSensorOnWedge.Zero();
    }

    

    LOG_DEBUG << "StFwdTrackMaker::Init" << endm;
    cout << "StFwdTrackMaker::Init" << endl;
    return kStOK;
};

TMatrixDSym makeSiCovMat(TVector3 hit, FwdTrackerConfig &xfg) {
    // we can calculate the CovMat since we know the det info, but in future we should probably keep this info in the hit itself

    float rSize = xfg.get<float>("SiRasterizer:r", 3.0);
    float phiSize = xfg.get<float>("SiRasterizer:phi", 0.004);

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

    TMatrixDSym tamvoc(3);
    tamvoc( 0, 0 ) = cm(0, 0); tamvoc( 0, 1 ) = cm(0, 1); tamvoc( 0, 2 ) = 0.0;
    tamvoc( 1, 0 ) = cm(1, 0); tamvoc( 1, 1 ) = cm(1, 1); tamvoc( 1, 2 ) = 0.0;
    tamvoc( 2, 0 ) = 0.0;      tamvoc( 2, 1 ) = 0.0; tamvoc( 2, 2 )      = 0.01*0.01;


    return tamvoc;
}

void StFwdTrackMaker::loadFttHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    LOG_INFO << "Loading FTT Hits" << endm;
    // Get the StEvent handle to see if the rndCollection is available
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    string fttFromSource = mFwdConfig.get<string>( "Source:ftt", "" );

    if (!event){ 
        LOG_ERROR << "No StEvent, cannot load Ftt Data" << endm;
        return;
    }

    // Load GEANT hits directly if requested
    if ( "GEANT" == fttFromSource  ) {
        LOG_INFO << "Loading sTGC hits directly from GEANT hits" << endm;
        loadFttHitsFromGEANT( mcTrackMap, hitMap, count );
        return;
    }
    
    if ( "FASTSIM" == fttFromSource ) {
        LOG_INFO << "Loading sTGC hits from StFttFastSim hits" << endm;
        loadFttHitsFromFastSim( mcTrackMap, hitMap, count );
        return;
    }

    StFttCollection *col = event->fttCollection();
    // From Data
    if ( col || "DATA" == fttFromSource ) {
        LOG_INFO << "Loading sTGC hits from StEvent hits" << endm;
        loadFttHitsFromStEvent( mcTrackMap, hitMap, count );
        return;
    }
} // loadFttHits

void StFwdTrackMaker::loadFttHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    LOG_INFO << "Loading FTT Hits from Data" << endm;
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    StFttCollection *col = event->fttCollection();
    
    mTreeData.fttN = 0;

    if ( col && col->numberOfPoints() > 0 ){
        LOG_INFO << "The Ftt Collection has " << col->numberOfPoints() << " points" << endm;
        TMatrixDSym hitCov3(3);
        const double sigXY = 0.2; // 
        hitCov3(0, 0) = sigXY * sigXY;
        hitCov3(1, 1) = sigXY * sigXY;
        hitCov3(2, 2) = 4; // unused since they are loaded as points on plane
        static const double mm_to_cm = 0.1;
        for ( auto point : col->points() ){
            
            mFwdHits.push_back(FwdHit(count++, point->xyz().x()*mm_to_cm, point->xyz().y()*mm_to_cm, point->xyz().z(), -point->plane(), 0, hitCov3, nullptr));
            auto hit = &(mFwdHits.back());
            cout << "Ftt hit->getLayer() = " << hit->getLayer() << endl;
            hit->setSensor(point->quadrant());
            mFttHits.push_back( TVector3( hit->getX(), hit->getY(), point->xyz().z() )  );
            if ( mGenHistograms ) {
                mHistograms[TString::Format("stgc%dHitMapSec", point->plane()).Data()]->Fill(hit->getX(), hit->getY());
            }
            // Add the hit to the hit map
            hitMap[hit->getSector()].push_back(hit);

            if (mGenTree && (unsigned)mTreeData.fttN < MAX_TREE_ELEMENTS) {
                LOG_DEBUG << "FttPoint( " << TString::Format( "[plane=%d, quad=%d, nClu=%d]", point->plane(), point->quadrant(), point->nClusters() ) << point->xyz().x()/10.0 << ", " << point->xyz().y()/10.0 << ", " << point->xyz().z() << " )" << endm;
                mTreeData.fttX.push_back( hit->getX() );
                mTreeData.fttY.push_back( hit->getY() );
                mTreeData.fttZ.push_back( point->xyz().z() );
                mTreeData.fttTrackId.push_back( 0 );
                mTreeData.fttVolumeId.push_back( point->plane() );
                mTreeData.fttPt.push_back( 0 );
                mTreeData.fttVertexId.push_back( 0 );
                mTreeData.fttN++;
            }
        }

        return;
    } else {
        LOG_INFO << "The Ftt Collection is EMPTY points" << endm;
    }
    LOG_INFO << "Number of FTT in TTree: " << mTreeData.fttN << endm;
}

void StFwdTrackMaker::loadFttHitsFromFastSim( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
  
    St_g2t_fts_hit *g2t_stg_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_stg_hit");

    //cout << "Loading FTT hits from GEANT" << endl;

    //this->mTreeData.fttN = 0;
    if (!g2t_stg_hits){
        //cout << "geant/g2t_stg_hit is empty" << endl; 
        LOG_WARN << "geant/g2t_stg_hit is empty" << endm;
        return;
    }

    // make the Covariance Matrix once and then reuse
    TMatrixDSym hitCov3(3);
    const double sigXY = 0.01;
    hitCov3(0, 0) = sigXY * sigXY;
    hitCov3(1, 1) = sigXY * sigXY;
    hitCov3(2, 2) = 0.001 * 0.001; // unused since they are loaded as points on plane

    int nstg = g2t_stg_hits->GetNRows();

    LOG_DEBUG << "This event has " << nstg << " stg hits in geant/g2t_stg_hit " << endm;
    //cout << "This event has " << nstg << " stg hits in geant/g2t_stg_hit " << endl;
    //if ( mGenHistograms ) {
    //    this->mHistograms["nHitsSTGC"]->Fill(nstg);
    //}


    bool filterGEANT = mFwdConfig.get<bool>( "Source:fttFilter", false );
    for (int i = 0; i < nstg; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_stg_hits->At(i);
        if (0 == git)
            continue; // geant hit
        int track_id = git->track_p;
        int volume_id = git->volume_id;
        //int plane_id = (volume_id - 1) / 4;           // from 1 - 16. four chambers per station
        // volume_id = (1 front | 2 back) + 10 * (quadrant 0-3) + 100 * (station 0-4)
        int plane_id = ((volume_id - 1) / 100) + 9 ;

        //LOG_DEBUG << "sTGC hit: volume_id = " << volume_id << " disk = " << plane_id << endm;
        cout << "sTGC hit: volume_id = " << volume_id << " disk = " << plane_id << endl;
        //if (disk < 9)
        //
        float x = git->x[0];// + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float y = git->x[1];// + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float z = git->x[2];

        int quadId = int((volume_id-1) % 100 / 10);
 
        cout << "plane_id = " << plane_id << ",    quadId = " << quadId << endl;

        //if(plane_id == 10 && quadId == 0) continue;

        cout << "sTGC HIT: x = " << git->x[0] << ", y = " << git->x[1] << ", z = " << git->x[2] << ",   track ID = " << track_id << endl;
        //cout << "global position: x = " << x << ", y = " << y << ", z = " << z << endl;      

        // Now that geometry has a front and back, we skip points on the back module for fast sim
        if (plane_id < 9 || volume_id % 2 == 0)
            continue;

        FwdHit *hit = new FwdHit(count++, x, y, z, -(plane_id-9), track_id, hitCov3, mcTrackMap[track_id]);
        hit->setSensor(quadId);
        // Add the hit to the hit map
        //hitMap[hit->getSector()].push_back(hit);
        //mFttHits.push_back( TVector3( x, y, z )  );

        // Add hit pointer to the track
        if (mcTrackMap[track_id])
            mcTrackMap[track_id]->addHit(hit);
    } // loop on hits

    LOG_INFO << "Loading FTT Hits from FastSim" << endm;
    StEvent *event = (StEvent *)this->GetDataSet("StEvent");

    StRnDHitCollection *rndCollection = event->rndHitCollection();
    const StSPtrVecRnDHit &hits = rndCollection->hits();

    if ( rndCollection && rndCollection->numberOfHits() > 0 ){
        LOG_INFO << "The RnD Collection has " << rndCollection->numberOfHits() << " points" << endm;
        //TMatrixDSym hitCov3(3);
        //const double sigXY = 0.01; // 
        //hitCov3(0, 0) = sigXY * sigXY;
        //hitCov3(1, 1) = sigXY * sigXY;
        //hitCov3(2, 2) = 0.001 * 0.001; // unused since they are loaded as points on plane
        for ( auto point : hits ){

            if(point->layer() < 9) continue;
            cout << "plane_id = " << point->layer() << ",    quadId = " << point->wafer() << endl;

            //if(point->layer() == 10 && point->wafer() == 0) continue;

            FwdHit *hit = new FwdHit(count++, point->position().x(), point->position().y(), point->position().z(), -(point->layer()-9), point->idTruth(), hitCov3, mcTrackMap[point->idTruth()]);
            hit->setSensor(point->wafer()); // quadrant 

            int plane_id = 0;
            mFttHits.push_back( TVector3( point->position().x(), point->position().y(), point->position().z())  );
            if ( mGenHistograms ) {
                this->mHistograms[TString::Format("stgc%dHitMapSec", point->layer()-9).Data()]->Fill(point->position().x(), point->position().y());
            }
            // Add the hit to the hit map
            hitMap[hit->getSector()].push_back(hit);

            if (mGenTree && mTreeData.fttN < MAX_TREE_ELEMENTS) {
                LOG_DEBUG << "Adding FTT Point to TTrees" << endm;
                LOG_DEBUG << "FttPoint( " << point->position().x() << ", " << point->position().y() << ", " << point->position().z() << " )" << endm;
                mTreeData.fttX.push_back( point->position().x() );
                mTreeData.fttY.push_back( point->position().y() );
                mTreeData.fttZ.push_back( point->position().z() );
                mTreeData.fttTrackId.push_back( 0 );
                mTreeData.fttVolumeId.push_back( point->layer() );
                mTreeData.fttPt.push_back( 0 );
                mTreeData.fttVertexId.push_back( 0 );
                mTreeData.fttN++;
            }
        }
        LOG_INFO << "Number of FTT hits in RnDHitCollection: " << mTreeData.fttN << endm;
        return;
    } else {
        LOG_INFO << "The Ftt Collection is EMPTY points" << endm;
    }
    LOG_INFO << "Number of FTT in TTree: " << mTreeData.fttN << endm;
}

void StFwdTrackMaker::loadFttHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    /************************************************************/
    // STGC Hits
    St_g2t_fts_hit *g2t_stg_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_stg_hit");


    mTreeData.fttN = 0;
    if (!g2t_stg_hits){
        LOG_WARN << "geant/g2t_stg_hit is empty" << endm; 
        return;
    }

    // make the Covariance Matrix once and then reuse
    TMatrixDSym hitCov3(3);
    const double sigXY = 0.01;
    hitCov3(0, 0) = sigXY * sigXY;
    hitCov3(1, 1) = sigXY * sigXY;
    hitCov3(2, 2) = 1.0; // unused since they are loaded as points on plane

    int nstg = g2t_stg_hits->GetNRows();

    LOG_DEBUG << "This event has " << nstg << " stg hits in geant/g2t_stg_hit " << endm;
    if ( mGenHistograms ) {
        mHistograms["nHitsSTGC"]->Fill(nstg);
    }
    

    bool filterGEANT = mFwdConfig.get<bool>( "Source:fttFilter", false );
    for (int i = 0; i < nstg; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_stg_hits->At(i);
        if (0 == git)
            continue; // geant hit
        int track_id = git->track_p;
        int volume_id = git->volume_id;
        int plane_id = (volume_id - 1) / 100;           // from 1 - 16. four chambers per station

        // only use the hits on the front modules
        if ( volume_id % 2 ==0 )
            continue;

        float x = git->x[0] + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float y = git->x[1] + gRandom->Gaus(0, sigXY); // 100 micron blur according to approx sTGC reso
        float z = git->x[2];


        
        if (mGenTree && (unsigned)mTreeData.fttN < MAX_TREE_ELEMENTS) {
            mTreeData.fttX.push_back( x );
            mTreeData.fttY.push_back( y );
            mTreeData.fttZ.push_back( z );
            mTreeData.fttTrackId.push_back( track_id );
            mTreeData.fttVolumeId.push_back( plane_id );
            mTreeData.fttPt.push_back( mcTrackMap[track_id]->mPt );
            mTreeData.fttVertexId.push_back( mcTrackMap[track_id]->mStartVertex );
            mTreeData.fttN++;
        } else if ( mGenTree ){
            LOG_WARN << "Truncating hits in TTree output" << endm;
        }

        if ( mGenHistograms ){
            mHistograms["stgc_volume_id"]->Fill(volume_id);
        }

        if (plane_id < 4 && plane_id >= 0) {
            if ( mGenHistograms ){
                mHistograms[TString::Format("stgc%dHitMap", plane_id).Data()]->Fill(x, y);
            }
        } else {
            continue;
        }

        // this rejects GEANT hits with eta -999 - do we understand this effect?
        if ( filterGEANT ) {
            if ( mcTrackMap[track_id] && fabs(mcTrackMap[track_id]->mEta) > 5.0 ){
                
                if ( mGenHistograms ) 
                    mHistograms[TString::Format("stgc%dHitMapSec", plane_id).Data()]->Fill(x, y);
                continue;
            } else if ( mcTrackMap[track_id] && fabs(mcTrackMap[track_id]->mEta) < 5.0 ){
                if ( mGenHistograms ) mHistograms[TString::Format("stgc%dHitMapPrim", plane_id).Data()]->Fill(x, y);
            }
        }

        // FwdHit *hit = new FwdHit(count++, x, y, z, -plane_id, track_id, hitCov3, mcTrackMap[track_id]);
        mFwdHits.push_back(FwdHit(count++, x, y, z, -plane_id, track_id, hitCov3, mcTrackMap[track_id]));
        FwdHit *hit = &(mFwdHits.back());

        // Add the hit to the hit map
        hitMap[mFwdHits.back().getSector()].push_back( hit );
        mFttHits.push_back( TVector3( x, y, z )  );

        // Add hit pointer to the track
        if (mcTrackMap[track_id]){
            mcTrackMap[track_id]->addHit(hit);
            
            if ( mGenTree ) mTreeData.mcNumFtt[track_id-1] ++;
        } else {
            LOG_ERROR << "Cannot find MC track for GEANT hit (FTT), track_id = " << track_id << endm;
        }
    } // loop on hits

    if (mGenTree){
        LOG_INFO << "Saving " << mTreeData.fttN << " hits in Tree" << endm;
    }
} // loadFttHits

/**
 * @brief Loads FST hits from various sources into the hitmap and McTrackMap (if availabale)
 * 
 * Order of precedence:
 * MuDst StMuFstCollection (Data)
 * StEvent StFstHitCollection (Data or slowsim)
 * StEvent StRndHitCollection (fast sim)
 * GEANT St_g2t_fts_hit (starsim only) - note if rasterizer is active this takes priority over FastSim
 * 
 * @param mcTrackMap : MC track map if running sim
 * @param hitMap : FST hitmap to populate
 * @param count  : number of hits loaded
 */
int StFwdTrackMaker::loadFstHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap ){
    
    int count = loadFstHitsFromMuDst(mcTrackMap, hitMap); 
    cout << "We have " << count << " Fst Hits From MuDst " << endl;
    if ( count > 0 ) return count; // only load from one source at a time
    cout << "Are we doing anything else? " << endl;    

    count += loadFstHitsFromStEvent(mcTrackMap, hitMap); 
    if ( count > 0 ) return count; // only load from one source at a time
    
    bool siRasterizer = mFwdConfig.get<bool>( "SiRasterizer:active", false );

    if ( !siRasterizer ) count += loadFstHitsFromStEventFastSim( mcTrackMap, hitMap );
    if ( count > 0 ) return count; // only load from one source at a time

    return loadFstHitsFromGEANT( mcTrackMap, hitMap );    
} // loadFstHits

int StFwdTrackMaker::loadFstHitsFromMuDst( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
    int count = 0;
    StMuDstMaker *mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
    if(!mMuDstMaker) {
        LOG_WARN << " No MuDstMaker ... bye-bye" << endm;
        return 0;
    }
    StMuDst *mMuDst = mMuDstMaker->muDst();
    if(!mMuDst) {
        LOG_WARN << " No MuDst ... bye-bye" << endm;
        return 0;
    }
    
    StMuFstCollection * fst = mMuDst->muFstCollection();
    if (!fst) {
        LOG_WARN << "No StMuFstCollection ... bye-bye" << endm;
        return 0;
    }

    LOG_INFO << "Loading " << fst->numberOfHits() << " StMuFstHits" << endm;

    double fstDefaultZ[12] = {150.008101+0.0009308,151.403100+0.000946,153.491899+0.0010220,152.096900+0.0010770,
                              166.989901+0.0010680,165.594901+0.001053,163.506101+0.0009766,164.901101+0.0009918,
                              177.039106+0.0009308,178.434106+0.000946,180.522905+0.0010220,179.127906+0.0010070};

    double PI = TMath::Pi();
    double PHI[] = {5.*PI/12.,  3.*PI/12.,  1.*PI/12.,  -1.*PI/12.,
                    -3.*PI/12., -5.*PI/12., -7.*PI/12., -9.*PI/12.,
                    -11.*PI/12.,-13.*PI/12.,-15.*PI/12.,-17.*PI/12.};

    TMatrixDSym hitCov3(3);
    for ( unsigned int index = 0; index < fst->numberOfHits(); index++){
        StMuFstHit * muFstHit = fst->getHit( index );

        float vR = muFstHit->localPosition(0);
        float vPhi = muFstHit->localPosition(1);
        float vZ = muFstHit->localPosition(2);

        float x0 = vR * cos( vPhi );
        float y0 = vR * sin( vPhi );

        cout << "vR = " << vR << endl;
        cout << "vPhi = " << vPhi << endl;

        //const float dz0 = fabs( vZ - 151.75 );
        //const float dz1 = fabs( vZ - 165.248 );
        //const float dz2 = fabs( vZ - 178.781 );

        //int d = 0 * ( dz0 < 1.0 ) + 1 * ( dz1 < 1.0 ) + 2 * ( dz2 < 1.0 );

        int iw = muFstHit->getWedge()-1;
        int is = muFstHit->getSensor();

        LOG_INFO << "iw = " << iw << ", is = " << is << endm;

        int d = iw/12;
        int ds = (is == 0)? 0 : 1; // +0 for inner, +1 for outer
        int defaultZidx = d * 4 + 2 * (iw % 2) + ds;
        vZ = fstDefaultZ[defaultZidx];
      
        cout << "defaultZidx = " << defaultZidx << endl;
        cout << "vZ = " << vZ << endl;

        int sensorIdx = iw*3 + is; 
        cout << "sensorIdx = " << sensorIdx << endl;
        
        //double phiAxisShift = 0.0;
        //if((d) == 0 || (d) == 2)
        //{
        //  if((iw)%2 == 0)
        //  {
        //    if((is)%3 == 1)      phiAxisShift = +8.0 * M_PI / 180.0;
        //    else if((is)%3 == 2) phiAxisShift = -8.0 * M_PI / 180.0;
        //  }
        //  else if((iw)%2 == 1)
        //  {
        //    if((is)%3 == 1)      phiAxisShift = -8.0 * M_PI / 180.0;
        //    else if((is)%3 == 2) phiAxisShift = +8.0 * M_PI / 180.0;
        //  } 
        //}
        //else if((d) == 1)
        //{
        //  if((iw)%2 == 0)
        //  {
        //    if((is)%3 == 1)      phiAxisShift = -8.0 * M_PI / 180.0;
        //    else if((is)%3 == 2) phiAxisShift = +8.0 * M_PI / 180.0;
        //  }
        //  else if((iw)%2 == 1)
        //  {
        //    if((is)%3 == 1)      phiAxisShift = +8.0 * M_PI / 180.0;
        //    else if((is)%3 == 2) phiAxisShift = -8.0 * M_PI / 180.0;
        //  } 
        //}
        // 
        //double xtemp = x0 *       TMath::Cos(PHI[iw] + phiAxisShift) + y0 * TMath::Sin(PHI[iw] + phiAxisShift);
        //double ytemp = x0 * -1. * TMath::Sin(PHI[iw] + phiAxisShift) + y0 * TMath::Cos(PHI[iw] + phiAxisShift);

        hitCov3 = makeSiCovMat( TVector3( x0, x0, vZ ), mFwdConfig );

        LOG_INFO << "FST HIT: d = " << d << ", x=" << x0 << ", y=" << y0 << ", z=" << vZ << endm;                
        mFstHits.push_back( TVector3( x0, y0, vZ)  );

        // we use d+4 so that both FTT and FST start at 4
        mFwdHits.push_back(FwdHit(count++, x0, y0, vZ, d+4, 0, hitCov3, nullptr));
        cout << "d+4 for setting layer = " << d+4 << endl;
        auto hit = &(mFwdHits.back());
        cout << "hit->getLayer() = " << hit->getLayer() << endl;
        hit->setSensor(sensorIdx);
        // Add the hit to the hit map
        hitMap[d+4].push_back(hit);

        mTreeData.fstX.push_back( x0 );
        mTreeData.fstY.push_back( y0 );
        mTreeData.fstZ.push_back( vZ );

        mTreeData.fstN++;
        //count++;
    } // index
    return count;
} // loadFstHitsFromMuDst

int StFwdTrackMaker::loadFstHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
    int count = 0;
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    if (!event) {
        LOG_WARN << "No StEvent, cannot load FST hits from StEvent StFstHitCollection" << endm;
        return 0;
    }
    LOG_DEBUG << "Got StEvent, loading Fst Hits" << endm;
    StFstHitCollection *fstHitCollection = event->fstHitCollection();

    double mFstDefaultZ[12] = {150.008101+0.0009308,151.403100+0.000946,153.491899+0.0010220,152.096900+0.0010770,
                               166.989901+0.0010680,165.594901+0.001053,163.506101+0.0009766,164.901101+0.0009918,
                               177.039106+0.0009308,178.434106+0.000946,180.522905+0.0010220,179.127906+0.0010070};

    if ( fstHitCollection && fstHitCollection->numberOfHits() > 0){
        // reuse this to store cov mat
        TMatrixDSym hitCov3(3);
        LOG_DEBUG << "StFstHitCollection is NOT NULL, loading hits" << endm;
        for ( unsigned int iw = 0; iw < kFstNumWedges; iw++ ){
            StFstWedgeHitCollection * wc = fstHitCollection->wedge( iw );
            if ( !wc ) continue;
            for ( unsigned int is = 0; is < kFstNumSensorsPerWedge; is++ ){
                StFstSensorHitCollection * sc = wc->sensor( is );
                if ( !sc ) continue;
                StSPtrVecFstHit fsthits = sc->hits();
                if ( mGenTree ) mTreeData.fstN = 0;
                for ( unsigned int ih = 0; ih < fsthits.size(); ih++ ){
                    float vR = fsthits[ih]->localPosition(0);
                    float vPhi = fsthits[ih]->localPosition(1);
                    float vZ = fsthits[ih]->localPosition(2);

                    //const float dz0 = fabs( vZ - 151.75 );
                    //const float dz1 = fabs( vZ - 165.248 );
                    //const float dz2 = fabs( vZ - 178.781 );

                    //int d = 0 * ( dz0 < 1.0 ) + 1 * ( dz1 < 1.0 ) + 2 * ( dz2 < 1.0 );
                    int d = iw/12;
                    int ds = (is == 0)? 0 : 1; // +0 for inner, +1 for outer
                    int defaultZidx = d * 4 + 2 * (iw % 2) + ds;
                    vZ = mFstDefaultZ[defaultZidx];

                    int sensorIdx = iw*3 + is; 

                    float x0 = vR * cos( vPhi );
                    float y0 = vR * sin( vPhi );
                    hitCov3 = makeSiCovMat( TVector3( x0, y0, vZ ), mFwdConfig );

                    LOG_DEBUG << "FST HIT: d = " << d << ", x=" << x0 << ", y=" << y0 << ", z=" << vZ << endm;                
                    mFstHits.push_back( TVector3( x0, y0, vZ)  );

                    // we use d+4 so that both FTT and FST start at 4
                    mFwdHits.push_back(FwdHit(count++, x0, y0, vZ, d+4, 0, hitCov3, nullptr));
                    auto hit = &(mFwdHits.back());
                    hit->setSensor(sensorIdx);
                    // Add the hit to the hit map
                    hitMap[d+4].push_back(hit);
                    cout << "Layer = " << hit->getLayer() << endl;

                    if ( mGenTree ){
                        mTreeData.fstX.push_back( x0 );
                        mTreeData.fstY.push_back( y0 );
                        mTreeData.fstZ.push_back( vZ );

                        mTreeData.fstN++;
                    }
                    count++;
                }
            } // loop is
        } // loop iw
        LOG_DEBUG << " FOUND " << mFstHits.size() << " FST HITS in StFstHitCollection" << endm;
    } // fstHitCollection
    return count;
} //loadFstHitsFromStEvent

int StFwdTrackMaker::loadFstHitsFromStEventFastSim( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap){
   {
        int countMc = 0;
        St_g2t_fts_hit *g2t_fsi_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_fsi_hit");

        if ( !g2t_fsi_hits )
            return countMc;

        int nfsi = g2t_fsi_hits->GetNRows();


        // reuse this to store cov mat
        TMatrixDSym hitCov3(3);

        if ( mGenHistograms ) this->mHistograms["nHitsFSI"]->Fill(nfsi);
        // LOG_INFO << "# fsi hits = " << nfsi << endm;

        mTreeData.fstN = 0;
        for (int i = 0; i < nfsi; i++) {

            g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_fsi_hits->At(i);

            if (0 == git)
                continue; // geant hit

            int track_id = git->track_p;
            int volume_id = git->volume_id;  // 4, 5, 6
            int d = volume_id / 1000;        // disk id
            int w = (volume_id % 1000) / 10; // wedge id
            int s = volume_id % 10;          // sensor id
            int globalSensorId = (d-4)*36 + (w-1)*3 + (s-1);
    
            int plane_id = d - 4;
            float x = git->x[0];
            float y = git->x[1];
            float z = git->x[2];

            cout << "FST MC HIT: x = " << git->x[0] << ", y = " << git->x[1] << ", z = " << git->x[2] << endl;

            hitCov3 = makeSiCovMat( TVector3( x, y, z ), mFwdConfig );
            FwdHit *hit = new FwdHit(countMc++, x, y, z, d, track_id, hitCov3, mcTrackMap[track_id]);
            hit->setSensor(globalSensorId);
            //mFstHits.push_back( TVector3( x, y, z )  );

            //mTreeData.fstX.push_back( x );
            //mTreeData.fstY.push_back( y );
            //mTreeData.fstZ.push_back( z );
            //mTreeData.fstTrackId.push_back( track_id );

            //mTreeData.fstN++;

            // Add the hit to the hit map
            if (mcTrackMap[track_id])
              mcTrackMap[track_id]->addFstHit(hit);
        }
    }

    int count = 0;
    // Get the StEvent handle
    StEvent *event = (StEvent *)GetDataSet("StEvent");
    if (!event) {
        LOG_DEBUG << "No StEvent, cannot load FST FastSim hits from StEvent StRndHitCollection" << endm;
        return 0;
    }

    StRnDHitCollection *rndCollection = event->rndHitCollection();
    if (!rndCollection) return 0;

    const StSPtrVecRnDHit &hits = rndCollection->hits();

    // we will reuse this to hold the cov mat
    TMatrixDSym hitCov3(3);
    
    if (mGenTree) mTreeData.fstN = 0;
    for (unsigned int fsthit_index = 0; fsthit_index < hits.size(); fsthit_index++) {
        StRnDHit *hit = hits[fsthit_index];
    
        //cout << "ABOUT TO ADD RC HITS  " <<endl;

        if ( hit->layer() > 6  || hit->layer() < 4){
            // skip sTGC hits here
            //cout << "JUST KIDDING DON'T  " <<endl;
            continue;
        }

        const StThreeVectorF pos = hit->position();

        StMatrixF covmat = hit->covariantMatrix();

        // copy covariance matrix element by element from StMatrixF
        hitCov3(0,0) = covmat[0][0]; hitCov3(0,1) = covmat[0][1]; hitCov3(0,2) = covmat[0][2];
        hitCov3(1,0) = covmat[1][0]; hitCov3(1,1) = covmat[1][1]; hitCov3(1,2) = covmat[1][2];
        hitCov3(2,0) = covmat[2][0]; hitCov3(2,1) = covmat[2][1]; hitCov3(2,2) = covmat[2][2];

        mFwdHits.push_back(FwdHit(count++, hit->position().x(), hit->position().y(), hit->position().z(), hit->layer(), hit->idTruth(), hitCov3, mcTrackMap[hit->idTruth()]));
        auto fhit = &(mFwdHits.back());
        fhit->setSensor(hit->volumeId()); 
        cout << "FST RC HIT: x = " << hit->position().x() << ", y = " << hit->position().y() << ", z = " << hit->position().z() << endl;
        
        size_t index = hit->layer()-4;
        if (mGenHistograms && index < 3 ){
            ((TH2*)mHistograms[TString::Format("fsi%luHitMapZ", index).Data()]) -> Fill( hit->position().x(), hit->position().y(), hit->position().z() );
        }

        // Add the hit to the hit map
        hitMap[fhit->getSector()].push_back(fhit);
        mFstHits.push_back( TVector3( hit->position().x(), hit->position().y(), hit->position().z())  );

        if ( mGenTree ){
            mTreeData.fstX.push_back( hit->position().x() );
            mTreeData.fstY.push_back( hit->position().y() );
            mTreeData.fstZ.push_back( hit->position().z() );
            mTreeData.fstTrackId.push_back( hit->idTruth() );

            mTreeData.fstN++;
        }
        count++;

    }
    return count;
} //loadFstHitsFromStEventFastSim

int StFwdTrackMaker::loadFstHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap ){
    int count = 0;
    LOG_DEBUG << "Looking for FST hits in geant struct" << endm;
    /************************************************************/
    // Load FSI Hits from GEANT
    St_g2t_fts_hit *g2t_fsi_hits = (St_g2t_fts_hit *)GetDataSet("geant/g2t_fsi_hit");

    if ( !g2t_fsi_hits ){
        LOG_DEBUG << "No g2t_fts_hits, cannot load FST hits from GEANT" << endm;
        return 0;
    }

    int nfsi = g2t_fsi_hits->GetNRows();
    

    // reuse this to store cov mat
    TMatrixDSym hitCov3(3);
    
    if ( mGenHistograms ) mHistograms["nHitsFSI"]->Fill(nfsi);

    if ( mGenTree ) mTreeData.fstN = 0;
    for (int i = 0; i < nfsi; i++) {

        g2t_fts_hit_st *git = (g2t_fts_hit_st *)g2t_fsi_hits->At(i);
        
        if (0 == git)
            continue; // geant hit
        
        int track_id = git->track_p;
        int volume_id = git->volume_id;  // 4, 5, 6
        int d = volume_id / 1000;        // disk id
        
        int plane_id = d - 4;
        float x = git->x[0];
        float y = git->x[1];
        float z = git->x[2];

        if (mSiRasterizer->active()) {
            TVector3 rastered = mSiRasterizer->raster(TVector3(git->x[0], git->x[1], git->x[2]));
            
            if ( mGenHistograms ) {
                mHistograms["fsiHitDeltaR"]->Fill(std::sqrt(x * x + y * y) - rastered.Perp());
                mHistograms["fsiHitDeltaPhi"]->Fill(std::atan2(y, x) - rastered.Phi());
            }
            x = rastered.X();
            y = rastered.Y();
        }


        if ( mGenHistograms ) mHistograms["fsi_volume_id"]->Fill(d);

        if (plane_id < 3 && plane_id >= 0) {

            if ( mGenHistograms ) {
                mHistograms[TString::Format("fsi%dHitMap", plane_id).Data()]->Fill(x, y);
                mHistograms[TString::Format("fsi%dHitMapR", plane_id).Data()]->Fill(std::sqrt(x * x + y * y));
                mHistograms[TString::Format("fsi%dHitMapPhi", plane_id).Data()]->Fill(std::atan2(y, x) + TMath::Pi());
            }
        } else {
            continue;
        }

        hitCov3 = makeSiCovMat( TVector3( x, y, z ), mFwdConfig );
        mFwdHits.push_back(FwdHit(count++, x, y, z, d, track_id, hitCov3, mcTrackMap[track_id]));
        auto hit = &(mFwdHits.back());

        mFstHits.push_back( TVector3( x, y, z )  );

        if ( mGenTree ){
            mTreeData.fstX.push_back( x );
            mTreeData.fstY.push_back( y );
            mTreeData.fstZ.push_back( z );
            mTreeData.fstTrackId.push_back( track_id );

            mTreeData.fstN++;
        }
        count++;

        // Add the hit to the hit map
        hitMap[hit->getSector()].push_back(hit);

        // Add hit pointer to the track
        if (mcTrackMap[track_id]){
            mcTrackMap[track_id]->addFstHit(hit);
            if ( mGenTree ) {
                mTreeData.mcNumFst[track_id-1] ++;
                LOG_DEBUG << TString::Format("Adding FST Hit to McTrack, mcPt=%f", mTreeData.mcPt[track_id-1] ) << endm;
            }
        } else {
            LOG_ERROR << "Cannot find MC track for GEANT hit (FTT), track_id = " << track_id << endm;
        }
    }
    return count;
} // loadFstHitsFromGEANT

size_t StFwdTrackMaker::loadMcTracks( FwdDataSource::McTrackMap_t &mcTrackMap ){
    

    LOG_DEBUG << "Looking for GEANT sim vertex info" << endm;
    St_g2t_vertex *g2t_vertex = (St_g2t_vertex *)GetDataSet("geant/g2t_vertex");

    if ( g2t_vertex != nullptr ) {
        // Set the MC Vertex for track fitting
        g2t_vertex_st *vert = (g2t_vertex_st*)g2t_vertex->At(0);
        mForwardTracker->setEventVertex( TVector3( vert->ge_x[0], vert->ge_x[1], vert->ge_x[2] ) );
    }

    // Get geant tracks
    St_g2t_track *g2t_track = (St_g2t_track *)GetDataSet("geant/g2t_track");

    if (!g2t_track)
        return 0;

    size_t nShowers = 0;

    mTreeData.mcN = 1;
    LOG_DEBUG << g2t_track->GetNRows() << " mc tracks in geant/g2t_track " << endm;
    if ( mGenHistograms ) mHistograms["nMcTracks"]->Fill(g2t_track->GetNRows());

    for (int irow = 0; irow < g2t_track->GetNRows(); irow++) {
        g2t_track_st *track = (g2t_track_st *)g2t_track->At(irow);

        if (0 == track)
            continue;

        int track_id = track->id;
        float pt2 = track->p[0] * track->p[0] + track->p[1] * track->p[1];
        float pt = std::sqrt(pt2);
        float eta = track->eta;
        float phi = std::atan2(track->p[1], track->p[0]); //track->phi;
        int q = track->charge;

        if (!mcTrackMap[track_id] ) 
            mcTrackMap[track_id] = shared_ptr<McTrack>(new McTrack(pt, eta, phi, q, track->start_vertex_p));
        
        if (mGenTree && (unsigned)mTreeData.mcN < MAX_TREE_ELEMENTS) {
            mTreeData.mcPt.push_back( pt );
            mTreeData.mcEta.push_back( eta );
            mTreeData.mcPhi.push_back( phi );
            mTreeData.mcCharge.push_back( q );
            mTreeData.mcVertexId.push_back( track->start_vertex_p );
            mTreeData.mcNumFtt.push_back(0);
            mTreeData.mcNumFst.push_back(0);

            if (track->is_shower)
                nShowers++;

            mTreeData.mcN++;
        } else if ( mGenTree ) {
            LOG_WARN << "Truncating Mc tracks in TTree output" << endm;
        }

    } // loop on track (irow)


    // now check the Mc tracks against the McEvent filter
    size_t nForwardTracks = 0;
    size_t nForwardTracksNoThreshold = 0;
    for (auto mctm : mcTrackMap ){
        if ( mctm.second == nullptr ) continue;

        if ( mGenHistograms ){
            mHistograms[ "McEventPt" ] ->Fill( mctm.second->mPt );
            mHistograms[ "McEventEta" ] ->Fill( mctm.second->mEta );
            mHistograms[ "McEventPhi" ] ->Fill( mctm.second->mPhi );
        }

        if ( mctm.second->mEta > 2.5 && mctm.second->mEta < 4.0 ){
            
            if ( mGenHistograms ){
                mHistograms[ "McEventFwdPt" ] ->Fill( mctm.second->mPt );
                mHistograms[ "McEventFwdEta" ] ->Fill( mctm.second->mEta );
                mHistograms[ "McEventFwdPhi" ] ->Fill( mctm.second->mPhi );
            }

            nForwardTracksNoThreshold++;
            if ( mctm.second->mPt > 0.05  )
                nForwardTracks++;
        }
    } // loop on mcTrackMap

    if ( mGenHistograms ) {
        mHistograms[ "nMcTracksFwd" ]->Fill( nForwardTracks );
        mHistograms[ "nMcTracksFwdNoThreshold" ]->Fill( nForwardTracksNoThreshold );
    }


    return nForwardTracks;
} // loadMcTracks

void StFwdTrackMaker::loadFcs( ) {
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    StFcsDb* fcsDb=static_cast<StFcsDb*>(GetDataSet("fcsDb"));
    if ( !stEvent || !fcsDb ){
        return;
    }
    StFcsCollection* fcsCol = stEvent->fcsCollection();
    if ( !fcsCol ){
        return;
    }

    StEpdGeom epdgeo;

    mTreeData.fcsN = 0;
    
    // LOAD ECAL / HCAL CLUSTERS
    for ( int idet = 0; idet  < 4; idet++ ){
        StSPtrVecFcsCluster& clusters = fcsCol->clusters(idet);
        int nc=fcsCol->numberOfClusters(idet);
        for ( int i = 0; i < nc; i++ ){
            StFcsCluster* clu = clusters[i];
            StThreeVectorD xyz = fcsDb->getStarXYZfromColumnRow(clu->detectorId(),clu->x(),clu->y());
            mFcsClusters.push_back( TVector3( xyz.x(), xyz.y(), xyz.z() - 2 ) );
            mFcsClusterEnergy.push_back( clu->energy() );

            mTreeData.fcsX.push_back( xyz.x() );
            mTreeData.fcsY.push_back( xyz.y() );
            mTreeData.fcsZ.push_back( xyz.z() - 2 );
            mTreeData.fcsE.push_back( clu->energy() );
            mTreeData.fcsDet.push_back( idet );
        } // i
    } // idet

    // LOAD PRESHOWER HITS (EPD)
    for ( int det = 4; det < 6; det ++ ) {

        StSPtrVecFcsHit& hits = stEvent->fcsCollection()->hits(det);
        int nh=fcsCol->numberOfHits(det);
        for ( int i = 0; i < nh; i++ ){
            StFcsHit* hit=hits[i];

            if(det==kFcsPresNorthDetId || det==kFcsPresSouthDetId){ //EPD
                double zepd=375.0;
                int pp,tt,n;
                double x[5],y[5];

                if ( hit->energy() < 0.2 ) continue;
                fcsDb->getEPDfromId(det,hit->id(),pp,tt);
                epdgeo.GetCorners(100*pp+tt,&n,x,y);
                double x0 = (x[0] + x[1] + x[2] + x[3]) / 4.0;
                double y0 = (y[0] + y[1] + y[2] + y[3]) / 4.0;
                mFcsPreHits.push_back( TVector3( x0, y0, zepd ) );

                mTreeData.epdX.push_back( x0 );
                mTreeData.epdY.push_back( y0 );
                mTreeData.epdZ.push_back( zepd );
                mTreeData.epdE.push_back( hit->energy() );
                // mTreeData.fcsDet.push_back( det );

            } // if det
        } // for i
    } // for det
} // loadFcs

TVector3 StFwdTrackMaker::GetEventPrimaryVertex(){
    TVector3 pv(0, 0, 0);

    // MuDst only for now
    int count = 0;
    StMuDstMaker *mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
    if(!mMuDstMaker) {
        return pv;
    }
    StMuDst *mMuDst = mMuDstMaker->muDst();
    if(!mMuDst) {
        return pv;
    }

    StMuPrimaryVertex *muPV = mMuDst->primaryVertex();
    if ( muPV ){
        pv.SetX(muPV->position().x());
        pv.SetY(muPV->position().y());
        pv.SetZ(muPV->position().z());
    } else {
        // No primary vertex available 
    }

    return pv;
}

//________________________________________________________________________
int StFwdTrackMaker::Make() {
    //reportMem();
    // START time for measuring tracking
    long long itStart = FwdTrackerUtils::nowNanoSecond();

    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));

    // Access forward Tracker maps
    FwdDataSource::McTrackMap_t &mcTrackMap = mForwardData->getMcTracks();
    FwdDataSource::HitMap_t &hitMap = mForwardData->getFttHits();
    FwdDataSource::HitMap_t &fsiHitMap = mForwardData->getFstHits();
    
    // clear vectors for visualization OBJ hits
    mFttHits.clear();
    mFstHits.clear();
    mFcsPreHits.clear();
    mFcsClusters.clear();
    mFwdTracks.clear();
    
    // default event vertex
    auto eventPV = GetEventPrimaryVertex();
    mForwardTracker->setEventVertex( eventPV );
    if ( fabs(eventPV.X()) < 0.001 && fabs(eventPV.Y()) < 0.001 && fabs(eventPV.Z()) < 0.001 ){
        LOG_DEBUG << "FWD Tracking on event without available Primary Vertex" << endm;
    } else {
        LOG_INFO << "Forward Tracking using primary vertex: " << TString::Format( "(%f, %f, %f)", eventPV.X(), eventPV.Y(), eventPV.Z() ) << endm;
    }


    /**********************************************************************/
    // Load MC tracks
    size_t nForwardTracks = loadMcTracks( mcTrackMap );
    size_t maxForwardTracks = mFwdConfig.get<size_t>( "McEvent.Mult:max", 10000 );
    if ( nForwardTracks > maxForwardTracks ){
        LOG_WARN << "Skipping event with more than " << maxForwardTracks << " forward tracks" << endm;
        return kStOk;
    }
    LOG_DEBUG << "We have " << nForwardTracks << " forward MC tracks" << endm;
    /**********************************************************************/
    // Load sTGC 
    LOG_DEBUG << ">>StFwdTrackMaker::loadFttHits" << endm;
    if ( IAttr("useFtt") ) {
        loadFttHits( mcTrackMap, hitMap );
    }
    

    /**********************************************************************/
    // Load FST
    if ( IAttr("useFst") ) {
        LOG_DEBUG << ">>StFwdTrackMaker::loadFstHits" << endm;
        int fstCount = loadFstHits( mcTrackMap, fsiHitMap );
        LOG_DEBUG << "Loaded " << fstCount << " FST hits" << endm;
    }

    /**********************************************************************/
    // Load FCS
    LOG_DEBUG << ">>StFwdTrackMaker::loadFcsHits" << endm;
    if ( IAttr("useFcs") ) {
        loadFcs();
    }

    /**********************************************************************/
    // Run Track finding + fitting
    LOG_DEBUG << ">>START Event Forward Tracking" << endm;
    mForwardTracker->doEvent();
    LOG_DEBUG << "<<FINISH Event Forward Tracking" << endm;
    LOG_DEBUG << "<<Found " << mForwardTracker -> getTrackSeeds().size() << " Track Seeds" << endm;
    LOG_DEBUG << "Fit " << mForwardTracker -> getTrackResults().size() << " GenFit Tracks" << endm;
    // perform vertex finding/fitting with Fwd tracks
    // FitVertex();

    /**********************************************************************/
    
    /**********************************************************************/
    // Output track visualization if configured to do so
    if ( mVisualize ){
        std::vector<genfit::Track *> genfitTracks;
        for ( auto gtr : mForwardTracker->getTrackResults() ) {
            if ( gtr.isFitConvergedFully == false ) continue;
            genfitTracks.push_back( gtr.track.get() );
        }

        if ( mVisualize && genfitTracks.size() > 0 && genfitTracks.size() < 400 && eventIndex < 50 ) {
            const auto &seed_tracks = mForwardTracker -> getTrackSeeds();

            ObjExporter woe;
            woe.output( 
                TString::Format( "ev%lu", eventIndex ).Data(), 
                stEvent,
                seed_tracks, genfitTracks, mRaveVertices, 
                mFttHits, mFstHits, mFcsPreHits, mFcsClusters, mFcsClusterEnergy );
            eventIndex++;
            LOG_DEBUG << "Done Writing OBJ " << endm;
        } else if (mVisualize && genfitTracks.size() == 0) {
            LOG_DEBUG << "Skipping visualization, no FWD tracks" << endm;
        } else if (mVisualize && genfitTracks.size() >= 400) {
            LOG_DEBUG << "Skipping visualization, too many FWD tracks" << endm;
        }
    }

    // LOG_DEBUG << "FillTrackDeltas" << endm;
    // Fill Track Deltas in ttree for helpful alignment info
    // FillTrackDeltas();

    LOG_DEBUG << "Forward tracking on this event took " << (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6 << " ms" << endm;

    if ( IAttr("fillEvent") ) {
        if (!stEvent) {
            LOG_WARN << "No StEvent found. Forward tracks will not be saved" << endm;
            return kStWarn;
        }
        FillEvent();
    } // IAttr FillEvent

    LOG_DEBUG << "Filling fwd Tree for event: " << GetIventNumber() << endm;
    // FillTTree();
    return kStOK;
} // Make

StFwdTrack * StFwdTrackMaker::makeStFwdTrack( GenfitTrackResult &gtr, size_t indexTrack ){
    LOG_DEBUG << "StFwdTrackMaker::makeStFwdTrack()" << endm;
    StFwdTrack *fwdTrack = new StFwdTrack(  );

    // Fit failed beyond use
    if ( gtr.track == nullptr  ){
        gtr.Clear();
        LOG_DEBUG << "GenfitTrack is nullptr, not making StFwdTrack" << endm;
        return nullptr;
    }    
    // Fill charge and quality info
    fwdTrack->setDidFitConverge( gtr.status->isFitConverged() );
    fwdTrack->setDidFitConvergeFully( gtr.status->isFitConvergedFully() );
    fwdTrack->setNumberOfFailedPoints( gtr.status->getNFailedPoints() );
    
    fwdTrack->setNumberOfFitPoints( gtr.track->getNumPoints() );
    fwdTrack->setChi2( gtr.status->getChi2() );
    fwdTrack->setNDF( gtr.status->getNdf() );
    fwdTrack->setPval( gtr.status->getPVal() );

    // auto cr = gtr.track->getCardinalRep();
    // charge at first point
    fwdTrack->setCharge( gtr.charge );

    TVector3 p = gtr.momentum;//cr->getMom( gtr.track->getFittedState( 0, cr ));
    fwdTrack->setPrimaryMomentum( StThreeVectorD( gtr.momentum.X(), gtr.momentum.Y(), gtr.momentum.Z() ) );
    LOG_DEBUG << "Making StFwdTrack with " << TString::Format( "p=(%f, %f, %f)",  fwdTrack->momentum().x(), fwdTrack->momentum().y(), fwdTrack->momentum().z() ) << endm;

    int nSeedPoints = 0;
    // store the seed points from FTT
    for ( auto s : gtr.fttSeed ){
        FwdHit * fh = static_cast<FwdHit*>( s );
        if (!fh) continue;
        float cov[9];
        cov[0] = fh->_covmat(0,0); cov[3] = fh->_covmat(1,0); cov[6] = fh->_covmat(2,0);
        cov[1] = fh->_covmat(0,1); cov[4] = fh->_covmat(1,1); cov[7] = fh->_covmat(2,1);
        cov[2] = fh->_covmat(0,2); cov[5] = fh->_covmat(1,2); cov[8] = fh->_covmat(2,2);

        StFwdTrackSeedPoint p( StThreeVectorD( fh->getX(), fh->getY(), fh->getZ() ), fh->getSector(), fh->getTrackId(), cov );
        fwdTrack->mFTTPoints.push_back( p );
        nSeedPoints++;
    }

    for ( auto s : gtr.fstSeed ){
        FwdHit * fh = static_cast<FwdHit*>( s );
        if (!fh) continue;
        float cov[9];
        cov[0] = fh->_covmat(0,0); cov[3] = fh->_covmat(1,0); cov[6] = fh->_covmat(2,0);
        cov[1] = fh->_covmat(0,1); cov[4] = fh->_covmat(1,1); cov[7] = fh->_covmat(2,1);
        cov[2] = fh->_covmat(0,2); cov[5] = fh->_covmat(1,2); cov[8] = fh->_covmat(2,2);

        StFwdTrackSeedPoint p( StThreeVectorD( fh->getX(), fh->getY(), fh->getZ() ), fh->getSector(), fh->getTrackId(), cov );
        fwdTrack->mFSTPoints.push_back( p );
        nSeedPoints++;
    }

    // set total number of seed points
    fwdTrack->setNumberOfSeedPoints( nSeedPoints );
    Seed_t combinedSeed;
    combinedSeed.insert( combinedSeed.begin(), gtr.fstSeed.begin(), gtr.fstSeed.end() ); // this is goofed but will fix
    combinedSeed.insert( combinedSeed.end(), gtr.fttSeed.begin(), gtr.fttSeed.end() );
    int idt = 0;
    double qual = 0;
    idt = MCTruthUtils::dominantContribution(combinedSeed, qual);
    fwdTrack->setMc( idt, qual );

    // compute projections to z-planes of various detectors
    vector<float> zPlanes = {
        0, // PV TODO, update with event vertex?
        mFstZFromGeom[0], mFstZFromGeom[1], mFstZFromGeom[2], // FST
        mFttZFromGeom[0], mFttZFromGeom[1], mFttZFromGeom[2], // FTT
        375.0, // EPD
        715.0, //ECAL
        807.0 // HCAL
    };
    
    // Match these to the z-planes above
    const int FST = kFstId;
    const int FTT = kFttId;
    vector<int> detMap = {
        kTpcId, 
        FST, FST, FST, 
        FTT, FTT, FTT, FTT, 
        kFcsPresId,
        kFcsWcalId,
        kFcsHcalId
    };

    size_t zIndex = 0;
    int detIndex = 0;
    TVector3 mom(0, 0, 0);
    float cov[9];
    TVector3 tv3(0, 0, 0);
    for ( float z : zPlanes ){
        detIndex = detMap[ zIndex];
        LOG_DEBUG << "Calculating Fwd Projection for detId=" << detIndex << " @ z=" << z << endm;
        tv3.SetXYZ(0, 0, 0);
        if ( detIndex != kFcsHcalId && detIndex != kFcsWcalId ){
            tv3 = ObjExporter::trackPosition( gtr.track.get(), z, cov, mom );
        } else {
            // use a straight line projection to HCAL since GenFit cannot handle long projections
            tv3 = ObjExporter::projectAsStraightLine( gtr.track.get(), 575.0, 625.0, z, cov, mom );
        }
        fwdTrack->mProjections.push_back( StFwdTrackProjection( detIndex, StThreeVectorF( tv3.X(), tv3.Y(), tv3.Z() ), StThreeVectorF( mom.X(), mom.Y(), mom.Z() ), cov) );

        // // Add Proj info to TTree
        if ( mGenTree ){
            if ( kFcsPresId == detIndex ){
                mTreeData.rcEpdX.push_back( tv3.X() ); 
                mTreeData.rcEpdY.push_back( tv3.Y() ); 
            } else if ( kFcsHcalId == detIndex ){

            } else if ( kFcsWcalId == detIndex ){

            }
            
        }

        zIndex++;
    }
    gtr.Clear();

    return fwdTrack;   
}

void StFwdTrackMaker::FillEvent() {
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent)
        return;
    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    if ( !ftc ){
        LOG_INFO << "Creating the StFwdTrackCollection" << endm;
        ftc = new StFwdTrackCollection();
        stEvent->setFwdTrackCollection( ftc );
    }

    if ( mGenTree ){
        mTreeData.tprojN = 0;
        mTreeData.tprojX.clear();
        mTreeData.tprojY.clear();
        mTreeData.tprojZ.clear();
        mTreeData.tprojIdD.clear();
        mTreeData.tprojIdT.clear();
    }

    size_t indexTrack = 0; 
    for ( auto gtr : mForwardTracker->getTrackResults() ) {
            StFwdTrack* fwdTrack = makeStFwdTrack( gtr, indexTrack );
            indexTrack++;
            if (nullptr == fwdTrack)
                continue;
            ftc->addTrack( fwdTrack );
    }


    if ( mGenTree ) mTreeData.tprojN = mTreeData.tprojX.size();
    LOG_INFO << "StFwdTrackCollection has " << ftc->numberOfTracks() << " tracks now" << endm;
    // ProcessFwdTracks();
}

void StFwdTrackMaker::FillTrackDeltas(){
    const auto &fittedTracks = mForwardTracker -> getTrackResults();

    for ( size_t i = 0; i < fittedTracks.size(); i++ ){
        auto gt = fittedTracks[i].track;

        Seed_t combinedSeed;
        combinedSeed.insert( combinedSeed.begin(), fittedTracks[i].fstSeed.begin(), fittedTracks[i].fstSeed.end() ); // this is goofed but will fix
        combinedSeed.insert( combinedSeed.end(), fittedTracks[i].fttSeed.begin(), fittedTracks[i].fttSeed.end() );

        if (fittedTracks[i].isFitConvergedFully == false){
            LOG_DEBUG << "Skipping track, failed fit" << endm;
            continue;
        }

        for ( KiTrack::IHit * hit : combinedSeed ){
            if ( !hit ){
                LOG_ERROR << "Null hit in unexpected place, check track creation" << endm;
                continue; // we should never have null hits, but guard 
            } 
            TVector3 htv3(hit->getX(), hit->getY(), hit->getZ());
            TLorentzVector lvHit, lvTrack;
            lvHit.SetPxPyPzE( hit->getX(), hit->getY(), hit->getZ(), 1 );
            
            auto ttv3 = ObjExporter::trackPosition( gt.get(), htv3.Z() );
            lvTrack.SetPxPyPzE( ttv3.X(), ttv3.Y(), ttv3.Z(), 1 );
            
            mTreeData.thdX.push_back( (ttv3.X() - htv3.X()) );
            mTreeData.thdY.push_back( (ttv3.Y() - htv3.Y()) );

            mTreeData.thdP.push_back( lvTrack.DeltaPhi( lvHit ) );
            mTreeData.thdR.push_back( lvTrack.Pt() - lvHit.Pt() );

            mTreeData.thaX.push_back( htv3.X() );
            mTreeData.thaY.push_back( htv3.Y() );
            mTreeData.thaZ.push_back( htv3.Z() );
            mTreeData.thdN++;
        } // hit
    } // i
} //FillTrackDeltas

void StFwdTrackMaker::FitVertex(){
    vector<genfit::Track *> genfitTracks;

    const auto &trackResults = mForwardTracker -> getTrackResults();
    // for ( auto gtr : trackResults ){
    //     // genfitTracks.push_back( gtr.track );
    // }
    if ( genfitTracks.size() >= 2 ){
        genfit::GFRaveVertexFactory gfrvf;

        TMatrixDSym bscm(3);
        const double bssXY = 2.0;
        bscm(0, 0) = bssXY*bssXY;
        bscm(1, 1) = bssXY*bssXY;
        bscm(2, 2) = 50.5 * 50.5;
        gfrvf.setBeamspot( TVector3( 0, 0, 0 ), bscm );
        
        mRaveVertices.clear();
        gfrvf.findVertices( &mRaveVertices, genfitTracks, false );
        
        LOG_DEBUG << "mRaveVertices.size() = " << mRaveVertices.size() << endm;
        for ( auto vert : mRaveVertices ){
            LOG_DEBUG << TString::Format( "RAVE vertex @(%f, %f, %f)\n\n", vert->getPos().X(), vert->getPos().Y(), vert->getPos().Z() ) << endm;
        }
    }
} // FitVertex

void StFwdTrackMaker::FillTTree(){

    St_g2t_vertex *g2t_vertex = (St_g2t_vertex *)GetDataSet("geant/g2t_vertex");
    if (mGenTree) {
        
        // VERTICES
        if ( g2t_vertex ){
            mTreeData.vmcN = g2t_vertex->GetNRows();
            if ( (unsigned)mTreeData.vmcN >= MAX_TREE_ELEMENTS ) mTreeData.vmcN = MAX_TREE_ELEMENTS;
            LOG_INFO << "Saving " << mTreeData.vmcN << " vertices in TTree" << endm;
            for ( int i = 0; i < mTreeData.vmcN; i++ ){
                g2t_vertex_st *vert = (g2t_vertex_st*)g2t_vertex->At(i);
                mTreeData.vmcX.push_back( vert->ge_x[0] );
                mTreeData.vmcY.push_back( vert->ge_x[1] );
                mTreeData.vmcZ.push_back( vert->ge_x[2] );
            }
        }

        // RAVE RECO VERTICES
        mTreeData.vrcN = mRaveVertices.size();
        if ( (unsigned)mTreeData.vrcN >= MAX_TREE_ELEMENTS ) mTreeData.vrcN = MAX_TREE_ELEMENTS;
        LOG_INFO << "Saving " << mTreeData.vrcN << " RAVE vertices in TTree" << endm;
        for ( int i = 0; i < mTreeData.vrcN; i++ ) {
            auto vert = mRaveVertices[i];
            mTreeData.vrcX.push_back( vert->getPos().X() );
            mTreeData.vrcY.push_back( vert->getPos().Y() );
            mTreeData.vrcZ.push_back( vert->getPos().Z() );
        }




        if (mForwardTracker->getSaveCriteriaValues() && mTreeData.saveCrit ) {
            for (auto crit : mForwardTracker->getTwoHitCriteria()) {
                string name = crit->getName();

                // special, save all hit info for this one
                

                if ( name == "Crit2_BDT" ){
                    mTreeData.Crits["Crit2_BDT_DeltaPhi"].clear(); 
                    mTreeData.Crits["Crit2_BDT_DeltaRho"].clear(); 
                    mTreeData.Crits["Crit2_BDT_RZRatio"].clear(); 
                    mTreeData.Crits["Crit2_BDT_StraightTrackRatio"].clear(); 

                    for (auto kv : mForwardTracker->getCriteriaAllValues(name)) {
                        mTreeData.Crits["Crit2_BDT_DeltaPhi"].push_back( kv["Crit2_BDT_DeltaPhi"] );
                        mTreeData.Crits["Crit2_BDT_DeltaRho"].push_back( kv["Crit2_BDT_DeltaRho"] );
                        mTreeData.Crits["Crit2_BDT_RZRatio"].push_back( kv["Crit2_BDT_RZRatio"] );
                        mTreeData.Crits["Crit2_BDT_StraightTrackRatio"].push_back( kv["Crit2_BDT_StraightTrackRatio"] );
                    }

                }

                if ( name == "Crit2_RZRatio" ){
                    LOG_INFO << "allValues.size() = " << mForwardTracker->getCriteriaAllValues(name).size() << " == " << mForwardTracker->getCriteriaTrackIds(name).size() << endm;
                    assert( mForwardTracker->getCriteriaAllValues(name).size() == mForwardTracker->getCriteriaTrackIds(name).size() && " Crit lengths must be equal" );
                    mTreeData.Crits["Crit2_RZRatio_x1"].clear();
                    mTreeData.Crits["Crit2_RZRatio_y1"].clear();
                    mTreeData.Crits["Crit2_RZRatio_z1"].clear();
                    mTreeData.Crits["Crit2_RZRatio_x2"].clear();
                    mTreeData.Crits["Crit2_RZRatio_y2"].clear();
                    mTreeData.Crits["Crit2_RZRatio_z2"].clear();

                    mTreeData.CritTrackIds["Crit2_RZRatio_h1"].clear();
                    mTreeData.CritTrackIds["Crit2_RZRatio_h2"].clear();
                    mTreeData.CritTrackIds["Crit2_RZRatio_h3"].clear();
                    

                    for (auto kv : mForwardTracker->getCriteriaAllValues(name)) {
                        mTreeData.Crits["Crit2_RZRatio_x1"].push_back( kv["Crit2_RZRatio_x1"] );
                        mTreeData.Crits["Crit2_RZRatio_y1"].push_back( kv["Crit2_RZRatio_y1"] );
                        mTreeData.Crits["Crit2_RZRatio_z1"].push_back( kv["Crit2_RZRatio_z1"] );

                        mTreeData.Crits["Crit2_RZRatio_x2"].push_back( kv["Crit2_RZRatio_x2"] );
                        mTreeData.Crits["Crit2_RZRatio_y2"].push_back( kv["Crit2_RZRatio_y2"] );
                        mTreeData.Crits["Crit2_RZRatio_z2"].push_back( kv["Crit2_RZRatio_z2"] );

                        mTreeData.CritTrackIds["Crit2_RZRatio_h1"].push_back( kv["Crit2_RZRatio_h1"] );
                        mTreeData.CritTrackIds["Crit2_RZRatio_h2"].push_back( kv["Crit2_RZRatio_h2"] );
                        mTreeData.CritTrackIds["Crit2_RZRatio_h3"].push_back( -1 );
                    }
                }


                LOG_DEBUG << "Saving Criteria values from " << name << " in TTree" << endm;
                mTreeData.Crits[name].clear();
                mTreeData.CritTrackIds[name].clear();
                // copy by value so ROOT doesnt get lost (uses pointer to vector)
                for (float v : mForwardTracker->getCriteriaValues(name)) {
                    mTreeData.Crits[name].push_back(v);
                }
                for (int v : mForwardTracker->getCriteriaTrackIds(name)) {
                    mTreeData.CritTrackIds[name].push_back(v);
                }
            }

            // three hit criteria
            for (auto crit : mForwardTracker->getThreeHitCriteria()) {
                string name = crit->getName();

                // special, save all hit info for this one
                if ( name == "Crit2_RZRatio" ){
                    LOG_INFO << "allValues.size() = " << mForwardTracker->getCriteriaAllValues(name).size() << " == " << mForwardTracker->getCriteriaTrackIds(name).size() << endm;
                    assert( mForwardTracker->getCriteriaAllValues(name).size() == mForwardTracker->getCriteriaTrackIds(name).size() && " Crit lengths must be equal" );

                    mTreeData.CritTrackIds["Crit2_RZRatio_h1"].clear();
                    mTreeData.CritTrackIds["Crit2_RZRatio_h2"].clear();
                    mTreeData.CritTrackIds["Crit2_RZRatio_h3"].clear();

                    for (auto kv : mForwardTracker->getCriteriaAllValues(name)) {
                        mTreeData.CritTrackIds["Crit2_RZRatio_h1"].push_back( kv["Crit2_RZRatio_h1"] );
                        mTreeData.CritTrackIds["Crit2_RZRatio_h2"].push_back( kv["Crit2_RZRatio_h2"] );
                        mTreeData.CritTrackIds["Crit2_RZRatio_h3"].push_back( kv["Crit2_RZRatio_h3"] );
                    }
                }


                LOG_DEBUG << "Saving Criteria values from " << name << " in TTree" << endm;
                mTreeData.Crits[name].clear();
                mTreeData.CritTrackIds[name].clear();
                // copy by value so ROOT doesnt get lost (uses pointer to vector)
                for (float v : mForwardTracker->getCriteriaValues(name)) {
                    mTreeData.Crits[name].push_back(v);
                }
                for (int v : mForwardTracker->getCriteriaTrackIds(name)) {
                    mTreeData.CritTrackIds[name].push_back(v);
                }
            }

            // clear them 
            mForwardTracker->clearSavedCriteriaValues();
        }

        // SAVE RECO tracks

        mTreeData.rcN = 0;
        const auto &fittedTracks = mForwardTracker -> getTrackResults();

        LOG_INFO << "There are " << fittedTracks.size() << " seed tracks to save" << endm;
        size_t maxToSave = fittedTracks.size();
        if (maxToSave >= 200) {
            maxToSave = 0;
            LOG_INFO << "More than 200 tracks , not saving unfit tracks" << endm;
        }
        
        for ( size_t i = 0; i < maxToSave; i++ ){
            if ( i >= MAX_TREE_ELEMENTS ){
                LOG_WARN << "Truncating Reco tracks in TTree output" << endm;
                break;
            }

            int idt = 0;
            double qual = 0;
            // TODO: combine these
            Seed_t combinedSeed;
            combinedSeed.insert( combinedSeed.begin(), fittedTracks[i].fstSeed.begin(), fittedTracks[i].fstSeed.end() ); // this is goofed but will fix
            combinedSeed.insert( combinedSeed.end(), fittedTracks[i].fttSeed.begin(), fittedTracks[i].fttSeed.end() );
            idt = MCTruthUtils::dominantContribution(combinedSeed, qual);
            if ( idt == -1 ){
                idt = MCTruthUtils::dominantContribution(fittedTracks[i].fttSeed, qual);
                if ( idt == -1 ){
                    idt = MCTruthUtils::dominantContribution(fittedTracks[i].fstSeed, qual);
                }
            }

            if ( fittedTracks[i].track == nullptr || fittedTracks[i].trackRep == nullptr ) {
                LOG_INFO << "Skip saving null track" << endm;
                continue;
            }

            if ( fittedTracks[i].isFitConverged == false ){
                LOG_INFO << "Skip saving track where fit did not converge" << endm;
                continue;
            }

            
            mTreeData.rcQuality.push_back( qual );
            mTreeData.rcTrackId.push_back( idt );

            mTreeData.rcCharge.push_back( fittedTracks[i].charge );
            mTreeData.rcPt.push_back( fittedTracks[i].momentum.Pt() );
            mTreeData.rcEta.push_back( fittedTracks[i].momentum.Eta() );
            mTreeData.rcPhi.push_back( fittedTracks[i].momentum.Phi() );

            mTreeData.rcNumPV.push_back( fittedTracks[i].nPV );
            mTreeData.rcNumFTT.push_back( fittedTracks[i].numFTT() );
            mTreeData.rcNumFST.push_back( fittedTracks[i].numFST() );

            mTreeData.rcN ++;
        }
        LOG_INFO << "Filling TTree" << endm;
        mTree->Fill();
    } // if mGenTree
} // FillTree

//________________________________________________________________________
void StFwdTrackMaker::Clear(const Option_t *opts) {
    LOG_DEBUG << "StFwdTrackMaker::CLEAR" << endm;
    mForwardData->clear();
    mForwardTracker->Clear();
    
    // clear fwd hits from fst and ftt
    mFwdHits.clear();

    if (mGenTree){
        mTreeData.clear();
    }
}
//________________________________________________________________________
void StFwdTrackMaker::ProcessFwdTracks(  ){
    // This is an example of how to process fwd track collection
    LOG_DEBUG << "StFwdTrackMaker::ProcessFwdTracks" << endm;
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    StFwdTrackCollection * ftc = stEvent->fwdTrackCollection();
    for ( auto fwdTrack : ftc->tracks() ){
        LOG_DEBUG << TString::Format("StFwdTrack[ nProjections=%lu, nFTTSeeds=%lu, nFSTSeeds=%lu, mPt=%f ]", fwdTrack->mProjections.size(), fwdTrack->mFTTPoints.size(), fwdTrack->mFSTPoints.size(), fwdTrack->momentum().perp()) << endm;
        for ( auto proj : fwdTrack->mProjections ) {
            LOG_DEBUG << TString::Format("Proj[ %d, %f, %f, %f ]", proj.mDetId, proj.mXYZ.x(), proj.mXYZ.y(), proj.mXYZ.z() ) << endm;
        }
    }
}


std::string StFwdTrackMaker::defaultConfigIdealSim = R"(
<?xml version="1.0" encoding="UTF-8"?>
<config>
    <Output url="fwdTrackMaker_ideal_sim.root" />
    <Source ftt="GEANT"  />

	<TrackFitter refit="false" mcSeed="true" active="true">
        <Vertex sigmaXY="0.001" sigmaZ="0.01" includeInFit="true" smearMcVertex="true" />
    </TrackFitter>
</config>
)";



std::string StFwdTrackMaker::defaultConfigData = R"(
<?xml version="1.0" encoding="UTF-8"?>
<config>
    <Output url="stfwdtrackmaker_data.root" />
    <Geometry misaligned="true">fGeom.root</Geometry>
    <Source ftt="DATA" />

    <SiRasterizer r="3.0" phi="0.004" />

    <TrackFinder nIterations="1">
        <Iteration nPhiSlices="2" > <!-- Options for first iteration -->
            <SegmentBuilder>
                <Criteria name="Crit2_RZRatio" min="0" max="1.20" />
                <Criteria name="Crit2_DeltaRho" min="-50" max="50.9"/>
                <Criteria name="Crit2_DeltaPhi" min="0" max="30.0" />
                <Criteria name="Crit2_StraightTrackRatio" min="0.01" max="5.85"/>
            </SegmentBuilder>

            <ThreeHitSegments>
				<Criteria name="Crit3_3DAngle" min="0" max="30" />
                <Criteria name="Crit3_PT" min="0" max="100" />
				<Criteria name="Crit3_ChangeRZRatio" min="0.8" max="1.21" />
				<Criteria name="Crit3_2DAngle" min="0" max="30" />
            </ThreeHitSegments>
        </Iteration>

        <Connector distance="1"/>

        <SubsetNN active="false" min-hits-on-track="3" >
            <!-- <InitialTemp>2.1</InitialTemp> -->
            <!-- <InfTemp>0.1</InfTemp> -->
            <Omega>0.99</Omega>
            <StableThreshold>0.001</StableThreshold>
        </SubsetNN> 

        <HitRemover active="false" />
    </TrackFinder>
    
	<TrackFitter refit="true" mcSeed="false" zeroB="true" refitGBL="true">
        <Vertex sigmaXY="0.01" sigmaZ="0.01" includeInFit="true" smearMcVertex="false" />
    </TrackFitter>
</config>
)";
