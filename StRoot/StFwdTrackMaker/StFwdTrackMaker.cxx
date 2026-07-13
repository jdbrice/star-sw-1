
#include "KiTrack/IHit.h"
#include "GenFit/Track.h"

#include "TMath.h"
#include "StBFChain/StBFChain.h"

#include "StEvent/StEvent.h"
#include "StEvent/StGlobalTrack.h"
#include "StEvent/StHelixModel.h"
#include "StEvent/StPrimaryTrack.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"
#include "StEvent/StTrack.h"
#include "StEvent/StBTofCollection.h"
#include "StEvent/StBTofHeader.h"
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

#include "StMcEvent/StMcEvent.hh"
#include "StMcEvent/StMcVertex.hh"


#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_g2t_event_Table.h"

#include "StarMagField/StarMagField.h"

#include "St_base/StMessMgr.h"
#include "StarClassLibrary/StPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"

#include <SystemOfUnits.h>

#include "TROOT.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StFcsDbMaker/StFcsDb.h"
#include "StFstUtil/StFstCollection.h"

#include "StEvent/StFwdTrack.h"
#include "GenFit/AbsMeasurement.h"
#include "GenFit/KalmanFitterInfo.h"
#include "GenFit/MeasurementOnPlane.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuFstCollection.h"
#include "StMuDSTMaker/COMMON/StMuFstHit.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"

#include <cstdlib>

#define LOG_DEBUG if(false) std::cerr
#define LOG_INFO if(false) std::cerr

#include "StFwdTrackMaker/StFwdTrackMaker.h"
#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StFwdTrackMaker/include/Tracker/FwdTracker.h"
#include "StFwdTrackMaker/include/Tracker/TrackFitter.h"
#include "StFwdTrackMaker/include/Tracker/FwdGeomUtils.h"
#include "StFwdTrackMaker/include/Tracker/ObjExporter.h"

#include <algorithm>
#include <cmath>

FwdSystem* FwdSystem::sInstance = nullptr;

namespace {
    constexpr float kInvalidAlignValue = -99999.0f;

    void setAlignmentVector(const TVectorD &source, float &x0, float &x1, float &x2) {
        x0 = source.GetNrows() > 0 ? source[0] : kInvalidAlignValue;
        x1 = source.GetNrows() > 1 ? source[1] : kInvalidAlignValue;
        x2 = source.GetNrows() > 2 ? source[2] : kInvalidAlignValue;
    }

    void setAlignmentPrediction(
        const TVectorD &measurement,
        const TVectorD &residual,
        float &x0, float &x1, float &x2
    ) {
        x0 = kInvalidAlignValue;
        x1 = kInvalidAlignValue;
        x2 = kInvalidAlignValue;

        float *prediction[] = {&x0, &x1, &x2};
        const int n = std::min(measurement.GetNrows(), residual.GetNrows());
        for (int i = 0; i < 3 && i < n; ++i) {
            // GenFit residuals are stored as measurement minus fitted prediction.
            *prediction[i] = measurement[i] - residual[i];
        }
    }

    void setAlignmentTrackKinematics(
        const TVector3 &momentum,
        float &px, float &py, float &pz,
        float &p, float &pt,
        float &eta
    ) {
        px = momentum.X();
        py = momentum.Y();
        pz = momentum.Z();
        pt = momentum.Perp();
        p = momentum.Mag();
        eta = kInvalidAlignValue;

        if (p > std::abs(static_cast<double>(pz))) {
            eta = 0.5 * std::log((p + pz) / (p - pz));
        }
    }

    void setAlignmentPulls(
        const genfit::MeasurementOnPlane &residual,
        float &sigma0, float &sigma1, float &sigma2,
        float &pull0, float &pull1, float &pull2
    ) {
        sigma0 = kInvalidAlignValue;
        sigma1 = kInvalidAlignValue;
        sigma2 = kInvalidAlignValue;
        pull0 = kInvalidAlignValue;
        pull1 = kInvalidAlignValue;
        pull2 = kInvalidAlignValue;

        const TVectorD &state = residual.getState();
        const TMatrixDSym &cov = residual.getCov();
        const int nState = state.GetNrows();
        const int nCov = cov.GetNrows();

        float *sigmas[] = {&sigma0, &sigma1, &sigma2};
        float *pulls[] = {&pull0, &pull1, &pull2};
        for (int i = 0; i < 3 && i < nState && i < nCov; ++i) {
            const double variance = cov(i, i);
            if (variance > 0 && std::isfinite(variance)) {
                const double sigma = std::sqrt(variance);
                *sigmas[i] = sigma;
                *pulls[i] = state[i] / sigma;
            }
        }
    }

    const FwdHit *findFstSeedHit(const GenfitTrackResult &gtr, int globalSensor) {
        if (globalSensor < 0)
            return nullptr;

        for (auto hit : gtr.mSeed) {
            const FwdHit *fwdHit = dynamic_cast<const FwdHit*>(hit);
            if (!fwdHit || !fwdHit->isFst())
                continue;
            if (static_cast<int>(fwdHit->_genfit_plane_index) == globalSensor)
                return fwdHit;
        }

        return nullptr;
    }
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

//  Wrapper class around the forward tracker
class ForwardTracker : public ForwardTrackMaker {
  public:
    // Replaces original initialization.  Config file and hitloader
    // will be provided by the maker.
    void initialize( TString geoCache, bool genHistograms ) {
        nEvents = 1; // only process single event

        // Create the forward system, cleaning up any previous instance
        if (FwdSystem::sInstance) {
            delete FwdSystem::sInstance;
        }
        FwdSystem::sInstance = new FwdSystem();

        // initialize the track fitter, cleaning up any previous instance
        if (mTrackFitter) {
            delete mTrackFitter;
        }
        mTrackFitter = new TrackFitter(mConfig, geoCache);
        mTrackFitter->setup();

        ForwardTrackMaker::initialize( geoCache, genHistograms );
    }

    void finish() {

        if (FwdSystem::sInstance){
            delete FwdSystem::sInstance;
            FwdSystem::sInstance = 0;
        }
        if (mTrackFitter){
            delete mTrackFitter;
            mTrackFitter= 0;
        }
    }
};

//________________________________________________________________________
StFwdTrackMaker::StFwdTrackMaker() : StMaker("fwdTrack"), mEventVertex(0,0,0), mForwardTracker(nullptr), mForwardData(nullptr), mGeoCache(""){
    LOG_DEBUG << "StFwdTrackMaker::StFwdTrackMaker()" << endm;
    mEventVertexCov.ResizeTo(3, 3);
    mEventVertexCov.Zero();
    
    SetAttr("useFtt",1);                 // Default Ftt on
    SetAttr("useFst",1);                 // Default Fst on
    SetAttr("useFcs",1);                 // Default Fcs on
    SetAttr("useEpd",1);                 // Default Epd on
    SetAttr("config", "config.xml");     // Default configuration file (user may override before Init())
    SetAttr("fillEvent",1); // fill StEvent
    SetAttr("fillAlignment",0); // Optional alignment diagnostics, off by default

    // Load the default configuration
    configLoaded = false;
    LoadConfiguration();

    // set additional default configuration values
    setOutputFilename( "stfwdtrackmaker_data.root" );
    LOG_DEBUG << "Done with StFwdTrackMaker::StFwdTrackMaker()" << endm;  
};

int StFwdTrackMaker::Finish() {
    if (mAlignmentFile) {
        mAlignmentFile->cd();
        if (mAlignmentTree) {
            mAlignmentTree->Write();
        }
        mAlignmentFile->Close();
        delete mAlignmentFile;
        mAlignmentFile = nullptr;
        mAlignmentTree = nullptr;
    }

    mForwardTracker->finish();
    return kStOk;
}

void StFwdTrackMaker::LoadConfiguration() {
    if (mConfigFile.length() < 5){
        // no config file specified, use default
        // 5 characters is the minimum length for a valid filename since we must have at least .xml
        mFwdConfig.load( defaultConfig, true );
        LOG_DEBUG << "Forward Tracker is using the default config" <<  mConfigFile << endm;
    } else {
        LOG_DEBUG << "Forward Tracker is using config from file : " <<  mConfigFile << endm;
        mFwdConfig.load( mConfigFile );
    }
    configLoaded = true;
}

//________________________________________________________________________
int StFwdTrackMaker::Init() {
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
    } else {
        LOG_INFO << "Using cached geometry file: " << mGeoCache << endm;
    }
    
    mForwardTracker = std::shared_ptr<ForwardTracker>(new ForwardTracker( ));
    mForwardTracker->setConfig(mFwdConfig);

    // in production we disable crit saving.
    mForwardTracker->setSaveCriteriaValues(false);

    mForwardData = std::shared_ptr<FwdDataSource>(new FwdDataSource());
    mForwardTracker->setData(mForwardData);
    mForwardTracker->initialize( mGeoCache, false );

    if ( IAttr("fillAlignment") ) {
        mAlignmentFile = new TFile(mAlignmentOutputFilename.c_str(), "RECREATE");
        mAlignmentTree = new TTree("fwdAlign", "Forward alignment diagnostics");
        mAlignmentTree->Branch("run", &mAlignRun, "run/I");
        mAlignmentTree->Branch("event", &mAlignEvent, "event/I");
        mAlignmentTree->Branch("trackIndex", &mAlignTrackIndex, "trackIndex/I");
        mAlignmentTree->Branch("pointIndex", &mAlignPointIndex, "pointIndex/I");
        mAlignmentTree->Branch("measurementIndex", &mAlignMeasurementIndex, "measurementIndex/I");
        mAlignmentTree->Branch("detId", &mAlignDetId, "detId/I");
        mAlignmentTree->Branch("hitId", &mAlignHitId, "hitId/I");
        mAlignmentTree->Branch("fstGlobalSensor", &mAlignFstGlobalSensor, "fstGlobalSensor/I");
        mAlignmentTree->Branch("fstDisk", &mAlignFstDisk, "fstDisk/I");
        mAlignmentTree->Branch("fstWedge", &mAlignFstWedge, "fstWedge/I");
        mAlignmentTree->Branch("fstSensor", &mAlignFstSensor, "fstSensor/I");
        mAlignmentTree->Branch("measurementDim", &mAlignMeasurementDim, "measurementDim/I");
        mAlignmentTree->Branch("residualDim", &mAlignResidualDim, "residualDim/I");
        mAlignmentTree->Branch("hasResidual", &mAlignHasResidual, "hasResidual/I");
        mAlignmentTree->Branch("nSeeds", &mAlignNSeeds, "nSeeds/I");
        mAlignmentTree->Branch("nFitTracks", &mAlignNFitTracks, "nFitTracks/I");
        mAlignmentTree->Branch("chi2", &mAlignChi2, "chi2/F");
        mAlignmentTree->Branch("ndf", &mAlignNdf, "ndf/I");
        mAlignmentTree->Branch("pval", &mAlignPval, "pval/F");
        mAlignmentTree->Branch("fitConverged", &mAlignFitConverged, "fitConverged/I");
        mAlignmentTree->Branch("fitConvergedFully", &mAlignFitConvergedFully, "fitConvergedFully/I");
        mAlignmentTree->Branch("fitConvergedPartially", &mAlignFitConvergedPartially, "fitConvergedPartially/I");
        mAlignmentTree->Branch("trackType", &mAlignTrackType, "trackType/I");
        mAlignmentTree->Branch("trackNHitsFit", &mAlignTrackNHitsFit, "trackNHitsFit/I");
        mAlignmentTree->Branch("trackNFstHits", &mAlignTrackNFstHits, "trackNFstHits/I");
        mAlignmentTree->Branch("trackPx", &mAlignTrackPx, "trackPx/F");
        mAlignmentTree->Branch("trackPy", &mAlignTrackPy, "trackPy/F");
        mAlignmentTree->Branch("trackPz", &mAlignTrackPz, "trackPz/F");
        mAlignmentTree->Branch("trackP", &mAlignTrackP, "trackP/F");
        mAlignmentTree->Branch("trackPt", &mAlignTrackPt, "trackPt/F");
        mAlignmentTree->Branch("trackEta", &mAlignTrackEta, "trackEta/F");
        mAlignmentTree->Branch("sorting", &mAlignSorting, "sorting/F");
        mAlignmentTree->Branch("meas0", &mAlignMeas0, "meas0/F");
        mAlignmentTree->Branch("meas1", &mAlignMeas1, "meas1/F");
        mAlignmentTree->Branch("meas2", &mAlignMeas2, "meas2/F");
        mAlignmentTree->Branch("trackPred0", &mAlignTrackPred0, "trackPred0/F");
        mAlignmentTree->Branch("trackPred1", &mAlignTrackPred1, "trackPred1/F");
        mAlignmentTree->Branch("trackPred2", &mAlignTrackPred2, "trackPred2/F");
        mAlignmentTree->Branch("fstRawR", &mAlignFstRawR, "fstRawR/F");
        mAlignmentTree->Branch("fstRawStripPhi", &mAlignFstRawStripPhi, "fstRawStripPhi/F");
        mAlignmentTree->Branch("fstMeanPhiStrip", &mAlignFstMeanPhiStrip, "fstMeanPhiStrip/F");
        mAlignmentTree->Branch("fstHitGlobalX", &mAlignFstHitGlobalX, "fstHitGlobalX/F");
        mAlignmentTree->Branch("fstHitGlobalY", &mAlignFstHitGlobalY, "fstHitGlobalY/F");
        mAlignmentTree->Branch("fstHitGlobalZ", &mAlignFstHitGlobalZ, "fstHitGlobalZ/F");
        mAlignmentTree->Branch("fstPlaneOriginX", &mAlignFstPlaneOriginX, "fstPlaneOriginX/F");
        mAlignmentTree->Branch("fstPlaneOriginY", &mAlignFstPlaneOriginY, "fstPlaneOriginY/F");
        mAlignmentTree->Branch("fstPlaneOriginZ", &mAlignFstPlaneOriginZ, "fstPlaneOriginZ/F");
        mAlignmentTree->Branch("fstPlaneUX", &mAlignFstPlaneUX, "fstPlaneUX/F");
        mAlignmentTree->Branch("fstPlaneUY", &mAlignFstPlaneUY, "fstPlaneUY/F");
        mAlignmentTree->Branch("fstPlaneUZ", &mAlignFstPlaneUZ, "fstPlaneUZ/F");
        mAlignmentTree->Branch("fstPlaneVX", &mAlignFstPlaneVX, "fstPlaneVX/F");
        mAlignmentTree->Branch("fstPlaneVY", &mAlignFstPlaneVY, "fstPlaneVY/F");
        mAlignmentTree->Branch("fstPlaneVZ", &mAlignFstPlaneVZ, "fstPlaneVZ/F");
        mAlignmentTree->Branch("fstMeasGlobalX", &mAlignFstMeasGlobalX, "fstMeasGlobalX/F");
        mAlignmentTree->Branch("fstMeasGlobalY", &mAlignFstMeasGlobalY, "fstMeasGlobalY/F");
        mAlignmentTree->Branch("fstMeasGlobalZ", &mAlignFstMeasGlobalZ, "fstMeasGlobalZ/F");
        mAlignmentTree->Branch("fstClosureX", &mAlignFstClosureX, "fstClosureX/F");
        mAlignmentTree->Branch("fstClosureY", &mAlignFstClosureY, "fstClosureY/F");
        mAlignmentTree->Branch("fstClosureZ", &mAlignFstClosureZ, "fstClosureZ/F");
        mAlignmentTree->Branch("fstClosureU", &mAlignFstClosureU, "fstClosureU/F");
        mAlignmentTree->Branch("fstClosureV", &mAlignFstClosureV, "fstClosureV/F");
        mAlignmentTree->Branch("fstClosureMag", &mAlignFstClosureMag, "fstClosureMag/F");
        mAlignmentTree->Branch("resBiased0", &mAlignResBiased0, "resBiased0/F");
        mAlignmentTree->Branch("resBiased1", &mAlignResBiased1, "resBiased1/F");
        mAlignmentTree->Branch("resBiased2", &mAlignResBiased2, "resBiased2/F");
        mAlignmentTree->Branch("resBiasedSigma0", &mAlignResBiasedSigma0, "resBiasedSigma0/F");
        mAlignmentTree->Branch("resBiasedSigma1", &mAlignResBiasedSigma1, "resBiasedSigma1/F");
        mAlignmentTree->Branch("resBiasedSigma2", &mAlignResBiasedSigma2, "resBiasedSigma2/F");
        mAlignmentTree->Branch("pullBiased0", &mAlignPullBiased0, "pullBiased0/F");
        mAlignmentTree->Branch("pullBiased1", &mAlignPullBiased1, "pullBiased1/F");
        mAlignmentTree->Branch("pullBiased2", &mAlignPullBiased2, "pullBiased2/F");
        mAlignmentTree->Branch("resUnbiased0", &mAlignResUnbiased0, "resUnbiased0/F");
        mAlignmentTree->Branch("resUnbiased1", &mAlignResUnbiased1, "resUnbiased1/F");
        mAlignmentTree->Branch("resUnbiased2", &mAlignResUnbiased2, "resUnbiased2/F");
        mAlignmentTree->Branch("resUnbiasedSigma0", &mAlignResUnbiasedSigma0, "resUnbiasedSigma0/F");
        mAlignmentTree->Branch("resUnbiasedSigma1", &mAlignResUnbiasedSigma1, "resUnbiasedSigma1/F");
        mAlignmentTree->Branch("resUnbiasedSigma2", &mAlignResUnbiasedSigma2, "resUnbiasedSigma2/F");
        mAlignmentTree->Branch("pullUnbiased0", &mAlignPullUnbiased0, "pullUnbiased0/F");
        mAlignmentTree->Branch("pullUnbiased1", &mAlignPullUnbiased1, "pullUnbiased1/F");
        mAlignmentTree->Branch("pullUnbiased2", &mAlignPullUnbiased2, "pullUnbiased2/F");
    }

    // Setup the mFwdHitLoader


    // geometry should be available from here (mForwardTracker will initialize cache if needed)
    if (gGeoManager) {
        FwdGeomUtils fwdGeoUtils( gGeoManager );
        // get the z-locations from geometry model and fallback to the defaults
        auto fstZ = fwdGeoUtils.fstZ( {151.750000, 165.248001, 178.781006} );
        mFstZFromGeom.assign( fstZ.begin(), fstZ.end() );
        auto fttZ = fwdGeoUtils.fttZ( {280.904999, 303.704987, 326.605011, 349.404999} );
        mFttZFromGeom.assign( fttZ.begin(), fttZ.end() );
    }
    return kStOK;
};

EventStats StFwdTrackMaker::GetEventStats() { 
    return mForwardTracker->getEventStats(); 
}

/**
 * Loads the Monte Carlo (MC) tracks from the GEANT simulation data.
 *
 * @param mcTrackMap A reference to the MC track map.
 *
 * @return The number of forward tracks.
 *
 * @throws None.
 */
size_t StFwdTrackMaker::loadMcTracks( FwdDataSource::McTrackMap_t &mcTrackMap ){

    LOG_DEBUG << "Looking for GEANT sim vertex info" << endm;
    St_g2t_vertex *g2t_vertex = (St_g2t_vertex *)GetDataSet("geant/g2t_vertex");

    if ( g2t_vertex != nullptr ) {
        // Set the MC Vertex for track fitting
        g2t_vertex_st *vert = (g2t_vertex_st*)g2t_vertex->At(0);
        TMatrixDSym cov;
        cov.ResizeTo(3, 3);
        cov(0, 0) = 0.001;
        cov(1, 1) = 0.001;
        cov(2, 2) = 0.001;
        mForwardTracker->setEventVertex( TVector3( vert->ge_x[0], vert->ge_x[1], vert->ge_x[2] ), cov );
    }
    // Get geant tracks
    St_g2t_track *g2t_track = (St_g2t_track *)GetDataSet("geant/g2t_track");

    if (!g2t_track)
        return 0;

    LOG_DEBUG << g2t_track->GetNRows() << " mc tracks in geant/g2t_track " << endm;

    for (int irow = 0; irow < g2t_track->GetNRows(); irow++) {
        g2t_track_st *track = (g2t_track_st *)g2t_track->At(irow);

        if (0 == track)
            continue;

        int track_id = track->id;
        TVector3 pp( track->p[0], track->p[1], track->p[2] );
        int q = track->charge;
        // sometimes the track->eta is wrong, pt, phi
        if (!mcTrackMap[track_id] )
            mcTrackMap[track_id] = shared_ptr<McTrack>(new McTrack(pp.Pt(), pp.Eta(), pp.Phi(), q, track->start_vertex_p));

    } // loop on track (irow)

    // now check the Mc tracks against the McEvent filter
    size_t nForwardTracks = 0;
    size_t nForwardTracksNoThreshold = 0;
    for (auto mctm : mcTrackMap ){
        if ( mctm.second == nullptr ) continue;
        if ( mctm.second->mEta > 2.5 && mctm.second->mEta < 4.0 ){
            nForwardTracksNoThreshold++;
            if ( mctm.second->mPt > 0.05  )
                nForwardTracks++;
        }
    } // loop on mcTrackMap
    return nForwardTracks;
} // loadMcTracks

TVector3 StFwdTrackMaker::GetEventPrimaryVertex(){
    if ( mFwdVertexSource != kFwdVertexSourceUnknown ){
        // This includes the case where we have already searched and found nothing
        return mEventVertex;
    }

    mEventVertexCov.ResizeTo(3, 3);
    mEventVertexCov.Zero();
    double sig2 = 1;// default resolution, overwritten if valid vtx found
    mEventVertexCov(0, 0) = sig2; 
    mEventVertexCov(1, 1) = sig2;
    mEventVertexCov(2, 2) = sig2;
    // if something is found it will overwrite this, if not
    // it will indicate that we have searched and found nothing
    mFwdVertexSource = kFwdVertexSourceNone;

    /*****************************************************
     * Add Primary Vertex to the track
     */
    St_g2t_vertex *g2t_vertex = (St_g2t_vertex *)GetDataSet("geant/g2t_vertex");
    LOG_DEBUG << "Searching for Event Vertex from geant/g2t_vertex: " << g2t_vertex << endm;
    if ( g2t_vertex != nullptr ) {
        // Set the MC Vertex for track fitting
        g2t_vertex_st *vert = (g2t_vertex_st*)g2t_vertex->At(0);
        LOG_INFO << "Setting Event Vertex from geant/g2t_vertex[0]: " << vert << endm;
        if ( vert ){
            mEventVertexCov.ResizeTo(3, 3);
            const double sigXY = 0.1; // TODO: read from MC vertex info?
            const double sigZ = 0.1;
            mEventVertexCov(0, 0) = pow(sigXY,2);
            mEventVertexCov(1, 1) = pow(sigXY,2);
            mEventVertexCov(2, 2) = pow(sigZ, 2);
            auto rhc = TVectorD( 3 );
            rhc[0] = vert->ge_x[0];
            rhc[1] = vert->ge_x[1];
            rhc[2] = vert->ge_x[2];
            mEventVertex.SetXYZ( vert->ge_x[0], vert->ge_x[1], vert->ge_x[2] );
            mFwdVertexSource = kFwdVertexSourceMc;
            return mEventVertex;
        }
    }

    // or try the McEvent
    StMcEvent *stMcEvent = static_cast<StMcEvent *>(GetInputDS("StMcEvent"));
    LOG_DEBUG << "Searching for Event Vertex from StMcEvent: " << stMcEvent << endm;
    if (stMcEvent && stMcEvent->primaryVertex() ) {
        StThreeVectorF vertex = stMcEvent->primaryVertex()->position();
        mEventVertex.SetXYZ( vertex.x(), vertex.y(), vertex.z() );
        mFwdVertexSource = kFwdVertexSourceMc;
        LOG_INFO << "FWD Tracking on event with MC Primary Vertex: " << mEventVertex.X() << ", " << mEventVertex.Y() << ", " << mEventVertex.Z() << endm;

        const double sigXY = 0.1; // TODO: read from MC vertex info?
        const double sigZ = 0.1;
        mEventVertexCov(0, 0) = pow(sigXY,2);
        mEventVertexCov(1, 1) = pow(sigXY,2);
        mEventVertexCov(2, 2) = pow(sigZ, 2);
        return mEventVertex;
    }

    StMuDstMaker *mMuDstMaker = (StMuDstMaker *)GetMaker("MuDst");
    if(mMuDstMaker && mMuDstMaker->muDst() && mMuDstMaker->muDst()->numberOfPrimaryVertices() > 0 && mMuDstMaker->muDst()->primaryVertex() ) {
        mEventVertex.SetX(mMuDstMaker->muDst()->primaryVertex()->position().x());
        mEventVertex.SetY(mMuDstMaker->muDst()->primaryVertex()->position().y());
        mEventVertex.SetZ(mMuDstMaker->muDst()->primaryVertex()->position().z());
        mFwdVertexSource = kFwdVertexSourceTpc;
        return mEventVertex;
    } 
    
    
    LOG_DEBUG << "FWD Tracking on event without available Mu Primary Vertex" << endm;
    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent) return mEventVertex; // if we get here and there is no StEvent, we are done
    
    StBTofCollection *btofC = stEvent->btofCollection();
    if (!btofC) {
        LOG_WARN << "Cannot get BTOF collections, Cannot use VPD vertex" << endm;
        return mEventVertex;
    }

    StBTofHeader * btofHeader = btofC->tofHeader();
    if (!btofHeader){
        LOG_WARN << "Cannot get BTOF Header, Cannot use VPD vertex" << endm;
        return mEventVertex;
    }

    if ( btofHeader->vpdVz() && fabs(btofHeader->vpdVz()) < 100 ){
        // default event vertex
        LOG_DEBUG << "FWD Tracking on event using VPD z vertex: (, 0, 0, " << btofHeader->vpdVz() << " )" << endm;
        mFwdVertexSource = kFwdVertexSourceVpd;
        const double sigXY = 1;
        const double sigZ = 6; // approximate resolution of VPD in p+p collisions
        mEventVertexCov(0, 0) = pow(sigXY,2);
        mEventVertexCov(1, 1) = pow(sigXY,2);
        mEventVertexCov(2, 2) = pow(sigZ, 2);
        mEventVertex.SetXYZ( 0, 0, btofHeader->vpdVz() );
        return mEventVertex;
    }
    
    // if we get here we failed to find a valid vtx
    return mEventVertex;
}

//________________________________________________________________________
int StFwdTrackMaker::Make() {
    // return kStOk;
    // START time for measuring tracking
    long long itStart = FwdTrackerUtils::nowNanoSecond();

    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    if (!stEvent) return kStOk;

    

    /**********************************************************************/
    // Access forward track and hit maps
    FwdDataSource::McTrackMap_t &mcTrackMap = mForwardData->getMcTracks();
    FwdDataSource::HitMap_t &hitMap         = mForwardData->getFttHits();
    FwdDataSource::HitMap_t &fsiHitMap      = mForwardData->getFstHits();
    FwdDataSource::HitMap_t &epdHitMap      = mForwardData->getEpdHits();

    mFwdHitLoader.setStEvent( stEvent );
    mFwdHitLoader.setMuDstMaker( (StMuDstMaker *)GetMaker("MuDst") );
    mFwdHitLoader.setTables(
        (St_g2t_fts_hit *)GetDataSet("geant/g2t_stg_hit"),
        (St_g2t_fts_hit *)GetDataSet("geant/g2t_fsi_hit"),
        nullptr
    );
    mFcsDb = dynamic_cast<StFcsDb*>(GetDataSet("fcsDb"));
    if ( !mFcsDb ) {
        LOG_WARN << "No FCS DB found, cannot load FCS hits" << endm;
    }

    /**********************************************************************/
    // get the primary vertex for use with FWD tracking
    mFwdVertexSource = StFwdTrackMaker::kFwdVertexSourceUnknown;
    GetEventPrimaryVertex();
    LOG_DEBUG << "FWD Vertex Source: " << mFwdVertexSource << endm;
    LOG_DEBUG << "Setting FWD event vertex to: " << TString::Format("mEventVertex=(%0.3f+/-%0.3f, %0.3f+/-%0.3f, %0.3f+/-%0.3f)", mEventVertex.X(), sqrt(mEventVertexCov(0, 0)), mEventVertex.Y(), sqrt(mEventVertexCov(1, 1)), mEventVertex.Z(), sqrt(mEventVertexCov(2, 2)) ) << endm;
    mForwardTracker->setEventVertex( mEventVertex, mEventVertexCov );

    
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
        mFwdHitLoader.loadFttHits( mcTrackMap, hitMap );
    }

    /**********************************************************************/
    // Load FST
    if ( IAttr("useFst") ) {
        LOG_DEBUG << ">>StFwdTrackMaker::loadFstHits" << endm;
        int fstCount = mFwdHitLoader.loadFstHits( mcTrackMap, fsiHitMap );
        LOG_DEBUG << "Loaded " << fstCount << " FST hits" << endm;
    }

    /**********************************************************************/
    // Load FCS
    LOG_DEBUG << ">>StFwdTrackMaker::loadFcsHits" << endm;
    if ( IAttr("useEpd") ) {
        LOG_DEBUG << ">>StFwdTrackMaker::loadEpdHits" << endm;
        int epdCount = mFwdHitLoader.loadEpdHits( mcTrackMap, epdHitMap, mFcsDb );
        LOG_DEBUG << "Loaded " << epdCount << " Epd hits" << endm;
    }

    /**********************************************************************/
    // Print out the MC tracks and their hit counts
    map<int, int> nFstMcTracks;
    map<int, int> nFttMcTracks;
    for ( auto kv : mcTrackMap ){
        if ( kv.second == nullptr ) continue;
        LOG_DEBUG << "MC Track: id=" << kv.first << ", nFTT=" << kv.second->mFttHits.size() << ", nFST=" << kv.second->mFstHits.size() << endm;
        nFstMcTracks[ kv.second->mFstHits.size() ]++;
        nFttMcTracks[ kv.second->mFttHits.size() ]++;
    }
    int idealNumberOfSeeds = (nFstMcTracks[3]);
    if (!mcTrackMap.empty()) {
        // Only print this if we have MC tracks
        LOG_INFO << "There are: " << Form( "%d with 0 FST, %d with 1 FST, %d with 2 FST, %d with 3 FST", nFstMcTracks[0], nFstMcTracks[1], nFstMcTracks[2], nFstMcTracks[3] ) << endm;
        LOG_INFO << "There are: " << Form( "%d with 0 FTT, %d with 1 FTT, %d with 2 FTT, %d with 3 FTT, %d with 4 FTT, %d with 5 FTT, %d with 6 FTT, %d with 7 FTT, %d with 8 FTT", nFttMcTracks[0], nFttMcTracks[1], nFttMcTracks[2], nFttMcTracks[3], nFttMcTracks[4], nFttMcTracks[5], nFttMcTracks[6], nFttMcTracks[7], nFttMcTracks[8] ) << endm;
        LOG_INFO << "There are " << Form( "%d McTracks with > 2 FST hits (#of possible seeds)", idealNumberOfSeeds  ) << endm;
    }

    /**********************************************************************/
    // Run Track finding + fitting
    LOG_DEBUG << ">>START Event Forward Tracking" << endm;
    LOG_INFO << "\tFinding FWD Track Seeds" << endm;
    mForwardTracker->findTrackSeeds();

    
    // Report the results of the seed finding, in the future we could provide more info about #hits etc.
    LOG_INFO << "<<Fwd Tracking Found : " << mForwardTracker -> getTrackSeeds().size() << " Track Seeds from " << fsiHitMap.size() << " FST hits and " << hitMap.size() << " sTGC hits"  << endm;
    if ( idealNumberOfSeeds > 0 ){
        float seedFindingEff = ( mForwardTracker -> getTrackSeeds().size() + 1e-5 ) / ( idealNumberOfSeeds + 1e-5 );
        LOG_INFO << "    (vs. " << idealNumberOfSeeds << " McTracks with FST>2, eff = " << seedFindingEff << ")" << endm;
    } 
    /**********************************************************************/
    
    /**********************************************************************/
    // Run Track fitting on the seeds we found
    LOG_INFO << "\tFitting FWD Track Seeds" << endm;
    mForwardTracker->doTrackFitting( mForwardTracker->getTrackSeeds() );
    if ( IAttr("fillAlignment") ) {
        FillAlignment();
    }
    LOG_INFO << "<<Fwd Tracking Fit :" << mForwardTracker -> getTrackResults().size() << " GenFit Tracks" << endm;
    LOG_DEBUG << "<<FINISH Event Forward Tracking" << endm;
    /**********************************************************************/


    LOG_DEBUG << "Forward tracking on this event took " << (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6 << " ms" << endm;
    if ( IAttr("fillEvent") ) {
        if (!stEvent) {
            LOG_WARN << "No StEvent found. Forward tracks will not be saved" << endm;
            return kStWarn;
        }
        FillEvent();
    } // IAttr FillEvent

    return kStOK;
} // Make

/**
 * Creates a StFwdTrack object from a GenfitTrackResult.
 *
 * @param gtr The GenfitTrackResult object containing the track information.
 * @param indexTrack The index of the track.
 *
 * @return A pointer to the created StFwdTrack object, or nullptr if the GenfitTrackResult is nullptr.
 *
 * @throws None.
 */
StFwdTrack * StFwdTrackMaker::makeStFwdTrack( GenfitTrackResult &gtr, size_t indexTrack ){
    LOG_DEBUG << "StFwdTrackMaker::makeStFwdTrack()" << endm;
    StFwdTrack *fwdTrack = new StFwdTrack( );
    TVector3 p = gtr.mMomentum;

    /*******************************************************************************/
    // store the seed points for the track
    int nSeedPoints = 0;
    float cov[9]; // reused covariance matrix for seed points
    for ( auto s : gtr.mSeed ){
        FwdHit * fh = static_cast<FwdHit*>( s );
        if (!fh) continue;
        cov[0] = fh->_covmat(0,0); cov[3] = fh->_covmat(1,0); cov[6] = fh->_covmat(2,0);
        cov[1] = fh->_covmat(0,1); cov[4] = fh->_covmat(1,1); cov[7] = fh->_covmat(2,1);
        cov[2] = fh->_covmat(0,2); cov[5] = fh->_covmat(1,2); cov[8] = fh->_covmat(2,2);

        StFwdTrackSeedPoint p(
            StThreeVectorD( fh->getX(), fh->getY(), fh->getZ() ),
            fh->_detid * 10 + fh->getSector(), // 10 * detid + sector
            fh->getTrackId(),
            cov
        );
        if ( fh->isFst() )
            fwdTrack->mFSTPoints.push_back( p );
        else if ( fh->isFtt() )
            fwdTrack->mFTTPoints.push_back( p );

        nSeedPoints++;
    }

    // set total number of seed points
    fwdTrack->setNumberOfSeedPoints( nSeedPoints ); 
    int idt = 0;
    double qual = 0;
    idt = MCTruthUtils::dominantContribution(gtr.mSeed, qual);
    fwdTrack->setMc( idt, qual*100 ); // QAtruth stored as UChar_t
    LOG_DEBUG << "Dominant contribution: " << idt << " with quality " << qual << endm;


    // for seed only, we save the seed charge and momentum computed from the seed points
    fwdTrack->setCharge( gtr.mCharge );
    fwdTrack->setPrimaryMomentum( StThreeVectorD( gtr.mMomentum.X(), gtr.mMomentum.Y(), gtr.mMomentum.Z() ) );
    fwdTrack->setVtxIndexAndTrackType( gtr.mVertexIndex, gtr.mTrackType );
    fwdTrack->setGlobalTrackIndex( gtr.mGlobalTrackIndex);

    // Fit failed beyond use
    if ( !gtr.mIsFitConvergedPartially|| gtr.mNumFitPoints == 0 ){
        // if num points == 0 then calling PVal seg faults :/
        fwdTrack->setDidFitConverge( false );
        fwdTrack->setDidFitConvergeFully( false );
        fwdTrack->setNumberOfFailedPoints( 99 );
        fwdTrack->setNumberOfFitPoints( 0 );
        fwdTrack->setChi2( 0 );
        fwdTrack->setNDF( 0 );
        fwdTrack->setPval( 0 );

        fwdTrack->setNumberOfFitPoints( 1 ); // setting this to 1 so that the charge is still saved as charge * n
        
        gtr.Clear();
        return fwdTrack;
    }
    // Fill fit quality info
    fwdTrack->setDidFitConverge( gtr.mIsFitConverged );
    fwdTrack->setDidFitConvergeFully( gtr.mIsFitConvergedFully );
    fwdTrack->setNumberOfFailedPoints( gtr.mNFailedPoints );

    fwdTrack->setNumberOfFitPoints( gtr.mNumFitPoints );
    fwdTrack->setChi2( gtr.mChi2 );
    fwdTrack->setNDF( gtr.mNdf );
    fwdTrack->setPval( gtr.mPval );

    // DCA and fitted momentum
    fwdTrack->setDCA( gtr.mDCA.X(), gtr.mDCA.Y(), gtr.mDCA.Z() );
    fwdTrack->setPrimaryMomentum( StThreeVectorD( gtr.mMomentum.X(), gtr.mMomentum.Y(), gtr.mMomentum.Z() ) );

    /*******************************************************************************/
    // if the track did not converged, do not try to project it
    if ( !gtr.mIsFitConvergedFully ){
        gtr.Clear();
        LOG_WARN << "Genfit track did not converge fully, skipping projections" << endm;

        return fwdTrack;
    }

    /*******************************************************************************/
    // compute projections to z-planes of various detectors
    // TODO: update FCS to use correct z + angle
    // Use vector<pair> instead of map to allow multiple entries per detector
    std::vector<std::pair<int, float>> detectorZPlanes;

    // Add TPC projection
    detectorZPlanes.push_back({ kTpcId, 0.0 });

    // Add FST projections (check vector size first)
    for (size_t i = 0; i < mFstZFromGeom.size() && i < 3; i++) {
        detectorZPlanes.push_back({ kFstId, mFstZFromGeom[i] });
    }

    // Add FTT projections (check vector size first)
    for (size_t i = 0; i < mFttZFromGeom.size() && i < 4; i++) {
        detectorZPlanes.push_back({ kFttId, mFttZFromGeom[i] });
    }

    // Add FCS projections
    detectorZPlanes.push_back({ kFcsPresId, 375.0 });
    detectorZPlanes.push_back({ kFcsWcalId, 715.0 });
    detectorZPlanes.push_back({ kFcsHcalId, 807.0 });

    size_t zIndex = 0;
    TVector3 mom(0, 0, 0);
    TVector3 tv3(0, 0, 0);
    for ( auto zp : detectorZPlanes ){
        int detIndex = zp.first;
        float z = zp.second;
        tv3.SetXYZ(0, 0, 0);
        std::fill(std::begin(cov), std::end(cov), 0.0f);
        LOG_DEBUG << "Projecting to: " << detIndex << " at z=" << z << endm;
        if ( detIndex != kFcsHcalId && detIndex != kFcsWcalId ){
            float detpos[3] = {0,0,z};
            float detnorm[3] = {0,0,1};
            tv3 = ObjExporter::trackPosition( gtr.mTrack.get(), detpos, detnorm, cov, mom );
        } else {
            // use a straight line projection to HCAL since GenFit cannot handle long projections
            int det=0;
            if( detIndex==kFcsWcalId ){
                det = 0;   // North side for negative px
                // South side for positive px, since px==0 does not hit detector choose south side for that case
                if( p[2]>=0 && p[0]>=0 ){ det=1; }
                if( p[2]<0  && p[0]<0  ){ det=1; }
            }
            //Since detIndex cannot be both don't need "else if"
            if( detIndex==kFcsHcalId ){
                det = 2;  // North side for negative px
                // South side for positive px, since px==0 does not hit detector choose south side for that case
                if( p[2]>=0 && p[0]>=0 ){ det=3; }
                if( p[2]<0  && p[0]<0  ){ det=3; }
            }
            if (!mFcsDb) {
                LOG_ERROR << "FCS database not initialized, cannot project to FCS" << endm;
                continue;
            }
            StThreeVectorD xyzoff = mFcsDb->getDetectorOffset(det);
            StThreeVectorD planenormal = mFcsDb->getNormal(det);
            float xyz0[3] = { 0, 0, 575.0 };
            float xyz1[3] = { 0, 0, 625.0 };
            float xyzdet[3] = { (float)xyzoff.x(), (float)xyzoff.y(), (float)xyzoff.z() };
            float detnorm[3] = { (float)planenormal.x(), (float)planenormal.y(), (float)planenormal.z() };
            LOG_DEBUG << "Projecting to: " << detIndex << endm;
            tv3 = ObjExporter::projectAsStraightLine( gtr.mTrack.get(), xyz0, xyz1, xyzdet, detnorm, cov, mom );
        }
        fwdTrack->mProjections.push_back( StFwdTrackProjection( detIndex, StThreeVectorF( tv3.X(), tv3.Y(), tv3.Z() ), StThreeVectorF( mom.X(), mom.Y(), mom.Z() ), cov) );
        // LOG_INFO << "Projection added for " << detIndex << " at z=" << z << endm;
        zIndex++;
    }
    /*******************************************************************************/

    /*******************************************************************************/
    // clear the GenfitTrackResult
    gtr.Clear();

    // return the StFwdTrack we made
    return fwdTrack;
}

void StFwdTrackMaker::FillAlignment() {
    if (!mAlignmentTree)
        return;

    StEvent *stEvent = static_cast<StEvent *>(GetInputDS("StEvent"));
    mAlignRun = stEvent ? stEvent->runId() : 0;
    mAlignEvent = stEvent ? stEvent->id() : 0;
    mAlignNSeeds = static_cast<int>(mForwardTracker->getTrackSeeds().size());
    mAlignNFitTracks = static_cast<int>(mForwardTracker->getTrackResults().size());

    FwdGeomUtils *alignmentGeo = nullptr;
    if (gGeoManager)
        alignmentGeo = new FwdGeomUtils(gGeoManager);

    int trackIndex = 0;
    for ( const auto &gtr : mForwardTracker->getTrackResults() ) {
        mAlignTrackIndex = trackIndex++;
        mAlignChi2 = gtr.mChi2;
        mAlignNdf = gtr.mNdf;
        mAlignPval = gtr.mPval;
        mAlignFitConverged = gtr.mIsFitConverged ? 1 : 0;
        mAlignFitConvergedFully = gtr.mIsFitConvergedFully ? 1 : 0;
        mAlignFitConvergedPartially = gtr.mIsFitConvergedPartially ? 1 : 0;
        mAlignTrackType = gtr.mTrackType;
        mAlignTrackNHitsFit = gtr.mNumFitPoints;
        mAlignTrackNFstHits = 0;
        setAlignmentTrackKinematics(
            gtr.mMomentum,
            mAlignTrackPx, mAlignTrackPy, mAlignTrackPz,
            mAlignTrackP, mAlignTrackPt,
            mAlignTrackEta
        );

        if (!gtr.mTrack)
            continue;

        auto rep = gtr.mTrack->getCardinalRep();
        for ( auto point : gtr.mTrack->getPoints() ) {
            if (!point)
                continue;
            const unsigned int nRawMeasurements = point->getNumRawMeasurements();
            for ( unsigned int iMeas = 0; iMeas < nRawMeasurements; ++iMeas ) {
                auto rawMeasurement = point->getRawMeasurement(iMeas);
                if (rawMeasurement && rawMeasurement->getDetId() == kFstId)
                    ++mAlignTrackNFstHits;
            }
        }

        int pointIndex = 0;
        for ( auto point : gtr.mTrack->getPoints() ) {
            mAlignPointIndex = pointIndex++;
            if (!point || point->getNumRawMeasurements() == 0)
                continue;

            genfit::KalmanFitterInfo *kfi = nullptr;
            if (rep && point->hasFitterInfo(rep)) {
                kfi = dynamic_cast<genfit::KalmanFitterInfo*>(point->getFitterInfo(rep));
            }

            const unsigned int nRawMeasurements = point->getNumRawMeasurements();
            for ( unsigned int iMeas = 0; iMeas < nRawMeasurements; ++iMeas ) {
                auto rawMeasurement = point->getRawMeasurement(iMeas);
                if (!rawMeasurement)
                    continue;

                mAlignMeasurementIndex = static_cast<int>(iMeas);
                mAlignDetId = rawMeasurement->getDetId();
                mAlignHitId = rawMeasurement->getHitId();
                mAlignMeasurementDim = static_cast<int>(rawMeasurement->getDim());
                mAlignResidualDim = 0;
                mAlignHasResidual = 0;
                const double sortingParameter = point->getSortingParameter();
                mAlignSorting = sortingParameter;
                mAlignFstGlobalSensor = -1;
                mAlignFstDisk = -1;
                mAlignFstWedge = -1;
                mAlignFstSensor = -1;
                if (mAlignDetId == kFstId) {
                    int globalSensor = static_cast<int>(sortingParameter) - 1;
                    if (globalSensor >= 0 && globalSensor < kFstNumSensors) {
                        mAlignFstGlobalSensor = globalSensor;
                        FwdHit::fstSensorWedgeDiskFromGlobalIndex(globalSensor, mAlignFstDisk, mAlignFstWedge, mAlignFstSensor);
                    }
                }
                mAlignResBiased0 = kInvalidAlignValue;
                mAlignResBiased1 = kInvalidAlignValue;
                mAlignResBiased2 = kInvalidAlignValue;
                mAlignResBiasedSigma0 = kInvalidAlignValue;
                mAlignResBiasedSigma1 = kInvalidAlignValue;
                mAlignResBiasedSigma2 = kInvalidAlignValue;
                mAlignPullBiased0 = kInvalidAlignValue;
                mAlignPullBiased1 = kInvalidAlignValue;
                mAlignPullBiased2 = kInvalidAlignValue;
                mAlignResUnbiased0 = kInvalidAlignValue;
                mAlignResUnbiased1 = kInvalidAlignValue;
                mAlignResUnbiased2 = kInvalidAlignValue;
                mAlignResUnbiasedSigma0 = kInvalidAlignValue;
                mAlignResUnbiasedSigma1 = kInvalidAlignValue;
                mAlignResUnbiasedSigma2 = kInvalidAlignValue;
                mAlignPullUnbiased0 = kInvalidAlignValue;
                mAlignPullUnbiased1 = kInvalidAlignValue;
                mAlignPullUnbiased2 = kInvalidAlignValue;
                mAlignTrackPred0 = kInvalidAlignValue;
                mAlignTrackPred1 = kInvalidAlignValue;
                mAlignTrackPred2 = kInvalidAlignValue;
                mAlignFstRawR = kInvalidAlignValue;
                mAlignFstRawStripPhi = kInvalidAlignValue;
                mAlignFstMeanPhiStrip = kInvalidAlignValue;
                mAlignFstHitGlobalX = kInvalidAlignValue;
                mAlignFstHitGlobalY = kInvalidAlignValue;
                mAlignFstHitGlobalZ = kInvalidAlignValue;
                mAlignFstPlaneOriginX = kInvalidAlignValue;
                mAlignFstPlaneOriginY = kInvalidAlignValue;
                mAlignFstPlaneOriginZ = kInvalidAlignValue;
                mAlignFstPlaneUX = kInvalidAlignValue;
                mAlignFstPlaneUY = kInvalidAlignValue;
                mAlignFstPlaneUZ = kInvalidAlignValue;
                mAlignFstPlaneVX = kInvalidAlignValue;
                mAlignFstPlaneVY = kInvalidAlignValue;
                mAlignFstPlaneVZ = kInvalidAlignValue;
                mAlignFstMeasGlobalX = kInvalidAlignValue;
                mAlignFstMeasGlobalY = kInvalidAlignValue;
                mAlignFstMeasGlobalZ = kInvalidAlignValue;
                mAlignFstClosureX = kInvalidAlignValue;
                mAlignFstClosureY = kInvalidAlignValue;
                mAlignFstClosureZ = kInvalidAlignValue;
                mAlignFstClosureU = kInvalidAlignValue;
                mAlignFstClosureV = kInvalidAlignValue;
                mAlignFstClosureMag = kInvalidAlignValue;

                setAlignmentVector(rawMeasurement->getRawHitCoords(), mAlignMeas0, mAlignMeas1, mAlignMeas2);
                if (mAlignDetId == kFstId && mAlignFstGlobalSensor >= 0) {
                    const FwdHit *fwdHit = findFstSeedHit(gtr, mAlignFstGlobalSensor);
                    if (fwdHit) {
                        mAlignFstHitGlobalX = fwdHit->getX();
                        mAlignFstHitGlobalY = fwdHit->getY();
                        mAlignFstHitGlobalZ = fwdHit->getZ();

                        if (fwdHit->_localPosition[0] >= 0.f && fwdHit->_localPosition[1] >= 0.f) {
                            mAlignFstRawR = fwdHit->_localPosition[0];
                            mAlignFstRawStripPhi = fwdHit->_localPosition[1];
                            mAlignFstMeanPhiStrip = fwdHit->_localPosition[1] / kFstStripPitchPhi;
                        }

                        if (alignmentGeo && mAlignMeasurementDim >= 2 &&
                            mAlignMeas0 > kInvalidAlignValue &&
                            mAlignMeas1 > kInvalidAlignValue) {
                            TVector3 u(1, 0, 0);
                            TVector3 v(0, 1, 0);
                            const int globalWedge =
                                mAlignFstDisk * kFstNumWedgePerDisk + mAlignFstWedge;
                            TVector3 o = alignmentGeo->getFstWedgeOrigin(globalWedge, u, v);
                            o.SetZ(alignmentGeo->getFstMeasurementZ(
                                mAlignFstGlobalSensor, o.Z()));
                            TVector3 measGlobal = o + mAlignMeas0 * u + mAlignMeas1 * v;
                            TVector3 hitGlobal(fwdHit->getX(), fwdHit->getY(), fwdHit->getZ());
                            TVector3 closure = measGlobal - hitGlobal;

                            mAlignFstPlaneOriginX = o.X();
                            mAlignFstPlaneOriginY = o.Y();
                            mAlignFstPlaneOriginZ = o.Z();
                            mAlignFstPlaneUX = u.X();
                            mAlignFstPlaneUY = u.Y();
                            mAlignFstPlaneUZ = u.Z();
                            mAlignFstPlaneVX = v.X();
                            mAlignFstPlaneVY = v.Y();
                            mAlignFstPlaneVZ = v.Z();
                            mAlignFstMeasGlobalX = measGlobal.X();
                            mAlignFstMeasGlobalY = measGlobal.Y();
                            mAlignFstMeasGlobalZ = measGlobal.Z();
                            mAlignFstClosureX = closure.X();
                            mAlignFstClosureY = closure.Y();
                            mAlignFstClosureZ = closure.Z();
                            mAlignFstClosureU = closure.Dot(u);
                            mAlignFstClosureV = closure.Dot(v);
                            mAlignFstClosureMag = closure.Mag();
                        }
                    }
                }

                if (kfi && iMeas < kfi->getNumMeasurements()) {
                    try {
                        genfit::MeasurementOnPlane biasedResidual = kfi->getResidual(iMeas, true, true);
                        genfit::MeasurementOnPlane unbiasedResidual = kfi->getResidual(iMeas, false, true);
                        mAlignResidualDim = static_cast<int>(unbiasedResidual.getState().GetNrows());
                        setAlignmentVector(biasedResidual.getState(), mAlignResBiased0, mAlignResBiased1, mAlignResBiased2);
                        setAlignmentVector(unbiasedResidual.getState(), mAlignResUnbiased0, mAlignResUnbiased1, mAlignResUnbiased2);
                        setAlignmentPrediction(
                            rawMeasurement->getRawHitCoords(),
                            unbiasedResidual.getState(),
                            mAlignTrackPred0, mAlignTrackPred1, mAlignTrackPred2
                        );
                        setAlignmentPulls(
                            biasedResidual,
                            mAlignResBiasedSigma0, mAlignResBiasedSigma1, mAlignResBiasedSigma2,
                            mAlignPullBiased0, mAlignPullBiased1, mAlignPullBiased2
                        );
                        setAlignmentPulls(
                            unbiasedResidual,
                            mAlignResUnbiasedSigma0, mAlignResUnbiasedSigma1, mAlignResUnbiasedSigma2,
                            mAlignPullUnbiased0, mAlignPullUnbiased1, mAlignPullUnbiased2
                        );
                        mAlignHasResidual = 1;
                    } catch (...) {
                        mAlignHasResidual = 0;
                    }
                }

                mAlignmentTree->Fill();
            }
        }
    }

    if (alignmentGeo)
        delete alignmentGeo;
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

    size_t indexTrack = 0;
    for ( auto &gtr : mForwardTracker->getTrackResults() ) {
            LOG_INFO << "Processing GenfitTrackResult(type=" << gtr.mTrackType << "): " << indexTrack << " mIsFitConverged=" << gtr.mIsFitConverged << ", mIsFitConvergedPartially=" << gtr.mIsFitConvergedPartially << ", mNumFitPoints=" << gtr.mNumFitPoints << endm;
            StFwdTrack* fwdTrack = makeStFwdTrack( gtr, indexTrack );
            indexTrack++;
            if (nullptr == fwdTrack)
                continue;
            ftc->addTrack( fwdTrack );
    }

    LOG_INFO << "StFwdTrackCollection has " << ftc->numberOfTracks() << " tracks now" << endm;


    // get the vertices from the forward tracker
    // and add them to the StEvent as Primary vertices
    auto fwdVertices = mForwardTracker->getVertices();
    for ( auto vert : fwdVertices ){
        StPrimaryVertex *pv = new StPrimaryVertex();
        pv->setPosition( StThreeVectorF( vert->getPos().X(), vert->getPos().Y(), vert->getPos().Z() ) );
        pv->setCovariantMatrix( vert->getCov().GetMatrixArray() );
        pv->setChiSquared( vert->getChi2() );
        pv->setNumTracksUsedInFinder( vert->getNTracks() );
        pv->setFwdVertex();
        stEvent->addPrimaryVertex( pv );
    }


    // Pico Dst requires a primary vertex,
    // if we have a PicoDst maker in the chain, we need to add a primary vertex
    // when one does not exist to get a "FWD" picoDst
    // auto mk = GetMaker("PicoDst");
    // LOG_INFO << "stEvent->numberOfPrimaryVertices() = " << stEvent->numberOfPrimaryVertices() << endm;
    // if ( mk && stEvent->numberOfPrimaryVertices() == 0 ){
    //     LOG_INFO << "Adding a primary vertex to StEvent since PicoDst maker was found in chain, but no vertices found" << endm;
    //     stEvent->addPrimaryVertex( new StPrimaryVertex() );
    //     LOG_INFO << "StPrimaryVertex::numberOfPrimaryVertices = " << stEvent->numberOfPrimaryVertices() << endm;
    // }



    LOG_INFO << "StFwdTrackCollection has " << ftc->numberOfTracks() << " tracks now" << endm;
}




//________________________________________________________________________
void StFwdTrackMaker::Clear(const Option_t *opts) {
    LOG_DEBUG << "StFwdTrackMaker::Clear" << endm;
    mForwardData->clear();
    mFwdHitLoader.clear();
    mForwardTracker->Clear();
}


std::string StFwdTrackMaker::defaultConfig = R"(
<?xml version="1.0" encoding="UTF-8"?>
<config>
    <TrackFinder nIterations="1">
        <Iteration nPhiSlices="1" > <!-- Options for first iteration -->
            <SegmentBuilder>
                <!-- <Criteria name="Crit2_RZRatio" min="0" max="1.20" /> -->
                <!-- <Criteria name="Crit2_DeltaRho" min="-50" max="50.9"/> -->
                <Criteria name="Crit2_DeltaPhi" min="0" max="2.0" />
                <!-- <Criteria name="Crit2_StraightTrackRatio" min="0.01" max="5.85"/> -->
            </SegmentBuilder>

            <ThreeHitSegments>
				<Criteria name="Crit3_3DAngle" min="0" max="1" />
                <!-- <Criteria name="Crit3_PT" min="0" max="100" /> -->
				<!-- <Criteria name="Crit3_ChangeRZRatio" min="0.8" max="1.21" /> -->
				<Criteria name="Crit3_2DAngle" min="0" max="1" />
            </ThreeHitSegments>

        </Iteration>

        <Connector distance="1"/>

        <SubsetNN active="true" min-hits-on-track="3" >
            <!-- <InitialTemp>2.1</InitialTemp> -->
            <!-- <InfTemp>0.1</InfTemp> -->
            <Omega>0.99</Omega>
            <StableThreshold>0.001</StableThreshold>
        </SubsetNN>

        <HitRemover active="false" />
    </TrackFinder>
</config>
)";


const std::vector<Seed_t> &StFwdTrackMaker::getTrackSeeds() const{
    return mForwardTracker->getTrackSeeds();
}

const std::vector<GenfitTrackResult> &StFwdTrackMaker::getFitResults()const{
    return mForwardTracker->getTrackResults();
}
