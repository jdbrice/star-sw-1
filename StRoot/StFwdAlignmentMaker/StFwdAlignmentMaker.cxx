#include "StFwdAlignmentMaker.h"

#include "StFwdTrackMaker/StFwdTrackMaker.h"
#include "StFwdTrackMaker/include/Tracker/FwdTracker.h"
#include "StFwdTrackMaker/include/Tracker/TrackFitter.h"
#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StFwdTrackMaker/include/Tracker/GenfitTrackResult.h"

#include "StEvent/StEvent.h"
#include "StEvent/StEventInfo.h"

#include "GenFit/Track.h"
#include "GenFit/MeasuredStateOnPlane.h"
#include "GenFit/Exception.h"
#include "GenFit/FitStatus.h"

#include <vector>

ClassImp(StFwdAlignmentMaker)

StFwdAlignmentMaker::StFwdAlignmentMaker(const char* outFile, const char* makerName)
    : StMaker(makerName), mOutFile(outFile) {}

StFwdAlignmentMaker::~StFwdAlignmentMaker() {}

Int_t StFwdAlignmentMaker::Init() {
    if (!mFwdTrackMaker) {
        LOG_ERROR << "StFwdAlignmentMaker::Init -- setTrackMaker() was never called, "
                     "this maker has nothing to read. Refusing to run." << endm;
        return kStFatal;
    }
    mFout = new TFile(mOutFile, "RECREATE");
    bookTree();
    return kStOK;
}

void StFwdAlignmentMaker::bookTree() {
    mFout->cd();
    mTree = new TTree("alignTree", "Unbiased (hit-removed) FST/FTT residuals");
    mTree->Branch("run",             &b_run,             "run/I");
    mTree->Branch("event",           &b_event,           "event/I");
    mTree->Branch("trackId",         &b_trackId,         "trackId/I");
    mTree->Branch("trackType",       &b_trackType,       "trackType/b");
    mTree->Branch("detType",         &b_detType,         "detType/I");
    mTree->Branch("genfitPlaneIndex",&b_genfitPlaneIndex,"genfitPlaneIndex/I");
    mTree->Branch("stripDir",        &b_stripDir,        "stripDir/I");
    mTree->Branch("hitX",  &b_hitX,  "hitX/F");
    mTree->Branch("hitY",  &b_hitY,  "hitY/F");
    mTree->Branch("hitZ",  &b_hitZ,  "hitZ/F");
    mTree->Branch("projX", &b_projX, "projX/F");
    mTree->Branch("projY", &b_projY, "projY/F");
    mTree->Branch("projZ", &b_projZ, "projZ/F");
    mTree->Branch("chi2ndf",     &b_chi2ndf,     "chi2ndf/F");
    mTree->Branch("nPointsUsed", &b_nPointsUsed, "nPointsUsed/I");
    mTree->Branch("converged",   &b_converged,   "converged/O");
}

Int_t StFwdAlignmentMaker::Make() {
    if (!mFwdTrackMaker) return kStFatal;

    std::shared_ptr<ForwardTrackMaker> ft = mFwdTrackMaker->GetForwardTrackerBase();
    if (!ft) {
        LOG_WARN << "StFwdAlignmentMaker: no ForwardTrackMaker available this event, skipping" << endm;
        return kStOK;
    }

    StEvent* evt = (StEvent*)GetDataSet("StEvent");
    int runId = 0, eventId = 0;
    if (evt) {
        runId = evt->runId();
        if (evt->info()) eventId = evt->info()->id();
    }

    int trackId = 0;
    for (auto& gtr : ft->getTrackResults()) {
        if (gtr.mTrackType != mTrackType) continue;
        if (!gtr.mIsFitConvergedFully) continue;
        mNTracksSeen++;
        processTrack(gtr, ft.get(), runId, eventId, trackId);
        trackId++;
    }

    return kStOK;
}

void StFwdAlignmentMaker::processTrack(GenfitTrackResult& gtr, ForwardTrackMaker* ft,
                                        int runId, int eventId, int trackId) {
    TrackFitter* trackFitter = ft->getTrackFitter();
    if (!trackFitter) return;

    Seed_t& fullSeed = gtr.mSeed;

    // Count FST+FTT hits up front (never remove the PV point, so it doesn't count
    // toward "removable" but does count toward what stays in the reduced fit).
    int nPlaneHits = 0;
    for (auto h : fullSeed) {
        FwdHit* fh = dynamic_cast<FwdHit*>(h);
        if (fh && (fh->isFst() || fh->isFtt())) nPlaneHits++;
    }
    if (nPlaneHits < mMinRemainingHits + 1) return; // can't remove even one hit and stay above the floor

    bool testedAny = false;
    for (size_t i = 0; i < fullSeed.size(); i++) {
        FwdHit* fh = dynamic_cast<FwdHit*>(fullSeed[i]);
        if (!fh) continue;
        if (!(fh->isFst() || fh->isFtt())) continue; // never remove the vertex/PV point
        if (nPlaneHits - 1 < mMinRemainingHits) continue;

        Seed_t reduced;
        reduced.reserve(fullSeed.size() - 1);
        for (size_t j = 0; j < fullSeed.size(); j++) {
            if (j != i) reduced.push_back(fullSeed[j]);
        }

        GenfitTrackResult refit;
        try {
            refit = ft->fitTrack(reduced, &gtr.mMomentum, gtr.mCharge);
        } catch (genfit::Exception& e) {
            LOG_WARN << "StFwdAlignmentMaker: genfit exception on leave-one-out refit: " << e.what() << endm;
            continue;
        } catch (std::exception& e) {
            LOG_WARN << "StFwdAlignmentMaker: std exception on leave-one-out refit: " << e.what() << endm;
            continue;
        }
        testedAny = true;

        // NOTE: deliberately NOT using trackFitter->getPlaneFor(fh) here.
        // getPlaneFor() looks up fh->_genfit_plane_index into the sensor-plane
        // vectors built at Init() time -- but that index is only reliably
        // propagated onto FST hits at hit-loading time (StFwdHitLoader.cxx);
        // it was found to read back as 0 for every hit (FST and FTT alike) by
        // the time a hit reaches gtr.mSeed here, which silently projected
        // every removed hit onto plane index 0 regardless of which plane it
        // actually came from (caught by comparing hitZ vs projZ in testing --
        // they should track together and did not). Building a plane directly
        // at the hit's own z instead sidesteps that indexing issue entirely
        // and matches an existing fallback pattern already used in this
        // codebase (FwdTracker.h's addFstHits()/addFttHits(), which project to
        // a plane built from the mean z of nearby hits rather than a fixed
        // sensor-index plane, for the same reason).
        genfit::SharedPlanePtr plane = genfit::SharedPlanePtr(
            new genfit::DetPlane(TVector3(0, 0, fh->getZ()), TVector3(0, 0, 1)));

        bool converged = refit.mIsFitConvergedFully;
        float chi2ndf = 0;
        if (converged && refit.mTrack) {
            try {
                auto status = refit.mTrack->getFitStatus();
                if (status && status->getNdf() > 0) chi2ndf = status->getChi2() / status->getNdf();
            } catch (...) { converged = false; }
        }
        if (!converged) continue;

        try {
            genfit::MeasuredStateOnPlane msp = trackFitter->projectToPlane(plane, refit.mTrack);
            TVector3 projPos = msp.getPos();

            b_run = runId;
            b_event = eventId;
            b_trackId = trackId;
            b_trackType = gtr.mTrackType;
            b_detType = fh->isFst() ? 0 : 1;
            b_genfitPlaneIndex = (Int_t)fh->_genfit_plane_index;
            // FTT hits are 1D (a single strip measures either x or y, never both --
            // see proposal_alignment_path.txt); which one is which axis has the
            // tight covariance, same test StFwdResidualMaker::processFttPoints()
            // already uses (sigX vs sigY from the hit's own 3x3 cov matrix).
            b_stripDir = 0;
            if (fh->isFtt()) {
                double sigX = sqrt(fabs(fh->_covmat(0,0)));
                double sigY = sqrt(fabs(fh->_covmat(1,1)));
                if (sigX > sigY + 0.01) b_stripDir = 2;      // H-strip: measures y tightly
                else if (sigY > sigX + 0.01) b_stripDir = 1; // V-strip: measures x tightly
            }
            b_hitX = fh->getX(); b_hitY = fh->getY(); b_hitZ = fh->getZ();
            b_projX = projPos.X(); b_projY = projPos.Y(); b_projZ = projPos.Z();
            b_chi2ndf = chi2ndf;
            b_nPointsUsed = refit.mNumFitPoints;
            b_converged = converged;
            mTree->Fill();
            mNRowsWritten++;
        } catch (genfit::Exception& e) {
            LOG_WARN << "StFwdAlignmentMaker: genfit exception projecting refit: " << e.what() << endm;
        } catch (std::exception& e) {
            LOG_WARN << "StFwdAlignmentMaker: std exception projecting refit: " << e.what() << endm;
        }
    }
    if (testedAny) mNTracksTested++;
}

Int_t StFwdAlignmentMaker::Finish() {
    if (!mFout) return kStOK;
    LOG_INFO << "StFwdAlignmentMaker: " << mNTracksSeen << " tracks seen, "
             << mNTracksTested << " tracks with >=1 hit tested, "
             << mNRowsWritten << " unbiased-residual rows written" << endm;
    mFout->cd();
    mTree->Write();
    mFout->Write();
    mFout->Close();
    mFout = nullptr;
    LOG_INFO << "StFwdAlignmentMaker: wrote " << mOutFile << endm;
    return kStOK;
}
