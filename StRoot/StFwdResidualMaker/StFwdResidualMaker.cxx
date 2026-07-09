#include "StFwdResidualMaker.h"

#include "StEvent/StEvent.h"
#include "StEvent/StFwdTrackCollection.h"
#include "StEvent/StFwdTrack.h"
#include "StEvent/StEnumerations.h"

#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrack.h"
#include "StMuDSTMaker/COMMON/StMuFwdTrackCollection.h"

#include "TDirectory.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <cstring>

ClassImp(StFwdResidualMaker)

StFwdResidualMaker::StFwdResidualMaker(const char* outFile, const char* makerName)
    : StMaker(makerName), mOutFile(outFile) {}

StFwdResidualMaker::~StFwdResidualMaker() {}

Int_t StFwdResidualMaker::Init() {
    mFout = new TFile(mOutFile, "RECREATE");
    bookHistos();
    return kStOK;
}

void StFwdResidualMaker::bookHistos() {
    const float FST_RRES = 1.0;   // FST RMS is ~0.03-0.05 cm; zoom in from the old +-5 cm
    const float RRES = 2.5;       // FTT RMS is ~0.6-1.2 cm; zoom in from the old +-5 cm
    const float RPOS = 65.0;
    const int   NR   = 100;
    const int   NP   = 80;

    const char* dimLbl[4] = {"x [cm]", "y [cm]", "r [cm]", "r#phi [cm]"};

    // FST per disk
    TDirectory* dFst = mFout->mkdir("FST");
    for (int d = 0; d < 3; d++) {
        TDirectory* dd = dFst->mkdir(Form("disk%d", d));
        dd->cd();
        hFst[d] = new TH1F(Form("h_fst_d%d_rdphi", d),
                            Form("FST disk %d;r#cdot#Delta#phi [cm];", d),
                            NR, -FST_RRES, FST_RRES);
        float xlo[4] = {-RPOS, -RPOS, 0,    -RPOS};
        float xhi[4] = { RPOS,  RPOS, RPOS,  RPOS};
        const char* rlbl = "r#cdot#Delta#phi [cm]";
        h2Fst[d][0] = new TH2F(Form("h2_fst_d%d_rdphi_vs_x",    d),
            Form("FST disk %d vs x;%s;%s",    d, dimLbl[0], rlbl), NP, xlo[0], xhi[0], NR, -FST_RRES, FST_RRES);
        h2Fst[d][1] = new TH2F(Form("h2_fst_d%d_rdphi_vs_y",    d),
            Form("FST disk %d vs y;%s;%s",    d, dimLbl[1], rlbl), NP, xlo[1], xhi[1], NR, -FST_RRES, FST_RRES);
        h2Fst[d][2] = new TH2F(Form("h2_fst_d%d_rdphi_vs_r",    d),
            Form("FST disk %d vs r;%s;%s",    d, dimLbl[2], rlbl), 60, xlo[2], xhi[2], NR, -FST_RRES, FST_RRES);
        h2Fst[d][3] = new TH2F(Form("h2_fst_d%d_rdphi_vs_rphi", d),
            Form("FST disk %d vs r#phi;%s;%s",d, dimLbl[3], rlbl), NP, xlo[3], xhi[3], NR, -FST_RRES, FST_RRES);
        h2Fst[d][4] = new TH2F(Form("h2_fst_d%d_rdphi_vs_phi", d),
            Form("FST disk %d vs #phi;#phi [rad];%s",d, rlbl), NP, -TMath::Pi(), TMath::Pi(), NR, -FST_RRES, FST_RRES);
    }

    // FTT per plane
    TDirectory* dFtt = mFout->mkdir("FTT");
    for (int p = 0; p < 4; p++) {
        TDirectory* dp = dFtt->mkdir(Form("plane%d", p));
        dp->cd();

        float xlo[4] = {-RPOS, -RPOS, 0,    -RPOS};
        float xhi[4] = { RPOS,  RPOS, RPOS,  RPOS};

        // Horizontal strips → dy
        hFttDy[p] = new TH1F(Form("h_ftt_p%d_dy", p),
                              Form("FTT plane %d H-strip;#Deltay [cm];", p),
                              NR, -RRES, RRES);
        h2FttDy[p][0] = new TH2F(Form("h2_ftt_p%d_dy_vs_x",    p),
            Form("FTT plane %d H-strip vs x;%s;#Deltay [cm]",    p, dimLbl[0]), NP, xlo[0], xhi[0], NR, -RRES, RRES);
        h2FttDy[p][1] = new TH2F(Form("h2_ftt_p%d_dy_vs_y",    p),
            Form("FTT plane %d H-strip vs y;%s;#Deltay [cm]",    p, dimLbl[1]), NP, xlo[1], xhi[1], NR, -RRES, RRES);
        h2FttDy[p][2] = new TH2F(Form("h2_ftt_p%d_dy_vs_r",    p),
            Form("FTT plane %d H-strip vs r;%s;#Deltay [cm]",    p, dimLbl[2]), 60, xlo[2], xhi[2], NR, -RRES, RRES);
        h2FttDy[p][3] = new TH2F(Form("h2_ftt_p%d_dy_vs_rphi", p),
            Form("FTT plane %d H-strip vs r#phi;%s;#Deltay [cm]",p, dimLbl[3]), NP, xlo[3], xhi[3], NR, -RRES, RRES);

        // Vertical strips → dx
        hFttDx[p] = new TH1F(Form("h_ftt_p%d_dx", p),
                              Form("FTT plane %d V-strip;#Deltax [cm];", p),
                              NR, -RRES, RRES);
        h2FttDx[p][0] = new TH2F(Form("h2_ftt_p%d_dx_vs_x",    p),
            Form("FTT plane %d V-strip vs x;%s;#Deltax [cm]",    p, dimLbl[0]), NP, xlo[0], xhi[0], NR, -RRES, RRES);
        h2FttDx[p][1] = new TH2F(Form("h2_ftt_p%d_dx_vs_y",    p),
            Form("FTT plane %d V-strip vs y;%s;#Deltax [cm]",    p, dimLbl[1]), NP, xlo[1], xhi[1], NR, -RRES, RRES);
        h2FttDx[p][2] = new TH2F(Form("h2_ftt_p%d_dx_vs_r",    p),
            Form("FTT plane %d V-strip vs r;%s;#Deltax [cm]",    p, dimLbl[2]), 60, xlo[2], xhi[2], NR, -RRES, RRES);
        h2FttDx[p][3] = new TH2F(Form("h2_ftt_p%d_dx_vs_rphi", p),
            Form("FTT plane %d V-strip vs r#phi;%s;#Deltax [cm]",p, dimLbl[3]), NP, xlo[3], xhi[3], NR, -RRES, RRES);
    }

    // Plane-usage histogram, one per track type
    const char* typeName[6] = {"Global", "BLC", "Primary", "FwdVtx", "BLCVtx", "FCSConstrained"};
    const char* binLbl[20] = {
        "Vtx", "FST1", "FST2", "FST3",
        "FTT1x", "FTT1y", "FTT1uv", "FTT2x", "FTT2y", "FTT2uv",
        "FTT3x", "FTT3y", "FTT3uv", "FTT4x", "FTT4y", "FTT4uv",
        "EPD", "FCSEcal", "FCSHcal", "AllTrk"
    };
    TDirectory* dUse = mFout->mkdir("PlaneUsage");
    dUse->cd();
    for (int t = 0; t < 6; t++) {
        hPlaneUsage[t] = new TH1F(Form("hPlaneUsage_%s", typeName[t]),
                                  Form("Plane usage, track type=%s;;tracks", typeName[t]),
                                  20, 0, 20);
        for (int b = 0; b < 20; b++) hPlaneUsage[t]->GetXaxis()->SetBinLabel(b+1, binLbl[b]);
    }

    mFout->cd();
}

void StFwdResidualMaker::fillFst(int disk, float res, float hx, float hy) {
    float r   = sqrt(hx*hx + hy*hy);
    float phi = atan2(hy, hx);
    hFst[disk]->Fill(res);
    h2Fst[disk][0]->Fill(hx,    res);
    h2Fst[disk][1]->Fill(hy,    res);
    h2Fst[disk][2]->Fill(r,     res);
    h2Fst[disk][3]->Fill(r*phi, res);
    h2Fst[disk][4]->Fill(phi,   res);
}

void StFwdResidualMaker::fillFttDy(int plane, float res, float hx, float hy) {
    float r   = sqrt(hx*hx + hy*hy);
    float phi = atan2(hy, hx);
    hFttDy[plane]->Fill(res);
    h2FttDy[plane][0]->Fill(hx,    res);
    h2FttDy[plane][1]->Fill(hy,    res);
    h2FttDy[plane][2]->Fill(r,     res);
    h2FttDy[plane][3]->Fill(r*phi, res);
}

void StFwdResidualMaker::fillFttDx(int plane, float res, float hx, float hy) {
    float r   = sqrt(hx*hx + hy*hy);
    float phi = atan2(hy, hx);
    hFttDx[plane]->Fill(res);
    h2FttDx[plane][0]->Fill(hx,    res);
    h2FttDx[plane][1]->Fill(hy,    res);
    h2FttDx[plane][2]->Fill(r,     res);
    h2FttDx[plane][3]->Fill(r*phi, res);
}

// ---- Source-agnostic physics core: operates only on the neutral FwdProj/FwdSeedPt
// structs, filled identically regardless of whether the track came from StEvent
// (afterburner refit) or MuDst (already-produced production file). ----

void StFwdResidualMaker::processFstPoints(
        const std::vector<FwdSeedPt>& fstPts,
        const std::vector<FwdProj>& fstProjs)
{
    for (unsigned int is = 0; is < fstPts.size(); is++) {
        const FwdSeedPt& sp = fstPts[is];
        float hz = sp.z;

        // find closest FST projection by z (within 15 cm)
        int best = -1;
        float bestDz = 15.0f;
        for (int ip = 0; ip < (int)fstProjs.size(); ip++) {
            float dz = fabs(fstProjs[ip].z - hz);
            if (dz < bestDz) { bestDz = dz; best = ip; }
        }
        if (best < 0 || best >= 3) continue;

        float hx  = sp.x;
        float hy  = sp.y;
        float r   = sqrt(hx*hx + hy*hy);
        if (r < 1.0f) continue;

        float hphi = atan2(hy, hx);
        float pphi = atan2(fstProjs[best].y, fstProjs[best].x);
        float dphi = hphi - pphi;
        while (dphi >  TMath::Pi()) dphi -= 2.f*TMath::Pi();
        while (dphi < -TMath::Pi()) dphi += 2.f*TMath::Pi();

        fillFst(best, r * dphi, hx, hy);
        mNFstRes++;
    }
}

void StFwdResidualMaker::processFttPoints(
        const std::vector<FwdSeedPt>& fttPts,
        const std::vector<FwdProj>& fttProjs)
{
    for (unsigned int is = 0; is < fttPts.size(); is++) {
        const FwdSeedPt& sp = fttPts[is];
        float hz = sp.z;

        int best = -1;
        float bestDz = 15.0f;
        for (int ip = 0; ip < (int)fttProjs.size(); ip++) {
            float dz = fabs(fttProjs[ip].z - hz);
            if (dz < bestDz) { bestDz = dz; best = ip; }
        }
        if (best < 0 || best >= 4) continue;

        float hx = sp.x;
        float hy = sp.y;
        float px = fttProjs[best].x;
        float py = fttProjs[best].y;

        float sigX = sqrt(fabs(sp.cov[0]));
        float sigY = sqrt(fabs(sp.cov[4]));

        if (sigX > sigY + 0.01f) {
            fillFttDy(best, hy - py, hx, hy);
            mNFttRes++;
        } else if (sigY > sigX + 0.01f) {
            fillFttDx(best, hx - px, hx, hy);
            mNFttRes++;
        }
        // diagonal: equal sigma, not in fit → skip
    }
}

int StFwdResidualMaker::findClosestZ(
        float hz, const std::vector<FwdProj>& projs, int maxIdx)
{
    int best = -1;
    float bestDz = 15.0f;
    for (int ip = 0; ip < (int)projs.size(); ip++) {
        float dz = fabs(projs[ip].z - hz);
        if (dz < bestDz) { bestDz = dz; best = ip; }
    }
    if (best >= maxIdx) best = -1;
    return best;
}

void StFwdResidualMaker::fillPlaneUsage(
        UChar_t trackType, bool hasEcal, bool hasHcal,
        const std::vector<FwdSeedPt>& fstPts,
        const std::vector<FwdSeedPt>& fttPts,
        const std::vector<FwdProj>& fstProjs,
        const std::vector<FwdProj>& fttProjs)
{
    if (trackType > 5) return;
    TH1F* h = hPlaneUsage[trackType];

    // trackType()==kGlobal(0) is unconstrained by definition; all other types
    // (BLC/Primary/FwdVtx/BLCVtx/FCSConstrained) are vertex-constrained.
    // (isPrimary()/mVtxIndex is bookkeeping, not a reliable signal of fit usage.)
    if (trackType != 0) h->Fill(0);

    for (unsigned int is = 0; is < fstPts.size(); is++) {
        int disk = findClosestZ(fstPts[is].z, fstProjs, 3);
        if (disk >= 0) h->Fill(1 + disk);
    }

    for (unsigned int is = 0; is < fttPts.size(); is++) {
        const FwdSeedPt& sp = fttPts[is];
        int plane = findClosestZ(sp.z, fttProjs, 4);
        if (plane < 0) continue;
        int base = 4 + 3*plane;   // x, y, uv
        float sigX = sqrt(fabs(sp.cov[0]));
        float sigY = sqrt(fabs(sp.cov[4]));
        if (sigY > sigX + 0.01f)      h->Fill(base + 0);   // vertical strip -> x
        else if (sigX > sigY + 0.01f) h->Fill(base + 1);   // horizontal strip -> y
        else                           h->Fill(base + 2);   // diagonal u/v (not used in fit yet)
    }

    // EPD: not yet tracked here — bin 16 left for future use
    if (hasEcal) h->Fill(17);
    if (hasHcal) h->Fill(18);

    h->Fill(19);
}

void StFwdResidualMaker::processTrack(
        UChar_t trackType, bool hasEcal, bool hasHcal,
        const std::vector<FwdSeedPt>& fstPts,
        const std::vector<FwdSeedPt>& fttPts,
        const std::vector<FwdProj>& fstProjs,
        const std::vector<FwdProj>& fttProjs)
{
    fillPlaneUsage(trackType, hasEcal, hasHcal, fstPts, fttPts, fstProjs, fttProjs);

    if (trackType != mResidualTrackType) return;  // residuals: selected track type only
    mNGoodTracks++;
    processFstPoints(fstPts, fstProjs);
    processFttPoints(fttPts, fttProjs);
}

// ---- Source-specific extraction: convert StEvent's StFwdTrack or MuDst's
// StMuFwdTrack into the neutral structs above, then hand off to processTrack(). ----

Int_t StFwdResidualMaker::makeFromStEvent() {
    StEvent* evt = (StEvent*)GetDataSet("StEvent");
    if (!evt) return kStWarn;

    StFwdTrackCollection* fwdCol = evt->fwdTrackCollection();
    if (!fwdCol || fwdCol->numberOfTracks() == 0) return kStOK;

    const StSPtrVecFwdTrack& tracks = fwdCol->tracks();
    for (unsigned int it = 0; it < tracks.size(); it++) {
        StFwdTrack* trk = tracks[it];
        if (!trk)                          continue;
        if (!trk->didFitConverge())        continue;
        if (trk->numberOfFitPoints() < 3)  continue;

        std::vector<FwdProj> fstProjs, fttProjs;
        for (unsigned int ip = 0; ip < trk->mProjections.size(); ip++) {
            const StFwdTrackProjection& p = trk->mProjections[ip];
            FwdProj fp; fp.x = p.mXYZ.x(); fp.y = p.mXYZ.y(); fp.z = p.mXYZ.z();
            if (p.mDetId == kFstId) fstProjs.push_back(fp);
            if (p.mDetId == kFttId) fttProjs.push_back(fp);
        }
        std::sort(fstProjs.begin(), fstProjs.end(),
            [](const FwdProj& a, const FwdProj& b){ return a.z < b.z; });
        std::sort(fttProjs.begin(), fttProjs.end(),
            [](const FwdProj& a, const FwdProj& b){ return a.z < b.z; });

        std::vector<FwdSeedPt> fstPts, fttPts;
        for (auto& sp : trk->mFSTPoints) {
            FwdSeedPt s; s.x = sp.mXYZ.x(); s.y = sp.mXYZ.y(); s.z = sp.mXYZ.z();
            memcpy(s.cov, sp.mCov, sizeof(s.cov));
            fstPts.push_back(s);
        }
        for (auto& sp : trk->mFTTPoints) {
            FwdSeedPt s; s.x = sp.mXYZ.x(); s.y = sp.mXYZ.y(); s.z = sp.mXYZ.z();
            memcpy(s.cov, sp.mCov, sizeof(s.cov));
            fttPts.push_back(s);
        }

        processTrack(trk->trackType(),
                     trk->ecalClusters().size() > 0,
                     trk->hcalClusters().size() > 0,
                     fstPts, fttPts, fstProjs, fttProjs);
    }
    return kStOK;
}

Int_t StFwdResidualMaker::makeFromMuDst() {
    StMuDst* muDst = (StMuDst*)GetInputDS("MuDst");
    if (!muDst) {
        StMuDstMaker* muMaker = (StMuDstMaker*)GetMaker("MuDst");
        if (muMaker) muDst = muMaker->muDst();
    }
    if (!muDst) return kStWarn;

    StMuFwdTrackCollection* fwdCol = muDst->muFwdTrackCollection();
    if (!fwdCol || fwdCol->numberOfFwdTracks() == 0) return kStOK;

    for (unsigned int it = 0; it < fwdCol->numberOfFwdTracks(); it++) {
        StMuFwdTrack* trk = fwdCol->getFwdTrack(it);
        if (!trk)                             continue;
        if (!trk->didFitConverge())            continue;
        if (trk->numberOfFitPoints() < 3)     continue;

        std::vector<FwdProj> fstProjs, fttProjs;
        for (unsigned int ip = 0; ip < trk->mProjections.size(); ip++) {
            const StMuFwdTrackProjection& p = trk->mProjections[ip];
            FwdProj fp; fp.x = p.mXYZ.X(); fp.y = p.mXYZ.Y(); fp.z = p.mXYZ.Z();
            if (p.mDetId == kFstId) fstProjs.push_back(fp);
            if (p.mDetId == kFttId) fttProjs.push_back(fp);
        }
        std::sort(fstProjs.begin(), fstProjs.end(),
            [](const FwdProj& a, const FwdProj& b){ return a.z < b.z; });
        std::sort(fttProjs.begin(), fttProjs.end(),
            [](const FwdProj& a, const FwdProj& b){ return a.z < b.z; });

        std::vector<FwdSeedPt> fstPts, fttPts;
        for (auto& sp : trk->mFSTPoints) {
            FwdSeedPt s; s.x = sp.mXYZ.X(); s.y = sp.mXYZ.Y(); s.z = sp.mXYZ.Z();
            memcpy(s.cov, sp.mCov, sizeof(s.cov));
            fstPts.push_back(s);
        }
        for (auto& sp : trk->mFTTPoints) {
            FwdSeedPt s; s.x = sp.mXYZ.X(); s.y = sp.mXYZ.Y(); s.z = sp.mXYZ.Z();
            memcpy(s.cov, sp.mCov, sizeof(s.cov));
            fttPts.push_back(s);
        }

        processTrack(trk->trackType(),
                     trk->mEcalClusters.GetEntriesFast() > 0,
                     trk->mHcalClusters.GetEntriesFast() > 0,
                     fstPts, fttPts, fstProjs, fttProjs);
    }
    return kStOK;
}

Int_t StFwdResidualMaker::Make() {
    return mUseMuDst ? makeFromMuDst() : makeFromStEvent();
}

Int_t StFwdResidualMaker::Finish() {
    if (!mFout) return kStOK;   // already Finish()'d (may be called both explicitly and at teardown)
    const char* typeName[6] = {"Global", "BLC", "Primary", "FwdVtx", "BLCVtx", "FCSConstrained"};
    const char* tn = (mResidualTrackType < 6) ? typeName[mResidualTrackType] : "?";
    LOG_INFO << "StFwdResidualMaker: " << mNGoodTracks << " " << tn << " tracks, "
             << mNFstRes << " FST residuals, " << mNFttRes << " FTT residuals" << endm;
    mFout->Write();
    mFout->Close();
    mFout = nullptr;
    LOG_INFO << "StFwdResidualMaker: wrote " << mOutFile << endm;
    return kStOK;
}
