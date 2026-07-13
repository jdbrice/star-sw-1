// Wedge-based ROOT QA macro for Forward STAR FST alignment.
//
// Usage:
//   root4star -l -b -q \
//     'fwd_alignment_residual_qa.C+("align_test.root","fwd_align_qa")'
//
// The trailing + compiles the macro with ACLiC. This is much faster than
// interpreting the explicit event loop with ROOT 5.
//
// Outputs:
//   fwd_align_qa.root
//   fwd_align_qa.pdf
//
// Coordinate convention used throughout this macro:
//   U = wedge radial axis
//   V = wedge counterclockwise azimuthal axis
//   r = sqrt(meas0^2 + meas1^2)
//   phi = atan2(meas1, meas0), relative to the wedge centerline
//   rphi = r * phi, with phi in radians
//
// The three FST readout regions in one wedge are deliberately not treated as
// separate alignment objects. Their rows are combined into one global wedge.

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
// USER EDIT SECTION

// Track-quality cuts applied to every residual, pull, prediction, and alignment
// profile. Uncut measurement maps are retained as geometry/acceptance checks.
const double kFwdQaTrackEtaMin = 2.5;
const double kFwdQaTrackEtaMax = 4.0;
const double kFwdQaTrackPMin = 1.0;  // GeV/c
const double kFwdQaTrackPtMin = 0.2; // GeV/c
const int kFwdQaDefaultMinTrackNHitsFit = 8;
const int kFwdQaDefaultMinTrackNFstHits = 3;

// -1 = all track types. 0 = Global, 1 = BLC, 2 = Primary, 3 = FwdVertex.
const int kFwdQaRequiredTrackType = 1;

// Wedge-local position binning. One bin is 1 cm in both r and rphi.
const int kFwdQaRBins = 25;
const double kFwdQaRMin = 5.0;
const double kFwdQaRMax = 30.0;
const int kFwdQaRPhiBins = 16;
const double kFwdQaRPhiMin = -8.0;
const double kFwdQaRPhiMax = 8.0;
const int kFwdQaPhiBins = 30;
const double kFwdQaPhiMinDeg = -15.0;
const double kFwdQaPhiMaxDeg = 15.0;

// Closure diagnostic ranges. The tree stores cm; plots use microns.
const int kFwdQaClosureBins = 160;
const double kFwdQaCmToMicron = 10000.0;
const double kFwdQaClosureSignedMinMicron = -20000.0;
const double kFwdQaClosureSignedMaxMicron = 20000.0;
const double kFwdQaClosureMagMaxMicron = 30000.0;

// Track-slope diagnostic binning.
const int kFwdQaEtaBins = 10;
const int kFwdQaPolarSlopeBins = 10;
const double kFwdQaPolarSlopeMin = 0.0;
const double kFwdQaPolarSlopeMax = 0.2;

///////////////////////////////////////////////////////////////////////////////
// Fixed detector constants.

const int kFwdQaNumFstDisks = 3;
const int kFwdQaNumWedgesPerDisk = 12;
const int kFwdQaNumGlobalWedges =
    kFwdQaNumFstDisks * kFwdQaNumWedgesPerDisk;
const int kFwdQaFstDetId = 45;
const float kFwdQaInvalid = -90000.f;

///////////////////////////////////////////////////////////////////////////////
// Small helpers.

bool fwdQaHasBranch(TTree *tree, const char *name) {
  return tree && tree->GetBranch(name);
}

bool fwdQaBindBranch(TTree *tree, const char *name, void *address,
                     bool required = true) {
  if (!fwdQaHasBranch(tree, name)) {
    if (required)
      std::cerr << "Missing required fwdAlign branch: " << name << std::endl;
    return false;
  }
  tree->SetBranchAddress(name, address);
  return true;
}

bool fwdQaValid(float value) { return value > kFwdQaInvalid; }

bool fwdQaValidDisk(int disk) {
  return disk >= 0 && disk < kFwdQaNumFstDisks;
}

bool fwdQaValidWedge(int wedge) {
  return wedge >= 0 && wedge < kFwdQaNumWedgesPerDisk;
}

int fwdQaGlobalWedge(int disk, int wedge) {
  if (!fwdQaValidDisk(disk) || !fwdQaValidWedge(wedge))
    return -1;
  return disk * kFwdQaNumWedgesPerDisk + wedge;
}

void fwdQaWedgePolar(double u, double v, double &r, double &phiRad,
                     double &phiDeg, double &rphi) {
  r = std::sqrt(u * u + v * v);
  phiRad = std::atan2(v, u);
  phiDeg = phiRad * TMath::RadToDeg();
  rphi = r * phiRad;
}

// Rotate a residual vector from wedge (U,V) into local radial/tangential
// components at the measured hit angle. resRPhi is a displacement in cm or um,
// matching the units supplied for resU/resV.
void fwdQaPolarResidual(double phiRad, double resU, double resV,
                        double &resR, double &resRPhi) {
  const double c = std::cos(phiRad);
  const double s = std::sin(phiRad);
  resR = c * resU + s * resV;
  resRPhi = -s * resU + c * resV;
}

void fwdQaDrawLine(double x1, double y1, double x2, double y2,
                   int color = kRed + 1, int style = 2) {
  TLine *line = new TLine(x1, y1, x2, y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->Draw();
}

void fwdQaDrawDiskBoundaries(double yMin, double yMax) {
  fwdQaDrawLine(11.5, yMin, 11.5, yMax, kGray + 2, 2);
  fwdQaDrawLine(23.5, yMin, 23.5, yMax, kGray + 2, 2);
}

void fwdQaLabelTrackTypeAxis(TH1D *hist) {
  if (!hist)
    return;
  hist->GetXaxis()->SetBinLabel(1, "Global");
  hist->GetXaxis()->SetBinLabel(2, "BLC");
  hist->GetXaxis()->SetBinLabel(3, "Primary");
  hist->GetXaxis()->SetBinLabel(4, "FwdVtx");
}

TH2D *fwdQaMakeUVMap(const TString &name, const TString &title) {
  TString fullTitle = title + ";wedge U [cm];wedge V [cm];rows";
  TH2D *hist = new TH2D(name.Data(), fullTitle.Data(), kFwdQaRBins,
                        kFwdQaRMin, kFwdQaRMax, kFwdQaRPhiBins,
                        kFwdQaRPhiMin, kFwdQaRPhiMax);
  hist->SetStats(false);
  return hist;
}

TH2D *fwdQaMakeRPhiVsR(const TString &name, const TString &title) {
  TString fullTitle = title + ";r [cm];r#phi [cm];rows";
  TH2D *hist = new TH2D(name.Data(), fullTitle.Data(), kFwdQaRBins,
                        kFwdQaRMin, kFwdQaRMax, kFwdQaRPhiBins,
                        kFwdQaRPhiMin, kFwdQaRPhiMax);
  hist->SetStats(false);
  return hist;
}

///////////////////////////////////////////////////////////////////////////////
// Branch bundle. Only branches used by this macro are bound.

struct FwdQaBranches {
  int detId;
  int fstDisk;
  int fstWedge;
  int measurementDim;
  int residualDim;
  int hasResidual;
  int fitConverged;
  int fitConvergedFully;
  int trackType;
  int trackNHitsFit;
  int trackNFstHits;

  float trackP;
  float trackPt;
  float trackPz;
  float trackEta;
  float meas0;
  float meas1;
  float trackPred0;
  float trackPred1;

  float resBiased0;
  float resBiased1;
  float resBiasedSigma0;
  float resBiasedSigma1;
  float pullBiased0;
  float pullBiased1;
  float resUnbiased0;
  float resUnbiased1;
  float resUnbiasedSigma0;
  float resUnbiasedSigma1;
  float pullUnbiased0;
  float pullUnbiased1;

  float fstClosureU;
  float fstClosureV;
  float fstClosureZ;
  float fstClosureMag;

  bool hasPullBranches;
  bool hasPredictionBranches;
  bool hasClosureBranches;
  bool hasTrackTypeBranch;
  bool hasTrackNHitsFitBranch;
  bool hasTrackNFstHitsBranch;
  bool hasTrackSlopeBranches;

  FwdQaBranches() {
    detId = 0;
    fstDisk = -1;
    fstWedge = -1;
    measurementDim = 0;
    residualDim = 0;
    hasResidual = 0;
    fitConverged = 0;
    fitConvergedFully = 0;
    trackType = -1;
    trackNHitsFit = 0;
    trackNFstHits = 0;
    trackP = 0;
    trackPt = 0;
    trackPz = 0;
    trackEta = 0;
    meas0 = kFwdQaInvalid;
    meas1 = kFwdQaInvalid;
    trackPred0 = kFwdQaInvalid;
    trackPred1 = kFwdQaInvalid;
    resBiased0 = kFwdQaInvalid;
    resBiased1 = kFwdQaInvalid;
    resBiasedSigma0 = kFwdQaInvalid;
    resBiasedSigma1 = kFwdQaInvalid;
    pullBiased0 = kFwdQaInvalid;
    pullBiased1 = kFwdQaInvalid;
    resUnbiased0 = kFwdQaInvalid;
    resUnbiased1 = kFwdQaInvalid;
    resUnbiasedSigma0 = kFwdQaInvalid;
    resUnbiasedSigma1 = kFwdQaInvalid;
    pullUnbiased0 = kFwdQaInvalid;
    pullUnbiased1 = kFwdQaInvalid;
    fstClosureU = kFwdQaInvalid;
    fstClosureV = kFwdQaInvalid;
    fstClosureZ = kFwdQaInvalid;
    fstClosureMag = kFwdQaInvalid;
    hasPullBranches = false;
    hasPredictionBranches = false;
    hasClosureBranches = false;
    hasTrackTypeBranch = false;
    hasTrackNHitsFitBranch = false;
    hasTrackNFstHitsBranch = false;
    hasTrackSlopeBranches = false;
  }

  bool bind(TTree *tree) {
    bool ok = true;
    ok &= fwdQaBindBranch(tree, "detId", &detId);
    ok &= fwdQaBindBranch(tree, "fstDisk", &fstDisk);
    ok &= fwdQaBindBranch(tree, "fstWedge", &fstWedge);
    ok &= fwdQaBindBranch(tree, "measurementDim", &measurementDim);
    ok &= fwdQaBindBranch(tree, "residualDim", &residualDim);
    ok &= fwdQaBindBranch(tree, "hasResidual", &hasResidual);
    ok &= fwdQaBindBranch(tree, "fitConverged", &fitConverged);
    ok &= fwdQaBindBranch(tree, "fitConvergedFully", &fitConvergedFully);
    ok &= fwdQaBindBranch(tree, "trackP", &trackP);
    ok &= fwdQaBindBranch(tree, "trackPt", &trackPt);
    ok &= fwdQaBindBranch(tree, "trackEta", &trackEta);
    ok &= fwdQaBindBranch(tree, "meas0", &meas0);
    ok &= fwdQaBindBranch(tree, "meas1", &meas1);
    ok &= fwdQaBindBranch(tree, "resBiased0", &resBiased0);
    ok &= fwdQaBindBranch(tree, "resBiased1", &resBiased1);
    ok &= fwdQaBindBranch(tree, "resUnbiased0", &resUnbiased0);
    ok &= fwdQaBindBranch(tree, "resUnbiased1", &resUnbiased1);

    hasTrackTypeBranch = fwdQaHasBranch(tree, "trackType");
    if (hasTrackTypeBranch)
      fwdQaBindBranch(tree, "trackType", &trackType, false);

    hasTrackNHitsFitBranch = fwdQaHasBranch(tree, "trackNHitsFit");
    if (hasTrackNHitsFitBranch)
      fwdQaBindBranch(tree, "trackNHitsFit", &trackNHitsFit, false);

    hasTrackNFstHitsBranch = fwdQaHasBranch(tree, "trackNFstHits");
    if (hasTrackNFstHitsBranch)
      fwdQaBindBranch(tree, "trackNFstHits", &trackNFstHits, false);

    hasPullBranches = fwdQaHasBranch(tree, "pullBiased0") &&
                      fwdQaHasBranch(tree, "pullBiased1") &&
                      fwdQaHasBranch(tree, "pullUnbiased0") &&
                      fwdQaHasBranch(tree, "pullUnbiased1") &&
                      fwdQaHasBranch(tree, "resBiasedSigma0") &&
                      fwdQaHasBranch(tree, "resBiasedSigma1") &&
                      fwdQaHasBranch(tree, "resUnbiasedSigma0") &&
                      fwdQaHasBranch(tree, "resUnbiasedSigma1");
    if (hasPullBranches) {
      fwdQaBindBranch(tree, "resBiasedSigma0", &resBiasedSigma0, false);
      fwdQaBindBranch(tree, "resBiasedSigma1", &resBiasedSigma1, false);
      fwdQaBindBranch(tree, "pullBiased0", &pullBiased0, false);
      fwdQaBindBranch(tree, "pullBiased1", &pullBiased1, false);
      fwdQaBindBranch(tree, "resUnbiasedSigma0", &resUnbiasedSigma0, false);
      fwdQaBindBranch(tree, "resUnbiasedSigma1", &resUnbiasedSigma1, false);
      fwdQaBindBranch(tree, "pullUnbiased0", &pullUnbiased0, false);
      fwdQaBindBranch(tree, "pullUnbiased1", &pullUnbiased1, false);
    }

    hasPredictionBranches = fwdQaHasBranch(tree, "trackPred0") &&
                            fwdQaHasBranch(tree, "trackPred1");
    if (hasPredictionBranches) {
      fwdQaBindBranch(tree, "trackPred0", &trackPred0, false);
      fwdQaBindBranch(tree, "trackPred1", &trackPred1, false);
    }

    hasClosureBranches = fwdQaHasBranch(tree, "fstClosureU") &&
                         fwdQaHasBranch(tree, "fstClosureV") &&
                         fwdQaHasBranch(tree, "fstClosureZ") &&
                         fwdQaHasBranch(tree, "fstClosureMag");
    if (hasClosureBranches) {
      fwdQaBindBranch(tree, "fstClosureU", &fstClosureU, false);
      fwdQaBindBranch(tree, "fstClosureV", &fstClosureV, false);
      fwdQaBindBranch(tree, "fstClosureZ", &fstClosureZ, false);
      fwdQaBindBranch(tree, "fstClosureMag", &fstClosureMag, false);
    }

    hasTrackSlopeBranches = fwdQaHasBranch(tree, "trackPz");
    if (hasTrackSlopeBranches)
      fwdQaBindBranch(tree, "trackPz", &trackPz, false);

    return ok;
  }
};

///////////////////////////////////////////////////////////////////////////////

void fwd_alignment_residual_qa(
    const char *inputFilename = "align_test.root",
    const char *outputPrefix = "fwd_align_qa",
    bool requireFullyConverged = false, double residualRangeMicron = 5000.0,
    int minTrackNFstHits = kFwdQaDefaultMinTrackNFstHits,
    int minTrackNHitsFit = kFwdQaDefaultMinTrackNHitsFit) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(1);

  TFile *inputFile = TFile::Open(inputFilename, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Cannot open input file: " << inputFilename << std::endl;
    return;
  }

  TTree *tree = dynamic_cast<TTree *>(inputFile->Get("fwdAlign"));
  if (!tree) {
    std::cerr << "Cannot find TTree fwdAlign in " << inputFilename << std::endl;
    inputFile->ls();
    return;
  }

  FwdQaBranches row;
  if (!row.bind(tree)) {
    std::cerr << "Stopping because required fwdAlign branches are missing."
              << std::endl;
    return;
  }

  TString rootOutput = TString::Format("%s.root", outputPrefix);
  TString pdfOutput = TString::Format("%s.pdf", outputPrefix);
  TFile *outputFile = TFile::Open(rootOutput, "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Cannot create output file: " << rootOutput.Data()
              << std::endl;
    return;
  }
  outputFile->cd();

  TString fitCut =
      requireFullyConverged ? "fitConvergedFully>0" : "fitConverged>0";
  TString trackCut = TString::Format(
      "trackEta>=%.1f&&trackEta<=%.1f&&trackP>=%.1f&&trackPt>=%.1f",
      kFwdQaTrackEtaMin, kFwdQaTrackEtaMax, kFwdQaTrackPMin,
      kFwdQaTrackPtMin);
  if (kFwdQaRequiredTrackType >= 0 && row.hasTrackTypeBranch)
    trackCut += TString::Format("&&trackType==%d", kFwdQaRequiredTrackType);
  if (minTrackNHitsFit > 0 && row.hasTrackNHitsFitBranch)
    trackCut += TString::Format("&&trackNHitsFit>=%d", minTrackNHitsFit);
  if (minTrackNFstHits > 0 && row.hasTrackNFstHitsBranch)
    trackCut += TString::Format("&&trackNFstHits>=%d", minTrackNFstHits);

  /////////////////////////////////////////////////////////////////////////////
  // Histogram booking.

  TH1D *hDetId =
      new TH1D("hDetId", "Alignment rows by detector;detId;rows", 80, 0, 80);
  TH2D *hDim = new TH2D(
      "hDim", "FST measurement and residual dimensions;measurementDim;"
              "residualDim;rows",
      6, -0.5, 5.5, 6, -0.5, 5.5);
  TH1D *hHasResidual = new TH1D(
      "hHasResidual", "FST hasResidual flag;hasResidual;rows", 3, -0.5, 2.5);
  TH1D *hTrackTypeFstRows = new TH1D(
      "hTrackTypeFstRows", "FST rows by track type;track type;rows", 4, -0.5,
      3.5);
  TH1D *hTrackTypeResidualRows = new TH1D(
      "hTrackTypeResidualRows",
      "Selected FST residual rows by track type;track type;rows", 4, -0.5,
      3.5);
  fwdQaLabelTrackTypeAxis(hTrackTypeFstRows);
  fwdQaLabelTrackTypeAxis(hTrackTypeResidualRows);

  TH2D *hUV = fwdQaMakeUVMap("hUV", "FST wedge-local measurement map");
  TH2D *hRPhiVsR =
      fwdQaMakeRPhiVsR("hRPhiVsR", "FST wedge-local r#phi vs r");
  TH2D *hRPhiVsPhi = new TH2D(
      "hRPhiVsPhi", "FST wedge-local r#phi vs #phi;#phi [deg];r#phi [cm];rows",
      kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg, kFwdQaRPhiBins,
      kFwdQaRPhiMin, kFwdQaRPhiMax);
  hRPhiVsPhi->SetStats(false);

  TH2D *hUVDisk[kFwdQaNumFstDisks] = {0};
  TH2D *hRPhiVsRDisk[kFwdQaNumFstDisks] = {0};
  TH2D *hRPhiVsPhiDisk[kFwdQaNumFstDisks] = {0};
  for (int disk = 0; disk < kFwdQaNumFstDisks; ++disk) {
    hUVDisk[disk] = fwdQaMakeUVMap(
        TString::Format("hUV_disk%d", disk),
        TString::Format("FST disk %d wedge-local measurement map", disk));
    hRPhiVsRDisk[disk] = fwdQaMakeRPhiVsR(
        TString::Format("hRPhiVsR_disk%d", disk),
        TString::Format("FST disk %d wedge-local r#phi vs r", disk));
    hRPhiVsPhiDisk[disk] = new TH2D(
        TString::Format("hRPhiVsPhi_disk%d", disk),
        TString::Format("FST disk %d wedge-local r#phi vs #phi;#phi [deg];"
                        "r#phi [cm];rows",
                        disk),
        kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg, kFwdQaRPhiBins,
        kFwdQaRPhiMin, kFwdQaRPhiMax);
    hRPhiVsPhiDisk[disk]->SetStats(false);
  }

  TH2D *hUVGlobalWedge[kFwdQaNumGlobalWedges] = {0};
  TH2D *hRPhiVsRGlobalWedge[kFwdQaNumGlobalWedges] = {0};
  TProfile2D *pResUMapGlobalWedge[kFwdQaNumGlobalWedges] = {0};
  TProfile2D *pResVMapGlobalWedge[kFwdQaNumGlobalWedges] = {0};
  for (int globalWedge = 0; globalWedge < kFwdQaNumGlobalWedges;
       ++globalWedge) {
    const int disk = globalWedge / kFwdQaNumWedgesPerDisk;
    const int wedge = globalWedge % kFwdQaNumWedgesPerDisk;
    hUVGlobalWedge[globalWedge] = fwdQaMakeUVMap(
        TString::Format("hUV_globalWedge%d", globalWedge),
        TString::Format("FST wedge-local map, global wedge %d (d%d w%d)",
                        globalWedge, disk, wedge));
    hRPhiVsRGlobalWedge[globalWedge] = fwdQaMakeRPhiVsR(
        TString::Format("hRPhiVsR_globalWedge%d", globalWedge),
        TString::Format("FST r#phi vs r, global wedge %d (d%d w%d)",
                        globalWedge, disk, wedge));
    pResUMapGlobalWedge[globalWedge] = new TProfile2D(
        TString::Format("pResUMap_globalWedge%d", globalWedge),
        TString::Format("Mean residual U, global wedge %d (d%d w%d);r [cm];"
                        "r#phi [cm];<resU> [um]",
                        globalWedge, disk, wedge),
        kFwdQaRBins, kFwdQaRMin, kFwdQaRMax, kFwdQaRPhiBins, kFwdQaRPhiMin,
        kFwdQaRPhiMax);
    pResVMapGlobalWedge[globalWedge] = new TProfile2D(
        TString::Format("pResVMap_globalWedge%d", globalWedge),
        TString::Format("Mean residual V, global wedge %d (d%d w%d);r [cm];"
                        "r#phi [cm];<resV> [um]",
                        globalWedge, disk, wedge),
        kFwdQaRBins, kFwdQaRMin, kFwdQaRMax, kFwdQaRPhiBins, kFwdQaRPhiMin,
        kFwdQaRPhiMax);
  }

  TH1D *hResidualU = new TH1D(
      "hResidualU", "FST unbiased wedge-U residual;resU [um];rows", 160,
      -residualRangeMicron, residualRangeMicron);
  TH1D *hResidualV = new TH1D(
      "hResidualV", "FST unbiased wedge-V residual;resV [um];rows", 160,
      -residualRangeMicron, residualRangeMicron);
  TH1D *hResidualR = new TH1D(
      "hResidualR", "FST unbiased radial residual;resR [um];rows", 160,
      -residualRangeMicron, residualRangeMicron);
  TH1D *hResidualRPhi = new TH1D(
      "hResidualRPhi", "FST unbiased tangential residual;resR#phi [um];rows",
      160, -residualRangeMicron, residualRangeMicron);
  TH1D *hBiasedU = new TH1D(
      "hBiasedU", "FST biased wedge-U residual;biased resU [um];rows", 160,
      -residualRangeMicron, residualRangeMicron);
  TH1D *hBiasedV = new TH1D(
      "hBiasedV", "FST biased wedge-V residual;biased resV [um];rows", 160,
      -residualRangeMicron, residualRangeMicron);
  TH1D *hBiasedR = new TH1D(
      "hBiasedR", "FST biased radial residual;biased resR [um];rows", 160,
      -residualRangeMicron, residualRangeMicron);
  TH1D *hBiasedRPhi = new TH1D(
      "hBiasedRPhi", "FST biased tangential residual;biased resR#phi [um];rows",
      160, -residualRangeMicron, residualRangeMicron);

  TH1D *hGlobalWedge = new TH1D(
      "hGlobalWedge", "Selected residual rows by global wedge;12*disk+wedge;rows",
      kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5);
  TH2D *hWedgeDisk = new TH2D(
      "hWedgeDisk", "Selected residual occupancy;disk;wedge;rows",
      kFwdQaNumFstDisks, -0.5, kFwdQaNumFstDisks - 0.5,
      kFwdQaNumWedgesPerDisk, -0.5, kFwdQaNumWedgesPerDisk - 0.5);

  TProfile *pResUByGlobalWedge = new TProfile(
      "pResUByGlobalWedge",
      "Mean unbiased wedge-U residual;global wedge = 12*disk+wedge;"
      "<resU> [um]",
      kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5);
  TProfile *pResVByGlobalWedge = new TProfile(
      "pResVByGlobalWedge",
      "Mean unbiased wedge-V residual;global wedge = 12*disk+wedge;"
      "<resV> [um]",
      kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5);
  TProfile *pResRByGlobalWedge = new TProfile(
      "pResRByGlobalWedge",
      "Mean unbiased radial residual;global wedge = 12*disk+wedge;"
      "<resR> [um]",
      kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5);
  TProfile *pResRPhiByGlobalWedge = new TProfile(
      "pResRPhiByGlobalWedge",
      "Mean unbiased tangential residual;global wedge = 12*disk+wedge;"
      "<resR#phi> [um]",
      kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5);

  TH2D *hResUByGlobalWedge = new TH2D(
      "hResUByGlobalWedge",
      "Wedge-U residual by global wedge;global wedge;resU [um];rows",
      kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5, 160,
      -residualRangeMicron, residualRangeMicron);
  TH2D *hResVByGlobalWedge = new TH2D(
      "hResVByGlobalWedge",
      "Wedge-V residual by global wedge;global wedge;resV [um];rows",
      kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5, 160,
      -residualRangeMicron, residualRangeMicron);

  TH2D *hResUVsU = new TH2D(
      "hResUVsU", "Wedge-U residual vs U;U [cm];resU [um];rows", kFwdQaRBins,
      kFwdQaRMin, kFwdQaRMax, 160, -residualRangeMicron,
      residualRangeMicron);
  TH2D *hResUVsV = new TH2D(
      "hResUVsV", "Wedge-U residual vs V;V [cm];resU [um];rows",
      kFwdQaRPhiBins, kFwdQaRPhiMin, kFwdQaRPhiMax, 160,
      -residualRangeMicron, residualRangeMicron);
  TH2D *hResVVsU = new TH2D(
      "hResVVsU", "Wedge-V residual vs U;U [cm];resV [um];rows", kFwdQaRBins,
      kFwdQaRMin, kFwdQaRMax, 160, -residualRangeMicron,
      residualRangeMicron);
  TH2D *hResVVsV = new TH2D(
      "hResVVsV", "Wedge-V residual vs V;V [cm];resV [um];rows",
      kFwdQaRPhiBins, kFwdQaRPhiMin, kFwdQaRPhiMax, 160,
      -residualRangeMicron, residualRangeMicron);

  TH2D *hResRVsR = new TH2D(
      "hResRVsR", "Radial residual vs r;r [cm];resR [um];rows", kFwdQaRBins,
      kFwdQaRMin, kFwdQaRMax, 160, -residualRangeMicron,
      residualRangeMicron);
  TH2D *hResRPhiVsR = new TH2D(
      "hResRPhiVsR", "Tangential residual vs r;r [cm];resR#phi [um];rows",
      kFwdQaRBins, kFwdQaRMin, kFwdQaRMax, 160, -residualRangeMicron,
      residualRangeMicron);
  TH2D *hResRVsPhi = new TH2D(
      "hResRVsPhi", "Radial residual vs wedge #phi;#phi [deg];resR [um];rows",
      kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg, 160,
      -residualRangeMicron, residualRangeMicron);
  TH2D *hResRPhiVsPhi = new TH2D(
      "hResRPhiVsPhi",
      "Tangential residual vs wedge #phi;#phi [deg];resR#phi [um];rows",
      kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg, 160,
      -residualRangeMicron, residualRangeMicron);

  TProfile *pResRVsR = new TProfile(
      "pResRVsR", "Mean radial residual vs r;r [cm];<resR> [um]", kFwdQaRBins,
      kFwdQaRMin, kFwdQaRMax);
  TProfile *pResRPhiVsR = new TProfile(
      "pResRPhiVsR",
      "Mean tangential residual vs r;r [cm];<resR#phi> [um]", kFwdQaRBins,
      kFwdQaRMin, kFwdQaRMax);
  TProfile *pResRVsPhi = new TProfile(
      "pResRVsPhi",
      "Mean radial residual vs wedge #phi;#phi [deg];<resR> [um]",
      kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg);
  TProfile *pResRPhiVsPhi = new TProfile(
      "pResRPhiVsPhi",
      "Mean tangential residual vs wedge #phi;#phi [deg];<resR#phi> [um]",
      kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg);

  TProfile *pResRVsRDisk[kFwdQaNumFstDisks] = {0};
  TProfile *pResRPhiVsRDisk[kFwdQaNumFstDisks] = {0};
  TProfile *pResRVsPhiDisk[kFwdQaNumFstDisks] = {0};
  TProfile *pResRPhiVsPhiDisk[kFwdQaNumFstDisks] = {0};
  for (int disk = 0; disk < kFwdQaNumFstDisks; ++disk) {
    pResRVsRDisk[disk] = new TProfile(
        TString::Format("pResRVsR_disk%d", disk),
        TString::Format("Disk %d mean radial residual vs r;r [cm];"
                        "<resR> [um]",
                        disk),
        kFwdQaRBins, kFwdQaRMin, kFwdQaRMax);
    pResRPhiVsRDisk[disk] = new TProfile(
        TString::Format("pResRPhiVsR_disk%d", disk),
        TString::Format("Disk %d mean tangential residual vs r;r [cm];"
                        "<resR#phi> [um]",
                        disk),
        kFwdQaRBins, kFwdQaRMin, kFwdQaRMax);
    pResRVsPhiDisk[disk] = new TProfile(
        TString::Format("pResRVsPhi_disk%d", disk),
        TString::Format("Disk %d mean radial residual vs #phi;#phi [deg];"
                        "<resR> [um]",
                        disk),
        kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg);
    pResRPhiVsPhiDisk[disk] = new TProfile(
        TString::Format("pResRPhiVsPhi_disk%d", disk),
        TString::Format("Disk %d mean tangential residual vs #phi;#phi [deg];"
                        "<resR#phi> [um]",
                        disk),
        kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg);
  }

  TH1D *hPullU = 0;
  TH1D *hPullV = 0;
  TH1D *hPullBiasedU = 0;
  TH1D *hPullBiasedV = 0;
  TProfile *pPullUByGlobalWedge = 0;
  TProfile *pPullVByGlobalWedge = 0;
  if (row.hasPullBranches) {
    hPullU = new TH1D("hPullU", "FST unbiased wedge-U pull;pullU;rows", 160,
                      -10.0, 10.0);
    hPullV = new TH1D("hPullV", "FST unbiased wedge-V pull;pullV;rows", 160,
                      -10.0, 10.0);
    hPullBiasedU = new TH1D(
        "hPullBiasedU", "FST biased wedge-U pull;biased pullU;rows", 160,
        -10.0, 10.0);
    hPullBiasedV = new TH1D(
        "hPullBiasedV", "FST biased wedge-V pull;biased pullV;rows", 160,
        -10.0, 10.0);
    pPullUByGlobalWedge = new TProfile(
        "pPullUByGlobalWedge",
        "Mean wedge-U pull;global wedge = 12*disk+wedge;<pullU>",
        kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5);
    pPullVByGlobalWedge = new TProfile(
        "pPullVByGlobalWedge",
        "Mean wedge-V pull;global wedge = 12*disk+wedge;<pullV>",
        kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5);
  }

  TH2D *hTrackPredUVsMeasU = 0;
  TH2D *hTrackPredVVsMeasV = 0;
  TH2D *hTrackPredRVsMeasR = 0;
  TH2D *hTrackPredRPhiVsMeasRPhi = 0;
  if (row.hasPredictionBranches) {
    hTrackPredUVsMeasU = new TH2D(
        "hTrackPredUVsMeasU", "Track prediction U vs measured U;measured U "
                              "[cm];predicted U [cm];rows",
        kFwdQaRBins, kFwdQaRMin, kFwdQaRMax, kFwdQaRBins, kFwdQaRMin,
        kFwdQaRMax);
    hTrackPredVVsMeasV = new TH2D(
        "hTrackPredVVsMeasV", "Track prediction V vs measured V;measured V "
                              "[cm];predicted V [cm];rows",
        kFwdQaRPhiBins, kFwdQaRPhiMin, kFwdQaRPhiMax, kFwdQaRPhiBins,
        kFwdQaRPhiMin, kFwdQaRPhiMax);
    hTrackPredRVsMeasR = new TH2D(
        "hTrackPredRVsMeasR", "Track prediction r vs measured r;measured r "
                              "[cm];predicted r [cm];rows",
        kFwdQaRBins, kFwdQaRMin, kFwdQaRMax, kFwdQaRBins, kFwdQaRMin,
        kFwdQaRMax);
    hTrackPredRPhiVsMeasRPhi = new TH2D(
        "hTrackPredRPhiVsMeasRPhi",
        "Track prediction r#phi vs measured r#phi;measured r#phi [cm];"
        "predicted r#phi [cm];rows",
        kFwdQaRPhiBins, kFwdQaRPhiMin, kFwdQaRPhiMax, kFwdQaRPhiBins,
        kFwdQaRPhiMin, kFwdQaRPhiMax);
  }

  TH1D *hClosureU = 0;
  TH1D *hClosureV = 0;
  TH1D *hClosureZ = 0;
  TH1D *hClosureMag = 0;
  TProfile *pClosureUByGlobalWedge = 0;
  TProfile *pClosureVByGlobalWedge = 0;
  TProfile *pClosureUVsR = 0;
  TProfile *pClosureVVsR = 0;
  TProfile *pClosureUVsPhi = 0;
  TProfile *pClosureVVsPhi = 0;
  if (row.hasClosureBranches) {
    hClosureU = new TH1D(
        "hClosureU", "Wedge-coordinate closure U;closureU [um];rows",
        kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    hClosureV = new TH1D(
        "hClosureV", "Wedge-coordinate closure V;closureV [um];rows",
        kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    hClosureZ = new TH1D(
        "hClosureZ",
        "Physical measurement z - flat FwdHit z;closureZ [um];rows",
        kFwdQaClosureBins, kFwdQaClosureSignedMinMicron,
        kFwdQaClosureSignedMaxMicron);
    hClosureMag = new TH1D(
        "hClosureMag", "3D closure magnitude;closure magnitude [um];rows",
        kFwdQaClosureBins, 0.0, kFwdQaClosureMagMaxMicron);
    pClosureUByGlobalWedge = new TProfile(
        "pClosureUByGlobalWedge",
        "Mean wedge closure U;global wedge;<closureU> [um]",
        kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5);
    pClosureVByGlobalWedge = new TProfile(
        "pClosureVByGlobalWedge",
        "Mean wedge closure V;global wedge;<closureV> [um]",
        kFwdQaNumGlobalWedges, -0.5, kFwdQaNumGlobalWedges - 0.5);
    pClosureUVsR = new TProfile(
        "pClosureUVsR", "Mean wedge closure U vs r;r [cm];<closureU> [um]",
        kFwdQaRBins, kFwdQaRMin, kFwdQaRMax);
    pClosureVVsR = new TProfile(
        "pClosureVVsR", "Mean wedge closure V vs r;r [cm];<closureV> [um]",
        kFwdQaRBins, kFwdQaRMin, kFwdQaRMax);
    pClosureUVsPhi = new TProfile(
        "pClosureUVsPhi",
        "Mean wedge closure U vs #phi;#phi [deg];<closureU> [um]",
        kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg);
    pClosureVVsPhi = new TProfile(
        "pClosureVVsPhi",
        "Mean wedge closure V vs #phi;#phi [deg];<closureV> [um]",
        kFwdQaPhiBins, kFwdQaPhiMinDeg, kFwdQaPhiMaxDeg);
  }

  TProfile *pResRVsEta = 0;
  TProfile *pResRPhiVsEta = 0;
  TProfile *pResRVsSlope = 0;
  TProfile *pResRPhiVsSlope = 0;
  if (row.hasTrackSlopeBranches) {
    pResRVsEta = new TProfile(
        "pResRVsEta", "Mean radial residual vs track #eta;track #eta;"
                       "<resR> [um]",
        kFwdQaEtaBins, kFwdQaTrackEtaMin, kFwdQaTrackEtaMax);
    pResRPhiVsEta = new TProfile(
        "pResRPhiVsEta", "Mean tangential residual vs track #eta;track #eta;"
                          "<resR#phi> [um]",
        kFwdQaEtaBins, kFwdQaTrackEtaMin, kFwdQaTrackEtaMax);
    pResRVsSlope = new TProfile(
        "pResRVsSlope", "Mean radial residual vs pT/pZ;pT/pZ;<resR> [um]",
        kFwdQaPolarSlopeBins, kFwdQaPolarSlopeMin, kFwdQaPolarSlopeMax);
    pResRPhiVsSlope = new TProfile(
        "pResRPhiVsSlope",
        "Mean tangential residual vs pT/pZ;pT/pZ;<resR#phi> [um]",
        kFwdQaPolarSlopeBins, kFwdQaPolarSlopeMin, kFwdQaPolarSlopeMax);
  }

  /////////////////////////////////////////////////////////////////////////////
  // Main tree loop.

  const Long64_t nEntries = tree->GetEntries();
  Long64_t nFstRows = 0;
  Long64_t nMeasurementRows = 0;
  Long64_t nSelectedResidualRows = 0;
  Long64_t nClosureRows = 0;

  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);

    const bool isFst = row.detId == kFwdQaFstDetId;
    const int globalWedge = fwdQaGlobalWedge(row.fstDisk, row.fstWedge);
    const bool validGlobalWedge =
        globalWedge >= 0 && globalWedge < kFwdQaNumGlobalWedges;
    const bool fstRow = isFst && validGlobalWedge;
    const bool measurementOk = fstRow && row.measurementDim == 2 &&
                               fwdQaValid(row.meas0) &&
                               fwdQaValid(row.meas1);

    const bool fitOk = requireFullyConverged ? row.fitConvergedFully > 0
                                             : row.fitConverged > 0;
    const bool trackTypeOk = kFwdQaRequiredTrackType < 0 ||
                             !row.hasTrackTypeBranch ||
                             row.trackType == kFwdQaRequiredTrackType;
    const bool trackOk =
        row.trackEta >= kFwdQaTrackEtaMin &&
        row.trackEta <= kFwdQaTrackEtaMax && row.trackP >= kFwdQaTrackPMin &&
        row.trackPt >= kFwdQaTrackPtMin && trackTypeOk &&
        (minTrackNHitsFit <= 0 || !row.hasTrackNHitsFitBranch ||
         row.trackNHitsFit >= minTrackNHitsFit) &&
        (minTrackNFstHits <= 0 || !row.hasTrackNFstHitsBranch ||
         row.trackNFstHits >= minTrackNFstHits);

    const bool residualOk =
        measurementOk && row.hasResidual > 0 && row.residualDim == 2 && fitOk &&
        trackOk && fwdQaValid(row.resBiased0) &&
        fwdQaValid(row.resBiased1) && fwdQaValid(row.resUnbiased0) &&
        fwdQaValid(row.resUnbiased1);
    const bool pullOk =
        residualOk && row.hasPullBranches && fwdQaValid(row.pullBiased0) &&
        fwdQaValid(row.pullBiased1) && fwdQaValid(row.pullUnbiased0) &&
        fwdQaValid(row.pullUnbiased1) && row.resBiasedSigma0 > 0 &&
        row.resBiasedSigma1 > 0 && row.resUnbiasedSigma0 > 0 &&
        row.resUnbiasedSigma1 > 0;
    const bool predictionOk =
        residualOk && row.hasPredictionBranches &&
        fwdQaValid(row.trackPred0) && fwdQaValid(row.trackPred1);
    const bool closureOk =
        measurementOk && row.hasClosureBranches &&
        fwdQaValid(row.fstClosureU) && fwdQaValid(row.fstClosureV) &&
        fwdQaValid(row.fstClosureZ) && fwdQaValid(row.fstClosureMag);

    hDetId->Fill(row.detId);
    if (fstRow) {
      ++nFstRows;
      hDim->Fill(row.measurementDim, row.residualDim);
      hHasResidual->Fill(row.hasResidual);
      if (row.hasTrackTypeBranch)
        hTrackTypeFstRows->Fill(row.trackType);
    }

    double r = 0.0;
    double phiRad = 0.0;
    double phiDeg = 0.0;
    double rphi = 0.0;
    if (measurementOk) {
      ++nMeasurementRows;
      fwdQaWedgePolar(row.meas0, row.meas1, r, phiRad, phiDeg, rphi);
      hUV->Fill(row.meas0, row.meas1);
      hRPhiVsR->Fill(r, rphi);
      hRPhiVsPhi->Fill(phiDeg, rphi);
      hUVDisk[row.fstDisk]->Fill(row.meas0, row.meas1);
      hRPhiVsRDisk[row.fstDisk]->Fill(r, rphi);
      hRPhiVsPhiDisk[row.fstDisk]->Fill(phiDeg, rphi);
      hUVGlobalWedge[globalWedge]->Fill(row.meas0, row.meas1);
      hRPhiVsRGlobalWedge[globalWedge]->Fill(r, rphi);
    }

    if (closureOk) {
      ++nClosureRows;
      const double closureUUm = row.fstClosureU * kFwdQaCmToMicron;
      const double closureVUm = row.fstClosureV * kFwdQaCmToMicron;
      const double closureZUm = row.fstClosureZ * kFwdQaCmToMicron;
      const double closureMagUm = row.fstClosureMag * kFwdQaCmToMicron;
      hClosureU->Fill(closureUUm);
      hClosureV->Fill(closureVUm);
      hClosureZ->Fill(closureZUm);
      hClosureMag->Fill(closureMagUm);
      pClosureUByGlobalWedge->Fill(globalWedge, closureUUm);
      pClosureVByGlobalWedge->Fill(globalWedge, closureVUm);
      pClosureUVsR->Fill(r, closureUUm);
      pClosureVVsR->Fill(r, closureVUm);
      pClosureUVsPhi->Fill(phiDeg, closureUUm);
      pClosureVVsPhi->Fill(phiDeg, closureVUm);
    }

    if (!residualOk)
      continue;

    ++nSelectedResidualRows;
    if (row.hasTrackTypeBranch)
      hTrackTypeResidualRows->Fill(row.trackType);

    const double resUUm = row.resUnbiased0 * kFwdQaCmToMicron;
    const double resVUm = row.resUnbiased1 * kFwdQaCmToMicron;
    const double biasedUUm = row.resBiased0 * kFwdQaCmToMicron;
    const double biasedVUm = row.resBiased1 * kFwdQaCmToMicron;
    double resRUm = 0.0;
    double resRPhiUm = 0.0;
    double biasedRUm = 0.0;
    double biasedRPhiUm = 0.0;
    fwdQaPolarResidual(phiRad, resUUm, resVUm, resRUm, resRPhiUm);
    fwdQaPolarResidual(phiRad, biasedUUm, biasedVUm, biasedRUm,
                       biasedRPhiUm);

    hResidualU->Fill(resUUm);
    hResidualV->Fill(resVUm);
    hResidualR->Fill(resRUm);
    hResidualRPhi->Fill(resRPhiUm);
    hBiasedU->Fill(biasedUUm);
    hBiasedV->Fill(biasedVUm);
    hBiasedR->Fill(biasedRUm);
    hBiasedRPhi->Fill(biasedRPhiUm);
    hGlobalWedge->Fill(globalWedge);
    hWedgeDisk->Fill(row.fstDisk, row.fstWedge);

    pResUByGlobalWedge->Fill(globalWedge, resUUm);
    pResVByGlobalWedge->Fill(globalWedge, resVUm);
    pResRByGlobalWedge->Fill(globalWedge, resRUm);
    pResRPhiByGlobalWedge->Fill(globalWedge, resRPhiUm);
    hResUByGlobalWedge->Fill(globalWedge, resUUm);
    hResVByGlobalWedge->Fill(globalWedge, resVUm);

    hResUVsU->Fill(row.meas0, resUUm);
    hResUVsV->Fill(row.meas1, resUUm);
    hResVVsU->Fill(row.meas0, resVUm);
    hResVVsV->Fill(row.meas1, resVUm);
    hResRVsR->Fill(r, resRUm);
    hResRPhiVsR->Fill(r, resRPhiUm);
    hResRVsPhi->Fill(phiDeg, resRUm);
    hResRPhiVsPhi->Fill(phiDeg, resRPhiUm);
    pResRVsR->Fill(r, resRUm);
    pResRPhiVsR->Fill(r, resRPhiUm);
    pResRVsPhi->Fill(phiDeg, resRUm);
    pResRPhiVsPhi->Fill(phiDeg, resRPhiUm);
    pResRVsRDisk[row.fstDisk]->Fill(r, resRUm);
    pResRPhiVsRDisk[row.fstDisk]->Fill(r, resRPhiUm);
    pResRVsPhiDisk[row.fstDisk]->Fill(phiDeg, resRUm);
    pResRPhiVsPhiDisk[row.fstDisk]->Fill(phiDeg, resRPhiUm);
    pResUMapGlobalWedge[globalWedge]->Fill(r, rphi, resUUm);
    pResVMapGlobalWedge[globalWedge]->Fill(r, rphi, resVUm);

    if (pullOk) {
      hPullU->Fill(row.pullUnbiased0);
      hPullV->Fill(row.pullUnbiased1);
      hPullBiasedU->Fill(row.pullBiased0);
      hPullBiasedV->Fill(row.pullBiased1);
      pPullUByGlobalWedge->Fill(globalWedge, row.pullUnbiased0);
      pPullVByGlobalWedge->Fill(globalWedge, row.pullUnbiased1);
    }

    if (predictionOk) {
      double predR = 0.0;
      double predPhiRad = 0.0;
      double predPhiDeg = 0.0;
      double predRPhi = 0.0;
      fwdQaWedgePolar(row.trackPred0, row.trackPred1, predR, predPhiRad,
                      predPhiDeg, predRPhi);
      hTrackPredUVsMeasU->Fill(row.meas0, row.trackPred0);
      hTrackPredVVsMeasV->Fill(row.meas1, row.trackPred1);
      hTrackPredRVsMeasR->Fill(r, predR);
      hTrackPredRPhiVsMeasRPhi->Fill(rphi, predRPhi);
    }

    if (row.hasTrackSlopeBranches && row.trackPz != 0.0f) {
      const double slope = row.trackPt / row.trackPz;
      pResRVsEta->Fill(row.trackEta, resRUm);
      pResRPhiVsEta->Fill(row.trackEta, resRPhiUm);
      pResRVsSlope->Fill(slope, resRUm);
      pResRPhiVsSlope->Fill(slope, resRPhiUm);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Terminal summary.

  std::cout << "Input file: " << inputFilename << std::endl;
  std::cout << "fwdAlign entries: " << nEntries << std::endl;
  std::cout << "FST rows with valid disk/wedge: " << nFstRows << std::endl;
  std::cout << "FST measurement rows: " << nMeasurementRows << std::endl;
  std::cout << "FST closure rows: " << nClosureRows << std::endl;
  std::cout << "Selected FST residual rows: " << nSelectedResidualRows
            << std::endl;
  std::cout << "Fit cut: " << fitCut.Data() << std::endl;
  std::cout << "Residual/pull track cut: " << trackCut.Data() << std::endl;
  std::cout << "Global wedge definition: 12*fstDisk + fstWedge (0..35)"
            << std::endl;
  std::cout << "Wedge polar definition: r=hypot(meas0,meas1), "
               "phi=atan2(meas1,meas0), rphi=r*phi"
            << std::endl;
  std::cout << "Residual units: microns (tree stores cm)" << std::endl;
  std::cout << "Pulls remain in wedge U/V; polar pulls are not calculated "
               "without the full covariance."
            << std::endl;
  std::cout << "Pull branches: "
            << (row.hasPullBranches ? "available" : "missing") << std::endl;
  std::cout << "Prediction branches: "
            << (row.hasPredictionBranches ? "available" : "missing")
            << std::endl;
  std::cout << "Closure branches: "
            << (row.hasClosureBranches ? "available" : "missing")
            << std::endl;
  if (kFwdQaRequiredTrackType >= 0 && !row.hasTrackTypeBranch)
    std::cout << "trackType cut requested but branch is missing; cut skipped."
              << std::endl;
  if (minTrackNHitsFit > 0 && !row.hasTrackNHitsFitBranch)
    std::cout << "trackNHitsFit cut requested but branch is missing; cut "
                 "skipped."
              << std::endl;
  if (minTrackNFstHits > 0 && !row.hasTrackNFstHitsBranch)
    std::cout << "trackNFstHits cut requested but branch is missing; cut "
                 "skipped."
              << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  // PDF pages.

  TCanvas *canvas =
      new TCanvas("cFwdAlignQa", "Wedge-based FST alignment QA", 1400, 950);
  canvas->Print(pdfOutput + "[");

  canvas->Clear();
  TPaveText *summary = new TPaveText(0.06, 0.08, 0.94, 0.92, "NDC");
  summary->SetFillColor(0);
  summary->SetTextAlign(12);
  summary->AddText("Forward STAR FST wedge-based alignment QA");
  summary->AddText(TString::Format("Input: %s", inputFilename));
  summary->AddText(TString::Format("Selected residual rows: %lld",
                                   nSelectedResidualRows));
  summary->AddText(TString::Format("Fit cut: %s", fitCut.Data()));
  summary->AddText(TString::Format("Track cut: %s", trackCut.Data()));
  summary->AddText("Alignment unit: global wedge = 12*fstDisk + fstWedge");
  summary->AddText("All three readout regions in a wedge share one U/V frame");
  summary->AddText("r = sqrt(U^2+V^2), phi = atan2(V,U), rphi = r*phi");
  summary->AddText("resR = cos(phi)*resU + sin(phi)*resV");
  summary->AddText("resRphi = -sin(phi)*resU + cos(phi)*resV");
  summary->AddText("Residuals and closures are shown in microns");
  summary->AddText("Measurement occupancy is uncut; alignment residual plots "
                   "use all listed cuts");
  summary->AddText("Pulls are kept in U/V because the full covariance is not "
                   "stored");
  summary->AddText("closureZ compares physical FTUS z with the flat FwdHit z "
                   "and need not be zero");
  summary->Draw();
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hDetId->Draw();
  canvas->cd(2);
  hDim->Draw("COLZ TEXT");
  canvas->cd(3);
  hHasResidual->Draw();
  canvas->cd(4);
  hWedgeDisk->Draw("COLZ TEXT");
  canvas->Print(pdfOutput);

  if (row.hasTrackTypeBranch) {
    canvas->Clear();
    canvas->Divide(2, 1);
    canvas->cd(1);
    hTrackTypeFstRows->Draw();
    canvas->cd(2);
    hTrackTypeResidualRows->Draw();
    canvas->Print(pdfOutput);
  }

  canvas->Clear();
  canvas->Divide(3, 1);
  canvas->cd(1);
  hUV->Draw("COLZ");
  canvas->cd(2);
  hRPhiVsR->Draw("COLZ");
  canvas->cd(3);
  hRPhiVsPhi->Draw("COLZ");
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(3, 3);
  for (int disk = 0; disk < kFwdQaNumFstDisks; ++disk) {
    canvas->cd(disk + 1);
    hUVDisk[disk]->Draw("COLZ");
    canvas->cd(disk + 4);
    hRPhiVsRDisk[disk]->Draw("COLZ");
    canvas->cd(disk + 7);
    hRPhiVsPhiDisk[disk]->Draw("COLZ");
  }
  canvas->Print(pdfOutput);

  for (int disk = 0; disk < kFwdQaNumFstDisks; ++disk) {
    canvas->Clear();
    canvas->Divide(4, 3);
    for (int wedge = 0; wedge < kFwdQaNumWedgesPerDisk; ++wedge) {
      const int globalWedge = disk * kFwdQaNumWedgesPerDisk + wedge;
      canvas->cd(wedge + 1);
      hRPhiVsRGlobalWedge[globalWedge]->Draw("COLZ");
    }
    canvas->Print(pdfOutput);
  }

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hResidualU->Draw();
  fwdQaDrawLine(0, 0, 0, hResidualU->GetMaximum());
  canvas->cd(2);
  hResidualV->Draw();
  fwdQaDrawLine(0, 0, 0, hResidualV->GetMaximum());
  canvas->cd(3);
  hResidualR->Draw();
  fwdQaDrawLine(0, 0, 0, hResidualR->GetMaximum());
  canvas->cd(4);
  hResidualRPhi->Draw();
  fwdQaDrawLine(0, 0, 0, hResidualRPhi->GetMaximum());
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hBiasedU->Draw();
  fwdQaDrawLine(0, 0, 0, hBiasedU->GetMaximum());
  canvas->cd(2);
  hBiasedV->Draw();
  fwdQaDrawLine(0, 0, 0, hBiasedV->GetMaximum());
  canvas->cd(3);
  hBiasedR->Draw();
  fwdQaDrawLine(0, 0, 0, hBiasedR->GetMaximum());
  canvas->cd(4);
  hBiasedRPhi->Draw();
  fwdQaDrawLine(0, 0, 0, hBiasedRPhi->GetMaximum());
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  pResUByGlobalWedge->Draw("E1");
  fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
  fwdQaDrawDiskBoundaries(pResUByGlobalWedge->GetMinimum(),
                          pResUByGlobalWedge->GetMaximum());
  canvas->cd(2);
  pResVByGlobalWedge->Draw("E1");
  fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
  fwdQaDrawDiskBoundaries(pResVByGlobalWedge->GetMinimum(),
                          pResVByGlobalWedge->GetMaximum());
  canvas->cd(3);
  pResRByGlobalWedge->Draw("E1");
  fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
  fwdQaDrawDiskBoundaries(pResRByGlobalWedge->GetMinimum(),
                          pResRByGlobalWedge->GetMaximum());
  canvas->cd(4);
  pResRPhiByGlobalWedge->Draw("E1");
  fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
  fwdQaDrawDiskBoundaries(pResRPhiByGlobalWedge->GetMinimum(),
                          pResRPhiByGlobalWedge->GetMaximum());
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(2, 1);
  canvas->cd(1);
  hResUByGlobalWedge->Draw("COLZ");
  fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
  canvas->cd(2);
  hResVByGlobalWedge->Draw("COLZ");
  fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hResUVsU->Draw("COLZ");
  canvas->cd(2);
  hResUVsV->Draw("COLZ");
  canvas->cd(3);
  hResVVsU->Draw("COLZ");
  canvas->cd(4);
  hResVVsV->Draw("COLZ");
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hResRVsR->Draw("COLZ");
  canvas->cd(2);
  hResRPhiVsR->Draw("COLZ");
  canvas->cd(3);
  hResRVsPhi->Draw("COLZ");
  canvas->cd(4);
  hResRPhiVsPhi->Draw("COLZ");
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  pResRVsR->Draw("E1");
  fwdQaDrawLine(kFwdQaRMin, 0.0, kFwdQaRMax, 0.0);
  canvas->cd(2);
  pResRPhiVsR->Draw("E1");
  fwdQaDrawLine(kFwdQaRMin, 0.0, kFwdQaRMax, 0.0);
  canvas->cd(3);
  pResRVsPhi->Draw("E1");
  fwdQaDrawLine(kFwdQaPhiMinDeg, 0.0, kFwdQaPhiMaxDeg, 0.0);
  canvas->cd(4);
  pResRPhiVsPhi->Draw("E1");
  fwdQaDrawLine(kFwdQaPhiMinDeg, 0.0, kFwdQaPhiMaxDeg, 0.0);
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(3, 2);
  for (int disk = 0; disk < kFwdQaNumFstDisks; ++disk) {
    canvas->cd(disk + 1);
    pResRVsRDisk[disk]->Draw("E1");
    fwdQaDrawLine(kFwdQaRMin, 0.0, kFwdQaRMax, 0.0);
    canvas->cd(disk + 4);
    pResRPhiVsRDisk[disk]->Draw("E1");
    fwdQaDrawLine(kFwdQaRMin, 0.0, kFwdQaRMax, 0.0);
  }
  canvas->Print(pdfOutput);

  canvas->Clear();
  canvas->Divide(3, 2);
  for (int disk = 0; disk < kFwdQaNumFstDisks; ++disk) {
    canvas->cd(disk + 1);
    pResRVsPhiDisk[disk]->Draw("E1");
    fwdQaDrawLine(kFwdQaPhiMinDeg, 0.0, kFwdQaPhiMaxDeg, 0.0);
    canvas->cd(disk + 4);
    pResRPhiVsPhiDisk[disk]->Draw("E1");
    fwdQaDrawLine(kFwdQaPhiMinDeg, 0.0, kFwdQaPhiMaxDeg, 0.0);
  }
  canvas->Print(pdfOutput);

  for (int disk = 0; disk < kFwdQaNumFstDisks; ++disk) {
    canvas->Clear();
    canvas->Divide(4, 3);
    for (int wedge = 0; wedge < kFwdQaNumWedgesPerDisk; ++wedge) {
      const int globalWedge = disk * kFwdQaNumWedgesPerDisk + wedge;
      canvas->cd(wedge + 1);
      pResUMapGlobalWedge[globalWedge]->Draw("COLZ");
    }
    canvas->Print(pdfOutput);

    canvas->Clear();
    canvas->Divide(4, 3);
    for (int wedge = 0; wedge < kFwdQaNumWedgesPerDisk; ++wedge) {
      const int globalWedge = disk * kFwdQaNumWedgesPerDisk + wedge;
      canvas->cd(wedge + 1);
      pResVMapGlobalWedge[globalWedge]->Draw("COLZ");
    }
    canvas->Print(pdfOutput);
  }

  if (row.hasPullBranches) {
    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    hPullU->Draw();
    fwdQaDrawLine(0, 0, 0, hPullU->GetMaximum());
    canvas->cd(2);
    hPullV->Draw();
    fwdQaDrawLine(0, 0, 0, hPullV->GetMaximum());
    canvas->cd(3);
    hPullBiasedU->Draw();
    fwdQaDrawLine(0, 0, 0, hPullBiasedU->GetMaximum());
    canvas->cd(4);
    hPullBiasedV->Draw();
    fwdQaDrawLine(0, 0, 0, hPullBiasedV->GetMaximum());
    canvas->Print(pdfOutput);

    canvas->Clear();
    canvas->Divide(2, 1);
    canvas->cd(1);
    pPullUByGlobalWedge->Draw("E1");
    fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
    canvas->cd(2);
    pPullVByGlobalWedge->Draw("E1");
    fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
    canvas->Print(pdfOutput);
  }

  if (row.hasPredictionBranches) {
    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    hTrackPredUVsMeasU->Draw("COLZ");
    fwdQaDrawLine(kFwdQaRMin, kFwdQaRMin, kFwdQaRMax, kFwdQaRMax);
    canvas->cd(2);
    hTrackPredVVsMeasV->Draw("COLZ");
    fwdQaDrawLine(kFwdQaRPhiMin, kFwdQaRPhiMin, kFwdQaRPhiMax,
                  kFwdQaRPhiMax);
    canvas->cd(3);
    hTrackPredRVsMeasR->Draw("COLZ");
    fwdQaDrawLine(kFwdQaRMin, kFwdQaRMin, kFwdQaRMax, kFwdQaRMax);
    canvas->cd(4);
    hTrackPredRPhiVsMeasRPhi->Draw("COLZ");
    fwdQaDrawLine(kFwdQaRPhiMin, kFwdQaRPhiMin, kFwdQaRPhiMax,
                  kFwdQaRPhiMax);
    canvas->Print(pdfOutput);
  }

  if (row.hasClosureBranches) {
    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    hClosureU->Draw();
    fwdQaDrawLine(0, 0, 0, hClosureU->GetMaximum());
    canvas->cd(2);
    hClosureV->Draw();
    fwdQaDrawLine(0, 0, 0, hClosureV->GetMaximum());
    canvas->cd(3);
    hClosureZ->Draw();
    fwdQaDrawLine(0, 0, 0, hClosureZ->GetMaximum());
    canvas->cd(4);
    hClosureMag->Draw();
    canvas->Print(pdfOutput);

    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    pClosureUByGlobalWedge->Draw("E1");
    fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
    canvas->cd(2);
    pClosureVByGlobalWedge->Draw("E1");
    fwdQaDrawLine(-0.5, 0.0, 35.5, 0.0);
    canvas->cd(3);
    pClosureUVsR->Draw("E1");
    fwdQaDrawLine(kFwdQaRMin, 0.0, kFwdQaRMax, 0.0);
    canvas->cd(4);
    pClosureVVsR->Draw("E1");
    fwdQaDrawLine(kFwdQaRMin, 0.0, kFwdQaRMax, 0.0);
    canvas->Print(pdfOutput);

    canvas->Clear();
    canvas->Divide(2, 1);
    canvas->cd(1);
    pClosureUVsPhi->Draw("E1");
    fwdQaDrawLine(kFwdQaPhiMinDeg, 0.0, kFwdQaPhiMaxDeg, 0.0);
    canvas->cd(2);
    pClosureVVsPhi->Draw("E1");
    fwdQaDrawLine(kFwdQaPhiMinDeg, 0.0, kFwdQaPhiMaxDeg, 0.0);
    canvas->Print(pdfOutput);
  }

  if (row.hasTrackSlopeBranches) {
    canvas->Clear();
    canvas->Divide(2, 2);
    canvas->cd(1);
    pResRVsEta->Draw("E1");
    fwdQaDrawLine(kFwdQaTrackEtaMin, 0.0, kFwdQaTrackEtaMax, 0.0);
    canvas->cd(2);
    pResRPhiVsEta->Draw("E1");
    fwdQaDrawLine(kFwdQaTrackEtaMin, 0.0, kFwdQaTrackEtaMax, 0.0);
    canvas->cd(3);
    pResRVsSlope->Draw("E1");
    fwdQaDrawLine(kFwdQaPolarSlopeMin, 0.0, kFwdQaPolarSlopeMax, 0.0);
    canvas->cd(4);
    pResRPhiVsSlope->Draw("E1");
    fwdQaDrawLine(kFwdQaPolarSlopeMin, 0.0, kFwdQaPolarSlopeMax, 0.0);
    canvas->Print(pdfOutput);
  }

  outputFile->Write();
  canvas->Print(pdfOutput + "]");
  outputFile->Close();
  inputFile->Close();

  std::cout << "Wrote " << rootOutput.Data() << std::endl;
  std::cout << "Wrote " << pdfOutput.Data() << std::endl;
}
