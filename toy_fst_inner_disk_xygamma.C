// Toy in-plane disk alignment solve for inner FST sensors only.
//
// Usage:
//   root -l -b -q
//   'toy_fst_inner_disk_xygamma.C("align_test.root","fGeom.root")'
//
// This uses only branches already written to fwdAlign, plus fGeom.root to
// recover each sensor plane's active-area origin and U/V axes.  It solves one
// rigid in-plane correction per FST disk:
//
//   residual0 ~= -U dot [(deltaX, deltaY, 0) + gammaZ * (zhat x P)]
//   residual1 ~= -V dot [(deltaX, deltaY, 0) + gammaZ * (zhat x P)]
//
// where P = O + meas0 * U + meas1 * V is the measured hit position using the
// active-area sensor origin O.  This intentionally does not solve deltaZ.

#include "TCanvas.h"
#include "TFile.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoShape.h"
#include "TGeoTube.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector3.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {
constexpr int kFstId = 45;
constexpr int kFstNumDisks = 3;
constexpr int kFstNumWedgePerDisk = 12;
constexpr int kFstNumSensorsPerWedge = 3;
constexpr int kFstNumSensors =
    kFstNumDisks * kFstNumWedgePerDisk * kFstNumSensorsPerWedge;
constexpr double kInvalidAlignValue = -99999.0;
constexpr double kCmToMicron = 10000.0;
constexpr double kRadToMrad = 1000.0;
constexpr double kTrackEtaMin = 2.5;
constexpr double kTrackEtaMax = 4.0;
constexpr double kTrackPMin = 1.0;
constexpr double kTrackPtMin = 0.2;
constexpr int kRequiredTrackType = 1;
constexpr int kMinTrackNHitsFit = 8;
constexpr int kMinTrackNFstHits = 3;

struct SensorGeom {
  bool ok = false;
  TVector3 origin;
  TVector3 u;
  TVector3 v;
};

bool invertSymmetric3x3(const double normal[3][3], double inverse[3][3]) {
  const double a00 = normal[0][0];
  const double a01 = normal[0][1];
  const double a02 = normal[0][2];
  const double a11 = normal[1][1];
  const double a12 = normal[1][2];
  const double a22 = normal[2][2];

  const double det = a00 * (a11 * a22 - a12 * a12) -
                     a01 * (a01 * a22 - a12 * a02) +
                     a02 * (a01 * a12 - a11 * a02);

  if (!std::isfinite(det) || std::abs(det) < 1e-20)
    return false;

  inverse[0][0] = (a11 * a22 - a12 * a12) / det;
  inverse[0][1] = (a02 * a12 - a01 * a22) / det;
  inverse[0][2] = (a01 * a12 - a02 * a11) / det;
  inverse[1][0] = inverse[0][1];
  inverse[1][1] = (a00 * a22 - a02 * a02) / det;
  inverse[1][2] = (a01 * a02 - a00 * a12) / det;
  inverse[2][0] = inverse[0][2];
  inverse[2][1] = inverse[1][2];
  inverse[2][2] = (a00 * a11 - a01 * a01) / det;

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (!std::isfinite(inverse[i][j]))
        return false;
    }
  }
  return true;
}

struct DiskSolveResult {
  bool valid = false;
  double deltaX = std::numeric_limits<double>::quiet_NaN();
  double deltaY = std::numeric_limits<double>::quiet_NaN();
  double gammaZ = std::numeric_limits<double>::quiet_NaN();
  double deltaXErr = std::numeric_limits<double>::quiet_NaN();
  double deltaYErr = std::numeric_limits<double>::quiet_NaN();
  double gammaZErr = std::numeric_limits<double>::quiet_NaN();
  double chi2 = std::numeric_limits<double>::quiet_NaN();
  int ndf = 0;
};

struct RigidAccumulator {
  int nRows = 0;
  int nTerms = 0;
  int nTerms0 = 0;
  int nTerms1 = 0;
  double normal[3][3] = {};
  double rhs[3] = {};
  double sumWeightedResidual2 = 0.0;

  static bool validValue(double value) {
    return std::isfinite(value) &&
           std::abs(value) < 0.9 * std::abs(kInvalidAlignValue);
  }

  bool add(double residual, double sigma, const double a[3], int dim) {
    if (!validValue(residual) || !std::isfinite(sigma) || sigma <= 0.0)
      return false;
    for (int i = 0; i < 3; ++i) {
      if (!std::isfinite(a[i]))
        return false;
    }

    const double weight = 1.0 / (sigma * sigma);
    for (int i = 0; i < 3; ++i) {
      rhs[i] += weight * a[i] * residual;
      for (int j = 0; j < 3; ++j)
        normal[i][j] += weight * a[i] * a[j];
    }
    sumWeightedResidual2 += weight * residual * residual;
    ++nTerms;
    if (dim == 0)
      ++nTerms0;
    if (dim == 1)
      ++nTerms1;
    return true;
  }

  DiskSolveResult solve() const {
    DiskSolveResult result;
    if (nTerms < 3)
      return result;

    double cov[3][3] = {};
    if (!invertSymmetric3x3(normal, cov))
      return result;

    double q[3] = {};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j)
        q[i] += cov[i][j] * rhs[j];
    }

    double qNormalQ = 0.0;
    double qRhs = 0.0;
    for (int i = 0; i < 3; ++i) {
      qRhs += q[i] * rhs[i];
      for (int j = 0; j < 3; ++j)
        qNormalQ += q[i] * normal[i][j] * q[j];
    }

    result.valid = true;
    result.deltaX = q[0];
    result.deltaY = q[1];
    result.gammaZ = q[2];
    result.deltaXErr = cov[0][0] > 0.0
                           ? std::sqrt(cov[0][0])
                           : std::numeric_limits<double>::quiet_NaN();
    result.deltaYErr = cov[1][1] > 0.0
                           ? std::sqrt(cov[1][1])
                           : std::numeric_limits<double>::quiet_NaN();
    result.gammaZErr = cov[2][2] > 0.0
                           ? std::sqrt(cov[2][2])
                           : std::numeric_limits<double>::quiet_NaN();
    result.chi2 = sumWeightedResidual2 - 2.0 * qRhs + qNormalQ;
    result.ndf = nTerms - 3;
    return result;
  }
};

struct SelectedResidual {
  int disk = -1;
  bool use0 = false;
  bool use1 = false;
  double residual0 = kInvalidAlignValue;
  double residual1 = kInvalidAlignValue;
  double a0[3] = {};
  double a1[3] = {};
};

bool hasBranch(TTree *tree, const char *name) {
  return tree && tree->GetBranch(name);
}

int globalSensorIndex(int disk, int wedge, int sensor) {
  return disk * kFstNumWedgePerDisk * kFstNumSensorsPerWedge +
         wedge * kFstNumSensorsPerWedge + sensor;
}

bool buildFstSensorGeometry(const char *geomFilename,
                            std::array<SensorGeom, kFstNumSensors> &sensors) {
  static const int kEventSensorToFtusCopy[kFstNumSensorsPerWedge] = {3, 1, 2};
  static const int kElecToGeantWedge[kFstNumWedgePerDisk] = {
      2, 7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8};
  static const int kElecToGeantWedgeDisk2[kFstNumWedgePerDisk] = {
      7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8, 2};

  TFile geomFile(geomFilename, "READ");
  if (geomFile.IsZombie()) {
    std::cerr << "Cannot open geometry file: " << geomFilename << std::endl;
    return false;
  }

  TGeoManager *geom = dynamic_cast<TGeoManager *>(geomFile.Get("dyson"));
  if (!geom) {
    std::cerr << "Cannot find TGeoManager named dyson in " << geomFilename
              << std::endl;
    return false;
  }

  for (int disk = 0; disk < kFstNumDisks; ++disk) {
    const int fstdCopy = disk + 4;
    const int *wedgeMap =
        (fstdCopy == 5) ? kElecToGeantWedgeDisk2 : kElecToGeantWedge;

    for (int wedge = 0; wedge < kFstNumWedgePerDisk; ++wedge) {
      const int geantFstwCopy = wedgeMap[wedge];

      for (int sensor = 0; sensor < kFstNumSensorsPerWedge; ++sensor) {
        const int ftusCopy = kEventSensorToFtusCopy[sensor];
        const int globalSensor = globalSensorIndex(disk, wedge, sensor);
        TString path =
            TString::Format("/HALL_1/CAVE_1/FSTM_1/FSTD_%d/FSTW_%d/FTUS_%d",
                            fstdCopy, geantFstwCopy, ftusCopy);

        if (!geom->cd(path)) {
          std::cerr << "Cannot cd to " << path.Data() << std::endl;
          continue;
        }

        TGeoMatrix *matrix = geom->GetCurrentMatrix();
        TGeoNode *node = geom->GetCurrentNode();
        if (!matrix || !node || !node->GetVolume())
          continue;

        const Double_t *translation = matrix->GetTranslation();
        TVector3 origin(translation[0], translation[1], translation[2]);

        // FTUS node translations are not the active silicon centers.  This
        // mirrors FwdGeomUtils/TrackFitter: use the active shape center and
        // transform it to global coordinates.
        TGeoShape *shape = node->GetVolume()->GetShape();
        Double_t activeLocal[3] = {0.0, 0.0, 0.0};
        bool haveActiveCenter = false;
        if (auto tube = dynamic_cast<TGeoTubeSeg *>(shape)) {
          const double rCenter = 0.5 * (tube->GetRmin() + tube->GetRmax());
          const double phiCenter =
              0.5 * (tube->GetPhi1() + tube->GetPhi2()) * TMath::DegToRad();
          activeLocal[0] = rCenter * std::cos(phiCenter);
          activeLocal[1] = rCenter * std::sin(phiCenter);
          activeLocal[2] = 0.0;
          haveActiveCenter = true;
        } else if (auto box = dynamic_cast<TGeoBBox *>(shape)) {
          const Double_t *boxOrigin = box->GetOrigin();
          activeLocal[0] = boxOrigin[0];
          activeLocal[1] = boxOrigin[1];
          activeLocal[2] = boxOrigin[2];
          haveActiveCenter = true;
        }
        if (haveActiveCenter) {
          Double_t activeMaster[3] = {0.0, 0.0, 0.0};
          matrix->LocalToMaster(activeLocal, activeMaster);
          origin.SetXYZ(activeMaster[0], activeMaster[1], activeMaster[2]);
        }

        const Double_t *rot = matrix->GetRotationMatrix();
        TVector3 u(rot[0], rot[3], rot[6]);
        TVector3 v(rot[1], rot[4], rot[7]);
        if (u.Cross(v).Z() < 0)
          v = -v;

        sensors[globalSensor].ok = true;
        sensors[globalSensor].origin = origin;
        sensors[globalSensor].u = u.Unit();
        sensors[globalSensor].v = v.Unit();
      }
    }
  }

  return true;
}

bool validAlignValue(double value) {
  return std::isfinite(value) &&
         std::abs(value) < 0.9 * std::abs(kInvalidAlignValue);
}

double dot3(const double a[3], const DiskSolveResult &result) {
  return a[0] * result.deltaX + a[1] * result.deltaY + a[2] * result.gammaZ;
}

void drawZeroLine(double x1, double x2) {
  TLine *line = new TLine(x1, 0.0, x2, 0.0);
  line->SetLineStyle(2);
  line->SetLineColor(kRed + 1);
  line->Draw();
}

void drawVerticalZeroLine(double yMax) {
  TLine *line = new TLine(0.0, 0.0, 0.0, yMax);
  line->SetLineStyle(2);
  line->SetLineColor(kRed + 1);
  line->Draw();
}
} // namespace

void toy_fst_inner_disk_xygamma(
    const char *alignFilename = "align_test.root",
    const char *geomFilename = "fGeom.root",
    const char *outputRoot = "toy_fst_inner_disk_xygamma.root",
    const char *outputPdf = "toy_fst_inner_disk_xygamma.pdf",
    const char *outputTxt = "toy_fst_inner_disk_xygamma.txt") {
  gStyle->SetOptStat(1110);

  std::array<SensorGeom, kFstNumSensors> sensorGeom;
  if (!buildFstSensorGeometry(geomFilename, sensorGeom))
    return;

  TFile alignFile(alignFilename, "READ");
  if (alignFile.IsZombie()) {
    std::cerr << "Cannot open alignment file: " << alignFilename << std::endl;
    return;
  }

  TTree *tree = dynamic_cast<TTree *>(alignFile.Get("fwdAlign"));
  if (!tree) {
    std::cerr << "Cannot find TTree fwdAlign in " << alignFilename << std::endl;
    return;
  }

  const char *requiredBranches[] = {"detId",
                                    "measurementDim",
                                    "fstGlobalSensor",
                                    "fstDisk",
                                    "fstWedge",
                                    "fstSensor",
                                    "hasResidual",
                                    "residualDim",
                                    "fitConverged",
                                    "trackP",
                                    "trackPt",
                                    "trackEta",
                                    "meas0",
                                    "meas1",
                                    "resUnbiased0",
                                    "resUnbiased1",
                                    "resUnbiasedSigma0",
                                    "resUnbiasedSigma1"};
  for (const char *branch : requiredBranches) {
    if (!hasBranch(tree, branch)) {
      std::cerr << "Missing required branch: " << branch << std::endl;
      return;
    }
  }

  int detId = -1;
  int measurementDim = 0;
  int fstGlobalSensor = -1;
  int fstDisk = -1;
  int fstWedge = -1;
  int fstSensor = -1;
  int hasResidual = 0;
  int residualDim = 0;
  int fitConverged = 0;
  int trackType = -1;
  int trackNHitsFit = 0;
  int trackNFstHits = 0;
  float trackP = 0.0f;
  float trackPt = 0.0f;
  float trackEta = 0.0f;
  float meas0 = kInvalidAlignValue;
  float meas1 = kInvalidAlignValue;
  float resUnbiased0 = kInvalidAlignValue;
  float resUnbiased1 = kInvalidAlignValue;
  float resUnbiasedSigma0 = kInvalidAlignValue;
  float resUnbiasedSigma1 = kInvalidAlignValue;

  tree->SetBranchAddress("detId", &detId);
  tree->SetBranchAddress("measurementDim", &measurementDim);
  tree->SetBranchAddress("fstGlobalSensor", &fstGlobalSensor);
  tree->SetBranchAddress("fstDisk", &fstDisk);
  tree->SetBranchAddress("fstWedge", &fstWedge);
  tree->SetBranchAddress("fstSensor", &fstSensor);
  tree->SetBranchAddress("hasResidual", &hasResidual);
  tree->SetBranchAddress("residualDim", &residualDim);
  tree->SetBranchAddress("fitConverged", &fitConverged);
  const bool hasTrackTypeBranch = hasBranch(tree, "trackType");
  const bool hasTrackNHitsFitBranch = hasBranch(tree, "trackNHitsFit");
  const bool hasTrackNFstHitsBranch = hasBranch(tree, "trackNFstHits");
  if (hasTrackTypeBranch)
    tree->SetBranchAddress("trackType", &trackType);
  if (hasTrackNHitsFitBranch)
    tree->SetBranchAddress("trackNHitsFit", &trackNHitsFit);
  if (hasTrackNFstHitsBranch)
    tree->SetBranchAddress("trackNFstHits", &trackNFstHits);
  tree->SetBranchAddress("trackP", &trackP);
  tree->SetBranchAddress("trackPt", &trackPt);
  tree->SetBranchAddress("trackEta", &trackEta);
  tree->SetBranchAddress("meas0", &meas0);
  tree->SetBranchAddress("meas1", &meas1);
  tree->SetBranchAddress("resUnbiased0", &resUnbiased0);
  tree->SetBranchAddress("resUnbiased1", &resUnbiased1);
  tree->SetBranchAddress("resUnbiasedSigma0", &resUnbiasedSigma0);
  tree->SetBranchAddress("resUnbiasedSigma1", &resUnbiasedSigma1);

  std::array<RigidAccumulator, kFstNumDisks> byDisk;
  std::vector<SelectedResidual> selectedResiduals;

  Long64_t nRead = 0;
  Long64_t nSelected = 0;
  Long64_t nNoGeom = 0;
  Long64_t nBadMeasurement = 0;
  Long64_t nNoTerms = 0;

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    ++nRead;

    if (detId != kFstId)
      continue;
    if (measurementDim != 2)
      continue;
    if (fstSensor != 0)
      continue;
    if (fstGlobalSensor < 0 || fstGlobalSensor >= kFstNumSensors)
      continue;
    if (fstDisk < 0 || fstDisk >= kFstNumDisks)
      continue;
    if (hasResidual != 1 || residualDim != 2 || fitConverged != 1)
      continue;
    if (!(trackEta >= kTrackEtaMin && trackEta <= kTrackEtaMax))
      continue;
    if (!(trackP >= kTrackPMin))
      continue;
    if (!(trackPt >= kTrackPtMin))
      continue;
    if (hasTrackTypeBranch && trackType != kRequiredTrackType)
      continue;
    if (hasTrackNHitsFitBranch && trackNHitsFit < kMinTrackNHitsFit)
      continue;
    if (hasTrackNFstHitsBranch && trackNFstHits < kMinTrackNFstHits)
      continue;
    if (!sensorGeom[fstGlobalSensor].ok) {
      ++nNoGeom;
      continue;
    }
    if (!validAlignValue(meas0) || !validAlignValue(meas1)) {
      ++nBadMeasurement;
      continue;
    }

    const SensorGeom &geom = sensorGeom[fstGlobalSensor];
    const TVector3 hitPosition = geom.origin + meas0 * geom.u + meas1 * geom.v;
    const TVector3 zCrossP(-hitPosition.Y(), hitPosition.X(), 0.0);

    SelectedResidual row;
    row.disk = fstDisk;
    row.residual0 = resUnbiased0;
    row.residual1 = resUnbiased1;
    row.a0[0] = -geom.u.X();
    row.a0[1] = -geom.u.Y();
    row.a0[2] = -geom.u.Dot(zCrossP);
    row.a1[0] = -geom.v.X();
    row.a1[1] = -geom.v.Y();
    row.a1[2] = -geom.v.Dot(zCrossP);

    bool usedAnyTerm = false;
    row.use0 = byDisk[fstDisk].add(resUnbiased0, resUnbiasedSigma0, row.a0, 0);
    row.use1 = byDisk[fstDisk].add(resUnbiased1, resUnbiasedSigma1, row.a1, 1);
    usedAnyTerm = row.use0 || row.use1;
    if (!usedAnyTerm) {
      ++nNoTerms;
      continue;
    }

    ++nSelected;
    ++byDisk[fstDisk].nRows;
    selectedResiduals.push_back(row);
  }

  std::array<DiskSolveResult, kFstNumDisks> results;
  for (int disk = 0; disk < kFstNumDisks; ++disk)
    results[disk] = byDisk[disk].solve();

  std::cout << "\nToy inner-FST disk deltaX/deltaY/gammaZ solve" << std::endl;
  std::cout << "  input tree entries: " << nRead << std::endl;
  std::cout << "  selected inner FST rows: " << nSelected << std::endl;
  std::cout << "  skipped no geometry: " << nNoGeom << std::endl;
  std::cout << "  skipped bad meas0/meas1: " << nBadMeasurement << std::endl;
  std::cout << "  skipped no valid residual terms: " << nNoTerms << std::endl;
  std::cout << "  selection: detId==45, measurementDim==2, fstSensor==0, "
               "hasResidual, residualDim==2, fitConverged, "
            << "2.5<=eta<=4.0, P>=1.0, pT>=0.2, trackType==1 if available, "
            << "trackNHitsFit>=8 if available, trackNFstHits>=3 if available"
            << std::endl;
  std::cout << "  origin convention: active FTUS shape center from fGeom.root, "
               "not raw FTUS node translation"
            << std::endl;
  std::cout << "  sign convention: residual = measurement - prediction; output "
               "parameters are geometry corrections"
            << std::endl;
  std::cout << "  model: residual0/1 ~= -axis dot [(deltaX,deltaY,0) + "
               "gammaZ*(zhat x P)]"
            << std::endl;

  for (int disk = 0; disk < kFstNumDisks; ++disk) {
    const DiskSolveResult &r = results[disk];
    if (!r.valid) {
      std::cout << Form("disk %d: no valid solve", disk) << std::endl;
      continue;
    }
    std::cout << Form("disk %d  rows=%6d terms=%6d (u=%5d v=%5d)  deltaX=% .6f "
                      "cm (% .2f um)  deltaY=% .6f cm (% .2f um)  gammaZ=% .6e "
                      "rad (% .3f mrad)  chi2/ndf=% .3f/%d",
                      disk, byDisk[disk].nRows, byDisk[disk].nTerms,
                      byDisk[disk].nTerms0, byDisk[disk].nTerms1, r.deltaX,
                      r.deltaX * kCmToMicron, r.deltaY, r.deltaY * kCmToMicron,
                      r.gammaZ, r.gammaZ * kRadToMrad, r.chi2, r.ndf)
              << std::endl;
  }

  TGraphErrors *gDeltaX = new TGraphErrors();
  gDeltaX->SetName("gDeltaXByDisk");
  gDeltaX->SetTitle("Inner-FST disk #Deltax solve;fstDisk;#Deltax [#mum]");
  TGraphErrors *gDeltaY = new TGraphErrors();
  gDeltaY->SetName("gDeltaYByDisk");
  gDeltaY->SetTitle("Inner-FST disk #Deltay solve;fstDisk;#Deltay [#mum]");
  TGraphErrors *gGammaZ = new TGraphErrors();
  gGammaZ->SetName("gGammaZByDisk");
  gGammaZ->SetTitle(
      "Inner-FST disk #gamma_{Z} solve;fstDisk;#gamma_{Z} [mrad]");
  TGraph *gChi2 = new TGraph();
  gChi2->SetName("gChi2PerNdfByDisk");
  gChi2->SetTitle("Inner-FST disk solve #chi^{2}/ndf;fstDisk;#chi^{2}/ndf");

  for (int disk = 0; disk < kFstNumDisks; ++disk) {
    const DiskSolveResult &r = results[disk];
    if (!r.valid)
      continue;
    int ip = gDeltaX->GetN();
    gDeltaX->SetPoint(ip, disk, r.deltaX * kCmToMicron);
    gDeltaX->SetPointError(ip, 0.0, r.deltaXErr * kCmToMicron);
    ip = gDeltaY->GetN();
    gDeltaY->SetPoint(ip, disk, r.deltaY * kCmToMicron);
    gDeltaY->SetPointError(ip, 0.0, r.deltaYErr * kCmToMicron);
    ip = gGammaZ->GetN();
    gGammaZ->SetPoint(ip, disk, r.gammaZ * kRadToMrad);
    gGammaZ->SetPointError(ip, 0.0, r.gammaZErr * kRadToMrad);
    if (r.ndf > 0) {
      ip = gChi2->GetN();
      gChi2->SetPoint(ip, disk, r.chi2 / r.ndf);
    }
  }

  TH1D *hRes0Before = new TH1D(
      "hRes0Before",
      "Inner FST residual 0 before disk correction;resUnbiased0 [#mum];terms",
      160, -2000, 2000);
  TH1D *hRes1Before = new TH1D(
      "hRes1Before",
      "Inner FST residual 1 before disk correction;resUnbiased1 [#mum];terms",
      160, -2000, 2000);
  TH1D *hRes0After = new TH1D("hRes0After",
                              "Inner FST residual 0 after disk "
                              "correction;corrected residual0 [#mum];terms",
                              160, -2000, 2000);
  TH1D *hRes1After = new TH1D("hRes1After",
                              "Inner FST residual 1 after disk "
                              "correction;corrected residual1 [#mum];terms",
                              160, -2000, 2000);

  for (const SelectedResidual &row : selectedResiduals) {
    if (row.disk < 0 || row.disk >= kFstNumDisks || !results[row.disk].valid)
      continue;
    const DiskSolveResult &r = results[row.disk];
    if (row.use0) {
      hRes0Before->Fill(row.residual0 * kCmToMicron);
      hRes0After->Fill((row.residual0 - dot3(row.a0, r)) * kCmToMicron);
    }
    if (row.use1) {
      hRes1Before->Fill(row.residual1 * kCmToMicron);
      hRes1After->Fill((row.residual1 - dot3(row.a1, r)) * kCmToMicron);
    }
  }

  TFile out(outputRoot, "RECREATE");
  TTree diskTree("diskXYGamma", "Toy FST disk deltaX/deltaY/gammaZ alignment "
                                "parameters from inner sensors");
  int outDisk = -1;
  int outRows = 0;
  int outTerms = 0;
  int outTerms0 = 0;
  int outTerms1 = 0;
  int outNdf = 0;
  int outValid = 0;
  double outDeltaX = std::numeric_limits<double>::quiet_NaN();
  double outDeltaY = std::numeric_limits<double>::quiet_NaN();
  double outGammaZ = std::numeric_limits<double>::quiet_NaN();
  double outDeltaXErr = std::numeric_limits<double>::quiet_NaN();
  double outDeltaYErr = std::numeric_limits<double>::quiet_NaN();
  double outGammaZErr = std::numeric_limits<double>::quiet_NaN();
  double outDeltaXMicron = std::numeric_limits<double>::quiet_NaN();
  double outDeltaYMicron = std::numeric_limits<double>::quiet_NaN();
  double outGammaZMrad = std::numeric_limits<double>::quiet_NaN();
  double outChi2 = std::numeric_limits<double>::quiet_NaN();
  double outChi2PerNdf = std::numeric_limits<double>::quiet_NaN();
  diskTree.Branch("fstDisk", &outDisk, "fstDisk/I");
  diskTree.Branch("rows", &outRows, "rows/I");
  diskTree.Branch("terms", &outTerms, "terms/I");
  diskTree.Branch("terms0", &outTerms0, "terms0/I");
  diskTree.Branch("terms1", &outTerms1, "terms1/I");
  diskTree.Branch("valid", &outValid, "valid/I");
  diskTree.Branch("deltaX", &outDeltaX, "deltaX/D");
  diskTree.Branch("deltaY", &outDeltaY, "deltaY/D");
  diskTree.Branch("gammaZ", &outGammaZ, "gammaZ/D");
  diskTree.Branch("deltaXErr", &outDeltaXErr, "deltaXErr/D");
  diskTree.Branch("deltaYErr", &outDeltaYErr, "deltaYErr/D");
  diskTree.Branch("gammaZErr", &outGammaZErr, "gammaZErr/D");
  diskTree.Branch("deltaXMicron", &outDeltaXMicron, "deltaXMicron/D");
  diskTree.Branch("deltaYMicron", &outDeltaYMicron, "deltaYMicron/D");
  diskTree.Branch("gammaZMrad", &outGammaZMrad, "gammaZMrad/D");
  diskTree.Branch("chi2", &outChi2, "chi2/D");
  diskTree.Branch("ndf", &outNdf, "ndf/I");
  diskTree.Branch("chi2PerNdf", &outChi2PerNdf, "chi2PerNdf/D");
  for (int disk = 0; disk < kFstNumDisks; ++disk) {
    const DiskSolveResult &r = results[disk];
    outDisk = disk;
    outRows = byDisk[disk].nRows;
    outTerms = byDisk[disk].nTerms;
    outTerms0 = byDisk[disk].nTerms0;
    outTerms1 = byDisk[disk].nTerms1;
    outValid = r.valid ? 1 : 0;
    outDeltaX = r.deltaX;
    outDeltaY = r.deltaY;
    outGammaZ = r.gammaZ;
    outDeltaXErr = r.deltaXErr;
    outDeltaYErr = r.deltaYErr;
    outGammaZErr = r.gammaZErr;
    outDeltaXMicron = r.deltaX * kCmToMicron;
    outDeltaYMicron = r.deltaY * kCmToMicron;
    outGammaZMrad = r.gammaZ * kRadToMrad;
    outChi2 = r.chi2;
    outNdf = r.ndf;
    outChi2PerNdf =
        r.ndf > 0 ? r.chi2 / r.ndf : std::numeric_limits<double>::quiet_NaN();
    diskTree.Fill();
  }
  diskTree.Write();
  hRes0Before->Write();
  hRes1Before->Write();
  hRes0After->Write();
  hRes1After->Write();
  gDeltaX->Write();
  gDeltaY->Write();
  gGammaZ->Write();
  gChi2->Write();
  out.Close();

  TString pdfOutput(outputPdf);
  TCanvas c("cToyFstInnerDiskXYGamma", "Toy inner FST disk x/y/gamma alignment",
            1400, 1000);
  c.Print(pdfOutput + "[");
  c.Divide(2, 2);
  c.cd(1);
  gDeltaX->SetMarkerStyle(20);
  gDeltaX->Draw("AP");
  drawZeroLine(-0.5, 2.5);
  c.cd(2);
  gDeltaY->SetMarkerStyle(20);
  gDeltaY->Draw("AP");
  drawZeroLine(-0.5, 2.5);
  c.cd(3);
  gGammaZ->SetMarkerStyle(20);
  gGammaZ->Draw("AP");
  drawZeroLine(-0.5, 2.5);
  c.cd(4);
  gChi2->SetMarkerStyle(20);
  gChi2->Draw("AP");
  c.Print(pdfOutput);

  TCanvas cResidual("cToyFstInnerDiskXYGammaResiduals",
                    "Toy inner FST disk x/y/gamma residual closure", 1400,
                    1000);
  cResidual.Divide(2, 2);
  cResidual.cd(1);
  hRes0Before->Draw();
  drawVerticalZeroLine(hRes0Before->GetMaximum());
  cResidual.cd(2);
  hRes0After->Draw();
  drawVerticalZeroLine(hRes0After->GetMaximum());
  cResidual.cd(3);
  hRes1Before->Draw();
  drawVerticalZeroLine(hRes1Before->GetMaximum());
  cResidual.cd(4);
  hRes1After->Draw();
  drawVerticalZeroLine(hRes1After->GetMaximum());
  cResidual.Print(pdfOutput);
  cResidual.Print(pdfOutput + "]");

  std::ofstream txt(outputTxt);
  if (txt.is_open()) {
    txt << "# Toy inner-FST disk deltaX/deltaY/gammaZ solve\n";
    txt << "# One rigid in-plane alignment correction per fstDisk, using all "
           "selected fstSensor==0 residual terms on that disk.\n";
    txt << "# This intentionally does not solve deltaZ; no local track slopes "
           "are used.\n";
    txt << "# selection: detId==45 measurementDim==2 fstSensor==0 hasResidual "
           "residualDim==2 fitConverged 2.5<=trackEta<=4 trackP>=1.0 "
           "trackPt>=0.2 trackType==1_if_available "
           "trackNHitsFit>=8_if_available trackNFstHits>=3_if_available\n";
    txt << "# origin: active FTUS shape center from fGeom.root, not raw FTUS "
           "node translation.\n";
    txt << "# residual = measurement - prediction; output parameters are "
           "geometry corrections.\n";
    txt << "# model: P=O+meas0*U+meas1*V; residual0=-U dot "
           "[(dx,dy,0)+gammaZ*(zhat x P)]; residual1=-V dot same.\n";
    txt << "fstDisk rows terms terms0 terms1 valid deltaX_cm deltaX_um "
           "deltaXErr_cm deltaY_cm deltaY_um deltaYErr_cm gammaZ_rad "
           "gammaZ_mrad gammaZErr_rad chi2 ndf chi2PerNdf\n";
    for (int disk = 0; disk < kFstNumDisks; ++disk) {
      const DiskSolveResult &r = results[disk];
      txt << disk << " " << byDisk[disk].nRows << " " << byDisk[disk].nTerms
          << " " << byDisk[disk].nTerms0 << " " << byDisk[disk].nTerms1 << " "
          << (r.valid ? 1 : 0) << " " << r.deltaX << " "
          << r.deltaX * kCmToMicron << " " << r.deltaXErr << " " << r.deltaY
          << " " << r.deltaY * kCmToMicron << " " << r.deltaYErr << " "
          << r.gammaZ << " " << r.gammaZ * kRadToMrad << " " << r.gammaZErr
          << " " << r.chi2 << " " << r.ndf << " "
          << (r.ndf > 0 ? r.chi2 / r.ndf
                        : std::numeric_limits<double>::quiet_NaN())
          << "\n";
    }
    txt.close();
  } else {
    std::cerr << "Could not write " << outputTxt << std::endl;
  }

  std::cout << "Wrote " << outputRoot << std::endl;
  std::cout << "Wrote " << outputPdf << std::endl;
  std::cout << "Wrote " << outputTxt << std::endl;
}
