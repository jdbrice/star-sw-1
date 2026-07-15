// Toy per-wedge in-plane FST alignment solve.
//
// Run from the repository root after producing a compact fwdAlign tree:
//   root4star -l -b -q \
//     'toy_fst_wedge_uvgamma.C+("align_test.root","toy_fst_wedge_uvgamma")'
//
// Outputs:
//   toy_fst_wedge_uvgamma.root
//   toy_fst_wedge_uvgamma.pdf
//   toy_fst_wedge_uvgamma.txt
//
// One independent (DeltaU, DeltaV, gamma) triplet is fitted for each of the
// 36 FST wedges. All three physical-z surfaces in a wedge contribute to the
// same alignment object. No fGeom.root input is needed because the compact
// tree is already expressed in the wedge-local frame.
//
// Residual convention:
//   resU = measured U - unbiased predicted U
//   resV = measured V - unbiased predicted V
//
// Linearized geometry model:
//   resU = -DeltaU + gamma * Vpred
//   resV = -DeltaV - gamma * Upred
//
// Positive gamma rotates the wedge U axis toward its V axis. The reported
// signs must be verified with an injected-misalignment refit before constants
// are applied to production geometry.
//
// The fit treats rows as independent. Correlations between different hits on
// the same fitted track are not stored in fwdAlign, so the formal parameter
// uncertainties can be optimistic.

#include "TCanvas.h"
#include "TDecompSVD.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMatrixD.h"
#include "TProfile.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVectorD.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// EDIT CUTS HERE. These track cuts match fwd_alignment_residual_qa.C.

const double kWedgeSolveEtaMin = 2.5;
const double kWedgeSolveEtaMax = 4.0;
const double kWedgeSolvePMin = 1.0;          // GeV/c
const double kWedgeSolvePtMin = 0.2;         // GeV/c
const int kWedgeSolveMinHitsFit = 8;
const int kWedgeSolveMinFstHits = 3;
const int kWedgeSolveTrackType = 1;          // -1=all, 0=Global, 1=BLC
const bool kWedgeSolveRequireFullyConverged = true;

// This is an additional solver-only core cut. Set it <= 0 to disable it.
const double kWedgeSolveMaxAbsResidualCm = 0.5;
const int kWedgeSolveMinRows = 20;
const double kWedgeSolveMaxCondition = 1e10;

///////////////////////////////////////////////////////////////////////////////

const int kWedgeSolveNumDisks = 3;
const int kWedgeSolveWedgesPerDisk = 12;
const int kWedgeSolveNumWedges =
    kWedgeSolveNumDisks * kWedgeSolveWedgesPerDisk;
const double kWedgeSolveCmToMicron = 10000.0;
const double kWedgeSolveRadToMrad = 1000.0;

struct WedgeSolveRow {
  int disk;
  int wedge;
  int trackType;
  int nHitsFit;
  int nFstHits;
  int fullyConverged;

  float trackP;
  float trackPt;
  float trackEta;
  float measU;
  float measV;
  float resU;
  float resV;
  float covUU;
  float covUV;
  float covVV;
};

struct WedgeSolveHit {
  int globalWedge;
  double predU;
  double predV;
  double resU;
  double resV;
};

struct WedgeSolveResult {
  bool valid;
  int nRows;
  int ndf;
  double deltaU;
  double deltaV;
  double gamma;
  double deltaUErr;
  double deltaVErr;
  double gammaErr;
  double chi2;
  double condition;

  WedgeSolveResult()
      : valid(false), nRows(0), ndf(0),
        deltaU(std::numeric_limits<double>::quiet_NaN()),
        deltaV(std::numeric_limits<double>::quiet_NaN()),
        gamma(std::numeric_limits<double>::quiet_NaN()),
        deltaUErr(std::numeric_limits<double>::quiet_NaN()),
        deltaVErr(std::numeric_limits<double>::quiet_NaN()),
        gammaErr(std::numeric_limits<double>::quiet_NaN()),
        chi2(std::numeric_limits<double>::quiet_NaN()),
        condition(std::numeric_limits<double>::quiet_NaN()) {}
};

bool wedgeSolveBind(TTree *tree, const char *name, void *address) {
  if (!tree || !tree->GetBranch(name)) {
    std::cerr << "Missing required fwdAlign branch: " << name << std::endl;
    return false;
  }
  tree->SetBranchAddress(name, address);
  return true;
}

bool wedgeSolveBindAll(TTree *tree, WedgeSolveRow &r) {
  bool ok = true;
  ok &= wedgeSolveBind(tree, "disk", &r.disk);
  ok &= wedgeSolveBind(tree, "wedge", &r.wedge);
  ok &= wedgeSolveBind(tree, "trackType", &r.trackType);
  ok &= wedgeSolveBind(tree, "nHitsFit", &r.nHitsFit);
  ok &= wedgeSolveBind(tree, "nFstHits", &r.nFstHits);
  ok &= wedgeSolveBind(tree, "fullyConverged", &r.fullyConverged);
  ok &= wedgeSolveBind(tree, "trackP", &r.trackP);
  ok &= wedgeSolveBind(tree, "trackPt", &r.trackPt);
  ok &= wedgeSolveBind(tree, "trackEta", &r.trackEta);
  ok &= wedgeSolveBind(tree, "measU", &r.measU);
  ok &= wedgeSolveBind(tree, "measV", &r.measV);
  ok &= wedgeSolveBind(tree, "resU", &r.resU);
  ok &= wedgeSolveBind(tree, "resV", &r.resV);
  ok &= wedgeSolveBind(tree, "covUU", &r.covUU);
  ok &= wedgeSolveBind(tree, "covUV", &r.covUV);
  ok &= wedgeSolveBind(tree, "covVV", &r.covVV);
  return ok;
}

bool wedgeSolvePassCuts(const WedgeSolveRow &r) {
  if (kWedgeSolveRequireFullyConverged && !r.fullyConverged)
    return false;
  if (kWedgeSolveTrackType >= 0 && r.trackType != kWedgeSolveTrackType)
    return false;
  if (r.trackEta < kWedgeSolveEtaMin || r.trackEta > kWedgeSolveEtaMax)
    return false;
  if (r.trackP < kWedgeSolvePMin || r.trackPt < kWedgeSolvePtMin)
    return false;
  if (r.nHitsFit < kWedgeSolveMinHitsFit ||
      r.nFstHits < kWedgeSolveMinFstHits)
    return false;
  if (kWedgeSolveMaxAbsResidualCm > 0 &&
      (std::abs(r.resU) > kWedgeSolveMaxAbsResidualCm ||
       std::abs(r.resV) > kWedgeSolveMaxAbsResidualCm))
    return false;
  return true;
}

struct WedgeAccumulator {
  int nRows;
  double normal[3][3];
  double rhs[3];
  double residualWeightResidual;

  WedgeAccumulator() : nRows(0), residualWeightResidual(0) {
    for (int i = 0; i < 3; ++i) {
      rhs[i] = 0;
      for (int j = 0; j < 3; ++j)
        normal[i][j] = 0;
    }
  }

  bool add(const WedgeSolveRow &r) {
    const double determinant = r.covUU * r.covVV - r.covUV * r.covUV;
    if (!std::isfinite(determinant) || r.covUU <= 0 || r.covVV <= 0 ||
        determinant <= 0)
      return false;

    const double weight[2][2] = {
        {r.covVV / determinant, -r.covUV / determinant},
        {-r.covUV / determinant, r.covUU / determinant}};
    const double predictedU = r.measU - r.resU;
    const double predictedV = r.measV - r.resV;
    const double design[2][3] = {
        {-1.0, 0.0, predictedV},
        {0.0, -1.0, -predictedU}};
    const double residual[2] = {r.resU, r.resV};

    for (int a = 0; a < 3; ++a) {
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          rhs[a] += design[i][a] * weight[i][j] * residual[j];
          for (int b = 0; b < 3; ++b) {
            normal[a][b] +=
                design[i][a] * weight[i][j] * design[j][b];
          }
        }
      }
    }

    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j)
        residualWeightResidual += residual[i] * weight[i][j] * residual[j];
    }
    ++nRows;
    return true;
  }

  WedgeSolveResult solve() const {
    WedgeSolveResult result;
    result.nRows = nRows;
    result.ndf = 2 * nRows - 3;
    if (nRows < kWedgeSolveMinRows)
      return result;

    TMatrixD normalMatrix(3, 3);
    TVectorD rhsVector(3);
    for (int i = 0; i < 3; ++i) {
      rhsVector[i] = rhs[i];
      for (int j = 0; j < 3; ++j)
        normalMatrix(i, j) = normal[i][j];
    }

    TDecompSVD decomposition(normalMatrix);
    Bool_t solved = kFALSE;
    const TVectorD parameter = decomposition.Solve(rhsVector, solved);
    result.condition = decomposition.Condition();
    if (!solved || !std::isfinite(result.condition) ||
        result.condition > kWedgeSolveMaxCondition)
      return result;

    Bool_t inverted = kFALSE;
    const TMatrixD covariance = decomposition.Invert(inverted);
    if (!inverted)
      return result;

    double parameterRhs = 0;
    for (int i = 0; i < 3; ++i)
      parameterRhs += parameter[i] * rhs[i];

    result.valid = true;
    result.deltaU = parameter[0];
    result.deltaV = parameter[1];
    result.gamma = parameter[2];
    result.deltaUErr = covariance[0][0] > 0
                           ? std::sqrt(covariance[0][0])
                           : std::numeric_limits<double>::quiet_NaN();
    result.deltaVErr = covariance[1][1] > 0
                           ? std::sqrt(covariance[1][1])
                           : std::numeric_limits<double>::quiet_NaN();
    result.gammaErr = covariance[2][2] > 0
                          ? std::sqrt(covariance[2][2])
                          : std::numeric_limits<double>::quiet_NaN();
    result.chi2 = residualWeightResidual - parameterRhs;
    return result;
  }
};

void wedgeSolveDrawZero(double xMin, double xMax) {
  TLine *line = new TLine(xMin, 0, xMax, 0);
  line->SetLineColor(kRed + 1);
  line->SetLineStyle(2);
  line->Draw();
}

void toy_fst_wedge_uvgamma(
    const char *inputFilename = "align_test.root",
    const char *outputPrefix = "toy_fst_wedge_uvgamma") {
  gStyle->SetOptStat(1110);

  TFile *inputFile = TFile::Open(inputFilename, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Cannot open " << inputFilename << std::endl;
    return;
  }
  TTree *inputTree = dynamic_cast<TTree *>(inputFile->Get("fwdAlign"));
  if (!inputTree) {
    std::cerr << "Cannot find TTree fwdAlign in " << inputFilename << std::endl;
    return;
  }

  WedgeSolveRow row;
  if (!wedgeSolveBindAll(inputTree, row)) {
    std::cerr << "This solver requires the compact planar fwdAlign schema."
              << std::endl;
    return;
  }

  WedgeAccumulator accumulators[kWedgeSolveNumWedges];
  std::vector<WedgeSolveHit> selectedHits;
  Long64_t selectedRows = 0;
  Long64_t rejectedCovariance = 0;

  const Long64_t nEntries = inputTree->GetEntries();
  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    inputTree->GetEntry(entry);
    if (!wedgeSolvePassCuts(row))
      continue;
    if (row.disk < 0 || row.disk >= kWedgeSolveNumDisks ||
        row.wedge < 0 || row.wedge >= kWedgeSolveWedgesPerDisk)
      continue;

    const int globalWedge =
        row.disk * kWedgeSolveWedgesPerDisk + row.wedge;
    if (!accumulators[globalWedge].add(row)) {
      ++rejectedCovariance;
      continue;
    }

    WedgeSolveHit hit;
    hit.globalWedge = globalWedge;
    hit.predU = row.measU - row.resU;
    hit.predV = row.measV - row.resV;
    hit.resU = row.resU;
    hit.resV = row.resV;
    selectedHits.push_back(hit);
    ++selectedRows;
  }

  WedgeSolveResult results[kWedgeSolveNumWedges];
  for (int globalWedge = 0; globalWedge < kWedgeSolveNumWedges;
       ++globalWedge) {
    results[globalWedge] = accumulators[globalWedge].solve();
  }

  const TString rootName = TString::Format("%s.root", outputPrefix);
  const TString pdfName = TString::Format("%s.pdf", outputPrefix);
  const TString textName = TString::Format("%s.txt", outputPrefix);
  TFile *outputFile = TFile::Open(rootName, "RECREATE");
  outputFile->cd();

  TGraphErrors *gDeltaU = new TGraphErrors();
  gDeltaU->SetName("gDeltaU");
  gDeltaU->SetTitle("Toy per-wedge #DeltaU;12#timesdisk+wedge;#DeltaU [#mum]");
  TGraphErrors *gDeltaV = new TGraphErrors();
  gDeltaV->SetName("gDeltaV");
  gDeltaV->SetTitle("Toy per-wedge #DeltaV;12#timesdisk+wedge;#DeltaV [#mum]");
  TGraphErrors *gGamma = new TGraphErrors();
  gGamma->SetName("gGamma");
  gGamma->SetTitle("Toy per-wedge #gamma;12#timesdisk+wedge;#gamma [mrad]");

  TH1D *hRows = new TH1D(
      "hRows", "Rows used per wedge;12#timesdisk+wedge;rows", 36, -0.5, 35.5);
  TH1D *hChi2Ndf = new TH1D(
      "hChi2Ndf", "Weighted fit #chi^{2}/ndf;12#timesdisk+wedge;#chi^{2}/ndf",
      36, -0.5, 35.5);
  TH1D *hLog10Condition = new TH1D(
      "hLog10Condition", "Normal-matrix condition;12#timesdisk+wedge;"
      "log_{10}(condition)", 36, -0.5, 35.5);
  TProfile *pResUBefore = new TProfile(
      "pResUBefore", "Mean U residual before toy correction;12#timesdisk+wedge;"
      "<resU> [#mum]", 36, -0.5, 35.5);
  TProfile *pResUAfter = new TProfile(
      "pResUAfter", "Mean U residual after algebraic correction;"
      "12#timesdisk+wedge;<resU-modelU> [#mum]", 36, -0.5, 35.5);
  TProfile *pResVBefore = new TProfile(
      "pResVBefore", "Mean V residual before toy correction;12#timesdisk+wedge;"
      "<resV> [#mum]", 36, -0.5, 35.5);
  TProfile *pResVAfter = new TProfile(
      "pResVAfter", "Mean V residual after algebraic correction;"
      "12#timesdisk+wedge;<resV-modelV> [#mum]", 36, -0.5, 35.5);

  int graphPoint = 0;
  for (int globalWedge = 0; globalWedge < kWedgeSolveNumWedges;
       ++globalWedge) {
    const WedgeSolveResult &result = results[globalWedge];
    hRows->SetBinContent(globalWedge + 1, result.nRows);
    if (std::isfinite(result.condition) && result.condition > 0)
      hLog10Condition->SetBinContent(globalWedge + 1,
                                     std::log10(result.condition));
    if (!result.valid)
      continue;
    hChi2Ndf->SetBinContent(
        globalWedge + 1, result.ndf > 0 ? result.chi2 / result.ndf : 0);
    gDeltaU->SetPoint(graphPoint, globalWedge,
                      result.deltaU * kWedgeSolveCmToMicron);
    gDeltaU->SetPointError(graphPoint, 0,
                           result.deltaUErr * kWedgeSolveCmToMicron);
    gDeltaV->SetPoint(graphPoint, globalWedge,
                      result.deltaV * kWedgeSolveCmToMicron);
    gDeltaV->SetPointError(graphPoint, 0,
                           result.deltaVErr * kWedgeSolveCmToMicron);
    gGamma->SetPoint(graphPoint, globalWedge,
                     result.gamma * kWedgeSolveRadToMrad);
    gGamma->SetPointError(graphPoint, 0,
                          result.gammaErr * kWedgeSolveRadToMrad);
    ++graphPoint;
  }

  for (std::vector<WedgeSolveHit>::const_iterator hit = selectedHits.begin();
       hit != selectedHits.end(); ++hit) {
    const WedgeSolveResult &result = results[hit->globalWedge];
    if (!result.valid)
      continue;
    const double modelU = -result.deltaU + result.gamma * hit->predV;
    const double modelV = -result.deltaV - result.gamma * hit->predU;
    pResUBefore->Fill(hit->globalWedge,
                      hit->resU * kWedgeSolveCmToMicron);
    pResUAfter->Fill(hit->globalWedge,
                     (hit->resU - modelU) * kWedgeSolveCmToMicron);
    pResVBefore->Fill(hit->globalWedge,
                      hit->resV * kWedgeSolveCmToMicron);
    pResVAfter->Fill(hit->globalWedge,
                     (hit->resV - modelV) * kWedgeSolveCmToMicron);
  }

  int outDisk = 0;
  int outWedge = 0;
  int outGlobalWedge = 0;
  int outRows = 0;
  int outNdf = 0;
  int outValid = 0;
  double outDeltaU = 0;
  double outDeltaV = 0;
  double outGamma = 0;
  double outDeltaUErr = 0;
  double outDeltaVErr = 0;
  double outGammaErr = 0;
  double outChi2 = 0;
  double outChi2Ndf = 0;
  double outCondition = 0;

  TTree outputTree("wedgeAlignment", "Toy per-wedge FST U/V/gamma solve");
  outputTree.Branch("disk", &outDisk, "disk/I");
  outputTree.Branch("wedge", &outWedge, "wedge/I");
  outputTree.Branch("globalWedge", &outGlobalWedge, "globalWedge/I");
  outputTree.Branch("nRows", &outRows, "nRows/I");
  outputTree.Branch("ndf", &outNdf, "ndf/I");
  outputTree.Branch("valid", &outValid, "valid/I");
  outputTree.Branch("deltaU", &outDeltaU, "deltaU/D");
  outputTree.Branch("deltaV", &outDeltaV, "deltaV/D");
  outputTree.Branch("gamma", &outGamma, "gamma/D");
  outputTree.Branch("deltaUErr", &outDeltaUErr, "deltaUErr/D");
  outputTree.Branch("deltaVErr", &outDeltaVErr, "deltaVErr/D");
  outputTree.Branch("gammaErr", &outGammaErr, "gammaErr/D");
  outputTree.Branch("chi2", &outChi2, "chi2/D");
  outputTree.Branch("chi2Ndf", &outChi2Ndf, "chi2Ndf/D");
  outputTree.Branch("condition", &outCondition, "condition/D");

  std::ofstream textOutput(textName.Data());
  textOutput << "# Toy per-wedge FST DeltaU/DeltaV/gamma solve\n";
  textOutput << "# resU=-DeltaU+gamma*Vpred, "
                "resV=-DeltaV-gamma*Upred\n";
  textOutput << "# Algebraic closure is not a tracking-refit closure.\n";
  textOutput << "# Inter-hit track correlations are ignored; formal errors "
                "may be optimistic.\n";
  textOutput << "# disk wedge globalWedge rows valid deltaU_cm deltaU_um "
                "deltaUErr_cm deltaV_cm deltaV_um deltaVErr_cm gamma_rad "
                "gamma_mrad gammaErr_rad chi2 ndf chi2Ndf condition\n";
  textOutput << std::setprecision(10);

  for (outGlobalWedge = 0; outGlobalWedge < kWedgeSolveNumWedges;
       ++outGlobalWedge) {
    const WedgeSolveResult &result = results[outGlobalWedge];
    outDisk = outGlobalWedge / kWedgeSolveWedgesPerDisk;
    outWedge = outGlobalWedge % kWedgeSolveWedgesPerDisk;
    outRows = result.nRows;
    outNdf = result.ndf;
    outValid = result.valid ? 1 : 0;
    outDeltaU = result.deltaU;
    outDeltaV = result.deltaV;
    outGamma = result.gamma;
    outDeltaUErr = result.deltaUErr;
    outDeltaVErr = result.deltaVErr;
    outGammaErr = result.gammaErr;
    outChi2 = result.chi2;
    outChi2Ndf = result.valid && result.ndf > 0
                     ? result.chi2 / result.ndf
                     : std::numeric_limits<double>::quiet_NaN();
    outCondition = result.condition;
    outputTree.Fill();

    textOutput << outDisk << " " << outWedge << " " << outGlobalWedge << " "
               << outRows << " " << outValid << " " << outDeltaU << " "
               << outDeltaU * kWedgeSolveCmToMicron << " " << outDeltaUErr
               << " " << outDeltaV << " "
               << outDeltaV * kWedgeSolveCmToMicron << " " << outDeltaVErr
               << " " << outGamma << " "
               << outGamma * kWedgeSolveRadToMrad << " " << outGammaErr << " "
               << outChi2 << " " << outNdf << " " << outChi2Ndf << " "
               << outCondition << "\n";
  }
  textOutput.close();

  TCanvas *canvas = new TCanvas(
      "cToyFstWedgeUVGamma", "Toy per-wedge FST alignment", 1400, 950);
  canvas->Print(pdfName + "[");

  canvas->Clear();
  canvas->Divide(1, 3);
  canvas->cd(1); gDeltaU->SetMarkerStyle(20); gDeltaU->Draw("AP");
  wedgeSolveDrawZero(-0.5, 35.5);
  canvas->cd(2); gDeltaV->SetMarkerStyle(20); gDeltaV->Draw("AP");
  wedgeSolveDrawZero(-0.5, 35.5);
  canvas->cd(3); gGamma->SetMarkerStyle(20); gGamma->Draw("AP");
  wedgeSolveDrawZero(-0.5, 35.5);
  canvas->Print(pdfName);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1); pResUBefore->Draw(); wedgeSolveDrawZero(-0.5, 35.5);
  canvas->cd(2); pResUAfter->Draw(); wedgeSolveDrawZero(-0.5, 35.5);
  canvas->cd(3); pResVBefore->Draw(); wedgeSolveDrawZero(-0.5, 35.5);
  canvas->cd(4); pResVAfter->Draw(); wedgeSolveDrawZero(-0.5, 35.5);
  canvas->Print(pdfName);

  canvas->Clear();
  canvas->Divide(3, 1);
  canvas->cd(1); hRows->Draw("hist text");
  canvas->cd(2); hChi2Ndf->Draw("hist text");
  canvas->cd(3); hLog10Condition->Draw("hist text");
  canvas->Print(pdfName);
  canvas->Print(pdfName + "]");

  outputFile->cd();
  gDeltaU->Write();
  gDeltaV->Write();
  gGamma->Write();
  outputFile->Write();
  outputFile->Close();
  inputFile->Close();

  std::cout << "Input rows: " << nEntries << std::endl;
  std::cout << "Rows used: " << selectedRows << std::endl;
  std::cout << "Rows rejected for invalid covariance: "
            << rejectedCovariance << std::endl;
  std::cout << "Valid wedge solves: " << graphPoint << " / "
            << kWedgeSolveNumWedges << std::endl;
  std::cout << "Wrote " << rootName << ", " << pdfName << ", and "
            << textName << std::endl;
  std::cout << "Reminder: the after plots are algebraic closure only. "
               "Apply constants and rerun tracking for real closure."
            << std::endl;
}
