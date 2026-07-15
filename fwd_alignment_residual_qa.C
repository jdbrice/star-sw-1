// Compact QA for the planar FST alignment tree written by StFwdTrackMaker.
//
// Run from the repository root:
//   root4star -l -b -q \
//     'fwd_alignment_residual_qa.C+("align_test.root","fwd_align_qa")'
//
// Outputs:
//   fwd_align_qa.root
//   fwd_align_qa.pdf
//
// Residual convention:
//   resU = measured U - unbiased track prediction U
//   resV = measured V - unbiased track prediction V
//
// FST alignment objects are wedges. The "surface" branch only identifies the
// three physical z surfaces inside a wedge; it is not a separate alignment unit.

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
// EDIT CUTS HERE

const double kQaEtaMin = 2.5;
const double kQaEtaMax = 4.0;
const double kQaPMin = 1.0;       // GeV/c
const double kQaPtMin = 0.2;      // GeV/c
const int kQaMinHitsFit = 8;
const int kQaMinFstHits = 3;
const int kQaTrackType = 1;       // -1=all, 0=Global, 1=BLC
const bool kQaRequireFullyConverged = true;

///////////////////////////////////////////////////////////////////////////////
// Plot settings

const int kQaNumDisks = 3;
const int kQaWedgesPerDisk = 12;
const int kQaNumWedges = kQaNumDisks * kQaWedgesPerDisk;
const double kQaCmToMicron = 10000.0;

const int kQaLocalUBins = 25;     // 1 cm bins
const double kQaLocalUMin = 5.0;
const double kQaLocalUMax = 30.0;
const int kQaLocalVBins = 16;     // 1 cm bins
const double kQaLocalVMin = -8.0;
const double kQaLocalVMax = 8.0;

///////////////////////////////////////////////////////////////////////////////
// Tree reader

struct FwdAlignRow {
  int run;
  int event;
  int trackIndex;
  int planeId;
  int disk;
  int wedge;
  int surface;
  int trackType;
  int nHitsFit;
  int nFstHits;
  int fullyConverged;

  float chi2Ndf;
  float trackP;
  float trackPt;
  float trackEta;
  float measU;
  float measV;
  float globalX;
  float globalY;
  float globalZ;
  float resU;
  float resV;
  float covUU;
  float covUV;
  float covVV;
  float slopeU;
  float slopeV;
};

bool fwdQaBind(TTree *tree, const char *name, void *address) {
  if (!tree || !tree->GetBranch(name)) {
    std::cerr << "Missing required fwdAlign branch: " << name << std::endl;
    return false;
  }
  tree->SetBranchAddress(name, address);
  return true;
}

bool fwdQaBindAll(TTree *tree, FwdAlignRow &r) {
  bool ok = true;
  ok &= fwdQaBind(tree, "run", &r.run);
  ok &= fwdQaBind(tree, "event", &r.event);
  ok &= fwdQaBind(tree, "trackIndex", &r.trackIndex);
  ok &= fwdQaBind(tree, "planeId", &r.planeId);
  ok &= fwdQaBind(tree, "disk", &r.disk);
  ok &= fwdQaBind(tree, "wedge", &r.wedge);
  ok &= fwdQaBind(tree, "surface", &r.surface);
  ok &= fwdQaBind(tree, "trackType", &r.trackType);
  ok &= fwdQaBind(tree, "nHitsFit", &r.nHitsFit);
  ok &= fwdQaBind(tree, "nFstHits", &r.nFstHits);
  ok &= fwdQaBind(tree, "fullyConverged", &r.fullyConverged);
  ok &= fwdQaBind(tree, "chi2Ndf", &r.chi2Ndf);
  ok &= fwdQaBind(tree, "trackP", &r.trackP);
  ok &= fwdQaBind(tree, "trackPt", &r.trackPt);
  ok &= fwdQaBind(tree, "trackEta", &r.trackEta);
  ok &= fwdQaBind(tree, "measU", &r.measU);
  ok &= fwdQaBind(tree, "measV", &r.measV);
  ok &= fwdQaBind(tree, "globalX", &r.globalX);
  ok &= fwdQaBind(tree, "globalY", &r.globalY);
  ok &= fwdQaBind(tree, "globalZ", &r.globalZ);
  ok &= fwdQaBind(tree, "resU", &r.resU);
  ok &= fwdQaBind(tree, "resV", &r.resV);
  ok &= fwdQaBind(tree, "covUU", &r.covUU);
  ok &= fwdQaBind(tree, "covUV", &r.covUV);
  ok &= fwdQaBind(tree, "covVV", &r.covVV);
  ok &= fwdQaBind(tree, "slopeU", &r.slopeU);
  ok &= fwdQaBind(tree, "slopeV", &r.slopeV);
  return ok;
}

bool fwdQaPassCuts(const FwdAlignRow &r) {
  if (kQaRequireFullyConverged && !r.fullyConverged)
    return false;
  if (kQaTrackType >= 0 && r.trackType != kQaTrackType)
    return false;
  if (r.trackEta < kQaEtaMin || r.trackEta > kQaEtaMax)
    return false;
  if (r.trackP < kQaPMin || r.trackPt < kQaPtMin)
    return false;
  if (r.nHitsFit < kQaMinHitsFit || r.nFstHits < kQaMinFstHits)
    return false;
  return true;
}

int fwdQaGlobalWedge(const FwdAlignRow &r) {
  if (r.disk < 0 || r.disk >= kQaNumDisks ||
      r.wedge < 0 || r.wedge >= kQaWedgesPerDisk)
    return -1;
  return r.disk * kQaWedgesPerDisk + r.wedge;
}

double fwdQaWrapPhi(double phi) {
  return std::atan2(std::sin(phi), std::cos(phi));
}

void fwdQaLocalPolar(double u, double v, double &radius, double &phi,
                     double &rphi) {
  radius = std::sqrt(u * u + v * v);
  phi = std::atan2(v, u);
  rphi = radius * phi;
}

void fwdQaPolarResidual(double phi, double resU, double resV, double &resR,
                        double &resT) {
  const double c = std::cos(phi);
  const double s = std::sin(phi);
  resR = c * resU + s * resV;
  resT = -s * resU + c * resV;
}

void fwdQaDrawZero(double xMin, double xMax) {
  TLine *line = new TLine(xMin, 0, xMax, 0);
  line->SetLineColor(kRed + 1);
  line->SetLineStyle(2);
  line->Draw();
}

void fwdQaDrawDiskBoundaries(double yMin, double yMax) {
  for (int disk = 1; disk < kQaNumDisks; ++disk) {
    TLine *line = new TLine(disk * kQaWedgesPerDisk - 0.5, yMin,
                            disk * kQaWedgesPerDisk - 0.5, yMax);
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(2);
    line->Draw();
  }
}

void fwdQaDrawPhiBoundaries(double yMin, double yMax) {
  for (int boundary = 1; boundary < kQaWedgesPerDisk; ++boundary) {
    const double phi = -TMath::Pi() + boundary * TMath::TwoPi() /
                                           kQaWedgesPerDisk;
    TLine *line = new TLine(phi, yMin, phi, yMax);
    line->SetLineColor(kGray + 2);
    line->SetLineStyle(3);
    line->Draw();
  }
}

///////////////////////////////////////////////////////////////////////////////

void fwd_alignment_residual_qa(
    const char *inputFilename = "align_test.root",
    const char *outputPrefix = "fwd_align_qa",
    double residualRangeMicron = 5000.0) {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  gStyle->SetPalette(1);

  TFile *inputFile = TFile::Open(inputFilename, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Cannot open " << inputFilename << std::endl;
    return;
  }

  TTree *tree = dynamic_cast<TTree *>(inputFile->Get("fwdAlign"));
  if (!tree) {
    std::cerr << "Cannot find TTree fwdAlign in " << inputFilename << std::endl;
    return;
  }

  FwdAlignRow row;
  if (!fwdQaBindAll(tree, row)) {
    std::cerr << "This file does not use the compact planar tree schema. "
                 "Rerun fwd_afterburner.C with the rebuilt StFwdTrackMaker."
              << std::endl;
    return;
  }

  TString rootName = TString::Format("%s.root", outputPrefix);
  TString pdfName = TString::Format("%s.pdf", outputPrefix);
  TFile *outputFile = TFile::Open(rootName, "RECREATE");
  outputFile->cd();

  // Track-selection QA. These are filled once per track, not once per hit.
  TH1D *hEtaAll = new TH1D("hEtaAll", "All tracks;track #eta;tracks", 60, 2, 5);
  TH1D *hEtaSelected =
      new TH1D("hEtaSelected", "Selected tracks;track #eta;tracks", 60, 2, 5);
  TH1D *hP = new TH1D("hP", "Selected tracks;track p [GeV/c];tracks", 60, 0, 12);
  TH1D *hPt =
      new TH1D("hPt", "Selected tracks;track p_{T} [GeV/c];tracks", 60, 0, 3);
  TH1D *hNHits =
      new TH1D("hNHits", "Selected tracks;fit points;tracks", 16, -0.5, 15.5);
  TH1D *hChi2Ndf =
      new TH1D("hChi2Ndf", "Selected tracks;#chi^{2}/ndf;tracks", 80, 0, 20);

  // Measurement coverage.
  TH2D *hDiskWedge = new TH2D(
      "hDiskWedge", "Selected FST rows;disk;wedge;rows", 3, -0.5, 2.5, 12,
      -0.5, 11.5);
  TH1D *hGlobalWedge = new TH1D(
      "hGlobalWedge", "Selected FST rows;12#timesdisk+wedge;rows", 36, -0.5,
      35.5);
  TH2D *hUV = new TH2D(
      "hUV", "Selected wedge-local measurements;U [cm];V [cm];rows",
      kQaLocalUBins, kQaLocalUMin, kQaLocalUMax, kQaLocalVBins, kQaLocalVMin,
      kQaLocalVMax);
  TH2D *hRPhiVsR = new TH2D(
      "hRPhiVsR", "Selected wedge-local measurements;r [cm];r#phi [cm];rows",
      kQaLocalUBins, kQaLocalUMin, kQaLocalUMax, kQaLocalVBins, kQaLocalVMin,
      kQaLocalVMax);

  // Residual and pull distributions.
  TH1D *hResU = new TH1D("hResU", "Unbiased wedge-U residual;resU [#mum];rows",
                         100, -residualRangeMicron, residualRangeMicron);
  TH1D *hResV = new TH1D("hResV", "Unbiased wedge-V residual;resV [#mum];rows",
                         100, -residualRangeMicron, residualRangeMicron);
  TH1D *hResR = new TH1D("hResR", "Unbiased radial residual;resR [#mum];rows",
                         100, -residualRangeMicron, residualRangeMicron);
  TH1D *hResT = new TH1D(
      "hResT", "Unbiased tangential residual;res(r#Delta#phi) [#mum];rows", 100,
      -residualRangeMicron, residualRangeMicron);
  TH1D *hPullU = new TH1D("hPullU", "Unbiased wedge-U pull;resU/#sigmaU;rows",
                          100, -10, 10);
  TH1D *hPullV = new TH1D("hPullV", "Unbiased wedge-V pull;resV/#sigmaV;rows",
                          100, -10, 10);

  // Wedge means.
  TProfile *pResUByWedge = new TProfile(
      "pResUByWedge", "Mean unbiased U residual by wedge;global wedge;<resU> [#mum]",
      36, -0.5, 35.5);
  TProfile *pResVByWedge = new TProfile(
      "pResVByWedge", "Mean unbiased V residual by wedge;global wedge;<resV> [#mum]",
      36, -0.5, 35.5);
  TProfile *pPullUByWedge = new TProfile(
      "pPullUByWedge", "Mean unbiased U pull by wedge;global wedge;<pullU>", 36,
      -0.5, 35.5);
  TProfile *pPullVByWedge = new TProfile(
      "pPullVByWedge", "Mean unbiased V pull by wedge;global wedge;<pullV>", 36,
      -0.5, 35.5);

  // Dependence on wedge-local hit position.
  TProfile *pResUVsU = new TProfile(
      "pResUVsU", "Mean U residual vs U;U [cm];<resU> [#mum]", kQaLocalUBins,
      kQaLocalUMin, kQaLocalUMax);
  TProfile *pResUVsV = new TProfile(
      "pResUVsV", "Mean U residual vs V;V [cm];<resU> [#mum]", kQaLocalVBins,
      kQaLocalVMin, kQaLocalVMax);
  TProfile *pResVVsU = new TProfile(
      "pResVVsU", "Mean V residual vs U;U [cm];<resV> [#mum]", kQaLocalUBins,
      kQaLocalUMin, kQaLocalUMax);
  TProfile *pResVVsV = new TProfile(
      "pResVVsV", "Mean V residual vs V;V [cm];<resV> [#mum]", kQaLocalVBins,
      kQaLocalVMin, kQaLocalVMax);

  TProfile *pResRVsR = new TProfile(
      "pResRVsR", "Mean radial residual vs r;r [cm];<resR> [#mum]",
      kQaLocalUBins, kQaLocalUMin, kQaLocalUMax);
  TProfile *pResTVsR = new TProfile(
      "pResTVsR", "Mean tangential residual vs r;r [cm];<res(r#Delta#phi)> [#mum]",
      kQaLocalUBins, kQaLocalUMin, kQaLocalUMax);
  TProfile *pResRVsPhi = new TProfile(
      "pResRVsPhi", "Mean radial residual vs local #phi;local #phi [deg];<resR> [#mum]",
      30, -15, 15);
  TProfile *pResTVsPhi = new TProfile(
      "pResTVsPhi", "Mean tangential residual vs local #phi;local #phi [deg];"
                    "<res(r#Delta#phi)> [#mum]",
      30, -15, 15);

  // Exact angular residual used by the STAR residual page.
  TH2D *hRDeltaPhiVsGlobalPhi = new TH2D(
      "hRDeltaPhiVsGlobalPhi",
      "Unbiased FST r#Delta#phi vs global #phi;global #phi [rad];r#Delta#phi [#mum];rows",
      72, -TMath::Pi(), TMath::Pi(), 100, -residualRangeMicron,
      residualRangeMicron);
  TProfile *pRDeltaPhiVsGlobalPhi = new TProfile(
      "pRDeltaPhiVsGlobalPhi",
      "Mean unbiased FST r#Delta#phi vs global #phi;global #phi [rad];"
      "<r#Delta#phi> [#mum]",
      72, -TMath::Pi(), TMath::Pi());
  TProfile *pRDeltaPhiDisk[kQaNumDisks];
  for (int disk = 0; disk < kQaNumDisks; ++disk) {
    pRDeltaPhiDisk[disk] = new TProfile(
        TString::Format("pRDeltaPhiDisk%d", disk),
        TString::Format("Disk %d mean r#Delta#phi vs global #phi;global #phi [rad];"
                        "<r#Delta#phi> [#mum]", disk),
        72, -TMath::Pi(), TMath::Pi());
  }

  // Local track slopes are the quantities needed for out-of-plane alignment.
  TProfile *pResUVsSlopeU = new TProfile(
      "pResUVsSlopeU", "Mean U residual vs du/dw;du/dw;<resU> [#mum]", 40,
      -0.25, 0.25);
  TProfile *pResVVsSlopeV = new TProfile(
      "pResVVsSlopeV", "Mean V residual vs dv/dw;dv/dw;<resV> [#mum]", 40,
      -0.25, 0.25);
  TProfile *pResUVsSlopeV = new TProfile(
      "pResUVsSlopeV", "Mean U residual vs dv/dw;dv/dw;<resU> [#mum]", 40,
      -0.25, 0.25);
  TProfile *pResVVsSlopeU = new TProfile(
      "pResVVsSlopeU", "Mean V residual vs du/dw;du/dw;<resV> [#mum]", 40,
      -0.25, 0.25);

  Long64_t selectedRows = 0;
  Long64_t tracksSeen = 0;
  Long64_t tracksSelected = 0;
  int previousRun = -1;
  int previousEvent = -1;
  int previousTrack = -1;

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    tree->GetEntry(entry);
    const bool newTrack = row.run != previousRun || row.event != previousEvent ||
                          row.trackIndex != previousTrack;
    const bool pass = fwdQaPassCuts(row);

    if (newTrack) {
      ++tracksSeen;
      hEtaAll->Fill(row.trackEta);
      if (pass) {
        ++tracksSelected;
        hEtaSelected->Fill(row.trackEta);
        hP->Fill(row.trackP);
        hPt->Fill(row.trackPt);
        hNHits->Fill(row.nHitsFit);
        hChi2Ndf->Fill(row.chi2Ndf);
      }
      previousRun = row.run;
      previousEvent = row.event;
      previousTrack = row.trackIndex;
    }

    if (!pass)
      continue;
    const int globalWedge = fwdQaGlobalWedge(row);
    if (globalWedge < 0 || row.covUU <= 0 || row.covVV <= 0)
      continue;

    ++selectedRows;
    double radius = 0;
    double localPhi = 0;
    double rphi = 0;
    fwdQaLocalPolar(row.measU, row.measV, radius, localPhi, rphi);

    double resR = 0;
    double resT = 0;
    fwdQaPolarResidual(localPhi, row.resU, row.resV, resR, resT);
    const double resUUm = row.resU * kQaCmToMicron;
    const double resVUm = row.resV * kQaCmToMicron;
    const double resRUm = resR * kQaCmToMicron;
    const double resTUm = resT * kQaCmToMicron;
    const double pullU = row.resU / std::sqrt(row.covUU);
    const double pullV = row.resV / std::sqrt(row.covVV);

    hDiskWedge->Fill(row.disk, row.wedge);
    hGlobalWedge->Fill(globalWedge);
    hUV->Fill(row.measU, row.measV);
    hRPhiVsR->Fill(radius, rphi);
    hResU->Fill(resUUm);
    hResV->Fill(resVUm);
    hResR->Fill(resRUm);
    hResT->Fill(resTUm);
    hPullU->Fill(pullU);
    hPullV->Fill(pullV);

    pResUByWedge->Fill(globalWedge, resUUm);
    pResVByWedge->Fill(globalWedge, resVUm);
    pPullUByWedge->Fill(globalWedge, pullU);
    pPullVByWedge->Fill(globalWedge, pullV);
    pResUVsU->Fill(row.measU, resUUm);
    pResUVsV->Fill(row.measV, resUUm);
    pResVVsU->Fill(row.measU, resVUm);
    pResVVsV->Fill(row.measV, resVUm);
    pResRVsR->Fill(radius, resRUm);
    pResTVsR->Fill(radius, resTUm);
    pResRVsPhi->Fill(localPhi * TMath::RadToDeg(), resRUm);
    pResTVsPhi->Fill(localPhi * TMath::RadToDeg(), resTUm);

    const double predU = row.measU - row.resU;
    const double predV = row.measV - row.resV;
    const double predPhi = std::atan2(predV, predU);
    const double rDeltaPhiUm =
        radius * fwdQaWrapPhi(localPhi - predPhi) * kQaCmToMicron;
    const double globalPhi = std::atan2(row.globalY, row.globalX);
    hRDeltaPhiVsGlobalPhi->Fill(globalPhi, rDeltaPhiUm);
    pRDeltaPhiVsGlobalPhi->Fill(globalPhi, rDeltaPhiUm);
    pRDeltaPhiDisk[row.disk]->Fill(globalPhi, rDeltaPhiUm);

    pResUVsSlopeU->Fill(row.slopeU, resUUm);
    pResVVsSlopeV->Fill(row.slopeV, resVUm);
    pResUVsSlopeV->Fill(row.slopeV, resUUm);
    pResVVsSlopeU->Fill(row.slopeU, resVUm);
  }

  std::cout << "Input rows: " << nEntries << std::endl;
  std::cout << "Tracks seen: " << tracksSeen << std::endl;
  std::cout << "Tracks selected: " << tracksSelected << std::endl;
  std::cout << "Rows selected: " << selectedRows << std::endl;
  std::cout << "Cuts: " << kQaEtaMin << " <= eta <= " << kQaEtaMax
            << ", p >= " << kQaPMin << ", pT >= " << kQaPtMin
            << ", nHitsFit >= " << kQaMinHitsFit
            << ", nFstHits >= " << kQaMinFstHits
            << ", trackType = " << kQaTrackType
            << ", fullyConverged = " << kQaRequireFullyConverged << std::endl;

  TCanvas *canvas = new TCanvas("cFwdAlignQa", "Planar FST alignment QA", 1400,
                                950);
  canvas->Print(pdfName + "[");

  canvas->Clear();
  canvas->Divide(3, 2);
  canvas->cd(1); hEtaAll->Draw();
  canvas->cd(2); hEtaSelected->Draw();
  canvas->cd(3); hP->Draw();
  canvas->cd(4); hPt->Draw();
  canvas->cd(5); hNHits->Draw();
  canvas->cd(6); hChi2Ndf->Draw();
  canvas->Print(pdfName);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1); hDiskWedge->Draw("colz text");
  canvas->cd(2); hGlobalWedge->Draw();
  canvas->cd(3); hUV->Draw("colz");
  canvas->cd(4); hRPhiVsR->Draw("colz");
  canvas->Print(pdfName);

  canvas->Clear();
  canvas->Divide(3, 2);
  canvas->cd(1); hResU->Draw(); fwdQaDrawZero(-residualRangeMicron, residualRangeMicron);
  canvas->cd(2); hResV->Draw(); fwdQaDrawZero(-residualRangeMicron, residualRangeMicron);
  canvas->cd(3); hResR->Draw(); fwdQaDrawZero(-residualRangeMicron, residualRangeMicron);
  canvas->cd(4); hResT->Draw(); fwdQaDrawZero(-residualRangeMicron, residualRangeMicron);
  canvas->cd(5); hPullU->Draw(); fwdQaDrawZero(-10, 10);
  canvas->cd(6); hPullV->Draw(); fwdQaDrawZero(-10, 10);
  canvas->Print(pdfName);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1); pResUByWedge->Draw(); fwdQaDrawZero(-0.5, 35.5); fwdQaDrawDiskBoundaries(pResUByWedge->GetMinimum(), pResUByWedge->GetMaximum());
  canvas->cd(2); pResVByWedge->Draw(); fwdQaDrawZero(-0.5, 35.5); fwdQaDrawDiskBoundaries(pResVByWedge->GetMinimum(), pResVByWedge->GetMaximum());
  canvas->cd(3); pPullUByWedge->Draw(); fwdQaDrawZero(-0.5, 35.5); fwdQaDrawDiskBoundaries(pPullUByWedge->GetMinimum(), pPullUByWedge->GetMaximum());
  canvas->cd(4); pPullVByWedge->Draw(); fwdQaDrawZero(-0.5, 35.5); fwdQaDrawDiskBoundaries(pPullVByWedge->GetMinimum(), pPullVByWedge->GetMaximum());
  canvas->Print(pdfName);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1); pResUVsU->Draw(); fwdQaDrawZero(kQaLocalUMin, kQaLocalUMax);
  canvas->cd(2); pResUVsV->Draw(); fwdQaDrawZero(kQaLocalVMin, kQaLocalVMax);
  canvas->cd(3); pResVVsU->Draw(); fwdQaDrawZero(kQaLocalUMin, kQaLocalUMax);
  canvas->cd(4); pResVVsV->Draw(); fwdQaDrawZero(kQaLocalVMin, kQaLocalVMax);
  canvas->Print(pdfName);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1); pResRVsR->Draw(); fwdQaDrawZero(kQaLocalUMin, kQaLocalUMax);
  canvas->cd(2); pResTVsR->Draw(); fwdQaDrawZero(kQaLocalUMin, kQaLocalUMax);
  canvas->cd(3); pResRVsPhi->Draw(); fwdQaDrawZero(-15, 15);
  canvas->cd(4); pResTVsPhi->Draw(); fwdQaDrawZero(-15, 15);
  canvas->Print(pdfName);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1);
  hRDeltaPhiVsGlobalPhi->Draw("colz");
  pRDeltaPhiVsGlobalPhi->SetMarkerStyle(20);
  pRDeltaPhiVsGlobalPhi->Draw("same");
  fwdQaDrawZero(-TMath::Pi(), TMath::Pi());
  fwdQaDrawPhiBoundaries(-residualRangeMicron, residualRangeMicron);
  for (int disk = 0; disk < kQaNumDisks; ++disk) {
    canvas->cd(disk + 2);
    pRDeltaPhiDisk[disk]->Draw();
    fwdQaDrawZero(-TMath::Pi(), TMath::Pi());
    fwdQaDrawPhiBoundaries(pRDeltaPhiDisk[disk]->GetMinimum(),
                           pRDeltaPhiDisk[disk]->GetMaximum());
  }
  canvas->Print(pdfName);

  canvas->Clear();
  canvas->Divide(2, 2);
  canvas->cd(1); pResUVsSlopeU->Draw(); fwdQaDrawZero(-0.25, 0.25);
  canvas->cd(2); pResVVsSlopeV->Draw(); fwdQaDrawZero(-0.25, 0.25);
  canvas->cd(3); pResUVsSlopeV->Draw(); fwdQaDrawZero(-0.25, 0.25);
  canvas->cd(4); pResVVsSlopeU->Draw(); fwdQaDrawZero(-0.25, 0.25);
  canvas->Print(pdfName);

  canvas->Print(pdfName + "]");
  outputFile->Write();
  outputFile->Close();
  inputFile->Close();

  std::cout << "Wrote " << rootName << " and " << pdfName << std::endl;
}
