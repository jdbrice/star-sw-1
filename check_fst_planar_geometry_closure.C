// Geometry closure for the real FST planar-measurement path.
//
// Generate a geometry cache for the data (for example y2022):
//   root4star -l -b -q 'StRoot/StFwdTrackMaker/macro/build_geom.C("y2022","fGeom.root")'
// Run from the repository root with StFwdTrackMaker built and its libraries on
// LD_LIBRARY_PATH, then load the STAR libraries and compile this test:
//   .L StRoot/StFwdTrackMaker/macro/mudst/fwd_afterburner.C
//   loadLibs();
//   .L check_fst_planar_geometry_closure.C+
//   gSystem->Exit(check_fst_planar_geometry_closure("fGeom.root"));
//
// The test covers every physical strip center on all 108 FST surfaces. It
// calls TrackFitter::createTrackPointFromPlanarMeasurement, compares the
// resulting wedge-local measurement with an independent AGML FTUS transform,
// checks both directions of the DetPlane local/global transform, and verifies
// that the former outer-gap signs fail by one degree.

#ifdef __CINT__

int check_fst_planar_geometry_closure(const char *geometryFile = "fGeom.root");

#else

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoTube.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TVector3.h"

#include "StEvent/StFstConsts.h"
#include "StFwdTrackMaker/FwdTrackerConfig.h"
#include "StFwdTrackMaker/include/Tracker/TrackFitter.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>

namespace {

double closureWrapPhi(double value) {
  return TMath::ATan2(TMath::Sin(value), TMath::Cos(value));
}

TVector2 formerOuterGapLocal(int globalSurface, double radius,
                             double stripPhi) {
  const int globalWedge = globalSurface / kFstNumSensorsPerWedge;
  const int disk = globalWedge / kFstNumWedgePerDisk;
  const int wedge = globalWedge % kFstNumWedgePerDisk;
  const int sensor = globalSurface % kFstNumSensorsPerWedge;
  const double sign = kFstzFilp[disk] * kFstzDirct[wedge];
  const double halfWedgePhi =
      0.5 * kFstNumPhiSegPerWedge * kFstStripPitchPhi;
  const double edgeToCenterPhi = halfWedgePhi - 0.5 * kFstStripPitchPhi;

  double dphi = sign * (stripPhi - edgeToCenterPhi);
  if (sensor != 0) {
    const double formerOffset =
        (sensor == 1) ? -0.5 * kFstStripGapPhi : 0.5 * kFstStripGapPhi;
    dphi = sign * (edgeToCenterPhi - stripPhi + formerOffset);
  }
  return TVector2(radius * TMath::Cos(dphi),
                  radius * TMath::Sin(dphi));
}

} // namespace

int check_fst_planar_geometry_closure(const char *geometryFile =
                                          "fGeom.root") {
  const int expectedSurfaces = 108;
  const long long expectedCenters = 36864;
  const double geometryToleranceCm = 1.0e-5;
  const double transformToleranceCm = 1.0e-10;

  if (!geometryFile || gSystem->AccessPathName(geometryFile)) {
    std::printf("GEOMETRY_IMPORT=FAIL file=%s\n",
                geometryFile ? geometryFile : "(null)");
    return 1;
  }
  TGeoManager::Import(geometryFile);
  TGeoManager *geometry = gGeoManager;
  if (!geometry) {
    std::printf("GEOMETRY_IMPORT=FAIL file=%s\n", geometryFile);
    return 1;
  }

  FwdTrackerConfig config;
  TrackFitter fitter(config, geometryFile);
  FwdGeomUtils geometryUtils(geometry);
  fitter.createAllFstPlanes(geometryUtils);

  static const int electronicToFstw[kFstNumWedgePerDisk] =
      {2, 7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8};
  static const int electronicToFstwDisk2[kFstNumWedgePerDisk] =
      {7, 1, 12, 6, 11, 5, 10, 4, 9, 3, 8, 2};
  static const int eventSensorToFtus[kFstNumSensorsPerWedge] = {3, 1, 2};

  int resolvedSurfaces = 0;
  long long checkedCenters = 0;
  long long measurementFailures = 0;
  long long geometryMappingFailures = 0;
  long long usedFwdHitGlobalFailures = 0;
  double maxLocalToAgmlCm = 0.0;
  double maxAgmlToLocalCm = 0.0;
  double maxDeltaPhiDeg = 0.0;
  double maxRoundTripCm = 0.0;
  double minFormerOuterDeltaPhiDeg =
      std::numeric_limits<double>::infinity();
  double maxFormerOuterDeltaPhiDeg = 0.0;

  // TrackFitter currently emits two INFO messages for every planar hit. Keep
  // the 36,864-point closure output readable without changing production code.
  const TString verboseLog = TString::Format(
      "/tmp/fst_planar_geometry_closure_%d.log", gSystem->GetPid());
  RedirectHandle_t redirectHandle;
  const bool redirected =
      gSystem->RedirectOutput(verboseLog.Data(), "w", &redirectHandle) == 0;

  for (int disk = 0; disk < kFstNumDisk; ++disk) {
    const int *wedgeMap =
        (disk == 1) ? electronicToFstwDisk2 : electronicToFstw;
    for (int wedge = 0; wedge < kFstNumWedgePerDisk; ++wedge) {
      const TString wedgePath = TString::Format(
          "/HALL_1/CAVE_1/FSTM_1/FSTD_%d/FSTW_%d", disk + 4,
          wedgeMap[wedge]);
      if (!geometry->cd(wedgePath)) {
        ++geometryMappingFailures;
        continue;
      }
      const TGeoMatrix *wedgeMatrix = geometry->GetCurrentMatrix();
      if (!wedgeMatrix) {
        ++geometryMappingFailures;
        continue;
      }
      const double expectedWedgePhi =
          0.5 * (kFstphiStart[wedge] + kFstphiStop[wedge]) *
          TMath::Pi() / 6.0;
      const Double_t *wedgeRotation = wedgeMatrix->GetRotationMatrix();
      const double mappedWedgePhi =
          TMath::ATan2(wedgeRotation[3], wedgeRotation[0]);
      if (TMath::Abs(closureWrapPhi(mappedWedgePhi - expectedWedgePhi)) >
          1.0e-10)
        ++geometryMappingFailures;

      for (int sensor = 0; sensor < kFstNumSensorsPerWedge; ++sensor) {
        const TString sensorPath = TString::Format(
            "%s/FTUS_%d", wedgePath.Data(), eventSensorToFtus[sensor]);
        if (!geometry->cd(sensorPath))
          continue;
        TGeoNode *sensorNode = geometry->GetCurrentNode();
        TGeoMatrix *currentMatrix = geometry->GetCurrentMatrix();
        if (!sensorNode || !sensorNode->GetVolume() || !currentMatrix)
          continue;
        TGeoTubeSeg *shape = dynamic_cast<TGeoTubeSeg *>(
            sensorNode->GetVolume()->GetShape());
        if (!shape)
          continue;
        const double expectedRmin = (sensor == 0) ? 5.0 : 16.5;
        const double expectedRmax = (sensor == 0) ? 16.5 : 28.0;
        if (TMath::Abs(shape->GetRmin() - expectedRmin) > 1.0e-10 ||
            TMath::Abs(shape->GetRmax() - expectedRmax) > 1.0e-10) {
          ++geometryMappingFailures;
          continue;
        }
        const TGeoHMatrix sensorMatrix(*currentMatrix);
        ++resolvedSurfaces;

        const int rBegin = (sensor == 0) ? 0 : 4;
        const int phiBegin = (sensor == 2) ? 64 : 0;
        const int phiCount = (sensor == 0) ? 128 : 64;
        const int globalSurface =
            (disk * kFstNumWedgePerDisk + wedge) *
                kFstNumSensorsPerWedge +
            sensor;

        for (int localRBin = 0; localRBin < 4; ++localRBin) {
          const int rStrip = rBegin + localRBin;
          const double radius =
              kFstrStart[rStrip] + 0.5 * kFstStripPitchR;
          const double agmlRadius =
              shape->GetRmin() +
              (localRBin + 0.5) * (shape->GetRmax() - shape->GetRmin()) /
                  4.0;

          for (int localPhiBin = 0; localPhiBin < phiCount;
               ++localPhiBin) {
            const int phiStrip = phiBegin + localPhiBin;
            const int agmlPhiBin =
                (sensor == 0) ? 127 - phiStrip : localPhiBin;
            const double agmlPhi =
                (shape->GetPhi1() +
                 (agmlPhiBin + 0.5) *
                     (shape->GetPhi2() - shape->GetPhi1()) / phiCount) *
                TMath::DegToRad();
            const Double_t local[3] = {
                agmlRadius * TMath::Cos(agmlPhi),
                agmlRadius * TMath::Sin(agmlPhi), 0.0};
            Double_t master[3] = {0.0, 0.0, 0.0};
            sensorMatrix.LocalToMaster(local, master);
            const TVector3 agmlGlobal(master[0], master[1], master[2]);

            TMatrixDSym covariance(3);
            covariance.Zero();
            covariance(0, 0) = 1.0e-4;
            covariance(1, 1) = 1.0e-4;
            covariance(2, 2) = 1.0e-4;
            // Deliberately wrong FwdHit XYZ: the FST planar coordinate path
            // must use the strip-native fields below, not this legacy global.
            const TVector3 sentinelGlobal =
                agmlGlobal + TVector3(3.0, 4.0, 5.0);
            FwdHit hit(checkedCenters, sentinelGlobal.X(), sentinelGlobal.Y(),
                       sentinelGlobal.Z(), disk + 4, kFstId, 0, covariance,
                       std::shared_ptr<McTrack>(), nullptr, globalSurface);
            hit._localPosition[0] = static_cast<float>(radius);
            hit._localPosition[1] = static_cast<float>(
                phiStrip * kFstStripPitchPhi);

            int hitId = 0;
            genfit::TrackPoint *point =
                fitter.createTrackPointFromPlanarMeasurement(
                    std::shared_ptr<genfit::Track>(), &hit, hitId);
            genfit::SharedPlanePtr plane = fitter.getPlaneFor(&hit);
            if (!point || !plane || point->getNumRawMeasurements() != 1) {
              ++measurementFailures;
              delete point;
              continue;
            }
            genfit::AbsMeasurement *raw = point->getRawMeasurement(0);
            if (!raw || raw->getRawHitCoords().GetNrows() != 2) {
              ++measurementFailures;
              delete point;
              continue;
            }

            const TVectorD &uv = raw->getRawHitCoords();
            const TVector2 measuredLocal(uv[0], uv[1]);
            const TVector3 measuredGlobal = plane->toLab(measuredLocal);
            const TVector2 agmlLocal = plane->LabToPlane(agmlGlobal);
            const TVector2 sentinelLocal = plane->LabToPlane(sentinelGlobal);
            const TVector2 roundTrip = plane->LabToPlane(measuredGlobal);

            maxLocalToAgmlCm = std::max(
                maxLocalToAgmlCm, (measuredGlobal - agmlGlobal).Mag());
            maxAgmlToLocalCm = std::max(
                maxAgmlToLocalCm, (agmlLocal - measuredLocal).Mod());
            maxDeltaPhiDeg = std::max(
                maxDeltaPhiDeg,
                TMath::Abs(closureWrapPhi(measuredGlobal.Phi() -
                                          agmlGlobal.Phi())) *
                    TMath::RadToDeg());
            maxRoundTripCm = std::max(
                maxRoundTripCm, (roundTrip - measuredLocal).Mod());
            if ((measuredLocal - sentinelLocal).Mod() < 1.0e-6)
              ++usedFwdHitGlobalFailures;

            if (sensor != 0) {
              const TVector2 formerLocal = formerOuterGapLocal(
                  globalSurface, radius,
                  phiStrip * kFstStripPitchPhi);
              const TVector3 formerGlobal = plane->toLab(formerLocal);
              const double formerDeltaPhiDeg = TMath::Abs(
                  closureWrapPhi(formerGlobal.Phi() - agmlGlobal.Phi())) *
                  TMath::RadToDeg();
              minFormerOuterDeltaPhiDeg = std::min(
                  minFormerOuterDeltaPhiDeg, formerDeltaPhiDeg);
              maxFormerOuterDeltaPhiDeg = std::max(
                  maxFormerOuterDeltaPhiDeg, formerDeltaPhiDeg);
            }

            delete point;
            ++checkedCenters;
          }
        }
      }
    }
  }

  if (redirected)
    gSystem->RedirectOutput(nullptr, nullptr, &redirectHandle);

  const bool coveragePass =
      resolvedSurfaces == expectedSurfaces &&
      checkedCenters == expectedCenters && measurementFailures == 0 &&
      geometryMappingFailures == 0 && usedFwdHitGlobalFailures == 0 &&
      !TrackFitter::kUseSpacePoints;
  const bool agmlPass =
      coveragePass && maxLocalToAgmlCm <= geometryToleranceCm &&
      maxAgmlToLocalCm <= geometryToleranceCm;
  const bool roundTripPass =
      coveragePass && maxRoundTripCm <= transformToleranceCm;
  const bool negativeControlPass =
      TMath::Finite(minFormerOuterDeltaPhiDeg) &&
      minFormerOuterDeltaPhiDeg > 0.99 &&
      maxFormerOuterDeltaPhiDeg < 1.01;
  const bool pass =
      coveragePass && agmlPass && roundTripPass && negativeControlPass;

  std::printf("FST_SURFACE_COVERAGE=%s resolved=%d expected=%d\n",
              resolvedSurfaces == expectedSurfaces ? "PASS" : "FAIL",
              resolvedSurfaces, expectedSurfaces);
  std::printf(
      "FST_STRIP_CENTER_COVERAGE=%s checked=%lld expected=%lld "
      "measurementFailures=%lld mappingFailures=%lld "
      "usedFwdHitGlobalFailures=%lld planarMode=%d\n",
      coveragePass ? "PASS" : "FAIL", checkedCenters, expectedCenters,
      measurementFailures, geometryMappingFailures,
      usedFwdHitGlobalFailures, TrackFitter::kUseSpacePoints ? 0 : 1);
  std::printf(
      "PLANAR_LOCAL_TO_AGML_GLOBAL=%s maxDistanceCm=%.12g "
      "maxDeltaPhiDeg=%.12g\n",
      agmlPass ? "PASS" : "FAIL", maxLocalToAgmlCm, maxDeltaPhiDeg);
  std::printf("AGML_GLOBAL_TO_PLANAR_LOCAL=%s maxDistanceCm=%.12g\n",
              agmlPass ? "PASS" : "FAIL", maxAgmlToLocalCm);
  std::printf("DETPLANE_ROUND_TRIP_SMOKE=%s maxDistanceCm=%.12g\n",
              roundTripPass ? "PASS" : "FAIL", maxRoundTripCm);
  std::printf(
      "FORMER_GAP_SIGN_NEGATIVE_CONTROL=%s minOuterDeltaPhiDeg=%.12g "
      "maxOuterDeltaPhiDeg=%.12g\n",
      negativeControlPass ? "PASS" : "FAIL",
      minFormerOuterDeltaPhiDeg, maxFormerOuterDeltaPhiDeg);
  std::printf("FST_PLANAR_GEOMETRY_CLOSURE=%s file=%s\n",
              pass ? "PASS" : "FAIL", geometryFile);

  if (redirected) {
    if (pass) {
      gSystem->Unlink(verboseLog.Data());
    } else {
      std::printf("TRACKFITTER_VERBOSE_LOG=%s\n", verboseLog.Data());
    }
  }
  return pass ? 0 : 1;
}

#endif
