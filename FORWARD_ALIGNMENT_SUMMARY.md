# Forward STAR Alignment: Current Audit and Plan

Last audited: 2026-07-13

Branch: `xihe-align`
Target detector: FST first; FTT is intentionally deferred.

## 1. Executive Summary

The project accumulated several generations of geometry experiments, tree branches,
and QA plots. The useful core is smaller than the repository currently suggests:

1. Read strip-native FST hit information from `StMuFstHit`.
2. Convert `(radial strip, phi strip)` to a two-dimensional measurement in a
   wedge-oriented `(U,V)` frame.
3. Fit those measurements with GENFIT `PlanarMeasurement` objects.
4. Ask GENFIT for the unbiased residual of each FST measurement.
5. Write one compact row per valid two-dimensional FST hit.
6. Diagnose residuals by disk and wedge before solving any constants.

The alignment object is now the **wedge**, not an individual silicon sensor.
There are 3 disks and 12 wedges per disk, for 36 potential alignment objects.
The three silicon pieces within a wedge are retained only as three physical-z
measurement surfaces. They do not receive independent alignment parameters.

The compact tree and rewritten QA macro are deliberately incompatible with old
`align_test.root` files. Rebuild `StFwdTrackMaker` and rerun the afterburner before
using the new QA.

## 2. Recommended Architecture

The primary alignment path should remain in planar mode:

```text
StMuFstHit strip data
  -> FwdHit strip-native local fields
  -> wedge-local PlanarMeasurement(U,V,C)
  -> GENFIT Kalman fit
  -> unbiased planar residual and covariance
  -> compact fwdAlign tree
  -> wedge QA
  -> constrained wedge alignment solve
  -> update geometry/constants
  -> rerun tracking for closure
```

This is preferable to treating global hit positions as the alignment
measurement. Alignment parameters move detector planes. The measured coordinates
should remain in the plane's local basis while the plane origin and orientation
come from geometry.

## 3. Current Coordinate Model

### 3.1 Raw FST information

`StFwdHitLoader.cxx` reads:

```text
mean radial strip -> radial center r [cm]
mean phi strip    -> stripPhi = meanPhiStrip * phi pitch [rad]
disk, wedge, sensor/surface identifiers
```

These strip-native values are stored in `FwdHit::_localPosition`:

```text
_localPosition[0] = r
_localPosition[1] = stripPhi
```

They are not global Cartesian coordinates.

### 3.2 Wedge-local measurement

`TrackFitter::createTrackPointFromPlanarMeasurement` converts the strip-native
values to wedge-oriented Cartesian coordinates:

```text
measU = r cos(dphi)
measV = r sin(dphi)
```

Here `dphi` is the hit's angular displacement from the wedge centerline. Its
sign depends on the disk and electronic wedge orientation. The outer regions
also require the half-wedge and central-gap convention encoded by the FST strip
numbering.

This conversion is detector readout decoding. It does not apply an alignment
correction. A future alignment changes the plane transform, not the raw strip
indices.

### 3.3 Plane transform

`FwdGeomUtils::getFstWedgeOrigin` obtains the nominal wedge origin and basis from
`fGeom.root`. `TrackFitter::createAllFstPlanes` then creates:

- 36 wedge reference planes, one per `(disk,wedge)`;
- 108 measurement planes, three physical-z surfaces per wedge.

All three surfaces in one wedge share the same wedge `(U,V)` orientation and
nominal wedge `x,y`. Their `z` values come from the corresponding active FTUS
geometry surface.

The measured global point is therefore:

```text
X_global = O_surface + measU * U_wedge + measV * V_wedge
```

GENFIT uses this plane and the two local coordinates during fitting.

### 3.4 Identifier convention

The planar FST plane ID is:

```text
planeId = disk * 36 + wedge * 3 + surface
```

with:

```text
disk    = 0, 1, 2
wedge   = 0 ... 11 within a disk
surface = 0, 1, 2 within a wedge
planeId = 0 ... 107
```

GENFIT sorting uses `planeId + 1` because sorting value zero is reserved for the
primary vertex.

## 4. How The Residual Is Calculated

For a planar measurement vector `m=(U,V)` and an unbiased predicted track state
`x_pred`, GENFIT forms:

```text
r = m - H x_pred
```

The stored convention is therefore:

```text
resU = measured U - unbiased predicted U
resV = measured V - unbiased predicted V
```

`StFwdTrackMaker::FillAlignment` calls:

```cpp
kfi->getResidual(iMeas, false, false)
```

The first `false` requests an unbiased residual. GENFIT combines the forward and
backward Kalman predictions, which do not include the current measurement. This
is a fast leave-one-measurement-out state without rerunning the full fit.

The second `false` requests the full residual covariance:

```text
C_residual = C_measurement + H C_prediction H^T
```

The previous code requested measurement errors only. That made the residual
values usable but made the reported pulls incomplete because track prediction
uncertainty was missing from the denominator.

The compact tree stores `covUU`, `covUV`, and `covVV`. The QA currently shows
marginal pulls:

```text
pullU = resU / sqrt(covUU)
pullV = resV / sqrt(covVV)
```

Because `covUV` can be nonzero, a later solver should use the full 2x2 covariance
or whiten the residual vector. The diagonal pulls are diagnostics, not the final
least-squares weighting implementation.

## 5. Compact Alignment Tree

The tree remains named `fwdAlign`. It contains one row for each valid 2D FST
planar measurement. It no longer writes FTT rows, primary-vertex rows, biased
residual duplicates, three-component placeholders, or `-99999` sentinels.

### 5.1 Event and track fields

| Branch | Meaning |
|---|---|
| `run`, `event` | STAR event identifiers |
| `trackIndex` | Track-result index within the event |
| `trackType` | Forward tracker fit category |
| `nHitsFit` | Number of points used in the fit |
| `nFstHits` | Number of FST measurements on the fitted track |
| `fullyConverged` | GENFIT full-convergence flag |
| `chi2Ndf` | Track chi-square divided by NDF |
| `trackP`, `trackPt`, `trackEta` | Fitted track kinematics |

### 5.2 Geometry fields

| Branch | Meaning |
|---|---|
| `planeId` | FST measurement plane, 0 through 107 |
| `disk` | Disk index, 0 through 2 |
| `wedge` | Electronic wedge index within the disk, 0 through 11 |
| `surface` | Physical-z surface within the wedge, 0 through 2 |
| `globalX/Y/Z` | Measured point reconstructed from plane plus `(measU,measV)` |

### 5.3 Measurement and fit fields

| Branch | Meaning |
|---|---|
| `measU`, `measV` | Wedge-local measured coordinates in cm |
| `resU`, `resV` | Unbiased local residuals in cm |
| `covUU`, `covUV`, `covVV` | Full unbiased residual covariance in cm2 |
| `slopeU`, `slopeV` | Local track slopes `(p.U)/(p.N)` and `(p.V)/(p.N)` |

The local slopes are necessary for out-of-plane translations and tilts. They are
computed using the same plane basis as the measurement and residual.

## 6. QA Macro

`fwd_alignment_residual_qa.C` is now a normal ROOT C++ event loop rather than a
large collection of `TTree::Draw` expressions.

All editable selection values are grouped at the top:

```text
2.5 <= eta <= 4.0
p >= 1.0 GeV/c
pT >= 0.2 GeV/c
nHitsFit >= 8
nFstHits >= 3
trackType == 1
fullyConverged == true
```

Every alignment residual plot uses the same selection function. This removes the
old risk that disk, wedge, sensor, pull, and trend pages silently used different
cuts.

The PDF has eight focused pages:

1. Track selection distributions.
2. Disk/wedge occupancy and local `(U,V)` / `(r,rphi)` coverage.
3. Local U/V and radial/tangential residual and pull distributions.
4. Mean residual and pull by the 36 wedge identifiers.
5. U/V residual trends versus local U and V.
6. Radial/tangential trends versus local radius and local phi.
7. Exact `r*Delta(phi)` versus global phi, for all disks and each disk.
8. Residual correlations with local track slopes.

The exact angular diagnostic uses measured and predicted local points:

```text
predU = measU - resU
predV = measV - resV
DeltaPhi = wrap(atan2(measV,measU) - atan2(predV,predU))
rDeltaPhi = sqrt(measU^2 + measV^2) * DeltaPhi
```

Global phi for the horizontal axis is read from:

```text
atan2(globalY, globalX)
```

This plot can reveal a 12-fold wedge-step pattern. It should not be confused with
the direct wedge-V residual, although they are approximately equal for small
angular differences.

## 7. Comparison With StFwdAlignmentMaker

The colleague implementation in `jdbrice/star-sw-1` takes a different and useful
approach:

| Topic | This branch | `StFwdAlignmentMaker` |
|---|---|---|
| Measurement model | 2D planar wedge `(U,V)` | Existing production space points |
| Unbiased method | GENFIT forward/backward Kalman prediction | Remove one hit and refit the track |
| Cost | One nominal track fit | Approximately one extra fit per removed hit |
| Output | Local residual, covariance, local slopes | Global hit and projected XYZ |
| Projection surface | Actual FST measurement plane | Simple horizontal plane at hit z |
| Main use | Mechanical planar alignment solve | Independent global residual/12-fold QA |
| Geometry sensitivity | Explicit wedge origin, basis, and surface z | Uses legacy hit xyz and z convention |

The explicit remove-and-refit method is a strong nonlinear cross-check. It is not
the only valid definition of an unbiased residual: GENFIT's forward/backward
prediction also excludes the current hit. The two methods should agree within
expected linearization and fitter differences on a controlled sample.

The colleague output is especially useful for reproducing the STAR global
`r*Delta(phi)` step plot and testing empirical wedge corrections. It is not yet a
replacement for a local planar alignment tree because it does not provide the
native local residual covariance or local slopes, and its removed hit is projected
to a horizontal plane rather than the wedge measurement plane.

There is also a configuration difference to control during comparisons:
the colleague afterburner commonly uses BLCVtx track type 4, while the current QA
default is BLC track type 1. Comparisons are not meaningful until track type,
cuts, geometry, magnetic field, and input sample are matched.

Recommended synthesis:

1. Use this branch's planar tree as the alignment input.
2. Run `StFwdAlignmentMaker` on a smaller matched sample as an independent check.
3. Compare global `r*Delta(phi)` after transforming the planar result to global
   coordinates.
4. If needed, adapt the separate maker to project onto the actual wedge plane,
   but do not replace the local alignment model with a hard-coded phi table.

## 8. File Ownership And Status

### 8.1 Core files to keep

| File | Responsibility |
|---|---|
| `StFwdHitLoader.cxx` | Load FST hit identity, covariance, and strip-native values |
| `include/Tracker/FwdGeomUtils.h` | Map detector identifiers to nominal geometry planes |
| `include/Tracker/TrackFitter.h` | Build planar measurements and run GENFIT |
| `StFwdTrackMaker.cxx/.h` | Write the compact alignment tree |
| `macro/mudst/fwd_afterburner.C` | Run tracking with alignment output enabled |
| `macro/build_geom.C` | Produce the geometry cache |
| `fwd_alignment_residual_qa.C` | Standard post-fit alignment QA |

### 8.2 Validation-only files

| File | Status |
|---|---|
| `inspect_fst_ftus_geometry.C` | Keep as a geometry debugging utility |
| `draw_fst_fgeom_active.C` | Useful scratch visualization; currently untracked |
| `debug_fst_raw_phi.C` | Useful strip-decoding scratch check; currently untracked |

### 8.3 Experimental or stale files

`toy_fst_inner_disk_xygamma.C` was useful while exploring sensor-level constants,
but its name, assumptions, and solve are no longer the target wedge model. Do not
use its output as an alignment result. Either rewrite it as a constrained
wedge-level solver after the QA is stable or archive it outside the production
workflow.

`toy_fst_wedge_uvgamma.C` is the current compact-tree toy. It independently fits
`DeltaU`, `DeltaV`, and `gamma` for each of the 36 wedges using both residual
components and their full 2x2 covariance. It combines all three physical-z
surfaces in a wedge and reports normal-matrix condition numbers. Its after plots
subtract the fitted linear model from the existing residuals; they are algebraic
closure only. The compact tree does not contain covariance between different
hits on the same track, so the toy treats rows as independent and its formal
parameter errors can be optimistic.

Generated `.root`, `.pdf`, `.png`, `.so`, `.d`, `.log`, and `paw.metafile` files
are run products, not source. They should remain untracked. No user-generated
outputs were deleted during this audit.

## 9. Important Remaining Risks

### 9.1 FST covariance model

`makeFstCovMat` currently defaults to:

```text
rSize   = 3.0 cm
phiSize = 0.0040906154 rad
```

and squares these values directly. A `sqrt(12)` variable is declared but not used.
If these inputs are full strip/bin widths rather than one-sigma resolutions, the
measurement weights and pulls are wrong. The covariance is also first constructed
using legacy global hit phi and then rotated to the plane basis. This deserves a
dedicated validation against the strip cluster definition and the corrected
wedge-local coordinate convention.

This issue can affect the fit itself, so it must be studied separately from tree
cleanup.

### 9.2 Legacy FwdHit z versus planar surface z

The legacy `FwdHit` position can carry a flat or hard-coded FST z convention.
The planar measurement now uses the active FTUS surface z from geometry. GENFIT
places the planar hit at the latter value, but seed finding and covariance setup
can still inspect the legacy FwdHit coordinates. Therefore the mismatch is not
just cosmetic and requires a controlled tracking comparison.

### 9.3 Strip-to-wedge decoding

The inner and outer strip formulas now pass the existing Cartesian closure tests,
but the logic remains encoded through constants and sign tables. Add deterministic
tests that cover every disk, wedge orientation, outer half, radial strip, and phi
edge. A passing occupancy plot alone is not enough.

### 9.4 Geometry mapping duplication

Electronic-to-GEANT wedge maps appear in more than one `FwdGeomUtils` method.
Centralize them after the current behavior is validated, so a future geometry
change cannot update one path and leave another stale.

### 9.5 Verbose debug output

`TrackFitter::kVerbose` is still set to 1 and prints per-hit plane checks. One
debug `deltaPhi` expression subtracts a value from itself and is always zero.
These should be removed or made configuration-driven after the compact workflow
has been rerun once.

### 9.6 No final alignment solve yet

The per-wedge U/V/gamma toy is now available, but it does not impose disk/global
gauge constraints, write geometry updates, or perform a tracking-refit closure.
Its signs still require an injected geometry displacement followed by a complete
refit. It is a solver-development tool, not yet a trusted constants producer.

## 10. Proposed Wedge Alignment Solve

Start with the best-constrained in-plane degrees of freedom per wedge:

```text
DeltaU, DeltaV, gamma
```

where `gamma` is a small rotation around the plane normal. For a local hit
`q=(U,V)`, the first-order residual response has the form:

```text
Delta(resU) approximately -DeltaU + gamma * V
Delta(resV) approximately -DeltaV - gamma * U
```

The exact sign should be verified with injected-geometry tests, not assumed from
notation alone.

Only after this closes should the solve add:

```text
DeltaW, alpha, beta
```

The out-of-plane parameters require `slopeU` and `slopeV`, because moving a plane
along its normal changes the track intersection according to the local track
slope. This is why a delta-z solve without slope information was underdetermined.

Use the full 2x2 residual covariance for every row. Solve all wedges in a disk
simultaneously and impose gauge constraints, for example zero mean translation
and rotation for each disk or one fixed reference wedge. Otherwise common disk
motion and track-parameter changes create weak modes.

## 11. Validation Sequence

Proceed in this order:

1. Rebuild `StFwdTrackMaker` and rerun a small afterburner sample.
2. Confirm the compact tree has only finite FST rows and the expected 27 branches.
3. Check occupancy for all 3 disks, 12 wedges, and 3 surfaces.
4. Check local-to-global closure for inner and outer regions at micron scale.
5. Validate the covariance model and pull widths before fitting constants.
6. Match cuts and compare fast GENFIT-unbiased residuals with explicit
   remove-and-refit residuals on a small sample.
7. Inject one known `DeltaU`, `DeltaV`, or `gamma` into geometry and verify the
   solver recovers its magnitude and sign.
8. Solve real data with gauge constraints.
9. Write updated geometry/alignment constants.
10. Rerun the entire tracking fit and require residual closure on independent
    events. Algebraically subtracting a correction from an existing tree is not
    a closure test.

## 12. Current Build And Run Commands

Build only the forward tracking package in the running STAR container:

```bash
docker exec fwd bash -lc 'cd /work && cons +StRoot/StFwdTrackMaker'
```

Run the afterburner using the repository macro and a MuDst input:

```bash
root4star -l -b -q 'StRoot/StFwdTrackMaker/macro/mudst/fwd_afterburner.C("input.MuDst.root",100)'
```

Run the compact QA after producing a new `align_test.root`:

```bash
root4star -l -b -q 'fwd_alignment_residual_qa.C+("align_test.root","fwd_align_qa")'
```

Outputs are:

```text
fwd_align_qa.root
fwd_align_qa.pdf
```

Run the toy per-wedge in-plane solve:

```bash
root4star -l -b -q 'toy_fst_wedge_uvgamma.C+("align_test.root","toy_fst_wedge_uvgamma")'
```

## 13. External References

- STAR residual QA: <https://www.star.bnl.gov/protected/spin/akio/fcs/residual/index.html>
- Colleague implementation: <https://github.com/jdbrice/star-sw-1/tree/dev/StRoot/StFwdAlignmentMaker>
- GENFIT residual implementation: <https://github.com/GenFit/GenFit/blob/master/fitters/src/KalmanFitterInfo.cc>

The STAR residual page requires protected-site authentication. The comparison in
this audit uses the supplied screenshot, the public colleague source, and the
local STAR/GENFIT code paths.

## 14. Immediate Next Step

Produce a fresh compact tree, inspect the eight-page QA, and run the per-wedge
toy. Treat wedges with high condition numbers or too few rows as unsolved. Before
applying any reported constants, the next physics change should be the isolated
FST covariance study, followed by a known geometry-injection and full-refit sign
test.
