#!/bin/bash
# One command to turn the finished DB-geometry campaign into numbers and plots.
#   merge -> per-station fit vs the hardcoded run -> per-quadrant table -> plots
# Usage: makeDbGeomResults.sh [outdir]
set -e
cd /direct/star+u/akio/fcstrk11/star-sw-fwd
O=${1:-/gpfs01/star/pwg_tasks/FwdCalib/akio/pico_dbgeom_20260929}
S=/tmp/claude-2546/-direct-star-u-akio-fcstrk11-star-sw-fwd/999e4fd5-b5d9-4cdc-bbed-3d4d48f7c08f/scratchpad

N=$(ls $O/blinddiag/*.root 2>/dev/null | wc -l)
echo "=== $N blinddiag files ==="
# Guard against the failure that invalidated the 20260928 run: confirm every job
# actually used the DB tables rather than silently falling back.
echo "jobs reporting DB tables : $(grep -ah 'quadrant offsets built from' $O/log/*.log 2>/dev/null | wc -l)"
echo "jobs reporting fallback  : $(grep -ah 'falling back to the hardcoded' $O/log/*.log 2>/dev/null | wc -l)"

hadd -f blinddiag_dbgeom_final.root $O/blinddiag/*.root > /dev/null 2>&1
echo "merged -> blinddiag_dbgeom_final.root"

# Each macro gets its OWN root process. Loading them into one CINT session collides:
# cmpBlindShift.C declares gMu as a scalar and stgcQuadOffsets.C declares gMu[4][4][2],
# and the result is silent garbage -- every cell reported the same number.
echo "--- per-station, hardcoded vs DB tables ---"
cat > $S/f1.C <<'M1'
void f1(){ gROOT->LoadMacro("script/cmpBlindShift.C");
  cmpBlindShift("blinddiag_stgcmis_20260925.root","blinddiag_dbgeom_final.root"); }
M1
root -l -b -q $S/f1.C 2>&1 | sed -n '/--- dx/,$p' | head -24

echo "--- per-quadrant residuals with the tables in ---"
cat > $S/f2.C <<'M2'
void f2(){ gROOT->LoadMacro("script/stgcQuadOffsets.C");
  stgcQuadOffsets("blinddiag_dbgeom_final.root"); }
M2
root -l -b -q $S/f2.C 2>&1 | grep "^>>>" | head -24

cat > $S/f3.C <<'M3'
void f3(){ gROOT->LoadMacro("script/plotStgcRadial.C");
  plotStgcRadial("blinddiag_dbgeom_final.root","FstFttFlipTest/plots/dbgeom_radial.png"); }
M3
root -l -b -q $S/f3.C > /dev/null 2>&1

cat > $S/f4.C <<'M4'
void f4(){ gROOT->LoadMacro("script/plotBlindQuadFit.C");
  plotBlindQuadFit("blinddiag_dbgeom_final.root","FstFttFlipTest/plots/dbgeom_quad_all.png",-1);
  for(int d=0;d<4;d++)
    plotBlindQuadFit("blinddiag_dbgeom_final.root",
                     Form("FstFttFlipTest/plots/dbgeom_quad_disk%d.png",d),d); }
M4
root -l -b -q $S/f4.C > /dev/null 2>&1

echo "=== plots ==="
ls -la FstFttFlipTest/plots/dbgeom_*.png 2>/dev/null
