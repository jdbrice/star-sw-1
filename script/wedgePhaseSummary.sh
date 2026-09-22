#!/bin/bash
# wedgePhaseSummary.sh -- merge the wedge-phase campaign and print the |delta| table.
# Usage: script/wedgePhaseSummary.sh [outdir]
#   g100 = mirror OFF (the question), g110 = mirror ON (positive control),
#   g010 = mirror ON without the gap fix (separates the two).
G=/gpfs01/star/pwg_tasks/FwdCalib/akio
OUT=${1:-/tmp/wedgephase}
mkdir -p "$OUT"
for c in g100 g110 g010; do
    d=$G/pico_wp_${c}_P_20260922/blinddiag
    n=$(ls "$d"/*.fwd_blind_diag.root 2>/dev/null | wc -l)
    echo "=== $c : $n / 5 diagnostic files"
    [ "$n" -eq 0 ] && continue
    hadd -f "$OUT/wp_${c}.root" "$d"/*.fwd_blind_diag.root > /dev/null 2>&1
    root4star -l -b -q "script/fitBlindWedgePhase.C(\"$OUT/wp_${c}.root\",\"$c\")" 2>&1 \
        | tr -d '\000' | grep ">>>"
done
echo "=== overlay plot"
[ -f "$OUT/wp_g100.root" ] && [ -f "$OUT/wp_g110.root" ] && \
  root4star -l -b -q "script/plotBlindWedgePhase.C(\"$OUT/wp_g100.root\",\"g100\",\"$OUT/wp_g110.root\",\"g110\",\"$OUT/wedge_phase.png\")" 2>&1 \
    | tr -d '\000' | grep -i "png has been created"
