#!/bin/sh
# Waits for the two fwd_afterburner_db gap-fix test runs to finish, copies their
# outputs out of the session scratchpad (which lives under /tmp and may not
# survive), then runs both check macros and leaves the results in gapfix_test/.
#
# Launched detached with nohup so it survives the interactive session ending.
# Progress/results: gapfix_test/STATUS
#
# See xihe_crosscheck_20260729.txt and fst_wedge_testA_results.txt for context.

SCR=/tmp/claude-2546/-direct-star-u-akio-fcstrk11-star-sw-fwd/999e4fd5-b5d9-4cdc-bbed-3d4d48f7c08f/scratchpad
REPO=/star/u/akio/fcstrk11/star-sw-fwd
OUT=$REPO/gapfix_test
BASE=st_fwd_23081008_raw_1000025.FwdAlignment_BLCVtx.root

mkdir -p $OUT
echo "waiting for afterburner runs to finish (started $(date))" > $OUT/STATUS

# The runs are `root4star ... macro/mudst/fwd_afterburner_db.C(...)`. This
# script's own cmdline is just "sh gapfix_finish.sh", so it cannot match itself.
while pgrep -f 'macro/mudst/fwd_afterburner_db' > /dev/null 2>&1 ; do
    sleep 60
done
echo "runs finished $(date); copying outputs" >> $OUT/STATUS

for d in gapoff gapon ; do
    if [ -f $SCR/$d/$BASE ] ; then
        cp -f $SCR/$d/$BASE $OUT/$d.FwdAlignment_BLCVtx.root
        echo "  copied $d" >> $OUT/STATUS
    else
        echo "  MISSING $d/$BASE" >> $OUT/STATUS
    fi
    [ -f $SCR/$d/run.log ] && tail -3000 $SCR/$d/run.log > $OUT/$d.run.tail.log
done

cd $REPO || exit 1

# 1. hit POSITIONS -- the decisive check, independent of any fit
echo "running checkFstGapFixEffect.C $(date)" >> $OUT/STATUS
root4star -l -b -q "script/checkFstGapFixEffect.C(\"$OUT/gapoff.FwdAlignment_BLCVtx.root\",\"$OUT/gapon.FwdAlignment_BLCVtx.root\",\"$OUT\")" \
    > $OUT/1_positions.txt 2>&1

# 2. centreline step in the residual, off vs on
echo "running checkFstOuterHalf.C x2 $(date)" >> $OUT/STATUS
root4star -l -b -q "script/checkFstOuterHalf.C(\"$OUT/gapoff.FwdAlignment_BLCVtx.root\",\"$OUT/outerhalf_off\")" \
    > $OUT/2_outerhalf_off.txt 2>&1
root4star -l -b -q "script/checkFstOuterHalf.C(\"$OUT/gapon.FwdAlignment_BLCVtx.root\",\"$OUT/outerhalf_on\")" \
    > $OUT/2_outerhalf_on.txt 2>&1

# short digest so the answer is readable without opening the big logs
{
    echo "=== 1. HIT POSITIONS (decisive; no fit involved) ==="
    grep -A 12 'Outer-sensor hits: occupancy' $OUT/1_positions.txt 2>/dev/null
    grep -A 5  'Inner-sensor control'         $OUT/1_positions.txt 2>/dev/null
    echo
    echo "=== 2. CENTRELINE STEP IN mean r*dphi -- gap fix OFF ==="
    grep -A 16 'disk parity  region' $OUT/2_outerhalf_off.txt 2>/dev/null
    echo
    echo "=== 2. CENTRELINE STEP IN mean r*dphi -- gap fix ON ==="
    grep -A 16 'disk parity  region' $OUT/2_outerhalf_on.txt 2>/dev/null
} > $OUT/SUMMARY.txt 2>&1

echo "ALL DONE $(date)" >> $OUT/STATUS
