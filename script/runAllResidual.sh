#!/bin/bash
# Regenerate all StFwdResidualMaker outputs: 2 datasets (real data, pythia) x
# 3 track types (Global, BLC, BLCVtx). Run from the repo top level.
set -e -o pipefail
cd "$(dirname "$0")/.."

DATA_FILE="/star/data20/reco/production_pp500_2022/ReversedFullField/DEV_fwd_test2/2022/081/23081015/st_fwd_23081015_raw_6000057.MuDst.root"
DATA_NEVENTS=300
PYTHIA_NEVENTS=500
LOGDIR="script/logs"
mkdir -p "$LOGDIR"

echo "############################# REAL DATA #############################"
for T in 0 1 4; do
  echo "=== data residualTrackType=$T ==="
  # StFwdTrackMaker's LOG_INFO is extremely verbose (GB-scale over hundreds of
  # events) and isn't gated by debug=0 -- keep only the tail so a failure is
  # still diagnosable without risking the disk quota again.
  root4star -b -q "fwd_afterburner_db.C(\"$DATA_FILE\", $DATA_NEVENTS, 0, $T)" 2>&1 \
    | tail -n 1000 > "$LOGDIR/data_run_t${T}.log"
  echo "=== data T=$T done, exit=${PIPESTATUS[0]} ==="
done

echo "############################# PYTHIA J/PSI ###########################"
for T in 0 1 4; do
  echo "=== pythia residualTrackType=$T ==="
  sed -i "s/^int residualTrackType = .*/int residualTrackType = $T;  \/\/ StFwdTrack::StFwdTrackType: 0=Global 1=BLC 2=Primary 3=FwdVtx 4=BLCVtx 5=FCSConstrained/" script/sim.C
  root4star -b -q "sim.C($PYTHIA_NEVENTS,1,\"JPsi\",0)" 2>&1 \
    | tail -n 1000 > "$LOGDIR/pythia_run_t${T}.log"
  echo "=== pythia T=$T done, exit=${PIPESTATUS[0]} ==="
done
# restore default
sed -i "s/^int residualTrackType = .*/int residualTrackType = 0;  \/\/ StFwdTrack::StFwdTrackType: 0=Global 1=BLC 2=Primary 3=FwdVtx 4=BLCVtx 5=FCSConstrained/" script/sim.C

echo "############################# ALL DONE ###############################"
ls -la st_fwd_23081015_raw_6000057.FwdDetResidual_*.root pythia.JPsi.vz0.run1.FwdDetResidual_*.root
