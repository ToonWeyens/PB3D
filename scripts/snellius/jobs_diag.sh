#!/bin/bash
# jobs_diag.sh -- W7-X check_vac with vac_diag_style 1 (indirect: row sum and identities, the
# default) and 2 (direct: corner term plus cell integral) on the grids of the alignment study,
# plus ctest Release as the build gate (issue #4, 6 November 2026, commit 225aa72).
# The small grids run with 1 and 4 processes to check the distributed direct diagonal.
#   nohup ~/pb3d/scripts/jobs_diag.sh > ~/pb3d/logs/jobs_diag.log 2>&1 &
set -uo pipefail
source "$HOME/pb3d/env.sh"
LOGS=$PB3D_ROOT/logs; D=$PB3D_ROOT/decks; FX=$PB3D_FIXTURE_DIR; R=$PB3D_ROOT/scripts/run_deck.sbatch
stamp() { echo "[$(date "+%F %T")] $*"; }
stamp "exe: $(ls -la $PB3D_ROOT/build/PB3D | cut -c25-) commit $(cat $PB3D_ROOT/src/.commit)"
SB="sbatch --parsable --partition=rome --nodes=1 --output=$LOGS/%x-%j.out"
$SB --ntasks=16 --time=00:40:00 --job-name=d8_ctest --wrap="source $PB3D_ROOT/env.sh; cd $PB3D_ROOT/build && ctest --output-on-failure --timeout 1400"
# run_deck.sbatch <deck> [np-list] [bc-style] [repeat] [exe] [keep] [eq-file]
for s in 1 2; do
    for g in na16_p61 na32_p61 na64_p61 na32_p121; do
        $SB --ntasks=4 --mem=60G --time=00:30:00 --job-name=d8s${s}_w7x_$g "$R" "$D/input_cv_w7x_${g}_d$s" "1 4" "" 1 "" "" "$FX/wout_w7xflip.nc"
    done
    for g in na128_p61 na64_p121 na128_p121; do
        $SB --ntasks=4 --mem=100G --time=01:30:00 --job-name=d8s${s}_w7x_$g "$R" "$D/input_cv_w7x_${g}_d$s" "4" "" 1 "" "" "$FX/wout_w7xflip.nc"
    done
    for g in na256_p121 na128_p241; do
        $SB --ntasks=4 --mem=200G --time=02:30:00 --job-name=d8s${s}_w7x_$g "$R" "$D/input_cv_w7x_${g}_d$s" "4" "" 1 "" "" "$FX/wout_w7xflip.nc"
    done
done
stamp "submitted; queued or running: $(squeue -u "$USER" -h | wc -l)"
