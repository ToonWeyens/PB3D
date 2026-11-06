#!/bin/bash
# pipeline.sh -- the whole Snellius bring-up, unattended and in order:
# dependencies (PETSc with MUMPS, SLEPc, PSPLINE, minimal LIBSTELL), the
# Release and Debug builds, then the regression decks and ctest as Slurm
# jobs, strictly one at a time (sbatch --wait), each one node, partition
# rome, 15 minutes, output under ~/pb3d/logs/%x-%j.out. Stages that are
# already installed are skipped, so re-running is cheap.
#
#   nohup ~/pb3d/scripts/pipeline.sh > ~/pb3d/logs/pipeline.log 2>&1 &
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/../env.sh"
LOGS=$PB3D_ROOT/logs
DECKS=$PB3D_ROOT/src/tests/regression
SB="sbatch --wait --parsable --partition=rome --nodes=1 --ntasks=16 --time=00:15:00 --output=$LOGS/%x-%j.out"
stamp() { echo "[$(date '+%F %T')] $*"; }
step() {   # step <name> <cmd...> ; stops the pipeline on failure
    stamp "START $1"; local name=$1; shift
    if "$@" > "$LOGS/$name.log" 2>&1; then stamp "OK    $name"; else stamp "FAIL  $name (see $LOGS/$name.log)"; exit 1; fi
}
job() {    # job <name> <sbatch args...> ; one job at a time, waits for it
    stamp "JOB   $1"; local name=$1; shift
    local id; id=$($SB --job-name="$name" "$@" | cut -d';' -f1)
    [ -n "$id" ] || { stamp "FAIL  sbatch $name gave no job id"; exit 1; }
    stamp "      $name = job $id finished: $(sacct -j "$id" -n -o State,Elapsed -X | head -1 | tr -s ' ')"
    grep -E "^EV |tests passed|tests failed|Total Test time" "$LOGS/$name-$id.out" | head -40
}

t0=$(date +%s)
[ -f "$PETSC_DIR/lib/pkgconfig/PETSc.pc" ]  || step petsc    "$SCRIPT_DIR/build_deps.sh" petsc
[ -f "$SLEPC_DIR/lib/pkgconfig/SLEPc.pc" ]  || step slepc    "$SCRIPT_DIR/build_deps.sh" slepc
[ -f "$PSPLINE_DIR/lib/libpspline.a" ]      || step pspline  "$SCRIPT_DIR/build_deps.sh" pspline
[ -f "$LIBSTELL_DIR/lib/libstell.a" ]       || step libstell "$SCRIPT_DIR/build_deps.sh" libstell
step pb3d_rel   "$SCRIPT_DIR/build_pb3d.sh" release
step pb3d_dbg   "$SCRIPT_DIR/build_pb3d.sh" debug
stamp "build time: $(( $(date +%s) - t0 )) s"

# run_deck.sbatch <deck> [np-list] [bc-style] [repeat] [exe]; options are
# positional because Snellius drops the submitting shell's environment
job pb3d-fixed      "$SCRIPT_DIR/run_deck.sbatch" "$DECKS/input_cbm18a_fixed_bnd" "0 1 2 4"
job pb3d-free-bc4   "$SCRIPT_DIR/run_deck.sbatch" "$DECKS/input_cbm18a_free_bnd"  "0 1 2 4"
job pb3d-free-bc2   "$SCRIPT_DIR/run_deck.sbatch" "$DECKS/input_cbm18a_free_bnd"  "0 1 2 4" "1 2"
job pb3d-fixed-x10  "$SCRIPT_DIR/run_deck.sbatch" "$DECKS/input_cbm18a_fixed_bnd" "0" "" 10
job pb3d-ctest      --wrap="source $PB3D_ROOT/env.sh; cd $PB3D_ROOT/build && ctest --output-on-failure --timeout 600"
job pb3d-ctest-dbg  --wrap="source $PB3D_ROOT/env.sh; cd $PB3D_ROOT/build-debug && ctest --output-on-failure --timeout 600 -L 'unit|smoke|fullstack'"

stamp "queue check: $(squeue -u $USER -h | wc -l) jobs of mine queued/running"
stamp "PIPELINE DONE in $(( $(date +%s) - t0 )) s"
