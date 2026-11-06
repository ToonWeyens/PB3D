#!/bin/bash
# run_deck.sbatch -- run one PB3D input deck on Snellius with a list of
# process counts, in one short single-node job on the rome partition.
#
# Submit (always from the login node, one job at a time):
#   sbatch ~/pb3d/scripts/run_deck.sbatch <deck> [np-list] [bc-style] [repeat] [exe] [keep] [eq-file]
#     <deck>     input deck, e.g. ~/pb3d/src/tests/regression/input_cbm18a_fixed_bnd
#     [np-list]  process counts, default "0 1 2 4"; 0 = MPI singleton (no mpirun)
#     [bc-style] e.g. "1 2" to override BC_style in the deck; "" keeps the deck's
#     [repeat]   run each process count this many times (default 1; crash statistics)
#     [exe]      executable (default ~/pb3d/build/PB3D; "" for the default)
#     [keep]     the word "keep" retains PB3D_out.h5 (about 270 MB per run, needed
#                for POST); by default it is deleted once the eigenvalue is read
#     [eq-file]  equilibrium file (default: the cbm18a fixture)
# All options are positional on purpose: Snellius submits with SBATCH_EXPORT=NONE,
# so exported shell variables do not reach the job.
#
# Each run gets its own directory ~/pb3d/runs/<job>-<jobid>/np<N>[_r<k>] with
# the deck copied to `input` and the equilibrium symlinked as `eq`, or under its
# own wout_*.nc name for a VMEC file, which PB3D recognises by the extension (PB3D has a
# 120-character limit on file paths). The first accepted eigenvalue of every
# run is printed as "EV np=<N> <value>" at the end of the job output.
#SBATCH --job-name=pb3d-deck
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --output=/home/%u/pb3d/logs/%x-%j.out
set -u
source "$HOME/pb3d/env.sh"

DECK=${1:?input deck}
NP_LIST=${2:-"0 1 2 4"}
BC_STYLE=${3:-}
REPEAT=${4:-1}
PB3D_EXE=${5:-$PB3D_ROOT/build/PB3D}
KEEP=${6:-}
EQ_FILE=${7:-$PB3D_FIXTURE_DIR/cbm18a_extended_vac.12}

JOBDIR=$PB3D_ROOT/runs/${SLURM_JOB_NAME:-manual}-${SLURM_JOB_ID:-$$}
mkdir -p "$JOBDIR"
echo "deck=$DECK exe=$PB3D_EXE eq=$EQ_FILE bc_style='${BC_STYLE}' np_list='$NP_LIST' repeat=$REPEAT"
echo "node=$(hostname) jobdir=$JOBDIR"
summary=()
for np in $NP_LIST; do
    for k in $(seq 1 "$REPEAT"); do
        rdir=$JOBDIR/np$np; [ "$REPEAT" -gt 1 ] && rdir=${rdir}_r$k
        mkdir -p "$rdir"; cd "$rdir"
        cp "$DECK" input
        [ -n "$BC_STYLE" ] && sed -i "s/^\( *BC_style *= *\).*/\1$BC_STYLE/" input
        EQ_LINK=eq; case "$EQ_FILE" in *.nc) EQ_LINK=$(basename "$EQ_FILE");; esac   # a VMEC file keeps its wout_*.nc name
        ln -sf "$EQ_FILE" "$EQ_LINK"
        t0=$(date +%s)
        if [ "$np" -eq 0 ]; then
            "$PB3D_EXE" input "$EQ_LINK" > run.log 2>&1; rc=$?
        else
            mpirun -np "$np" "$PB3D_EXE" input "$EQ_LINK" > run.log 2>&1; rc=$?
        fi
        dt=$(( $(date +%s) - t0 ))
        if [ -f PB3D_out_EV_R_1.txt ]; then
            ev=$(awk '$1 == 1 {print $2}' PB3D_out_EV_R_1.txt)
            line="EV np=$np run=$k rc=$rc time=${dt}s $ev"
        else
            line="EV np=$np run=$k rc=$rc time=${dt}s NO_EIGENVALUE_FILE (log tail: $(tail -3 run.log | tr '\n' '|'))"
        fi
        echo "$line"; summary+=("$line")
        [ "$KEEP" = keep ] || rm -f PB3D_out.h5 fort.*
    done
done
echo "=== summary ==="
printf '%s\n' "${summary[@]}"
