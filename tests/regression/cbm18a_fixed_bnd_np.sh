#!/bin/bash
# Multi-process consistency: the cbm18a fixed-boundary deck run with one
# process and under an MPI launcher must give the same eigenvalue.
#
# Usage: cbm18a_fixed_bnd_np.sh <PB3D-executable> <input-deck> <equilibrium-file> <launcher>
#   <launcher>  the launcher with its process count, e.g. "mpirun -n 2"
#               (one argument, split on spaces on purpose)
#
# Why 1e-6: until 30 October 2026 two and four processes shifted every
# eigenvalue by 2e-5 relative, inside the 1e-3 anchor tolerance of
# cbm18a_fixed_bnd.sh, because redistribute_var (MPI_utilities) filled the
# overlap of two process ranges from the lower rank's ghost region (issue #4,
# step 0 open items). With the owner rule the shift is 1.5e-8 on Snellius
# (PETSc 3.25.3 with MUMPS, gfortran 14), and the SLEPc tolerance series put
# the reproducibility of the one-process eigenvalue itself at 6e-9.
#
# More than one process needs a distributed LU inside PETSc (MUMPS or
# SuperLU_DIST); without one EPSSolve aborts with "Could not locate a solver
# type for factorization type LU and matrix type mpibaij". The row is
# therefore opt-in: -DPB3D_TEST_MPI_REGRESSION=ON (see Documentation/testing.md).
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3
LAUNCHER=$4

RTOL=1e-6

# fresh scratch dir with one subdirectory per run, as PB3D writes several
# output files
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cp "$INPUT" "$WORKDIR/input_deck"
cd "$WORKDIR"

run_one() {   # run_one <subdir> [launcher words...]; leaves the eigenvalue in <subdir>/ev.txt
    local sub=$1; shift
    mkdir -p "$sub"
    cp input_deck "$sub/"
    cd "$sub"
    if ! "$@" "$PB3D_EXE" input_deck "$EQ_FILE" > run.log 2>&1; then
        echo "FAIL: PB3D ($sub) exited with status $?; log tail:"
        tail -20 run.log
        exit 1
    fi
    if [ ! -f PB3D_out_EV_R_1.txt ]; then
        echo "FAIL: no eigenvalue output file ($sub); log tail:"
        tail -20 run.log
        exit 1
    fi
    awk '$1 == 1 {print $2}' PB3D_out_EV_R_1.txt > ev.txt
    cd ..
}

run_one single
run_one multi $LAUNCHER          # unquoted on purpose
EV1=$(cat single/ev.txt)
EV2=$(cat multi/ev.txt)
if [ -z "$EV1" ] || [ -z "$EV2" ]; then
    echo "FAIL: no accepted eigenvalue (single: '$EV1', under '$LAUNCHER': '$EV2')"
    exit 1
fi

echo "eigenvalue with one process: $EV1, under '$LAUNCHER': $EV2 (rtol: $RTOL)"
awk -v a="$EV1" -v b="$EV2" -v rtol="$RTOL" 'BEGIN {
    rel = (b - a) / a
    if (rel < 0) rel = -rel
    if (rel <= rtol) { print "PASS: relative deviation", rel; exit 0 }
    else { print "FAIL: relative deviation", rel; exit 1 }
}'
