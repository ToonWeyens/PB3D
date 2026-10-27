#!/bin/bash
# Physics regression test: the HELENA shear S of the cbm18a equilibrium
# against the direct formula from the metric derivatives (issue #4, step 2).
#
# Usage: cbm18a_shear_HEL.sh <PB3D-executable> <input-deck> <equilibrium-file> <python3>
#
# Runs PB3D once (the fixed-boundary deck; the equilibrium quantities do not
# depend on the boundary condition) and hands PB3D_out.h5 to
# check_shear_HEL.py, which compares eq_2/S with
#   S = -(D_theta h^{12}/h^{22} - D_theta h^{22} h^{12}/(h^{22})^2)/J
# on the whole equilibrium grid and requires a non-zero outboard-midplane
# value. With the even spline condition on the odd Theta^alpha (the bug of
# June 2018) the deviation is 0.79 and S(theta = 0) is 1e-9 of max |S|.
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3
PYTHON=$4
CHECK=$(cd "$(dirname "$0")" && pwd)/check_shear_HEL.py

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cp "$INPUT" "$WORKDIR/input_deck"
cd "$WORKDIR"

# note: singleton MPI (no mpirun); one launch, a crash is a failure
if ! "$PB3D_EXE" input_deck "$EQ_FILE" > run.log 2>&1; then
    echo "FAIL: PB3D exited with status $?; log tail:"
    tail -20 run.log
    exit 1
fi
if [ ! -f PB3D_out.h5 ]; then
    echo "FAIL: no PB3D_out.h5 produced; log tail:"
    tail -20 run.log
    exit 1
fi

"$PYTHON" "$CHECK" PB3D_out.h5 1e-4
