#!/bin/bash
# Physics regression test: FREE-boundary peeling-ballooning eigenvalue for
# the cbm18a case, using the vacuum response (BC_style(2) = 4) instead of a
# fixed boundary.
#
# Usage: cbm18a_free_bnd.sh <PB3D-executable> <input-deck> <equilibrium-file> [python3]
#
# With a python3 (with h5py) as fourth argument, the vacuum response matrix
# of the style-2 run is also checked for Hermiticity and positive
# semi-definiteness (check_vac_res.py; issue #4 step 3).
#
# The natural companion to cbm18a_fixed_bnd.sh: same deck window (psi 0.45 to
# 0.75, where the n = 10 mode still has an edge amplitude of 0.1 of its
# peak), but the boundary is free and the vacuum response is computed by the
# boundary element method of vac_ops. The vacuum can only add energy, so the
# free-boundary eigenvalue must lie between the no-vacuum value (-4.6041e-2,
# measured with a diagnostic build in issue #4 step 3) and the fixed-boundary
# one (-4.3573e-2, cbm18a_fixed_bnd.sh); the test asserts the upper bound.
#
# Anchors (MISHKA normalization, single process, established 2026-10-27 after
# issue #4 steps 2 and 3, asserted to 1e-3 relative):
#   BC_style(2) = 2, variational (Hermitian) imposition: -4.5019e-2. This is
#     the validated style: on the circular external-kink fixture it agrees
#     with MISHKA to 1 % (circ_eps005_free_bnd.sh).
#   BC_style(2) = 4, explicit natural-boundary-condition row: -4.4712e-2.
#     NOT validated: it differs from style 2 by 0.7 % here and by much more on
#     the circular case, independently of the radial resolution (issue #4,
#     step 4); it is anchored only so that a change shows up. Use style 2.
#
# The deck specifies style 4; the script runs it as given and then with
# style 2. With a python3 (with h5py) as fourth argument, the vacuum response
# matrix of the style-2 run is also checked for Hermiticity and positive
# semi-definiteness (check_vac_res.py; issue #4 step 3).
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3
PYTHON=${4:-}
CHECK=$(cd "$(dirname "$0")" && pwd)/check_vac_res.py

ANCHOR_2=-4.5019e-2
ANCHOR_4=-4.4712e-2
RTOL=1e-3
FIXED=-4.3573e-2

# fresh scratch dir, as PB3D writes several output files
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cp "$INPUT" "$WORKDIR/input_deck"
cd "$WORKDIR"

# run one deck and extract the eigenvalue into EV
# note: singleton MPI (no mpirun); one launch, no retries: a crash is a
# failure (see cbm18a_fixed_bnd.sh)
run_deck() {
    local deck=$1
    rm -f PB3D_out*
    if ! "$PB3D_EXE" "$deck" "$EQ_FILE" > run.log 2>&1; then
        echo "FAIL: PB3D exited with status $?; log tail:"
        tail -20 run.log
        exit 1
    fi
    if [ ! -f PB3D_out_EV_R_1.txt ]; then
        echo "FAIL: no eigenvalue output file produced; log tail:"
        tail -20 run.log
        exit 1
    fi
    EV=$(awk '$1 == 1 {print $2}' PB3D_out_EV_R_1.txt)
    if [ -z "$EV" ]; then
        echo "FAIL: no accepted eigenvalue in output:"
        cat PB3D_out_EV_R_1.txt
        exit 1
    fi
}

# style as specified in the deck (4)
run_deck input_deck
EV_4=$EV

# variational/Hermitian style 2
sed 's/BC_style                = 1 4/BC_style                = 1 2/' input_deck > input_deck_bc2
run_deck input_deck_bc2
EV_2=$EV
if [ -n "$PYTHON" ]; then
    "$PYTHON" "$CHECK" PB3D_out.h5 || exit 1
fi

echo "eigenvalues: style 4: $EV_4 (anchor: $ANCHOR_4), style 2: $EV_2 (anchor: $ANCHOR_2, rtol: $RTOL, fixed: $FIXED)"
awk -v ev4="$EV_4" -v ev2="$EV_2" -v a4="$ANCHOR_4" -v a2="$ANCHOR_2" \
    -v rtol="$RTOL" -v fixed="$FIXED" 'BEGIN {
    ok = 1
    rel4 = (ev4 - a4) / a4; if (rel4 < 0) rel4 = -rel4
    rel2 = (ev2 - a2) / a2; if (rel2 < 0) rel2 = -rel2
    if (rel4 > rtol) { print "FAIL: style 4 deviation", rel4; ok = 0 }
    if (rel2 > rtol) { print "FAIL: style 2 deviation", rel2; ok = 0 }
    if (ev2 > fixed) { print "FAIL: free boundary more stable than fixed:", ev2, ">", fixed; ok = 0 }
    if (ok) { print "PASS: deviations", rel4, rel2, "free below fixed"; exit 0 }
    exit 1
}'
