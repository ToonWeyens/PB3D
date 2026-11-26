#!/bin/bash
# Physics regression test: FREE-boundary peeling-ballooning eigenvalue for
# the cbm18a case, with the vacuum response (BC_style(2) = 2).
#
# Usage: cbm18a_free_bnd.sh <PB3D-executable> <input-deck> <equilibrium-file> [python3]
#
# The companion to cbm18a_fixed_bnd.sh: the same deck window (psi 0.45 to
# 0.75, where the n = 10 mode still has an edge amplitude of 0.1 of its
# peak), but the boundary is free and the vacuum response is computed by the
# boundary element method of vac_ops. The vacuum can only add energy, so the
# free-boundary eigenvalue must lie between the no-vacuum value (-4.6041e-2)
# and the fixed-boundary one (-4.3570e-2, cbm18a_fixed_bnd.sh); the test
# asserts the upper bound.
#
# Anchor (MISHKA normalization, single process, asserted to 1e-3 relative):
# -4.5049e-2 with the 16-digit fixture (23 November 2026, issue #65;
# -4.5051e-2 with the 8-digit one). On the circular external-kink fixtures the same free-boundary
# treatment agrees with MISHKA to 1.2 % or better (circ_eps005_free_bnd.sh).
# The 200-point value is 0.5 % from the radially converged -4.5263e-2. That
# error falls as 1/n_r_sol and does not come from the outer edge: the fixed
# case converges the same way.
#
# With a python3 (with h5py) as fourth argument, the vacuum response matrix
# is also checked for Hermiticity and positive semi-definiteness
# (check_vac_res.py).
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3
PYTHON=${4:-}
CHECK=$(cd "$(dirname "$0")" && pwd)/check_vac_res.py

ANCHOR=-4.5049e-2
RTOL=1e-3
FIXED=-4.3570e-2

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
    "$PB3D_EXE" "$deck" "$EQ_FILE" > run.log 2>&1
    STATUS=$?
    if [ "$STATUS" -ne 0 ]; then
        echo "FAIL: PB3D exited with status $STATUS; log tail:"
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

run_deck input_deck
if [ -n "$PYTHON" ]; then
    "$PYTHON" "$CHECK" PB3D_out.h5 || exit 1
fi

echo "eigenvalue: $EV (anchor: $ANCHOR, rtol: $RTOL, fixed: $FIXED)"
awk -v ev="$EV" -v a="$ANCHOR" -v rtol="$RTOL" -v fixed="$FIXED" 'BEGIN {
    ok = 1
    rel = (ev - a) / a; if (rel < 0) rel = -rel
    if (rel > rtol) { print "FAIL: deviation", rel; ok = 0 }
    if (ev > fixed) { print "FAIL: free boundary more stable than fixed:", ev, ">", fixed; ok = 0 }
    if (ok) { print "PASS: deviation", rel, "free below fixed"; exit 0 }
    exit 1
}'
