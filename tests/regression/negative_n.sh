#!/bin/bash
# Physics regression test: a negative toroidal mode number (issue #22).
#
# Usage: negative_n.sh <PB3D-executable> <input-deck> <equilibrium-file>
#
# With prim_X = -n every Fourier mode is the complex conjugate of the mode
# with prim_X = n and the opposite secondary mode number, so the eigenvalue
# must be that of n. The test runs the free-boundary cbm18a deck (HELENA,
# n = 10, 20 poloidal mode numbers) with n = 10 and with n = -10, both with
# tol_SLEPC = 1e-12:
#   - at fixed boundary with max_njq_change = 0.05. The perturbation grid
#     (X_grid_style 3) gets floor(|prim_X dq| / max_njq_change) points more
#     than the equilibrium grid in every interval of the equilibrium grid, dq
#     being the change of q over the interval: 368 points with 0.05, none
#     with the default 0.49. The test asks for the same number, larger than
#     zero, in both runs.
#   - at free boundary with the deck as it is, through the axisymmetric
#     vacuum.
# Measured on 25 November 2026: 3.0e-12 relative at fixed and 4e-13 at free
# boundary; the test asks for 1e-9. The build of 2.55 stops on every run
# with a negative prim_X, in calc_eqd_grid_3D ("The angular array has to
# have a length of at least 1").
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3

RTOL=1e-9

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

sed 's/tol_SLEPC               = 1.0E-6/tol_SLEPC               = 1.0E-12/' "$INPUT" > free_pos
sed 's/BC_style                = 1 2/BC_style                = 1 1\
    max_njq_change          = 0.05/' free_pos > fixed_pos
for deck in free fixed; do
    sed 's/prim_X                  = 10/prim_X                  = -10/' ${deck}_pos > ${deck}_neg
    if ! grep -q 'prim_X *= -10' ${deck}_neg; then
        echo "FAIL: the deck has no line prim_X = 10 to change"
        exit 1
    fi
done

# run one deck: the eigenvalue into EV, the number of added points into ADDED
run_deck() {
    rm -f PB3D_out*
    "$PB3D_EXE" "$1" "$EQ_FILE" > run.log 2>&1
    STATUS=$?
    if [ "$STATUS" -ne 0 ]; then
        echo "FAIL: PB3D exited with status $STATUS on $1; log tail:"
        tail -20 run.log
        exit 1
    fi
    EV=$(awk '$1 == 1 {print $2}' PB3D_out_EV_R_1.txt 2>/dev/null)
    if [ -z "$EV" ]; then
        echo "FAIL: no accepted eigenvalue for $1:"
        cat PB3D_out_EV_R_1.txt 2>/dev/null
        exit 1
    fi
    ADDED=$(awk '/ Added [0-9]+ points$/ {print $(NF-1)}' run.log)
}

run_deck fixed_pos
EV_FIXED=$EV; ADDED_POS=${ADDED:-0}
run_deck fixed_neg
EV_FIXED_NEG=$EV; ADDED_NEG=${ADDED:-0}
run_deck free_pos
EV_FREE=$EV
run_deck free_neg
EV_FREE_NEG=$EV

echo "fixed boundary: n = 10: $EV_FIXED ($ADDED_POS points added), n = -10: $EV_FIXED_NEG ($ADDED_NEG points added)"
echo "free boundary: n = 10: $EV_FREE, n = -10: $EV_FREE_NEG (rtol: $RTOL)"
awk -v ev="$EV_FIXED" -v evn="$EV_FIXED_NEG" -v evf="$EV_FREE" -v evfn="$EV_FREE_NEG" \
    -v ap="$ADDED_POS" -v an="$ADDED_NEG" -v rtol="$RTOL" 'BEGIN {
    ok = 1
    if (ap <= 0) { print "FAIL: no points added to the perturbation grid for n = 10"; ok = 0 }
    if (an != ap) { print "FAIL: points added for n = -10 and n = 10 differ:", an, ap; ok = 0 }
    rel = (evn - ev) / ev; if (rel < 0) rel = -rel
    relf = (evfn - evf) / evf; if (relf < 0) relf = -relf
    if (rel > rtol) { print "FAIL: fixed boundary, n = -10 and n = 10 differ by", rel; ok = 0 }
    if (relf > rtol) { print "FAIL: free boundary, n = -10 and n = 10 differ by", relf; ok = 0 }
    if (ok) { print "PASS: n = -10 and n = 10 differ by", rel, "(fixed)", relf, "(free)"; exit 0 }
    exit 1
}'
