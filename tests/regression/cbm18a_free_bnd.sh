#!/bin/bash
# Physics regression test: FREE-boundary peeling-ballooning eigenvalue for
# the cbm18a case, using the vacuum response (BC_style(2) = 4) instead of a
# fixed boundary.
#
# Usage: cbm18a_free_bnd.sh <PB3D-executable> <input-deck> <equilibrium-file>
#
# The natural companion to cbm18a_fixed_bnd.sh: there, the *extended
# pseudo-vacuum* HELENA equilibrium mimics the vacuum with a cold-plasma
# layer and the boundary is held fixed; here, the plasma boundary is free
# and the vacuum response is computed by the boundary element method of
# vac_ops. The two growth rates should be comparable (the pseudo-vacuum
# construction was designed for exactly this correspondence).
#
# Anchor: TBD_ANCHOR - to be established on the first verified run with the
# HELENA fixture; see Documentation/testing.md.
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3

ANCHOR=TBD_ANCHOR
RTOL=0.10

if [ "$ANCHOR" = "TBD_ANCHOR" ]; then
    echo "FAIL: anchor eigenvalue not yet established for this test."
    echo "Run manually, verify the physics, then substitute the anchor."
    exit 1
fi

# fresh scratch dir, as PB3D writes several output files
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cp "$INPUT" "$WORKDIR/input_deck"
cd "$WORKDIR"

# note: singleton MPI (no mpirun); retries for the known macOS MPI-I/O
# flakiness, see cbm18a_fixed_bnd.sh
ATTEMPTS=6
ok=0
for attempt in $(seq 1 $ATTEMPTS); do
    rm -f PB3D_out*
    "$PB3D_EXE" input_deck "$EQ_FILE" > run.log 2>&1
    if [ -f PB3D_out_EV_R_1.txt ]; then
        ok=1
        [ "$attempt" -gt 1 ] && echo "NOTE: needed $attempt attempts (known macOS MPI flakiness)"
        break
    fi
    sleep 1
done
if [ "$ok" -ne 1 ]; then
    echo "FAIL: no eigenvalue output file produced in $ATTEMPTS attempts; log tail:"
    tail -20 run.log
    exit 1
fi

EV=$(awk '$1 == 1 {print $2}' PB3D_out_EV_R_1.txt)
if [ -z "$EV" ]; then
    echo "FAIL: no accepted eigenvalue in output:"
    cat PB3D_out_EV_R_1.txt
    exit 1
fi

echo "eigenvalue: $EV (anchor: $ANCHOR, rtol: $RTOL)"
awk -v ev="$EV" -v anchor="$ANCHOR" -v rtol="$RTOL" 'BEGIN {
    rel = (ev - anchor) / anchor
    if (rel < 0) rel = -rel
    if (rel <= rtol) { print "PASS: relative deviation", rel; exit 0 }
    else { print "FAIL: relative deviation", rel; exit 1 }
}'
