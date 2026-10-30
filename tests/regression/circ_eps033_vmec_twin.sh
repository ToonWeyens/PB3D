#!/bin/bash
# Physics regression test: the HELENA fixture circ_eps033_q104 (circular,
# eps = 0.33, flat current, zero pressure, q from 1.041 to 1.038) and its VMEC
# twin wout_circ_eps033_q104.nc (issue #4, step 6e; made from the HELENA file
# with tests/reference/make_vmec_twin.sh: PB3D --export_HEL, then VMEC++) must
# give the same fixed-boundary eigenvalue for the n = 10 external kink, m = 11,
# which no vacuum enters. The two equilibrium codes differ in the Fourier fit of
# the boundary and in the radial coordinate (VMEC's toroidal flux against
# HELENA's poloidal flux), so the agreement is 1.6 % (HELENA +7.8430e-2, VMEC
# +7.9716e-2 on 30 October 2026); the test asks for 2.5 %. The VMEC value is
# also anchored at 1e-3, as the regression layer's only fixed-boundary VMEC
# number, so that a change in the VMEC path shows.
#
# Usage: circ_eps033_vmec_twin.sh <PB3D-executable> <input-deck> <HELENA-file> <wout-file>
set -u

PB3D_EXE=$1
INPUT=$2
HEL_FILE=$3
WOUT_FILE=$4

ANCHOR_VMEC=7.9716e-2
RTOL=1e-3
TWIN_RTOL=0.025

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cp "$INPUT" "$WORKDIR/input_deck"
cd "$WORKDIR"
sed 's/BC_style                = 1 2/BC_style                = 1 1/' input_deck > input_fixed
# PB3D recognises a VMEC file by its .nc extension and limits path lengths,
# so both equilibria are reached through short links
ln -s "$HEL_FILE" eq_HEL
ln -s "$WOUT_FILE" wout_eq.nc

run_deck() {   # run_deck <deck> <equilibrium>
    rm -f PB3D_out*
    if ! "$PB3D_EXE" "$1" "$2" > run.log 2>&1; then
        echo "FAIL: PB3D exited with status $? on $2; log tail:"
        tail -20 run.log
        exit 1
    fi
    EV=$(awk '$1 == 1 {print $2}' PB3D_out_EV_R_1.txt 2>/dev/null)
    if [ -z "$EV" ]; then
        echo "FAIL: no accepted eigenvalue for $2:"
        cat PB3D_out_EV_R_1.txt 2>/dev/null
        exit 1
    fi
}

run_deck input_fixed eq_HEL
EV_HEL=$EV
run_deck input_fixed wout_eq.nc
EV_VMEC=$EV

echo "fixed-boundary eigenvalues: HELENA: $EV_HEL, VMEC twin: $EV_VMEC (anchor: $ANCHOR_VMEC, rtol: $RTOL; twins within $TWIN_RTOL)"
awk -v evh="$EV_HEL" -v evv="$EV_VMEC" -v anchor="$ANCHOR_VMEC" -v rtol="$RTOL" -v trtol="$TWIN_RTOL" 'BEGIN {
    ok = 1
    if (evh <= 0) { print "FAIL: HELENA fixed boundary not stable:", evh; ok = 0 }
    if (evv <= 0) { print "FAIL: VMEC fixed boundary not stable:", evv; ok = 0 }
    rel = (evv - anchor) / anchor; if (rel < 0) rel = -rel
    relt = (evv - evh) / evh; if (relt < 0) relt = -relt
    if (rel > rtol) { print "FAIL: VMEC deviation from anchor", rel; ok = 0 }
    if (relt > trtol) { print "FAIL: HELENA and VMEC twins differ by", relt; ok = 0 }
    if (ok) { print "PASS: deviations", rel, "(anchor)", relt, "(twins)"; exit 0 }
    exit 1
}'
