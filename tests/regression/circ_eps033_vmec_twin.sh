#!/bin/bash
# Physics regression test: the HELENA fixture circ_eps033_q104 (circular,
# eps = 0.33, flat current, zero pressure, q from 1.041 to 1.038) and its VMEC
# twin wout_circ_eps033_q104.nc (issue #4, step 6e; made from the HELENA file
# with tests/reference/make_vmec_twin.sh: PB3D --export_HEL, then VMEC++) must
# give the same fixed-boundary eigenvalue for the n = 10 external kink, m = 11,
# which no vacuum enters: HELENA +7.84304e-2, VMEC +7.84306e-2 (2.4e-6 apart) on
# 30 October 2026, once PB3D's reference field B_0 is the (0,0) mode of |B| on
# the axis; with the old B_0 (all modes summed at theta = zeta = 0) the VMEC
# value was +7.9716e-2, 1.6 % off. The test asks for 1e-4 between the twins and
# anchors the VMEC value at 1e-3, as the regression layer's only fixed-boundary
# VMEC number, so that a change in the VMEC path shows.
#
# The twin is then run free-boundary through the 3-D vacuum (style 1:
# alpha_style 2, 64 field lines, the deck's 100 parallel points). The response
# of this first-order boundary-element method is 0.80 of the axisymmetric one
# at this resolution (converging to 1), so the eigenvalue, -0.51892 on
# 30 October 2026, lies 17 % below the HELENA style-2 value -0.44383; the test
# anchors it at 1e-3, asks for 25 % and for free below fixed. This pins the
# orientation of the normal that store_vac_VMEC hands to the kernels: with the
# outward normal the response was 0.27 of the axisymmetric one and the
# eigenvalue -0.083 (issue #4, step 6e).
#
# Usage: circ_eps033_vmec_twin.sh <PB3D-executable> <input-deck> <HELENA-file> <wout-file>
set -u

PB3D_EXE=$1
INPUT=$2
HEL_FILE=$3
WOUT_FILE=$4

ANCHOR_VMEC=7.8431e-2
RTOL=1e-3
TWIN_RTOL=1e-4
ANCHOR_FREE=-5.1892e-1
HEL_FREE=-4.4383e-1
FREE_RTOL=0.25
N_ALPHA=64

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

# free boundary through the 3-D vacuum: the deck with the field lines added
{
    echo '&inputdata_PB3D'
    echo '    alpha_style             = 2'
    echo "    n_alpha                 = $N_ALPHA"
    sed '1d; s/BC_style                = 1 1/BC_style                = 1 2/' input_fixed
} > input_free
run_deck input_free wout_eq.nc
EV_FREE=$EV

echo "fixed-boundary eigenvalues: HELENA: $EV_HEL, VMEC twin: $EV_VMEC (anchor: $ANCHOR_VMEC, rtol: $RTOL; twins within $TWIN_RTOL)"
echo "free-boundary eigenvalue of the twin through the 3-D vacuum, $N_ALPHA field lines: $EV_FREE (anchor: $ANCHOR_FREE, rtol: $RTOL; HELENA style 2: $HEL_FREE, rtol: $FREE_RTOL)"
awk -v evh="$EV_HEL" -v evv="$EV_VMEC" -v anchor="$ANCHOR_VMEC" -v rtol="$RTOL" -v trtol="$TWIN_RTOL" \
    -v evf="$EV_FREE" -v fanchor="$ANCHOR_FREE" -v hfree="$HEL_FREE" -v frtol="$FREE_RTOL" 'BEGIN {
    ok = 1
    if (evh <= 0) { print "FAIL: HELENA fixed boundary not stable:", evh; ok = 0 }
    if (evv <= 0) { print "FAIL: VMEC fixed boundary not stable:", evv; ok = 0 }
    rel = (evv - anchor) / anchor; if (rel < 0) rel = -rel
    relt = (evv - evh) / evh; if (relt < 0) relt = -relt
    if (rel > rtol) { print "FAIL: VMEC deviation from anchor", rel; ok = 0 }
    if (relt > trtol) { print "FAIL: HELENA and VMEC twins differ by", relt; ok = 0 }
    relf = (evf - fanchor) / fanchor; if (relf < 0) relf = -relf
    relh = (evf - hfree) / hfree; if (relh < 0) relh = -relh
    if (relf > rtol) { print "FAIL: free-boundary deviation from anchor", relf; ok = 0 }
    if (relh > frtol) { print "FAIL: 3-D vacuum result off the HELENA style-2 value by", relh; ok = 0 }
    if (evf >= evv) { print "FAIL: free boundary not below fixed:", evf, evv; ok = 0 }
    if (ok) { print "PASS: deviations", rel, "(anchor)", relt, "(twins)", relf, "(free anchor)", relh, "(3-D vacuum vs HELENA)"; exit 0 }
    exit 1
}'
