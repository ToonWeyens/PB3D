#!/bin/bash
# Physics regression test: the HELENA fixture circ_eps033_q104 (circular,
# eps = 0.33, flat current, zero pressure, q from 1.041 to 1.038) and its VMEC
# twin wout_circ_eps033_q104.nc (issue #4, step 6e; made from the HELENA file
# with tests/reference/make_vmec_twin.sh: PB3D --export_HEL, then VMEC++) must
# give the same fixed-boundary eigenvalue for the n = 10 external kink, m = 11,
# which no vacuum enters: HELENA +7.84304e-2, VMEC +7.84306e-2 (2.4e-6 apart) on
# 30 October 2026, once PB3D's reference field B_0 is the (0,0) mode of |B| on
# the axis; with the old B_0 (all modes summed at theta = zeta = 0) the VMEC
# value was +7.9716e-2, 1.6 % off. With VMEC's half-mesh values (lambda, the
# Jacobian, B_u, B_v) at VMEC's positions (j - 3/2)/(ns - 1) (issue #16), the
# VMEC value is +7.84303e-2, 1.1e-6 below HELENA's. The test asks for 1e-4
# between the twins and anchors the VMEC value at 1e-3, as the regression
# layer's only fixed-boundary VMEC number, so that a change in the VMEC path
# shows.
#
# Between the fixed VMEC run and the free one, check_vmec_read.py compares what
# PB3D stored in PB3D_out.h5 with the wout file (issue #16): the file's mode
# lists against VMEC's conventions, PB3D's mode lists and sizes, R, Z and
# B_s, and the half-mesh coefficients against a cubic spline from VMEC's
# positions, to 1e-12 away from the two surfaces at each end. It fails when
# read_VMEC places the half mesh elsewhere or leaves B_s unset. It needs h5py
# and scipy, and is skipped with a note when they are missing.
#
# Then the same equilibrium with the magnetic field reversed (issue #19): in a
# copy of the wout file, reverse_field.py changes the sign of the fluxes and of
# the components of B and J, which changes neither the equilibrium nor its
# stability. The twin has a negative toroidal flux and a negative rotational
# transform; in the copy the two have opposite signs, so that the poloidal flux
# decreases outwards and the Jacobian of PB3D's flux coordinates is negative.
# With the poloidal and with the toroidal flux as normal coordinate, the copy
# must give the twin's fixed-boundary eigenvalue to the last digit, and
# check_vmec_read.py must find the copy stored as it is. These runs need the
# same python as the check and are skipped with it.
#
# The twin is then run free-boundary through the 3-D vacuum (style 1:
# alpha_style 2, 64 field lines, the deck's 100 parallel points). With the
# kernels integrated along the field lines in the near field and the G
# diagonal set by the linear-potential identity (31 October 2026), the
# response is 0.98 to 0.99 of the axisymmetric one at this resolution and the
# eigenvalue, -0.43914, lies 1.2 % below the HELENA style-2 value -0.43380
# (both with the trapezoidal weight of the free edge, 31 October 2026, later;
# before it -0.44910 and -0.44383, the same 1.2 % apart); the test anchors it
# at 1e-3, asks for 3 % and for free below fixed. Since 2 November 2026 the
# corner term of the axisymmetric vacuum (vac%style 2, calc_GH_2) carries the
# plasma-side angle, which moves the HELENA style-2 value at the fixture's 257
# poloidal points to -0.43679, 0.54 % from the twin's; the twin's value does
# not change (the 3-D path sets its diagonal by the row sum). Since 17
# November 2026 the last points of the free edge carry the weights of
# calc_end_weights in the normal integral, which moves the twin's value to
# -0.43901. (Before the kernel changes
# the response was 0.85 of the axisymmetric one and the eigenvalue
# -0.51892, 17 % below.) This pins the orientation of the normal that
# store_vac_VMEC hands to the kernels: with the outward normal the response
# was 0.27 of the axisymmetric one and the eigenvalue -0.083 (issue #4,
# step 6e).
#
# Last, the fixed run of the twin with n = -10 (issue #22): the modes are
# then the complex conjugates of those with n = 10, with m -> -m, so the
# eigenvalue must be that of n = 10. The modes come in the opposite order
# and SLEPc stops at another iterate, so the two agree only to the solver's
# tolerance: with the deck's tol_SLEPC of 1e-6 the eigenvalues differ by
# 9e-8 relative, with 1e-12 by 5e-15. The test asks for 1e-6. negative_n.sh
# has the same comparison on the HELENA path, with points added to the
# perturbation grid and at free boundary.
#
# Usage: circ_eps033_vmec_twin.sh <PB3D-executable> <input-deck> <HELENA-file> <wout-file> [python3]
set -u

PB3D_EXE=$1
INPUT=$2
HEL_FILE=$3
WOUT_FILE=$4
PYTHON=${5:-}
CHECK=$(cd "$(dirname "$0")" && pwd)/check_vmec_read.py
REVERSE=$(cd "$(dirname "$0")" && pwd)/reverse_field.py

ANCHOR_VMEC=7.8431e-2
RTOL=1e-3
TWIN_RTOL=1e-4
ANCHOR_FREE=-4.3901e-1
HEL_FREE=-4.3679e-1
FREE_RTOL=0.03
NEG_RTOL=1e-6
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
    "$PB3D_EXE" "$1" "$2" > run.log 2>&1
    STATUS=$?
    if [ "$STATUS" -ne 0 ]; then
        echo "FAIL: PB3D exited with status $STATUS on $2; log tail:"
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

same_ev() {   # same_ev <label> <eigenvalue> <eigenvalue>
    if [ "$2" != "$3" ]; then
        echo "FAIL: $1: eigenvalues differ: $2 and $3"
        exit 1
    fi
}

# what PB3D stored from the wout file, before the free run overwrites it, and
# the same equilibrium with the field reversed
if [ -z "$PYTHON" ]; then
    echo "NOTE: check_vmec_read.py and the reversed field skipped: no python3 with h5py"
elif ! "$PYTHON" -c 'import scipy' > /dev/null 2>&1; then
    echo "NOTE: check_vmec_read.py and the reversed field skipped: $PYTHON has no scipy"
else
    "$PYTHON" "$CHECK" PB3D_out.h5 wout_eq.nc || exit 1
    "$PYTHON" "$REVERSE" wout_eq.nc wout_rev.nc > /dev/null || exit 1
    run_deck input_fixed wout_rev.nc
    "$PYTHON" "$CHECK" PB3D_out.h5 wout_rev.nc > check_rev.log || { cat check_rev.log; exit 1; }
    same_ev "twin and the copy with the field reversed" "$EV_VMEC" "$EV"
    sed 's/use_pol_flux_F          = .true./use_pol_flux_F          = .false./' input_fixed > input_tor
    run_deck input_tor wout_eq.nc
    EV_TOR=$EV
    run_deck input_tor wout_rev.nc
    same_ev "toroidal flux as normal coordinate, twin and the copy with the field reversed" "$EV_TOR" "$EV"
    echo "reversed field: the copy gives the twin's eigenvalue $EV_VMEC, and with the toroidal flux $EV_TOR, to the last digit"
fi

# free boundary through the 3-D vacuum: the deck with the field lines added
{
    echo '&inputdata_PB3D'
    echo '    alpha_style             = 2'
    echo "    n_alpha                 = $N_ALPHA"
    sed '1d; s/BC_style                = 1 1/BC_style                = 1 2/' input_fixed
} > input_free
run_deck input_free wout_eq.nc
EV_FREE=$EV

# the fixed run with n = -10
sed 's/prim_X                  = 10/prim_X                  = -10/' input_fixed > input_fixed_neg
run_deck input_fixed_neg wout_eq.nc
EV_NEG=$EV

echo "fixed-boundary eigenvalues: HELENA: $EV_HEL, VMEC twin: $EV_VMEC (anchor: $ANCHOR_VMEC, rtol: $RTOL; twins within $TWIN_RTOL)"
echo "free-boundary eigenvalue of the twin through the 3-D vacuum, $N_ALPHA field lines: $EV_FREE (anchor: $ANCHOR_FREE, rtol: $RTOL; HELENA style 2: $HEL_FREE, rtol: $FREE_RTOL)"
echo "fixed-boundary eigenvalue of the twin with n = -10: $EV_NEG (equal to n = 10 within $NEG_RTOL)"
awk -v evh="$EV_HEL" -v evv="$EV_VMEC" -v anchor="$ANCHOR_VMEC" -v rtol="$RTOL" -v trtol="$TWIN_RTOL" \
    -v evf="$EV_FREE" -v fanchor="$ANCHOR_FREE" -v hfree="$HEL_FREE" -v frtol="$FREE_RTOL" \
    -v evn="$EV_NEG" -v nrtol="$NEG_RTOL" 'BEGIN {
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
    reln = (evn - evv) / evv; if (reln < 0) reln = -reln
    if (reln > nrtol) { print "FAIL: fixed boundary, n = -10 and n = 10 differ by", reln; ok = 0 }
    if (ok) { print "PASS: deviations", rel, "(anchor)", relt, "(twins)", relf, "(free anchor)", relh, "(3-D vacuum vs HELENA)", reln, "(n = -10)"; exit 0 }
    exit 1
}'
