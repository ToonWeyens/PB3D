#!/bin/bash
# Physics regression test: the external kink of a circular, flat-current,
# zero-pressure, large-aspect-ratio (eps = 0.05) HELENA equilibrium with
# q_a = 1.038 (issue #4, step 0 fixture circ_eps005_q104), n = 10, mode
# m = 11: no resonance inside, so the fixed boundary is stable and the free
# boundary (vacuum response, BC_style(2) = 2, no wall) is unstable.
#
# Usage: circ_eps005_free_bnd.sh <PB3D-executable> <input-deck> <equilibrium-file> [python3]
#
# References (tests/reference/mishka_circular.txt): MISHKA-1 on the same
# equilibrium without a wall gives -omega^2/omega_A^2 = 0.42679, and the
# cylinder formula of tests/reference/cyl_kink.py 0.4366 (uniform current;
# the fixture has a 1 % current peaking). PB3D's own value is the anchor at
# 1e-3: -4.3034e-1 from 2026-10-27 (step 3, vacuum factor -2 pi) until the
# free edge got its trapezoidal quadrature weight on 31 October 2026, then
# -4.2168e-1, which is 1.2 % below MISHKA, and since 2 November 2026, when the
# corner term of the axisymmetric vacuum (vac%style 2, calc_GH_2) got the
# plasma-side angle, -4.2452e-1, 0.53 % below MISHKA; at n = 15 and 20 the gap
# is 0.24 and 0.13 % (0.55 and 0.26 % before), falling as 1/n^2 like the
# finite-n terms the high-n ordering drops. The MISHKA value must hold to 2 %.
# The deck is also run with 400 radial points: the eigenvalue must move by
# less than 5e-4 relative, which holds for the second-order convergence the
# edge weight gives (1.5e-4) and fails for the first-order one it replaced
# (1.2e-2).
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3
PYTHON=${4:-}
CHECK=$(cd "$(dirname "$0")" && pwd)/check_vac_res.py

ANCHOR=-4.2452e-1
RTOL=1e-3
MISHKA=-0.42679
MISHKA_RTOL=0.02
ORDER_RTOL=5e-4

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cp "$INPUT" "$WORKDIR/input_deck"
cd "$WORKDIR"

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

# fixed boundary: must be stable
sed 's/BC_style                = 1 2/BC_style                = 1 1/' input_deck > input_deck_fixed
run_deck input_deck_fixed
EV_FIXED=$EV

# free boundary with the vacuum response
run_deck input_deck
EV_FREE=$EV
if [ -n "$PYTHON" ]; then
    "$PYTHON" "$CHECK" PB3D_out.h5 || exit 1
fi

# the same with twice the radial points: second-order convergence of the free edge
sed -E 's/^( *n_r_sol *= *)200/\1400/' input_deck > input_deck_r400
run_deck input_deck_r400
EV_FREE_400=$EV

echo "eigenvalues: fixed: $EV_FIXED, free: $EV_FREE (anchor: $ANCHOR, rtol: $RTOL; MISHKA: $MISHKA, rtol: $MISHKA_RTOL), free at 400 points: $EV_FREE_400 (rtol: $ORDER_RTOL)"
awk -v evf="$EV_FIXED" -v ev="$EV_FREE" -v anchor="$ANCHOR" -v rtol="$RTOL" \
    -v mishka="$MISHKA" -v mrtol="$MISHKA_RTOL" -v ev4="$EV_FREE_400" -v ortol="$ORDER_RTOL" 'BEGIN {
    ok = 1
    relo = (ev4 - ev) / ev; if (relo < 0) relo = -relo
    if (relo > ortol) { print "FAIL: 200- and 400-point eigenvalues differ by", relo, "(first-order edge quadrature?)"; ok = 0 }
    if (evf <= 0) { print "FAIL: fixed boundary not stable:", evf; ok = 0 }
    rel = (ev - anchor) / anchor; if (rel < 0) rel = -rel
    relm = (ev - mishka) / mishka; if (relm < 0) relm = -relm
    if (rel > rtol) { print "FAIL: deviation from anchor", rel; ok = 0 }
    if (relm > mrtol) { print "FAIL: deviation from MISHKA", relm; ok = 0 }
    if (ok) { print "PASS: deviations", rel, "(anchor)", relm, "(MISHKA)"; exit 0 }
    exit 1
}'
