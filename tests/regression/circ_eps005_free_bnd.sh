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
# the fixture has a 1 % current peaking). PB3D's own value, established
# 2026-10-27 after step 3 (vacuum factor -2 pi), is the anchor at 1e-3; the
# MISHKA value must hold to 2 % (PB3D's radial discretization is about 1 %,
# the boundary-element response 1 to 2 %).
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3
PYTHON=${4:-}
CHECK=$(cd "$(dirname "$0")" && pwd)/check_vac_res.py

ANCHOR=-4.3034e-1
RTOL=1e-3
MISHKA=-0.42679
MISHKA_RTOL=0.02

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

echo "eigenvalues: fixed: $EV_FIXED, free: $EV_FREE (anchor: $ANCHOR, rtol: $RTOL; MISHKA: $MISHKA, rtol: $MISHKA_RTOL)"
awk -v evf="$EV_FIXED" -v ev="$EV_FREE" -v anchor="$ANCHOR" -v rtol="$RTOL" \
    -v mishka="$MISHKA" -v mrtol="$MISHKA_RTOL" 'BEGIN {
    ok = 1
    if (evf <= 0) { print "FAIL: fixed boundary not stable:", evf; ok = 0 }
    rel = (ev - anchor) / anchor; if (rel < 0) rel = -rel
    relm = (ev - mishka) / mishka; if (relm < 0) relm = -relm
    if (rel > rtol) { print "FAIL: deviation from anchor", rel; ok = 0 }
    if (relm > mrtol) { print "FAIL: deviation from MISHKA", relm; ok = 0 }
    if (ok) { print "PASS: deviations", rel, "(anchor)", relm, "(MISHKA)"; exit 0 }
    exit 1
}'
