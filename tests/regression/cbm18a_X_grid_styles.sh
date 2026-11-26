#!/bin/bash
# Regression test: the perturbation grid and the solution grid hold the same
# matrices, also when a Fourier mode has another index in their tables.
#
# Usage: cbm18a_X_grid_styles.sh <PB3D-executable> <input-deck> <equilibrium-file>
#
# The integrated matrices PV and KV are calculated in the perturbation grid
# and used in the solution grid. With X_grid_style 2 the two grids are the
# same. With X_grid_style 3, the default for X_style 2, the perturbation grid
# is the equilibrium grid with points added, and interp_V interpolates the
# matrices to the solution grid.
#
# The deck is the cbm18a fixed-boundary deck with three changes:
#  - min_r_sol = 0.4552: the window of Fourier modes at the first equilibrium
#    point of the run starts at m = 5 and the one at the first solution point
#    at m = 6, so a Fourier mode has another index in the tables of the two
#    grids;
#  - a parallel range of 0.2 pi to 1.5 pi, which is not symmetric and makes
#    the matrices complex. The symmetric PV_0, PV_2, KV_0 and KV_2 have one
#    entry for the two complex conjugate elements of a pair of Fourier modes,
#    and the result then depends on which of the two every routine takes it
#    to be (the one whose row has the higher mode number, is_conj_X). With
#    real matrices, as in the other cbm18a decks, it does not;
#  - 100 normal points instead of 200, to keep the test short.
# The case is a numerical check and not a physical one, because the parallel
# range is not a whole turn.
#
# Asserted:
#  - the eigenvalues of X_grid_style 2 and 3 agree to 1e-5 relative. They
#    agree to 2e-7; before version 2.56 they differed by 1.2e-3, because the
#    two grids held different elements of some pairs of Fourier modes;
#  - X_grid_style 2 gives the anchor omega^2/omega_A^2 = -8.8610e-2 to 1e-3
#    relative (single process, established 2026-11-25 with PETSc 3.25.3 /
#    SLEPc 3.25.1, gfortran 16.1, macOS arm64; master 8ad111e gives the same
#    value to all digits for X_grid_style 2).
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3

ANCHOR=-8.8610e-2
RTOL=1e-3
RTOL_STYLES=1e-5

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

# The deck of one X_grid_style.
for style in 2 3; do
    sed -E -e 's/^( *min_par_X *= *).*/\10.2/' \
        -e 's/^( *max_par_X *= *).*/\11.5/' \
        -e 's/^( *min_r_sol *= *).*/\10.4552/' \
        -e 's/^( *n_r_sol *= *).*/\1100/' "$INPUT" | \
        awk -v s="$style" '{print} /^&inputdata_PB3D/ {print "    X_grid_style            = " s}' \
        > input_deck_$style
    for key in "min_par_X *= *0.2" "max_par_X *= *1.5" "min_r_sol *= *0.4552" \
        "n_r_sol *= *100" "X_grid_style *= *$style"; do
        if ! grep -Eq "^ *$key *$" input_deck_$style; then
            echo "FAIL: the deck for X_grid_style $style has no line \"$key\""
            exit 1
        fi
    done
done

# One launch per style, no retries: a crash is a failure.
for style in 2 3; do
    mkdir run_$style
    (cd run_$style && "$PB3D_EXE" ../input_deck_$style "$EQ_FILE" > run.log 2>&1)
    STATUS=$?
    if [ "$STATUS" -ne 0 ]; then
        echo "FAIL: PB3D exited with status $STATUS for X_grid_style $style; log tail:"
        tail -20 run_$style/run.log
        exit 1
    fi
    if [ ! -f run_$style/PB3D_out_EV_R_1.txt ]; then
        echo "FAIL: no eigenvalue output file for X_grid_style $style; log tail:"
        tail -20 run_$style/run.log
        exit 1
    fi
done

EV2=$(awk '$1 == 1 {print $2}' run_2/PB3D_out_EV_R_1.txt)
EV3=$(awk '$1 == 1 {print $2}' run_3/PB3D_out_EV_R_1.txt)
if [ -z "$EV2" ] || [ -z "$EV3" ]; then
    echo "FAIL: no accepted eigenvalue (X_grid_style 2: \"$EV2\", 3: \"$EV3\")"
    exit 1
fi

echo "eigenvalue with X_grid_style 2: $EV2 (anchor: $ANCHOR, rtol: $RTOL)"
echo "eigenvalue with X_grid_style 3: $EV3 (against X_grid_style 2, rtol: $RTOL_STYLES)"
awk -v e2="$EV2" -v e3="$EV3" -v anchor="$ANCHOR" -v rtol="$RTOL" -v rtol_s="$RTOL_STYLES" 'BEGIN {
    rel = (e2 - anchor) / anchor; if (rel < 0) rel = -rel
    rel_s = (e3 - e2) / e2; if (rel_s < 0) rel_s = -rel_s
    if (rel > rtol) { print "FAIL: X_grid_style 2 deviates from the anchor by", rel; exit 1 }
    if (rel_s > rtol_s) { print "FAIL: X_grid_style 2 and 3 differ by", rel_s; exit 1 }
    print "PASS: deviation from the anchor", rel, "and between the styles", rel_s
}'
