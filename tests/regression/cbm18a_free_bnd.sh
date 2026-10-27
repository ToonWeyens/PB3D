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
# The natural companion to cbm18a_fixed_bnd.sh: there, the *extended
# pseudo-vacuum* HELENA equilibrium mimics the vacuum with a cold-plasma
# layer and the boundary is held fixed; here, the plasma boundary is free
# and the vacuum response is computed by the boundary element method of
# vac_ops. The two growth rates should be comparable (the pseudo-vacuum
# construction was designed for exactly this correspondence).
#
# Anchor: omega^2/omega_A^2 = -4.5144e-2 (MISHKA normalization), single
# process, established 2026-10-27 after the HELENA shear fix of issue #4
# step 2 (before it: -4.2349e-2). Asserted to 1e-3 relative. NOTE: at
# max_r_sol = 1.0 the edge displacement is 1e-8 of the peak, so this deck
# cannot see the vacuum: fixed and free agree to 1e-8 here (issue #4,
# report section 9); step 4 moves the deck to the pedestal foot.
#
# The deck is run twice: with the BC style it specifies (4, explicit
# surface-energy minimization) and with BC style 2 (variational/Hermitian
# imposition of the same natural boundary condition). The two eigenvalues
# are asserted against the anchor and, mutually, to 1e-5 relative (they
# agree to ~3e-8 when established: two independent discretizations of the
# same physics).
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3
PYTHON=${4:-}
CHECK=$(cd "$(dirname "$0")" && pwd)/check_vac_res.py

ANCHOR=-4.5144e-2
RTOL=1e-3
MUTUAL_RTOL=1e-5

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

echo "eigenvalues: style 4: $EV_4, style 2: $EV_2 (anchor: $ANCHOR, rtol: $RTOL)"
awk -v ev4="$EV_4" -v ev2="$EV_2" -v anchor="$ANCHOR" -v rtol="$RTOL" \
    -v mrtol="$MUTUAL_RTOL" 'BEGIN {
    ok = 1
    rel4 = (ev4 - anchor) / anchor; if (rel4 < 0) rel4 = -rel4
    rel2 = (ev2 - anchor) / anchor; if (rel2 < 0) rel2 = -rel2
    relm = (ev4 - ev2) / anchor; if (relm < 0) relm = -relm
    if (rel4 > rtol) { print "FAIL: style 4 deviation", rel4; ok = 0 }
    if (rel2 > rtol) { print "FAIL: style 2 deviation", rel2; ok = 0 }
    if (relm > mrtol) { print "FAIL: mutual deviation", relm; ok = 0 }
    if (ok) { print "PASS: deviations", rel4, rel2, "mutual", relm; exit 0 }
    exit 1
}'
