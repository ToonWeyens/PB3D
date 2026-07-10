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
# Anchor: omega^2/omega_A^2 = -4.2349e-2 (MISHKA normalization), single
# process, established 2026-07-09 with distro PETSc 3.19/SLEPc 3.19
# (complex), gfortran 13, Ubuntu 24.04. Bit-identical between the STRUMPACK
# and ScaLAPACK vacuum solver paths; ~4% more unstable than the
# fixed-boundary anchor of the same deck window, as expected from freeing
# the boundary. Moves by ~1% when n_r_sol goes 200 -> 300 (normal radial
# convergence; the anchor is defined at the deck's resolution). Asserted to
# +-10% like the fixed-boundary anchor.
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

ANCHOR=-4.2349e-2
RTOL=0.10
MUTUAL_RTOL=1e-5

# fresh scratch dir, as PB3D writes several output files
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cp "$INPUT" "$WORKDIR/input_deck"
cd "$WORKDIR"

# run one deck and extract the eigenvalue into EV
# note: singleton MPI (no mpirun); retries for the known macOS MPI-I/O
# flakiness, see cbm18a_fixed_bnd.sh
run_deck() {
    local deck=$1
    local attempts=6
    local ok=0
    for attempt in $(seq 1 $attempts); do
        rm -f PB3D_out*
        "$PB3D_EXE" "$deck" "$EQ_FILE" > run.log 2>&1
        if [ -f PB3D_out_EV_R_1.txt ]; then
            ok=1
            [ "$attempt" -gt 1 ] && echo "NOTE: needed $attempt attempts (known macOS MPI flakiness)"
            break
        fi
        sleep 1
    done
    if [ "$ok" -ne 1 ]; then
        echo "FAIL: no eigenvalue output file produced in $attempts attempts; log tail:"
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
