#!/bin/bash
# Physics regression test: fixed-boundary peeling-ballooning eigenvalue for
# the cbm18a case (extended pseudo-vacuum HELENA equilibrium).
#
# Usage: cbm18a_fixed_bnd.sh <PB3D-executable> <input-deck> <equilibrium-file>
#
# The equilibrium fixture (cbm18a_extended_vac.12, ~25MB, generated with the
# PB3D-patched HELENA from the original cbm18a input) is committed
# xz-compressed under tests/fixtures and decompressed at configure time.
#
# Anchor: omega^2/omega_A^2 = -4.074e-2 (MISHKA normalization), single
# process, established 2026-07-08 with PETSc 3.25.3/SLEPc 3.25.1, gfortran
# 16.1, macOS arm64. Asserted to +-10% to allow for platform/solver noise;
# tighten once more platforms have been sampled.
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3

ANCHOR=-4.074e-2
RTOL=0.10

# fresh scratch dir, as PB3D writes several output files
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cp "$INPUT" "$WORKDIR/input_deck"
cd "$WORKDIR"

# note: singleton MPI (no mpirun); see Documentation/testing.md
# One launch, no retries: a crash is a failure. (The intermittent macOS
# crashes of 2026-07 were a heap overrun in PSPLINE's r8genxpkg for 3-point
# splines, fixed by Libraries/pspline-r8genxpkg-nx3.patch; issue #4.)
if ! "$PB3D_EXE" input_deck "$EQ_FILE" > run.log 2>&1; then
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

echo "eigenvalue: $EV (anchor: $ANCHOR, rtol: $RTOL)"
awk -v ev="$EV" -v anchor="$ANCHOR" -v rtol="$RTOL" 'BEGIN {
    rel = (ev - anchor) / anchor
    if (rel < 0) rel = -rel
    if (rel <= rtol) { print "PASS: relative deviation", rel; exit 0 }
    else { print "FAIL: relative deviation", rel; exit 1 }
}'
