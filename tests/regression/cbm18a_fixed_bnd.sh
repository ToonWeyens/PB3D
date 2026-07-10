#!/bin/bash
# Physics regression test: fixed-boundary peeling-ballooning eigenvalue for
# the cbm18a case (extended pseudo-vacuum HELENA equilibrium).
#
# Usage: cbm18a_fixed_bnd.sh <PB3D-executable> <input-deck> <equilibrium-file>
#
# The equilibrium fixture (cbm18a_extended_vac.12, ~25MB, generated with the
# PB3D-patched HELENA from the original cbm18a input) is not committed to the
# repository; the test is only registered when PB3D_FIXTURE_DIR is passed to
# CMake and the file is present there.
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
#
# KNOWN FLAKY on macOS with OpenMPI 5 + parallel HDF5 1.14: runs
# intermittently die (SIGTRAP / failed H5Fcreate) in the parallel-I/O layer,
# at ~50% per launch on the machine where this was characterized (2026-07-08).
# The physics is unaffected: successful runs reproduce the eigenvalue to all
# digits. Retried here until the I/O survives; root cause is tracked as
# follow-up work (suspect: OpenMPI 5 I/O path on macOS; try MPICH).
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
