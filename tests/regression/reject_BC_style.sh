#!/bin/bash
# Regression test: PB3D stops on an input value that it does not accept, says
# why, and returns an exit status that is not zero.
#
# Usage: reject_BC_style.sh <PB3D-executable> <input-deck> <equilibrium-file> [launcher]
#   <input-deck>  a deck with the line "BC_style = 1 2", such as
#                 input_cbm18a_free_bnd
#   [launcher]    optional: the launcher with its process count, e.g.
#                 "mpirun -n 2" (one argument, split on spaces on purpose)
#
# The test runs PB3D on a copy of the deck in which BC_style is 1 4. The
# second number is the boundary condition at the plasma edge: 1 is the fixed
# boundary, 2 the free boundary with vacuum, and no other value exists
# (style 4 was removed in November 2026, issue #36). The input check of
# adapt_run_PB3D (input_ops) refuses the deck with error code 1.
#
# Three things are required:
#  - the exit status is not zero (issue #40: it was 0 after every error,
#    because abort_MPI passed 0 to MPI_Abort);
#  - the output contains the message of the input check;
#  - no eigenvalue file exists, which means that PB3D did not go on.
# The output is redirected to a file, as a batch system does. PB3D has to
# write out its buffered output before it aborts, or the message is lost.
set -u

PB3D_EXE=$1
INPUT=$2
EQ_FILE=$3
LAUNCHER=${4:-}

MESSAGE="BC_style(2) has to be 1 (fixed boundary) or 2 (free boundary with vacuum)"

# fresh scratch dir, as PB3D writes several output files
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
sed 's/^\( *BC_style *= *\)1 2 *$/\11 4/' "$INPUT" > "$WORKDIR/input_deck"
cd "$WORKDIR"
if ! grep -q '^ *BC_style *= *1 4$' input_deck; then
    echo "FAIL: the deck $INPUT has no line 'BC_style = 1 2' to change"
    exit 1
fi

# one launch; without a launcher this is singleton MPI (see cbm18a_fixed_bnd.sh)
$LAUNCHER "$PB3D_EXE" input_deck "$EQ_FILE" > run.log 2>&1    # $LAUNCHER unquoted on purpose
STATUS=$?

OK=1
if [ "$STATUS" -eq 0 ]; then
    echo "FAIL: PB3D returned exit status 0 for a deck that it has to refuse"
    OK=0
fi
if ! grep -qF "$MESSAGE" run.log; then
    echo "FAIL: the output does not contain the message '$MESSAGE'"
    OK=0
fi
if [ -f PB3D_out_EV_R_1.txt ]; then
    echo "FAIL: PB3D wrote an eigenvalue file for a deck that it has to refuse"
    OK=0
fi
if [ "$OK" -ne 1 ]; then
    echo "log tail:"
    tail -20 run.log
    exit 1
fi
echo "PASS: PB3D refused BC_style = 1 4 with exit status $STATUS and the message of the input check"
