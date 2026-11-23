#!/bin/bash
# Smoke test: PB3D or POST stops when the input file on its command line does
# not exist, says so, and returns an exit status that is not zero.
#
# Usage: missing_input.sh <executable> <program name: PB3D or POST>
#
# There is no default input (issue #40). Before, the programs printed
# "No input file found. Default used", then read their input from a unit that
# was not open, so that the Fortran library created an empty file fort.50, and
# stopped with "Cannot read user-provided file".
#
# Four things are required:
#  - the exit status is not zero;
#  - the output names the missing file and links to the input section of the
#    tutorial;
#  - the line that names the calling routine names the program itself (POST
#    named itself "PB3D (main)");
#  - no file fort.50 exists.
# The second argument of the program does not exist either: the program has
# to stop at the input file, before it looks at the second file.
set -u

EXE=$1
PROG=$2

# fresh scratch dir, as the programs write output files
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

"$EXE" no_such_input no_such_file > run.log 2>&1
STATUS=$?

OK=1
if [ "$STATUS" -eq 0 ]; then
    echo "FAIL: $PROG returned exit status 0 for an input file that does not exist"
    OK=0
fi
for MESSAGE in 'Cannot open input file "no_such_input"' \
        'page_tutorial.html#tutorial_inputs' \
        ">> calling routine: $PROG (main)"; do
    if ! grep -qF "$MESSAGE" run.log; then
        echo "FAIL: the output does not contain '$MESSAGE'"
        OK=0
    fi
done
if [ -e fort.50 ]; then
    echo "FAIL: $PROG created fort.50"
    OK=0
fi
if [ "$OK" -ne 1 ]; then
    echo "log tail:"
    tail -20 run.log
    exit 1
fi
echo "PASS: $PROG stopped on the missing input file with exit status $STATUS"
