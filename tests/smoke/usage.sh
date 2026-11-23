#!/bin/bash
# Smoke test: PB3D or POST prints the usage message when it is called without
# arguments and the help message when it is called with --help, and then every
# process stops with the expected exit status.
#
# Usage: usage.sh <executable> [launcher]
#   [launcher]  optional: the launcher with its process count, e.g.
#               "mpirun -n 2" (one argument, split on spaces on purpose)
#
# Required for both calls: the output contains "Usage:", and the exit status
# is 0 after --help, which the user asked for, and 2 without arguments, which
# is a mistake that a script has to notice (option A of issue #40).
#
# Under a launcher with more than one process, all processes have to stop.
# Only rank 0 reads the command line; the other processes learn from its
# broadcast that they have to stop (issue #40). A run in which they do not
# learn it never ends: ctest's TIMEOUT then ends the test, and it fails.
set -u

EXE=$1
LAUNCHER=${2:-}

# fresh scratch dir, as the programs write output files
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

OK=1
# run_case <description> <expected exit status> [arguments of the program]
run_case() {
    local what=$1 expected=$2
    shift 2
    $LAUNCHER "$EXE" "$@" > run.log 2>&1    # $LAUNCHER unquoted on purpose
    local status=$?
    local case_ok=1
    if [ "$status" -ne "$expected" ]; then
        echo "FAIL: $what returned exit status $status instead of $expected"
        case_ok=0
    fi
    if ! grep -qF "Usage:" run.log; then
        echo "FAIL: $what did not print the usage message"
        case_ok=0
    fi
    if [ "$case_ok" -ne 1 ]; then
        echo "log tail of $what:"
        tail -20 run.log
        OK=0
    fi
}

run_case "the call with --help" 0 --help
run_case "the call without arguments" 2

if [ "$OK" -ne 1 ]; then
    exit 1
fi
echo "PASS: usage and help message, exit status 0 after --help and 2 without arguments"
