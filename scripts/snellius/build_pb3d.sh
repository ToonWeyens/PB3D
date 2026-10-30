#!/bin/bash
# build_pb3d.sh -- configure and build PB3D on Snellius against the
# dependencies from build_deps.sh.
#
# Usage: build_pb3d.sh [release|debug|all]      (default: release)
#   release -> ~/pb3d/build        (CMAKE_BUILD_TYPE=Release)
#   debug   -> ~/pb3d/build-debug  (CMAKE_BUILD_TYPE=Debug, PB3D_ENABLE_DEBUG=ON)
#
# The source tree is expected in ~/pb3d/src (a clean `git archive` of the
# PB3D worktree). Builds on the login node with nice and at most 16 jobs.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/../env.sh"

SRC=$PB3D_ROOT/src
NP=16
NICE="nice -n 10"

configure_and_build() {   # configure_and_build <build-dir> <extra cmake args...>
    local bdir=$1; shift
    $NICE cmake -S "$SRC" -B "$bdir" \
        -DPSPLINE_DIR="$PSPLINE_DIR" \
        -DLIBSTELL_DIR="$LIBSTELL_DIR" \
        -DPB3D_FIXTURE_DIR="$PB3D_FIXTURE_DIR" \
        "$@"
    $NICE cmake --build "$bdir" -j $NP
    ls -la "$bdir/PB3D" "$bdir/POST"
}

case ${1:-release} in
    release) configure_and_build "$PB3D_ROOT/build" -DCMAKE_BUILD_TYPE=Release ;;
    debug)   configure_and_build "$PB3D_ROOT/build-debug" -DCMAKE_BUILD_TYPE=Debug -DPB3D_ENABLE_DEBUG=ON ;;
    all)     "$0" release; "$0" debug ;;
    *) echo "usage: $0 [release|debug|all]"; exit 1 ;;
esac
