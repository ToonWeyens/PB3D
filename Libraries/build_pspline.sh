#!/bin/bash
# build_pspline.sh
# ----------------
# Build PSPLINE v2.5.0 (https://github.com/PrincetonUniversity/pspline) as the
# static library libpspline.a, without netCDF, and install it with its .mod
# files into <prefix>/lib and <prefix>/include, the layout that
# cmake/FindPSPLINE.cmake expects via PSPLINE_DIR.
#
# Usage: build_pspline.sh <pspline checkout> [<prefix>]   (default prefix:
#        <pspline checkout>/install)
#
#   git clone --depth 1 --branch v2.5.0 \
#       https://github.com/PrincetonUniversity/pspline.git ~/Code/pspline
#   Libraries/build_pspline.sh ~/Code/pspline
#
# PSPLINE's Makefile compiles its netCDF save and load routines when nf-config
# is on the PATH; PB3D does not use them, so the netCDF variables of the
# Makefile are set empty. Only the target libpspline.a is built: the default
# target also links the C++ test program czspline_test, which fails on macOS
# because clang++ does not find libgfortran.
set -e
[ "$(uname)" = Darwin ] && export PATH="/opt/homebrew/bin:$PATH"

SRC="$(cd "${1:?usage: build_pspline.sh <pspline checkout> [<prefix>]}" && pwd)"
PREFIX="${2:-$SRC/install}"
NP=${NP:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}   # parallel compile jobs

echo "=== PSPLINE $(git -C "$SRC" describe --tags 2>/dev/null || echo '(no git tag)') in $SRC ==="
make -C "$SRC" -j"$NP" NETCDF_FLAGS= FINCL= FLIBS= libpspline.a
make -C "$SRC" PREFIX="$PREFIX" install
echo "=== installed in $PREFIX ==="
ls "$PREFIX/lib" "$PREFIX/include"
