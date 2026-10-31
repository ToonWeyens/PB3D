#!/bin/bash
# build_deps.sh -- build PB3D's non-module dependencies on Snellius under
# ~/pb3d/deps: PETSc (complex), SLEPc, PSPLINE, minimal LIBSTELL.
#
# Usage: build_deps.sh [petsc|slepc|pspline|libstell|all]   (default: all)
#
# Runs on the login node with nice and at most 16 parallel jobs. Tarballs and
# build trees are removed after a successful install; only the install
# prefixes stay. Downloads are anonymous (generic User-Agent, no identity).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/../env.sh"

DEPS=$PB3D_ROOT/deps
LOGS=$PB3D_ROOT/logs
mkdir -p "$DEPS" "$LOGS"
NP=16
NICE="nice -n 10"
UA="Mozilla/5.0"

PETSC_VER=3.25.3
SLEPC_VER=3.25.1
PETSC_URL=https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-$PETSC_VER.tar.gz
SLEPC_URL=https://slepc.upv.es/download/distrib/slepc-$SLEPC_VER.tar.gz
PSPLINE_URL=https://github.com/ilonster/pspline.git

stage=${1:-all}

fetch() {  # fetch <url> <file>
    [ -f "$2" ] || curl -sSL -A "$UA" -o "$2" "$1"
}

build_petsc() {
    # MUMPS (with ParMETIS/METIS orderings) is required: PB3D's shift-invert
    # eigensolve factorizes a distributed (mpibaij) matrix with LU, and PETSc's
    # built-in LU is sequential only -- without MUMPS every run with more than
    # one MPI process fails with "Could not locate a solver type for
    # factorization type LU and matrix type mpibaij".
    echo "=== PETSc $PETSC_VER (complex, double, optimized, MUMPS) -> $PETSC_DIR ==="
    local prefix=$PETSC_DIR
    cd "$DEPS"
    fetch "$PETSC_URL" petsc-$PETSC_VER.tar.gz
    rm -rf petsc-$PETSC_VER
    tar xzf petsc-$PETSC_VER.tar.gz
    (   # configure/make need PETSC_DIR = the build tree, not the prefix
        cd petsc-$PETSC_VER
        export PETSC_DIR=$PWD PETSC_ARCH=arch-complex-opt
        $NICE ./configure \
            --prefix="$prefix" \
            --with-scalar-type=complex \
            --with-precision=double \
            --with-debugging=0 \
            --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            --with-blaslapack-lib="-L$EBROOTFLEXIBLAS/lib -lflexiblas" \
            --with-scalapack-lib="-L$EBROOTSCALAPACK/lib -lscalapack" \
            --download-mumps --download-parmetis --download-metis \
            --with-fortran-bindings=1 \
            --with-shared-libraries=1 \
            --with-x=0 \
            --with-make-np=$NP \
            COPTFLAGS="-O2" CXXOPTFLAGS="-O2" FOPTFLAGS="-O2"
        $NICE make all
        make install
    )
    cd "$DEPS"
    rm -rf petsc-$PETSC_VER petsc-$PETSC_VER.tar.gz
    pkg-config --libs PETSc
}

build_slepc() {
    echo "=== SLEPc $SLEPC_VER -> $SLEPC_DIR ==="
    local prefix=$SLEPC_DIR
    cd "$DEPS"
    fetch "$SLEPC_URL" slepc-$SLEPC_VER.tar.gz
    rm -rf slepc-$SLEPC_VER
    tar xzf slepc-$SLEPC_VER.tar.gz
    (   # configure/make need SLEPC_DIR = the build tree; PETSC_DIR = the installed prefix
        cd slepc-$SLEPC_VER
        export SLEPC_DIR=$PWD
        $NICE ./configure --prefix="$prefix"
        $NICE make MAKE_NP=$NP
        make install
    )
    cd "$DEPS"
    rm -rf slepc-$SLEPC_VER slepc-$SLEPC_VER.tar.gz
    pkg-config --libs SLEPc
}

build_pspline() {
    echo "=== PSPLINE -> $PSPLINE_DIR ==="
    cd "$DEPS"
    [ -d pspline ] || git -c http.userAgent="$UA" clone -q --depth 1 "$PSPLINE_URL" pspline
    # the nx = 3 heap-overrun guard from the PB3D tree (Libraries/), when present
    PSPLINE_PATCH=${PSPLINE_PATCH:-$PB3D_ROOT/src/Libraries/pspline-r8genxpkg-nx3.patch} \
        $NICE "$SCRIPT_DIR/build_pspline_linux.sh" "$DEPS/pspline"
    rm -rf "$DEPS/pspline/build_linux"
}

build_libstell() {
    echo "=== minimal LIBSTELL -> $LIBSTELL_DIR ==="
    # STELLOPT's LIBSTELL/Sources is expected at $DEPS/stellopt/LIBSTELL/Sources
    # (rsynced from a STELLOPT checkout; a full clone is not needed)
    STELLOPT_SRC=$DEPS/stellopt NETCDF_FORTRAN_DIR=$EBROOTNETCDFMINFORTRAN \
        $NICE "$SCRIPT_DIR/build_libstell_min_linux.sh" "$LIBSTELL_DIR"
}

case $stage in
    petsc)    build_petsc ;;
    slepc)    build_slepc ;;
    pspline)  build_pspline ;;
    libstell) build_libstell ;;
    all)      build_petsc; build_slepc; build_pspline; build_libstell ;;
    *) echo "unknown stage $stage"; exit 1 ;;
esac
echo "=== build_deps.sh $stage: done ==="
