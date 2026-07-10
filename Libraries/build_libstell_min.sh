#!/bin/bash
# build_libstell_min.sh
# ---------------------
# Build the minimal LIBSTELL subset that PB3D needs (read_wout_mod + its
# transitive dependencies) as a small static library, laid out so that
# cmake/FindLIBSTELL.cmake picks it up:
#
#   $PREFIX/lib/libstell.a
#   $PREFIX/include/*.mod      (includes the marker files vsvd0.mod, vmec_input.mod)
#
# Validated on macOS (Apple Silicon), gfortran/mpif90 from Homebrew
# (GCC 16.1 + OpenMPI 5), netcdf-fortran from spack.
#
# Usage:
#   STELLOPT_SRC=~/Code/stellopt NETCDF_FORTRAN_DIR=<prefix> ./build_libstell_min.sh [PREFIX]
#
# Defaults below match the current machine.
set -euo pipefail
export PATH="/opt/homebrew/bin:$PATH"

STELLOPT_SRC=${STELLOPT_SRC:-$HOME/Code/stellopt}
SRC=$STELLOPT_SRC/LIBSTELL/Sources
# netcdf-fortran prefix: containing include/netcdf.inc and lib/libnetcdff.*
NETCDF_FORTRAN_DIR=${NETCDF_FORTRAN_DIR:-$(ls -d $HOME/spack/opt/spack/*/netcdf-fortran-* 2>/dev/null | head -1)}
PREFIX=${1:-$HOME/Code/PB3D/Libraries/libstell-min}

[ -f "$SRC/Modules/read_wout_mod.f90" ] || { echo "STELLOPT sources not found at $SRC"; exit 1; }
[ -f "$NETCDF_FORTRAN_DIR/include/netcdf.inc" ] || { echo "netcdf.inc not found under $NETCDF_FORTRAN_DIR/include"; exit 1; }

BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT
mkdir -p "$BUILD/mod" "$BUILD/obj" "$PREFIX/lib" "$PREFIX/include"
cd "$BUILD/obj"

# -DNETCDF: enable the netcdf wout reader (the only one PB3D needs)
# -DMPI_OPT: required -- mpi_params.f does not compile without it (upstream
#            always builds LIBSTELL with MPI); PB3D is an MPI code anyway.
# -fallow-argument-mismatch: F77-style netcdf calls in Ezcdf (GCC >= 10).
COMMON="-c -cpp -DNETCDF -DMPI_OPT -O2 -fPIC -fallow-argument-mismatch \
        -I$NETCDF_FORTRAN_DIR/include -I../mod -J../mod"
# .f sources are dual-format: fixed-form, and the DEFAULT 72-column limit is
# REQUIRED (free-form '&' continuations deliberately sit beyond column 72).
# Do NOT add -ffixed-line-length-none.
FIX="mpif90 $COMMON -ffixed-form"
FRE="mpif90 $COMMON -ffree-form -ffree-line-length-none"

# --- layer 0: kinds/constants/params (order matters) ---
$FIX $SRC/Modules/stel_kinds.f
$FIX $SRC/Modules/stel_constants.f
$FIX $SRC/Modules/vsvd0.f
$FIX $SRC/Modules/vparams.f
$FIX $SRC/Modules/mpi_params.f
$FIX $SRC/Modules/mpi_inc.f
$FRE $SRC/Modules/mpi_sharmem.f90
$FIX $SRC/Modules/v3_utilities.f
$FIX $SRC/Modules/vmec_input.f
$FIX $SRC/Modules/system_mod.f
$FRE $SRC/Modules/safe_open_mod.f90

# --- layer 1: ezcdf netcdf wrapper ---
$FRE $SRC/Ezcdf/handle_err.f90
$FRE $SRC/Ezcdf/ezcdf_inqvar.f90
$FRE $SRC/Ezcdf/ezcdf_opncls.f90
$FRE $SRC/Ezcdf/ezcdf_attrib.f90
$FRE $SRC/Ezcdf/ezcdf_GenGet.f90
$FRE $SRC/Ezcdf/ezcdf_GenPut.f90
$FRE $SRC/Ezcdf/ezcdf.f90

# --- layer 2: mgrid + read_wout ---
$FIX $SRC/Modules/mgrid_mod.f
$FRE $SRC/Modules/read_wout_mod.f90

# --- layer 3: external procedures referenced at link time ---
$FRE $SRC/Miscel/vmec_getenv.f90      # getenv interface used by mgrid_mod
$FIX $SRC/Miscel/parse_extension.f    # called by read_wout_mod/readw_and_open

ar -crs "$PREFIX/lib/libstell.a" ./*.o
cp ../mod/*.mod "$PREFIX/include/"

echo
echo "libstell (minimal) installed to: $PREFIX"
echo "  lib/libstell.a  ($(du -h "$PREFIX/lib/libstell.a" | cut -f1))"
echo "  include/: $(ls "$PREFIX/include" | wc -l | tr -d ' ') .mod files"
echo
echo "Configure PB3D with:  -DLIBSTELL_DIR=$PREFIX"
echo "Consumers must also link netcdf-fortran and LAPACK, e.g.:"
echo "  -L$NETCDF_FORTRAN_DIR/lib -lnetcdff -framework Accelerate"
