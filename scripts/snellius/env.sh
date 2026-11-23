#!/bin/bash
# env.sh -- PB3D build/run environment on Snellius (source this file).
# Everything lives under ~/pb3d; nothing is installed outside it, and
# `rm -rf ~/pb3d` removes all of it.
export PB3D_ROOT=$HOME/pb3d

module purge
module load 2025
module load foss/2025a                       # GCC 14.2.0, OpenMPI 5.0.7, FlexiBLAS 3.4.5 (OpenBLAS 0.3.29), ScaLAPACK 2.2.2
module load HDF5/1.14.6-gompi-2025a          # parallel HDF5 with Fortran bindings
module load netCDF-Fortran/4.6.2-gompi-2025a
module load CMake/3.31.3-GCCcore-14.2.0
module load pkgconf/2.3.0-GCCcore-14.2.0

# prefix installs (no PETSC_ARCH); the build trees are deleted after install.
# PB3D_PETSC_DIR and PB3D_SLEPC_DIR, when set before this file is sourced,
# select another pair of prefixes, for example a newer release built next to
# the working one (see build_deps.sh).
export PETSC_DIR=${PB3D_PETSC_DIR:-$PB3D_ROOT/deps/petsc}
unset PETSC_ARCH
export SLEPC_DIR=${PB3D_SLEPC_DIR:-$PB3D_ROOT/deps/slepc}
export PSPLINE_VERSION=v2.5.0                 # release of github.com/PrincetonUniversity/pspline
export PSPLINE_DIR=$PB3D_ROOT/deps/pspline-$PSPLINE_VERSION/install
export LIBSTELL_DIR=$PB3D_ROOT/deps/libstell-min
export PB3D_FIXTURE_DIR=$PB3D_ROOT/fixtures

# pkg-config: PETSc.pc / SLEPc.pc from the prefixes; netcdf-fortran.pc and
# scalapack.pc come from the modules' own PKG_CONFIG_PATH entries
export PKG_CONFIG_PATH=$PETSC_DIR/lib/pkgconfig:$SLEPC_DIR/lib/pkgconfig:${PKG_CONFIG_PATH:-}

# CMake >= 4 from a venv under ~/pb3d (pip wheel): the site CMake 3.31's
# Fortran dependency scanner loops forever on PETSc 3.25's finclude headers
# (petscsys.h <-> petscmath.h); the Mac builds use CMake 4.3 without trouble.
# Created with: python3 -m venv ~/pb3d/deps/venv && ~/pb3d/deps/venv/bin/pip install cmake
[ -x "$PB3D_ROOT/deps/venv/bin/cmake" ] && export PATH=$PB3D_ROOT/deps/venv/bin:$PATH

export OMP_NUM_THREADS=1
