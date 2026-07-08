#!/bin/bash
# build_pspline.sh
# ----------------
# Build PSPLINE (NTCC) with gfortran, without a real netcdf, from a checkout
# of https://github.com/ilonster/pspline. Copy this script into the pspline
# source root (or run it from there) -- it locates the sources relative to
# its own path. Produces install/lib/libpspline.a + install/include/*.mod,
# the layout expected by cmake/FindPSPLINE.cmake via PSPLINE_DIR.
# The ezcdf modules are compiled against the dummy netcdf.inc from cdf_dummy
# so that ezspline_save/load interfaces compile; netcdf symbols are only
# referenced by archive members that PB3D never pulls in.
set -e
export PATH="/opt/homebrew/bin:$PATH"

ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD="$ROOT/build_macos"
INSTALL="$ROOT/install"
FC=gfortran
FFLAGS="-O2 -fPIC -fallow-argument-mismatch -fallow-invalid-boz -std=legacy"

rm -rf "$BUILD"
mkdir -p "$BUILD" "$INSTALL/lib" "$INSTALL/include"
cd "$BUILD"

echo "=== 1. ezcdf modules (with dummy netcdf.inc) ==="
for f in ezcdf_inqvar ezcdf_opncls ezcdf_attrib ezcdf_genget ezcdf_genput ezcdf handle_err; do
    echo "  FC $f.f90"
    $FC $FFLAGS -I"$ROOT/cdf_dummy" -c "$ROOT/ezcdf/$f.f90" -o "$f.o"
done

echo "=== 2. pspline fixed-form sources ==="
EXCLUDE_F="pspltest.f pspltsub.f lookup_test.f r8lookup_test.f f2test.f f3test.f"
for p in "$ROOT"/pspline/*.f; do
    f=$(basename "$p")
    skip=0
    for e in $EXCLUDE_F; do [ "$f" = "$e" ] && skip=1; done
    [ $skip = 1 ] && continue
    echo "  FC $f"
    $FC $FFLAGS -c "$p" -o "$BUILD/${f%.f}.o"
done

echo "=== 3. ezspline modules and routines ==="
# module file first (needs ezcdf.mod for interface bodies)
$FC $FFLAGS -I"$BUILD" -c "$ROOT/pspline/ezspline_mod.f90" -o ezspline_mod.o
# remaining f90 files (all depend at most on ezspline_obj/ezspline/ezcdf/pspline_calls)
EXCLUDE_F90="ezspline_mod.f90 ezspline_test_r4.f90 ezspline_test_r8.f90 ezspline_perf_r4.f90 ezspline_perf_r8.f90 ezspline_io_test.f90 qk_pspline.f90"
for p in "$ROOT"/pspline/*.f90; do
    f=$(basename "$p")
    skip=0
    for e in $EXCLUDE_F90; do [ "$f" = "$e" ] && skip=1; done
    [ $skip = 1 ] && continue
    echo "  FC $f"
    $FC $FFLAGS -I"$BUILD" -c "$p" -o "$BUILD/${f%.f90}.o"
done

echo "=== 4. czspline C-API wrappers (.F90, need cpp + headers) ==="
$FC $FFLAGS -I"$BUILD" -I"$ROOT/include" -c "$ROOT/pspline/czspline_pointer_types.F90" -o "$BUILD/czspline_pointer_types.o"
for p in "$ROOT"/pspline/czspline_*.F90; do
    f=$(basename "$p")
    [ "$f" = "czspline_pointer_types.F90" ] && continue
    echo "  FC $f"
    $FC $FFLAGS -I"$BUILD" -I"$ROOT/include" -c "$p" -o "$BUILD/${f%.F90}.o"
done

echo "=== 5. archive ==="
rm -f libpspline.a
ar rcs libpspline.a *.o
ranlib libpspline.a

echo "=== 6. install ==="
cp libpspline.a "$INSTALL/lib/"
cp *.mod "$INSTALL/include/"

echo "=== DONE ==="
ls -la "$INSTALL/lib" "$INSTALL/include"
