# Build System & Dependencies

## Two build systems, one preferred

As of **version 2.47** (commit `5c2fe76`, "Modernized build system with CMake and Spack
support"), PB3D builds with **CMake**. The previous hand-written `Makefile` was renamed to
[`/Makefile.legacy`](/Makefile.legacy) and kept only for reference — it hardcodes
machine-specific library paths (e.g. `LIBSTELL_DIR=/opt/stellinstall/bin # 1. XPS 9360`,
commented-out alternates for an "ITER" cluster) and should not be treated as the source of truth
for how to build the project today.

Prefer:
- [`/CMakeLists.txt`](/CMakeLists.txt) for the actual build graph, or
- [`/Documentation/spack-setup.md`](/Documentation/spack-setup.md) for a full reproducible
  environment setup (Spack + manually built PSPLINE/LIBSTELL/STRUMPACK-Dense) on macOS or Linux.

## CMake build graph (`CMakeLists.txt`)

- Requires CMake ≥ 3.18, out-of-source builds only (in-source builds fail fast with a
  `FATAL_ERROR`).
- Detects GNU vs. Intel Fortran and sets compiler flags per build type (`Debug`, `Release`,
  `RelWithDebInfo`); default build type is `Release` (`-O3`).
- Two build-time options:
  - `PB3D_ENABLE_DEBUG` → defines the `ldebug` preprocessor macro, which gates a large amount of
    debug/plotting/testing code throughout `Modules/*.f90` (search for `#if ldebug`).
  - `PB3D_ENABLE_INFINIBAND` → defines `lIB`.
  - The compiler family is also exposed to Fortran code as `lwith_intel` / `lwith_gnu` (used by
    [`/include/wrappers.h`](/include/wrappers.h) to pick `real`/`aimag` vs. `realpart`/`imagpart`
    intrinsics).
- Two internal static libraries are built from vendored sources in `Libraries/`:
  - `dfftpack` from `Libraries/dfft.f` (FFT routines).
  - `foul` from `Libraries/foul.f90` (formatted terminal output).
- All physics/numerics modules are compiled once into an **object library** `pb3d_modules`
  (`Modules/*.f90`, listed in dependency order — see
  [Module Map](module-map.md)) and then linked into **both** executables:
  - `PB3D` ← `PB3D.f90` + `pb3d_modules`
  - `POST` ← `POST.f90` + `pb3d_modules` (POST's own CMake target is defined further down in the
    same file, following the same pattern as `PB3D`)

## External dependencies

CMake locates these via `find_package`/`pkg_check_modules`; several require environment
variables or manual builds (see the Spack guide for exact steps):

| Dependency | How it's found | Why PB3D needs it |
|---|---|---|
| MPI (Fortran/C/CXX) | `find_package(MPI REQUIRED ...)` | All parallelism; `start_MPI()`/`stop_MPI()` wrap the whole program |
| PETSc | `pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)` | Linear algebra backend underlying SLEPc |
| SLEPc | `pkg_check_modules(SLEPC REQUIRED IMPORTED_TARGET SLEPc)` | Solves the generalized eigenvalue problem for stability modes (`SLEPC_ops.f90`) |
| HDF5 (Fortran) | `find_package(HDF5 REQUIRED COMPONENTS Fortran)` | All PB3D/POST output (`HDF5_ops.f90`) plus XDMF metadata for visualization |
| NetCDF-Fortran | `pkg_check_modules(NETCDF_FORTRAN REQUIRED ...)` | Reading VMEC equilibrium files |
| ScaLAPACK | `pkg_check_modules` then a manual `find_library` fallback | Optional; a `WARNING` (not a hard failure) is emitted if missing |
| STRUMPACK-Dense **1.1.1** (old, not the modern STRUMPACK) | Custom `find_library`/`find_path` against `STRUMPACK_DIR` (defaults to `$HOME/Code/STRUMPACK-Dense-1.1.1`) | **Optional** compressed (HSS) solver for the **vacuum** boundary-element system (`vac_ops.f90`); without it the same system is solved with ScaLAPACK LU (`pdgesv`) — CMake sets the `PB3D_WITH_STRUMPACK` preprocessor flag when found |
| PSPLINE | `find_package(PSPLINE REQUIRED)` via [`/cmake/FindPSPLINE.cmake`](/cmake/FindPSPLINE.cmake) | Spline interpolation library (Princeton) used throughout equilibrium/grid interpolation |
| LIBSTELL | `find_package(LIBSTELL REQUIRED)` via [`/cmake/FindLIBSTELL.cmake`](/cmake/FindLIBSTELL.cmake) | Part of the STELLOPT suite; provides the `read_wout_mod` module used by `VMEC_ops.f90` to read VMEC NetCDF output |
| LAPACK / BLAS | `find_package`, optional | General linear algebra |

Note that **STRUMPACK-Dense 1.1.1** is deliberately pinned to an old version — the build uses
`NO_DEFAULT_PATH` specifically to avoid accidentally picking up a modern STRUMPACK install (which
has an incompatible API for this code).

## Setting up dependencies with Spack

[`/Documentation/spack-setup.md`](/Documentation/spack-setup.md) walks through:
1. Installing Spack and bootstrapping a Fortran-capable GCC (Homebrew `gcc` on macOS, since Apple
   Clang has no `gfortran`).
2. Creating and concretizing a Spack environment (`spack.yaml` in the project root — check
   whether it exists before assuming it does; the guide references it as the dependency list to
   edit for cluster-specific MPI/math libraries).
3. Manually building the three dependencies Spack does **not** package: PSPLINE, LIBSTELL (from
   `STELLOPT`, with a hand-written `make_spack.inc`), and STRUMPACK-Dense 1.1.1.

If you are changing build configuration, update both `CMakeLists.txt` and
`Documentation/spack-setup.md` together so they stay consistent — the doc is the only place that
explains *how* to obtain the non-Spack-packaged libraries the CMake `find_*` calls expect.

## Compile-time macros used across the codebase

- [`/include/PB3D_macros.h`](/include/PB3D_macros.h) defines `CHCKERR(s)` (and `CHCKSTT`), the
  standard error-check-and-return-early macro used after nearly every function call that returns
  an `ierr`. Understanding this macro is essential before editing control flow in any module.
- [`/include/wrappers.h`](/include/wrappers.h) abstracts `real`/`aimag` (Intel) vs.
  `realpart`/`imagpart` (GNU) complex-number accessors behind `rp`/`ip`.
