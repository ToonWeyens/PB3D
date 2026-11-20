---
type: "Reference"
title: "Build & Dependencies"
description: "How PB3D and POST are built with CMake: build options, the dependency-light core library, the external libraries (PETSc/SLEPc 3.25 with complex scalars, HDF5, NetCDF-Fortran, ScaLAPACK, optional STRUMPACK-Dense 1.1.1, PSPLINE, minimal LIBSTELL), the scripts that build the libraries Spack does not carry, and the Python and lint tooling."
tags: [build, cmake, dependencies, spack, petsc, slepc, snellius, tooling]
verified:
  - by: openwiki/0.9.0
    at: 2026-11-20T18:32:33.626Z
sources:
  - id: openwiki-source-d44494ef3e497fea81240ef8
    resource: repo://CMakeLists.txt
  - id: openwiki-source-c00e5a43ef412324f20dad00
    resource: repo://fortitude.toml
  - id: openwiki-source-f49115d9b2b34576c4a81974
    resource: repo://Libraries/build_pspline.sh
  - id: openwiki-source-2cf1be2deef2573fd19af4cf
    resource: repo://Modules/vac_ops.f90
  - id: openwiki-source-d716790afdbcef047c07c7a1
    resource: repo://scripts/snellius/build_deps.sh
  - id: openwiki-source-283b096e1a4a818b5ef65b2e
    resource: repo://spack.yaml
  - id: openwiki-source-fa68239bf614d837d7e5522c
    resource: repo://tests/CMakeLists.txt
generated: { by: "claude-code", at: "2026-11-20T18:32:33.626Z" }
---

# Build & Dependencies

<!-- openwiki: broken internal link [/CMakeLists.txt] link "/CMakeLists.txt" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/Makefile.legacy] link "/Makefile.legacy" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/Documentation/spack-setup.md] link "/Documentation/spack-setup.md" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/Documentation/snellius-setup.md] link "/Documentation/snellius-setup.md" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
PB3D builds with CMake ([`/CMakeLists.txt`](/CMakeLists.txt)). The old hand-written makefile is kept as [`/Makefile.legacy`](/Makefile.legacy), together with the object and dependency lists `ObjectList` and `PB3D.dep`; it hard-codes machine-specific library paths and is not used by CMake. For a full environment from scratch, follow [`/Documentation/spack-setup.md`](/Documentation/spack-setup.md); for the national supercomputer Snellius, follow [`/Documentation/snellius-setup.md`](/Documentation/snellius-setup.md).

## CMake build graph

- CMake 3.18 or newer is required, and in-source builds stop with a fatal error. The project version in `project()` is the program version (2.50).
- The default build type is `Release` (`-O3` for gfortran). `Debug` adds `-Wall -Wextra -finit-real=snan -ffpe-trap=invalid -fcheck=all`. For gfortran 14.x the Debug build drops the bounds check (`-fcheck=all,no-bounds`), because those compiler versions miscompile the bounds check of a subscript that calls a function with an array-constructor argument; the reproducer is `tests/reference/gfortran14_bounds_repro.f90`.
- Options:
  - `PB3D_ENABLE_DEBUG` defines the preprocessor macro `ldebug`, which switches on the debug, plotting and test code behind `#if ldebug` throughout `Modules/`.
  - `PB3D_ENABLE_INFINIBAND` defines `lIB`.
  - `PB3D_BUILD_EXECUTABLES` (default `ON`) builds `PB3D` and `POST` and therefore needs the full dependency stack. With `OFF`, only the dependency-light core library and its unit tests are built.
  - `BUILD_TESTING` (from `include(CTest)`, default `ON`) adds `tests/`; see [Testing & Extension](../testing-and-extension.md) for the test options `PB3D_FIXTURE_DIR` and `PB3D_TEST_MPI_REGRESSION`.
<!-- openwiki: broken internal link [/include/wrappers.h] link "/include/wrappers.h" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- The compiler family reaches the Fortran code as `lwith_gnu` or `lwith_intel`; [`/include/wrappers.h`](/include/wrappers.h) uses it to choose the complex-number accessors.

### Targets

| Target | Sources | Purpose |
|---|---|---|
| `dfftpack` (static) | `Libraries/dfft.f` | vendored FFT; compiled with `-w` and `-fcheck=no-bounds`, because its F77 work arrays are declared with length 1 |
| `foul` (static) | `Libraries/foul.f90` | vendored formatted terminal output |
| `pb3d_core` (object) | `num_vars`, `str_utilities`, `messages`, `files_utilities`, `dtorh`, `flat_vars`, `grid_vars`, `X_vars`, `PB3D_utilities`; plus `num_utilities` and `X_utilities` when LAPACK is found, and `num_ops` when also not in debug mode | modules that compile with only a Fortran compiler (and LAPACK); the unit tests link against it |
| `pb3d_modules` (object) | all of `Modules/*.f90` in dependency order (see [Module Map](module-map.md)) | shared by both executables |
| `PB3D`, `POST` | `PB3D.f90`, `POST.f90` + `pb3d_modules` | the two programs; `install` puts them in `bin/` |

The core library writes its `.mod` files to `mod/core` and the full library to `mod/full`, because the core sources are compiled into both.

## External dependencies

| Dependency | How CMake finds it | What PB3D uses it for |
|---|---|---|
| MPI (Fortran, C, C++) | `find_package(MPI)` | all parallelism; also needed by the core library in debug mode, where `messages.f90` monitors memory |
| PETSc with complex scalars | `pkg_check_modules(PETSC)` | matrices and linear algebra under SLEPc; the code uses the PETSc/SLEPc 3.25 API |
| SLEPc | `pkg_check_modules(SLEPC)` | the generalised eigenvalue problem (`SLEPC_ops.f90`) |
| HDF5 with Fortran bindings | `find_package(HDF5)` | all PB3D and POST output (`HDF5_ops.f90`) |
| NetCDF-Fortran | `pkg_check_modules(NETCDF_FORTRAN)` | reading VMEC `wout` files through LIBSTELL |
| ScaLAPACK | pkg-config, then `find_library`; a warning if missing | the LU solve of the vacuum boundary-element system |
| STRUMPACK-Dense 1.1.1 | `find_library`/`find_path` in `STRUMPACK_DIR` only (`NO_DEFAULT_PATH`, default `$HOME/Code/STRUMPACK-Dense-1.1.1`) | optional compressed solver for the vacuum system; sets `PB3D_WITH_STRUMPACK` |
<!-- openwiki: broken internal link [/cmake/FindPSPLINE.cmake] link "/cmake/FindPSPLINE.cmake" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
| PSPLINE | [`/cmake/FindPSPLINE.cmake`](/cmake/FindPSPLINE.cmake) | splines, wrapped by `spline_utilities.f90` |
<!-- openwiki: broken internal link [/cmake/FindLIBSTELL.cmake] link "/cmake/FindLIBSTELL.cmake" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
| LIBSTELL | [`/cmake/FindLIBSTELL.cmake`](/cmake/FindLIBSTELL.cmake) | `read_wout_mod`, which `VMEC_ops.f90` uses to read VMEC output |
| LAPACK, BLAS | `find_package`, optional | dense linear algebra; LAPACK also enables `num_utilities` in the core library |

STRUMPACK-Dense is the old dense-only package with a Fortran interface, not the modern STRUMPACK, whose interface is different; the `NO_DEFAULT_PATH` search keeps CMake from picking up a modern install. Without it, `solve_Phi_BEM` in `vac_ops.f90` solves the vacuum system with ScaLAPACK LU (`pdgesv`). With it, the compressed solution is checked against the uncompressed operator, and the solve falls back to ScaLAPACK LU if the relative residual is too large.

Runs on more than one MPI process need a distributed LU in PETSc, such as MUMPS (this is also what the option `PB3D_TEST_MPI_REGRESSION` requires).

## Building the dependencies

<!-- openwiki: broken internal link [/spack.yaml] link "/spack.yaml" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
[`/spack.yaml`](/spack.yaml) pins the Spack part of the stack: `petsc@3.25 +complex`, `slepc@3.25`, parallel HDF5 with Fortran bindings, `netcdf-fortran`, `netlib-scalapack` and `pkgconf`, with OpenMPI and GCC preferred. Three libraries are not in Spack and have their own recipes:

<!-- openwiki: broken internal link [/Libraries/build_pspline.sh] link "/Libraries/build_pspline.sh" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/Libraries/pspline-r8genxpkg-nx3.patch] link "/Libraries/pspline-r8genxpkg-nx3.patch" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- **PSPLINE**: [`/Libraries/build_pspline.sh`](/Libraries/build_pspline.sh) builds it with gfortran from a PSPLINE checkout, against a dummy NetCDF include. Before compiling, it applies [`/Libraries/pspline-r8genxpkg-nx3.patch`](/Libraries/pspline-r8genxpkg-nx3.patch), a one-line guard against a write past the end of an array in `r8genxpkg` for three unevenly spaced points.
<!-- openwiki: broken internal link [/Libraries/build_libstell_min.sh] link "/Libraries/build_libstell_min.sh" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- **LIBSTELL**: [`/Libraries/build_libstell_min.sh`](/Libraries/build_libstell_min.sh) builds only `read_wout_mod` and what it needs from the STELLOPT sources, in the layout that `FindLIBSTELL.cmake` expects. Section 5 of the Spack guide describes the full STELLOPT build instead.
- **STRUMPACK-Dense 1.1.1**: section 6 of the Spack guide.

<!-- openwiki: broken internal link [/scripts/snellius/build_deps.sh] link "/scripts/snellius/build_deps.sh" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
On Snellius, [`/scripts/snellius/build_deps.sh`](/scripts/snellius/build_deps.sh) builds PETSc 3.25.3 (complex, with MUMPS) and SLEPc 3.25.1 against the system's GCC 14 and OpenMPI modules, `build_pspline_linux.sh` and `build_libstell_min_linux.sh` build the other two, and `build_pb3d.sh` configures the Release build (and the Debug build with `PB3D_ENABLE_DEBUG=ON`) with `PSPLINE_DIR`, `LIBSTELL_DIR` and `PB3D_FIXTURE_DIR`. The same directory holds the Slurm job scripts and the comparison helpers; the setup guide lists them.

When you change the build, change `CMakeLists.txt` and the Spack guide together: the guide is the only place that says how to obtain the libraries the `find_*` calls expect.

## Python and lint tooling

<!-- openwiki: broken internal link [/pyproject.toml] link "/pyproject.toml" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- [`/pyproject.toml`](/pyproject.toml) is a `uv` project, not a package: `uv sync` makes `.venv` with `h5py`, `numpy` and `scipy` (used by the regression checks in `tests/regression/check_*.py`), the script group (`matplotlib`, `netCDF4`, `sympy`) and the dev group (`fortls`, and `fortitude-lint` at the version CI uses). `vmecpp` is in an optional group for the VMEC twin recipe.
<!-- openwiki: broken internal link [/fortitude.toml] link "/fortitude.toml" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- [`/fortitude.toml`](/fortitude.toml) configures the Fortran linter: line length 132, rule `E001` ignored because the linter misreads the preprocessor lines, vendored code (`Libraries/`, `tests/vendor/`) excluded.
<!-- openwiki: broken internal link [/.fortlsrc] link "/.fortlsrc" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- [`/.fortlsrc`](/.fortlsrc) configures the Fortran language server: source folders `Modules`, `Libraries` and `tests`, include folder `include`, and the macros `lwith_gnu` and `PB3D_WITH_STRUMPACK` defined.

## Preprocessor headers

<!-- openwiki: broken internal link [/include/PB3D_macros.h] link "/include/PB3D_macros.h" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- [`/include/PB3D_macros.h`](/include/PB3D_macros.h) defines `CHCKERR` and `CHCKSTT`, the check-and-return macros used after almost every call that returns `ierr`. Read it before you change control flow in any module.
<!-- openwiki: broken internal link [/include/wrappers.h] link "/include/wrappers.h" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- [`/include/wrappers.h`](/include/wrappers.h) defines `rp` and `ip`, the real and imaginary parts, for gfortran and Intel.
