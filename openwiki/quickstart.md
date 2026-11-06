# PB3D — OpenWiki Quickstart

## What this repository is

**PB3D** ("Peeling Ballooning in 3D") is a scientific/HPC **Fortran** code that computes the
linear ideal-MHD **peeling-ballooning stability** of 3-D (tokamak or stellarator) plasma
equilibria. It was written by Toon Weyens (ITER Organization / TU Eindhoven / Universidad
Carlos III de Madrid, 2012–present); see [`/README.md`](/README.md) for authorship and a
version-by-version changelog (the README is almost entirely changelog — there is no separate
narrative "about" section in-repo, which is why this wiki exists).

PB3D takes a pre-computed 3-D MHD equilibrium — from either **VMEC** (NetCDF output, via the
LIBSTELL library) or **HELENA** (plain-text output) — sets up field-aligned perturbation modes,
optionally computes the vacuum response outside the plasma, and solves a generalized eigenvalue
problem with **SLEPc/PETSc** to find unstable peeling-ballooning normal modes. A companion program,
**POST**, reads PB3D's HDF5 output and produces derived quantities and plots (as HDF5/XDMF for
external visualization tools such as ParaView).

This is not a web application: there is no server, API, database, or UI in the usual sense. The
two "entry points" are the compiled executables `PB3D` and `POST`, run from the command line
against namelist input files and equilibrium files.

## Repository layout

| Path | Contents |
|---|---|
| `PB3D.f90`, `POST.f90` | Top-level programs (main entry points) |
| `Modules/` | ~44 Fortran modules implementing all physics, numerics, I/O, and drivers |
| `Libraries/` | Vendored third-party Fortran sources built as static libs (`dfft.f` FFT, `foul.f90` formatted output) |
| `include/` | C-preprocessor headers used across modules (`PB3D_macros.h` error-check macro, `wrappers.h` compiler-dependent real/imag macros) |
| `Examples/` | Sample namelist input decks and perturbation data files |
| `Test/` | Small standalone Fortran scratch programs used to test isolated features (not part of the CMake/CTest build — see [Testing & Extension](testing-and-extension.md)) |
| `Documentation/` | `spack-setup.md` build guide plus an offline PDF; full Doxygen API docs are hosted externally at <https://pb3d.github.io/Doxygen/html/index.html> |
| `cmake/` | Custom CMake `Find*` modules for PSPLINE and LIBSTELL |
| `CMakeLists.txt`, `Makefile.legacy` | Build system (see [Build & Dependencies](architecture/build-and-dependencies.md)) |
| `ObjectList`, `PB3D.dep` | Leftover object-file/dependency lists from the legacy Makefile build; not used by CMake |

## How to navigate this wiki

- **[Build & Dependencies](architecture/build-and-dependencies.md)** — how to compile PB3D/POST with CMake or Spack, the legacy Makefile, and the heavy scientific-library dependency chain (PETSc, SLEPc, HDF5, NetCDF-Fortran, ScaLAPACK, STRUMPACK-Dense, PSPLINE, LIBSTELL).
- **[Module Map](architecture/module-map.md)** — what every module in `Modules/` does and how they depend on each other.
- **[Simulation Pipeline](workflows/simulation-pipeline.md)** — what actually happens when you run `PB3D` or `POST`: the Richardson-extrapolation loop, equilibrium/perturbation/solution drivers, and command-line usage.
- **[Input & Configuration](workflows/input-and-configuration.md)** — the namelist input files, runtime "style" switches, and example input decks in `Examples/`.
- **[Physics Domains](domain/physics-domains.md)** — the scientific concepts behind each module group: equilibrium (VMEC/HELENA), perturbation modes, vacuum response, eigenvalue solution, Richardson extrapolation.
- **[Testing & Extension](testing-and-extension.md)** — how correctness is checked today (ad-hoc `Test/` programs and `#if ldebug` hooks), versioning conventions, and where to start when extending the code.

## Key facts to keep in mind

- **Two executables, one shared module library**: `PB3D` and `POST` are separate CMake targets that both link against the same `pb3d_modules` object library built from `Modules/*.f90` (see `CMakeLists.txt`).
- **Everything is MPI-parallel** from the start (`start_MPI()` is the first call in both `PB3D.f90` and `POST.f90`); output uses parallel HDF5.
- **VMEC vs. HELENA** equilibrium input is selected at runtime via `eq_style` (1 = VMEC, 2 = HELENA), not at compile time.
- **The vacuum module (`vac_ops.f90`) is verified by the full-stack test suite** (`tests/fullstack`, see `Documentation/testing.md`): axisymmetric kernels/matrices/solve/response, the free-boundary chain end-to-end (`regression_cbm18a_free_bnd`), and the 3-D field-line machinery (cross-checked against the axisymmetric response). Style-1 needs `n_alpha > 1`, and the response has only been validated on axisymmetric boundaries so far — extend the tests when touching the vacuum.
- **Companion tooling lives outside this repo**: cluster run scripts and parameter-scan helpers were moved to the separate `PB3D_tools` repository (see README changelog, version 2.34); this repo's `Examples/` only holds sample input files.
- License: GNU GPLv3 (see [`/LICENSE`](/LICENSE)).
