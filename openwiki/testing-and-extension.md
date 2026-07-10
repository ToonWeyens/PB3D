# Testing & Extension

## What "testing" means in this repository today

Since v2.48 (branch `dev/testing-infra`) there is a **CTest-driven test suite** (`tests/`,
unit tests in the vendored [test-drive](https://github.com/fortran-lang/test-drive) framework,
plus a GitHub Actions workflow) — see [`/Documentation/testing.md`](/Documentation/testing.md)
for how to build and run it, including without the heavy dependency stack
(`-DPB3D_BUILD_EXECUTABLES=OFF`). Currently covered: `dtorh` (toroidal harmonics, validated
against 30-digit mpmath references) and `str_utilities`.

In addition, the two historical ad-hoc mechanisms remain:

### 1. Interactive in-program tests (`Modules/test.f90`)

`generic_tests()` in [`/Modules/test.f90`](/Modules/test.f90) is only compiled when the build
defines `ldebug` (i.e. `-DPB3D_ENABLE_DEBUG=ON` / `PB3D_ENABLE_DEBUG` CMake option) and is only
invoked when `PB3D` is run with the `--test`/`-t` command-line flag (see
`files_ops.f90::init_files`). Its own doc comment recommends also passing
`--do_execute_command_line` so it can generate plots interactively while running. This is a
developer/debugging tool, not an automated regression suite — treat it as such, and don't assume
it exercises correctness of the full pipeline.

Throughout the physics modules, many routines have paired `#if ldebug` debug blocks (e.g.
`debug_calc_derived_q` in `eq_ops.f90`, `debug_run_driver_X_1/2` in `driver_X.f90`,
`debug_setup_mats`/`debug_set_BC` in `SLEPC_ops.f90`) that are toggled by a local `logical`
flag defaulting to `.false.` inside each module. **When debugging a specific routine, search for
its module's `debug_*` flags and flip them rather than adding new print statements** — this is
the established pattern.

### 2. Standalone scratch programs (`Test/`)

[`/Test/`](/Test/) contains small, independent Fortran programs with their **own** makefiles,
compiled and run outside of the main CMake build:
- `array_index/` — tests array-indexing logic in isolation (`run_test_array_index.sh`).
- `Compare_PV_KV/` — a GNUPlot-based comparison of two saved datasets (`test.gnu`, `test.dat`).
- `HDF5/` — checks whether HDF5 usage causes issues under Valgrind (see `readme.txt`).
- `MPI_IO/` — MPI file I/O experiments.
- `shell_mats/` — experiments with PETSc "shell" (matrix-free) matrices for the SLEPc solver,
  with a `readme.txt` full of links to PETSc/SLEPc documentation on shell-matrix preconditioning
  — useful background reading before modifying `SLEPC_ops.f90`'s matrix assembly.
- `test_mem/` — memory-usage measurement experiments.

These are historical scratch/spike code, not a maintained regression suite. If you add a new
isolated experiment, follow this pattern (own directory, own tiny makefile/script) rather than
wiring it into the main build.

## Practical verification for changes

Since there's no automated suite, verify changes the way the existing changelog implies the
author does:
1. Build with `PB3D_ENABLE_DEBUG=ON` and run the relevant `debug_*`-gated plots/checks for the
   module you changed.
2. Run `PB3D`/`POST` end-to-end against one of the `Examples/` input decks
   (`input_cbm18a`, `input_Hmode`, `input_qps`, `input_cdxu`, etc. — each pairs with a specific
   equilibrium type/shape) and compare HDF5/plot output before and after your change.
3. For anything touching mode setup (`X_ops::setup_modes`), equilibrium derived quantities
   (`eq_ops::calc_derived_q`), or the SLEPc matrix assembly (`SLEPC_ops::setup_mats`/`set_BC`),
   read the relevant README changelog entries first — these three areas account for a large
   fraction of historical "important bug" fixes, meaning subtle regressions are easy to
   reintroduce.

## Versioning convention

`prog_version` in `num_vars.f90` (currently `2.47`) is the single source of truth for the
program's version and is checked against `min_PB3D_version` when POST reads a PB3D output file
compatibility check. When bumping the version:
1. Update `prog_version` in `num_vars.f90`.
2. Update the CMake `project(PB3D VERSION ...)` line in `CMakeLists.txt`.
3. Add a new `## X.YY:` entry at the top of `README.md`'s changelog describing the change — this
   is the project's de facto release-notes mechanism and is treated as required practice (every
   historical commit touches `README.md`).

## Where to start when extending PB3D

- **New runtime option**: follow the checklist in
  [Input & Configuration](workflows/input-and-configuration.md#where-to-look-when-adding-a-new-option).
- **New physics quantity derived from equilibrium**: extend `eq_ops.f90::calc_derived_q` and the
  corresponding type in `eq_vars.f90`; add HDF5 output in `eq_ops.f90::print_output_eq` and, if it
  should be visualized, a `_plot` routine plus a `plot_*` flag.
- **New POST-only diagnostic**: add to `driver_POST.f90`, gated by a new `plot_*` namelist flag
  the same way existing diagnostics are.
- **Build/dependency changes**: update `CMakeLists.txt` and
  `Documentation/spack-setup.md` together (see
  [Build & Dependencies](architecture/build-and-dependencies.md)).
