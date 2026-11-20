---
type: "Guide"
title: "Testing & Extension"
description: "How PB3D is tested: the CTest layers (unit, smoke, full-stack, physics regression) and how to run them, the committed fixtures, the seven CI jobs and their gates, the debug hooks and the old scratch programs, exit status, the versioning rules (prog_version, min_PB3D_version), and where to start when extending the code."
tags: [testing, ctest, ci, regression, test-drive, versioning, debugging]
verified:
  - by: openwiki/0.9.0
    at: 2026-11-20T18:32:33.626Z
sources:
  - id: openwiki-source-80bf3c2c3aef51f0e51d3b0c
    resource: repo://.github/apt-mirror.sh
  - id: openwiki-source-d44494ef3e497fea81240ef8
    resource: repo://CMakeLists.txt
  - id: openwiki-source-c6d176c1d717145bec6bbe53
    resource: repo://Doxyfile
  - id: openwiki-source-2b951757a1f7314b5a33d0bf
    resource: repo://Modules/MPI_ops.f90
  - id: openwiki-source-0e5296407709fa9254158dd8
    resource: repo://Modules/num_vars.f90
  - id: openwiki-source-fa68239bf614d837d7e5522c
    resource: repo://tests/CMakeLists.txt
  - id: openwiki-source-865160cfc030414588d92099
    resource: repo://tests/regression/post_checks.sh
  - id: openwiki-source-8d6a34c174c9363d878d9c95
    resource: repo://tests/unit/test_num_utilities.f90
  - id: openwiki-source-6572c54ba134d228e86232d8
    resource: repo://tests/unit/test_PB3D_utilities.f90
generated: { by: "claude-code", at: "2026-11-20T18:32:33.626Z" }
---

# Testing & Extension

<!-- openwiki: broken internal link [/tests/] link "/tests/" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/.github/workflows/tests.yml] link "/.github/workflows/tests.yml" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/Documentation/testing.md] link "/Documentation/testing.md" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/Documentation/quality-and-ci.md] link "/Documentation/quality-and-ci.md" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
PB3D has a CTest suite in [`/tests/`](/tests/), with unit tests in the vendored fortran-lang `test-drive` framework (`tests/vendor/`), and a GitHub Actions pipeline in [`/.github/workflows/tests.yml`](/.github/workflows/tests.yml). [`/Documentation/testing.md`](/Documentation/testing.md) explains how to run and extend it; [`/Documentation/quality-and-ci.md`](/Documentation/quality-and-ci.md) describes the CI strategy and the quality gates.

## Test layers

| Layer | CTest label | Needs | Contents |
|---|---|---|---|
| Unit | `unit` | a Fortran compiler (LAPACK optional) | suites against the core library `pb3d_core`: `dtorh`, `str_utilities`, `files_utilities`, `num_utilities` (including the end weights), `num_ops`, `X_utilities` (mode-index arithmetic) and `PB3D_utilities` (Richardson index arithmetic) |
| Smoke | `smoke` | the full build | the usage messages of `PB3D` and `POST` |
| Full-stack | `fullstack` | the full build | suites against all of `pb3d_modules`: the vacuum kernels, Green identities, exterior solves and responses of both styles, the splines, the volume integral and the HELENA parser on the committed fixture; the suites with distributed linear algebra also run under `mpirun -n 2` and `-n 4` |
| Physics regression | `regression` | the full build and the fixtures | complete `PB3D` runs whose eigenvalues are compared with recorded anchors |

The full-stack suites assert on quantities gathered over all ranks, so the same assertions hold on every process. [Vacuum Response](domain/vacuum-response.md#tests) lists the vacuum tests one by one.

### Running them

Without the dependency stack, build only the core library and its unit tests:

```bash
cmake -S . -B build-tests -DPB3D_BUILD_EXECUTABLES=OFF
cmake --build build-tests -j
ctest --test-dir build-tests --output-on-failure
```

In a full build, `ctest` runs every layer that the build supports. Select a layer with `ctest -L unit` (or `fullstack`, `regression`), one suite with `ctest -R unit_dtorh`, or exclude the MPI rows with `ctest -LE mpi`.

### Regression decks and fixtures

The equilibrium fixtures are committed xz-compressed in `tests/fixtures/` (the cbm18a HELENA case, two circular HELENA cases at ε = 0.05 and 0.33, and the VMEC twin of the ε = 0.33 case). CMake decompresses them into the build tree at configure time, so the regression layer runs anywhere `xz` is installed; `-DPB3D_FIXTURE_DIR=<dir>` points to other copies. The decks and scripts are in `tests/regression/`:

- `cbm18a_fixed_bnd`, `cbm18a_free_bnd` (with the response-matrix check `check_vac_res.py`) and `cbm18a_shear_HEL` (`check_shear_HEL.py`);
- `circ_eps005_free_bnd`, the circular external kink, which also requires the 200- and 400-point eigenvalues to agree;
- `circ_eps033_vmec_twin`, the HELENA case against its VMEC twin, fixed and free through the 3-D vacuum, with `check_vmec_read.py` comparing what PB3D stored with the `wout` file;
- `reject_BC_style` (also under `mpirun -n 2`), which requires a deck with a removed boundary style to be refused with a status that is not zero;
- `post_checks`, the first test that runs POST: the mode-range input checks, and POST on a run with three solutions;
- with `-DPB3D_TEST_MPI_REGRESSION=ON`, `cbm18a_fixed_bnd_np2` and `_np4`, which require the 2- and 4-process eigenvalues to agree with the 1-process one to 1e-6 (this needs a distributed LU in PETSc, such as MUMPS).

A regression script fails on a crash, on a status that is not zero, and on the word FAIL in its output. The Python checks need `h5py`, `numpy` and `scipy` (`uv sync`; see [Build & Dependencies](architecture/build-and-dependencies.md#python-and-lint-tooling)).

## CI jobs

All jobs run on every push to `master` and `dev/**`, on pull requests and on manual dispatch:

| Job | What it does | Gate |
|---|---|---|
| `unit` | the unit layer, Release and Debug | blocking |
| `sanitize` | the unit layer under the address and undefined-behaviour sanitizers | blocking |
<!-- openwiki: broken internal link [/.github/warnings-baseline.txt] link "/.github/warnings-baseline.txt" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
| `fullstack` | the full stack in Release, all layers including the regression and MPI rows, and the compiler-warning count against [`/.github/warnings-baseline.txt`](/.github/warnings-baseline.txt) | blocking |
| `fullstack-debug` | Debug with `PB3D_ENABLE_DEBUG`, all layers except regression | blocking |
| `coverage` | line coverage of `Modules/` with gcovr | report |
| `lint` | fortitude; the cleaned categories are blocking, the rest a report | gate and report |
<!-- openwiki: broken internal link [/.github/doxygen-warnings-baseline.txt] link "/.github/doxygen-warnings-baseline.txt" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
| `docs` | the Doxygen build, with the warning count against [`/.github/doxygen-warnings-baseline.txt`](/.github/doxygen-warnings-baseline.txt) | ratchet |

<!-- openwiki: broken internal link [/.github/apt-mirror.sh] link "/.github/apt-mirror.sh" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
On AWS-hosted runners, [`/.github/apt-mirror.sh`](/.github/apt-mirror.sh) switches the apt sources to `archive.ubuntu.com`, because the regional mirror is slow.

## Exit status

`PB3D` and `POST` end with a status that is not zero after an error: `abort_MPI` passes the error code to `MPI_Abort` when it lies between 1 and 125, and 1 for every other code. `sudden_stop` writes out the buffered output first, so that the error message reaches a redirected log. Scripts, batch systems and the regression tests rely on this.

## Debug hooks and old test programs

<!-- openwiki: broken internal link [/Modules/test.f90] link "/Modules/test.f90" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- `generic_tests()` in [`/Modules/test.f90`](/Modules/test.f90) is compiled only with `ldebug` (`-DPB3D_ENABLE_DEBUG=ON`) and runs only with the command-line flag `--test` (`-t`), usually together with `--do_execute_command_line` so that it can make plots. It is an interactive developer tool, not an automated test.
- Many routines have `#if ldebug` blocks controlled by a module-level logical such as `debug_calc_derived_q` in `eq_ops.f90`, `debug_run_driver_X_1`/`_2` in `driver_X.f90`, or `debug_setup_mats` and `debug_set_BC` in `SLEPC_ops.f90`, all `.false.` by default. To debug a routine, find its module's `debug_*` flag and set it, rather than adding print statements.
<!-- openwiki: broken internal link [/Test/] link "/Test/" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- [`/Test/`](/Test/) holds old standalone experiments with their own makefiles (array indexing, HDF5 under Valgrind, MPI I/O, PETSc shell matrices, memory use). They are not part of the CMake build or of CTest.

## Versioning

`prog_version` in `num_vars.f90` (2.50) is the program version. PB3D writes it to its output. POST, and PB3D when it continues from earlier output (`rich_restart_lvl` above 1), stop with a message if that output is older than `min_PB3D_version` (also 2.50). When you change the version:

1. Change `prog_version` in `num_vars.f90`, and `min_PB3D_version` as well when the layout of the output changes.
2. Change `project(PB3D VERSION ...)` in `CMakeLists.txt` and `PROJECT_NUMBER` in `Doxyfile`.
3. Add a `## X.YY:` section at the top of the changelog in `README.md`. The changelog is the project's release notes, and most commits add to it.

## Where to start when extending PB3D

- **A new runtime option:** follow the checklist in [Input & Configuration](workflows/input-and-configuration.md#where-to-look-when-adding-a-new-option).
- **A new quantity derived from the equilibrium:** extend `calc_derived_q` in `eq_ops.f90` and the type in `eq_vars.f90`; add the HDF5 output in `print_output_eq`, and, if it should be plotted, a `_plot` routine and a `plot_*` switch.
- **A new POST diagnostic:** add it to `driver_POST.f90`, behind a new `plot_*` switch of the POST namelist.
- **A change to the vacuum:** read [Vacuum Response](domain/vacuum-response.md) and extend the full-stack tests.
- **A change to mode set-up, derived equilibrium quantities, or the matrix assembly** (`setup_modes`, `calc_derived_q`, `setup_mats`/`set_BC`): read the README changelog for that area first and run the regression layer; most of the past serious bugs were there.
- **A build change:** change `CMakeLists.txt` and `Documentation/spack-setup.md` together (see [Build & Dependencies](architecture/build-and-dependencies.md)).
- **A new test:** add a suite under `tests/unit` or `tests/fullstack` (registration is described in the section "Adding a test suite" of `Documentation/testing.md`), or a deck and script under `tests/regression` with an anchor.
