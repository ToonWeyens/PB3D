# Testing PB3D

Since version 2.48, PB3D has an automated test suite driven by **CTest**, with
unit tests written in the [test-drive](https://github.com/fortran-lang/test-drive)
framework (vendored under `tests/vendor/test-drive/`, so no extra dependency
is needed).

## Quick start (no dependency stack needed)

The unit tests build against a dependency-light subset of the code
(`pb3d_core`: `num_vars`, `str_utilities`, `messages`, `files_utilities`,
`dtorh`), so they compile with nothing but a Fortran compiler and CMake:

```bash
cmake -S . -B build-tests -DPB3D_BUILD_EXECUTABLES=OFF
cmake --build build-tests -j
ctest --test-dir build-tests --output-on-failure
```

This is also what the GitHub Actions workflow (`.github/workflows/tests.yml`)
runs on every push.

In a full build (with the PETSc/SLEPc/HDF5/... stack available, see
`spack-setup.md`), the same tests plus executable smoke tests are available
from the regular build directory.

## Test layers

| Layer | Label | Needs | Contents |
|---|---|---|---|
| Unit | `unit` | Fortran compiler only | test-drive suites against `pb3d_core` |
| Smoke | `smoke` | full executable build | `PB3D`/`POST` usage-message checks |
| Physics regression | *(planned)* | full build + equilibrium fixtures | end-to-end eigenvalue comparisons on `Examples/` decks |

Select layers with labels: `ctest -L unit`, a single suite with
`ctest -R unit_dtorh`, or run the test binary directly
(`./tests/pb3d_unit_tests [suite [test]]`).

## What is currently covered

- **`dtorh` (toroidal harmonics)** — the mathematical backbone of the
  axisymmetric vacuum response (`vac_ops::calc_GH_2`):
  - values of \(P_{n-1/2}(z)\), \(Q_{n-1/2}(z)\) against 30-digit mpmath
    references (54 points, including the BEM-relevant near-singular regime
    \(z \to 1^+\)), at 1e-11 relative tolerance;
  - the three-term recurrence in the degree;
  - rejection of invalid arguments (\(z \le 1\)).
- **`str_utilities`** — all public conversion/case/merge routines, pinning
  the exact output formats other modules rely upon.

## Adding a test suite

1. Create `tests/unit/test_<name>.f90` defining a module with a
   `collect_<name>` subroutine returning `unittest_type` entries
   (copy `test_dtorh.f90` as a template).
2. Register the module in `tests/unit/main.f90` (one `use`, one
   `new_testsuite` line).
3. Add the source file to `pb3d_unit_tests` and the suite name to
   `PB3D_UNIT_SUITES` in `tests/CMakeLists.txt`.

If the module under test needs more of PB3D than `pb3d_core` provides, first
check whether its module can be added to `PB3D_CORE_SOURCES` in the top-level
`CMakeLists.txt` (only possible if it does not pull in the heavy external
libraries). Otherwise the test belongs in a full-stack layer (to be created —
link against `pb3d_modules` and gate on `PB3D_BUILD_EXECUTABLES`).

Reference-value provenance: values hard-coded in tests should state their
source in a comment (e.g. mpmath version and call), so they can be
regenerated. The generator scripts do not need to be committed, but the
command should be reproducible from the comment.

## Relation to the older test mechanisms

- `Test/` (capital T) holds historical standalone scratch programs with their
  own makefiles — spikes, not regression tests. New isolated experiments can
  still go there, but anything that should *keep* passing belongs in `tests/`.
- The `--test`/`-t` flag of a debug-built `PB3D` runs interactive developer
  checks (`Modules/test.f90`) and is unaffected by this infrastructure.

## Roadmap

- Full-stack unit tests (grid, equilibrium and vacuum quantities) once the
  dependency stack builds on developer machines again — in particular
  `vac_utilities::calc_GH_int_1/2` (singular Green's-function integrals)
  against brute-force quadrature, as groundwork for completing the vacuum
  module.
- Executable smoke tests are registered but unverified until the full stack
  builds (`tests/CMakeLists.txt`, label `smoke`).
- Physics regression anchors: fixed-boundary eigenvalues for
  `Examples/input_cbm18a` (HELENA) and one VMEC case, compared across
  versions.
