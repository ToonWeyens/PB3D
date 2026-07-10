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
| Full-stack | `fullstack` | full executable build | test-drive suites against the complete `pb3d_modules` (currently: the vacuum module) |
| Physics regression | `regression` | full build + `-DPB3D_FIXTURE_DIR=<dir>` | end-to-end eigenvalue anchors on local equilibrium fixtures |

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
- **`vac_kernels` (full-stack)** — the vacuum Green's function interval
  kernels (`vac_utilities::calc_GH_int_1/2`) against
  implementation-independent references: direct toroidal-harmonic
  evaluation for regular G, numerical directional derivatives along the
  source normal for H (validates the analytical `Aij` algebra), and
  brute-force panel quadrature with exact log subtraction for the
  (near-)singular analytical integrals, at two interval sizes. Also pins
  the asymptote \(Q_{n-1/2}(1+x) = -\tfrac{1}{2}\ln(x/32) - b_n\) against
  `dtorh`.
- **`vac_greens` (full-stack)** — the assembled axisymmetric \(G\), \(H\)
  matrices on a synthetic circular boundary (built like `store_vac_HEL`'s,
  seam point duplicated) through the two jump relations
  \(H\phi = G\,\text{d}\phi\) for interior-harmonic \(\phi\) and
  \((H + 4\pi I)\phi = G\,\text{d}\phi\) for exterior-harmonic decaying
  \(\phi\), their first-order convergence with resolution, and the
  `solve_Phi_BEM` round trip (Neumann data of a known exterior harmonic
  returns its boundary trace), both with STRUMPACK and with the ScaLAPACK
  fallback.
- **`vac_3d` (full-stack)** — the field-line 3-D (style 1) vacuum on an
  analytical circular torus covered by field lines \(\zeta = \alpha +
  q\theta\):
  - the singular interval kernel (`vac_utilities::calc_GH_int_1`) against
    brute-force quadrature of \(-1/d\) over the metric half-cell (the exact
    primitive is documented in the kernel);
  - Green's identities for the assembled matrices: with the row-sum
    construction of the H diagonal and the inward normal \(J\nabla\psi\),
    interior harmonics satisfy \((H+4\pi I)\phi = G\,\text{d}\phi\) and
    exterior ones \(H\phi = G\,\text{d}\phi\) (opposite roles to style 2!),
    converging with resolution;
  - the vacuum response of the axisymmetric boundary computed with the
    full 3-D machinery against the (independently validated) axisymmetric
    result: diagonals agree to 0.2 % (m = 1) through 9 % (m = 5, 12 points
    per poloidal wavelength) at 61×60.

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

## Physics regression layer

Regression tests run the real `PB3D` executable on equilibrium fixtures and
compare eigenvalues against recorded anchors. The fixtures are large binary
files kept *outside* the repository; pass `-DPB3D_FIXTURE_DIR=<dir>` at
configure time to enable the layer (tests are skipped silently otherwise).

Current anchors (see `tests/regression/`):

- `cbm18a_fixed_bnd`: fixed-boundary peeling-ballooning mode of the cbm18a
  extended-pseudo-vacuum HELENA equilibrium, n = 10, pedestal window,
  \(\omega^2/\omega_A^2 = -4.074\times10^{-2}\) (MISHKA normalization),
  asserted to ±10%. Fixture `cbm18a_extended_vac.12` is generated by running
  the PB3D-patched HELENA (writes the extended mapping format including
  `RAXIS, B0` and the full R/Z maps) on the original cbm18a input.
  Cross-platform: reproduced on Linux (PETSc/SLEPc 3.19, gfortran 13) to
  within \(10^{-4}\) relative of the macOS (PETSc 3.25, gfortran 16) anchor.
- `cbm18a_free_bnd`: the same deck but **free-boundary** (`BC_style = 1 4`,
  solution grid up to the plasma edge), exercising the whole vacuum chain:
  `store_vac_HEL` → `calc_GH` → `calc_vac_res` (exterior BEM solve) →
  `set_BC_4`. Anchor \(\omega^2/\omega_A^2 = -4.235\times10^{-2}\), ±10%,
  established 2026-07-09 on Linux; bit-identical between the STRUMPACK and
  ScaLAPACK vacuum solver paths, and ~4% more unstable than the
  fixed-boundary result, as expected from freeing the boundary.

**Known flakiness (macOS)**: with Homebrew OpenMPI 5 and parallel HDF5 1.14,
singleton-MPI PB3D runs intermittently die in the I/O layer (SIGTRAP or a
failed `H5Fcreate`, roughly every other launch); successful runs reproduce
eigenvalues to all digits, and the issue is insensitive to OMPI io/btl/osc
component overrides, malloc debugging, and disappears under a debugger. The
regression runner retries up to 6 times as mitigation. Root-causing (e.g.
comparing against an MPICH-based stack) is open follow-up work; Linux/HPC
stacks are expected to be unaffected.

## Relation to the older test mechanisms

- `Test/` (capital T) holds historical standalone scratch programs with their
  own makefiles — spikes, not regression tests. New isolated experiments can
  still go there, but anything that should *keep* passing belongs in `tests/`.
- The `--test`/`-t` flag of a debug-built `PB3D` runs interactive developer
  checks (`Modules/test.f90`) and is unaffected by this infrastructure.

## Roadmap

- Extend the full-stack vacuum coverage to the response matrix itself:
  `calc_vac_res` on the circular boundary against the analytical
  large-aspect-ratio (cylinder) vacuum response, validating sign and
  magnitude of `vac%res` as it enters the SLEPc boundary condition
  (`set_BC_4`).
- Full-stack tests for grid and equilibrium quantities.
- More regression anchors: a VMEC fixed-boundary case and multi-process runs
  (blocked on an `mpirun` launcher crash in Homebrew OpenMPI 5's prte on
  macOS). The free-boundary anchor exists (`cbm18a_free_bnd`); comparing it
  against an *external* code (e.g. MISHKA with vacuum) or against the
  fixed-boundary run on a differently-extended domain would further harden
  it.
- Root-cause the macOS singleton-MPI I/O flakiness (see above).
