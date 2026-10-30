# Testing PB3D

Since version 2.48, PB3D has an automated test suite driven by **CTest**, with
unit tests written in the [test-drive](https://github.com/fortran-lang/test-drive)
framework (vendored under `tests/vendor/test-drive/`, so no extra dependency
is needed).

## Quick start (no dependency stack needed)

The unit tests build against a dependency-light subset of the code
(`pb3d_core`: `num_vars`, `str_utilities`, `messages`, `files_utilities`,
`dtorh`, the light type modules `var_1D_vars`/`grid_vars`/`X_vars`,
`PB3D_utilities`, and - when LAPACK is available - `num_utilities`,
`X_utilities` and `num_ops`; the PSPLINE dependency of `num_utilities` was
split off into `spline_utilities` and the HDF5 dependency of
`PB3D_utilities` into `var_1D_vars`), so they compile with nothing but a
Fortran compiler, CMake and optionally LAPACK:

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
| Full-stack | `fullstack` | full executable build | test-drive suites against the complete `pb3d_modules` (currently: the vacuum module); the distributed suites also run under `mpirun -n 2` and `-n 4` |
| Physics regression | `regression` | full build + `-DPB3D_FIXTURE_DIR=<dir>` | end-to-end eigenvalue anchors on local equilibrium fixtures |

Select layers with labels: `ctest -L unit`, a single suite with
`ctest -R unit_dtorh`, or run the test binary directly
(`./tests/pb3d_unit_tests [suite [test]]`).

The full-stack suites are **rank-agnostic**: they assert on globally
gathered quantities (`tests/fullstack/fullstack_utils.f90` provides a
distributed matrix-vector product through `pdgemv` and gathers through
`vec_dis2loc`/`MPI_Bcast`), so the same assertions run identically on every
process. The suites with distributed linear algebra (`vac_greens`,
`vac_3d`) are registered three times: in MPI singleton mode (where the
whole matrix is one ScaLAPACK block) and as `fullstack_<suite>_np2`/`_np4`
under `mpirun`, where G and H genuinely live in the 2-D block-cyclic
distribution (blocksize 16, BLACS grids 1×2 and 2×2). This covers the
distribution machinery — descriptor setup, `lims_r`/`lims_c` index
bookkeeping, `dgsum2d` gathers, multi-process STRUMPACK/`pdgesv` solves —
that single-process runs bypass. Rank 0 reports to stderr; other ranks
write to `pb3d_fullstack_tests_rank<r>.log`, and the failure count is
combined over all ranks.

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
- **`num_utilities`** — trapezoidal integration vs analytic, finite-
  difference weights vs the classical values, the Björck-Pereyra Vandermonde
  solver, symmetric-storage indexing, sorting, LAPACK determinant/inverse,
  polynomial extrapolation, GCD/LCM/factorial.
- **`num_ops`** — the Householder (orders 1-3) and Zhang zero finders on
  functions with known roots.
- **`X_utilities`** — the pure mode-index logic: local-to-total secondary-
  index translation, the symmetric-storage necessity criterion, the
  contiguous HDF5 ranges of tensorial perturbation variables (pinning the
  worked example from the `get_sec_X_range` documentation) and the mode-
  range trimming.
- **`PB3D_utilities`** — the Richardson/parallel index arithmetic:
  `setup_par_id` checked against its docstring formulas, including the
  partition property (the per-level index sets of an interlaced parallel
  grid tile it exactly), windows (`par_lim`) and per-level memory indices;
  `setup_rich_id`; and the 1-D/n-D storage conversion `conv_1D2ND`.
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
- **`splines` (full-stack)** — the spline interpolation wrapper
  (`spline_utilities::spline`, backed by PSPLINE/EZspline; converted from
  the interactive legacy check): exact reproduction of polynomials in the
  interpolation space (linear/order 1, cubic/order 3 with prescribed
  endpoint derivatives, all derivatives 0-3), the quadratic-Taylor
  extrapolation convention, convergence at the expected rates on
  \(\sin 2\pi x\) for all three orders, periodic boundary conditions, and
  the refusal to extrapolate when not allowed.
- **`calc_int_vol` (full-stack)** — the volume integral
  (`grid_utilities::calc_int_vol`; converted from the interactive legacy
  check) against the analytic torus integral
  \(\int f J = R_0\pi^2 + i\,2\pi^2/3\) for
  \(f = 1 - r^2 + i\cos\theta\), its second-order convergence, and the
  singleton-dimension convention (a missing angular dimension contributes
  a full turn \(2\pi\)).
- **`read_HEL` (full-stack, needs the fixture)** — golden-file test of the
  HELENA equilibrium parser on the committed cbm18a fixture: grid sizes,
  profile values, MISHKA-normalization factors and flux-surface geometry
  pinned; exact-by-construction relations (\(\iota = 1/q\),
  \(\Phi' = q\,\Psi'\), \(\Psi' = 2\pi\)) and a trapezoidal cross-check of
  \(\Phi = \int q \,\text{d}\Psi\) as implementation-independent guards.
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
compare eigenvalues against recorded anchors. The fixtures are committed
xz-compressed under `tests/fixtures/` (HELENA `*.12.xz`, cbm18a: 26 MB →
5.2 MB; the VMEC twin `*.nc.xz`) and decompressed into the build tree at configure time, so the layer runs from
a fresh clone and in CI. Passing `-DPB3D_FIXTURE_DIR=<dir>` overrides this
with a local fixture directory.

Current anchors (see `tests/regression/`; the values moved in October 2026
with the fixes of issue #4, steps 2 to 4, and are asserted to 1e-3 since):

- `cbm18a_fixed_bnd`: fixed-boundary peeling-ballooning mode of the cbm18a
  extended-pseudo-vacuum HELENA equilibrium, n = 10, ψ_N from 0.45 to 0.75,
  \(\omega^2/\omega_A^2 = -4.3573\times10^{-2}\) (MISHKA normalization).
  Fixture `cbm18a_extended_vac.12` is generated by running the PB3D-patched
  HELENA (writes the extended mapping format including `RAXIS, B0` and the
  full R/Z maps) on the original cbm18a input. Cross-platform, one process:
  macOS (gfortran 16, PETSc 3.25), Ubuntu CI (gfortran 13, PETSc 3.19) and
  Snellius (gfortran 14, PETSc 3.25 with MUMPS) agree to 1e-11 relative;
  with 2 and 4 processes the eigenvalue shifts by 2e-5 relative (issue #4,
  open item).
- `cbm18a_free_bnd`: the same deck **free-boundary** with the vacuum
  response, exercising the whole vacuum chain `store_vac_HEL` → `calc_GH` →
  `calc_vac_res` (exterior BEM solve) → `set_BC`. Style 2 (variational,
  Hermitian imposition) \(-4.5019\times10^{-2}\) is the validated value;
  style 4 (explicit natural-boundary-condition row) \(-4.4712\times10^{-2}\)
  is anchored only so that a change shows, it is not a cross-check of style 2.
  The test also asserts free ≤ fixed and that the stored response matrix is
  Hermitian and positive semi-definite (`check_vac_res.py`).
- `cbm18a_shear_HEL`: the HELENA shear S on the whole equilibrium grid
  against the direct formula from the metric derivatives
  (`check_shear_HEL.py`), the check that pins the spline boundary condition
  of step 2.
- `circ_eps005_free_bnd`: the external kink of the circular, flat-current,
  zero-pressure HELENA fixture `circ_eps005_q104` (ε = 0.05, q_a = 1.038),
  n = 10, m = 11: fixed boundary stable, free boundary
  \(-4.3034\times10^{-1}\) against MISHKA's 0.42679 within 2 %
  (`tests/reference/mishka_circular.txt`).
- `circ_eps033_vmec_twin`: the circular ε = 0.33 HELENA fixture and its VMEC
  twin `wout_circ_eps033_q104.nc` (made from the HELENA file by
  `tests/reference/make_vmec_twin.sh`: PB3D `--export_HEL`, then VMEC++) give
  the same fixed-boundary eigenvalue at n = 10 within 2.5 % (measured 1.6 %,
  the boundary fit and the radial coordinate differing between the codes);
  the VMEC value \(+7.9716\times10^{-2}\) is anchored at 1e-3 as the
  layer's only VMEC number. The twin is then run free-boundary through the
  3-D vacuum (`alpha_style 2`, 64 field lines): \(-5.2743\times10^{-1}\),
  anchored at 1e-3, below the fixed value and within 25 % of the HELENA
  style-2 result \(-4.4383\times10^{-1}\) (the first-order BEM error at this
  resolution), which pins the orientation of the normal that
  `store_vac_VMEC` hands to the kernels (issue #4, step 6e).

**The macOS crashes of July 2026** (singleton runs dying with SIGTRAP or a
failed `H5Fcreate` on about every other launch) were a heap overrun in
PSPLINE: `r8genxpkg.f` writes `xpkg(4,4)` for three unevenly spaced points,
8 bytes past the array, which glibc hides in malloc padding and the macOS
allocator traps at a later allocation. `Libraries/build_pspline.sh` now
applies the one-line guard `Libraries/pspline-r8genxpkg-nx3.patch`; the
regression runners launch PB3D once and treat any crash as a failure
(issue #4, step 1).

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
- Multi-process *end-to-end* regression runs: the full-stack suites run at 2
  and 4 processes, the regression decks are single-process in CI (macOS has
  an `mpirun` launcher crash in Homebrew OpenMPI 5's prte). On Snellius the
  decks run at 1, 2 and 4 processes and the 2e-5 shift above is the open
  question (issue #4). The external-code comparison of the free boundary
  exists since step 3 (`circ_eps005_free_bnd` against MISHKA).
- Report the PSPLINE `r8genxpkg` overrun upstream (github.com/ilonster/pspline).
