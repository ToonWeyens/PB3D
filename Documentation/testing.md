# Testing PB3D

Since version 2.48, PB3D has an automated test suite driven by **CTest**, with
unit tests written in the [test-drive](https://github.com/fortran-lang/test-drive)
framework (vendored under `tests/vendor/test-drive/`, so no extra dependency
is needed).

## Quick start (no dependency stack needed)

The unit tests build against a dependency-light subset of the code
(`pb3d_core`: `num_vars`, `str_utilities`, `messages`, `files_utilities`,
`dtorh`, the light type modules `flat_vars`/`grid_vars`/`X_vars`,
`PB3D_utilities`, and - when LAPACK is available - `num_utilities`,
`X_utilities` and `num_ops`; the PSPLINE dependency of `num_utilities` was
split off into `spline_utilities` and the HDF5 dependency of
`PB3D_utilities` into `flat_vars`), so they compile with nothing but a
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
| Smoke | `smoke` | full executable build | `PB3D`/`POST` usage and help message with their exit status, also under `mpirun -n 2`; stop on a missing input file |
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
`vac_3d`) are registered three times: in MPI singleton mode (where one
process holds the whole matrix) and as `fullstack_<suite>_np2`/`_np4`
under `mpirun`, where G and H genuinely live in the 2-D block-cyclic
distribution (block size 64, BLACS grids 1×2 and 2×2). This covers the
distribution machinery — descriptor setup, `lims_r`/`lims_c` index
bookkeeping, `dgsum2d` gathers, multi-process ScaLAPACK LU solves —
that single-process runs bypass. Rank 0 reports to stderr; other ranks
write to `pb3d_fullstack_tests_rank<r>.log`, and the failure count is
combined over all ranks.

## What is currently covered

- **`dtorh` (toroidal harmonics)** — the mathematical backbone of the
  axisymmetric vacuum response (`vac_ops::calc_GH_2`):
  - values of \(P_{n-1/2}(z)\), \(Q_{n-1/2}(z)\) against 30-digit mpmath
    references (54 points, including the BEM-relevant near-singular regime
    \(z \to 1^+\)), at 1e-11 relative tolerance, and 32 more at the degrees
    of peeling-ballooning modes, n = 20 to 80 (catalogue item T11 of issue
    #4; measured worst case 1.4e-11, tolerance 5e-11);
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
  worked example from the `get_sec_X_range` documentation) and the check
  that two grids have the same rows of modes.
- **`PB3D_utilities`** — the Richardson/parallel index arithmetic:
  `setup_par_id` checked against its docstring formulas, including the
  partition property (the per-level index sets of an interlaced parallel
  grid tile it exactly), windows (`par_lim`) and per-level memory indices;
  `setup_rich_id`; and the flat/n-D storage conversion `conv_flat2ND`.
- **`vac_kernels` (full-stack)** — the vacuum Green's function interval
  kernels (`vac_utilities::calc_GH_int_1/2`) against
  implementation-independent references: direct toroidal-harmonic
  evaluation for regular G, numerical directional derivatives along the
  source normal for H (validates the analytical `Aij` algebra), and
  brute-force panel quadrature with exact log subtraction for the
  (near-)singular analytical integrals, at two interval sizes. Also pins
  the asymptote \(Q_{n-1/2}(1+x) = -\tfrac{1}{2}\ln(x/32) - b_n\) against
  `dtorh`.
- **`vac_3d` additions (5 November 2026)** — `solid_angle_node` (the exact
  solid angle at a faceted node against the axisymmetric corner and Girard's
  theorem), `H_sing_node` (the cell integral of the H kernel against a
  64-point brute force, to 1e-5), `diag_direct_3d` (the row-sum diagonal
  against the smooth direct one, converging), `point_source_direct_3d` and
  `exact_exterior_2d_options` (the point-source solve of the field-line
  style and the exact toroidal harmonics of the axisymmetric style with the
  direct and the indirect diagonals, `vac_diag_style` 2 and 1, at n = 10
  and n = 30, printed for the benchmark; the solution error must fall from
  301 to 601 points, for the direct diagonals at n = 10 only).
- **`fit_diag_pair` (full-stack, issue 29)** — `fit_diag` on synthetic
  \(H\), \(G\) with known diagonals and two potentials, the first \(10^{-9}\)
  of the second at every row: the exact solve of the two equations per row
  returns the diagonals to rounding; with the two equations of one row made
  dependent, the routine has to stop with an error instead of returning NaN.
- **`vac_greens` (full-stack)** — the assembled axisymmetric \(G\), \(H\)
  matrices on a synthetic circular boundary (built like `store_vac_HEL`'s,
  seam point duplicated) through the two jump relations
  \(H\phi - G\,\text{d}\phi = -4\pi\phi\) for interior-harmonic \(\phi\) and
  \(H\phi = G\,\text{d}\phi\) for exterior-harmonic decaying \(\phi\), with
  \(\text{d}\phi\) the derivative along the stored normal \(-J\nabla\psi\) (one
  convention for both vacuum styles), their first-order convergence with
  resolution, and the
  `solve_Phi_BEM` round trip (Neumann data of a known exterior harmonic
  returns its boundary trace) with ScaLAPACK's LU.
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
- **`modes` (full-stack)** — the mode tables of `X_style` 2
  (`X_ops::setup_modes`) for windows on the equilibrium points that move
  up, down, up and then down (a minimum of the safety factor), down and
  then up, up and back down at one equilibrium point, and up and then down
  on a normal coordinate that decreases: at every point of
  every grid the rows of `sec` hold the window once, each mode in its own
  slot of `m`; every grid has the rows of the equilibrium grid, in the same
  order; and a grid on which the window moves over its full width between
  two points is refused.
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
    primitive is documented in the kernel), and the half-interval integrals
    at the regular end of the same interval against closed forms;
  - `GH_int_1_near`: the near-field branch of the kernel (Gauss–Legendre
    along the cubic interpolant through four points of the line) on an arc
    with the radius of curvature of the W7-X bean tip and the spacing of
    the 121-point field lines, the influence point a quarter interval off
    the arc, against brute force over the exact arc (6e-4 for G, 1.2e-3 for
    H, the error of the cubic interpolant);
  - Green's identities for the assembled matrices: with the row-sum
    construction of the H diagonal, the G diagonal from the linear
    potentials and the inward normal \(J\nabla\psi\), interior harmonics
    satisfy \((H+4\pi I)\phi = G\,\text{d}\phi\) and exterior ones
    \(H\phi = G\,\text{d}\phi\) (opposite roles to style 2!); the residuals
    fall from 1.0e-3 and 5.0e-3 at 41×40 to 1.2e-4 and 1.8e-3 at 81×80;
  - a manufactured exterior potential (point source on the magnetic axis)
    through the whole style-1 chain, on the circular and on a rippled
    torus: the solve returns the boundary trace to discretisation accuracy
    at 81×80 and the error falls with resolution;
  - the vacuum response of the axisymmetric boundary computed with the
    full 3-D machinery against the (independently validated) axisymmetric
    result: the five diagonal ratios at 81×80 are 0.994 (m = 1) to 0.982
    (m = 5) and are held to 0.5 % of their recorded values (0.965 to 0.886
    before the near-field quadrature and the G diagonal of 31 October
    2026);
  - `exact_exterior_2d` (catalogue item T4 of issue #4): the exact toroidal
    harmonics \(\sqrt{\cosh\eta - \cos\theta}\,F^n_{k-1/2}(\cosh\eta)
    \cos k\theta\) of the circular torus, \(F = P\) exterior and \(F = Q\)
    interior, through the axisymmetric (style 2) chain: identity residuals
    1e-3 to 1.2e-2 and solution errors 0.9 to 2.4 % at 301 boundary points
    for n = 2, k = 1 to 5 and n = 10, k = 6 to 15, halving at 601 (first
    order). Style 2 takes the normal data as \(-J\nabla\psi\cdot\nabla\phi\)
    with the inward \(J\nabla\psi\) of the fixture, the opposite sign of
    style 1.

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
  \(\omega^2/\omega_A^2 = -4.3570\times10^{-2}\) (MISHKA normalization).
  Fixture `cbm18a_extended_vac.12` is generated by running the PB3D-patched
  HELENA (writes the extended mapping format including `RAXIS, B0` and the
  full R/Z maps) on the original cbm18a input. Like the circle fixtures, it
  is written with 16 digits since 23 November 2026 (issue #65): the same
  HELENA with the formats of its mapping output changed from `E16.8` to
  `E24.16` (see "Patching HELENA" in the tutorial). With the original
  formats that HELENA reproduces the earlier 8-digit fixtures byte for byte.
  The three anchors moved by 6e-5 or less. Cross-platform (measured on the
  8-digit fixture), one process:
  macOS (gfortran 16, PETSc 3.25), Ubuntu CI (gfortran 13, PETSc 3.19) and
  Snellius (gfortran 14, PETSc 3.25 with MUMPS) agree to 1e-11 relative.
  With 2, 3 and 4 processes the eigenvalue agrees to 1.5e-8 since the owner
  rule in `redistribute_var` (30 October 2026); before, it shifted by 2e-5
  because the overlap of two process ranges was filled from the lower rank's
  ghost region.
- `cbm18a_fixed_bnd_np2`, `_np4` (opt-in, `-DPB3D_TEST_MPI_REGRESSION=ON`):
  the same deck with one process and under `mpirun -n 2` / `-n 4`, eigenvalues
  required to agree to 1e-6. Needs a distributed LU in PETSc (MUMPS or
  SuperLU_DIST), which the Snellius build has and the CI and macOS builds do
  not; the 1e-3 anchor rows cannot see a decomposition-dependent shift.
- `cbm18a_free_bnd`: the same deck **free-boundary** with the vacuum
  response, exercising the whole vacuum chain `store_vac_HEL` → `calc_GH` →
  `calc_vac_res` (exterior BEM solve) → `set_BC`. The free boundary is
  `BC_style(2) = 2`: the vacuum energy is added to the edge block and the
  matrices stay Hermitian. Anchor \(-4.5049\times10^{-2}\).
  The test also asserts free ≤ fixed and that the stored response matrix is
  Hermitian and positive semi-definite (`check_vac_res.py`).
- A cross-check that is no longer in the code. Until November 2026 PB3D had
  a second free-boundary treatment, `BC_style(2) = 4`: it replaced the last
  row of the matrices by the natural boundary condition, written with a left
  finite difference, which made the matrices non-Hermitian. With the weights
  of `calc_end_weights` it converged to the value of style 2, with an error
  proportional to 1/`n_r_sol`: on the circular \(\epsilon = 0.33\) case
  \(-0.43266\), \(-0.43531\), \(-0.43661\), \(-0.43725\) at 200, 400, 800
  and 1600 radial points, against \(-0.437883\) for style 2 from 400 points
  on; on this deck \(-4.5088\times10^{-2}\) at 200 points, 0.08 % from style
  2. This is independent evidence that the last row of style 2 tends to the
  natural boundary condition as the step size goes to zero. At a finite step
  size that row also holds the contribution of the last point to the Euler
  equation, which is smaller by a factor of the order of the step size; as
  equations, the last rows of the two styles differed by this term only. With
  half a weight at the last point instead, style 4 went to \(-0.469\) on the
  circular case, 7 % off. The code of style 4 was removed in pull request #38
  (issue #36).
- `cbm18a_shear_HEL`: the HELENA shear S on the whole equilibrium grid
  against the direct formula from the metric derivatives
  (`check_shear_HEL.py`), the check that pins the spline boundary condition
  of step 2.
- `circ_eps005_free_bnd`: the external kink of the circular, flat-current,
  zero-pressure HELENA fixture `circ_eps005_q104` (ε = 0.05, q_a = 1.038),
  n = 10, m = 11: fixed boundary stable, free boundary
  \(-4.2168\times10^{-1}\) against MISHKA's 0.42679 within 2 % (measured
  1.2 %; at n = 15 and 20 the gap is 0.55 and 0.26 %, falling as 1/n²;
  `tests/reference/mishka_circular.txt`), and the same deck at 400 radial
  points within 5e-5 of the 200-point value. This guards the weights of the
  last points of the free edge in the normal integral
  (`num_utilities::calc_end_weights`, 17 November 2026): the move is 1.3e-5
  with them, 1.4e-4 with the trapezoidal weight h/2 of 31 October 2026 and
  1.2e-2 with the full weight before that (when the anchor was
  \(-4.3034\times10^{-1}\), 2.3 % from the converged value). The anchor is
  \(-4.2552\times10^{-1}\) now.
- `circ_eps033_vmec_twin`: the circular ε = 0.33 HELENA fixture and its VMEC
  twin `wout_circ_eps033_q104.nc` (made from the HELENA file by
  `tests/reference/make_vmec_twin.sh`: PB3D `--export_HEL`, then VMEC++) give
  the same fixed-boundary eigenvalue at n = 10 within 1e-4 (measured 2.4e-6
  once B_0 is the toroidal average of |B| on the axis, 1.6 % with the old
  B_0, and 1.1e-6 since `read_VMEC` places VMEC's half-mesh values at VMEC's
  positions, issue #16); the VMEC value \(+7.8431\times10^{-2}\) is anchored at 1e-3 as the
  layer's only VMEC number. The twin is then run free-boundary through the
  3-D vacuum (`alpha_style 2`, 64 field lines): \(-4.3914\times10^{-1}\),
  anchored at 1e-3, below the fixed value and within 3 % of the HELENA
  style-2 result \(-4.3380\times10^{-1}\) (measured 1.2 %, the end-correction
  error of the field-line BEM at 64 lines; both values with the trapezoidal
  weight of the free edge of 31 October 2026), which pins the orientation of the normal that
  `store_vac_VMEC` hands to the kernels (issue #4, steps 6e and 6f).
  Between the fixed VMEC run and the free one, `check_vmec_read.py` compares
  what PB3D stored in `PB3D_out.h5` with the wout file (issue #16). It checks
  the file's mode lists against VMEC's conventions and PB3D's mode lists and
  sizes exactly, R, Z and B_s to 1e-13, and the half-mesh coefficients (λ, the
  Jacobian, B_u, B_v) against a cubic spline from VMEC's positions to 1e-12,
  away from the two surfaces at each end. With the half-mesh positions that
  `read_VMEC` used before 17 November 2026 it fails by 5.4e-6 in λ and 1.1e-6
  in the Jacobian. It needs scipy as well as h5py, and is skipped with a note
  without it. With the same python, the test then runs a copy of the twin
  with the magnetic field reversed (issue #19; `reverse_field.py` changes the
  sign of the fluxes and of the components of B and J). In the copy the
  rotational transform and the toroidal flux have opposite signs, so that
  the poloidal flux decreases outwards and the Jacobian of the flux
  coordinates is negative. The copy must give the twin's fixed-boundary
  eigenvalue to the last digit, with the poloidal and with the toroidal flux
  as normal coordinate, and `check_vmec_read.py` must find it stored as it
  is. The build of 2.55 stops on the copy with exit status 3.
- `reject_BC_style`: PB3D on the cbm18a free-boundary deck with
  `BC_style = 1 4`, a value that the input check refuses (issue #40). The
  test requires an exit status that is not zero, the message of the input
  check in the output, which is redirected to a file, and no eigenvalue
  file. `reject_BC_style_np2` runs the same under `mpirun -n 2` and checks
  that the launcher returns the status. PB3D stops before any linear
  algebra, so this row needs no distributed LU.
- `smoke_PB3D_usage` and `smoke_POST_usage` (label `smoke`, no fixture):
  the program without arguments and with `--help`. The test requires the
  usage message and the exit status of each call. `smoke_PB3D_usage_np2`
  and `smoke_POST_usage_np2` run the same under `mpirun -n 2`, where all
  processes have to stop: only rank 0 reads the command line, and a run in
  which the other process does not learn that it has to stop never ends
  (issue #40); the TIMEOUT of 60 s then fails the test.
- `smoke_PB3D_missing_input` and `smoke_POST_missing_input` (label
  `smoke`, no fixture): the program with an input file that does not exist.
  There is no default input (issue #40). The test requires an exit status
  that is not zero, the name of the file and the link to the tutorial in
  the output, the program's own name in the line that names the calling
  routine, and no file `fort.50`.

**Exit status.** `PB3D` and `POST` return exit status 0 after a complete
run and after `--help`, and 2 after the usage message, that is after a call
with too few arguments. After an error they return
the error code that they print ("ERROR CODE ...") if it lies between 1 and
125, and 1 for every other code, such as the negative codes of HDF5 and
NetCDF (`get_exit_status` in `MPI_ops`). The regression runners test this
status first and look for the output files second.

**The macOS crashes of July 2026** (singleton runs dying with SIGTRAP or a
failed `H5Fcreate` on about every other launch) were a heap overrun in
PSPLINE: `r8genxpkg.f` of the 2011 copy PB3D used then writes `xpkg(4,4)`
for three unevenly spaced points, 8 bytes past the array, which glibc hides
in malloc padding and the macOS allocator traps at a later allocation. Up to
2.53 `Libraries/build_pspline.sh` applied a one-line guard as a patch (issue
#4, step 1); since 2.54 PB3D uses PSPLINE v2.5.0, which has the guard
itself (issue #58). The regression runners launch PB3D once and treat any
crash as a failure.

## Python helpers

The HDF5 checks of the regression layer (`tests/regression/check_*.py`), the
reference derivations in `tests/reference/`, the figure scripts and the
Snellius and TERPSICHORE scripts need numpy, h5py, scipy, sympy, matplotlib
and netCDF4. `pyproject.toml` and `uv.lock` describe that environment for
[uv](https://docs.astral.sh/uv/); `.venv/` itself is not committed:

```bash
uv sync                                   # .venv with everything except vmecpp
uv sync --group vmec                      # also vmecpp, for tests/reference/vmec_twin.py
uv run python tests/reference/cyl_kink.py
```

CMake finds a Python by itself. To make the regression layer use this one,
configure with `-DPython3_EXECUTABLE=$PWD/.venv/bin/python`. CI does not use
uv: its jobs install the pinned `h5py` with pip, and the full-stack job also
`scipy` for `check_vmec_read.py` (`.github/workflows/tests.yml`).

The default group `dev` also installs `fortls`, the Fortran language server of
VS Code's Modern Fortran extension (point its setting `fortran.fortls.path` at
the absolute path of `.venv/bin/fortls`), and `fortitude` at the version of
the CI lint job, so the blocking gate can be run before a push:

```bash
uv run fortitude check --select=OB Modules/ tests/unit tests/fullstack PB3D.f90 POST.f90
```

## Relation to the older test mechanisms

- Anything that should *keep* passing belongs in `tests/`.
- The `--test`/`-t` flag of a debug-built `PB3D` runs interactive developer
  checks (`Modules/test.f90`) and is unaffected by this infrastructure.

## Roadmap

- Extend the full-stack vacuum coverage to the response matrix itself:
  `calc_vac_res` on the circular boundary against the analytical
  large-aspect-ratio (cylinder) vacuum response, validating sign and
  magnitude of `vac%res` as it enters the SLEPc boundary condition
  (`set_BC_2`).
- Full-stack tests for grid and equilibrium quantities.
- Multi-process *end-to-end* regression runs in CI: the full-stack suites
  run at 2 and 4 processes everywhere, the regression decks are
  single-process in CI (its PETSc has no distributed LU; macOS has an
  `mpirun` launcher crash in Homebrew OpenMPI 5's prte). The opt-in
  `cbm18a_fixed_bnd_np2`/`_np4` rows run on Snellius. The external-code
  comparison of the free boundary exists since step 3 (`circ_eps005_free_bnd`
  against MISHKA).
