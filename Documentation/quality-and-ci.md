# Code quality & CI strategy

This document is the outcome of a systematic investigation (2026-07) into how
to raise and *hold* the quality of the PB3D code base with automated testing
and continuous integration. It records what exists, what was found, what is
now automated, and the ordered roadmap of what to automate next. Companion
document: `testing.md` (how to run and extend the test suite).

## 1. Guiding principles

1. **Every claim of correctness is a test.** Physics results are pinned by
   implementation-independent references (analytical limits, brute-force
   quadrature, cross-discretization agreement), not by "it looks right".
   The vacuum-module campaign (see `testing.md`) is the template: kernel →
   assembled operator → solve → response → end-to-end eigenvalue, each layer
   asserted in CI.
2. **Everything runs from a fresh clone.** No externally hosted fixtures, no
   manual steps: the regression equilibrium is committed xz-compressed
   (26 MB → 5.2 MB) and decompressed at configure time.
3. **Gates ratchet, they never flap.** Blocking checks are only introduced at
   the current state (warning baseline = today's count; linter = report-only
   until a rule category is clean) and are tightened monotonically.
4. **Fast feedback first.** The dependency-light unit layer runs in seconds
   without the PETSc/SLEPc/HDF5 stack; the expensive layers run in parallel
   CI jobs with cached third-party builds.

## 2. Test architecture (current state)

| Layer | Label | Needs | Runtime | Contents |
|---|---|---|---|---|
| Unit | `unit` | gfortran + LAPACK | < 1 s | `str_utilities`, `files_utilities`, `dtorh` (mpmath references), `num_utilities`, `num_ops` |
| Smoke | `smoke` | full stack | ~1 s | executable usage checks |
| Full-stack | `fullstack` | full stack | ~40 s | vacuum module: kernels vs brute force, Green identities, solves, responses (axisymmetric + 3-D) |
| Regression | `regression` | full stack (+ committed fixture) | ~50 s | end-to-end eigenvalue anchors: cbm18a fixed-boundary and free-boundary (both BC styles, mutual agreement 1e-5) |

Key structural enabler added during this investigation: the PSPLINE-dependent
spline wrappers were split out of `num_utilities` into `spline_utilities`
(34 use-sites updated), and `num_ops`'s HDF5 dependency was confined to
`ldebug`. This lets `num_utilities` (integration, finite-difference weights,
Vandermonde solver, index helpers, sorting, extrapolation, ...) and `num_ops`
(Householder and Zhang zero finders) compile into the light `pb3d_core` and
be unit-tested without any scientific stack. LAPACK became an optional core
dependency (the numerics suites are skipped without it).

## 3. CI pipeline (`.github/workflows/tests.yml`)

All jobs run on every push to `master`/`dev/**`, on pull requests, and on
manual dispatch; superseded runs are cancelled. Third-party builds (PSPLINE,
minimal LIBSTELL) are cached on the hashes of their build scripts.

| Job | Purpose | Blocking? |
|---|---|---|
| `unit` (Release + Debug) | fast unit layer; the Debug leg runs with `-fcheck=all -ffpe-trap=invalid -finit-real=snan` | yes |
| `fullstack` (Release) | complete stack, **all** test layers including both physics regressions; enforces the compiler-warning baseline | yes |
| `fullstack-debug` | Debug + `PB3D_ENABLE_DEBUG` (ldebug code paths compiled in) with full runtime checking; all layers except regression | yes |
| `coverage` | line coverage of `Modules/` from all layers (gcovr); summary + HTML/XML artifact | reporting |
| `lint` | fortitude static analysis, findings histogram + full report artifact | report-only (see §5) |
| `docs` | Doxygen reference build (`Doxyfile`), warning count, HTML artifact | reporting |

Notes:
- Pushes made by the Claude GitHub App do not trigger workflows (GitHub
  anti-loop safeguard); user pushes, merges and manual dispatch do.
- The `PB3D_TEST_MPI_ENV` mechanism in `tests/CMakeLists.txt` handles
  sandboxed-container OpenMPI quirks automatically.

## 4. What the investigation found (inventory)

### 4.1 Unit-testable routines (dependency analysis of `Modules/`)

- **Now tested** (after the spline split): `num_utilities`
  (GCD/LCM/factorial, trapezoidal integration vs analytic, finite-difference
  weights vs classical values, Björck-Pereyra Vandermonde solver, symmetric-
  storage indexing `c`/`is_sym`, sorting with pivots, LAPACK determinant/
  inverse, polynomial extrapolation) and `num_ops` (Householder orders 1-3
  and Zhang zero finders on known roots).
- **Next candidates, cheap** (class "a/b" — light or minor decoupling):
  - `X_utilities` (`sec_ind_loc2tot`, `get_sec_X_range`, `is_necessary_X`,
    `trim_modes`): pure integer mode-index logic; only needs the light
    `X_vars`/`grid_vars` type modules in the core.
  - `PB3D_utilities::setup_par_id`/`setup_rich_id`: pure Richardson/parallel
    index arithmetic with a ~50-line formula docstring — ideal spec tests;
    blocked only by a module-level `use HDF5_vars` (move `var_1D_type` to a
    light module).
  - `num_utilities` remainder: `con`, `conv_mat`, `calc_mult`,
    `add_arr_mult`, `con2dis`/`dis2con`, `round_with_tol`, `derivs`,
    `shift_F`, `order_per_fun`.
- **Medium** (need `use output_ops` demoted to `#if ldebug`, as done for
  `num_ops`): `eq_utilities` (`calc_g`, `calc_inv_met` — which has a genuine
  X·Y=I self-check under ldebug —, `transf_deriv`), `sol_utilities`
  (`calc_tot_sol_vec`/`calc_loc_sol_vec`), `grid_utilities`
  (`calc_eqd_grid`, `nufft`, `calc_int_vol`, `find_compr_range`).
- **Golden-file parser test**: `HELENA_ops::read_HEL` is deterministic given
  a fixture; the committed cbm18a file enables a parsing test that pins the
  derived quantities (fluxes, safety factor, ellipticity) — full-stack layer
  because of the PSPLINE post-processing.
- **Inherently heavy**: `MPI_utilities` (all collectives), the SLEPc/HDF5
  layers, the drivers — covered by fullstack/regression layers instead.

### 4.2 Legacy test machinery (`--test` flag, `Modules/test.f90`, `Test/`)

- Converted to automated suites (self-contained analytic references):
  `test_splines` → `tests/fullstack/test_splines.f90` (polynomial
  exactness, convergence rates, periodic BCs, extrapolation convention)
  and `test_calc_int_vol` → `tests/fullstack/test_calc_int_vol.f90`
  (analytic torus volume integral `R₀π² + i·2π²/3`, convergence,
  singleton-dimension convention). Still convertible:
  `test_calc_D2_smooth` (Holoborodko derivative vs analytic; ldebug-only
  symbol).
- Already superseded: `test_tor_fun` → `tests/unit/test_dtorh.f90`; the
  vacuum debug toggles → `tests/fullstack/test_vac_*`.
- Not automatable: `test_lock` (timing/concurrency stress), ~25 plot-only
  `debug_*` toggles.
- The capital-`Test/` directory is a scratch graveyard (one program
  self-labeled "NOT WORKING", stale committed binaries, an obsolete PETSc
  API demo). **Recommendation: delete it** (history preserves it); nothing
  is worth salvaging.

### 4.3 Compiler warnings

With `-Wall -Wextra` the full build emits **321** genuine warnings (the
baseline in `.github/warnings-baseline.txt`; cpp quote-noise from
apostrophes in comments and vendored `Libraries/` excluded). The dominant
categories are unused imported symbols, real-equality comparisons and unused
dummy arguments. The `fullstack` CI job fails if the count *rises* and asks
for the baseline to be lowered when it falls — a pure ratchet.

### 4.4 Static analysis (fortitude 0.9)

~15k findings across `Modules/`, dominated by style rules (5.9k trailing
whitespace, 3.6k line lengths, 1.1k missing `end module` names). Correctness-
adjacent categories are small: 140 `use` without `only`, 45
`implicit none` without `external`, 24 obsolescent-feature findings.
`fortitude.toml` excludes the E001 false positives on preprocessor lines and
vendored code. Adoption is staged: the CI job is report-only; the ratchet is
to clean one rule category at a time (start with `OB` obsolescent, 24
findings, then `C003`), then add it to a blocking `--select` gate.

### 4.5 Multi-rank (MPI) coverage gap — **closed for the fullstack layer**

Originally, all CI tests ran as a single MPI process: the test harnesses
assumed local = global (plain `matmul` on distributed arrays), leaving the
block-cyclic distribution paths unexercised. The harnesses have since been
made rank-agnostic (`tests/fullstack/fullstack_utils.f90`: distributed
matrix-vector products through `pdgemv`, global gathers, response broadcast
from the last process), and the `vac_greens`/`vac_3d` suites now also run
under `mpirun -n 2` and `-n 4` in every CI leg — covering BLACS grids 1×2
and 2×2, blocksize-16 block-cyclic G/H, `vec_dis2loc` gathers and
multi-process STRUMPACK and `pdgesv` solves, with results identical to the
single-process runs. Remaining gap: the *regression* decks (full PB3D
executable) are still single-process.

## 5. Roadmap (ordered)

1. ~~**Rank-agnostic fullstack tests + multi-rank CI rows.**~~ **Done** (see
   §4.5): the vacuum test harnesses assert on globally gathered quantities
   and the distributed suites are registered at `-n 1, 2, 4` via `mpirun`
   in CTest. Follow-up: a multi-process regression deck.
2. ~~**Legacy-test conversion**~~ **Done**: `test_splines` and
   `test_calc_int_vol` are fullstack suites (`splines`, `calc_int_vol`);
   `test_calc_D2_smooth` remains (ldebug-only symbol).
3. **Warning burn-down**: fix the ~10 core-module warnings, then chip at the
   321 baseline per module; when a module reaches zero, consider
   `-Werror`-listing it.
4. **Fortitude ratchet**: clean `OB` (24), gate it; then `C003` (45), gate;
   style categories via `fortitude check --fix` in one mechanical commit
   each (S101 trailing whitespace is auto-fixable).
5. **Second-stage core expansion**: `X_utilities` and
   `PB3D_utilities::setup_par_id`/`setup_rich_id` unit suites (move
   `var_1D_type` to a light module); demote `use output_ops` to ldebug in
   `eq_utilities`/`sol_utilities`/`grid_utilities` and add their pure
   routines.
6. **`read_HEL` golden-file test** against the committed cbm18a fixture.
7. **Sanitizer job**: `-fsanitize=address,undefined` on the unit layer
   (no MPI involved) is cheap; full-stack ASAN under MPI is noisy — nightly
   at most.
8. **Docs gate**: baseline the 374 Doxygen warnings and ratchet like the
   compiler warnings; optionally publish the HTML artifact to GitHub Pages.
9. **PETSc 3.25 matrix leg** (scheduled weekly, spack build cache) to guard
   the 3.19/3.25 version guards from both sides.
10. **Delete `Test/`** and the stale committed binaries after sign-off.
11. **Valgrind nightly** on the fixed-boundary regression (the old
    `Test/HDF5` scratch program suggests HDF5+Valgrind history) — catches
    leaks the sanitizers miss in MPI context.

## 6. How to keep it healthy

- New code lands with tests in the appropriate layer (see `testing.md` for
  the how-to); reviewers should treat an eigenvalue-affecting change without
  a regression-anchor update as incomplete.
- When a physics anchor legitimately moves, update the anchor *and its
  provenance comment* (platform, solver versions, tolerance rationale).
- Never raise a baseline (warnings, lint, docs) to make CI pass; baselines
  only go down.
- Prefer implementation-independent references (analytic limits, brute
  force, cross-discretization) over stored outputs of the code itself:
  golden values pin behavior, references pin *correctness*.
