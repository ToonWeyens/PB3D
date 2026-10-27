# PB3D: bug fixes and validation before any 3-D vacuum work

Written 27 October 2026. This plan went through three critique rounds; section 7 records what each round changed.

**Starting point: commit 17b712c.** That is the tip of `origin/dev/ci-quality`, the third of the three stacked draft PRs on github.com/ToonWeyens/PB3D, pushed by a cloud session on 10 July 2026 (message: "tests: unit suites for the mode-index and Richardson index arithmetic"). It is the newest PB3D code anywhere: all of Toon's 2.48 work on `dev/testing-infra`, `dev/petsc-port` and `dev/vacuum` (up to 78fb49e, 8 July) plus the 21 cloud commits (vacuum repairs, BC style 2, 3-D kernel rewrite, test suite, CI). The other branches lack the vacuum work: `master` is 5c2fe76 (version 2.47, January 2026) and `dev/vacuum` stops at 78fb49e. `~/Code/PB3D` has 17b712c checked out.

**Bug list:** `~/Code/MHD`, branch `claude/pb3d-report`, file `pb3d_code_analysis.md` (8 October 2026). It was written against 17b712c, so every file and line number in the report and in this plan refers to that commit.

## 0. Ground rules

1. **Smallest change that fixes the bug.** No refactors, no renames, no new abstractions, no new styles. Every fix is a few lines at the place where the bug is, with a one-line comment that says why. The physics code stays readable as Toon wrote it. Section 6 counts the lines.
2. **One bug, one commit, one test.** The test is written first and must fail on the unfixed code with a specific wrong number. If a bug "found by reading" cannot be made to fail a test, it is not fixed until it can.
3. **Predict before you run.** Each physics step states the number it expects (from the report's runs) before the run. A surprise stops the step.
4. **Nothing moves unless the gate passes.** Infrastructure steps must leave every eigenvalue bit-identical. Physics steps must change them by the predicted amount and nothing else.
5. **Tolerances are honest.** Regression anchors are asserted to 1e-3 relative, not 10 %. Crashed runs fail the test; no retries.
6. **Every step passes CI** (Linux, 1, 2 and 4 processes) before the next starts. The Mac cannot run `mpirun`, so CI is the only multi-process signal.
7. **Out of scope** (section 8): negative-flux VMEC files, toroidal-flux vacuum, family coupling, the radial-convergence cause, Richardson improvements, STRUMPACK compression on 3-D kernels.
8. **Local before external.** Every check that runs on this Mac or on CI comes first: the theory, the logical conditions, MISHKA, COBRA, the two public GPEC vacuum codes, Bempp-cl and the published numbers. Requests to other groups (ELITE, CASTOR3D, CAS3D, TERPSICHORE) are sent only in step 7, as confirmation, once the local catalogue is exhausted and every local gate has passed. Nothing waits on an answer from outside.

## 1. Branch

- New branch `dev/vacuum-fixes` from `origin/dev/ci-quality` at 17b712c. This worktree's branch (`claude/pb3d-vacuum-3d-030b8d`, now at `master` with no commits of its own) is reset to it.
- The cloud session's changes outside `vac_*` are mechanical (the spline routines moved to `spline_utilities`, `#if ldebug` guards, PETSc version guards, style 3 removed, default changed). They are kept. Its two vacuum commits (429f003 exterior operator, d6d49e6 3-D kernel) are physics and get read line by line in step 6 before the 3-D fixes.
- This plan is posted as the tracking issue github.com/ToonWeyens/PB3D/issues/4 (27 October 2026). Every step's commit message references the issue, each step's gate result is reported there as a comment with its numbers, and the final PR closes it. The plan file itself is committed on `dev/vacuum-fixes` as `Documentation/plan_vacuum_bugfix.md` and updated when a step changes it.
- At the end: one PR `dev/vacuum-fixes` to `master`, version 2.49, whose description replaces the claims of PRs 1 to 3. Toon merges.

## 2. Step 0: references, before any code change

Nothing in PB3D changes in this step. It produces the numbers every later step is measured against.

| Reference | How | Used by |
|---|---|---|
| R1 Baseline table | Build Release and Debug (`-fcheck=all`). Run every deck of this plan on 17b712c and record ω² to all digits: cbm18a fixed at ψ_max 0.75 and 1.0; free at 0.75 and 1.0 with styles 2 and 4; vacuum off (scale 0) at 0.75; the cylinder `res` errors; the style 1/2 ratios of `test_vac_3d`. Same on CI. | every step |
| R2 Circular external-kink case | Regenerate the PhD case `Simulations/TEST_VAC/2.25` with the rebuilt HELENA (`~/Code/helena/hel12a`): circular, ε = 0.33, flat current, B = 0.001, XIAB = 2.25, q_edge = 1.038. Modes with m/n just above q_edge have no resonance inside, so the fixed boundary is stable and the free boundary is unstable. Make the VMEC twin with `hel2vmec.sh` (PB3D's HELENA-to-VMEC export). Also a large-aspect-ratio variant, ε = 0.05, for the cylinder formula. | steps 3, 4, 5, 6 |
| R3 MISHKA | Compile MISHKA-1 from `Documents/PHD/MISHKA_simulations/Sources` (recipe in its README). Run R2 for n = 10 to 20 with `RWALL` ≤ 1 (fixed boundary) and `RWALL` = 2, 3, 5, 10, 20 (conformal wall, vacuum on). Extrapolate to no wall. Record both ω²_fixed and ω²_free, because the vacuum effect is the difference. | steps 3, 4, 6 |
| R4 Analytic cylinder | Shafranov's uniform-current external kink with no wall: δW and the growth rate in PB3D's normalization, as a function of nq_a − m, derived with sympy and checked against MISHKA on the ε = 0.05 variant. | step 3 |
| R5 PhD numbers | The 2017 PB3D versus MISHKA (`RWALL` = 3) tables in `TEST_VAC/OLD`: PB3D −0.298 versus MISHKA 0.434 at n = 10, −0.430 versus 0.429 at n = 16, −0.330 versus 0.323 at n = 20 for the 2.25 case; two to four times apart for the 0.77 case. The folder's README says "sometimes positive, sometimes not". These numbers are only history; they show the comparison was never on equal footing (MISHKA had a wall, PB3D had the sign and factor errors). | context |

Gate: R1 reproduces the report's numbers (−4.0736e−2 fixed at 0.75, −4.3467e−2 free style 2 at 0.75, −4.5574e−2 style 4 at 0.75). R3 converges in `RWALL` and its wall sequence follows the cylinder wall factor (T2) on the ε = 0.05 case before MISHKA is trusted. R4 agrees with R3 at ε = 0.05 to the accuracy MISHKA's resolution allows. Step 0 also runs the radial-order diagnosis (O4) and measures every convergence order (L12), so each later comparison knows its own accuracy.

Risk: MISHKA's recipe links `-fdefault-integer-8` against LAPACK. On macOS that needs an ILP64 LAPACK, or dropping the flag if `ppplib` tolerates it. If MISHKA cannot be built in a day, steps 3 and 4 fall back on R4 plus the bounds, and MISHKA moves to CI (Linux) later.

### 2.1 Step 0 results (27 October 2026)

All runs on the macOS machine (gfortran 16.1, PETSc 3.25.3 complex, SLEPc 3.25.1, singleton MPI), code 17b712c unchanged, Release build unless stated. ω² is ω²/ω_A² in the MISHKA normalisation, as PB3D prints it. The Snellius build (Toon: "this should be on Snellius") is reported separately when it is done.

**R1 baseline table.** Regression deck cbm18a, n = 10, 20 modes, 200 radial and 100 parallel points (Simpson 3/8), `tol_SLEPC` 1e−6.

| Run | ω² | Report |
|---|---|---|
| fixed, ψ 0.45 to 0.75 | −4.0735999821e−2 | −4.0736e−2 |
| free style 2, to 0.75 | −4.3467138270e−2 | −4.3467e−2 |
| free style 4, to 0.75 | −4.5574082987e−2 | −4.5574e−2 |
| fixed, 0.45 to 1.0 | −4.2348778332e−2 | anchor −4.2349e−2 |
| free style 2, to 1.0 | −4.2348778421e−2 | |
| free style 4, to 1.0 | −4.2348776924e−2 | |
| vacuum × 0 at 0.75 (diagnostic patch) | −4.3142687497e−2 | −4.3143e−2 |
| vacuum × (−2π) at 0.75 | −4.2150080544e−2 | −4.2150e−2 |

Edge diagnostics from the patched scratch build: |X_edge|/max|X| = 0.100 (× 1), 0.094 (× 0), 0.068 (× −2π); X_e† res X_e = −2.6476e−2 (× 1), +3.9853e−2 (× −2π); max|res − res†|/max|res| = 2.65e−5. The three ψ = 1.0 runs agree to 4e−8, so that deck cannot see the vacuum (report section 9). Gate R1 passed: every number of the report reproduces to five digits.

Crashes: 9 of 15 launches of the regression decks on the unpatched build died (SIGTRAP or segfault) before the eigenvalue was written; the 800-point radial deck died 8 of 8 times. A scratch build of the same source linked against PSPLINE with the `r8genxpkg` guard applied ran 23 of 23 launches clean, with eigenvalues identical for the fixed deck and within 4e−13 relative for the free decks (SLEPc noise). That is the step 1 fix, confirmed in advance; step 1's gate tolerance should read 1e−10 relative, not bit-identical, for free-boundary runs.

Unit-level numbers (`pb3d_fullstack_tests`): cylinder response `res` against the analytic value, m = 1 to 5: 1.0, 1.3, 1.6, 1.8, 1.9 % (report 1.0 to 1.9 %); Green-identity residuals halve from 51 to 101 boundary points (6.7e−2 → 3.5e−2, 3.5e−2 → 1.8e−2, 2.5e−2 → 1.3e−2: order 1); style 1 over style 2 response ratios m = 1 to 5: 1.002, 0.981, 0.957, 0.931, 0.906 (report 1.00 → 0.91). STRUMPACK's HSS solve fell back to ScaLAPACK on the 3-D test (residual 0.15).

**O4 and L12, convergence orders** (fixed boundary, ψ 0.45 to 0.75, unpatched shear).

| Direction | Values | Order |
|---|---|---|
| parallel points 97, 100, 193, 385 (Simpson 3/8) | −4.0824e−2, −4.0736e−2, −4.2082e−2, −4.2800e−2 | differences 1.26e−3 then 7.2e−4: about 0.8 (the end-interval error of the shear bug, report section 6) |
| parallel 97 with `magn_int_style = 1` | −4.0118424e−2 | equals the report's table: that table used the trapezoid rule, the regression deck uses Simpson 3/8 |
| radial points 100, 200, 400, 800, 1600 (left stencil, patched build) | −4.14807e−2, −4.07360e−2, −4.08179e−2, −4.08868e−2, −4.09221e−2 | differences 7.4e−4, −8.2e−5, −6.9e−5, −3.5e−5: not monotonic, at best order 1 |
| radial 800, centred stencil | −4.20465e−2 | 2.8 % from the left stencil (report 2.7 %) |
| HELENA `NRMAP` 201, 401, 801 at 200 radial points | −4.07128e−2, −4.07142e−2, −4.07360e−2 | 0.06 %: the equilibrium grid is not what limits the radial order |

So the slow radial convergence lives in the solution discretisation or in the equilibrium-to-solution interpolation, not in HELENA's grid. Still out of scope (section 8), but every comparison below carries about 1 % of radial uncertainty on PB3D's side.

**R2, the circular external-kink cases.** HELENA `hel12a` (PB3D-patched) on the PhD input `TEST_VAC/2.25/HELENA/fort.10`: ε = 0.33, flat Γ, B = 0.001, XIAB = 2.25, q₀ = 1.0411, q_a = 1.0382. The large-aspect-ratio variant: ε = 0.05, XIAB = 0.3035. A perfectly flat q (q₀ = q_a = 1.0378) is rejected by PB3D ("The code has to be adapted for reversed-shear regions still", `calc_eq_1`, error 3), so the committed variant has Γ = 1 − 0.02ψ (AGA = −0.02): q₀ = 1.0274, q_a = 1.0378. Both equilibria are committed as `tests/fixtures/circ_eps033_q104.12.xz` and `circ_eps005_q104.12.xz` with their HELENA inputs. PB3D on them, 17b712c, deck ψ 0.2 to 1.0, 200 radial, 100 parallel points, 10 modes, `BC_style 1 1` (fixed) and `1 2` (free):

| Case | n | modes | fixed | free (style 2) |
|---|---|---|---|---|
| ε = 0.33 | 10 | 6 to 15 | +7.8446e−2 | +7.8190e−2 |
| | 15 | 11 to 20 | +8.6149e−2 | −1.7221e−1 |
| | 20 | 16 to 25 | +1.2056e−2 | −3.4585e−1 |
| ε = 0.05 | 10 | | +7.7650e−2 | +7.7650e−2 |
| | 15 | | +1.6146e−1 | −6.4661e−1 |
| | 20 | | +4.8685e−2 | −4.0010e−1 |

Fixed boundary is stable in every case, as it must be (no resonance inside). The free-boundary numbers are what the unfixed code gives; at n = 10 the vacuum changes nothing although |X_edge|/max|X| is 0.81 to 0.86 in the free runs. These are step 3's starting values. The VMEC twins were not made: there is no VMEC binary on this machine (VMEC++ exists in `~/Code/desc-vs-vmec`); they are first needed in step 5 and are made there.

**R3, MISHKA.** The source that reads the 2017 HELENA format (`JS0, IAS` on the first line, `RAXIS, B0` after the metric) is `Simulations/MISHKA/Sources` (November 2017, MANZ = 21), not `Documents/PHD/MISHKA_simulations/Sources` (2015, old format, fails with "Bad integer for item 409"). Built with `gfortran -std=legacy -w -O2 -fdefault-integer-8 -fdefault-real-8` and linked against an ILP64 OpenBLAS (`spack install openblas+ilp64 symbol_suffix=none`); with 4-byte integers it segfaults at once. It reads the namelist from stdin and the mapping from `fort.12`, writes `fort.20`; the line `INSTABILITY FOUND : n ... q₀ q_a Re Im` gives Re = γ² = −ω²/ω_A². Validation: the archived PhD run (ε = 0.33, n = 10, `RWALL` = 3) reproduces to all printed digits, 0.433847, with the archived `fort.12`, and gives 0.433764 with the regenerated one. `RWALL` = 1 (fixed boundary) does not converge in the inverse vector iteration (99 iterations), so MISHKA is used for the free-boundary values only.

| Case | n | RWALL 1.05 | 1.1 | 1.2 | 1.5 | 2 | 3 | 10 |
|---|---|---|---|---|---|---|---|---|
| ε = 0.05, flat q | 10 | 0.06443 | 0.33592 | 0.42327 | 0.43652 | 0.43662 | 0.43662 | 0.43663 |
| | 15 | 0.36337 | 0.43842 | 0.45472 | 0.45576 | 0.45576 | 0.45576 | 0.45577 |
| | 20 | 0.32596 | 0.34023 | 0.34223 | 0.34229 | 0.34229 | 0.34229 | 0.34229 |
| ε = 0.05, Γ = 1 − 0.02ψ | 10 | 0.05473 | 0.32609 | 0.41342 | 0.42667 | | | 0.42679 |
| | 15 | 0.35606 | 0.43109 | 0.44739 | 0.44843 | | | 0.44844 |
| | 20 | 0.32160 | 0.33586 | 0.33786 | 0.33792 | | | 0.33792 |
| ε = 0.33 | 10 | 0.10405 | 0.33363 | 0.41715 | 0.43343 | 0.43376 | 0.43376 | 0.43378 |
| | 15 | 0.35788 | 0.42539 | 0.44123 | 0.44281 | 0.44282 | 0.44282 | 0.44283 |
| | 20 | 0.30787 | 0.32065 | 0.32253 | 0.32262 | 0.32262 | 0.32262 | 0.32262 |

Values are −ω²/ω_A² (positive = unstable). Converged in `RWALL` from 1.5 on: a wall at 1.5a or further is no wall for m ≥ 11.

**R4, the cylinder.** Derived with sympy (`tests/reference/cyl_kink.py`): reduced-MHD energy with the exterior Laplace solution written out. Uniform current, zero pressure, uniform density, no wall:

ω²/ω_A² = −(2/q_a²)(m − nq_a)(nq_a − m + 1), unstable for m − 1 < nq_a < m;

with a conformal wall at b: −(m − nq_a)[(1 + λ)(nq_a − m) + 2]/q_a², λ = (1 + (a/b)^{2m})/(1 − (a/b)^{2m}), which closes the window to m − 2/(1 + λ) < nq_a < m. Against MISHKA on the flat ε = 0.05 case without a wall: 0.43660 vs 0.43663 (n = 10), 0.45590 vs 0.45577 (15), 0.34254 vs 0.34229 (20), i.e. 1e−4 to 7e−4. The wall factor (T2): MISHKA over formula between 0.9992 and 1.0011 for every `RWALL` ≥ 1.05 and every n. On ε = 0.33 the ratio is 0.990, 0.9755, 0.965 for n = 10, 15, 20, independent of the wall: the toroidal correction alone. On the Γ = 1 − 0.02ψ variant the ratio is 0.977, 0.984, 0.987, also wall-independent: the current peaking alone. Gates R3 and R4 passed.

**R5.** The PhD MISHKA numbers (0.434, 0.429, 0.323 at n = 10, 16, 20, `RWALL` = 3) are reproduced by the rebuilt code and are within 1 to 4 % of the cylinder formula, so in 2017 MISHKA was right and PB3D's −0.298, −0.430, −0.330 were wrong.

**Reference values for step 3** (free boundary, no wall, −ω²/ω_A²; PB3D must reach them within its own discretisation error, about 1 %, with the fixed boundary staying positive):

| Case | n = 10 | n = 15 | n = 20 |
|---|---|---|---|
| ε = 0.33 (`circ_eps033_q104`) | 0.43378 | 0.44283 | 0.32262 |
| ε = 0.05, Γ = 1 − 0.02ψ (`circ_eps005_q104`) | 0.42679 | 0.44844 | 0.33792 |

MISHKA input for these: the archived `TEST_VAC/OLD/2.25/MISHKA/10/fort.10` with `NTOR` and `RWALL` changed (`IFAST = 1` shifts the m window with nq by itself; `DSURF = 1` is a flat density). Table kept as `tests/reference/mishka_circular.txt`.

**Still open from step 0:** the same table on CI (CI is red for other reasons; step 1), the Snellius build and its 1, 2, 4-process runs (in progress), the VMEC twins (step 5).

## 3. Steps 1 to 8

| Step | Bug (report #) | Change | Physics lines | Gate |
|---|---|---|---|---|
| 1 | Crashes on macOS (1); CI red (12) | PSPLINE guard applied by `Libraries/build_pspline.sh`; retries removed from the regression scripts; CI warning count de-duplicated, `if: always()` on the test step, gcovr flag | 0 | 20 of 20 clean runs of each regression deck; eigenvalues identical to R1 for the fixed deck and within 1e−10 relative for the free decks (SLEPc noise, measured in step 0); CI green |
| 2 | HELENA shear (2) | `bcs(:,1)` to `bcs(:,2)` at `eq_ops.f90:4409`, comment at 4354 corrected | 1 | predicted −4.0118e−2 to −4.3568e−2 at 97 parallel points with `magn_int_style = 1` (the report's deck); the regression deck (`magn_int_style = 2`) gives −4.0824e−2 at 97 points and must reach the same converged value; Δθ-independence; identity test; VMEC twin; POST energy 0.07 % |
| 3a | Vacuum sign (3) | none yet: the test | 0 | ordering test fails on current code: free (−4.3467e−2) lies outside [no-vacuum −4.3143e−2, fixed −4.0736e−2] |
| 3b | Vacuum sign and factor (3) | the energy factor applied once, at the end of `calc_vac_res`, with the derivation in the comment; `res` symmetrised; POST's `if (1.eq.2` guard removed | 3 | ordering holds; circular case: fixed stable, free unstable, marginal point at the right nq_a; growth rates against R3 difference and R4; POST reproduces ω² |
| 4 | Style 4 drops the vacuum (4); default (4) | (0,0) block assembled once before insertion; default back to `[1,2]` | about 6 | style 2 = style 4 at ψ = 0.75 to 1e-5; regression deck moved to 0.75; anchors re-set with 1e-3 |
| 5 | `setup_modes` slots (5) | slot computed from the mode number, same formula as `nm_X`; consistency check of the two tables | about 8 | simulation script: 0 wrong entries in all 8 sequences; HELENA and VMEC twins of R2 and cbm18a agree fixed-boundary |
| 6a | review | read 429f003 and d6d49e6 against the vacuum paper and its known errors | 0 | written verdict per hunk: keep, or revert |
| 6b | 3-D α spacing (9); normal orientation (13) | `dalpha` with `n_alpha` at `vac_ops.f90:839, 1567`; test grid built like PB3D's | 2 | manufactured exterior potential: BEM error falls with resolution on a torus and on a rippled torus; style 1/2 ratio → 1 with measured order |
| 6c | Richardson double count (10) | style 1 recomputes G, H each level (skip the copy) or zeroes before accumulation | 1 to 3 | level-2 `res` equals a fresh run at level-2 resolution |
| 6d | aliasing (11) | warning when `n_alpha` divides the field period | about 5 | input check fires |
| 6e | end to end at ε = 0 | none | 0 | circular VMEC twin through the style-1 vacuum equals the HELENA style-2 result and R3 |
| 6f | continuity in δ | none | 0 | rippled circular VMEC: ω²(δ) − ω²(0) ∝ δ², smooth; ordering holds on W7-X |
| 7 | external confirmation, changelog, PR | requests to ELITE, CASTOR3D and CAS3D groups go out now; 2.49 entry, one line per fix, the cloud's silent default changes documented | 0 | external results agree within their known blind spots; PR open; Toon merges |
| 8 | 3-D equilibrium path | none in PB3D: DESC-written `wout` of a Landreman case; compare ι, S, κ, σ with the exact formulas | 0 | agreement to the resolution of the DESC fit |

### Step 1: stop the crashes, make CI mean something

- **PSPLINE.** `r8genxpkg.f` writes `xpkg(4,4)` for 3 unevenly spaced points, 8 bytes past the array. PB3D hits it in `interp_V_spline` through `spline_real` for mode ranges with exactly 3 points. The report's two-line guard goes into a patch file in `Libraries/` that `build_pspline.sh` applies after cloning, and the bug gets reported upstream. PB3D's own code is not touched; the alternative (route n ≤ 3 around EZspline in `spline_real`, like the 2.46 fix for n = 2) stays in reserve for users with an unpatched PSPLINE.
- **Regression scripts.** Remove the six-attempt loop. A nonzero exit is a failure.
- **CI.** Count unique `file:line` warnings and exclude `tests/`; run the test step with `if: always()`; add `--gcov-ignore-parse-errors=suspicious_hits.warn`; pin gcovr and fortitude.
- **Gate.** 20 consecutive clean runs of both regression decks on the Mac in Release and Debug; eigenvalues bit-identical to R1; all CI jobs green including `_np2` and `_np4`.

**Step 1 results (27 October 2026).** This step's commit. No physics line changed.

- *PSPLINE.* `Libraries/pspline-r8genxpkg-nx3.patch` (the two-line guard, with the explanation in its header) is applied by `Libraries/build_pspline.sh` after the clone; CI copies the patch along with the script and keys its PSPLINE cache on both files. The Mac's `~/Code/pspline/install` was rebuilt the same way. Text for the upstream report is in the step 1 comment on issue 4; posting it on github.com/ilonster/pspline is Toon's call.
- *Regression scripts.* One launch, a crash is a failure; the stale "fixture not committed" header and the flakiness paragraph of `Documentation/testing.md` replaced by the cause.
- *Vendored FFT in Debug.* The free-boundary deck had never run on a Debug build: `Libraries/dfft.f` declares its work arrays with the F77 idiom `DIMENSION WSAVE(1)`, so `-fcheck=all` stops the first FFT call in `store_vac` with a false bounds error (20 of 20 launches). The `dfftpack` target now compiles with `-fcheck=no-bounds` next to its existing `-w` (one option in `CMakeLists.txt`). This was not in the report because its Debug runs were fixed-boundary only.
- *CI.* Warning ratchet counts unique source locations in `Modules/`, `PB3D.f90`, `POST.f90` from gfortran's location headers (the old counter counted every repeated diagnostic and the duplicates of sources compiled into two targets, and could not see `tests/`); baseline 74 (63 maybe-uninitialized, 11 unused dummy arguments) measured on the Linux build of 2058627. Test step runs with `if: always()`. gcovr pinned to 8.8 with `--gcov-ignore-parse-errors=suspicious_hits.warn`; fortitude pinned to 0.9.2 (0.10.0, picked up unpinned on 27 October, adds rule OB091 and fails the OB gate on two labelled `DO` loops in `dtorh.f90`; the code is not changed for a linter).
- *Gate.* Release: 20 of 20 clean launches of each deck, fixed eigenvalue identical to R1 (−4.0735999821004454e−2), free style 4 identical (−4.2348776924299336e−2), free style 2 spread 5e−12 relative (SLEPc). Debug: fixed 20 of 20, −4.0735999821102098e−2 (2e−12 from Release); free 20 of 20, style 4 −4.2348776924365450e−2, style 2 spread 1e−11. CI: the run on this commit is the first one with the new ratchet and the pins; its result is appended below once it has run.
- *Learning for the plan.* The PSPLINE overrun was the whole "MPI flake"; nothing is wrong with OpenMPI 5 or parallel HDF5 on macOS. Debug builds of the vacuum path need the dfftpack exemption, so step 3's Debug sweeps (O5) are possible now. Free-boundary eigenvalues move at 1e−12 between launches, so "bit-identical" gates apply to fixed-boundary runs only.

---

### Step 2: the HELENA shear

- **Change.** One line. `Θ^α` is odd in θ, so the spline on a top-bottom symmetric equilibrium needs the odd condition (zero second derivative), `bcs(:,2)`. The comment on line 4354 says "zero first derivative" for `[2,2]`; it should say second.
- **Test first.** A parity test on cbm18a: S at the outboard midplane is not zero and does not change by more than 1e-4 between 97 and 193 parallel points. On the unfixed code it changes by 3.7 %.
- **Predicted.** −4.0118e−2 becomes −4.3568e−2 at 97 points, −4.3571e−2 at 1537. Step 0 established that the report's parallel-points table used `magn_int_style = 1` (trapezoid); with the regression deck's `magn_int_style = 2` (Simpson 3/8) the unpatched values are −4.0824e−2, −4.2082e−2, −4.2800e−2 at 97, 193, 385 points. The converged value is the gate, not the 97-point number.
- **Independent checks.** (i) The identity J|∇ψ|²S + μ0 J B² σ = K (`test_S_HEL`, already in the code under `ldebug`) to 1e-6. (ii) The VMEC twin of cbm18a computes S from Fourier derivatives and never had the bug: HELENA and VMEC fixed-boundary ω² must agree to the VMEC fit accuracy (expected 1e-3). This check also depends on step 5 for VMEC; if they disagree, step 5 is pulled forward. (iii) POST energy reconstruction: E_pot/E_kin misses ω² by 5.0 % now and by 0.07 % after the fix.
- **Anchors.** Fixed regression re-set to about −4.357e−2, tolerance 1e-3.
- **More checks (section 4).** L8: the same equilibrium through the periodic branch (`ias = 1`) gives the same ω² as the symmetric branch. O1: the JCP 2017 PB3D curve, which predates the bug, is reproduced with the data-file settings. T8: the 1/n extrapolation reaches COBRA again. L3: ω² falls monotonically as `max_r_sol` and `n_mod_X` grow.

**Step 2 results (27 October 2026).** This step's commit. One physics line changed (`eq_ops.f90`: `bcs(:,2)` instead of `bcs(:,1)` in the call of `calc_derived_S_from_deriv_HEL`), three comments corrected ("zero second derivative" for `[2,2]`).

- *Test first.* `tests/regression/cbm18a_shear_HEL.sh` runs the fixed deck once and `check_shear_HEL.py` compares `eq_2/S` of `PB3D_out.h5` with the direct formula −(D_θh¹²/h²² − D_θh²² h¹²/h²²²)/J from the stored metric derivatives, over the whole HELENA grid, tolerance 1e−4 on max|S|, plus a non-zero outboard-midplane value. Before the fix: deviation 0.79 (a boundary layer from both midplanes, decaying by 3.7 per interval; at θ = 0 the stored S is 1e−9 of max|S| where the formula gives 3.4e−2), FAIL. After: 1.47e−5 on the whole grid, midplane value 4.33e−2 from both routes, PASS. Registered in CMake when a python3 with h5py is found (`-DPython3_EXECUTABLE`); CI installs h5py.
- *Prediction met.* Fixed boundary ψ 0.45 to 0.75: 97 points with `magn_int_style = 1`: −4.3568189e−2 (predicted −4.35682e−2); Simpson 3/8 at 97, 193, 385, 100 points: −4.35691e−2, −4.35699e−2, −4.35714e−2, −4.35729e−2: converged to 1e−4 at 97 points where the unfixed code still moved by 5 % (−4.08e−2 to −4.28e−2).
- *POST energy balance* (`plot_E_rec`, POST_style 1, 101 poloidal points): E_pot/E_kin = −4.35343e−2 against ω² = −4.35682e−2, 0.078 % (plan: 0.07 %). On the unfixed run the same POST gives 0.17 %, not the 5 % of the report; the report's POST settings are unknown, so only the fixed value is a gate number.
- *L8, periodic branch.* HELENA cbm18a with `IAS = 1` needs `NCHI = 2^n`; at the same poloidal spacing (1024 points over 2π against 513 over π) PB3D gives −4.35487e−2 (periodic) against −4.35878e−2 (symmetric), 9e−4 apart, unchanged when HELENA's own poloidal resolution is doubled (NP 101 to 201: −4.35463e−2 to −4.35487e−2). The symmetric branch itself moves by 3.4e−4 between NCHI 401 and 513. So the two branches agree to 1e−3, the level of the discretisation, and no sharper; the residual is not the shear (both routes agree to 1e−5 now) and is left as an open item.
- *New anchors.* Fixed deck −4.3573e−2, free deck (ψ = 1.0, styles 2 and 4) −4.5144e−2, both at 1e−3 relative. The free deck at ψ = 1.0 still cannot see the vacuum (styles agree to 1e−9); step 4 moves it. The free style 2 run at ψ = 0.75 now gives −4.6377e−2 and style 4 −4.7684e−2: the vacuum sign error is unchanged and is step 3's business.
- *Not done.* The VMEC twin (ii) waits for step 5 as planned. MISHKA on the extended cbm18a fixture (q_a = 25, zero pressure beyond the pedestal) returns NaN from the first iteration with every density option tried, so C1/O1 on cbm18a is open; the circular cases work. COBRA (T8) is not built.
- *Learning for the plan.* PB3D and POST exit silently with status 0 when a path exceeds 120 characters (POST from a long scratch directory wrote nothing at all); run from short paths. POST needs a `Data/` directory for its HDF5 plots. The regression anchors are now at 1e−3, as rule 5 asked.

---

### Step 3: the vacuum sign and factor (axisymmetric)

Split in two because the sign is certain and the factor is not.

- **3a, the test that fails.** A vacuum can only add energy, so at the same boundary ω²(no vacuum) ≤ ω²(free) ≤ ω²(fixed). On cbm18a at ψ = 0.75 (edge displacement 0.1 of the maximum) the current code gives −4.3467e−2 for free, below the no-vacuum value −4.3143e−2. The test asserts the ordering and fails.
- **3b, the derivation.** On paper, with the actual definitions: A sums ∫PV with `alpha_fac` = 2π/n_alpha per line (so A corresponds to 2δW_p including the toroidal integral); `calc_vac_res` builds `res` from `EP^T Φ` with trapezoid weights `I_loc` and divides by μ0; the vacuum paper's surface form gives 2δW_v in terms of that `res`. The report's result is 2δW_v = −2π X_e^† res X_e, with the −1 certain (the paper's jump-term sign) and the 2π derived once. Write the derivation into the comment where the factor is applied.
- **Where the factor goes.** Once, at the end of `calc_vac_res`, so that `res` means "the surface-energy Hessian as it enters A". Then `set_BC_2`, `set_BC_4` and POST all use `res` unchanged, and the cylinder test's reference formula is multiplied by the same factor, which makes that test a check of the derivation. One line of physics code, one line in the test reference.
- **Hermiticity.** `res` differs from its conjugate transpose by 1e-5 relative, enough to make SLEPc reject the eigenvalue as complex at tight tolerance. The exact operator is Hermitian, so `res` is replaced by its Hermitian part, after printing the asymmetry as a diagnostic. Two lines.
- **POST.** The vacuum energy term has been inside `if (1.eq.2 .and. ...)` since February 2018, with a "NO VACUUM" print. Delete the guard and the print (two lines). POST then reproduces a free-boundary ω² from the eigenvector, with no shared code path for the surface term beyond `res` itself.
- **Checks, three independent of each other.**
  1. Ordering on cbm18a at 0.75. Predicted free ≈ −4.215e−2 with −2π (the report's rerun), inside the bounds. The pseudo-vacuum HELENA run (wall at ψ = 1) gives −4.198e−2 at 800 radial points; the no-wall result must not be more stable than that at the same resolution.
  2. The circular case R2. Fixed boundary stable (ω² > 0), free boundary unstable for m − 1 < nq_a < m. The marginal nq_a and the growth rates against R4 (ε = 0.05 variant) and against R3: compare ω²_free − ω²_fixed between PB3D and MISHKA, not the absolute values, because PB3D's high-n ordering and MISHKA's full MHD differ by a few per cent at n = 10 to 20 independently of the vacuum (the JCP comparison showed that gap).
  3. Energy balance. POST's E_pot/E_kin equals ω² for the free-boundary run to the fixed-boundary accuracy (0.07 %), and POST's vacuum term equals the cylinder formula's δW_v evaluated with the computed edge displacement on the ε = 0.05 case.
  If the three checks want different factors, the derivation is wrong and the step stops.
- **More checks (section 4).** L2: every eigenvalue of the final `res` is non-negative on cbm18a and on the circular case. L5: the Hellmann-Feynman slope. L9 and T2: MISHKA's wall sequence is monotone and follows the cylinder wall factor, and PB3D is at least as unstable as MISHKA at any wall. T3: ξ_r ∝ r^{m−1} and Φ ∝ r^{−m} on the ε = 0.05 case. T4: `res` on the circular torus against the mpmath toroidal-harmonic solution at finite ε, including the m ↔ m ± 1 coupling. T9: the high-n peeling criterion on the `PEELING` equilibrium. C4: the published cbm18 curves. C5: the public VACUUM matrix on the same boundary. C3, ELITE, only in step 7. O7: the dominant m per n against MISHKA.
- **Optional, measured, 3c.** The report found that using the vacuum-side solid angle instead of the plasma-side one lowers the cylinder errors for m = 2 to 5 from 1.3 to 1.9 % to 0.1 to 0.7 %. That is a one-sign change in the diagonal. Do it only if the error falls as predicted and the order of convergence does not change for the worse.

### Step 4: style 4 and the default

- **Bug.** `set_BC_4` inserts the (0,0) block with V1† + vac using INSERT_VALUES, then the V2 loop reaches jd = 0 and inserts the same block again with INSERT_VALUES, which replaces it. The vacuum and V1† are lost. At ψ = 0.75 style 4 gives −4.5574e−2 for any vacuum scale, 4.8 % from style 2.
- **Change.** Assemble the whole (0,0) block in `loc_block` first (V1† + vac, plus the jd = 0 part of V2), insert it once, and run the V2 loop from −st_size to −1. About six lines moved, no new logic. Keeping style 4 costs little and gives a second discretization of the same boundary condition as a permanent cross-check; deleting it would remove 80 lines of the original code.
- **Default.** Back to `[1,2]`, Toon's choice since 2017. One line.
- **Checks.** Style 2 and style 4 agree to 1e-5 at ψ = 0.75 with a nonzero edge displacement (they agree to 3.5e−8 today only because the edge displacement at ψ = 1 is 2e−8). The free-boundary regression deck moves to `max_r_sol = 0.75`, runs both styles, asserts the ordering of step 3a, and is anchored to 1e-3. The circular HELENA fixture (xz) is committed with a second regression deck anchored to R3 and R4. L6: the natural-boundary-condition residual of the style-2 eigenvector at the edge falls with resolution.

### Step 5: `setup_modes`

- **Bug.** With `X_style = 2` the mode window follows nq. `nm_X` puts mode number m in slot modulo(m − m_min(1), M) + 1, which depends only on m. `ind_tot` computes a new mode's slot from the slot of the mode created last plus the sub-step, which equals the same thing only while the window moves up. For nq < 0 (every VMEC case, because PB3D flips VMEC's q), reversed shear, or jitter in flat q, two active modes share a slot and the coupling coefficients are mixed without a message.
- **Change.** Two expressions: the new mode's slot is modulo(m_new − nm_X(1,1), M) + 1 in both branches, the same formula `nm_X` uses, so the two tables agree by construction. The 2020 rewrite on `origin/issue_vac_not_used_VMEC` is not used; it is unfinished and wrong for jumps of 2 or more.
- **Consistency check.** After the loop, assert for every normal point that the slot of each active mode in `ind_tot` is the slot where `nm_X` holds its number; CHCKERR otherwise. About six lines, O(M n_r), always on. This is the test that keeps the fix honest without extracting the routine into a unit test.
- **Checks.** (i) `pb3d/scripts/setup_modes_sim.py` from the report branch gets the new formula as a fourth variant: 0 wrong entries in all 8 sequences (the current code has 8, 1, 8, 4, 6 wrong entries in five of them). (ii) VMEC twins: the cbm18a and circular VMEC runs agree with their HELENA runs fixed-boundary to the fit accuracy. The size of the effect on eigenvalues is unknown today; the consistency check tells immediately whether a given run was affected.
- **More checks (section 4).** L3 on the VMEC twin: ω² must fall monotonically as `n_mod_X` grows, and mixed-up slots break that. L7: the poloidal and toroidal flux labels give the same fixed-boundary ω² (today's W7-X runs differ by 30 % at 60 radial points, which is either resolution or this bug), and shifting `min_alpha` changes nothing on an axisymmetric VMEC file. T8: COBRA from the same VMEC file as the n → ∞ limit.

### Step 6: the 3-D vacuum

- **6a, review.** The cloud commits 429f003 and d6d49e6 rewrote the singular integral (now the exact primitive, checked by the report to 8e−10), the far-field kernel sign, the row-sum diagonal and the line bridging. Read every hunk against the vacuum paper and the paper's known errors (singular-cell prefactor, jump-term sign, missing shell term, half solid angle). Verdict per hunk, written down. Anything that is not needed for correctness is reverted.
- **6b, α spacing and orientation.** PB3D's field lines use `calc_eqd_grid(..., excl_last=.true.)`: n_alpha lines with spacing (max − min)π/n_alpha. `calc_GH_1` and `calc_vac_res` use n_alpha − 1. The weights are too large by n_alpha/(n_alpha − 1) and line j's phase is off by n(j − 1)2π/(n_alpha(n_alpha − 1)). Change the two lines. `test_vac_3d` builds its own grid with both endpoints, which hides this; it gets PB3D's convention instead.
  The test that decides orientation and weights at once: a manufactured exterior potential. Put a point source inside the surface, Φ = 1/|x − x0|, harmonic outside and decaying. Its normal derivative on the surface is exact, so the BEM solve must return Φ on the surface to discretization accuracy, and the error must fall with resolution at the BEM's order. Run it on the circular torus and on a rippled torus a(ζ) = a0(1 + δ cos Nζ) with field-line-like coordinates ζ = α + qθ. The surface need not be an equilibrium; only G and H are under test. This settles the style-1 normal orientation (report 13, "inferred") and the empirical factor `−res/((max_alpha − min_alpha)π)` in `calc_vac_res`, which was fitted to style 2 with the wrong grid.
- **6c, Richardson levels.** `store_vac_VMEC` copies the old G and H into the refined arrays and `calc_GH_1` adds with `+=`. For style 1, skip the copy (one line) or zero before accumulation. Test: a two-level run's level-2 `res` equals a one-level run at the level-2 resolution, bit for bit.
- **6d, aliasing.** The n_alpha-point sum aliases every toroidal harmonic that is a multiple of n_alpha into the average. A warning when n_alpha divides the field period. Five lines in the input checks.
- **6e, end to end at ε = 0.** The circular VMEC twin, free boundary through the style-1 3-D vacuum (alpha_style 2, n_alpha ≥ 4), against the HELENA style-2 result and R3. This is the only full-path check of the 3-D vacuum with a known answer, and it exercises `setup_modes` with nq < 0, the VMEC equilibrium path and the field-line quadrature together.
- **6f, continuity.** Ripple the circular VMEC boundary with PB3D's own perturbation tools (`pert_ripple*.dat`, `create_VMEC_input`), rerun VMEC, and check ω²(δ) − ω²(0) ∝ δ² with no jump, which is what finding 2 of the papers review predicts. Then W7-X (`wout_w7x.nc`, B-reversed copy) as a smoke test with the ordering of step 3a, style 2 only, and the STRUMPACK fallback logged.
- **More checks (section 4).** T5: the full family of Green identities with the order measured. L2: `res` positive semi-definite on the rippled torus and on W7-X. L7: rotating a 3-D boundary by 2π/N_fp leaves `res` unchanged. L10: n_alpha convergence and the aliasing jump. L11: PB3D at large n at most COBRA's worst field line on W7-X and LHD. C6 first: the open-source 3-D VACUUM code [R3] on the same boundary, element by element. C9: an independent BEM on the same surfaces. C5: the axisymmetric VACUUM matrix on the circular case as the ε = 0 anchor of C6. C7 and C8, CASTOR3D and CAS3D on the rippled circular case and W7-X, free minus fixed, only in step 7. T7: the δ² coefficient from the rippled cylinder, optional.

### Step 7: external confirmation, changelog and PR

Only now go the requests to other groups (section 4.5): one ELITE run of cbm18a (C3), CASTOR3D on the rippled circular case (C7), CAS3D on W7-X (C8). They confirm what the local catalogue has already pinned; a disagreement beyond the stated blind spots reopens the step concerned. The PR does not wait for them.

A 2.49 entry with one line per fix and the numbers that changed. The cloud session's undocumented changes (default `[1,4]`, style 3 removed, sign of the response flipped) get their lines too. PR to `master` with a description that corrects PRs 1 to 3 ("4 % more unstable than fixed" measured the longer domain; "style 2 agrees with style 4" held only at zero edge displacement; "3-D verified end to end" was the synthetic torus). Before the PR: O3, the mutation check, so that the suite is shown to catch each bug it was written for. Toon merges.

### Step 8: where the analytic 3-D equilibria come in

Not a bug fix, so it comes after. The Landreman and Issan-Citrin-Landreman-Tracey equilibria give exact B, ψ, ι, p, curvature and shear on a strongly 3-D torus. PB3D reads only VMEC and HELENA, and VMEC does not converge on the exact-pressure cases (Mercier-unstable; `desc-vs-vmec/stability/REPORT.md`). DESC does, and DESC can write a `wout` file. So: DESC solution of family 2 at the paper's parameters, written as `wout`, read by PB3D, and PB3D's ι, S, κ_n, κ_g, σ compared with the closed forms from `desc-vs-vmec/scripts/landreman.py`. That is a test of the 3-D equilibrium path (`eq_ops`, `VMEC_ops`, `grid_utilities`) that nothing else gives. It is not a vacuum test: for the vacuum, step 6b's manufactured potential on any surface is stronger, because it has an exact answer. T10 applies here: the exact-pressure cases are Mercier-unstable on every surface, so PB3D must find them unstable at large n.

## 4. Validation catalogue

Each check has an id, says what it can see and what it is blind to, what it costs, and which step uses it. "Local" means it runs on this Mac or on CI with what is in the archive; "ask" means another group runs a code. Section 3 cites the ids.

### 4.1 Theory: exact or asymptotic answers

| Id | Check | Sees | Blind to | Cost | Steps |
|---|---|---|---|---|---|
| T1 | Shafranov's uniform-current external kink without a wall: δW, growth rate and the window m − 1 < nq_a < m, on the ε = 0.05 HELENA case | vacuum sign and factor, edge physics, inertia normalisation | toroidal effects | sympy, half a day | 3 |
| T2 | The same with a conformal wall at radius b: the vacuum term carries the factor (1 + (a/b)^{2m})/(1 − (a/b)^{2m}). Applied to MISHKA's `RWALL` sequence before MISHKA is used as a reference | the reference itself | | hours | 0, 3 |
| T3 | Eigenfunction shapes in the cylinder: ξ_r ∝ r^{m−1} inside, Φ ∝ r^{−m} outside (re-enable `vac_pot_plot` in POST, one line) | the eigenvector and the exterior field, not only ω² | | hours | 3 |
| T4 | Exterior Laplace solution on a circular torus by toroidal harmonics in mpmath at 30 digits: `res` at finite ε and finite n, including the m ↔ m ± 1 coupling. The exterior Dirichlet-to-Neumann map as a series in toroidal harmonics is written out by Ashtab, Morais and Porter 2022 [R6]; Jones, Chacón et al. 2026 [R7] use U = √(2 cosh ρ − 2 cos θ) P^n_{m−1/2}(cosh ρ) e^{i(mθ + nφ)} on a torus with R = 1.65, r = 1 as an exact exterior test and see order 2 | style 2 at finite aspect ratio, `dtorh`, the solid-angle side | non-circular surfaces | a day | 3, 6 |
| T11 | `dtorh` at high n: Chance, Turnbull and Snyder 2006 [R5] found that building the toroidal functions by upward recursion in n loses precision already at the moderately high n of peeling-ballooning modes, and reach 12 to 16 digits by direct integration with contour deformation. Check `dtorh` against mpmath at n = 20 to 80 and at the near-singular arguments the BEM uses | the n range PB3D is for | | hours | 3 |
| T5 | Green identities for a family of potentials: interior harmonics (x, y, z, Re and Im of (x + iy)^k) and exterior multipoles centred inside (1/r and its derivatives): H φ = G ∂φ/∂n and (H + 4π) φ = G ∂φ/∂n, with the convergence order measured | G, H, diagonal, orientation, weights, on any surface | the projection onto modes | partly exists; a day | 6 |
| T6 | Manufactured exterior potential (point source inside) on a torus and a rippled torus: the BEM solve returns Φ on the surface to discretization accuracy | the whole BEM solve and the normal orientation | the equilibrium | in the plan | 6 |
| T7 | Rippled cylinder, r = a(1 + δ cos(Nζ − lθ)): boundary perturbation theory of the exterior Neumann problem, carried to order δ², gives the δ² coefficient of the α-averaged vacuum energy | the only analytic 3-D vacuum number | finite aspect ratio | sympy, two days; optional | 6f |
| T8 | n → ∞ limits: COBRA (C2) and HELENA's own infinite-n ballooning flag per surface. PB3D extrapolated in 1/n must reach COBRA on the axisymmetric case again (JCP 2017: 0.454 to 0.459 against 0.453) | fixed-boundary physics after the shear fix; the whole VMEC path | the vacuum | hours once COBRA runs | 2, 5 |
| T9 | The high-n peeling criterion (Connor, Hastie, Wilson, Miller 1998) on an equilibrium with edge current (`Simulations/PEELING/HELENA`): the marginal edge current from the criterion against PB3D's marginal point | free-boundary edge physics in PB3D's own ordering | | a day of derivation plus a HELENA scan | 3, 4 |
| T10 | Mercier: every surface with DMerc < 0 in VMEC or DESC must be unstable in PB3D at large n (a necessary condition) | the 3-D equilibrium path and the high-n limit | | hours | 8 |

### 4.2 Logical conditions: inequalities and invariances

| Id | Condition | Catches | Cost | Steps |
|---|---|---|---|---|
| L1 | ω²(res = 0) ≤ ω²(free) ≤ ω²(fixed) at the same boundary, and ω²(free) ≤ the pseudo-vacuum fixed result at ψ = 1 | the vacuum sign | in the plan | 3, 4, 6 |
| L2 | The vacuum energy is ∫\|∇Φ\|² ≥ 0, so the final `res` is Hermitian positive semi-definite: every eigenvalue ≥ 0, on cbm18a, the circular case, the rippled torus and W7-X | sign, Hermiticity, orientation in 3-D (the 2019 referee's unanswered question 5) | minutes per case | 3, 6 |
| L3 | Rayleigh-Ritz monotonicity: ω² cannot rise when the trial space grows: in `max_r_sol` for a fixed boundary, in `n_mod_X`, in radial points on nested grids | mixed-up couplings (`setup_modes`), broken Hermiticity | hours | 2, 5 |
| L4 | A Hermitian, B positive definite, ω² real: an eigenvalue rejected as complex is a bug signal, never noise | any broken symmetry (style 4 on W7-X was one) | free, every run | all |
| L5 | Hellmann-Feynman: with `res` scaled by s, dω²/ds at s = 1 equals X_e† res X_e / (X† B X) from the eigenvector. The energy-balance check in differential form | the factor's placement and the eigenvector normalisation | two runs | 3 |
| L6 | The natural boundary condition (V1† X + V2 X′ + res X) evaluated at the edge from the converged style-2 eigenvector falls with resolution; style 4 satisfies it by construction | the variational imposition | a script on the HDF5 output | 4 |
| L7 | Invariances, each a bit-identical or 1e-10 pair: n → −n; poloidal against toroidal flux label (fixed boundary); shifts of `min_alpha` and of `min_par_X` on axisymmetric equilibria; rigid toroidal rotation of a 3-D boundary by 2π/N_fp (`res` unchanged) and by any angle (`res` changes by a phase only); scaling of R0, B0 and ρ0 (normalised ω² unchanged); 1, 2 and 4 processes; STRUMPACK against ScaLAPACK; Richardson restart against a fresh run | sign conventions, index bugs, parallel bugs. Today's W7-X runs already fail the flux-label pair | hours; becomes a CI suite | 1, 5, 6 |
| L8 | An up-down symmetric HELENA equilibrium run once with `ias = 0` (symmetric branch, where the shear bug lives) and once with `ias = 1` (periodic branch, never affected) gives the same ω² | step 2, independently of VMEC | one extra HELENA run | 2 |
| L9 | Wall monotonicity: MISHKA's ω² falls as `RWALL` grows; PB3D without a wall is at least as unstable as MISHKA at any wall | gross vacuum errors | with R3 | 3 |
| L10 | n_alpha convergence: results converge as n_alpha grows when n_alpha is coprime to N_fp, and jump when n_alpha divides N_fp | aliasing (report 11) | hours | 6d |
| L11 | PB3D's n → ∞ growth rate on a 3-D equilibrium is at most COBRA's worst-field-line value, and equals it at ε = 0 | the α-averaging as implemented | with C2 | 6f, 8 |
| L12 | The measured convergence order equals the order of the discretization: faster than any power in parallel points, order 3 radially (measured 1 today), order 1 or 2 for the BEM | silent errors such as the radial one | hours per direction | 0, 6 |

### 4.3 Other codes

| Id | Code | Access | Validates | Blind spot | Ask |
|---|---|---|---|---|---|
| C1 | MISHKA-1 | source in the archive, local build | axisymmetric fixed and free (conformal wall) eigenvalues, dominant m | full MHD against PB3D's high-n ordering: compare free − fixed and the marginal points, not absolute values | none |
| C2 | COBRAVMEC | STELLOPT clone at `~/Code/stellopt` | n = ∞ ballooning on VMEC files, 2-D and 3-D; the upper bound L11 | no vacuum, no finite n | none |
| C3 | ELITE | Wilson (York), Snyder | the same high-n peeling-ballooning ordering as PB3D, free boundary, no wall: the closest axisymmetric reference that exists; cbm18 is its standard case | 3-D | one cbm18a run, n = 5 to 30, free boundary, no wall |
| C4 | Published cbm18 benchmarks: Ferraro, Jardin and Snyder 2010 [R1] give ELITE, GATO and M3D-C1 growth rates against n (Fig. 5, normalised to ω_A built from B0 and R0 on axis and the density at the pedestal top) for cbm18: circular, pedestal foot at the plasma-vacuum interface. ELITE assumes an infinite vacuum, so its curve is a no-wall reference like PB3D; GATO and M3D-C1 use a conformal conducting wall with R = 3 m, a = 2 m, Ψ_w = 10/7, which gives the wall effect for free. The paper warns that slight differences in mapping the cbm18 equilibrium change the growth rates noticeably. Snyder et al. 2002 [R2] is the original ELITE and GATO cbm18 work. ELITE has itself been benchmarked against GATO, MISHKA, DCON, CASTOR, MARG2D, ILSA, MARS, BOUT++ and NIMROD [R1] | the papers | free-boundary growth rates against n at zero cost | PB3D has no wall, so MISHKA at the same wall flux is the bridge; the normalisation must be pinned (the 2017 referee exchange had a factor 1.5 in R/B) and the equilibrium confirmed to be the same cbm18 | none |
| C5 | VACUUM (Chance) inside GPEC, public: github.com/PrincetonUniversity/GPEC (Fortran, release 1.5.6, June 2025). Method: Chance 1997 [R4], Green's second identity with collocation, the same formulation as style 2; high-n Green's function in [R5] | public | the vacuum response matrix W_V in the Fourier basis e^{i(mθ − nζ)} on an axisymmetric boundary: element by element against `res`, far sharper than any eigenvalue | 3-D | none; build it |
| C6 | 3-D VACUUM: Halpern, Logan, Nührenberg, Paul and Paz-Soldan, arXiv:2610.00536, 30 September 2026 [R3]: the VACUUM method extended to non-axisymmetric boundaries in the open-source Julia GPEC 2.0 at github.com/OpenFUSIONToolkit/GPEC. Response matrix in the basis e^{i(mθ − nζ)} on straight-field-line angles, no wall or a conformal wall, order 4 convergence; agrees with CAS3D on the W7-X standard configuration to 1 % in the eigenvalues and 2 % in the Frobenius norm, and with VACUUM on a Solov'ev tokamak to 1e−5. STARWALL (Merkel and Strumberger 2015 [R8], IPP) is the fallback if the Julia release lags | public (release pending); STARWALL by request | the 3-D vacuum response matrix for a given boundary without a wall: the matrix-level 3-D reference. PB3D's field-line modes exp i(nα + (nq − m)θ) = exp i(nζ − mθ) are the same basis, so `res` compares element by element once the normalisation is settled | the plasma side | none for the Julia code; one boundary plus the basis definition for STARWALL |
| C7 | CASTOR3D (Strumberger and Günter 2017 [R9]) | Strumberger, Günter (IPP) | 3-D free-boundary eigenvalues with the vacuum, mode families included | PB3D has no family coupling, so agreement only at O(δ⁰) and O(δ¹); compare free − fixed | the rippled circular VMEC case at n = 10, 15, 20, fixed and free, no wall |
| C8 | CAS3D (Nührenberg, IPP Greifswald), the reference used in [R3]; TERPSICHORE (Lazerson, Cooper) second, because it treats the vacuum as a pseudo-plasma with a conformal wall, which is also what Toon's extended HELENA equilibria do | ask | 3-D free-boundary eigenvalues with the vacuum; TERPSICHORE is the GIGA paper's own tool | wall and pseudo-plasma for TERPSICHORE; families for both | CAS3D: the W7-X standard case, no wall, the n = 10 family. TERPSICHORE: the same with the wall at 30 % and at 300 % |
| C9 | An independent BEM: an open-source Laplace boundary-element library on triangulated surfaces (bempp-cl, or a hand-written Nyström code) | pip, local | the Neumann-to-Dirichlet map for PB3D's actual normal data on the rippled torus and the W7-X boundary; style 1's Φ and `res` against it | nothing beyond the BEM | two days |

The matrix-level references (C5, C6, C9) compare `res` itself and come first; C5 and C6 are public, so no strings need pulling for them. The eigenvalue codes (C1, C3, C7, C8) also see the ordering differences between the codes, so they confirm rather than pin.

### 4.4 Other ways

| Id | Way | Steps |
|---|---|---|
| O1 | Toon's own published numbers: the JCP 2017 PB3D-against-MISHKA curve for cbm18a (fixed boundary, ψ 0.2 to 0.9, 500 radial points) predates the June 2018 shear bug and must be reproduced after step 2 with the data-file settings | 2 |
| O2 | Golden intermediate arrays: `res`, G and H on the fixtures at fixed resolution, stored and compared at every step, not only ω² | all |
| O3 | Mutation check of the test suite: at the end, re-introduce each fixed bug one at a time; every one must fail at least one test, or that test is not a test | 7 |
| O4 | Radial-order diagnosis: vary HELENA's `NRMAP` (201, 401, 801) at fixed `n_r_sol`, and `n_r_sol` at fixed `NRMAP`, to find which grid limits the order; the report's candidate is the 158-point equilibrium grid | 0 |
| O5 | Memory sweeps on every regression deck in Debug: `-fcheck=all`, the CI sanitizer job, Guard Malloc once per step | 1, all |
| O6 | The `setup_modes` simulation script from the report branch as the index-logic test, extended with the new variant | 5 |
| O7 | Dominant-mode structure: PB3D's dominant m per n against MISHKA's, and the radial position of the eigenfunction peak | 3 |

### 4.5 What to ask whom, and when

All of these are step 7 items, sent after steps 1 to 6 have passed their local gates. The W7-X and Solov'ev inputs of [R3] are the one exception, because they only make the local C6 check easier and can be requested at any time.

| Group | Ask | Gives |
|---|---|---|
| IPP Garching (Strumberger, Hölzl, Günter) | CASTOR3D free and fixed at three n on the rippled circular case; STARWALL only if the Julia 3-D VACUUM is not released | C7, C6 |
| IPP Greifswald (Nührenberg) | CAS3D on the W7-X standard case without a wall, the n = 10 family | C8 |
| York / GA (Wilson, Snyder) | one ELITE run of cbm18a, free boundary, no wall, n = 5 to 30 | C3 |
| Gauss Fusion (Lazerson) | TERPSICHORE on W7-X with the wall at 30 % and at 300 %; they get the PB3D comparison they asked for in return | C8 |
| Columbia / PPPL (Halpern, Logan) | the W7-X and Solov'ev benchmark inputs of [R3], offered on request in the paper | C6 |

### 4.6 References found in the literature search of 27 October 2026

- [R1] N. M. Ferraro, S. C. Jardin, P. B. Snyder, "Ideal and resistive edge stability calculations with M3D-C1", Phys. Plasmas 17, 102508 (2010). cbm18, dbm18 and Meudas1 benchmarks; ELITE growth rates against n; wall and interface sensitivity. Local copy of the text in the scratchpad of this session.
- [R2] P. B. Snyder et al., "Edge localized modes and the pedestal: a model based on coupled peeling-ballooning modes", Phys. Plasmas 9, 2037 (2002). ELITE formalism companion: H. R. Wilson et al., Phys. Plasmas 9, 1277 (2002).
- [R3] J. Halpern, N. Logan, C. Nührenberg, E. Paul, C. Paz-Soldan, "Rapid vacuum response calculations for non-axisymmetric plasma geometries", arXiv:2610.00536 (30 Sep 2026). Code: github.com/OpenFUSIONToolkit/GPEC. Singular quadrature from Malhotra et al., J. Comput. Phys. 2019 and Plasma Phys. Control. Fusion 2019.
- [R4] M. S. Chance, "Vacuum calculations in azimuthally symmetric geometry", PPPL-3216 (1997); Phys. Plasmas 4, 2161.
- [R5] M. S. Chance, A. D. Turnbull, P. B. Snyder, "Calculation of the vacuum Green's function valid even for high toroidal mode number in tokamaks", PPPL-4152 (2006).
- [R6] Z. Ashtab, J. Morais, R. M. Porter, "Neumann problem on a torus", arXiv:2212.00224 (2022).
- [R7] S. Jones, L. Chacón, J. Hamilton, D. Barnes, A. Yagüe-López, "Implementation and verification of toroidal resistive wall boundary conditions in the PIXIE3D MHD code using a boundary integral method", arXiv:2606.05446 (2026). Section 6.1: the exact exterior toroidal-harmonic test.
- [R8] P. Merkel, E. Strumberger, "Linear MHD stability studies with the STARWALL code", arXiv:1508.04911 (2015).
- [R9] E. Strumberger, S. Günter, "CASTOR3D: linear stability studies for 2D and 3D tokamak equilibria", Nucl. Fusion 57, 016032 (2017).
- [R10] J. W. Connor, R. J. Hastie, H. R. Wilson, R. L. Miller, "Magnetohydrodynamic stability of tokamak edge plasmas", Phys. Plasmas 5, 2687 (1998). The high-n peeling criterion of T9.
- [R11] Bempp-cl, open-source boundary element library with Laplace single- and double-layer operators, pypi.org/project/Bempp-cl (C9).
- Not found: a published table of cbm18 growth rates with numbers; a published analytic expansion of the toroidal corrections to the external-kink vacuum energy. Fig. 5 of [R1] has to be digitised, and T7 has to be derived.

## 5. What each step must not do

- Step 1 must not change `spline_real` unless the build-script patch proves insufficient.
- Step 2 must not touch the VMEC shear or the asymmetric HELENA branch (periodic conditions, not affected).
- Step 3 must not move the factor into `set_BC_2`/`set_BC_4` (two places to keep in sync, and POST would need a third).
- Step 4 must not delete style 4 or add a style 5.
- Step 5 must not merge the 2020 rewrite.
- Step 6 must not touch the axisymmetric style 2, which is anchored by then.

## 6. Line count

| Step | Physics lines changed | Where |
|---|---|---|
| 1 | 0 | build script, test scripts, CI |
| 2 | 1 (+1 comment) | `eq_ops.f90` |
| 3 | 3 (+ derivation comment) | `vac_ops.f90` (factor, symmetrise), `sol_ops.f90` (guard removed) |
| 4 | about 6 moved, 1 default | `SLEPC_ops.f90`, `input_ops.f90` |
| 5 | 2 expressions + about 6 | `X_ops.f90` |
| 6 | 2 + 1 to 3 + about 5 | `vac_ops.f90`, `vac_utilities.f90`, `input_ops.f90` |
| total | under 40 | |

Tests grow by a few hundred lines, all under `tests/`, plus two committed HELENA fixtures (circular 2.25 and its ε = 0.05 variant, xz).

## 7. Critique log

**Round 1, against the report's order (crashes, shear, sign, modes, CI, radial, 3-D).**
- No external reference until late, so the sign fix would rest on bounds alone. Changed: step 0 builds MISHKA, the circular case and the cylinder formula before any code changes.
- CI was repaired in step 5 of the report's list, but the Mac cannot run `mpirun`, so until CI is green nothing multi-process is tested. Changed: CI goes into step 1.
- Sign and factor were one step. The sign follows from bounds; the factor needs a derivation and two numerical pins. Changed: 3a and 3b, with three independent checks that must agree.
- "Repair or remove style 4" was open. Changed: repair, because the fix is six lines and the agreement of two discretizations is a free permanent test; removing it deletes original code.
- The 3-D step had no exact reference. Changed: the manufactured exterior potential on a torus and a rippled torus (6b), and the ε = 0 end-to-end comparison through the VMEC path (6e).
- Richardson exponent "fix" dropped: the parallel trapezoid rule on a periodic integrand converges faster than any power, so there is no exponent to fix. Richardson stays off by default and documented as such. The double count (6c) is still fixed because it is one line and has a bit-exact test.
- Radial convergence (order 1, two stencils 2.7 % apart) was a gate. Changed: not a gate, but every benchmark comparison runs 200, 400 and 800 radial points and reports the extrapolated value with its spread. If the cause turns out to be a bug it becomes its own step.

**Round 2, against round 1.**
- MISHKA versus PB3D absolute values differ by a few per cent from the ordering alone (JCP 2017, gap growing with n). Changed: compare ω²_free − ω²_fixed between the codes, and the marginal nq_a, not absolute values.
- The circular case at ε = 0.33 is not a cylinder. Changed: an ε = 0.05 variant for the analytic formula; the 0.33 case stays for MISHKA.
- The Landreman surface as a BEM test geometry was proposed. Dropped for the vacuum: G and H do not care whether the surface is an equilibrium, and a rippled torus is simpler and equally exact. The analytic equilibria move to step 8, the 3-D equilibrium path, where they are the only exact reference.
- Predictions were missing. Changed: every physics step states its expected number from the report's runs before running.
- The `res` Hermiticity was in "smaller items". It makes SLEPc reject eigenvalues at tight tolerance, which hides inside every free-boundary convergence study. Moved into step 3.
- The `setup_modes` fix was "the corrected logic of section 10" (search the extreme again). Replaced by a smaller and more obviously right change: use the slot formula `nm_X` already uses, so the two tables agree by construction, plus a runtime assertion instead of a unit test that would need grids and equilibria.

**Round 3, against round 2.**
- Bugs found by reading (6b, 6c, 13) could be "fixed" without being bugs. Changed: rule 2, the test must fail first. For 6b the manufactured potential fails on the current grid convention or the bug is reconsidered; for 6c the bit-exact level test; for the orientation, the manufactured potential with the sign as coded.
- The empirical factor `−res/((max_alpha − min_alpha)π)` in `calc_vac_res` was fitted to style 2 using the test's grid that includes both endpoints. After 6b it may change. Changed: 6b re-derives it from the manufactured potential rather than keeping the fitted value.
- The MISHKA build is a real risk (ILP64). Changed: named as the risk of step 0 with the fallback (cylinder formula plus bounds now, MISHKA on Linux CI later).
- The VMEC twin check in step 2 depends on `setup_modes` (step 5) being right for nq < 0. Changed: noted in step 2; a disagreement pulls step 5 forward.
- Fixtures: the circular HELENA files must be committed (xz) or the regression cannot run on CI. Added to step 4.

**Round 4, Toon's comment: more validation, by theory, by logic, by other codes and by other means.**
- The validation matrix became a catalogue (section 4): 38 checks of four kinds, each with what it sees, what it is blind to and its cost. Every step gained checks from it.
- Three checks that are both cheap and sharp were missing: `res` positive semi-definite (L2), the Hellmann-Feynman slope (L5), and the symmetric-against-periodic HELENA pair (L8), which tests the shear fix with HELENA alone.
- Matrix-level references (C5 VACUUM, C6 STARWALL, C9 an independent BEM) compare `res` itself instead of eigenvalues. They are now the preferred external checks for the vacuum, with the eigenvalue codes (C1, C3, C7, C8) second, because those also see the ordering differences between the codes.
- Toon's own 2017 numbers (O1) and the published cbm18 curves (C4) cost nothing and were not in the plan.
- A mutation check of the test suite (O3) was added, because a suite that has never failed proves nothing.

**Round 5, after the literature search.**
- Chance's VACUUM code is public inside GPEC, and a September 2026 preprint [R3] extends it to non-axisymmetric boundaries as open source, in the same Fourier basis PB3D's field-line modes span. The strongest 3-D vacuum reference therefore needs no strings pulled. C5 and C6 were rewritten and now come before C9 in step 6.
- TERPSICHORE treats the vacuum as a pseudo-plasma with a conformal wall, so it cannot anchor a no-wall vacuum. CAS3D takes its place in C8; TERPSICHORE stays as the GIGA paper's own tool.
- The cbm18 literature reference is Ferraro, Jardin and Snyder 2010 [R1], not the Burke paper the plan had guessed; it uses a conformal wall at Ψ_w = 10/7, so MISHKA at that wall flux is the bridge to PB3D's no-wall result.
- Chance 2006 [R5] warns that toroidal functions built by upward recursion in n lose precision at the n of peeling-ballooning modes. A check of `dtorh` at n = 20 to 80 (T11) was added to step 3.
- Two exterior-torus references make T4 concrete: the toroidal-harmonic Dirichlet-to-Neumann series [R6] and PIXIE3D's exact exterior test [R7].

**Round 6, Toon's comment: everything local before CAS3D and the other asks.**
- Ground rule 8 added. The requests to other groups moved into step 7 as confirmation; no step waits on them. With VACUUM and its 3-D extension public, the local catalogue already contains a matrix-level reference for both the axisymmetric and the 3-D vacuum, so the asks confirm rather than decide.

## 8. Out of scope, in order of value afterwards

1. Negative-flux VMEC files (W7-X, GIGA): `origin/flip_angles` from 2018, or the B-reversal workaround made official. Needed before any real stellarator work.
2. Step 8, the 3-D equilibrium path against the analytic equilibria.
3. Family coupling n ± N_fp (the GIGA question).
4. The radial-convergence cause.
5. Toroidal flux in `store_vac`; G and H persisted for restarts.
6. STRUMPACK HSS compression failing on 3-D kernels (LU fallback works; performance only).

## 9. What I am still not happy with

- (added after step 2) The symmetric (`ias = 0`) and periodic (`ias = 1`) HELENA branches give ω² 9e−4 apart at equal poloidal spacing on cbm18a; the shear is not the cause. Not pursued in this plan.
- (added after step 2) MISHKA returns NaN on the extended cbm18a fixture (zero pressure beyond the pedestal), so the cbm18a fixed-boundary cross-check C1/O1 is open; the circular cases are fine.


- Step 3's factor rests on a derivation I have not done yet and on two numerical pins. If the pins disagree with the derivation and with each other, the step stops and the vacuum stays off until the derivation is settled; there is no shortcut.
- The 3-D vacuum has no exact end-to-end reference beyond ε = 0 and smoothness in δ. That is a property of the problem, not of the plan.
- The radial convergence caps every benchmark at about 3 % until its cause is known.
- MISHKA may not build on arm64 macOS in a day. ELITE (C3) and the published cbm18 curves (C4) now cover the axisymmetric vacuum if it does not.
