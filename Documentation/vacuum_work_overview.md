# PB3D vacuum work: overview

Branch `dev/vacuum-fixes`, draft PR #5 to `master`, tracking issue #4. Seventy-three commits in two stretches: the modernisation of July 2026 and the vacuum work of 27 October to 1 November 2026. The detailed log is `Documentation/plan_vacuum_bugfix.md`; this page is the outcome only.

Two different "styles" appear below. `vac%style` is the vacuum method, set by the equilibrium code: `vac%style` 1 is the 3-D boundary-element method on field lines (VMEC), `vac%style` 2 the axisymmetric one with toroidal harmonics (HELENA). `BC_style` is the input for the boundary conditions; its second value is the condition at the plasma edge: 1 fixed boundary, 2 the free boundary through the variational (Hermitian) surface energy, 4 the free boundary through an explicit natural-boundary-condition row.

## 1. Where we started, where we are

| | Before (PB3D 2.47, 2018 to 2020 code) | Now (2.49, this branch) |
|---|---|---|
| Build | PETSc 3.9 API, no longer buildable on a current machine | PETSc/SLEPc 3.25 on macOS arm64, Snellius and GitHub Actions |
| Tests | none that ran | unit, full-stack and regression layers under CTest; eight CI jobs, all green |
| Axisymmetric vacuum (HELENA, `vac%style` 2) | solved the interior problem; result never checked | exterior problem, Hermitian response, corner term on the right side of the boundary (2 November); free-boundary eigenvalues within 0.4 % of TERPSICHORE on the circular kink at the fixtures' 257 poloidal points and 0.2 % at 1025, and within 0.5 % on the cbm18a pedestal mode at n = 10 |
| 3-D vacuum (VMEC, field lines, `vac%style` 1) | code present, never gave a right answer | converges; reproduces TERPSICHORE's free-boundary eigenvalue to 2e−5 on the circular twin and runs on a rippled torus and on W7-X |
| Fixed boundary | the one result that was trusted, but 6.5 % too stable through the HELENA shear bug and 2e−5 process-count dependent | HELENA and VMEC paths agree to 2e−6, process counts to 1.5e−8, TERPSICHORE to 0.1 % on the circular case and 0.5 % on the cbm18a pedestal mode at n = 10 |

## 2. Changes to the PB3D and POST code

Fixed means the code was wrong and now does what it was meant to do; new means the capability did not exist.

| Where | Fixed | New |
|---|---|---|
| `SLEPC_ops`, `SLEPC_utilities` | port from the PETSc 3.9 API to 3.25 (July) | |
| `vac_ops` `solve_Phi_BEM` | the axisymmetric vacuum (`vac%style` 2) solved the interior system; it now solves the exterior one, (H + 4πI)Φ = G R (July) | |
| `vac_ops` `calc_GH_2` | the corner term of the H diagonal used the boundary's angle on the vacuum side, β, where the interior operator that `calc_GH_2` assembles needs the plasma side's, 2π − β; the diagonal was off by about 8π/N for N boundary points, and the response `vac%res`(m, m) 0.7 % too large for every m at the fixtures' 257 points (2 November) | |
| `vac_utilities` `calc_GH_int_2`, `vac_ops` `calc_GH_2`, `solve_Phi_BEM`, `calc_vac_res` | the two styles stored H in opposite sign conventions (the axisymmetric one with the plasma's normal in the kernel and a +4π shift at solve time, the field-line one with the vacuum's) and fed plasma-side data, compensated by a sign flip of the 3-D response; now one vacuum-side convention: kernels, corner term −2β and data −∂X/∂θ all along the stored normal −J∇ψ, no shift, response factor +2π; results unchanged, two style branches gone (5 November) | |
| test layer, `tests/CMakeLists.txt` | the regression tests passed on any "PASS" in their output, so the one `check_vac_res.py` prints hid a failed anchor; they now also fail on "FAIL" (2 November) | |
| `vac_utilities` `calc_GH_int_1`, `vac_ops` `calc_GH_1` | singular interval integral, quadrature weights, sign and a `**(−3)` far-field factor in the 3-D kernels (July) | near-field integration along the field line, Gauss–Legendre on the cubic interpolant through four points (31 October) |
| `vac_ops` `calc_GH_1` | α spacing 2π/n_alpha instead of 2π/(n_alpha − 1); G and H no longer counted twice across Richardson levels | G diagonal from the identity for the linear interior harmonics, as the row sum sets the H diagonal |
| `vac_ops` `store_vac_VMEC` | handed the kernels the outward normal where they take the inward one; the 3-D response was 3 to 4 times too small | |
| `vac_ops` | | `check_vac_exact` and the input `check_vac`: Green identities and exact exterior solves on the real surface, per-point files, condition estimate |
| `SLEPC_ops` energy assembly | the vacuum entered with the wrong sign and without the factor 2π, so it destabilised; now −2π X_e^H res X_e with the Hermitian part of the response | |
| `SLEPC_ops` `fill_mat` | the free end point carried the full radial weight h instead of h/2; free-boundary eigenvalues converged at first order in `n_r_sol` | |
| `SLEPC_utilities` boundary conditions | `BC_style` 4 overwrote the vacuum block; `BC_style` 3 removed; default back to `BC_style = 1 2` | `BC_style` 2, the variational (Hermitian) free-boundary condition (July) |
| `eq_ops` `calc_derived_q`, HELENA branch | the shear used the even spline condition on the odd Θ^α since June 2018; cbm18a was 6.5 % too stable | |
| `X_ops`, `X_utilities` `setup_modes` | the slot of a mode was wrong for nq < 0; now from the mode number, with a consistency check | |
| `VMEC_ops` `read_VMEC` | Jacobian, |B| and covariant B taken from the Nyquist arrays as if they followed the main mode list (37 % error on VMEC++ files with ntor > 0); matched by (m, n) | |
| `VMEC_ops` normalisation | B_0 was |B| on the axis at θ = ζ = 0, which varies with the ripple; now the toroidal average | |
| `VMEC_utilities` `calc_trigon_factors`, `fourier2real` | | one-angle basis, full Nyquist set, one code path for both uses; half the equilibrium time |
| `MPI_utilities` `redistribute_var` | the overlap of two process ranges came from the wrong rank; eigenvalues shifted 2e−5 with 2 or 4 processes, now 1.5e−8 | |
| `vac_utilities` | | distributed vector helpers in one place (consolidation) |
| POST (`sol_ops`, `driver_POST`) | the vacuum energy term and the vacuum potential plot were guarded off; the energy term is back | |
| `Libraries/build_pspline.sh` | PSPLINE `r8genxpkg` wrote past its array for three points, the macOS "MPI flake"; one-line guard applied | |
| CMake | gfortran 14.2 and 14.3 miscompile `-fcheck=bounds` on one pattern; the check is dropped for those versions | |
| Build system | | STRUMPACK optional with a ScaLAPACK fallback; `PB3D_BUILD_EXECUTABLES` for a dependency-light core |

## 3. Helper material outside the code

| Kind | What |
|---|---|
| Tests, unit | `dtorh` against 30-digit references up to degree 80, string and numerical utilities, mode-index and Richardson arithmetic |
| Tests, full-stack | vacuum kernels, Green identities, 3-D vacuum (`vac%style` 1) on an analytic torus against the axisymmetric one (`vac%style` 2) (`response_1_vs_2`, recorded ratios), `GH_int_1_near`, `exact_exterior_2d` (exact toroidal harmonics through `vac%style` 2), `point_source_3d`, splines, volume integral, HELENA golden file |
| Tests, regression | cbm18a fixed (one, two and four processes), free with `BC_style` 2 and 4 and the response check, HELENA shear, circular ε = 0.05 kink against MISHKA, circular ε = 0.33 HELENA against its VMEC twin fixed and free through the 3-D vacuum; anchors reproducible to 1e−12 on macOS, Ubuntu and Snellius |
| CI | eight GitHub Actions jobs: lint, unit Debug and Release, full stack, full stack Debug, sanitiser, docs, coverage; compiler-warning and Doxygen ratchets |
| Fixtures | cbm18a (xz), circular ε = 0.05 and 0.33 HELENA cases with their inputs, the VMEC twin with its recipe (`make_vmec_twin.sh`, `vmec_twin.py`), MISHKA and cylinder references |
| Equilibria made here | the circular ε = 0.33 and 0.05 HELENA cases, their VMEC twins at 159 to 1279 surfaces, the rippled twins, the cbm18a plasma cut at the pedestal foot with HELENA (`tests/reference/cbm18a_cut/`, scaled from the extended input without fitting) and its VMEC twin |
| Scripts | `scripts/snellius/` (build, `run_deck.sbatch`, response and check comparisons), `scripts/terpsichore/` (netCDF to `fort.18` with `--flip` and `--nfp`; card generator with the parallel current from the geometry, the pressure gradient from the file, a poloidal offset and a sized toroidal table; recipe and traps), `Documentation/figures/step0/make_figures.py` |
| Documentation | the plan and log, `benchmarking_prep.md` with three literature surveys, `testing.md`, `quality-and-ci.md`, `snellius-setup.md`, `spack-setup.md`, the openwiki pages |
| External codes set up | MISHKA-1 rebuilt on arm64 with ILP64 BLAS; HELENA `hel12a`; TERPSICHORE 1.2 built on Snellius, with a corrected copy for its kinetic-energy parity error; VMEC++ for the twins; a Snellius stack with PETSc, SLEPc, MUMPS, PSPLINE and LIBSTELL under `~/pb3d` |

## 4. Agreement with other codes

ω² in PB3D's Alfvén units, n = 10 unless stated.

| Case | PB3D | Other code | Distance |
|---|---|---|---|
| Circular ε = 0.33, fixed boundary (HELENA and VMEC twin agree to 2e−6) | +0.07843 | TERPSICHORE +0.07851 | 0.1 % |
| Circular ε = 0.33, external kink, 3-D vacuum at 128 lines × 400 points | −0.43849 | TERPSICHORE −0.43848 (converged over 159 to 639 surfaces) | 2e−5 |
| same, axisymmetric vacuum, HELENA mapping with 257 → 513 → 1025 points | −0.43679 → −0.43731 → −0.43769 since the corner-term fix of 2 November (−0.43380 → −0.43582 → −0.43695 before) | TERPSICHORE −0.43848 | 0.39 % at 257 points, 0.18 % at 1025 (1.07 % and 0.35 % before) |
| same | about −0.438 (PB3D and TERPSICHORE) | MISHKA-1 0.43377, converged in its grids | 1.1 %, probably the parallel inertia that MISHKA-1 keeps and the other two drop; not verified |
| Circular ε = 0.05 kink, n = 10, 15, 20 | | MISHKA-1 | 0.53, 0.24, 0.13 % since 2 November (1.2, 0.55, 0.26 % before), falling as 1/n² (the high-n ordering) |
| Rippled twin δ = 0.04, five periods, external kink, 3-D vacuum | −0.42100 at 128 × 200; −0.41822 at 128 × 400 | TERPSICHORE n = 0..20 family −0.41626 in PB3D's units | 1.1 % at 128 × 200, 0.47 % at 128 × 400; the unrippled twin is at 2e−5 at that resolution, so half a percent sits in the treatment of the ripple itself (PB3D's field-line average of the sideband coupling, or TERPSICHORE's truncation) |
| Rippled twin, sideband coupling | carried by PB3D's field-line representation | TERPSICHORE single row n = 10 against the family: 1.6 % | consistent |
| cbm18a cut at the pedestal foot, n = 10, 15, 20, fixed boundary | −0.04012, −0.08128, −0.11120 | TERPSICHORE −0.03991, −0.07951, −0.10701 | 0.5, 2, 4 % |
| same, free boundary | −0.05002, −0.08577, −0.11209 | TERPSICHORE (wall at 2a) −0.04979, −0.08441, −0.10753 | 0.5, 1.6, 4 % |
| same, both boundaries | | MISHKA-1 0.03637, 0.07009, 0.09290 fixed; 0.04492, 0.07407, 0.09453 free | MISHKA-1 is 9 to 16 % less unstable than PB3D and TERPSICHORE; cause not identified |
| W7-X standard case, n = 10 | stable, dense spectrum | TERPSICHORE, no eigenvalue below −1e−4 | no discrete mode to compare at this pressure |
| TERPSICHORE's own benchmark | | reproduced to the last digit | its release has a parity error in the physical kinetic-energy model; corrected in a source copy, report drafted |

## 5. Tests that pass, tests that do not

Pass:

- CI eight of eight jobs; Release CTest 27 of 27 on macOS and Snellius, Debug suite on both.
- Every regression anchor to 1e−3 and reproducible to 1e−12 across three machines; one to four processes agree to 1.5e−8.
- Exact-solution checks on the real surfaces: Green identities 1e−5 to 1e−4 and exterior solves 3e−4 to 4e−3 at 128 lines × 200 points, second order or better in the points, on the circular twin and on W7-X.
- Invariances: rigid rotation of a 3-D boundary written as a non-symmetric VMEC file leaves the response unchanged to 1e−16; scaling of R and B leaves the dimensionless eigenvalues unchanged to 1e−13; n → the ripple amplitude enters only at second order.
- Convergence: second order in the radial grid for the circular cases (0.03 % at 200 points); first order in field lines and converged by 128; the 3-D response reaches the axisymmetric one as the parallel points per poloidal wavelength grow (0.981, 0.988, 0.996 at 9, 18, 36).

Do not pass, or are not validated:

- `BC_style` 4 is 7 % off `BC_style` 2 on the circular kink and first order in the radial grid; kept as a change detector only.
- The axisymmetric vacuum, `vac%style` 2, was first order in HELENA's poloidal points (1 % at the fixtures' 257 and 401 points). Most of that was the corner term fixed on 2 November; the circular kink is now 0.39 % from TERPSICHORE at 257 points and 0.18 % at 1025, but still converges slower than second order, the interior identity (not used) got worse by as much as the exterior one improved, and the sequence still stops being monotone at 2049 points. The free-boundary anchors therefore pin reproducibility, not converged values.
- cbm18a converges only at first order in the radial grid (its 22 rational surfaces); the 200-point anchors are 0.5 % from the limit, and the cut case moves 1 % between 200 and 400 points at n = 20.
- The mutation test for the 2π/n_alpha spacing is inside the tolerances of the 3-D tests; the full-stack tests check orders, not values, there.
- Sixteen MPI processes hang in the parallel HDF5 write on the W7-X decks; four work.
- The Debug build on gfortran 14 runs without bounds checking.
- The toroidal flux label solves a different problem (`prim_X` changes meaning) and is untested; negative n is rejected at input.
- POST's energy reconstruction is 7 % off on the circular edge mode (a plot-grid problem, not pursued).

## 6. Where uncertainty remains

- MISHKA-1 sits 1.1 % less unstable than PB3D and TERPSICHORE on the converged circular kink. The inertia model is the likely cause and has not been checked.
- MISHKA-1 is 9 to 16 % less unstable than PB3D and TERPSICHORE on the cbm18a pedestal mode (the plasma cut at the foot, `tests/reference/cbm18a_cut/`), fixed boundary included, and 1.1 % on the circular kink. Two codes agree against one; the cause in MISHKA-1 (inertia model, its equation set) is not identified, and ELITE would be the natural third opinion on the pedestal mode.
- The 3-D vacuum is expensive at the accuracy shown: 128 lines × 400 parallel points is 51 000 boundary points, 220 GB and 39 minutes on eight cores for 0.1 %; at 64 lines × 100 points expect about 1 %. Which part of the method sets this cost is not settled: the projection of the modes on the field-line grid and the solve for oscillating data both contribute (section 11.9 of the report). Options to explore are in the work bucket.
- W7-X has no benchmark yet: the standard case is stable at n = 10 in both codes, and the W7-X vacuum response converges only when lines and points are refined together (61 points per line are too few). A higher-pressure W7-X is the way to a 3-D stellarator comparison.
- The fast mode window (`X_style` 2) is not a nested trial space (6e−4 on cbm18a); double `n_mod_X` once when it matters.
- The TERPSICHORE correction rests on the parity argument and two cases, and its default parallel current from charge conservation loses a pedestal's edge current (one eighth of the n = 10 eigenvalue of the cut cbm18a, about a third of the growth rate); both belong in the report to EPFL, which is Toon's to send, as are the GCC bug report and the merge of PR #5.

## 7. Work bucket

- POST plotting into the vacuum: check `plot_vac_pot` on the 3-D vacuum (`vac%style` 1), draw the perturbation and the potential across the interface, add the vacuum energy density |∇φ|²/2μ₀ to the energy reconstruction.
- **Exploration only, not decided and not required:** ways to make the 3-D vacuum cheaper. None replaces the present method unless it beats it on the `check_vac` exact potentials (circular twin and W7-X) and on the twin's eigenvalue; Toon decides whether any is tried.
  - A doubly periodic vacuum grid in the straight-field-line angles (θ, ζ), separate from the plasma's field-line grid, with `res` still handed over in mode space. The (α, θ) field-line grid is not doubly periodic (one poloidal turn at fixed α lands on another line), which is where the seam and bridge come from.
  - On that grid, a locally corrected trapezoidal rule: the trapezoidal sum everywhere, with new weights on a small stencil around each singular point (Wu–Martinsson zeta corrections, Alpert or Duan–Rokhlin corrections, Malhotra's polar patches as in BIEST and GPEC's Julia `Vacuum`). The literature gives 4th to 10th order; not tried here.
  - A vacuum in a Fourier basis, as in CAS3D (Merkel, Nührenberg and Cooper, "Free-boundary ideal MHD modes in W7-AS", CRPP Lausanne report LRP 555/96, pp. 51–58 (1996); the NESTOR method of Merkel 1986): Φ solved as Fourier coefficients on one field period with stellarator symmetry, the singular kernels regularised analytically. Same physics and edge-block coupling as now.
  - A solver whose memory and work grow like N or N log N in the number of boundary points N, instead of the dense LU's N² memory and N³ work: GMRES with a fast-multipole or hierarchical-matrix product (the operator's condition number is 3.5 to 3.8 on every grid), or a HODLR or H² direct solve with the points ordered geometrically. STRUMPACK's compression may have failed because the field-line ordering scatters neighbouring points and its sampling parameters were hard-coded; not checked.
  - Richardson would either refine the vacuum grid with each level or compute the vacuum once, accurately enough for every level.
  - Guesses, not measurements: the corrected quadrature might need ten times fewer points, and with a fast solver 220 GB might become a few GB.
- A Claude skill in the repository that runs PB3D for other users: equilibrium in, deck choices made with the convergence rules above, run locally or through Slurm, eigenvalue and POST plots out.
- From the benchmarking plan: GPEC's 3-D VACUUM (C6) or CASTOR3D (C7) on the rippled case as a fourth vacuum method; DESC's singular integrals on the kernels (C9); COBRAVMEC at large n (L11); ELITE on the cut cbm18a as the third opinion against MISHKA-1; a higher-pressure W7-X; TERPSICHORE's radial convergence at n = 20 on the cut case (its input stops at 997 surfaces).
