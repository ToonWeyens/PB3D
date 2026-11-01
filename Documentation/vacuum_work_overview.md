# PB3D vacuum work: overview

Branch `dev/vacuum-fixes`, draft PR #5 to `master`, tracking issue #4. Seventy commits in two stretches: the modernisation of July 2026 and the vacuum work of 27 October to 1 November 2026. The detailed log is `Documentation/plan_vacuum_bugfix.md`; this page is the outcome only.

## 1. Where we started, where we are

| | Before (PB3D 2.47, 2018 to 2020 code) | Now (2.49, this branch) |
|---|---|---|
| Build | PETSc 3.9 API, no longer buildable on a current machine | PETSc/SLEPc 3.25 on macOS arm64, Snellius and GitHub Actions |
| Tests | none that ran | unit, full-stack and regression layers under CTest; eight CI jobs, all green |
| Axisymmetric vacuum (HELENA, style 2) | solved the interior problem; result never checked | exterior problem, Hermitian response; free-boundary eigenvalues within 0.3 % of TERPSICHORE when the HELENA mapping is fine enough |
| 3-D vacuum (VMEC, field-line style 1) | code present, never gave a right answer | converges; reproduces TERPSICHORE's free-boundary eigenvalue to 2e−5 on the circular twin and runs on a rippled torus and on W7-X |
| Fixed boundary | the one result that was trusted, but 6.5 % too stable through the HELENA shear bug and 2e−5 process-count dependent | HELENA and VMEC paths agree to 2e−6, process counts to 1.5e−8, TERPSICHORE to 0.1 % |

## 2. Changes to the PB3D and POST code

Fixed means the code was wrong and now does what it was meant to do; new means the capability did not exist.

| Where | Fixed | New |
|---|---|---|
| `SLEPC_ops`, `SLEPC_utilities` | port from the PETSc 3.9 API to 3.25 (July) | |
| `vac_ops` `solve_Phi_BEM` | the axisymmetric style solved the interior system; it now solves the exterior one, (H + 4πI)Φ = G R (July) | |
| `vac_ops` `calc_GH_int_1` | singular interval integral, quadrature weights, sign and a `**(−3)` far-field factor in the 3-D kernels (July) | near-field integration along the field line, Gauss–Legendre on the cubic interpolant through four points (31 October) |
| `vac_ops` `calc_GH_1` | α spacing 2π/n_alpha instead of 2π/(n_alpha − 1); G and H no longer counted twice across Richardson levels | G diagonal from the identity for the linear interior harmonics, as the row sum sets the H diagonal |
| `vac_ops` `store_vac_VMEC` | handed the kernels the outward normal where they take the inward one; the 3-D response was 3 to 4 times too small | |
| `vac_ops` | | `check_vac_exact` and the input `check_vac`: Green identities and exact exterior solves on the real surface, per-point files, condition estimate |
| `SLEPC_ops` energy assembly | the vacuum entered with the wrong sign and without the factor 2π, so it destabilised; now −2π X_e^H res X_e with the Hermitian part of the response | |
| `SLEPC_ops` `fill_mat` | the free end point carried the full radial weight h instead of h/2; free-boundary eigenvalues converged at first order in `n_r_sol` | |
| `SLEPC_utilities` boundary conditions | style 4 overwrote the vacuum block; style 3 removed; default back to `BC_style = 1 2` | style 2, the variational (Hermitian) free-boundary condition (July) |
| `HELENA_ops` | the shear used the even spline condition on the odd Θ^α since June 2018; cbm18a was 6.5 % too stable | |
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
| Tests, full-stack | vacuum kernels, Green identities, 3-D vacuum on an analytic torus against the axisymmetric style (`response_1_vs_2`, recorded ratios), `GH_int_1_near`, `exact_exterior_2d` (exact toroidal harmonics through style 2), `point_source_3d`, splines, volume integral, HELENA golden file |
| Tests, regression | cbm18a fixed (one, two and four processes), free with styles 2 and 4 and the response check, HELENA shear, circular ε = 0.05 kink against MISHKA, circular ε = 0.33 HELENA against its VMEC twin fixed and free through the 3-D vacuum; anchors reproducible to 1e−12 on macOS, Ubuntu and Snellius |
| CI | eight GitHub Actions jobs: lint, unit Debug and Release, full stack, full stack Debug, sanitiser, docs, coverage; compiler-warning and Doxygen ratchets |
| Fixtures | cbm18a (xz), circular ε = 0.05 and 0.33 HELENA cases with their inputs, the VMEC twin with its recipe (`make_vmec_twin.sh`, `vmec_twin.py`), MISHKA and cylinder references |
| Equilibria made here | the circular ε = 0.33 and 0.05 HELENA cases, their VMEC twins at 159 to 1279 surfaces, the rippled twins, the cbm18a plasma cut at the pedestal foot with HELENA (`tests/reference/cbm18a_cut/`, scaled from the extended input without fitting) and its VMEC twin |
| Scripts | `scripts/snellius/` (build, `run_deck.sbatch`, response and check comparisons), `scripts/terpsichore/` (netCDF to `fort.18` with `--flip` and `--nfp`, card generator, recipe), `Documentation/figures/step0/make_figures.py` |
| Documentation | the plan and log, `benchmarking_prep.md` with three literature surveys, `testing.md`, `quality-and-ci.md`, `snellius-setup.md`, `spack-setup.md`, the openwiki pages |
| External codes set up | MISHKA-1 rebuilt on arm64 with ILP64 BLAS; HELENA `hel12a`; TERPSICHORE 1.2 built on Snellius, with a corrected copy for its kinetic-energy parity error; VMEC++ for the twins; a Snellius stack with PETSc, SLEPc, MUMPS, PSPLINE and LIBSTELL under `~/pb3d` |

## 4. Agreement with other codes

ω² in PB3D's Alfvén units, n = 10 unless stated.

| Case | PB3D | Other code | Distance |
|---|---|---|---|
| Circular ε = 0.33, fixed boundary (HELENA and VMEC twin agree to 2e−6) | +0.07843 | TERPSICHORE +0.07851 | 0.1 % |
| Circular ε = 0.33, external kink, 3-D vacuum at 128 lines × 400 points | −0.43849 | TERPSICHORE −0.43848 (converged over 159 to 639 surfaces) | 2e−5 |
| same, axisymmetric vacuum, HELENA mapping with 257 → 513 → 1025 points | −0.43380 → −0.43582 → −0.43695, limit near −0.438 | TERPSICHORE −0.43848 | about 0.1 to 0.3 % after extrapolation |
| same | about −0.438 (PB3D and TERPSICHORE) | MISHKA-1 0.43377, converged in its grids | 1.1 %, probably the parallel inertia that MISHKA-1 keeps and the other two drop; not verified |
| Circular ε = 0.05 kink, n = 10, 15, 20 | | MISHKA-1 | 1.2, 0.55, 0.26 %, falling as 1/n² (the high-n ordering) |
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

- Boundary-condition style 4 is 7 % off style 2 on the circular kink and first order in the radial grid; kept as a change detector only.
- The axisymmetric style 2 is first order in HELENA's poloidal points (1 % at the fixtures' 257 and 401 points), and its sequence stops being monotone at 2049 points; the free-boundary anchors therefore pin reproducibility, not converged values.
- cbm18a converges only at first order in the radial grid (its 22 rational surfaces); the 200-point anchors are 0.5 % from the limit.
- The mutation test for the 2π/n_alpha spacing is inside the tolerances of the 3-D tests; the full-stack tests check orders, not values, there.
- Sixteen MPI processes hang in the parallel HDF5 write on the W7-X decks; four work.
- The Debug build on gfortran 14 runs without bounds checking.
- The toroidal flux label solves a different problem (`prim_X` changes meaning) and is untested; negative n is rejected at input.
- POST's energy reconstruction is 7 % off on the circular edge mode (a plot-grid problem, not pursued).

## 6. Where uncertainty remains

- MISHKA-1 sits 1.1 % less unstable than PB3D and TERPSICHORE on the converged circular kink. The inertia model is the likely cause and has not been checked.
- MISHKA-1 is 9 to 16 % less unstable than PB3D and TERPSICHORE on the cbm18a pedestal mode (the plasma cut at the foot, `tests/reference/cbm18a_cut/`), fixed boundary included, and 1.1 % on the circular kink. Two codes agree against one; the cause in MISHKA-1 (inertia model, its equation set) is not identified, and ELITE would be the natural third opinion on the pedestal mode.
- The 3-D vacuum is expensive at the accuracy shown: 128 lines × 400 parallel points is 51 000 boundary points, 220 GB and 39 minutes on eight cores for 0.1 %; at 64 lines × 100 points expect about 1 %. A higher-order projection of the modes along the lines would cut this.
- W7-X has no benchmark yet: the standard case is stable at n = 10 in both codes, and the W7-X vacuum response converges only when lines and points are refined together (61 points per line are too few). A higher-pressure W7-X is the way to a 3-D stellarator comparison.
- The fast mode window (`X_style` 2) is not a nested trial space (6e−4 on cbm18a); double `n_mod_X` once when it matters.
- The TERPSICHORE correction rests on the parity argument and two cases; the report to EPFL is Toon's to send, as are the GCC bug report and the merge of PR #5.

## 7. Work bucket

- POST plotting into the vacuum: check `plot_vac_pot` on the 3-D style, draw the perturbation and the potential across the interface, add the vacuum energy density |∇φ|²/2μ₀ to the energy reconstruction.
- A Claude skill in the repository that runs PB3D for other users: equilibrium in, deck choices made with the convergence rules above, run locally or through Slurm, eigenvalue and POST plots out.
- From the benchmarking plan: GPEC's 3-D VACUUM (C6) or CASTOR3D (C7) on the rippled case as a fourth vacuum method; DESC's singular integrals on the kernels (C9); COBRAVMEC at large n (L11); the cbm18a cut equilibrium; a higher-pressure W7-X.
