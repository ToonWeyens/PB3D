# PB3D vacuum work: detailed report

Branch `dev/vacuum-fixes`, draft PR #5 to `master`, tracking issue #4. Written on 1 November 2026 from the sources listed in section 1.4.

## 1. Purpose, scope, how to read this report, and its sources

### 1.1 Purpose

This report is the detailed companion of `Documentation/vacuum_work_overview.md`. The overview states the outcome on one page. This report gives the evidence behind each line of the overview, ordered by topic.

The chronological log of the same work is `Documentation/plan_vacuum_bugfix.md`. This report cites its dated blocks but does not follow their order.

### 1.2 Scope

The branch holds 73 commits in two stretches. The first stretch, 6 to 10 July 2026, ported PB3D to the current PETSc and SLEPc, built a test suite and CI, and made a first repair of the vacuum code. The second stretch, 27 October to 1 November 2026, followed the plan of issue #4.

That plan fixed the known bugs one at a time. It then validated the vacuum against exact solutions and against other codes, and measured what the 3-D vacuum costs. Toon has not merged the branch; PR #5 is still a draft.

### 1.3 How to read it

- ω² means ω²/ω_A² in PB3D's normalisation (section 3.3). A negative value is unstable.
- MISHKA-1 prints γ² = −ω²/ω_A², so its unstable values are positive. Tables that hold both say so.
- A date in brackets marks when a result was obtained. Dates without a year are in 2026.
- "Fixed" in section 4 means the code was wrong and now does what it was meant to do. "New" in section 5 means the capability did not exist. This is the overview's split.
- Two different "styles" appear. `vac%style` is the vacuum method, set by the equilibrium code, not an input: `vac%style` 1 is the 3-D boundary-element method on field lines (VMEC), `vac%style` 2 the axisymmetric one with toroidal harmonics (HELENA). `BC_style` is the input for the radial boundary conditions; this report always means its second value, the plasma edge: 1 fixed boundary, 2 the free boundary through the variational surface energy, 4 the free boundary through an explicit natural-boundary-condition row. Other inputs with "style" in their name (`alpha_style`, `X_style`, `norm_style` and so on) are always named in full (Appendix C).
- A statement that no source supports is marked "[not in the sources]".
- Paths are relative to the repository root unless they start with `~`. Paths with `~` are on Toon's Mac or, where stated, on Snellius under `~/pb3d`.

### 1.4 Sources

| Rank | Source | What this report takes from it |
|---|---|---|
| 1 | `Documentation/vacuum_work_overview.md` | the verdicts and the summary numbers |
| 2 | `Documentation/plan_vacuum_bugfix.md` | steps 0 to 9, the dated blocks of 27 October to 1 November, the validation catalogue of its section 4, the open points of its section 9 |
| 3 | `README.md`, sections 2.49 and 2.48 | the changelog |
| 4 | `Documentation/testing.md`, `Documentation/quality-and-ci.md` | the test layers, the anchors, the CI |
| 5 | `Documentation/benchmarking_prep.md` and `Documentation/benchmarking/lit_*.md` | the literature, the normalisations, the priority list |
| 6 | `scripts/terpsichore/README.md`, `tests/reference/cbm18a_cut/README.md` | the TERPSICHORE recipe and traps; the cut cbm18a |
| 7 | `Documentation/snellius-setup.md` | the Snellius stack and the job rules |
| 8 | `tests/regression/*.sh` | what each anchor is and pins |
| 9 | `git log`, `git diff --stat master...HEAD`, the commit messages | the file map and the commit list |
| 10 | the code, read where a section names what a routine does | `vac_ops.f90`, `vac_utilities.f90`, `SLEPC_ops.f90`, `VMEC_ops.f90`, `MPI_utilities.f90`, `grid_ops.f90`, `X_ops.f90`, `eq_ops.f90`, `input_ops.f90`, `CMakeLists.txt` |

When two sources disagree, the later date wins, and the overview wins over all. The text says where this rule was used.

The known-bugs report of 8 October 2026 is not in this list. It lives in `~/Code/MHD`, branch `claude/pb3d-report`, file `pb3d_code_analysis.md`. This report cites it only through the plan.

## 2. Starting point

### 2.1 PB3D 2.47

`master` is at 5c2fe76, version 2.47, January 2026. Version 2.47 moved the build from a Makefile to CMake, added a Spack environment and kept the old Makefile as `Makefile.legacy`.

The physics code was that of 2018 to 2020. It used the PETSc 3.9 API, which no longer builds on a current machine. No test ran.

### 2.2 PB3D 2.48 and the July stretch

Toon's 2.48 work came first, 6 to 8 July, in nine commits from a8c07b6 to 78fb49e. The 2.48 changelog lists the port to PETSc 3.25.3 and SLEPc 3.25.1, confined to `SLEPC_ops` and `SLEPC_utilities`. Fixed-boundary runs skip the vacuum cleanly, so HELENA runs work for the first time since 2.44. The cbm18a fixed-boundary eigenvalue became a CTest anchor. Build scripts for PSPLINE and a minimal LIBSTELL, a CTest suite in the vendored `test-drive` framework, the option `PB3D_BUILD_EXECUTABLES`, `dtorh` tests against 30-digit mpmath values, a GitHub Actions workflow and the committed `spack.yaml` complete the list.

The 2.48 changelog ends with "Free-boundary vacuum remains under construction".

A cloud session then added 22 commits, 8 to 10 July, from 1862df0 to 17b712c. The plan counts 21; the git log shows 22 commits by that author. By topic:

| Topic | Commits | What they did |
|---|---|---|
| build | 1862df0 | guards for PETSc 3.19 and later; fixes for MPI in containers |
| vacuum | 429f003, 3679ee7, d6d49e6 | exterior operator, 3-D kernel repairs, negative n in the axisymmetric kernels, optional STRUMPACK, HELENA vacuum storage |
| boundary conditions | cc4c96d, 41b5274, 1d74b8f | `alpha_style` 2 required for 3-D free runs; `BC_style` 3 removed and the default set to `BC_style = [1,4]`; `BC_style` 2 implemented |
| tests | 31fe4a7, 89996a1, 9f9a69c, 0f5b27e, 9db1d9e, fa3560b, 1c7283c, 88153e4, 17b712c | full-stack vacuum layer, cylinder test, free-boundary regression, committed cbm18a fixture, multi-process suites, converted legacy checks, HELENA golden file, index-arithmetic units |
| quality | c16015c, 8ac37d7, 6f0f8d5, 6ce5029, 17716c4, a53869c | CI on the full stack, PSPLINE split from `num_utilities`, six-job pipeline, warnings from 321 to 124, fortitude OB gate, sanitiser job, Doxygen ratchet |

On 27 October the tip 17b712c was the newest PB3D code anywhere. It sat at the top of three stacked draft PRs, numbered 1 to 3. Their descriptions made three claims that the October work found wrong (plan, step 7):

- "4 % more unstable than fixed" measured a longer domain, not the vacuum;
- "style 2 agrees with style 4" (`BC_style` 2 and `BC_style` 4) held only where the edge displacement is zero;
- "3-D verified end to end" referred to the synthetic torus of the test suite.

### 2.3 The known-bugs report of October 2026

The report `pb3d_code_analysis.md` of 8 October was written against 17b712c. The plan cites its items by number.

| Report item | Subject | Plan step |
|---|---|---|
| 1 | crashes on macOS | 1 |
| 2 | HELENA shear | 2 |
| 3 | vacuum sign and factor | 3 |
| 4 | `BC_style` 4 drops the vacuum; the default `BC_style` | 4 |
| 5 | `setup_modes` slots | 5 |
| 9 | 3-D field-line spacing α | 6b |
| 10 | Richardson double count of G and H | 6c |
| 11 | aliasing of the field-line average | 6d |
| 12 | CI red | 1 |
| 13 | orientation of the 3-D normal, marked "inferred" | 6b; found in 6e |

Report section 9 noted that the ψ_N = 1.0 decks cannot see the vacuum. Report section 12 noted that the input path for negative n is broken. One claim of the report was wrong: not every VMEC case has nq < 0 (step 5).

### 2.4 The numbers at the start

Step 0 (27 October) ran every deck on 17b712c, unchanged. The runs used a Release build on the Mac: gfortran 16.1, PETSc 3.25.3 with complex scalars, SLEPc 3.25.1, singleton MPI. The deck is cbm18a at n = 10 with 20 modes, 200 radial and 100 parallel points, Simpson 3/8 and `tol_SLEPC` 1e−6.

| Run | ω² at 17b712c |
|---|---|
| fixed, ψ_N 0.45 to 0.75 | −4.0735999821e−2 |
| free, `BC_style` 2, to 0.75 | −4.3467138270e−2 |
| free, `BC_style` 4, to 0.75 | −4.5574082987e−2 |
| fixed, 0.45 to 1.0 | −4.2348778332e−2 |
| free, `BC_style` 2, to 1.0 | −4.2348778421e−2 |
| free, `BC_style` 4, to 1.0 | −4.2348776924e−2 |
| vacuum × 0 at 0.75, diagnostic patch | −4.3142687497e−2 |
| vacuum × (−2π) at 0.75 | −4.2150080544e−2 |

Every number of the bugs report reproduced to five digits. The three ψ_N = 1.0 runs agreed to 4e−8, so that deck could not see the vacuum. Nine of fifteen launches of the regression decks crashed before they wrote the eigenvalue, and the 800-point radial deck crashed eight times of eight.

### 2.5 The plan's method

The plan of 27 October lists eight ground rules (plan, section 0):

1. Make the smallest change that fixes the bug, with a one-line comment at the place.
2. One bug, one commit, one test. The test is written first and must fail on the unfixed code with a specific wrong number.
3. State the expected number before the run. A surprise stops the step.
4. Infrastructure steps leave every eigenvalue unchanged. Physics steps change them by the predicted amount only.
5. Regression anchors hold to 1e−3 relative. A crashed run is a failure, with no retries.
6. Every step passes CI before the next step starts.
7. Some topics are out of scope (section 12.3).
8. Local checks come before any request to another group.

The plan adds two working rules. Every commit references issue #4. Only Toon merges: the work stops at a pushed branch and a draft PR.

The plan went through six critique rounds (plan, section 7). Round 4 turned the checks into a catalogue of 38 items in four kinds: theory (T), logical conditions (L), other codes (C) and other means (O). Rounds 5 and 6 added the public vacuum codes and rule 8.

| Step | Subject | Gate | Outcome |
|---|---|---|---|
| 0 | references: baseline, circular cases, MISHKA-1, cylinder formula | the bugs report reproduces | passed (27 October) |
| 1 | crashes and CI | 20 of 20 clean runs, eigenvalues unchanged | passed (27 October) |
| 2 | HELENA shear | predicted value reached, shear test passes | passed (27 October) |
| 3 | vacuum sign and factor | no vacuum ≤ free ≤ fixed; MISHKA-1 and cylinder | passed (27 October) |
| 4 | `BC_style` 4 and the default | `BC_style` 2 equals `BC_style` 4 to 1e−5 | failed (27 October); `BC_style` 4 kept as change detector |
| 5 | `setup_modes` | no wrong slot in the simulation | passed (27 October) |
| 6a to 6d | 3-D review, α spacing, Richardson, aliasing | point-source test; level-2 equals fresh run | passed (27 October) |
| 6e | circular VMEC twin through the 3-D vacuum | twin equals HELENA through `vac%style` 2 | failed, then passed after the normal fix (30 October) |
| 6f | rippled twin | ω²(δ) − ω²(0) ∝ δ² | passed after two fixes (30 October) |
| 7 | mutation check, changelog, PR | each fix caught by its test | three of four caught (30 October) |
| 8 | 3-D equilibrium path against analytic equilibria | DESC-written wout of a Landreman case | not started |
| 9 | trigonometric factors in basis form | anchors unchanged | passed (30 October) |

The plan budgeted under 40 physics lines for steps 1 to 6 (plan, section 6). The external requests of step 7 were not sent. TERPSICHORE, open source since 2025, needed no request and was run on Snellius (section 8.6).

## 3. Reference cases, codes and normalisations

### 3.1 Equilibria

| Case | File | Made with | Defining parameters | Used for |
|---|---|---|---|---|
| cbm18a extended | `tests/fixtures/cbm18a_extended_vac.12.xz` | PB3D-patched HELENA on the original cbm18a input | pedestal ends at ψ_N = 0.70; extension without pressure or current to ψ_N = 1; q = 1.05 on axis, 3.02 at the plasma edge, 25 at ψ_N = 1 | regression decks; MISHKA-1 |
| cbm18a cut | `tests/reference/cbm18a_cut/` | hel12a from the extended input, ending at the foot | EPS = 0.344262; q from 1.052 to 3.040 | pedestal benchmark |
| circular ε = 0.33 | `tests/fixtures/circ_eps033_q104.12.xz` | hel12a on the PhD input `TEST_VAC/2.25` | flat Γ, B = 0.001, XIAB = 2.25, q_0 = 1.0411, q_a = 1.0382 | external kink; source of the twins |
| circular ε = 0.05 | `tests/fixtures/circ_eps005_q104.12.xz` | hel12a | XIAB = 0.3035, Γ = 1 − 0.02ψ, q_0 = 1.0274, q_a = 1.0378 | cylinder formula; MISHKA-1 |
| circular VMEC twin | `tests/fixtures/wout_circ_eps033_q104.nc.xz` | `--export_HEL`, then VMEC++ | as ε = 0.33; R_0 = B_0 = 1 | HELENA against VMEC; 3-D vacuum at ε = 0 |
| rippled twins | Snellius fixtures, e.g. `wout_circ_eps033_rip04.nc` | VMEC++ with NTOR = 6 | boundary a(1 + δ cos 5ζ), δ = 0.01, 0.02, 0.04 | δ² law; TERPSICHORE |
| W7-X | `wout_w7x.nc`, `wout_w7xflip.nc` | STELLOPT file; copy with the field reversed | N_fp = 5; 288 VMEC modes; q = 1.1313 at the edge | smoke test; response convergence; TERPSICHORE |

**cbm18a extended.** The fixture is HELENA's extended pseudo-vacuum equilibrium. The pressure gradient vanishes at ψ_N = 0.70. From there to the HELENA boundary at ψ_N = 1 the domain carries neither pressure gradient nor current. Lengths are in units of R_mag.

The plasma edge has half-width 0.333 R_mag about R = 0.969 R_mag. The HELENA boundary has half-width 0.555 R_mag about R = 0.832 R_mag. The magnetic axis sits 0.17 R_mag outside the centre of the HELENA boundary.

| ψ_N | p/p(0) | q | shear s |
|---|---|---|---|
| 0 | 1 | 1.05 | |
| 0.45 | 0.92 | | 1.5 |
| 0.60 | 0.33 | | flat stretch |
| 0.65 | 0.03 | | |
| 0.70, plasma edge | 0 | 3.02 | 5.4 |
| 0.75, deck boundary | 0 | 3.71 | 6.4 |
| 1.0, HELENA boundary | 0 | 25 | |

The regression decks solve from ψ_N = 0.45 to 0.75. Their boundary lies 7 % of the plasma minor radius outside the plasma edge, in the pseudo-vacuum. The rational surfaces q = m/10 with m = 15 to 36 lie in this window, 22 of them. The n = 10 mode peaks at ψ_N = 0.58 with m = 20 and 21 and balloons on the outboard side. The fixture is 26 MB, 5.2 MB with xz, and has NCHI = 401 poloidal points. MISHKA-1 needs a 513-point mapping (section 3.2).

**cbm18a cut at the pedestal foot (1 November).** The extended input's p′ and FF′ tables were truncated at the foot. HELENA's tables live on a uniform grid in √ψ_N, so the foot at table index 217 of 259 is ψ_N = (216/258)² = 0.7009 of the extended flux. The foot flux surface of the extended mapping became the boundary, as 24 cosine coefficients with `ISHAPE = 2`.

The amplitude parameters were scaled, not fitted. `B` scales by (R0_new/R0_old)², `XIAB` by (a_old R0_new)/(a_new R0_old). The cut plasma then matches the extended one:

| Quantity | Agreement with the extended run |
|---|---|
| q on axis and at the foot | 6e−6 |
| q profile inside | 3e−4 |
| p profile inside | 7e−6 |
| axis position | 3e−7 |

A first attempt read the foot off a uniform-ψ grid. It needed a fitted B 50 % off the scaling and put the axis 1.5 % out of place, which is how the √ψ_N grid was found. The HELENA run uses NR 251, NP 201, NRMAP 401, NPMAP 513 and NCHI 513. The mappings `cbm18a_cut.12` and `cbm18a_cut_nchi1025.12` sit in Snellius' fixtures. The VMEC twin has 639 surfaces and reproduces q to five digits. The pedestal of the cut plasma lies at ψ_N 0.75 to 1.

**Circular cases.** The ε = 0.33 case is the PhD case `Simulations/TEST_VAC/2.25`, regenerated with the rebuilt HELENA `~/Code/helena/hel12a`. Its current is flat and its pressure zero. The lowest rational surfaces of the decks are q = 11/10, 16/15 and 21/20, all outside the plasma. No mode resonates inside, so the fixed boundary must be stable, and every instability is an external kink: m = 11, 16, 21 for n = 10, 15, 20.

PB3D rejects a perfectly flat q, q_0 = q_a = 1.0378, with "The code has to be adapted for reversed-shear regions still" (`calc_eq_1`, error 3). The committed ε = 0.05 variant therefore has Γ = 1 − 0.02ψ, set by `AGA = −0.02`. Both fixtures come with their HELENA inputs `*.fort.10` and have 257 poloidal points. On 1 November the ε = 0.33 case was remapped at 513, 1025 and 2049 poloidal points.

**The VMEC twins.** `tests/reference/make_vmec_twin.sh` makes a twin in two steps. PB3D's `--export_HEL` writes a VMEC INDATA file from the HELENA mapping: the boundary as Fourier coefficients in the geometrical angle, pressure and iota as cubic splines against the normalised toroidal flux, and R_0 = B_0 = 1. VMEC++, the `vmecpp` package in a Python venv, then solves the fixed-boundary equilibrium. `tests/reference/vmec_twin.py` sets five multigrid stages up to 319 surfaces and a force tolerance of 1e−15.

| Twin | Surfaces | q on axis, at the edge | VMEC++ time | Date |
|---|---|---|---|---|
| circular ε = 0.33 (`wout_circ_eps033_q104.nc`, a5c8c41) | 319 | 1.04112, 1.03821 | 3 s | 30 October |
| circular, for TERPSICHORE | 159, 319, 639, 1279 | | | 1 November |
| circular, NTOR = 6 copy | | | | 30 October |
| cbm18a extended | 639 | 1.05212, 25.032 | 31 s | 31 October |
| cbm18a cut (`input.cbm18a_cut`) | 639 | q to five digits | | 1 November |

The cbm18a export needed one change. `--export_HEL` refuses a region of reversed shear, and HELENA's q dips by 5e−5 at its first three points, s < 0.004. A copy with those three values made monotone was exported. PB3D on this twin gives −4.3461e−2 fixed in the regression window against HELENA's −4.3573e−2, 0.26 % apart.

**Rippled twins (30 October).** The boundary of the circular twin is rippled as a(1 + δ cos 5ζ), with δ = 0, 0.01, 0.02, 0.04. VMEC++ runs with NTOR = 6 and writes N_fp = 1. The unrippled run at NTOR = 6 is the reference. VMEC's n = 0 content changes by 1.1e−5, 4.5e−5 and 3.0e−4 in R, as δ². Two derived copies exist on Snellius: `wout_circ_eps033_rip04rot.nc`, rotated by 3·2π/32 and written non-symmetric, and copies scaled in length and field (section 9.5).

**W7-X.** `wout_w7x.nc` is the STELLOPT standard-configuration file with 288 modes. Its two mode lists coincide: `mnmax_nyq = mnmax`. PB3D cannot run negative-flux VMEC files (section 11.7), so the runs use the copy with the field reversed, `wout_w7xflip.nc`. The decks of step 6 take ψ_N 0.85 to 1, n = 10 and 8 modes, m = −14 to −7, with 61 to 241 parallel points and 16 to 256 field lines. ι = 10/11 lies at s = 0.63, inside the plasma. VMEC's DMerc is negative only for s < 0.07.

**Other files.** `wout_cbm18a_small_640.nc` is a cbm18a VMEC file with q from 1.48 to 3.71. It is not a twin of the HELENA fixture: it has no extension. Its fixed-boundary run with the regression window does not converge, and SLEPc rejects its eigenvalue as complex (step 5).

### 3.2 External codes

**MISHKA-1.** Built on the Mac (27 October) from `Simulations/MISHKA/Sources`, the November 2017 source with MANZ = 21. That source reads the 2017 HELENA format, with `JS0, IAS` on the first line and `RAXIS, B0` after the metric. The 2015 source in `Documents/PHD/MISHKA_simulations/Sources` fails with "Bad integer for item 409".

- Build: `gfortran -std=legacy -w -O2 -fdefault-integer-8 -fdefault-real-8`, linked against an ILP64 OpenBLAS from `spack install openblas+ilp64 symbol_suffix=none`. With 4-byte integers it segfaults at once.
- Input and output: the namelist on stdin, the mapping in `fort.12`, the result in `fort.20`. The line "INSTABILITY FOUND" gives Re = γ².
- Deck: the archived `TEST_VAC/OLD/2.25/MISHKA/10/fort.10` with `NTOR` and `RWALL` changed; `IFAST = 1` shifts the m window with nq; `DSURF = 1` is a flat density.
- Validation: the archived PhD run, ε = 0.33, n = 10, `RWALL` = 3, reproduces 0.433847 to all printed digits with the archived `fort.12`, and gives 0.433764 with the regenerated one.

Pitfalls found:

| Pitfall | Date | What helped |
|---|---|---|
| `RWALL` = 1, a wall on the boundary, does not converge on the circular cases (99 iterations) | 27 October | PB3D's own fixed-boundary runs; MISHKA-1 used for free values |
| NaN from the first iteration on the extended cbm18a | 27 October | found on 31 October: `SQRT(2P/ρ)` with P a few 1e−20 below zero in the extension; clamped |
| its FFT needs 2(NCHI − 1) to be a power of two; the cbm18a mappings have NCHI = 401 | 31 October | a 513-point HELENA mapping, 33 MB, 801 radial points |
| 41 harmonics overflow the in-core workspace | 31 October | workspace enlarged |
| `-ffpe-trap=invalid` | 31 October | located the failing places |

MISHKA-1 is incompressible: with γ = 0 its results are identical to all digits (31 October). It is converged in `RWALL` from 1.5 on for m ≥ 11, in its radial grid NG from 201 to 401, and in its harmonics. On 1 November it gave 0.433764, 0.433765 and 0.433765 on the 257, 513 and 1025-point mappings of the circular case, and 0.433775 with 401 radial points. On the cut cbm18a the sources give fixed-boundary values with the wall on the boundary (1 November).

**HELENA `hel12a`.** The PB3D-patched HELENA at `~/Code/helena/hel12a` writes the extended mapping format, with `RAXIS, B0` and the full R and Z maps. It runs on the Mac as `hel12a < /dev/null`; the plot prompt at the end is harmless. Known pitfalls:

- the periodic branch, `IAS = 1`, needs NCHI = 2ⁿ (step 2);
- MISHKA-1 needs NCHI = 2ᵏ + 1;
- HELENA crashed writing its output at NRMAP = 1201 and 1601 (30 October);
- the 2049-point mapping of the circular case gives a `vac%style` 2 value that breaks the monotone sequence (1 November, section 7.3);
- the p′ and FF′ tables are on a uniform grid in √ψ_N (1 November).

**TERPSICHORE 1.2 (31 October).** TERPSICHORE is open source since 2025: EPFL GitLab, Apache 2.0. It was cloned on Snellius into `~/pb3d/ext/terpsichore`. The shipped Makefile builds it with gfortran 14.2 in 22 s without a change. Its QAS3 n = 1 benchmark, with 191 modes and 64 + 16 radial intervals, reproduces the reference eigenvalue −4.87061836e−5 to the last digit in 20 s and 1.2 GB.

TERPSICHORE reads neither netCDF nor namelists. Its equilibrium file `fort.18` comes from a VMEC 8 text wout through its own converter. Its control input is fixed-column cards. Two scripts bridge the gap: `scripts/terpsichore/wout2terps.py` writes `fort.18` from a netCDF wout, and `scripts/terpsichore/terps_data.py` writes the cards. The model has no compressibility term and no parallel inertia, like PB3D's ordering. Its vacuum is a pseudo-plasma inside a wall. Each run returns the one eigenvalue nearest the shift. Section 8.6 lists the traps.

**VMEC++.** VMEC++ runs in a Python venv on the Mac, from `~/Code/desc-vs-vmec`. The Mac has no VMEC binary. A PB3D run needs the VMEC file name to keep its `wout_*.nc` form (`Documentation/snellius-setup.md`).

**The cylinder formula.** `tests/reference/cyl_kink.py` derives with sympy the reduced-MHD energy of a straight cylinder with uniform current, zero pressure and uniform density, with the exterior Laplace solution written out. Without a wall:

ω²/ω_A² = −(2/q_a²)(m − nq_a)(nq_a − m + 1), unstable for m − 1 < nq_a < m.

With a conformal wall at radius b the factor λ = (1 + (a/b)^{2m})/(1 − (a/b)^{2m}) enters: ω²/ω_A² = −(m − nq_a)[(1 + λ)(nq_a − m) + 2]/q_a². The wall closes the window to m − 2/(1 + λ) < nq_a < m.

On the flat-q ε = 0.05 case without a wall the formula and MISHKA-1 agree to 1e−4 to 7e−4. With a wall, MISHKA-1 over the formula lies between 0.9992 and 1.0011 for every `RWALL` ≥ 1.05 and every n (27 October).

**Codes not built or not run.** COBRAVMEC (a STELLOPT clone exists at `~/Code/stellopt`), ELITE (not public), CAS3D and CASTOR3D (no public release found), GPEC's VACUUM and its Julia 3-D version, DESC's singular integrals, BIEST and bempp-cl. Section 12.2 lists what each would check.

### 3.3 Normalisations

| Code | Quantity printed | Relation to PB3D's ω²/ω_A² |
|---|---|---|
| PB3D | ω²/ω_A², MISHKA normalisation (`norm_style = 1`) | itself |
| MISHKA-1 | γ² = −ω²/ω_A² | sign change |
| TERPSICHORE, MODELK 2 after the parity fix | λ = WP/WK | λ = ω² (B_0/R_0)² |
| TERPSICHORE growth rate | GROWTH = √(\|λ\|/WALF) | (γ/ω_A)² = GROWTH² · WALF · R_0²/B_0² |
| cylinder formula | ω²/ω_A² with ω_A² = B_z²/(μ0 ρ R0²) | same |
| COBRA | τ_A = (μ0 ρ R0²/B0²)^½ | `norm_style = 2` states EV_COBRA = EV_PB3D · β/2 |

PB3D's ω_A uses B_0 and R_0 on the magnetic axis and a flat density. For VMEC files B_0 is the (m, n) = (0, 0) Fourier mode of |B| on the axis, the toroidal average, since 30 October (section 4.8). R_0 is the major radius of the axis. PB3D keeps the perpendicular kinetic energy only (`K_style = 1`), uses the incompressible ordering, has no wall and does not couple different n.

For twins written with B_0 = R_0 = 1, TERPSICHORE's λ equals ω². For the rippled twin, PB3D's B_0 is 1.003373 and R_0 is 1.000053, so the factor (B_0/R_0)² is 1.006652.

Before the parity fix of section 8.6, TERPSICHORE gave λ = +0.15099 fixed and −1.7248 with the wall at 2a on the circular twin. These are 1.93 and 3.93 times PB3D's values: no single factor relates them.

The published cbm18 curves use other choices of R and B. `Documentation/benchmarking_prep.md` gives the conversion (γ/ω_A)_published = (γ/ω_A)_PB3D × (R_pub/B_pub)/(R_mag/B_mag), with R_mag/B_mag = 1.3986 R_geo/B_geo for cbm18a. The choice of R and B alone spreads PB3D's converted n = 10 value from −8 % to +33 % against ELITE, so no absolute comparison has been made.

### 3.4 Sign conventions that mattered

| Convention | What holds now | Where it mattered |
|---|---|---|
| vacuum term in the energy | 2δW_vac = −2π X_e^H res X_e; the exterior Neumann-to-Dirichlet map of the BEM is negative definite | step 3: the vacuum destabilised |
| normal handed to the kernels | inward; `store_vac_HEL` stores −J∇ψ; `store_vac_VMEC` stored +J∇ψ until 30 October | step 6e: 3-D response 3 to 4 times too small |
| normal data of the potential | `vac%style` 1 takes them along J∇ψ; `vac%style` 2 takes −J∇ψ·∇φ, the opposite sign | July: `vac%style` 1 response flipped and divided by the α span |
| side of the boundary | `vac%style` 2 solves (H + 4πI)Φ = G R; `vac%style` 1 needs no shift, because its row-sum diagonal and the inward normal already give the exterior identity | July: interior problem solved |
| VMEC angles | boundary points follow VMEC's ζ_E; metric and mode phases use the flux frame, ζ_F = −ζ_E | checked on 30 October; not a fault |
| sign of iota | a copy of the twin with positive iota gives the same eigenvalues to 1e−9 | ruled out on 30 October |
| sign of the toroidal flux | PB3D's HELENA twins carry a negative flux; `wout2terps.py --flip` reverses it for TERPSICHORE | 31 October |
| mode numbers | the W7-X deck has m = −14 to −7; the report's claim that every VMEC case has nq < 0 is false for the cbm18a VMEC file | step 5 |
| MISHKA-1 toroidal number | `NTOR = −n` in its namelist | `tests/reference/mishka_circular.txt` |
| SLEPc target | PB3D asks for the eigenvalue nearest `EV_guess` = −0.3, not the most unstable one | step 6e: eigenvalue series hop between branches |

## 4. Defects found and fixed in the PB3D and POST code

Each subsection gives the symptom, the cause, the fix with routine and file, the effect on the numbers, and the test that now pins it.

### 4.1 The PETSc port (July)

- **Symptom.** PB3D 2.47 used the PETSc 3.9 API and did not build on a current machine.
- **Fix.** Commit 775219a (7 July) ported `SLEPC_ops` and `SLEPC_utilities` to PETSc 3.25.3 and SLEPc 3.25.1. Commit 1862df0 (8 July) added guards for PETSc 3.19, the version on CI's Ubuntu.
- **Effect.** PB3D builds and runs on macOS arm64 with gfortran 16, on Snellius with gfortran 14.2 and on Ubuntu 24.04 with gfortran 13 and PETSc 3.19. On one process the four regression decks agree across the three to 1e−12 to 6e−12 relative (30 October; table in section 9.3).
- **Test.** Every layer of the suite. CI has no PETSc 3.25 leg; `Documentation/quality-and-ci.md` proposes a weekly one.

### 4.2 The interior problem instead of the exterior one (July)

- **Symptom.** The axisymmetric vacuum returned the interior Neumann-to-Dirichlet map. In the cylinder limit that flips the response from +|m| to −|m|, so stabilising and destabilising swap (commit message of 429f003).
- **Cause.** `solve_Phi_BEM` in `vac_ops.f90` solved H Φ = G R. The assembled H of `vac%style` 2 satisfies the interior identity H φ = G ∂φ.
- **Fix.** 429f003 (8 July): `vac%style` 2 solves (H + 4πI) Φ = G R. The 4π term also makes the system nonsingular where the first and last points of a closed HELENA boundary coincide; there H alone has two identical rows.
- **Effect.** The sources give no eigenvalue before and after this change alone.
- **Test.** `vac_greens` checks both jump relations, their convergence, and the round trip of a known exterior harmonic, with STRUMPACK and with ScaLAPACK. `exact_exterior_2d` (31 October) checks the exterior solve against exact toroidal harmonics (section 7.2).

### 4.3 The 3-D kernel integrals of July

- **Symptom.** The field-line vacuum, `vac%style` 1, "was unusable" (commit message of d6d49e6, 9 July).
- **Causes and fixes**, in `calc_GH_int_1` (`vac_utilities.f90`), `calc_GH_1` and `calc_vac_res` (`vac_ops.f90`):

| Defect | Fix | Commit |
|---|---|---|
| singular interval integral had the same logarithm twice, diverged for coordinates at right angles and had the wrong sign | exact primitive of −1/distance over the metric half-cell (formula in section 6.1) | d6d49e6 |
| the regular end of the singular interval got the same analytic value | it keeps the regular kernel | d6d49e6 |
| far-field H used (−G)\*\*(−3), that is \|r\|³, instead of the dipole kernel n·r/\|r\|³: off by \|r\|⁶ | dipole kernel | 429f003 |
| interval weights and the −4π row sum used different units, which broke every jump relation at O(1) | each interval carries its surface weight Δp Δα/2 | d6d49e6 |
| the first point of each field line got zero weight, the last double | the bridge interval between two lines is excluded | d6d49e6 |
| `calc_vac_res` halved the weights of whole field lines through an integer division | weight halved by position along the line | d6d49e6 |
| `vac%style` 1 response on a different convention from `vac%style` 2 | sign flip and division by the α span, −res/((max_alpha − min_alpha)π) | d6d49e6 |
| `n_alpha` = 1 silently zeroed the singular integrals | rejected with an error | d6d49e6 |

- **Effect.** In July the `vac%style` 1 response of an axisymmetric torus matched `vac%style` 2 to 0.2 % at m = 1 and 9 % at m = 5, at 61 × 60. The October steps later found three more defects in the same path: the α spacing, the Richardson double count and the normal orientation. The 0.2 % at m = 1 turned out to be a cancellation (section 4.10).
- **Review.** Step 6a (27 October) read every hunk of 429f003 and d6d49e6 against the vacuum paper and its known errors. Verdict: keep all, revert nothing. The one fitted quantity, the factor −1/((max_alpha − min_alpha)π), is the α integral over one toroidal turn against the per-mode normalisation of `vac%style` 2, with the sign from the opposite conventions for the normal data.
- **Test.** `vac_3d`: `GH_int_1_singular` against brute-force quadrature, `greens_identity_3d`, `point_source_3d`, `response_1_vs_2`.

### 4.4 The vacuum sign and the factor 2π (step 3, 27 October)

- **Symptom.** A vacuum can only add energy, so ω²(no vacuum) ≤ ω²(free) ≤ ω²(fixed) at the same boundary. The code after step 2 broke this:

| cbm18a at ψ_N = 0.75, code after step 2 | ω² |
|---|---|
| no vacuum, diagnostic knob | −4.6041e−2 |
| free, `BC_style` 2 | −4.6377e−2 |
| fixed | −4.3573e−2 |

  The stored response had every eigenvalue negative: 20 of 20 on cbm18a, 10 of 10 on the circular case. Its Hermitian asymmetry was 2.7e−5.
- **Cause.** The vacuum term entered the matrix A with the wrong sign and without the factor 2π. A corresponds to 2δW with the toroidal integral included. In the large-aspect-ratio cylinder the vacuum energy of an edge displacement equals the line-bending energy of a plasma layer of thickness a/|m| at the edge. In A that energy is 4π²(nq − m)² |X_e|²/(μ0 R_0 |m|), while `res` up to that point is −2π(nq − m)²/(μ0 R_0 |m|). Hence 2δW_vac = −2π X_e^H res X_e. The derivation sits in the comment where the factor is applied.
- **Measurement before the fix.** A scratch build scaled the vacuum term on the ε = 0.05 fixture, free boundary, no wall. Values are ω²:

| Scale of the vacuum term | n = 10 | n = 15 | n = 20 |
|---|---|---|---|
| −2π | −0.43034 | −0.45161 | −0.33888 |
| −1 | −0.155 | −0.431 | −0.383 |
| −4π | −0.080 | −0.282 | −0.285 |
| 0, no vacuum | −0.545 | −0.615 | −0.392 |
| MISHKA-1, γ² | 0.42679 | 0.44844 | 0.33792 |

  Only −2π lands near MISHKA-1, at +0.8, +0.7 and +0.3 %. The scan is not monotone in the scale because SLEPc returns the eigenvalue nearest the target, not the most unstable one (section 6.2).
- **Fix.** Two lines at the end of `calc_vac_res` (`vac_ops.f90`): `vac%res = -2 pi vac%res`, then its Hermitian part. The Hermitian part matters because the discrete asymmetry of about 1e−5 makes SLEPc reject tightly converged eigenvalues as complex. The cylinder test's reference gained the same factor. POST's guard went too (section 4.16).
- **Effect.** cbm18a free at ψ_N = 0.75 became −4.5019e−2, between no vacuum, −4.6041e−2, and fixed, −4.3573e−2. Together with the shear fix this moved the free value from −4.3467e−2 at 17b712c to −4.5019e−2 (changelog). Near marginality on a circular fixture, n = 26 with nq_a = 26.98 and m = 27, PB3D gave −0.02969 against MISHKA-1's 0.03066 and the cylinder's 0.03118.
- **Test.** `check_vac_res.py` in `regression_cbm18a_free_bnd` requires the stored response to be Hermitian and positive semi-definite; it failed before the fix and passes after. `regression_circ_eps005_free_bnd` anchors the circular kink. `vac_response_cylinder` uses the reference times −2π. In the mutation check (30 October) the undone factor gave −0.1547 against the anchor −0.4303 and was caught; this needed a clean build, because an incremental scratch build once passed.

### 4.5 The HELENA shear (step 2, 27 October)

- **Symptom.** The shear S was zero at both midplanes. At θ = 0 the stored S was 1e−9 of max|S|, where the direct formula gives 3.4e−2. A boundary layer from both midplanes decayed by a factor 3.7 per interval. The cbm18a fixed-boundary eigenvalue was 6.5 % too stable and converged slowly in the parallel points:

| Parallel points, Simpson 3/8, before the fix | 97 | 100 | 193 | 385 |
|---|---|---|---|---|
| ω² | −4.0824e−2 | −4.0736e−2 | −4.2082e−2 | −4.2800e−2 |

- **Cause.** `calc_derived_q` in `eq_ops.f90`, HELENA branch, splines Θ^α in θ to get S through `calc_derived_S_from_deriv_HEL`. Θ^α is odd in θ, so a top-bottom symmetric equilibrium needs the odd condition `bcs(:,2) = [2,2]`, zero second derivative. Since June 2018 the call passed the even condition `bcs(:,1)`. The comments called `[2,2]` "zero first derivative".
- **Fix.** One argument changed and three comments corrected (28f8dab).
- **Effect.** The plan predicted −4.35682e−2 at 97 points with the trapezoidal rule; the run gave −4.3568189e−2. With the deck's Simpson 3/8 the value is converged to 1e−4 at 97 points:

| Parallel points, after the fix | 97 | 193 | 385 | 100 |
|---|---|---|---|---|
| ω² | −4.35691e−2 | −4.35699e−2 | −4.35714e−2 | −4.35729e−2 |

  The fixed anchor moved from −4.0736e−2 to −4.3573e−2. POST's energy balance closes to 0.078 % against the plan's 0.07 %.
- **Test.** `regression_cbm18a_shear_HEL` runs the fixed deck and `check_shear_HEL.py` compares the stored S with −(D_θh¹²/h²² − D_θh²² h¹²/h²²²)/J over the whole HELENA grid, tolerance 1e−4 on max|S|. Deviation 0.79 before the fix, 1.47e−5 after; the outboard-midplane value is 4.33e−2 from both routes. The mutation check caught the undone fix.
- **Side check.** The periodic HELENA branch, `IAS = 1`, gives −4.35487e−2 against −4.35878e−2 for the symmetric branch at equal poloidal spacing: 9e−4 apart. The shear is not the cause; both routes agree to 1e−5. This stays open (plan, section 9).

### 4.6 The outward normal in `store_vac_VMEC` (step 6e, 30 October)

- **Symptom.** The circular VMEC twin through the 3-D vacuum did not approach the HELENA value through `vac%style` 2, −0.443831, as the field lines grew (n = 10, 30 October):

| Field lines | 4 | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|---|
| ω² | +0.050011 | −0.356288 | −0.671667 | −0.194866 | −0.083333 | −0.049255 |

  The response matrix showed why. Its diagonal over that of the HELENA fixture, same case:

| Lines × parallel points | m = 6 | 8 | 10 | 11 | 12 | 14 | 15 |
|---|---|---|---|---|---|---|---|
| 8 × 100 | 0.95 | 4.37 | 1.76 | 1.30 | 1.14 | 1.35 | 1.73 |
| 16 × 100 | 0.24 | 0.27 | 0.37 | 0.46 | 0.60 | 1.12 | 1.58 |
| 32 × 100 | 0.249 | 0.229 | 0.215 | 0.209 | 0.206 | 0.205 | 0.208 |
| 32 × 200 | 0.267 | 0.252 | 0.243 | 0.240 | 0.239 | 0.242 | 0.247 |
| 64 × 100 | 0.293 | 0.280 | 0.271 | 0.267 | 0.264 | 0.259 | 0.257 |
| 64 × 200 | 0.289 | 0.275 | 0.263 | 0.257 | 0.251 | 0.241 | 0.236 |
| 128 × 100 | 0.306 | 0.296 | 0.288 | 0.284 | 0.281 | 0.275 | 0.273 |

  From 32 lines on the ratio converged at order about 1.8 to 0.29 ± 0.01. The 3-D response was Hermitian and positive definite, and three times too small.
- **How it was found.** The boundary points and normals of the two paths agreed to five digits. A copy of the twin with positive iota gave the same eigenvalues to 1e−9. The unit test `response_1_vs_2` run with the twin's geometry gave 0.865 to 0.769, so the method was sound and the real path was not. Two scratch switches settled it: mirroring the points into the flux frame changed nothing; flipping the normal gave 0.865 to 0.768, the unit test's numbers to three digits.
- **Cause.** `store_vac_VMEC` (`vac_ops.f90`) stored +J∇ψ, the outward normal. `calc_GH` takes the inward normal, which `store_vac_HEL` and the analytic torus of the test supply. With the outward normal the off-diagonal H has the wrong sign and the row-sum identity then sets a wrong diagonal.
- **Fix.** One minus sign, with the convention written at the place (e3b4847).
- **Effect.** Confirmed on Snellius to 11 digits the same evening:

| Field lines, 100 parallel points | ratio at m = 11, before | after | free ω², after | HELENA, `vac%style` 2 |
|---|---|---|---|---|
| 32 | 0.209 | 0.629 | −0.597 | −0.4438 |
| 64 | 0.267 | 0.799 | −0.527 | |
| 128 | 0.284 | 0.853 | −0.506 | |

  The W7-X smoke case barely moved, free +6.2297e−2 against +6.2295e−2, because its edge displacement is small. Sections 4.8, 4.11 and 6.3 brought the twin further.
- **Test.** `regression_circ_eps033_vmec_twin` runs the twin free through the 3-D vacuum. Its free anchor moved four times: −0.52743 within 25 % of HELENA (30 October), −0.51892 after the B_0 fix (30 October), −0.44910 within 5 % after the kernel changes (31 October), −0.43914 within 3 % after the edge weight (31 October). The mutation check caught a flipped normal.

### 4.7 The Nyquist mode matching in `read_VMEC` (step 6f, 30 October)

- **Symptom.** The unrippled twin rerun by VMEC++ with NTOR = 6 gave a fixed-boundary ω² of +0.108918 against +0.079716 for the NTOR = 0 file. The two files agree to 1e−13 in every n = 0 coefficient and have exactly zero n ≠ 0 coefficients. `alpha_style` was ruled out.
- **Cause.** VMEC stores the Jacobian, |B| and the covariant B components on its Nyquist mode set, which has its own (m, n) list. `read_VMEC` (`VMEC_ops.f90`) copied their first `mnmax` entries as if they followed the main list. That is right only when the lists coincide: for ntor = 0, and for files with `mnmax_nyq = mnmax` such as STELLOPT's W7-X file. The VMEC++ file with ntor = 6 had 143 of 150 coefficients mislabelled. The changelog puts the eigenvalue 37 % off.
- **Fix.** 0d49f19 matched every mode to its place in the Nyquist list. Step 9 (3b99109, the same night) replaced the matching: the full Nyquist set is kept with its own list `mn_V_nyq`.
- **Effect.** NTOR = 6: +7.9716389579421565e−2; NTOR = 0: +7.9716389579472524e−2; 6e−13 apart. W7-X unchanged to all digits.
- **Test.** The committed twin has ntor = 0, where the old code was right. The sources name no test that would fail if this fix were undone.

### 4.8 The B_0 normalisation (step 6f, 30 October)

- **Symptom.** The first pass of the rippled runs gave shifts linear in δ: fixed +2.7 % per 0.01, free −0.013, −0.024, −0.028 at 64 lines. That cannot be. δ → −δ is a toroidal rotation by π/5, so every α-averaged quantity is even in δ. A field-line scan was flat to all digits, because the `alpha` input does nothing with `alpha_style` 2. Five lines, commensurate with the ripple, gave a different number from six or more. The run logs showed the cause:

| δ | 0 | 0.01 | 0.02 | 0.04 |
|---|---|---|---|---|
| old B_0 [T] | 0.99190 | 0.97847 | 0.96573 | 0.94277 |
| new B_0 [T] | 1.000 | 1.0002 | 1.0008 | 1.0034 |

- **Cause.** `read_VMEC` summed all |B| Fourier modes at the axis: |B| at θ = ζ = 0. On a rippled axis that depends on the toroidal origin. The sum also picked up the m ≠ 0 modes extrapolated to the axis, so even the axisymmetric twin's B_0 was 0.8 % low.
- **Fix.** B_0 is the (0, 0) Nyquist mode of `bmnc`, extrapolated from the half mesh to the axis as 1.5 b₂ − 0.5 b₃ (0d49f19).
- **Effect.** The fixed-boundary twins, HELENA against VMEC, went from 1.6 % apart to 2.4e−6: +7.84304e−2 against +7.84306e−2. The 1.6 % had been blamed on the boundary fit all afternoon. The twin's free value at 64 lines became −0.518921. W7-X rescaled: B_0 from 2.506 to 2.429 T, fixed +6.6283e−2, free +6.6263e−2 at 16 lines and 61 points.
- **Test.** `regression_circ_eps033_vmec_twin` requires the HELENA and VMEC fixed values within 1e−4. On the old anchors it failed by 1.6 %, as it should; with the new ones it passes with deviations of 5e−6 and 3e−6. The scaling invariance of section 9.5 checks the normalisation as a whole.

### 4.9 `redistribute_var` and the multi-process shift (night of 30 October)

- **Symptom.** Two and four processes moved every eigenvalue towards stability, relative to one process:

| Deck | 2 processes | 4 processes |
|---|---|---|
| cbm18a fixed, 200 radial points | +2.4e−5 | +1.8e−5 |
| cbm18a fixed, 400 and 800 points | +2.6e−5, +2.6e−5 | |
| cbm18a free, `BC_style` 2 | +2.4e−5 | +1.8e−5 |
| cbm18a free, `BC_style` 4 | +2.4e−5 | +1.8e−5 |
| circular ε = 0.05 free | +1.6e−6 | +4.9e−6 |

  The shift did not depend on the radial resolution. It fell slowly with the parallel resolution: 1.05e−6, 0.86e−6 and 0.62e−6 absolute at 100, 202 and 400 points with two processes. Three processes gave only +6e−7, because their interfaces avoid the eigenfunction's peak at ψ_N = 0.58. `tol_SLEPC` from 1e−8 to 1e−12 moved the one-process value by 6e−9 and left the shifts unchanged.
- **How it was found.** The integrated perturbation matrices agreed for 1 to 4 processes to 1e−11. The assembled PETSc matrices, dumped with `PETSC_OPTIONS="-mat_view binary:…"`, differed only in the 12 to 16 rows after each process interface. The differences reached 6e−6 of the largest entry of A and 2e−5 of B.
- **Cause.** `redistribute_var` (`MPI_utilities.f90`) moves distributed variables between decompositions. It assumed that the source ranges do not overlap and filled each destination from the sources in rank order. Every process carries a ghost region of 4 × `norm_disc_prec_eq` = 12 equilibrium points. The overlap after an interface therefore came from the lower rank's ghost region, whose values near the end of its range are the less accurate ones.
- **Fix.** Each destination point is taken from its owner; the overlap of two neighbours is split in the middle, as `trim_grid` does (798fe50). A second, smaller part, 1e−6 relative, came from the cubic spline of `interp_V` over the local range without margin. `redistribute_output_grid` in `grid_ops.f90` now adds a 20-point margin.
- **Effect.**

| cbm18a fixed deck | 1 process | 2 | 3 | 4 |
|---|---|---|---|---|
| ω² | −4.3572898102e−2 | −4.3572897434e−2 | −4.3572897397e−2 | −4.3572898258e−2 |
| shift | | 1.5e−8 | 1.6e−8 | 3.6e−9 |

  The free deck shifts by 1.5e−8 and 1.3e−8 with 2 and 4 processes, the 400-point deck by 1.4e−8 with two, the circular free deck by 1.4e−12 and 1.8e−12. That is the level at which the eigenvalue itself reproduces.
- **Test.** The opt-in rows `regression_cbm18a_fixed_bnd_np2` and `_np4` require one process and `mpirun` to agree to 1e−6. They need a distributed LU in PETSc, which the Snellius build has and CI's does not. They pass on Snellius. The default rows use the 1e−3 anchor and never saw the shift.

### 4.10 The α spacing and the Richardson double count (steps 6b and 6c, 27 October)

**α spacing.** PB3D lays out the field lines with `calc_eqd_grid(..., excl_last=.true.)`: n_alpha lines with spacing (max_alpha − min_alpha)π/n_alpha. `calc_GH_1` and `calc_vac_res` used n_alpha − 1. The weights were too large by n_alpha/(n_alpha − 1), and line j's phase was off by n(j − 1)2π/(n_alpha(n_alpha − 1)). The test `test_vac_3d` built its own grid with both α end points, which hid the error.

The fix changed two lines (50c9511). The test grid now follows PB3D's layout. The old grid's duplicated line gave the system two identical rows; that is why STRUMPACK had fallen back to LU with a residual of 0.15 in the old test.

The new test `point_source_3d` puts a point source on the magnetic axis. Its exterior potential has exact normal data, and `solve_Phi_BEM` must return its boundary trace:

| Torus | 41 × 40 | 81 × 80 | ratio |
|---|---|---|---|
| circular | 4.15e−2 | 2.09e−2 | 0.50 |
| rippled, a(1 + 0.1 cos 5ζ) | 3.70e−2 | 1.86e−2 | 0.50 |

The errors halve with the resolution: first order. The normal orientation, the α weights and the exterior solve are right together on the analytic tori. The Green-identity residuals improved from 1.8e−2 to 7.5e−3 for the interior potential and from 6.2e−2 to 3.0e−2 for the exterior one, at 1640 points.

The `vac%style` 1 over `vac%style` 2 response ratios then converged towards 1 at first order:

| Grid | m = 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| 61 × 60 | 0.952 | 0.927 | 0.901 | 0.875 | 0.849 |
| 101 × 100 | 0.972 | 0.957 | 0.941 | 0.924 | 0.908 |

The old test's 1.002 at m = 1 was the duplicated line's extra weight cancelling against the near-singular solve.

**Richardson double count.** `store_vac_VMEC` copies the G and H of the previous Richardson level into the refined arrays, and `calc_GH_1` adds with `+=`. The copied entries were counted twice. `calc_GH` now zeroes G and H before either vacuum method (`vac%style` 1 or 2) adds to them (two lines). On W7-X with `max_it_rich = 2`, the level-2 response at 121 parallel points equals a fresh single-level run at 121 points to the last bit. The level-1 response equals the fresh 61-point run. The eigenvalues differ by 4.6e−6, inside `tol_SLEPC`.

**Aliasing warning (step 6d).** The sum over n_alpha lines aliases every toroidal harmonic that is a multiple of n_alpha into the average. After `read_VMEC`, `input_ops.f90` warns when LCM(n_alpha, N_fp) < n_alpha N_fp. On W7-X, n_alpha = 10 prints the warning and n_alpha = 16 does not. The condition is a shared factor, not "n_alpha divides N_fp" as the plan first said: 10 does not divide 5 but aliases the second field-period harmonic.

**Mutation check (30 October).** Undoing the spacing changes the surface weights by 1 to 2.5 % at 40 to 80 lines. The 3-D full-stack tests checked convergence orders, not values, and did not catch it. The same night `response_1_vs_2` gained its five recorded ratios at 81 × 80 with a 0.5 % tolerance, sized to catch the 1.25 % change. The sources do not report a rerun of this mutation, and the overview still lists the spacing mutation as not caught.

### 4.11 The free-edge trapezoidal weight (evening of 31 October)

- **Symptom.** The free-boundary decks converged only at first order in the radial points, while the fixed circular case converged at second order (night of 30 October):

| Radial points, before | cbm18a free, `BC_style` 2 | circular ε = 0.33 free | twin free, 64 lines |
|---|---|---|---|
| 200 | −4.50190e−2 | −0.443831 | −0.449100 |
| 400 | −4.51134e−2 | −0.438768 | −0.444073 |
| 800 | −4.51873e−2 | −0.436219 | −0.441544 |
| 1600 | −4.52253e−2 | −0.434941 | −0.440275 |
| 3200 | −4.52436e−2 | | |
| successive ratios | 1.3, 1.9, 2.1 | 1.99, 1.99 | 1.99, 1.99 |
| first-order limit | −4.5262e−2 | −0.43366 | −0.43901 |

- **Cause.** `fill_mat` (`SLEPC_ops.f90`) sums the energy over the radial points with the uniform weight h, also at the free end point, where the trapezoidal rule gives h/2. The discrete surface term then differs from the continuous natural boundary condition at O(h).
- **Fix.** A weight `w_q` = 1/2 at the last point for `BC_style` 2 and 4; three lines (4915d4b).
- **Effect.**

| Radial points, after | cbm18a free, `BC_style` 2 | circular ε = 0.33 free | twin free, 64 lines | circular ε = 0.05 free |
|---|---|---|---|---|
| 200 | −4.50394e−2 | −0.433799 | −0.439144 | −0.421677 |
| 400 | −4.51237e−2 | −0.433714 | −0.439058 | −0.421616 |
| 800 | −4.51925e−2 | −0.433683 | −0.439027 | −0.421585 |
| 1600 | −4.52279e−2 | −0.433670 | | −0.421573 |
| successive ratios | 1.2, 1.9 | 2.7, 2.4 | 2.8 | 2.0, 2.6 |

  The circular cases are second order and their 200-point values sit within 0.03 % of the limits. cbm18a stays first order with the same limit, because of its 22 rational surfaces. The 2.3 % gap to MISHKA-1 that the plan had taken for the toroidal correction of the high-n ordering was this quadrature (section 7.6). `BC_style` 4 on the circular case moved from −0.138 to −0.461. Fixed-boundary values and the vacuum identities are untouched: circular fixed +7.84304e−2, W7-X free at 16 × 61 +6.6271e−2.
- **Anchors re-recorded.** cbm18a free `BC_style` 2 −4.5039e−2, a move of 4.5e−4; `BC_style` 4 −4.5124e−2, a move of 0.9 %; ε = 0.05 free −0.42168; twin free −0.43914 against HELENA −0.43380.
- **Test.** `circ_eps005_free_bnd.sh` also runs its deck at 400 radial points and requires a move below 5e−4. The move is 1.5e−4 with the weight and was 1.2e−2 without it.

### 4.12 `BC_style` 4 and the default (step 4, 27 October)

- **Symptom.** At ψ_N = 0.75 `BC_style` 4 gave −4.5574e−2 for any vacuum scale, 4.8 % from `BC_style` 2.
- **Cause.** `set_BC_4` (`SLEPC_ops.f90`) inserted the edge (0,0) block with V1† + vac using `INSERT_VALUES`. The V2 loop then reached jd = 0 and inserted the same block again with `INSERT_VALUES`, which replaced it. The vacuum and V1† were lost.
- **Fix.** The whole (0,0) block is assembled in `loc_block` first, from V1† + vac and the jd = 0 part of V2, and inserted once. The V2 loop runs from −st_size to −1. About eight lines moved (83c02df).
- **Effect.** The overwrite is gone, but `BC_style` 4 is not equivalent to `BC_style` 2:

| Case | `BC_style` 4 | `BC_style` 2 |
|---|---|---|
| cbm18a, ψ_N 0.75, 200 points | −4.4712e−2 | −4.5019e−2 |
| absolute gap at 200, 400, 800 points | 3.1e−4, 3.3e−4, 3.4e−4 | |
| cbm18a after the edge weight (31 October) | −4.5124e−2 | −4.5039e−2 |
| circular ε = 0.05, n = 10, 15, 20 | −0.106, −0.292, −0.288 | −0.430, −0.452, −0.339 |
| circular ε = 0.05, n = 10, 400 points | −0.112 | |
| circular ε = 0.33 after the edge weight (31 October) | −0.461 | −0.434 |

  The gap does not shrink with the radial resolution. The gate "`BC_style` 2 = `BC_style` 4 to 1e−5" failed. A one-sided derivative in a replaced matrix row is not a second discretisation of the variational condition. The backward stencil gives the edge row of the variational scheme a factor step·ndc_e = 11/6 on the bulk terms relative to the surface term, which the explicit row ignores. Toon's 2017 note said "boundary style 4 does not work well yet" (`BC_style` 4).
- **Decision.** `BC_style` 4 stays in the code with the overwrite fixed, marked "not validated, use style 2" (that is, `BC_style` 2). It is anchored only as a change detector.
- **The default.** In July 41b5274 removed `BC_style` 3, a stub, and set the default to `BC_style = [1,4]`, because `[1,2]` pointed at another stub. 1d74b8f then implemented `BC_style` 2, but the default stayed at 4. Step 4 restored `BC_style = [1,2]`, Toon's choice from 2017 to July 2026.
- **Test.** `regression_cbm18a_free_bnd` anchors `BC_style` 4 at −4.5124e−2.

### 4.13 The `setup_modes` slots (step 5, 27 October)

- **Bug.** With `X_style` 2 the mode window follows nq. `nm_X` puts mode number m in slot modulo(m − m_min, M) + 1, which depends only on m. `ind_tot` computed a new mode's slot from the slot of the mode created last plus the sub-step. The two agree only while the window moves up. For nq < 0, reversed shear or jitter in a flat q, two active modes share a slot and the coupling coefficients mix without a message.
- **Fix.** Both branches use modulo(m_new − nm_X(1,1), n_mod_X) + 1, the formula of `nm_X`, in `setup_modes` (`X_ops.f90`, e4e3fdd). A consistency check after the loop requires every mode of `ind_tot` to find its own number in its slot of `nm_X` over its normal range; otherwise `CHCKERR`. It is always on and costs O(M n_r). The 2020 rewrite on `origin/issue_vac_not_used_VMEC` was not used: it is unfinished and wrong for jumps of 2 or more.
- **Test.** `tests/reference/setup_modes_sim.py`, registered as `unit_setup_modes_sim`. The old formula gives 8, 1, 8, 4 and 6 wrong entries on the five sequences that move down, reverse or jitter. The new one gives 0 on all eight, including jumps of 2 and 3.
- **Effect.** None measurable. Every equilibrium at hand moves the window up, so the check stays silent and every eigenvalue is unchanged. Negative n would move the window down on a HELENA fixture, but stops earlier (section 11.7).

### 4.14 PSPLINE (step 1, 27 October)

- **Symptom.** On macOS about every other launch died with SIGTRAP, a segfault or a failed `H5Fcreate`. This was the "MPI flake" of July.
- **Cause.** PSPLINE's `r8genxpkg.f` writes `xpkg(4,4)` for three unevenly spaced points, 8 bytes past the array. PB3D reaches it in `interp_V_spline` through `spline_real` for mode ranges with exactly three points. glibc hides the overrun in malloc padding; the macOS allocator traps at a later allocation. Ten singleton runs on Snellius with the unpatched library all succeeded.
- **Fix.** `Libraries/pspline-r8genxpkg-nx3.patch` holds the guard and its explanation. `Libraries/build_pspline.sh` applies it after the clone. CI keys its PSPLINE cache on both files. PB3D's own code is not touched.
- **Effect.** A patched build ran 23 of 23 launches clean in step 0. The step-1 gate ran 20 of 20 clean launches of each deck in Release and in Debug, with the fixed eigenvalue identical to the baseline and free values within 1e−11. Free-boundary eigenvalues move at 1e−12 between launches, so bit-identical gates apply to fixed-boundary runs only.
- **Test.** The regression runners launch PB3D once and count a crash as a failure; the six-attempt retry loop is gone. The upstream report text is in the step-1 comment on issue #4; posting it is Toon's call.
- **Found on the way.** The free-boundary deck had never run in a Debug build. `Libraries/dfft.f` declares its work arrays as `DIMENSION WSAVE(1)`, so `-fcheck=all` stopped the first FFT in `store_vac` with a false bounds error, 20 of 20 times. The `dfftpack` target now compiles with `-fcheck=no-bounds`.

### 4.15 The gfortran 14 bounds-check miscompile (30 October)

- **Symptom.** On Snellius the Debug suite passed 19 of 22 at first. The three regression decks segfaulted in gfortran 14.2's bounds-checking code.
- **Cause.** gfortran 14.2 and 14.3 miscompile `-fcheck=bounds` for the pattern `T_EF(:,:,:,c([1,1],.false.),...)`: a function with an array-constructor argument inside a subscript of an actual argument. gfortran 13.3 and 16.1 are fine.
- **Fix.** `CMakeLists.txt` replaces `-fcheck=all` by `-fcheck=all,no-bounds` for gfortran 14 and prints why. The reproducer is `tests/reference/gfortran14_bounds_repro.f90` (dfe439e).
- **Effect.** The Debug fixed deck then gives the Release value to all digits, −4.3572898102173241e−2. Snellius Debug ctest passed 24 of 24 on 31 October.
- **Open.** The Debug build on gfortran 14 has no bounds checking. The GCC bug report is Toon's to send.

### 4.16 The POST vacuum guards (step 3, 27 October)

- **Symptom.** POST never used the vacuum. In `sol_ops.f90` the vacuum energy term sat inside `if (1.eq.2 .and. ...)` since February 2018, with a print "NO VACUUM". In `driver_POST.f90` the vacuum read stopped with "Vacuum has not been implemented yet!".
- **Fix.** The guard and the print are deleted. `sol_ops.f90` now checks `allocated(vac%res)`, and `driver_POST.f90` reads the vacuum only when `BC_style(2)` ≠ 1, so fixed-boundary runs still work.
- **Effect.** On cbm18a free at ψ_N = 0.75 POST gives E_pot/E_kin = −4.49894e−2 against ω² = −4.50190e−2, 0.066 %. The vacuum term is 1.3 % of E_pot. On the circular case POST's vacuum term is 1.4045, equal to res₁₁,₁₁ |X_e|² = 1.4062 × 0.999, and to the cylinder formula within the 1.3 % error of the BEM. POST's total balance there is 2.1 % off, and 7.4 % off for the fixed circular run: that is POST's plot grid on an m = 11 edge mode, not the vacuum term.
- **Not done.** The overview says the energy term is back. The vacuum potential plot, `plot_vac_pot`, has not been checked on the 3-D vacuum, `vac%style` 1 (section 12.1).

### 4.17 Smaller repairs

| Repair | Where | Date |
|---|---|---|
| fixed-boundary runs skip the vacuum instead of erroring (HELENA) or computing it (VMEC) | drivers, `vac_ops` | July, ad2313b |
| free-boundary VMEC runs with `alpha_style` 1 stop at input with a clear message | `input_ops.f90` | July, cc4c96d |
| a 3-D vacuum with a single field line stops with a message | `calc_GH_1` | July, d6d49e6 |
| a 3-D vacuum with fewer than 4 points per line stops with a message, as the near-field stencil needs four | `calc_GH_1` | 31 October, 0bbda47 |
| CI's warning counter counts unique source locations; baseline 74 | `.github/workflows/tests.yml` | 27 October, step 1 |
| four `-Wcompare-reals` warnings in `fourier2real` removed | `VMEC_utilities.f90` | 30 October, 329b7a4 |
| a lowercase `\note` keeps the Doxygen count at its baseline | `vac_ops.f90` | 31 October, d8a1c6c |

PB3D and POST still exit silently with status 0 when a path exceeds 120 characters. This bit twice in steps 2 and 3. It is not fixed; the scripts use short paths and symbolic links.

### 4.18 The corner term of the axisymmetric vacuum (2 November)

- **Symptom.** `vac%style` 2 converged at first order in HELENA's poloidal points (section 7.3). On the circular ε = 0.33 kink at n = 10, free ω² was −0.43380, −0.43582 and −0.43695 at 257, 513 and 1025 points against TERPSICHORE's −0.43848, while the fixed-boundary value moved by 1e−5.
- **Cause.** `calc_GH_2` adds the corner term −2β to the diagonal of H, with β the boundary's angle on the vacuum side: cos(β − π) is the scalar product of the two unit chords at the point, and β > π on a convex boundary (eqs. 145–146 of the 2019 manuscript on the vacuum energy, T. Weyens, "Expressions for perturbed vacuum potential energy for 3D linear MHD stability", 21 March 2019). The manuscript derives its boundary equation, eq. (22), from Green's identity (21) over the vacuum volume, but the bracket of (22) is the negative of the one (21) gives with the vacuum's outward element −J∇ψ dα dθ. As printed, (22) holds for potentials harmonic inside the plasma. The H assembled from it is the interior operator, as July's jump-relation measurements found (section 4.2), and July's +4π in `solve_Phi_BEM` moved it to the vacuum side. The shift leaves the corner term as it was: the interior operator needs the plasma-side angle 2π − β, so every diagonal element was off by 4(β − π). On the HELENA mappings θ runs counterclockwise in (R, Z), the normal −J∇ψ points inward, and β − π, the angle through which the boundary turns at a point, is 0.008 to 0.019 at 257 points (mean 2π/512). The error is therefore about 8π/N for N boundary points around: first order.
- **Fix.** `- 2*(2*pi-beta(rd))` instead of `- 2*beta(rd)` in `calc_GH_2` (`vac_ops.f90`), one line; the comments of `solve_Phi_BEM` and `calc_GH_int_2` say which angle is used.
- **Effect on the circular ε = 0.33 kink,** n = 10, free boundary (`BC_style` 1 2), 200 radial points:

| HELENA mapping points | 257 | 513 | 1025 | 2049 |
|---|---|---|---|---|
| free ω², β (before) | −0.433799 | −0.435824 | −0.436949 | −0.436699 |
| free ω², 2π − β (after) | −0.436786 | −0.437313 | −0.437692 | −0.437070 |
| change | −2.99e−3 | −1.49e−3 | −7.4e−4 | −3.7e−4 |
| fixed ω², unchanged | +7.84304e−2 | | +7.84296e−2 | +7.84271e−2 |

  The change halves with each doubling of the points. The diagonal elements `vac%res`(m, m) of the response matrix, for the ten poloidal harmonics m = 6 to 15 at the edge, fall by the same fraction for every m: after over before is 0.993 at 257 points and 0.998 at 1025, against the estimate 4/N of 0.78 % and 0.20 % (N = 512 and 2048). Against TERPSICHORE's −0.43848, the 257-point value is now 0.39 % off (1.07 % before) and the 1025-point value 0.18 % (0.35 %).
- **Effect on the full-stack tests:**

| Test | Before | After |
|---|---|---|
| `response_1_vs_2`, `vac%style` 1 over `vac%style` 2, m = 1 to 5, 300 boundary points for `vac%style` 2 | 0.994, 0.992, 0.988, 0.985, 0.982 | 1.006, 1.005, 1.002, 0.999, 0.996 |
| `vac_response_cylinder`, R_0/a = 20, 200 points, `vac%res`(m, m) against the cylinder limit, m = 1 to 5 | −1.0, +1.3, +1.6, +1.8, +1.9 % | −3.0, −0.7, −0.4, −0.2, −0.1 % |
| `exact_exterior_2d`, 301 points: exterior solve error | 0.89 to 2.44 % | 0.19 to 1.07 % |
| same, residual of the P-family (exterior) identity | 4.5e−3 to 1.2e−2 | 1.0e−3 to 5.4e−3 |
| same, residual of the Q-family (interior) identity | 1.0e−3 to 5.3e−3 | 4.2e−3 to 1.2e−2 |
| `vac_greens`, 101 points, interior potentials 1 and 2; exterior potential 3 | 3.5e−2, 1.8e−2; 1.3e−2 | 4.6e−2, 2.8e−2; 7.3e−3 |

  In `response_1_vs_2` the `vac%style` 2 response falls by 1.2 to 1.4 % (4/N = 1.3 %), in `vac_response_cylinder` by 2.0 % for every m (4/N = 2.0 %). The cylinder limit leaves out corrections of order a/R_0 = 5 %, which the test allows (25 %), so its m = 1 value is not an accuracy measure.
- **What the identities show.** The exterior identity, the exterior solve and every quantity PB3D uses improve. The interior identity, which PB3D does not use, worsens by the same amount: at n = 10, k = 15 the interior residual goes from 5.3e−3 to 1.2e−2 and the exterior one from 1.2e−2 to 5.3e−3, a change of 4(β − π)/(4π) = 6.7e−3 in the normalisation by max |4πφ|, with β − π = 2π/300. The first-order error that remains has opposite signs in the two identities, so it is not in the diagonal. Its source is not identified; all three columns of `exact_exterior_2d` still halve from 301 to 601 points.
- **Effect on the MISHKA-1 comparisons,** 257-point fixtures, 200 radial points, no wall; gap = (|ω²| − γ²)/γ², negative where PB3D is less unstable:

| Case | n | MISHKA-1 γ² | PB3D before | PB3D after | Gap before | Gap after |
|---|---|---|---|---|---|---|
| ε = 0.05 | 10 | 0.42679 | −0.421677 | −0.424515 | −1.20 % | −0.53 % |
| | 15 | 0.44844 | −0.445993 | −0.447368 | −0.55 % | −0.24 % |
| | 20 | 0.33792 | −0.337036 | −0.337474 | −0.26 % | −0.13 % |
| ε = 0.33 | 10 | 0.43378 | −0.433799 | −0.436786 | +0.004 % | +0.69 % |
| | 15 | 0.44283 | −0.445624 | −0.447062 | +0.63 % | +0.96 % |
| | 20 | 0.32262 | −0.325481 | −0.325993 | +0.89 % | +1.05 % |

  On ε = 0.05 the gap now falls as 1/n² (factor 2.2 from n = 10 to 15, 1.8 from 15 to 20; (3/2)² = 2.25 and (4/3)² = 1.78). On ε = 0.33 PB3D moves away from MISHKA-1, towards the converged value it shares with TERPSICHORE, from which MISHKA-1 is 1.1 % apart (section 11.1).
- **Effect on cbm18a.** The free regression deck (401-point mapping, 800 boundary points) moves by 1.5e−4 with `BC_style` 2 and 1.3e−4 with `BC_style` 4: this pedestal mode depends little on the vacuum.
- **Anchors re-recorded.** ε = 0.05 free −4.2168e−1 → −4.2452e−1; cbm18a free `BC_style` 2 −4.5039e−2 → −4.5046e−2, `BC_style` 4 −4.5124e−2 → −4.5130e−2; the HELENA reference of the twin test −4.3380e−1 → −4.3679e−1, against which the twin's 3-D value −0.43914 is now 0.54 % off (1.2 % before); the `response_1_vs_2` ratios as in the table.
- **Test.** The ε = 0.05 anchor (moves 0.67 %) and the `response_1_vs_2` ratios (move 1.2 to 1.4 %, tolerance 0.5 %) catch the old corner term. The ε = 0.05 deck first passed ctest 6.7e−3 off its anchor: ctest decided on the regular expression "PASS", and `check_vac_res.py` prints "PASS" for the response matrix before the script's anchor verdict. The regression tests now also carry `FAIL_REGULAR_EXPRESSION "FAIL"`.
- **Still open.** The `vac%style` 2 sequence is still not second order (steps of 5.3e−4 and 3.8e−4 from 257 to 513 to 1025 points), and the 2049-point value still breaks it, in both builds (section 11.2).

## 5. New capabilities

### 5.1 `BC_style` 2 (July)

`BC_style` 2 imposes the free boundary variationally. `set_bulk_lims` already includes the edge row in the bulk assembly for the free-boundary values of `BC_style`, and `fill_mat` assembles that row like any bulk row, with one-sided stencils near the edge. The only boundary term of the energy is X_edge^H res X_edge. `set_BC_2` adds its Hessian, `vac%res`, to the edge diagonal block of A. The natural boundary condition (V_1^H X + V_2 X′ + res X) = 0 at the edge then follows from the stationarity of the discrete energy.

Commit 1d74b8f (10 July) implemented it by removing the guard of a 2.44-era stub. Unlike `BC_style` 4, `BC_style` 2 keeps A and B Hermitian. In July it gave a spurious imaginary part about 180 times smaller than `BC_style` 4, 1.1e−10 against 2.0e−8. `BC_style` 2 has been the default since step 4 and is the validated boundary condition (sections 4.12 and 7).

### 5.2 Near-field integration along the field lines and the G diagonal (31 October)

Two changes to the field-line boundary elements (0bbda47) removed the grid dependence of the W7-X response. Section 6.1 describes the method and section 6.2 the evidence.

- `calc_GH_int_1` (`vac_utilities.f90`) takes a stencil of four consecutive points of a line. For every interval within four interval lengths of the influence point it integrates both kernels over each half-interval with 8-point Gauss-Legendre on the cubic interpolant of position and normal. Farther intervals keep the endpoint values; the singular endpoint keeps its analytic half-cell. The bridge interval between two lines is skipped before any kernel is evaluated.
- `calc_GH_1` (`vac_ops.f90`) sets the diagonal of G from the identity for the linear interior harmonics, H x − G n_x = −4π x and the same for y and z, projected on the normal. This takes three distributed matrix-vector products with H and three with G, the machinery of the row sum. The two diagonal loops became one internal routine, `set_diag`, and `vac_ops` lost 24 lines net.

The new full-stack test `GH_int_1_near` puts an influence point a quarter interval off an arc with the radius of curvature of the W7-X bean tip and the spacing of the 121-point lines. The routine agrees with brute force over the exact arc to 6e−4 for G and 1.2e−3 for H, the error of the cubic interpolant. The cost of a run changed by less than 10 %: the near field is a few per cent of the kernel evaluations.

### 5.3 The `check_vac` input (31 October)

`check_vac = .true.` makes `check_vac_exact` (`vac_ops.f90`) test the 3-D boundary operators on the actual boundary of the run, after `calc_GH`. It works for VMEC equilibria. It is the catalogue's T5 and T6 on the real surface instead of the analytic torus.

| Test potential | Harmonic | Check |
|---|---|---|
| 1, x, y, z, x² − y² | inside | Green identity residual, with the free-term offset c = −4π |
| point source on the magnetic axis at ζ = 0 | outside, decaying | identity residual; exterior solve against the exact trace |
| ring of 256 sources along the axis | outside | same |
| the same ring with weights cos nζ, n the primary mode number | outside | same; the toroidal structure of the modes |

The results print as `CHECKVAC` lines (6894199, 612fc1f). Since 83724cc the routine also writes per-point files `check_vac_identity_<potential>.dat` and `check_vac_solve_<potential>.dat`, with point index, parallel index, field line, position and value. Since 13bdef4 `solve_Phi_BEM` factors with `pdgetrf` and solves with `pdgetrs`, and under `check_vac` estimates the reciprocal 1-norm condition number with `pdgecon` in between. `scripts/snellius/checkvac_table.py` turns the log lines into one table row per run. A run's eigenvalue is unchanged by the check.

### 5.4 Trigonometric factors in basis form and the full Nyquist set (step 9, 30 October)

Before step 9, `calc_trigon_factors` stored cos(mθ − nζ) and sin(mθ − nζ) for every mode of the main VMEC list and every grid point. Every Fourier reconstruction, `fourier2real`, summed against that table. The Jacobian, |B| and B_i could therefore only be rebuilt on the main modes, and their extra Nyquist content was dropped.

Step 9 (3b99109) stores the one-angle basis instead: cos mθ and sin mθ for m = 0 to `mpol_tf` − 1, and cos nN_fp ζ and sin nN_fp ζ for n = −`ntor_tf` to `ntor_tf`. `fourier2real` forms the per-mode products for whatever mode list it is handed. `read_VMEC` keeps J, |B| and B_i on the full Nyquist list `mn_V_nyq`, written to the output file with `mpol_tf` and `ntor_tf`. POST refuses older output files with a message.

| File | Modes before | Basis entries after |
|---|---|---|
| circular twin | 12 | 17 |
| NTOR = 6 twins | 150 | 33 |
| W7-X | 288 | 37 |

The plan expected the extra product per mode and point to cost nothing. The first build tripled the equilibrium phase of the NTOR = 6 twin, 40 → 137 s per equilibrium job, because the basis with the mode index first was read with a stride through memory. The table now has the grid points first, and `fourier2real` loops over the normal points outside and the modes inside. The equilibrium phase then took half its old time (section 10.4).

A consolidation followed (1a53f50, 31 October). The angle version of `fourier2real` builds the basis for its angles and calls the table version: one code path, 90 lines fewer, eigenvalues within 6e−13.

The gate on Snellius, before and after:

| Case | Before | After | Change |
|---|---|---|---|
| circular twin, fixed | +7.84305776064e−2 | +7.84305776064e−2 | 4.5e−13 |
| circular twin, free, 64 lines | −0.518921354424 | −0.518921354431 | 1.4e−11 |
| rippled δ = 0.01, fixed | +7.83564105571e−2 | +7.83564105572e−2 | 1.2e−12 |
| rippled δ = 0.02, fixed | +7.81344313123e−2 | +7.81344309607e−2 | 4.5e−9 |
| rippled δ = 0.04, fixed | +7.72831270067e−2 | +7.72832418107e−2 | 1.5e−6 |
| rippled δ = 0.04, free, 64 lines | −0.496300637 | −0.496301131 | 1.0e−6 |
| W7-X, fixed and free | | | identical to 17 digits |

The rippled changes grow as δ², as the Nyquist content beyond (12, 6) does. W7-X is identical because the STELLOPT file's two lists coincide. The truncation that step 9 removed therefore costs nothing on the three files in use; a VMEC++ file with strong shaping is where it would show.

### 5.5 Consolidations

| Consolidation | Where | Lines | Date |
|---|---|---|---|
| PSPLINE wrappers split from `num_utilities` into `spline_utilities`; 34 use sites | `Modules/spline_utilities.f90` | | July, 8ac37d7 |
| `var_1D_type` moved to the light `var_1D_vars` | `Modules/var_1D_vars.f90` | | July |
| distributed vector helpers `n_col1_loc`, `vec_glob2dis`, `dis_matvec` moved from the test helpers | `vac_utilities.f90` | 57 fewer | 31 October, fe305fd |
| one routine `set_diag` for both diagonals | `vac_ops.f90` | 24 fewer net | 31 October, 0bbda47 |
| the angle version of `fourier2real` calls the table version | `VMEC_utilities.f90` | 90 fewer | 31 October, 1a53f50 |

### 5.6 STRUMPACK optional, ScaLAPACK fallback (July)

STRUMPACK-Dense 1.1.1 has become hard to obtain. Since 429f003 PB3D builds without it: `solve_Phi_BEM` then uses ScaLAPACK LU on the same distributed system, and CMake sets `PB3D_WITH_STRUMPACK` when it finds the library. Since d6d49e6 a STRUMPACK solution is checked against the uncompressed operator and falls back to LU when the randomized compression fails. It fails on the 3-D kernels: residual 1.75e2 at n = 3660 in July, and residuals 1e3 to 1e13 on every W7-X solve on 27 October. Snellius has no STRUMPACK; every number of October comes from the ScaLAPACK path. Making STRUMPACK's compression work on 3-D kernels is out of scope (plan, section 8).

### 5.7 A dependency-light core build (July)

`-DPB3D_BUILD_EXECUTABLES=OFF` builds only the core library `pb3d_core` and its unit tests. The core holds `num_vars`, `str_utilities`, `messages`, `files_utilities`, `dtorh`, the light type modules `var_1D_vars`, `grid_vars` and `X_vars`, `PB3D_utilities`, and with LAPACK `num_utilities`, `X_utilities` and `num_ops`. It needs a Fortran compiler and CMake, LAPACK optional. CI's unit and sanitiser jobs use it.

## 6. The 3-D vacuum in depth

### 6.1 The method: boundary elements on field lines

**Points.** For a VMEC equilibrium the boundary is covered by `n_alpha` field lines ζ = α + qθ, equally spaced in α over (max_alpha − min_alpha)π with the last line left out. Each line carries `n_par_X` points equally spaced in the parallel angle, over one poloidal turn by default. `store_vac_VMEC` stores for each of the n_bnd points its position, the inward normal −J∇ψ and the metric of the field-line coordinates, g_αα, g_αθ and g_θθ.

**Kernels.** G is the single-layer kernel −1/|x_s − x|. H is the double-layer kernel n_s·(x_s − x)/|x_s − x|³, the derivative of G along the source normal. `calc_GH_1` loops over the intervals between consecutive points of each line. Each interval gives its two endpoints a contribution weighted with the half-cell element Δp Δα/2. Seen from an influence point x, three zones exist:

| Zone | Intervals | Integration | Since |
|---|---|---|---|
| singular | x is an endpoint | exact primitive of −1/distance over the metric half-cell | July |
| near | an endpoint within 4 interval lengths of x | 8-point Gauss-Legendre per half-interval on the cubic interpolant through four points of the line | 31 October |
| far | all others | kernel value at the endpoint: the trapezoidal rule along the line | July |

The interval that bridges two lines is skipped. Across the lines, each line enters with one sample per point: a one-point rule in α.

**The singular cell.** With U = √g_θθ Δp/2, V = √g_αα Δα/2, s = g_αθ/√(g_αα g_θθ) and D± = √(U² ± 2sUV + V²), the half-cell integral of 1/distance is

F = U ln((V + sU + D₊)/(−V + sU + D₋)) + V ln((U + sV + D₊)/(−U + sV + D₋)).

The singular endpoint then contributes −2F/(Δp Δα √(g_αα g_θθ)) to G. The test `GH_int_1_singular` checks this against brute-force quadrature.

**The diagonals.** The diagonal of H is set from the row sum H·1 = −4π·1, the identity of the interior harmonic 1. It holds the free term and the principal value together. Since 31 October the diagonal of G is set from the identities of the interior harmonics x, y and z, projected on the normal. The analytic half-cell integral is exact for the singular cell. But the one-point rule across neighbouring lines leaves an error that depends on the lattice of points; the identity puts that error into the diagonal, as the row sum does for H.

**The side of the boundary.** With the inward normal and the row-sum diagonal, potentials harmonic outside and decaying satisfy H φ = G ∂φ directly. Interior ones satisfy (H + 4πI) φ = G ∂φ. `vac%style` 1 therefore solves without a shift, the opposite of `vac%style` 2.

**The response.** `calc_vac_res` builds the normal data of each mode at every boundary point: (nq − m) times the cosine and sine of the mode phase. `solve_Phi_BEM` returns Φ for all 2 `n_mod_X` right-hand sides at once. The normal data are then weighted with the trapezoidal rule along each line, Δp Δα, halved at the first and last point of a line, and res = (I·EP)ᵀ Φ. For `vac%style` 1 the result is multiplied by −1/((max_alpha − min_alpha)π). The factor −2π and the Hermitian part follow for both values of `vac%style` (section 4.4).

**What the identities do not fix.** The off-diagonal far field stays the trapezoidal rule along open lines. The plan therefore called the method first order in the end corrections at the line ends (31 October). The projection of the modes onto the points stays the trapezoidal rule too. The 1 November study found that this projection, not the solve and not the line ends, sets the point count (section 6.3).

### 6.2 W7-X: the grid-dependence story

**27 October: the smoke test.** `wout_w7xflip.nc`, ψ_N 0.85 to 1, n = 10, 8 modes, 61 parallel points, 16 lines; `BC_style` 2 as the boundary condition, `vac%style` 1 as the vacuum. Fixed +6.2316e−2, free +6.2295e−2: free ≤ fixed held. The response was Hermitian positive definite, eigenvalues +0.066 to +27. At 121 points the response changed by 29 % and 4 of its 8 eigenvalues were negative, −3.1 to +32. STRUMPACK's compression failed on every W7-X solve. The verdict: sixteen lines do not resolve the W7-X boundary. These runs still had the outward normal and the old B_0.

**30 October, evening: the eigenvalue hides the response.** With the normal fixed, on Snellius:

| Parallel points | fixed | free, 16 lines | 32 lines | 64 lines |
|---|---|---|---|---|
| 61 | +6.231568e−2 | +6.229730e−2 | +6.229835e−2 | +6.230321e−2 |
| 121 | +6.231656e−2 | +6.229142e−2 | +6.229627e−2 | +6.230018e−2 |

Free lies below fixed everywhere. The vacuum lowers ω² by only 1.3e−5 to 2.5e−5, because this mode's edge displacement is small, and that difference moves by about 1e−5 across the grids. The response, not the eigenvalue, must be converged.

**Night of 30 October: no convergence seen.** Diagonal of the response at m = −7, smallest eigenvalue of its Hermitian part in brackets:

| Field lines | 61 parallel points | 121 parallel points |
|---|---|---|
| 16 | 13.3 (+0.16) | −12.4 (−14.2) |
| 32 | 47.0 (+0.13) | 19.5 (+0.11) |
| 64 | 92.4 (+0.98) | 46.3 (+0.49) |
| 128 | 44.4 (+0.46) | 72.9 (+0.88) |

The diagonal moved by factors of two between neighbouring grids, in either direction. At 16 × 121 the response was indefinite. The eigenvalue stayed between +6.6257e−2 and +6.6269e−2. The night's verdict: the first-order method cannot resolve the W7-X vacuum within a dense solve.

**31 October, weekend: a convergence that was not one.** `check_vac` put exact potentials on the W7-X surface (table in section 6.4). The operators were right everywhere; the identities held to 1e−2 and fell at first order. Two runs of about 31 000 points, 256 × 121 and 128 × 241, agreed to 2.6 % on the response. The whole diagonal moved together:

| m | −14 | −13 | −12 | −11 | −10 | −9 | −8 | −7 |
|---|---|---|---|---|---|---|---|---|
| 256 × 121 | 23.8 | 13.6 | 5.6 | 0.81 | 0.55 | 6.8 | 22.6 | 52.8 |
| 128 × 241 | 24.8 | 14.1 | 5.8 | 0.84 | 0.56 | 7.0 | 23.1 | 53.8 |
| 64 × 121 | 22.6 | 12.9 | 5.3 | 0.75 | 0.49 | 6.0 | 19.9 | 46.3 |

The weekend block read this as convergence at first order from 64 × 121 on, with 128 × 121 as an outlier 57 % off. A map around 128 × 121 found more bad grids: 112 × 121 and 128 × 141 were 160 % and 108 % off; 128 × 101 and 144 × 121 were within 20 % and 7 %. The condition number was 3.5 to 3.8 on every grid. The identity residual of x was 2 to 5 times larger on the bad grids.

**The wrong explanation.** The per-point files showed a spike at the inboard tip of the bean in the five symmetry planes, R ≈ 1.02 R_0, Z ≈ 0, where the curvature is highest and the axis passes 0.02 R_0 away. On 128 × 121 the point-source identity residual there was 0.137 against a median of 2e−4; on 144 × 121 it was 0.018. The plan blamed the curvature at the tips for a day.

**31 October, later: the 62 000-point run.** 256 × 241 lines and points, 61 696 boundary points, four processes with 200 GB, 72 minutes. Its diagonal at m = −7 was 67.5, against 52.8 and 53.8 for the two 31 000-point runs: 28 % off, 65 % at m = −14. Nothing was converged.

**The measured mechanism.** The per-point files of three grids showed that the surface is resolved. Along a line the points are 0.06 R_0 apart and the smallest radius of curvature is 0.27 R_0, a ratio of at most 0.22. What is not resolved is the kernel across the lines. At the inboard tips the lines run almost toroidally, so neighbouring lines pass 0.007 R_0 apart; that is the median over the surface at 128 lines, with 0.054 at most on the outboard side. From a point on one line, the kernel 1/r along the next line is a peak of width 0.007 sampled every 0.06.

The one-point rule sums those samples. The result depends on where the samples fall relative to the peak: on the fractional part of q Δθ/Δα, the offset along the line between points of neighbouring lines. At the edge of `wout_w7xflip.nc` q = 1.1313. That one number sorts every W7-X grid of the tables:

| Grid | Fractional part of q Δθ/Δα | Diagonal at m = −7 |
|---|---|---|
| 128 × 141 | 0.03 | 100.2 |
| 112 × 121 | 0.06 | 125.9 |
| 64 × 61 | 0.21 | 92.4 |
| 128 × 121 | 0.21 | 72.9 |
| 256 × 241 | 0.21 | 67.5 |
| 144 × 121 | 0.36 | 49.1 |
| 128 × 241 | 0.40 | 53.8 |
| 64 × 121 | 0.40 | 46.3 |
| 128 × 61 | 0.41 | 44.4 |
| 256 × 121 | 0.41 | 52.8 |
| 128 × 101 | 0.45 | 44.0 |

Points of neighbouring lines that line up put a sample on the peak and overestimate the near-line contributions. Staggered points underestimate them. On the circular twin the lines are 0.03 apart against 0.06 along them. There the Fourier component of the peak at the sampling wavenumber is 1 % of its integral, against 30 % on W7-X, which is why the twin converged and W7-X did not. The analytic self-cell is not involved.

**The two changes and the result.** The near-field integration alone brought six grids, 112 × 121, 128 × 101, 128 × 121, 128 × 61, 144 × 121 and 128 × 141, to 61.4 to 63.8, a 3 % spread against 44 to 126 before. It left 64 × 121 at 50.7. With both changes:

| Lines × points | n_bnd | diag m = −7 | eig min, max | identity x | identity z | ring | ring cos 10ζ |
|---|---|---|---|---|---|---|---|
| 64 × 61 | 3904 | 68.1 | 0.74, 68.1 | 1.8e−3 | 3.1e−3 | 0.0088 | 0.047 |
| 128 × 61 | 7808 | 64.1 | 0.73, 64.1 | 2.4e−3 | 3.0e−3 | 0.017 | 0.096 |
| 64 × 121 | 7744 | 66.1 | 0.75, 66.1 | 1.4e−3 | 2.0e−3 | 0.0066 | 0.017 |
| 112 × 121 | 13552 | 65.4 | 0.74, 65.4 | 6.5e−4 | 1.7e−3 | 0.0058 | 0.032 |
| 128 × 121 | 15488 | 65.9 | 0.74, 66.0 | 1.2e−3 | 1.3e−3 | 0.0059 | 0.016 |
| 144 × 121 | 17424 | 65.6 | 0.74, 65.7 | 1.5e−3 | 1.9e−3 | 0.0056 | 0.022 |
| 128 × 141 | 18048 | 65.3 | 0.74, 65.3 | 1.6e−3 | 1.9e−3 | 0.0031 | 0.020 |
| 256 × 121 | 30976 | 65.4 | 0.74, 65.4 | 7.8e−4 | 8.7e−4 | 0.0028 | 0.027 |

The six grids with 121 points or more agree to 1.2 % on the whole diagonal: 65.3 to 66.1 at m = −7, 35.0 to 35.9 at m = −14. The same grids differed by 160 % the night before. The first change alone converges to the same limit from below, 50.7, 63.2 and 64.5 at 64, 128 and 256 lines × 121 points, so the second change removes the first-order term without moving the limit.

The 61-point lines, 12 points per field period, stay apart by 3 to 6 %. The old 31 000-point value, 52.8, was 20 % low for the same reason as the twin's: the G diagonal. The W7-X eigenvalue of the step-6 mode is +6.62681e−2 on all five grids checked, identical to six digits.

**The lesson the plan drew.** The night's "no convergence" and the morning's "converges at 31 000 points" were both read off the response alone. The curvature explanation fitted the location of the error and survived a day. The spacing between lines against the spacing along them, measured on the grid itself, explained every grid. When a BEM result depends erratically on the grid, measure the lattice of points before the kernel (plan, section 7).

### 6.3 The convergence study on the circular twin

The circular twin has a known answer: its HELENA twin through `vac%style` 2, and since 31 October TERPSICHORE. Sections 4.6, 4.8 and 4.11 trace its free eigenvalue at 64 lines through the fixes, from −0.083333 with the outward normal to −0.43914 with all fixes.

The two quadrature changes on the twin, against HELENA through `vac%style` 2 at −0.44383 before the edge weight; the eigenvalue is in brackets:

| Field lines | before | along-line quadrature | both changes |
|---|---|---|---|
| 32 | (−0.5974) | 0.78 to 0.72 (−0.5555) | 1.00 to 1.20 (−0.4185) |
| 64 | about 0.85 (−0.5189) | 0.90 to 0.82 (−0.4971) | 0.981 to 0.993 (−0.4491) |
| 128 | (−0.5064) | 0.968 to 0.936 (−0.4607) | 0.980 to 0.962 (−0.4528) |

On the analytic torus at 81 × 80 the `vac%style` 1 over `vac%style` 2 ratios became 0.994 to 0.982, against 0.990 to 0.967 with the first change alone and 0.965 to 0.886 before.

**The study of 1 November.** PB3D ran the twin and the rippled twin free at n = 10 on 64 to 256 lines and 100 to 400 parallel points, with the response kept. The HELENA fixture was remapped at 513, 1025 and 2049 poloidal points. 157 lines make q n_alpha an integer, 1.03821 × 157 = 162.999, so the lines close on each other and the seam between one line's end and the next line's start disappears.

| 3-D path, lines × points | over `vac%style` 2 on 257 points, m = 6, 11, 15 | over `vac%style` 2 on 1025 points | twin ω² | rippled ω² |
|---|---|---|---|---|
| 64 × 100 | 0.981, 0.985, 0.993 | | −0.43914 | −0.418839 |
| 64 × 200 | 0.990, 0.990, 1.000 | | −0.43771 | −0.417492 |
| 128 × 100 | 0.980, 0.974, 0.962 | | −0.44293 | −0.422510 |
| 128 × 200 | 0.987, 0.980, 0.977 | 0.993, 0.988, 0.986 | −0.44133 | −0.420996 |
| 128 × 400 | | 0.997, 0.996, 0.994 | −0.43849 | −0.418224 |
| 256 × 100 | 0.978, 0.974, 0.965 | 0.984, 0.981, 0.974 | −0.44323 | −0.422795 |
| 157 × 100, closing lines | 0.979, 0.974, 0.963 | | −0.44306 | −0.422645 |

| Reference on the circular case | Value |
|---|---|
| TERPSICHORE, 159, 319, 639 surfaces | −0.438422, −0.438468, −0.438480 |
| TERPSICHORE on the rippled twin, 16 and 32 vacuum intervals | −0.419032, −0.419045 |
| PB3D `vac%style` 2 on 257, 513, 1025, 2049 poloidal points | −0.43380, −0.43582, −0.43695, −0.43670; with the plasma-side corner angle (2 November, section 4.18) −0.43679, −0.43731, −0.43769, −0.43707 |
| MISHKA-1 on 257, 513, 1025 points; 401 radial points | 0.433764, 0.433765, 0.433765; 0.433775 |

Readings of the plan:

1. The response converges in field lines by 128. 128, 157 and 256 lines agree to 3e−4, so the seam between line ends is not a source of error, and the 64-line values are under-resolved.
2. The parallel points set the accuracy. Against the 1025-point axisymmetric response the m = 11 ratio is 0.981, 0.988 and 0.996 at 9, 18 and 36 points per poloidal wavelength. The free eigenvalue at 128 lines goes −0.44293, −0.44133, −0.43849.
3. The boundary-element solve is not the slow part: `check_vac` at 128 lines is second order or better and accurate to 0.4 % at 200 points (section 6.4).
4. The axisymmetric reference itself was not converged (section 7.3).
5. Three methods meet: TERPSICHORE −0.43848, PB3D's 3-D path −0.43849 at 128 × 400, PB3D's `vac%style` 2 about −0.438 by extrapolation. This is the first time PB3D's 3-D vacuum reproduces another code's free-boundary eigenvalue to better than 0.1 %.

### 6.4 The exact-solution checks

**On the analytic torus** (`vac_3d`, after 31 October):

| Test | 41 × 40 | 81 × 80 |
|---|---|---|
| `greens_identity_3d`, interior residual | 1.0e−3 | 1.2e−4 |
| `greens_identity_3d`, exterior residual | 5.0e−3 | 1.8e−3 |
| `point_source_3d`, circular torus, 27 October | 4.15e−2 | 2.09e−2 |
| `point_source_3d`, rippled torus, 27 October | 3.70e−2 | 1.86e−2 |
| `response_1_vs_2`, ratios m = 1 to 5 | | 0.994 to 0.982; 1.006 to 0.996 since 2 November (section 4.18) |

**On the real surfaces, before the two changes (weekend of 31 October).** Error of the exterior solve, max |Φ − φ|/max |φ|; interior identity residuals for x and z in the last column:

| Surface, lines × points | n_bnd | point source | ring | ring cos 10ζ | identity x, z |
|---|---|---|---|---|---|
| circular twin 16 × 100 | 1600 | 0.109 | 0.031 | | 1.6e−2, 3.4e−2 |
| 32 × 100 | 3200 | 0.042 | 0.012 | 0.155 | 7.3e−3, 2.2e−2 |
| 64 × 100 | 6400 | 0.035 | 0.010 | 0.114 | 5.2e−3, 1.0e−2 |
| 128 × 100 | 12800 | 0.025 | | | 2.4e−3, 7.3e−3 |
| rippled δ = 0.04, 32 × 100 | 3200 | 0.040 | | | 7.2e−3, 2.2e−2 |
| 64 × 100 | 6400 | 0.034 | | | 5.2e−3, 1.0e−2 |
| W7-X 16 × 61 | 976 | 0.62 | 0.128 | | 1.6e−2, 9.2e−2 |
| 32 × 61 | 1952 | 0.56 | 0.110 | | 8.5e−3, 4.3e−2 |
| 64 × 61 | 3904 | 0.48 | | | 7.3e−3, 4.0e−2 |
| 128 × 61 | 7808 | 0.061 | 0.022 | 0.051 | 4.1e−3, 1.8e−2 |
| 16 × 121 | 1936 | 0.55 | 0.188 | | 1.4e−2, 7.2e−2 |
| 32 × 121 | 3872 | 0.30 | 0.058 | 0.274 | 7.0e−3, 3.8e−2 |
| 64 × 121 | 7744 | 0.068 | 0.020 | 0.046 | 3.7e−3, 1.7e−2 |
| 128 × 121 | 15488 | 0.30 | 0.070 | 0.221 | 3.5e−3, 2.0e−2 |
| 128 × 241 | 30848 | 0.086 | 0.0084 | 0.024 | 1.6e−3, 8.4e−3 |
| 256 × 121 | 30976 | 0.088 | 0.0086 | 0.027 | 1.6e−3, 9.1e−3 |

The source on the axis is 0.31 R_0 from the twin's surface but only 0.02 to 0.05 R_0 from W7-X's, so on W7-X the point source is a near-field test and the rings are the fair comparison. On the twin the cos 10ζ ring was only 11 % accurate at 64 lines and improved by 0.74 per doubling, the size and rate of the response deficit at the time. Its error was 2.5 times larger at the line ends than in the interior, 0.051 against 0.021.

**On W7-X after the two changes.** At 128 × 121 the exterior solves reach 0.6 % for the ring and 1.6 % for the cos 10ζ ring, against 7.0 % and 22 % before, and the z identity 1.3e−3 against 2.0e−2 (table in section 6.2).

**On the twin at 128 lines, 1 November.**

| Check | 100 points | 200 points |
|---|---|---|
| Green identities | 7.3e−5 to 2.0e−4 | 9.5e−6 to 3.5e−5 |
| exterior solve, point source | 1.4e−3 | 3.3e−4 |
| exterior solve, axis ring | 2.5e−4 | 3.6e−5 |
| exterior solve, ring cos 10ζ | 2.0e−2 | 3.8e−3 |

The overview summarises these as Green identities of 1e−5 to 1e−4 and exterior solves of 3e−4 to 4e−3 at 128 × 200, second order or better in the points.

**Invariances.** The rippled twin rotated by 3·2π/32 and written as a non-symmetric VMEC file gives a response identical to 1.5e−16 and ω² = −0.56285066726668 against −0.56285066726690 at 32 lines (31 October). The asymmetric path through the sine coefficients, the Nyquist arrays and the trigonometric basis is consistent with the symmetric one. Section 9.5 has the scalings.

### 6.5 The rippled twin and the sideband coupling

**The δ² law (step 6f, 30 October).** With B_0 corrected, n = 10 through the 3-D vacuum:

| δ | fixed | free, 32 lines | free, 64 lines |
|---|---|---|---|
| 0 | +0.0784306 | −0.587820 | −0.518921 |
| 0.01 | +0.0783564 | −0.586580 | −0.517794 |
| 0.02 | +0.0781344 | −0.582717 | −0.514300 |
| 0.04 | +0.0772831 | −0.562865 | −0.496301 |

The fixed and 64-line columns are direct runs; the 32-line column is the earlier runs rescaled by (B_0,old/B_0,new)². The shifts from δ = 0 are −7.4e−5, −3.0e−4 and −1.15e−3 fixed, ratios 4.0 and 3.9. At 64 lines they are +1.1e−3, +4.6e−3 and +2.3e−2, ratios 4.1 and 4.9; the last ratio carries the δ⁴ term and a less converged VMEC run at δ = 0.04. Gate 6f passed: smooth, no jump, δ². These free values predate the two quadrature changes and the edge weight.

**Against TERPSICHORE.** TERPSICHORE solves the family n = 0, 5, 10, 15, 20 of the five-fold ripple. Its external kink with the family is −0.419032 as λ, −0.41626 in PB3D's units; with the single row n = 10 it is −0.412459. The sidebands n = 10 ± 5 and ± 10 move its kink by 1.6 %. PB3D's field-line representation carries that coupling: on 31 October PB3D was 0.6 % more unstable than the family value at 64 lines and 1.5 % at 128, with the gap growing with the lines.

The 1 November study resolved the growth. The rippled twin follows the unrippled one point for point. At 128 × 200 it is 1.1 % from TERPSICHORE's family value, the same distance as the unrippled twin at that resolution. At 128 × 400 it is −0.418224 against −0.416263, 0.47 %, while the unrippled twin is at 2e−5. So about half a percent sits in the treatment of the ripple itself: PB3D's field-line average of the sideband coupling, or TERPSICHORE's truncation in modes and Boozer harmonics. The sources do not separate the two.

**Mercier (T10).** VMEC's DMerc says the δ = 0.04 twin is Mercier-unstable on every surface, DGeod −1.1e−3 against DShear +2.7e−5 at the edge. PB3D fixed at n = 20 and 40, 400 points, gives +2.056e−2 and +7.963e−2, within 0.2 % of the unrippled twin's +2.059e−2 and +7.975e−2. That is no contradiction at finite n, but T10 cannot serve as a check at the n PB3D runs (31 October).

### 6.6 The resulting user rule

| Rule | Evidence |
|---|---|
| at least 128 field lines | 128, 157 and 256 lines agree to 3e−4 (1 November) |
| about 36 parallel points per poloidal wavelength of the mode: 400 points over one turn at m = 11 | ratio 0.996 at 128 × 400; −0.43849 against TERPSICHORE's −0.43848 |
| 64 lines × 100 points gives about 1 % | overview |
| on W7-X at least 121 points per line; 61 points, 12 per field period, are too coarse | spread of 3 to 6 % at 61 points (31 October) |
| refine lines and points together and compare the two directions | 31 October |
| converge the response matrix, not the eigenvalue | the W7-X eigenvalue stayed between +6.6257e−2 and +6.6269e−2 while the diagonal moved from 44 to 126 |
| choose `n_alpha` without a common factor with N_fp | aliasing warning (section 4.10) |
| check `EV_guess` when the vacuum term changes a lot | SLEPc returns the eigenvalue nearest the target |
| `check_vac = .true.` once per new geometry | exact potentials on the real surface |

## 7. The axisymmetric vacuum

### 7.1 `vac%style` 2

HELENA equilibria run with `alpha_style` 1: one field line followed over the turns. `store_vac_HEL` stores the boundary points of the poloidal mapping, with the seam point duplicated, and the inward normal −J∇ψ. `calc_GH_2` builds G and H for one toroidal mode number from the toroidal harmonics P and Q of half-integer degree, evaluated by `dtorh`, with |n| for negative n. Near the singular point `calc_GH_int_2` uses the analytic asymptote Q_{n−1/2}(1 + x) = −½ ln(x/32) − b_n. The corner term of the diagonal of H is −2(2π − β), with β the boundary's angle on the vacuum side; it was −2β until 2 November (section 4.18). The exterior problem is (H + 4πI)Φ = G R, and the normal data are −J∇ψ·∇φ.

### 7.2 Exact and analytic checks

| Check | Result | Date |
|---|---|---|
| Green identities in `vac_greens`, 51 to 101 boundary points | residuals halve, e.g. 6.7e−2 → 3.5e−2: first order | 27 October |
| `vac_response_cylinder`, m = 1 to 5 | 1.0, 1.3, 1.6, 1.8, 1.9 % from the cylinder | 27 October |
| response diagonal of the ε = 0.05 fixture against 4π²(nq_a − m)²/(μ0 R_0 \|m\|) | 1.007 to 1.016 for m = 6 to 15 | 27 October |
| `exact_exterior_2d`: n = 2, k = 1 to 5; n = 10, k = 6 to 15 | identities 1e−3 to 1.2e−2; solve 0.9 to 2.4 % at 301 points, halving at 601 | 31 October |
| `unit_dtorh` at n = 20 to 80, z = 1.0001 to 100 | 1.4e−11 relative | 31 October |
| the same, with the plasma-side corner angle (section 4.18) | `vac_greens` at 101 points: interior 4.6e−2 and 2.8e−2, exterior 7.3e−3, still halving from 51; `vac_response_cylinder` −3.0, −0.7, −0.4, −0.2, −0.1 %; `exact_exterior_2d` at 301 points: solve 0.19 to 1.07 %, exterior identity 1.0e−3 to 5.4e−3, interior 4.2e−3 to 1.2e−2, all halving at 601 | 2 November |

`exact_exterior_2d` uses √(cosh η − cos θ) F^n_{k−1/2}(cosh η) cos kθ on the circular torus with R_0 = 3 and a = 1, with F = P outside and F = Q inside. It is the finite-aspect-ratio anchor of `vac%style` 2 that MISHKA-1 could not give. Fitting a constant to the identity, κ = −1.006 and −0.993, showed that `vac%style` 2 takes the normal data with the opposite sign of `vac%style` 1. Both are consistent, and each is pinned by its own exact solution. The precision loss of upward recursion in n that Chance, Turnbull and Snyder report does not affect `dtorh1`.

### 7.3 First order in HELENA's poloidal points (1 November)

| Poloidal points of the mapping | 257 | 513 | 1025 | 2049 |
|---|---|---|---|---|
| free ω², `vac%style` 2, circular ε = 0.33 | −0.43380 | −0.43582 | −0.43695 | −0.43670 |
| the same with the plasma-side corner angle (2 November) | −0.43679 | −0.43731 | −0.43769 | −0.43707 |
| fixed ω² | +7.84304e−2 | +7.84296e−2 | +7.84296e−2 | +7.84271e−2 (2 November) |

The response on the 1025-point mapping over the 257-point one is 0.994, 0.992 and 0.990 at m = 6, 11 and 15. The fixed value moves by 1e−5, so the first-order error sits in the vacuum. The 2049-point value breaks the monotone sequence; either HELENA's mapping at NPMAP = 2049 or the evaluation of the toroidal harmonics there is suspect. The limit lies near −0.438. Against TERPSICHORE's −0.43848 that is "about 0.1 to 0.3 % after extrapolation" (overview).

Most of this first-order error was the corner term of `calc_GH_2`, on the wrong side of the boundary (section 4.18). With the plasma-side angle the diagonal elements `vac%res`(m, m) on the 1025-point mapping over the 257-point one are 0.999, 0.997 and 0.996 at m = 6, 11 and 15, the steps from 257 to 513 to 1025 points are 5.3e−4 and 3.8e−4, and the 257-point value is 0.39 % from TERPSICHORE. The rest still converges slower than second order. The 2049-point value still breaks the sequence; the fixed value moves by 3e−5 relative there, and the response diagonal by at most 0.1 % in the direction that would lower ω², so the jump is neither in the plasma part nor in the response diagonal.

### 7.4 What the regression anchors therefore mean

The circular free anchors, −0.42168 for ε = 0.05 and −0.43380 for the `vac%style` 2 run of ε = 0.33, sit on 257-point mappings. The cbm18a free anchors sit on a 401-point mapping and carry the same first-order error of about 1 %. All free-boundary anchors pin reproducibility, not converged values. The cbm18a anchors are also 0.5 % from their radial limit (section 9.2).

Since 2 November the anchors are −0.42452 for ε = 0.05 and −0.43679 for the HELENA reference of ε = 0.33, 0.39 % from TERPSICHORE. The cbm18a free anchors moved by 1.5e−4 only, because this pedestal mode depends little on the vacuum (section 4.18).

### 7.5 The circular kink: four values

| Method | ω² or γ² | Converged in |
|---|---|---|
| MISHKA-1 | γ² = 0.43377 | mapping, radial grid, wall |
| TERPSICHORE | λ = −0.43848 | surfaces 159 to 639; vacuum intervals |
| PB3D, 3-D path, 128 × 400 | −0.43849 | lines; points to 0.4 % in the response |
| PB3D, `vac%style` 2 | −0.43695 at 1025 points; limit near −0.438; −0.43769 with the plasma-side corner angle (2 November) | slower than second order in poloidal points |

PB3D's two paths and TERPSICHORE agree at −0.4385. MISHKA-1 is 1.1 % less unstable. PB3D and TERPSICHORE with MODELK 2 weight the kinetic energy with the perpendicular displacement only. MISHKA-1 keeps the parallel displacement through incompressibility. The plan names that model difference as the likely cause; it is not verified (section 11.1).

### 7.6 How the "0.005 % against MISHKA-1" came about

| Date | PB3D `vac%style` 2, 257 points, 200 radial points | MISHKA-1 | Gap | Reading at the time |
|---|---|---|---|---|
| 2017, PhD | −0.298 | 0.434, `RWALL` 3 | 31 % | sign and factor of the vacuum wrong |
| 27 to 30 October | −0.443831 | 0.433779 | 2.3 % | the toroidal correction that the high-n ordering drops |
| 31 October, evening | −0.433799 | 0.433779 | 0.005 % | agreement |
| 1 November | −0.43695 at 1025 points, limit near −0.438 | 0.43377 | 1.1 % | coincidence |
| 2 November | −0.436786, with the plasma-side corner angle | 0.433779 | 0.69 % | most of the first-order error of `vac%style` 2 was the corner term (section 4.18) |

The 2.3 % was the free-edge weight (section 4.11), not physics. Once the weight was fixed, the first-order error of `vac%style` 2 at 257 points happened to land on MISHKA-1's number. MISHKA-1 itself sits 1.1 % from the two other codes. The changelog entry of 31 October still gives the 0.005 %, with a note that points to 1 November and calls it a coincidence.

## 8. Comparisons with external codes

### 8.1 Circular ε = 0.05 against MISHKA-1 and the cylinder formula

MISHKA-1, γ² against the wall radius in units of a:

| Case | n | `RWALL` 1.05 | 1.1 | 1.2 | 1.5 | 10 |
|---|---|---|---|---|---|---|
| flat q, not committed | 10 | 0.06443 | 0.33592 | 0.42327 | 0.43652 | 0.43663 |
| | 15 | 0.36337 | 0.43842 | 0.45472 | 0.45576 | 0.45577 |
| | 20 | 0.32596 | 0.34023 | 0.34223 | 0.34229 | 0.34229 |
| Γ = 1 − 0.02ψ, the fixture | 10 | 0.05473 | 0.32609 | 0.41342 | 0.42667 | 0.42679 |
| | 15 | 0.35606 | 0.43109 | 0.44739 | 0.44843 | 0.44844 |
| | 20 | 0.32160 | 0.33586 | 0.33786 | 0.33792 | 0.33792 |

The cylinder formula against MISHKA-1 on the flat case without a wall: 0.43660 against 0.43663, 0.45590 against 0.45577, 0.34254 against 0.34229. On the fixture MISHKA-1 over the formula is 0.977, 0.984, 0.987 at every wall: the current peaking alone. Gates R3 and R4 passed on 27 October.

| n | PB3D, 27 October | PB3D, after the edge weight | MISHKA-1 | Gap after | PB3D, plasma-side corner angle (2 November) | Gap |
|---|---|---|---|---|---|---|
| 10 | −0.43034 | −0.421677 | 0.42679 | −1.2 % | −0.424515 | −0.53 % |
| 15 | −0.45161 | −0.445993 | 0.44844 | −0.55 % | −0.447368 | −0.24 % |
| 20 | −0.33888 | −0.337036 | 0.33792 | −0.26 % | −0.337474 | −0.13 % |

The gap falls as 1/n². That is the finite-n correction that PB3D's high-n ordering drops, now visible under the discretisation error. With the corner term of section 4.18 the gaps are about half as large and fall as 1/n² more closely (factors 2.2 and 1.8 against 2.25 and 1.78). The fixed boundary is stable, +7.7650e−2 at n = 10. `regression_circ_eps005_free_bnd` anchors −0.42452 (−0.42168 until 2 November) and requires MISHKA-1's value within 2 %.

### 8.2 Circular ε = 0.33, fixed and free

**Fixed boundary, n = 10.**

| PB3D, HELENA | PB3D, VMEC twin | TERPSICHORE, corrected | Distance |
|---|---|---|---|
| +7.84304e−2 | +7.84306e−2 | +0.0785072 | 2.4e−6 between twins; 0.10 % to TERPSICHORE |

The fixed values at n = 16 and 20 are +7.5200e−2 and +1.2068e−2 on the HELENA path. Their VMEC counterparts were only computed with the old B_0.

**Free boundary.** MISHKA-1 prints γ², PB3D ω²:

| n | MISHKA-1 2017, `RWALL` 3 | MISHKA-1 now, no wall | PB3D 2017 | PB3D 30 October | PB3D after the edge weight | Gap after | PB3D, plasma-side corner angle (2 November) | Gap |
|---|---|---|---|---|---|---|---|---|
| 10 | 0.434 | 0.433779 | −0.298 | −0.443831 | −0.433799 | 0.005 % | −0.436786 | 0.69 % |
| 15 | | 0.442830 | | −0.452884 | −0.445624 | 0.63 % | −0.447062 | 0.96 % |
| 16 | 0.429 | 0.429043 | −0.430 | −0.438604 | | | | |
| 20 | 0.323 | 0.322624 | −0.330 | −0.329400 | −0.325481 | 0.89 % | −0.325993 | 1.05 % |

The 2017 PB3D values were wrong because of the sign and factor of section 4.4. Their near-agreement at n = 16 and 20 was luck; the n = 10 value was 31 % low. MISHKA-1 over the cylinder formula on this case is 0.990, 0.9755 and 0.965 at n = 10, 15 and 20, at every wall radius: the toroidal correction, not the vacuum. The gap to MISHKA-1 grows with n here, while it falls on ε = 0.05. Running both codes at 400 radial points at n = 15 and 20 would show whether that growth is resolution (`Documentation/benchmarking_prep.md`, item 3); it was not run. The converged values at n = 10 are in section 7.5.

### 8.3 The rippled twin, δ = 0.04, against TERPSICHORE

TERPSICHORE (31 October, late): `wout2terps.py --flip --nfp 5` keeps the 35 of 150 modes whose n is a multiple of 5 and declares five field periods. The Boozer table has m ≤ 12 and n/5 = −2 to 2. The stability table is either the row n = 10 with m = 1 to 20, or the family n = 0, 5, 10, 15, 20 with m = 1 to 25. MODELK 2 with the parity fix, 318 radial intervals, a conformal wall at 2a with 16 vacuum intervals. PB3D: n = 10, 20 modes, 200 radial points from ψ_N = 0.2.

| Run, PB3D units | PB3D ω² | TERPSICHORE, family | Distance |
|---|---|---|---|
| fixed | +7.72821e−2 | no discrete value | |
| free, 64 × 100 | −0.418839 | −0.41626 | 0.62 % |
| free, 128 × 100 | −0.422510 | −0.41626 | 1.5 % |
| free, 128 × 200 (1 November) | −0.420996 | −0.41626 | 1.1 % |
| free, 128 × 400 (1 November) | −0.418224 | −0.41626 | 0.47 % |

TERPSICHORE's single row n = 10 gives −0.412459 as λ free and +0.078853 fixed, 1.4 % above PB3D's fixed value as λ, +0.077796. The single row lacks the sideband coupling that PB3D carries, and its 1.4 % matches the size of that coupling in the free case. The fixed family spectrum is dense: the shifts 0.070, 0.078 and 0.086 return 0.07766, 0.07797 and 0.08601, with 634, 637 and 689 eigenvalues below them. Fifty-five eigenvalues lie between 0.070 and 0.086, the continua of the n = 5 and n = 20 rows. The fixed comparison is therefore ill-posed. Section 6.5 has the reading of the free values.

### 8.4 The cbm18a pedestal mode

**On the extended equilibrium (evening of 31 October).** MISHKA-1 ran on a 513-point mapping of the same input:

| n | MISHKA-1, γ² | PB3D, plasma to ψ_N = 1, ω² | Gap in ω² |
|---|---|---|---|
| 10 | 0.042118 at NG 201; 0.042134 at NG 401; 0.042323 with 41 harmonics | −4.514e−2, −4.473e−2, −4.482e−2, −4.488e−2 at 200, 400, 800, 1600 points | 6 % |
| 15 | 0.072840 | −8.353e−2 | 15 % |
| 20 | 0.094028 at NG 401; 0.094174 with 41 harmonics | −1.1125e−1 at 800 points | 18 % |

MISHKA-1 did not depend on the wall radius from 1.2 to 10 or on the density beyond ψ_N = 0.75. PB3D did not depend on its modes to 0.16 %. `Documentation/benchmarking_prep.md` gives the same comparison as 6.6, 15 and 18 %. PB3D's run with the pseudo-vacuum equals its true-vacuum run, −0.11125 against −0.11124 at n = 20. The suspected cause that evening was the rational surfaces of the extension, about 230 for n = 10, which MISHKA-1's grid resolves and PB3D does not see.

TERPSICHORE on the extended twin gave no usable spectrum. With the release code its converged eigenvalues sat at q = 1.09 to 1.11, an m = 11 resonance in the core that PB3D's window excludes. With the corrected kinetic energy and the shift at −0.045 it found no eigenvalue below −0.045, and landed on −0.01119 at q = 2.04 inside the extension, where ideal MHD allows no negative energy. On TERPSICHORE's toroidal-flux grid the extension takes three quarters of the surfaces.

**On the cut equilibrium (1 November).**

| n = 10, 15, 20 | PB3D, HELENA path | PB3D, VMEC twin | TERPSICHORE, twin | MISHKA-1, γ² |
|---|---|---|---|---|
| fixed | −0.04012, −0.08128, −0.11120 | −0.04020, not run, −0.10999 | −0.03991, −0.07951, −0.10701 | 0.03637, 0.07009, 0.09290, wall on the boundary |
| free | −0.05002, −0.08577, −0.11209 | | −0.04979, −0.08441, −0.10753, wall at 2a | 0.04492, 0.07407, 0.09453, wall at 3a |

The PB3D n = 10 values used 400 radial points on the 1025-point mapping; the other HELENA-path values used 200 points. MISHKA-1 with the wall at 10a matches 3a to 1e−5.

| Pair | n = 10 | n = 15 | n = 20 |
|---|---|---|---|
| PB3D against TERPSICHORE, fixed and free | 0.5 % | 2 %, 1.6 % | 4 % |
| PB3D, HELENA path against VMEC path | 0.2 % | | 1 % |
| MISHKA-1 against the two others | 9 to 16 % less unstable, fixed included | | |

PB3D's radial grid is first order on cbm18a and moves 1 % between 200 and 400 points at n = 20. TERPSICHORE's radial convergence at n = 20 is untested. The n = 20 distance is therefore a resolution figure. MISHKA-1 is converged in its radial grid, its harmonics, 21 against 41 at 0.08 %, the wall and the shift.

**Resolution.** The 6 to 18 % of 31 October was never the vacuum or the extension. It is a difference between MISHKA-1 and the two other codes in the plasma part, growing with n, with the same sign as the 1.1 % on the circular kink. Its cause is not identified: the inertia model or MISHKA-1's `MODE = 4` equations. ELITE would be the natural third opinion.

**The published cbm18 curves.** They have not been compared. Their ω_A is not pinned (section 3.3), their plasma edge is ψ_N = 0.701, and their compressible curves need Γ = 0 for PB3D.

### 8.5 W7-X

| Quantity | Value |
|---|---|
| PB3D fixed, edge window, 61 and 121 points | +6.6282557206921389e−2, +6.6283501023612579e−2 |
| PB3D free, 16 lines × 61 points, after step 9 | +6.6263019318233807e−2 |
| PB3D free after the two quadrature changes, five grids | +6.62681e−2 |
| TERPSICHORE, family n = 0 to 20, m ≤ 30, n = 10, fixed | 281 eigenvalues below 0.0124, 55 below 0.001, none below −1e−4; lowest +4.2e−6 |
| PB3D, full range ψ_N 0.2 to 1, 200 points, 20 modes, `input_w7x_fixed_full_n10` | no eigenvalue in 15 s; not pursued |

The lowest TERPSICHORE eigenvalue sits at the ι = 10/11 surface, s = 0.63, inside the plasma. W7-X at this pressure is stable at n = 10 with a dense spectrum in both codes, so there is no discrete mode to compare. TERPSICHORE needed its Boozer table at n/5 = −12 to 12 and m ≤ 14 to hold the 288 VMEC modes (section 8.6). A W7-X benchmark needs an equilibrium whose pressure is raised until a family goes unstable in both codes (section 11.4).

### 8.6 TERPSICHORE itself

**The benchmark.** The QAS3 n = 1 case reproduced to the last digit (section 3.2).

**Its vacuum.** The open release 1.2 treats the vacuum as a pressureless, currentless pseudo-plasma on extra surfaces out to a conformal wall; every free-boundary number in this report used a wall at 1.5a or 2a, which for these modes acts as no wall. The input `IGREEN` allocates a Green's-function vacuum matrix `wvac`, but nothing in the release fills it. The Green's-function TERPSICHORE results of Merkel, Nührenberg and Cooper, "Free-boundary ideal MHD modes in W7-AS", CRPP Lausanne report LRP 555/96, pp. 51–58 (1996) are therefore not reproducible with the open code. The agreement with PB3D is between two different vacuum methods: boundary elements without a wall, and a pseudo-plasma with a distant wall.

**The kinetic-energy parity error (31 October, night).** On the circular twin the release gave λ = +0.15099 fixed and −1.7248 with the wall at 2a; −1.7253 with the wall at 1.5a. A source reading found the cause in the MODELK 2 kinetic energy, `srcs/tprgl5_block_ap.f` lines 776 and 779. The ηη term, even in parity, is summed against the sine table, and the ξη term, odd, against the cosine table. For an up-down symmetric equilibrium both vanish, and the poloidal displacement carries no mass. A copy of the source with the two tables swapped, `srcs_fix`, rebuilt in 22 s:

| Circular twin, n = 10 | release | corrected | PB3D | Distance |
|---|---|---|---|---|
| fixed | +0.15099 | +0.0785072 | +0.078431 | 0.10 % |
| wall at 2a | −1.7248 | −0.438468 | −0.43914 at 64 lines, −0.43901 radially converged | 0.15 % |

The correction rests on the parity argument and on these two agreements, not on the authors. The report to EPFL is drafted and is Toon's to send (`external_asks.md`).

**Traps found.**

| Trap | Fix in the scripts | Date |
|---|---|---|
| the default parallel current from charge conservation, `LCURRF = 1`, loses a pedestal's edge current: −0.00502 instead of −0.03991 at n = 10 on the cut cbm18a | card default `LCURRF = 9`, current from the geometry; the circular kink does not care, 4e−5 | 1 November |
| `LPRESS = 2`, pressure gradient from the VMEC file | default; changed nothing on the cut case | 1 November |
| the toroidal convolution table needs `NMAXDF` ≥ 2n + 2; the default 38 fails at n = 20 | sized as 2 max\|N\| + 2 | 1 November |
| the VMEC input may hold at most 997 surfaces; a 1279-surface twin stops at once | none; limits radial convergence | 1 November |
| m up to 65 does not fit in 56 columns | optional `MPINIT` shifts the columns | 1 November |
| the Boozer table must hold at least as many modes as the VMEC input; `fvjac` is sized by one count and filled by the other, `tprgl0_block_ap.f:666` | W7-X needs n/5 = −12 to 12, m ≤ 14 | 31 October |
| a conformal wall on an axisymmetric equilibrium needs the Boozer rows n = −1 to 1 and a negative `RPLMIN` | in the card generator | 31 October |
| the kinetic energy is physical only for MODELK = 2 | MODELK 2 | 31 October |
| one eigenvalue per run, nearest the shift; a stable shift in a family run returns a continuum point | shift scans with the count of eigenvalues below each shift | 31 October |
| HELENA twins carry a negative toroidal flux | `wout2terps.py --flip` | 31 October |

With `LCURRF = 1` the n = 10 eigenvalue λ is one eighth of the right one, −0.00502 against −0.03991. λ scales as the growth rate squared, so the growth rate is about a third of the right one.

## 9. Tests and CI

### 9.1 The layers and their suites

`tests/CMakeLists.txt` registers 25 CTest entries, and 27 with `-DPB3D_TEST_MPI_REGRESSION=ON`. Labels select the layers: `unit`, `smoke`, `fullstack`, `regression`.

| Layer | CTest entry | What it checks |
|---|---|---|
| unit | `unit_str_utilities` | conversion, case and merge routines; the exact output formats |
| unit | `unit_files_utilities` | the file utilities (8a85ac2) |
| unit | `unit_dtorh` | P and Q of half-integer degree against 30-digit mpmath values, 54 points at 1e−11 and 32 more at n = 20 to 80 at 5e−11; the recurrence; rejection of z ≤ 1 |
| unit | `unit_num_utilities` | trapezoidal integration, finite-difference weights, the Björck-Pereyra solver, symmetric-storage indexing, sorting, LAPACK determinant and inverse, extrapolation, GCD and LCM |
| unit | `unit_num_ops` | Householder orders 1 to 3 and Zhang zero finders on known roots |
| unit | `unit_X_utilities` | mode-index logic, including the documented example of `get_sec_X_range` |
| unit | `unit_PB3D_utilities` | `setup_par_id`, `setup_rich_id`, `conv_1D2ND`; the partition property of the Richardson levels |
| unit | `unit_setup_modes_sim` | the slot simulation of section 4.13; needs python3 |
| smoke | `smoke_PB3D_usage`, `smoke_POST_usage` | usage messages |
| full-stack | `fullstack_vac_kernels` | interval kernels of `vac%style` 2 against direct evaluation, numerical derivatives and brute force; the Q asymptote |
| full-stack | `fullstack_vac_greens` | axisymmetric G and H: both jump relations, first-order convergence, the solve round trip, `vac_response_cylinder` |
| full-stack | `fullstack_vac_3d` | `GH_int_1_singular`, `GH_int_1_near`, `greens_identity_3d`, `point_source_3d`, `response_1_vs_2` with five recorded ratios at 0.5 %, `exact_exterior_2d` |
| full-stack | `fullstack_splines` | polynomial exactness, convergence rates, periodic ends, extrapolation |
| full-stack | `fullstack_calc_int_vol` | the analytic torus integral R_0π² + i·2π²/3 and its second-order convergence |
| full-stack | `fullstack_read_HEL` | golden-file test of the HELENA parser on the cbm18a fixture |
| full-stack | `fullstack_vac_greens_np2`, `_np4`, `fullstack_vac_3d_np2`, `_np4` | the distributed suites under `mpirun`, BLACS grids 1×2 and 2×2, block size 16 |
| regression | five entries, section 9.2 | eigenvalue anchors |
| regression, opt-in | `regression_cbm18a_fixed_bnd_np2`, `_np4` | one process against `mpirun`, 1e−6 |

The full-stack suites assert on globally gathered quantities, so the same assertions run on every process. The fixtures are committed with xz and decompressed at configure time, so the regression layer runs from a fresh clone.

### 9.2 The regression anchors

| Test | Deck | Anchor | Tolerance | Other assertions | Recorded |
|---|---|---|---|---|---|
| `cbm18a_fixed_bnd` | cbm18a, ψ_N 0.45 to 0.75, n = 10, 20 modes, 200 radial and 100 parallel points | −4.3573e−2 | 1e−3 | | 27 October |
| `cbm18a_free_bnd` | same deck, free | `BC_style` 2: −4.5046e−2; `BC_style` 4: −4.5130e−2 | 1e−3 each | free ≤ fixed; response Hermitian and positive semi-definite | 2 November |
| `cbm18a_shear_HEL` | fixed deck | S against the direct formula | 1e−4 of max\|S\| | outboard-midplane S non-zero | 27 October |
| `circ_eps005_free_bnd` | ε = 0.05, ψ_N 0.2 to 1.0, n = 10, 10 modes, 200 radial points | −4.2452e−1 | 1e−3 | fixed stable; MISHKA-1 within 2 %; 400-point value within 5e−4; response check | 2 November |
| `circ_eps033_vmec_twin` | ε = 0.33 and its twin, n = 10; twin free at 64 lines | VMEC fixed +7.8431e−2; twin free −4.3914e−1 | 1e−3 | twins within 1e−4; twin free within 3 % of HELENA's −4.3679e−1 (−4.3380e−1 until 2 November) and below fixed | 31 October; HELENA reference 2 November |
| `cbm18a_fixed_bnd_np2`, `_np4` | fixed deck | one-process value | 1e−6 | | 30 October |

The anchors reproduce to 1e−12 across three machines on one process. They are reproducibility checks. The cbm18a anchors are 0.5 % from their radial limit, because the cbm18a case converges only at first order in the radial grid. The radial limit of the fixed deck is about −4.3800e−2 by Richardson extrapolation from 200 to 3200 points:

| Radial points | 200 | 400 | 800 | 1600 | 3200 |
|---|---|---|---|---|---|
| cbm18a fixed ω² | −4.35729e−2 | −4.36594e−2 | −4.37296e−2 | −4.37652e−2 | −4.37827e−2 |
| circular ε = 0.33 fixed, HELENA | +7.84304e−2 | +7.83971e−2 | +7.83892e−2 | | |

The circular case converges at second order on the same HELENA grid, ratio 4.2. The cbm18a case has 22 rational surfaces in its window, each a point where the eigenfunction's derivative jumps. A benchmark at the 1 % level needs 800 points or a Richardson pair.

### 9.3 Reproducibility across platforms (30 October)

One process, the branch tip before the edge weight:

| Deck | macOS, gfortran 16, PETSc 3.25 | Snellius, gfortran 14.2, PETSc 3.25 with MUMPS | CI, Ubuntu 24.04, gfortran 13, PETSc 3.19 |
|---|---|---|---|
| cbm18a fixed | −4.3572898102123836e−2 | −4.3572898102173241e−2 | −4.3572898102188728e−2 |
| cbm18a free, `BC_style` 2 | −4.5018998767843765e−2 | −4.5018998767941465e−2 | −4.5018998768032115e−2 |
| cbm18a free, `BC_style` 4 | −4.4712139717450083e−2 | −4.4712139717485277e−2 | −4.4712139717494936e−2 |
| circular ε = 0.05 free | −4.3033715881358026e−1 | −4.3033715881628887e−1 | −4.3033715881637524e−1 |

The spread is 1e−12 to 6e−12 relative. The 17b712c baseline never ran on a green CI: step 1 repaired CI and step 2 changed the numbers in the same hour.

### 9.4 Multi-process rows

The full-stack suites `vac_greens` and `vac_3d` run at 2 and 4 processes in every CI leg, with results identical to one process. The regression decks run on one process in CI, because CI's PETSc has no distributed LU. The Mac cannot launch `mpirun`: Homebrew OpenMPI 5 crashes at launch for any program. The opt-in rows of section 9.2 therefore run on Snellius only. Sixteen MPI processes hang in the parallel HDF5 write of the flux equilibrium on the W7-X decks with 60 radial points; four processes work. This was noted and not pursued (31 October).

### 9.5 Invariance and monotonicity checks

| Check | Result | Date |
|---|---|---|
| 1, 2, 3 and 4 processes | 1.5e−8 | 30 October |
| `tol_SLEPC` from 1e−8 to 1e−12 | one-process value moves 6e−9 | 30 October |
| rigid rotation by 3·2π/32, non-symmetric file | response to 1.5e−16; ω² −0.56285066726668 against −0.56285066726690 | 31 October |
| rotation by 2π/N_fp | the VMEC file is literally unchanged; the pair is empty | 31 October |
| every length × 2 | fixed and free at 32 lines to all digits | 31 October |
| every field × 3, pressure × 9 | fixed 5e−14; free 5.5e−13 | 31 October |
| ripple δ → −δ | ω² even in δ, shifts as δ² | 30 October |
| `max_r_sol` 0.75, 0.85, 0.95, 1.0 | −4.3573e−2, −4.4753e−2, −4.4966e−2, −4.5144e−2: monotone | 31 October |
| prescribed windows [12,27], [10,29], [5,34], [0,39] | −2.5154e−2, −3.3723e−2, −4.1894e−2, −4.3478e−2: monotone | 31 October |
| fast window, 20 to 24 modes | rises by 6e−4 relative; section 11.5 | 31 October |
| HELENA symmetric against periodic branch | 9e−4 apart | 27 October |
| poloidal against toroidal flux label | 28 % apart; struck as ill-posed, section 11.7 | 31 October |
| n → −n | cannot run, section 11.7 | |
| Richardson restart against a fresh run; STRUMPACK against ScaLAPACK | not run | |

### 9.6 Mutation checks (30 October, Snellius)

Each fix was undone in a scratch copy of the source, rebuilt, and its test run.

| Fix undone | Test | Caught |
|---|---|---|
| HELENA shear, step 2 | `regression_cbm18a_shear_HEL` | yes |
| −2π factor, step 3 | `regression_circ_eps005_free_bnd`: −0.1547 against −0.4303 | yes, on a clean build only |
| normal orientation, 6e | `regression_circ_eps033_vmec_twin` | yes |
| α spacing, 6b | the 3-D full-stack tests | no; recorded ratios added afterwards, mutation not rerun |

The later fixes were not part of the mutation check. For some the test was shown to fail on the old code: the twin test on the old B_0 by 1.6 %, the circular 400-point check on the old edge weight by 1.2e−2 against 5e−4.

### 9.7 Status by date

| Date | Where | Result |
|---|---|---|
| 27 October | Mac, Release | 20 of 24; the four `mpirun` rows fail on the Mac |
| 28 October | Snellius, Release and Debug | 21 of 21 in 69.5 s; 18 of 18 without the regression layer |
| 30 October | Snellius, branch tip | Release 23 of 23 in 272 s; Debug 19 of 22, then the regression layer 5 of 5 with the gfortran 14 workaround |
| 30 October, evening | Snellius | Release 25 of 25 in 432 s |
| 31 October | Snellius | Release 27 of 27; Debug 24 of 24 |
| 1 November, overview | macOS and Snellius; CI | Release 27 of 27 on both; Debug suite on both; CI 8 of 8 |

The overview's "27 of 27 on macOS" is the later statement and stands; the plan's Mac numbers of 27 October predate it.

### 9.8 What fails or is not covered

- `BC_style` 4 is 7 % off `BC_style` 2 on the circular kink and first order in the radial grid; kept as a change detector.
- `vac%style` 2 is first order in HELENA's poloidal points; its sequence breaks at 2049 points. The free-boundary anchors pin reproducibility. Since 2 November most of the first-order error is gone (the corner term, section 4.18); the rest converges slower than second order, and the 2049-point break remains.
- Until 2 November a regression runner could fail its anchor and still pass ctest: ctest matched "PASS" anywhere in the output, and `check_vac_res.py` prints one before the anchor verdict. The ε = 0.05 deck passed 6.7e−3 off its anchor that way. The regression tests now also fail on "FAIL".
- cbm18a converges at first order in the radial grid; its anchors are 0.5 % from the limit, and the cut case moves 1 % between 200 and 400 points at n = 20.
- The overview lists the α-spacing mutation as inside the tolerances of the 3-D tests. The recorded ratios of 30 October were sized to catch it, but the mutation was not rerun.
- Sixteen processes hang on the W7-X decks.
- The Debug build on gfortran 14 has no bounds checking.
- The toroidal flux label is untested; negative n is rejected at input.
- POST's energy reconstruction is 7 % off on the circular edge mode, a plot-grid problem not pursued.
- No committed test uses a VMEC file whose two mode lists differ, so the Nyquist handling of sections 4.7 and 5.4 has no CTest guard. The sources name no CTest entry that runs `check_vac`.

### 9.9 CI jobs and ratchets

| Job | Purpose | Blocking |
|---|---|---|
| `unit`, Release | unit layer | yes |
| `unit`, Debug | unit layer with `-fcheck=all -ffpe-trap=invalid -finit-real=snan` | yes |
| `sanitize` | unit layer under `-fsanitize=address,undefined` | yes |
| `fullstack` | all layers, the regression decks and the 2- and 4-process rows; the warning ratchet; uploads the ctest log since 2d1363b | yes |
| `fullstack-debug` | Debug with `PB3D_ENABLE_DEBUG`; all layers but regression | yes |
| `coverage` | line coverage of `Modules/` with gcovr; uploads the ctest log | reporting |
| `lint` | fortitude; the OB category is a blocking gate | gate and report |
| `docs` | Doxygen build with a warning ratchet | ratchet |

| Ratchet | Baseline | Note |
|---|---|---|
| compiler warnings | 74 unique source locations in `Modules/`, `PB3D.f90`, `POST.f90`: 63 maybe-uninitialized, 11 unused dummy arguments | counted per location since step 1; the July count of 124 counted duplicates |
| Doxygen warnings | 374 | `.github/doxygen-warnings-baseline.txt` |
| fortitude | OB clean and gated; pinned to 0.9.2 | 0.10.0 adds rule OB091 and would fail two labelled loops in `dtorh.f90` |
| gcovr | pinned to 8.8 with `--gcov-ignore-parse-errors=suspicious_hits.warn` | step 1 |

Baselines only go down. Pushes made by the Claude GitHub App do not trigger workflows; user pushes, merges and manual dispatch do.

## 10. Infrastructure and cost

### 10.1 The Snellius stack

Set up on 27 and 28 October under `~/pb3d` on Snellius; `rm -rf ~/pb3d` removes all of it.

| Item | Version or path |
|---|---|
| modules | `foss/2025a`: GCC 14.2.0, OpenMPI 5.0.7, FlexiBLAS 3.4.5 over OpenBLAS 0.3.29, ScaLAPACK 2.2.2; parallel HDF5 1.14.6; netCDF-Fortran 4.6.2 |
| CMake | 4.4.4 from a venv; the site 3.31 loops for a quarter of an hour on PETSc 3.25's Fortran headers |
| PETSc | 3.25.3, complex scalars, with MUMPS 5.8.2, ParMETIS and METIS, built from source; the site module has real scalars |
| SLEPc | 3.25.1 |
| PSPLINE | patched with `pspline-r8genxpkg-nx3.patch` |
| LIBSTELL | minimal, `read_wout_mod`, from rsynced STELLOPT sources |
| STRUMPACK | not built |
| TERPSICHORE | `~/pb3d/ext/terpsichore`, release and `srcs_fix` |
| size | `~/pb3d/deps` 241 MB |

On the login node PETSc with MUMPS builds in 35 minutes and PB3D in 144 s for Release; `Documentation/snellius-setup.md` gives the other build times.

More than one MPI process needs MUMPS: the shift-and-invert solve needs a distributed LU, and without one `EPSSolve` aborts. The downloads of `build_deps.sh` are anonymous.

### 10.2 Scripts

Appendix A.4 lists the scripts. Two of them carry the job rules. `run_deck.sbatch <deck> [np-list] [bc-style] [repeat] [exe] [keep] [eq-file]` (`bc-style` is the edge value of `BC_style`) runs a deck once per process count with 16 tasks and 15 minutes, prints `EV np=<N>` lines, and deletes the 270 MB `PB3D_out.h5` unless the sixth argument is `keep`. `mem_run.sbatch` runs one executable on one deck under `/usr/bin/time -v` and prints wall time, peak memory and eigenvalue.

The job options are positional because Slurm on Snellius drops the submitting shell's environment. The first day's runs piled up 19 GB of output before `keep` existed. `sacct`'s MaxRSS is empty or wrong for jobs shorter than about a minute, so `mem_run.sbatch` measures memory itself. PB3D limits file paths to 120 characters, so `run_deck.sbatch` symlinks the equilibrium as `eq`.

### 10.3 Decks and fixtures

| Item | Where | Content |
|---|---|---|
| regression decks | `tests/regression/input_cbm18a_fixed_bnd`, `input_cbm18a_free_bnd`, `input_circ_eps005`, `input_circ_eps033` | the decks of section 9.2 |
| run decks | `~/pb3d/decks/` on Snellius, e.g. `input_twin_fixed`, `input_w7x_fixed`, `input_w7x_fixed_full_n10` | decks of the October runs; not in the repository |
| committed fixtures | `tests/fixtures/` | cbm18a extended, both circular cases, the VMEC twin; section 3.1 |
| Snellius fixtures | `~/pb3d/fixtures/` | `cbm18a_extended_vac.12`, `wout_circ_eps033_q104.nc`, `wout_w7x.nc`, `wout_w7xflip.nc`, the rippled and rotated twins, `cbm18a_cut.12`, `cbm18a_cut_nchi1025.12` |
| TERPSICHORE runs | `~/pb3d/ext/terpsichore/runs/` | one directory per case, e.g. `rip04_*` |
| run output and logs | `~/pb3d/runs/<job>-<jobid>/np<N>/`, `~/pb3d/logs/` | Slurm logs named after the jobs, e.g. `t3-rip04-*`, `t4-*` to `t9-*`, `terps-ns`, `terps-ivac` |

### 10.4 Run times and memory

| Run | Where | Wall time | Memory |
|---|---|---|---|
| cbm18a fixed regression deck | one Rome core | 4 to 5 s | |
| cbm18a free regression deck | one Rome core | 8 to 10 s | |
| Release ctest | Snellius | 70 s; 432 s for 25 entries on 30 October | |
| W7-X fixed, 121 parallel points, after step 9 | one core | 4.1 s | 821 MB |
| W7-X free, 64 lines × 61 points, after step 9 | one core | 19.2 s | 1574 MB |
| circular twin free, 32 lines, after step 9 | one core | 35.2 s | 9147 MB |
| NTOR = 6 twin, fixed, after step 9 | one core | 16 s | |
| NTOR = 6 twin, free, 64 lines, per equilibrium job | one core | 20 s | |
| W7-X 128 × 121, 15 488 points | one core | 13 min | |
| W7-X maps around 128 × 121 | four processes | | 100 GB |
| W7-X 256 × 121 and 128 × 241, about 31 000 points | four processes | 26 min each | 200 GB; 7.7 GB per dense matrix |
| W7-X 256 × 241, 61 696 points | four processes | 72 min | 200 GB |
| circular twin 128 × 400, about 51 000 points | eight cores | 39 min | 220 GB |
| TERPSICHORE build; QAS3 benchmark | Snellius | 22 s; 20 s | 1.2 GB |
| VMEC++, circular twin; cbm18a twin with 639 surfaces | Mac | 3 s; 31 s | |
| MISHKA-1 on the 513-point cbm18a mapping | Mac | 1 to 2 s | |
| HELENA cbm18a mapping with 513 poloidal and 801 radial points | Mac | | 33 MB file |

Step 9 saved 22 % of the W7-X run's memory, not the 3.6 GB the plan had estimated: the W7-X decks cover only ψ_N 0.85 to 1, so the old table was 230 to 460 MB.

### 10.5 The budget rules Toon set

| Rule | Detail |
|---|---|
| Snellius jobs | one node per job, partition `rome`, `--time` of 15 minutes or less, 30 for `mem_run.sbatch`; output under `~/pb3d/logs`; nothing left queued or running |
| parallel jobs | allowed (Toon, 30 October); up to six queued that night |
| where to run | PB3D on Snellius (Toon: "this should be on Snellius"); the Mac compiled and ran HELENA, MISHKA-1 and VMEC++ |
| Fable weekly use | under an 80 % ceiling: 38 % → 47 % (weekend), 56 %, 58 % → 63 %, 70 %, then 79 % on the night of 31 October, the rule's limit |
| from 1 November | 30 % of the all-models window; Fable may run out |
| outside requests | anonymous, with no name, address or token, and only to services that allow it; downloads of papers and data need Toon's approval |

## 11. Open questions and uncertainty

### 11.1 MISHKA-1's offset

- **Evidence.** MISHKA-1 is 1.1 % less unstable than PB3D and TERPSICHORE on the converged circular kink, and 9 to 16 % less unstable on the cut cbm18a, fixed boundary included. The gap grows with n and has the same sign on both cases. MISHKA-1 is converged in its mapping, radial grid, harmonics, wall and shift. Compressibility is ruled out: with γ = 0 its results do not change.
- **Candidate cause.** PB3D and TERPSICHORE with MODELK 2 weight the kinetic energy with the perpendicular displacement only. MISHKA-1 keeps the parallel displacement through incompressibility. Its `MODE = 4` equation set is the other suspect. Neither is checked.
- **Next experiment.** A MISHKA-1 run with the parallel inertia switched off, or the answer from the literature (plan, 1 November). ELITE on the cut cbm18a is the third opinion on the pedestal mode (section 12.2).

### 11.2 The `vac%style` 2 sequence beyond 1025 points

- **Evidence.** Free ω² of −0.43380, −0.43582, −0.43695 and −0.43670 at 257, 513, 1025 and 2049 poloidal points; the fixed value does not move (section 7.3). The 2049-point value breaks the monotone sequence. The plan names two suspects: HELENA's mapping accuracy at NPMAP = 2049, and the evaluation of the toroidal harmonics there. With the plasma-side corner angle (2 November, section 4.18) the sequence is −0.43679, −0.43731, −0.43769 and −0.43707: the corner term removed most of the first-order error, the steps 5.3e−4 and 3.8e−4 are still slower than second order, and the 2049-point value still rises, by 6.2e−4. The fixed value moves by 3e−5 relative from 1025 to 2049 points and the response diagonal by at most 0.1 %, in the direction that would lower ω²; the jump comes from neither.
- **Next experiment.** The sources name none. `exact_exterior_2d` runs `vac%style` 2 on an analytic torus without HELENA; running it at 1025 and 2049 boundary points would separate the two suspects [not in the sources]. Survey C points to product quadratures of the logarithmic singularity, as in Chance 1997 and Kress, as the higher-order route for the axisymmetric path.

### 11.3 The cost of the 3-D vacuum

- **Evidence.** 128 lines × 400 parallel points is 51 000 boundary points, 220 GB and 39 minutes on eight cores, for 0.1 % on the circular kink. At 64 × 100 the error is about 1 %. Each dense n_bnd × n_bnd matrix is 7.7 GB at 31 000 points. STRUMPACK's compression fails on the 3-D kernels, so the solve is a dense LU.
- **Next experiment.** None decided. Section 11.9 discusses a higher-order projection of the modes; section 12.5 lists other options to explore: a doubly periodic vacuum grid with a locally corrected quadrature, a Fourier basis as in CAS3D, and a solver that scales like the number of points. `Documentation/benchmarking_prep.md`, item 7, lists end corrections or a closed field-line lattice, with the Wu-Martinsson diagonal as a check. Toon decides whether any is tried.

### 11.4 W7-X needs an unstable case

- **Evidence.** At n = 10 the standard case is stable in PB3D and TERPSICHORE, with a dense spectrum and ι = 10/11 inside the plasma (section 8.5). The 3-D vacuum response on W7-X converges to 1.2 % since 31 October; the eigenvalue of the edge mode does not test it.
- **Next experiment.** Raise the pressure until a mode family goes unstable in both codes, then compare free boundaries. Candidates from the survey: the SIMSOPT test file `input.W7-X_without_coil_ripple_beta0p05_d23p4_tm` at ⟨β⟩ = 5 %, a β scan with VMEC++, and COBRA's marginal β.

### 11.5 The fast mode window

- **Evidence.** With `X_style` 2 the fixed cbm18a eigenvalue rises from −4.35729e−2 at 20 modes to −4.35450e−2 at 24, 6e−4 relative, and stays at 30 and 40 modes. The rise is the same at 100 to 400 parallel and 800 radial points. Prescribed windows behave monotonically. At n = 20, free, 800 points, 20, 30 and 40 modes give −0.111241, −0.111426 and −0.111427.
- **Mechanism.** `insert_block_mat` calls `setup_local_block` in `SLEPC_utilities`, which shifts each block by the change of mode numbers between the two radial points it couples and drops the entries outside the window. Wherever the window moves by one, the coupling between the mode that leaves and the one that enters is zero. The routine's comment says so. A trial space with dropped couplings is not nested, so Rayleigh-Ritz monotonicity does not apply.
- **Rule and fix.** Double `n_mod_X` once; if the eigenvalue moves by less than the accuracy wanted, the window is wide enough. The fix, a perturbation window widened by the stencil's mode change, is a design change and was not made.

### 11.6 `BC_style` 4

- **Evidence.** 7 % from `BC_style` 2 on the circular kink and first order in the radial grid. On cbm18a 0.7 % off before the edge weight, with a gap that did not shrink with the radial points, and 0.2 % after it (section 4.12).
- **What would settle it.** Its own derivation against the discrete energy, including the factor 11/6 that the backward stencil gives the bulk terms of the edge row. Out of scope; the plan forbids deleting `BC_style` 4 or adding a `BC_style` 5.

### 11.7 The toroidal flux label and negative n

- **Evidence.** With the toroidal flux label the parallel angle is ζ and `prim_X` is the poloidal mode number, so `prim_X = 10` means m = 10 with ten n's. On the twin, fixed, the two labels give +7.8431e−2 and +1.0046e−1, 28 % apart: different problems. The 30 % of the W7-X runs of 27 October has the same explanation. PB3D refuses the toroidal label for free-boundary runs, "TOR. FLUX HAS NOT BEEN TESTED" in `store_vac`, and for decreasing flux in `calc_eq_1`. Negative n stops in `calc_norm_range_PB3D_X` with "angular array has to have a length of at least 1".
- **Next step.** The input path for negative n needs a repair before the n → −n invariance can run. Negative-flux VMEC files are item 1 of the out-of-scope list: `origin/flip_angles` from 2018, or the field-reversal workaround made official. Both are needed before real stellarator work.

### 11.8 TERPSICHORE's radial convergence at n = 20

- **Evidence.** On the cut cbm18a PB3D and TERPSICHORE differ by 4 % at n = 20. PB3D moves 1 % between 200 and 400 radial points there. TERPSICHORE's VMEC input stops at 997 surfaces; the twin used 639. On the circular kink TERPSICHORE moved from −0.438422 to −0.438480 between 159 and 639 surfaces.
- **Next experiment.** TERPSICHORE on cut twins with increasing surface counts up to 997, and PB3D at 800 radial points or a Richardson pair at n = 20 [the surface counts are not in the sources].

### 11.9 Would a higher-order projection along the lines cut the point count?

Toon asked for this explanation on 1 November. The plan's claim: a higher-order projection of the modes along the field lines would cut the point count of the 3-D vacuum by about four for the same accuracy. The evidence behind it, all from the 1 November study:

1. The boundary-element solve is accurate. At 128 lines `check_vac` gives identity residuals of 9.5e−6 to 3.5e−5 and solve errors of 3.3e−4, 3.6e−5 and 3.8e−3 at 200 points, second order or better.
2. The field lines are converged. 128, 157 and 256 lines agree to 3e−4.
3. The response still moves with the parallel points:

| Lines × points | Points per poloidal wavelength at m = 11 | Response over `vac%style` 2 on 1025 points, m = 11 | Twin ω² |
|---|---|---|---|
| 256 × 100 | 9 | 0.981 | −0.44323 |
| 128 × 200 | 18 | 0.988 | −0.44133 |
| 128 × 400 | 36 | 0.996 | −0.43849 |

4. Between the solve and the response sits the projection, res = (I·EP)ᵀΦ: the trapezoidal rule along each line over the mode's normal data times Φ (section 6.1).

From these the plan concludes that the projection of an m = 11 mode on the field-line grid is the slow part, and that a rule of higher order would need far fewer than 36 points per wavelength. A dense matrix holds n_bnd² entries, 7.7 GB at 31 000 points, so a fourfold cut in the points cuts each matrix sixteenfold [computed from the 7.7 GB figure].

What complicates the claim:

- The 157 closing lines remove the line ends and the seam, and at 100 points they give the same response as 128 and 256 open lines. The error that the points control is therefore not at the line ends (plan, 1 November). End corrections of the Gregory or Alpert kind, which survey C lists, would then not reduce it [not in the sources].
- The sources do not state the order of the remaining error. The deficits 1 − ratio of the table above are 0.019, 0.012 and 0.004, which fall by 1.6 and then 3 per doubling [computed from the table].
- No run with a higher-order projection exists. The factor of four is an estimate, not a measurement.
- (2 November) Point 1 overstates the solve for oscillating data. The cos 10ζ ring, which has the modes' toroidal structure, is solved to 2.0e−2 at 100 points and 3.8e−3 at 200, at 128 lines (section 6.4). At 100 points that is the size of the response deficit, 1.9 % at 256 × 100; at 200 points it is a third of it, 0.38 % against 1.2 %. The solve and the projection therefore both contribute, and a higher-order projection alone may not give the factor of four. The quadrature options of section 12.5 act on the solve.

**Next experiment.** The plan names the exact potentials of `check_vac` on the real surface as the way to measure the order (31 October). Projecting the exact boundary trace of the cos 10ζ ring onto the modes with the present rule and with a higher-order one, at 100, 200 and 400 points, would measure both orders without a stability run [not in the sources].

## 12. Work bucket and roadmap

### 12.1 Toon's items (31 October, late; not started)

- **POST plotting into the vacuum.** Check `plot_vac_pot`, which computes the potential at external points through `calc_GH`, on the 3-D field-line vacuum, `vac%style` 1. Draw the perturbation and the potential across the plasma-vacuum interface into the vacuum. Add the vacuum energy density |∇φ|²/2μ₀ to the energy reconstruction `plot_E_rec`, so that δW_vac is drawn next to the plasma terms and its volume integral can be checked against the surface term of the free-boundary condition.
- **A Claude skill in the repository that runs PB3D for other users.** Equilibrium in; deck choices made with the convergence rules of sections 6.6, 9.2 and 11.5; run locally or through Slurm; eigenvalue and POST plots out, with a report of what converged.

### 12.2 Catalogue items still open

| Item | What it would check | Status and need | Priority in `benchmarking_prep.md` |
|---|---|---|---|
| C6, GPEC's 3-D VACUUM: Halpern, Logan, Nührenberg, Paul, Paz-Soldan, arXiv:2610.00536; Julia GPEC | the 3-D response matrix in the basis e^{i(mθ − nζ)}, element by element; order 4; 1 % from CAS3D on W7-X | basis, angle convention and (2π)² factors must be matched; circular torus first, then W7-X; 2 to 4 days | 8 |
| C7, CASTOR3D | 3-D free-boundary eigenvalues on the rippled case; a fourth vacuum method | request to IPP Garching after the local checks; TERPSICHORE now gives the number to stand next to | 10 |
| C9, DESC `singular_integral` | G·σ, H·φ, the Neumann-to-Dirichlet map and ∫Φ ∂Φ/∂n on PB3D's own normal data; 6e−7 on the W7-X constant identity | Python, `pip install desc-opt`; 1 to 2 days | 4 |
| L11, COBRAVMEC | PB3D's n → ∞ limit at most COBRA's worst field line in 3-D, equal in axisymmetry; fit ω²(n) = c1 + c2/n with n = 10 to 160 | COBRAVMEC from the STELLOPT clone, not checked to be built; 2 to 3 days | 5 |
| ELITE on the cut cbm18a | third opinion on the pedestal mode against MISHKA-1 | request text ready in `benchmarking_prep.md`, section 2c; Toon's to send | 9 |
| higher-pressure W7-X | a 3-D stellarator benchmark with a discrete mode | section 11.4 | 6 |
| TERPSICHORE at n = 20 on the cut case | radial convergence | section 11.8 | |
| published cbm18 curves, C4, and Pamela 2020's Zenodo set | γ(n) shape and first unstable n; absolute values after ω_A is pinned | PDF and data downloads need Toon's approval | 2, 3 |

Other catalogue items not run: T3 eigenfunction shapes, T7 rippled-cylinder δ² coefficient, T8 and L11 COBRA limits, T9 peeling criterion, L5 Hellmann-Feynman slope, L6 natural-boundary residual, O1 the JCP 2017 curve, O7 dominant m per n, C5 the axisymmetric VACUUM matrix, and the vacuum-side solid-angle option 3c of step 3.

### 12.3 Out of scope, in the plan's order of value

1. Negative-flux VMEC files such as W7-X and GIGA.
2. Step 8: the 3-D equilibrium path against the analytic equilibria of Landreman and of Issan, Citrin, Landreman and Tracey, through a DESC-written wout.
3. Coupling of n ± N_fp, the question of the GIGA inquiry.
4. The cause of the slow radial convergence; section 9.2 attributes it to the rational surfaces.
5. The toroidal flux in `store_vac`; G and H kept for restarts.
6. STRUMPACK's compression on 3-D kernels, a matter of speed only.

### 12.4 Decisions and actions that are Toon's

- Review and merge PR #5, with the merge strategy of his choice.
- Send the TERPSICHORE parity report to EPFL, the GCC bug report on gfortran 14, and the PSPLINE overrun report upstream.
- Send the ELITE and IPP requests once the local checks are done, by rule 8.
- Approve the downloads of the cbm18 papers and the Pamela 2020 data set.
- Decide whether any of the exploration options for the 3-D vacuum's cost (sections 11.9 and 12.5) is worth trying. None is planned or required.

### 12.5 Exploration: making the 3-D vacuum cheaper (2 November; options, nothing decided)

These are directions to explore, not planned work, and none is needed for the branch. Toon decides whether any is tried. An option would replace the present method only if it beats it on the exact potentials of `check_vac`, on the circular twin and on W7-X, and then on the twin's eigenvalue against TERPSICHORE.

**Where the error sits.** The surface integrals of `vac%style` 1 have three parts:

1. the singular cell, integrated analytically (section 6.1);
2. the near field, the intervals within four interval lengths of the influence point, where 1/r is large but finite: since 31 October Gauss–Legendre along the line, but not across lines (section 5.2);
3. the far field, a trapezoidal sum.

On a smooth closed surface the trapezoidal rule on a doubly periodic grid converges faster than any power of the spacing for smooth integrands, so the order of such a method is set by how the singular and near-singular parts are corrected [general property of the rule; not measured here]. PB3D's (α, θ) field-line grid is not doubly periodic: one poloidal turn at fixed α ends on another line, shifted by 2πq. That is why the method needs the seam and the bridge between lines. Section 11.9 shows that the solve for oscillating data and the projection of the modes both contribute to the present error.

**Options.**

- **A. A doubly periodic vacuum grid.** Assemble G and H on a regular grid in the straight-field-line angles (θ, ζ), separate from the plasma's field-line grid, and still hand `res` over in mode space. The Green's function does not depend on the coordinates, so the grid need not follow the field lines; the normal data of each mode are evaluated at its points. Stellarator symmetry and the field periods could then be used, which the field-line grid cannot. Gate: the `check_vac` exact potentials at least as accurate as now for the same number of points.
- **B. A locally corrected trapezoidal rule on grid A.** Keep the trapezoidal sum everywhere and change only the weights of a small stencil around each singular point, so that the rule is exact for the singularity's known expansion. Examples: the zeta-corrected trapezoidal rules of Wu and Martinsson, the corrections of Alpert and of Duan and Rokhlin, and the polar patches of Malhotra used in BIEST and in GPEC's Julia `Vacuum`. The literature gives 4th to 10th order at the cost of the trapezoidal sum; none was tried here. Gate: the exact potentials converge at high order on the twin and on W7-X.
- **C. A Fourier basis, as in CAS3D.** Merkel, Nührenberg and Cooper, "Free-boundary ideal MHD modes in W7-AS", CRPP Lausanne report LRP 555/96, pp. 51–58 (1996) describe the free-boundary CAS3D. Its vacuum solves Green's identity on the plasma surface, 2πΦ + ∮Φ ∂G/∂n′ = ∮G ∂Φ/∂n′, for the coefficients Φ_mn of a truncated double Fourier series cos 2π(mθ + nφ) in magnetic coordinates, using stellarator symmetry. The singular kernels are regularised by subtracting analytic functions with the same singularity, whose Fourier transforms are known (Merkel 1986, the NESTOR method). The physics and the coupling to the plasma are those of PB3D: ½∮Φ ∂ₙΦ with ∂ₙΦ = |∇s|⁻¹ B·∇(ξ·∇s), an infinite vacuum, the surface term dropped by pressure balance, and a real, symmetric, non-negative matrix added to the edge block. The unknowns are a few hundred to a few thousand Fourier coefficients, against Φ at each of up to 51 000 field-line points now. NESTOR is second order on the constant-potential test (Malhotra et al. 2020), as `check_vac` is now, so the gain would be in the size of the system rather than the order. Sources: the C++ NESTOR in VMEC++, with Schilling 2026 for its singularity treatment; GPEC's Julia `Vacuum` (C6) as the element-by-element check. C is an alternative to A with B.
- **D. A solver that scales like the number of points.** The dense LU needs N² memory and N³ work for N boundary points. The operator is a second-kind integral equation, and `check_vac` measured its 1-norm condition number at 3.5 to 3.8 on every grid, W7-X included, so GMRES should converge in a number of iterations that does not grow with N [estimate]. Each iteration needs one matrix-vector product, which a fast multipole method (PVFMM, for example) or a hierarchical matrix does in O(N) to O(N log N) work and memory, without storing G and H. There are 2 `n_mod_X`, about 40, right-hand sides, so a block or recycled solver would help. The direct alternative: order the points by a quadtree on the surface and compress with HODLR or H² matrices (STRUMPACK's HODLR through ButterflyPACK, or H2Lib), roughly O(N log² N). STRUMPACK's compression failed in July (relative residual 175 at 3660 points, section 5.6); the field-line ordering, which scatters neighbouring points, and the hard-coded sampling parameters may be the reasons [not checked]. Gate: the same `res` with memory that grows like N.

**Richardson.** A separate vacuum grid need not give up the Richardson levels. Either the vacuum grid doubles with each level (a doubly periodic grid nests as well as the field-line grid, so the interlaced reuse still works, and a compressed or fast-multipole operator is rebuilt), or the vacuum is made accurate enough that its error stays far below the plasma's at every level and is computed once. Which fits the Richardson scheme better is Toon's call.

**Order, if any is tried.** A, then B, then D, each gated by `check_vac` before any stability run; or C in place of A and B.

**Guesses, not measurements.** B alone might cut the point count tenfold. B with D might bring 51 000 points and 220 GB down to 5 000 to 10 000 points and a few GB.

## Appendix A. File map of the branch

From `git diff --stat master...HEAD`: 140 files, 17 668 lines added and 1255 removed. One line per file; files that changed together for one purpose share a line.

### A.1 Code

| File | What changed |
|---|---|
| `Modules/vac_ops.f90` | exterior solve with a shift that depends on `vac%style`; `vac%style` 1 assembly, weights, bridge, row-sum H diagonal, G diagonal from the linear harmonics; α spacing; zeroing before accumulation; inward normal in `store_vac_VMEC`; −2π factor and Hermitian part in `calc_vac_res`; `check_vac_exact`; condition estimate; ScaLAPACK fallback |
| `Modules/vac_utilities.f90` | `calc_GH_int_1` rewritten: exact singular half-cell, dipole far field, near-field Gauss-Legendre on a four-point stencil; distributed vector helpers moved in |
| `Modules/vac_vars.f90` | STRUMPACK under `PB3D_WITH_STRUMPACK`; `BLACSCTXTSIZE` public |
| `Modules/SLEPC_ops.f90` | PETSc 3.25 port; `BC_style` 3 removed; `BC_style` 2 documented; `BC_style` 4 edge block assembled once; free-edge weight `w_q` in `fill_mat` |
| `Modules/SLEPC_utilities.f90` | PETSc 3.25 port of `insert_block_mat`: logical flags, `InsertMode`, flattened blocks for `MatSetValuesBlocked` |
| `Modules/VMEC_ops.f90` | B_0 as the (0,0) mode of \|B\| on the axis; Jacobian, \|B\| and B_i on the full Nyquist list; `mpol_tf`, `ntor_tf` |
| `Modules/VMEC_utilities.f90` | `calc_trigon_factors` in basis form, grid points first; `fourier2real` for any mode list; angle version calls the table version |
| `Modules/VMEC_vars.f90` | `mn_V_nyq`, `mnmax_nyq_V`, `mpol_tf`, `ntor_tf` |
| `Modules/eq_ops.f90` | odd spline condition for the HELENA shear; J, \|B\| and B_i reconstructed on the Nyquist list; spline imports |
| `Modules/X_ops.f90` | `setup_modes` slots from the mode number and the consistency check; spline imports |
| `Modules/MPI_utilities.f90` | `redistribute_var` takes every point from its owner; `lock_return_acc` evaluates its impure call separately |
| `Modules/grid_ops.f90` | 20-point margin in `redistribute_output_grid`; spline imports |
| `Modules/grid_utilities.f90` | index order of `trigon_factors`; spline imports |
| `Modules/grid_vars.f90` | `trigon_factors` documented as the trigonometric basis |
| `Modules/input_ops.f90` | input `check_vac`; default `BC_style = [1,2]`; free 3-D runs need `alpha_style` 2; aliasing warning; `BC_style` 3 removed |
| `Modules/num_vars.f90` | `check_vac`; `prog_version` 2.47 → 2.48; `atan` for `datan` |
| `Modules/MPI_ops.f90` | broadcasts `check_vac` |
| `Modules/PB3D_ops.f90` | reads the Nyquist list and the basis extents back for POST; refuses output files older than 2.49 |
| `Modules/sol_ops.f90` | POST's vacuum energy guard and "NO VACUUM" print removed; `allocated(vac%res)` guard; index order of `trigon_factors` |
| `Modules/driver_POST.f90` | reads the vacuum only for free-boundary runs |
| `Modules/driver_eq.f90` | stores the vacuum only for free-boundary runs; `ltest` under `ldebug` |
| `Modules/driver_sol.f90` | skips the vacuum operations for a fixed boundary |
| `Modules/HELENA_ops.f90` | spline imports; `invert_top_bottom_H` under `ldebug` |
| `Modules/spline_utilities.f90` | new: the PSPLINE wrappers split out of `num_utilities` |
| `Modules/num_utilities.f90` | spline wrappers moved out, 390 lines fewer |
| `Modules/var_1D_vars.f90` | new: `var_1D_type` without HDF5 |
| `Modules/HDF5_vars.f90` | `var_1D_type` re-exported from `var_1D_vars` |
| `Modules/HDF5_ops.f90` | a debug-only use under `#if ldebug` |
| `Modules/HDF5_utilities.f90` | one declaration line reformatted |
| `Modules/PB3D_utilities.f90` | `var_1D_type` from `var_1D_vars` |
| `Modules/sol_utilities.f90` | spline imports; `ldebug` guards |
| `Modules/dtorh.f90` | type-specific intrinsics such as `DSQRT` replaced by generic ones |
| `Modules/num_ops.f90` | `output_ops` only under `ldebug` |
| `Modules/messages.f90`, `Modules/files_ops.f90` | imports of memory variables moved |
| `Modules/test.f90` | spline import |
| `Libraries/build_pspline.sh` | new: builds PSPLINE and applies the overrun patch |
| `Libraries/pspline-r8genxpkg-nx3.patch` | new: the guard for three-point splines |
| `Libraries/build_libstell_min.sh` | new: minimal LIBSTELL with `read_wout_mod` |

`prog_version` in `num_vars.f90` reads 2.48 while the README's newest entry is 2.49.

### A.2 Tests

| File | What changed |
|---|---|
| `tests/CMakeLists.txt` | new: all layers, MPI rows, opt-in MPI regression rows, fixture decompression |
| `tests/fullstack/main.f90`, `fullstack_utils.f90` | new: driver and the gathered result helpers |
| `tests/fullstack/test_vac_kernels.f90` | new: `vac%style` 2 interval kernels |
| `tests/fullstack/test_vac_greens.f90` | new: axisymmetric G and H, solve round trip, cylinder response |
| `tests/fullstack/test_vac_3d.f90` | new: the six 3-D and exact-solution tests of section 9.1 |
| `tests/fullstack/test_splines.f90` | new: converted legacy spline check |
| `tests/fullstack/test_calc_int_vol.f90` | new: converted legacy volume-integral check |
| `tests/fullstack/test_read_HEL.f90` | new: HELENA golden file |
| `tests/unit/main.f90` | new: unit driver |
| `tests/unit/test_dtorh.f90` | new: toroidal harmonics, including n = 20 to 80 |
| `tests/unit/test_str_utilities.f90`, `test_files_utilities.f90` | new: string and file utilities |
| `tests/unit/test_num_utilities.f90`, `test_num_ops.f90` | new: numerics and zero finders |
| `tests/unit/test_X_utilities.f90`, `test_PB3D_utilities.f90` | new: mode-index and Richardson index arithmetic |
| `tests/regression/cbm18a_fixed_bnd.sh`, `input_cbm18a_fixed_bnd` | new: fixed anchor and deck |
| `tests/regression/cbm18a_fixed_bnd_np.sh` | new: one process against `mpirun` |
| `tests/regression/cbm18a_free_bnd.sh`, `input_cbm18a_free_bnd` | new: free anchors for `BC_style` 2 and 4 |
| `tests/regression/cbm18a_shear_HEL.sh`, `check_shear_HEL.py` | new: shear against the direct formula |
| `tests/regression/check_vac_res.py` | new: response Hermitian and positive semi-definite |
| `tests/regression/circ_eps005_free_bnd.sh`, `input_circ_eps005` | new: circular kink against MISHKA-1, 400-point order check |
| `tests/regression/circ_eps033_vmec_twin.sh`, `input_circ_eps033` | new: HELENA against the VMEC twin, fixed and free |
| `tests/vendor/test-drive/testdrive.F90`, `LICENSE-MIT`, `README.md` | new: vendored test framework |

### A.3 Fixtures and references

| File | What it holds |
|---|---|
| `tests/fixtures/cbm18a_extended_vac.12.xz` | cbm18a extended mapping, 5.2 MB |
| `tests/fixtures/circ_eps033_q104.12.xz`, `circ_eps033_q104.fort.10` | circular ε = 0.33 mapping and HELENA input |
| `tests/fixtures/circ_eps005_q104.12.xz`, `circ_eps005_q104.fort.10` | circular ε = 0.05 mapping and HELENA input |
| `tests/fixtures/circ_eps033_q104.vmec.input` | VMEC INDATA of the twin |
| `tests/fixtures/wout_circ_eps033_q104.nc.xz` | the VMEC twin |
| `tests/reference/mishka_circular.txt` | MISHKA-1 values against `RWALL` |
| `tests/reference/cyl_kink.py` | the cylinder formula, sympy |
| `tests/reference/setup_modes_sim.py` | the slot simulation |
| `tests/reference/make_vmec_twin.sh`, `vmec_twin.py` | the twin recipe |
| `tests/reference/gfortran14_bounds_repro.f90` | the compiler-bug reproducer |
| `tests/reference/cbm18a_cut/README.md` | the cut equilibrium and its results |
| `tests/reference/cbm18a_cut/cbm18a_cut.fort.10` | its HELENA input |
| `tests/reference/cbm18a_cut/input.cbm18a_cut` | its VMEC INDATA |
| `tests/reference/cbm18a_cut/foot_surface.py`, `make_cut_input.py`, `parse_hel12.py` | the scripts that made it |

### A.4 Scripts

| File | What it does |
|---|---|
| `scripts/snellius/env.sh` | environment |
| `scripts/snellius/build_deps.sh`, `build_pb3d.sh` | dependency and PB3D builds |
| `scripts/snellius/build_pspline_linux.sh`, `build_libstell_min_linux.sh` | Linux variants of the library scripts |
| `scripts/snellius/pipeline.sh` | unattended build and run sequence |
| `scripts/snellius/run_deck.sbatch`, `mem_run.sbatch` | deck runs; timed runs |
| `scripts/snellius/compare_eq.py`, `compare_np.py`, `compare_petsc_mats.py` | equilibrium, process-count and matrix comparisons |
| `scripts/snellius/compare_res.py`, `compare_res_series.py`, `checkvac_table.py` | response and `check_vac` tables |
| `scripts/terpsichore/wout2terps.py`, `terps_data.py`, `README.md` | TERPSICHORE input and recipe |
| `Documentation/figures/step0/make_figures.py` | the step-0 figures |

### A.5 Documentation, build and CI

| File | What changed |
|---|---|
| `Documentation/plan_vacuum_bugfix.md` | new: plan and log of issue #4 |
| `Documentation/vacuum_work_overview.md` | new: the outcome |
| `Documentation/benchmarking_prep.md` | new: benchmark plan and priorities |
| `Documentation/benchmarking/lit_A_axisymmetric.md`, `lit_B_3d_codes.md`, `lit_C_vacuum_bem.md` | new: the three literature surveys |
| `Documentation/testing.md` | new: test layers and anchors |
| `Documentation/quality-and-ci.md` | new: quality and CI strategy |
| `Documentation/snellius-setup.md` | new: Snellius recipe |
| `Documentation/spack-setup.md` | STRUMPACK-Dense 1.1.1 recipe updated |
| `Documentation/figures/step0/*.png` | four step-0 figures |
| `README.md` | changelog entries 2.48 and 2.49 |
| `AGENTS.md`, `openwiki/` | new: pointer to and pages of the OpenWiki documentation |
| `CMakeLists.txt` | `PB3D_BUILD_EXECUTABLES`, optional STRUMPACK, gfortran 14 bounds workaround, `-fcheck=no-bounds` for `dfftpack`, core library |
| `.github/workflows/tests.yml` | new: the eight CI jobs |
| `.github/warnings-baseline.txt`, `.github/doxygen-warnings-baseline.txt` | new: 74 and 374 |
| `Doxyfile`, `fortitude.toml` | new: Doxygen and linter configuration |
| `spack.yaml` | new: pinned PETSc and SLEPc |
| `.gitignore` | `spack.yaml` allowed; build directories and the LIBSTELL build ignored |

## Appendix B. Commits by step

Seventy-three commits. Dates from `git log`.

| Step or topic | Commits | Dates |
|---|---|---|
| 2.48: testing infrastructure, PETSc port, build scripts, vacuum skip, regression layer, docs | a8c07b6, 8a85ac2, d32ef27, 775219a, d4812f7, ad2313b, da1a74f, db4427d, 78fb49e | 6 to 8 July |
| cloud session: build and vacuum | 1862df0, 429f003, 31fe4a7, 3679ee7, 89996a1, d6d49e6 | 8 and 9 July |
| cloud session: CI and regression | c16015c, 9f9a69c, 0f5b27e | 9 July |
| cloud session: boundary conditions | cc4c96d, 41b5274, 1d74b8f | 10 July |
| cloud session: quality and tests | 8ac37d7, 9db1d9e, 6f0f8d5, 1c7283c, fa3560b, 6ce5029, 17716c4, a53869c, 88153e4, 17b712c | 10 July |
| plan | 70e6de8, 2058627 | 27 October |
| step 0: references | 6d2ee51 | 27 October |
| steps 1 to 5 | 28241e1, 28f8dab, 212e3ec, 83c02df, e4e3fdd | 27 October |
| step 6a to 6d | 50c9511 | 27 October |
| step 0: figures; CI log artifact | d4eceb7; 2d1363b | 29 and 30 October |
| step 6e: twin, gate failure, normal fix | a5c8c41, dfe439e, e3b4847 | 30 October |
| steps 6f and 7: Nyquist matching, B_0, changelog | 0d49f19 | 30 October |
| plan for step 9 | a7a8f6a | 30 October |
| open items: `redistribute_var`, MPI rows, recorded ratios | 798fe50 | 30 October |
| step 9 and its follow-up | 3b99109, 329b7a4 | 30 October |
| log of the night; Snellius scripts | 6b818ad; 9b91567, d6186ef | 30 October |
| catalogue weekend: `check_vac`, T4, T11, helpers, condition estimate, per-point files | 6894199, 612fc1f, fe305fd, 13bdef4, a4cb73d, 83724cc, 10b1291 | 31 October |
| `fourier2real` consolidation | 1a53f50, 96c140e | 31 October |
| field-line BEM: near field and G diagonal | 0bbda47, 5f82da8, c9f597e, d8a1c6c | 31 October |
| free-edge weight, MISHKA-1 on cbm18a, literature | 4915d4b | 31 October |
| TERPSICHORE: build, parity error, rippled twin | be13baa, b19d462 | 31 October |
| convergence study | 07bd613 | 31 October in git; the log dates it after midnight |
| overview; cut cbm18a; overview update | 12fdec9, 3073a7f, 795c695 | 1 November |

## Appendix C. Glossary of PB3D inputs and options

| Input | Meaning |
|---|---|
| `vac%style` | Not an input: the vacuum method, set from `eq_style`. 1: 3-D boundary elements on the field-line grid (VMEC, `calc_GH_1`). 2: axisymmetric boundary elements with toroidal harmonics for one n (HELENA, `calc_GH_2`). |
| `BC_style` | Two integers, left and right boundary. 1 is a fixed boundary. 2 is the variational free boundary, the default since step 4. 4 is the explicit natural-boundary-condition row, not validated. 3 was removed in July. |
| `X_style` | 1: prescribed mode window from `min_sec_X` to `max_sec_X`. 2, the default: a fast window of `n_mod_X` modes that follows nq radially; not a nested trial space (section 11.5). |
| `prim_X` | The primary mode number: n with the poloidal flux label, m with the toroidal one. |
| `n_mod_X` | Number of secondary modes, m with the poloidal flux label. |
| `min_sec_X`, `max_sec_X` | Limits of a prescribed window, `X_style` 1. |
| `alpha_style` | 1: one field line over several turns, the HELENA default. 2: several field lines over one turn, the VMEC default with 10 lines; required for free-boundary VMEC runs. |
| `n_alpha` | Number of field lines for `alpha_style` 2. The 3-D vacuum needs at least 2; 128 lines with about 36 points per wavelength give 0.1 % (section 6.6). |
| `min_alpha`, `max_alpha` | Range of the field-line label in units of π; default 0 to 2. The last line is excluded. |
| `alpha` | Field-line label for `alpha_style` 1; has no effect with `alpha_style` 2. |
| `min_n_par_X` | Number of points along a field line on the first Richardson level; 100 in the regression decks. |
| `min_par_X`, `max_par_X` | Parallel range in units of π; default 0 to 2, one poloidal turn. |
| `n_r_sol` | Number of radial points of the solution grid; 200 in the regression decks. |
| `min_r_sol`, `max_r_sol` | Radial window of the solution in normalised flux. |
| `use_pol_flux_F` | Poloidal flux as the normal label; the toroidal label is untested. |
| `check_vac` | Logical, default off. Tests the 3-D vacuum operators on the real surface with exact potentials; VMEC equilibria only (section 5.3). |
| `EV_guess` | Target eigenvalue for SLEPc's shift-and-invert, default −0.3. PB3D returns the eigenvalue nearest to it. |
| `n_sol_requested` | Number of eigenvalues requested; 1 in the decks. |
| `tol_SLEPC` | Solver tolerance; 1e−6 in the decks. |
| `max_it_rich` | Number of Richardson levels; 1 means no extrapolation. |
| `magn_int_style` | Parallel integration rule: 1 trapezoidal, 2 Simpson 3/8; the decks use 2. |
| `norm_disc_style_sol` | Radial discretisation: 1 central, 2 left differences, the default. |
| `norm_style` | 1: MISHKA normalisation, used throughout. 2: a modified B_0 and ψ_0, with EV_COBRA = EV_PB3D · β/2. |
| `K_style` | 1: perpendicular kinetic energy only. |
| `X_grid_style` | Normal grid of the perturbation: 1 equilibrium, 2 solution, 3 enriched. |
| `plot_vac_pot` | POST: plot the vacuum potential of the computed eigenmodes outside the plasma, both vacuum styles (6 November 2026, checked against the exact exterior potentials of `check_vac`); box `min_Rvac_plot` … `max_Zvac_plot` in metres (default from the boundary), `n_vac_plot`, planes as the plasma plots. |
| `--export_HEL` | Command-line option: PB3D writes a VMEC INDATA file from a HELENA mapping; used for the twins. |
| `keep` | Sixth argument of `run_deck.sbatch`: keep `PB3D_out.h5` for the response tools. |
| `PB3D_TEST_MPI_REGRESSION` | CMake option for the opt-in 2- and 4-process regression rows. |
| `PB3D_BUILD_EXECUTABLES` | CMake option; OFF builds only the dependency-light core and its unit tests. |
