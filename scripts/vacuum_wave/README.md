# Wave-amplitude vacuum: a prototype (issue #4, exploration only)

Nothing here is decided or required. This tests one idea from issue #4: let a Richardson level of the 3-D vacuum keep the previous level's integrals, as the plasma does in `calc_magn_ints`, instead of N x N matrices over the boundary points.

## The idea

PB3D's field-line vacuum (`vac%style` 1) solves for the potential at each of the N boundary points (N = field lines x parallel points) with two dense N x N matrices, G and H. The prototype keeps the same rows of G and H, but writes the potential as a sum of M' waves exp(i(n φ − m θ)) and tests the boundary equation against the same waves (Galerkin). In the notation of section V of the 2019 manuscript:

  (E′† I H E′) a = E′† I G q,  Φ = E′ a,

where E′ holds the waves at the points (like E of eq. 106, with more modes) and I the quadrature weights (eq. 100). Every entry of the two small matrices is a sum over pairs of points. H and G are built one block of rows at a time and added straight into the small matrices, so no N x N matrix exists.

A Richardson level (parallel points n_par → 2 n_par − 1) halves the weight of every old point, so a pair of old points counts ¼ of what it counted before:

  new sums = ¼ × old sums + the pairs that contain a new point.

The self terms come from the interior identities, as in `calc_GH_1`: H's from the constant potential, G's from the linear potentials x, y, z. They need the row sums of each point, which carry over in the same way (½ × old + new columns), a few numbers per point.

## Test case

- Circular torus, R0 = 3, a = 1 (aspect ratio 3, like the circular VMEC twin).
- Field lines φ = α + q θ, one poloidal turn per line, q = 1.04296875. This q is chosen so that no line end meets a line start for 32, 64 or 128 lines.
- Exact exterior harmonics with n = 10 and m = 6, 11, 15 in toroidal coordinates (η, u, φ): sqrt(cosh η − cos u) P^n_{|m|−1/2}(cosh η) exp(i(n φ − m u)).
- Both methods get the exact normal derivative q = ∂Φ/∂n on the torus and solve for Φ. Errors are measured against:
  - the exact Φ, as a relative weighted L2 norm;
  - the exact vacuum energy −∫ q̄ Φ dS, computed on a periodic (θ, φ) grid.
- The waves are exp(i(n φ − m θ)) with θ the geometric poloidal angle. m runs over the test window 6 to 15, widened by k sidebands on each side.
- No near-field rule: the kernels use the plain trapezoidal rule along and across the lines, with the identity self terms. Both methods use exactly the same rows, so the comparison isolates the change of unknowns.

## Success criteria, fixed before the runs

1. The carried-over sums equal the ones computed from scratch to 1e−12.
2. With enough waves, the wave form is as accurate as the point-value solve on the same grid (error at most 1.5 times larger).
3. The number of waves stays a small multiple of the plasma's mode window.
4. The peak memory of the wave form stays far below one N x N matrix.

## Run

    uv sync
    .venv/bin/python scripts/vacuum_wave/prototype.py waves     # sidebands needed
    .venv/bin/python scripts/vacuum_wave/prototype.py grid      # parallel points, carry-over, point-value solve
    .venv/bin/python scripts/vacuum_wave/prototype.py lines     # field lines
    .venv/bin/python scripts/vacuum_wave/prototype.py memory 801

A short report in the layout of the 2019 vacuum manuscript, with figures, is `Documentation/vacuum_wave_report/report.pdf`. Its numbers come from `figdata.json`, written by `figures.py data` (about 6 minutes); `figures.py plot` draws the figures, and `latexmk -pdf report.tex` in that directory builds the PDF. The report adds two results to the ones below:
- with single waves as data, as PB3D's plasma modes are, the potential stays on the data wave (neighbours 0.1 to 6 %, three waves away 0.1 % or less);
- the wave-form matrix is nearly tridiagonal: its first neighbours are 0.17 of the diagonal, close to a/(2 R0) = 1/6.

## Results (3 November 2026, Mac, one process)

All four criteria are met. The exact harmonics are harmonic to the accuracy of the finite-difference check: relative Laplacian 2e−6 or less. Errors below are listed for m = 6, 11, 15; "Φ" is the relative weighted L2 error of the potential, "energy" the relative error of −∫ q̄ Φ dS.

**Waves needed (criterion 3).** 64 lines × 201 points. With k sidebands on each side of the window 6 to 15:

| k | waves | Φ | energy |
|---|---|---|---|
| 0 | 10 | 0.71, 0.23, 0.67 | 0.40, 0.063, 0.55 |
| 2 | 14 | 0.066, 0.054, 0.56 | 1.0e−3, 9.3e−3, 0.40 |
| 6 | 22 | 3.2e−3, 1.4e−2, 0.18 | 3.1e−3, 1.4e−2, 1.6e−2 |
| 12 | 34 | 3.2e−3, 1.4e−2, 3.4e−2 | 3.1e−3, 1.4e−2, 3.3e−2 |
| 20 | 50 | 3.2e−3, 1.4e−2, 3.4e−2 | 3.1e−3, 1.4e−2, 3.3e−2 |

m = 6 and 11 are settled at 6 sidebands, m = 15 at 12. That is 34 waves for a window of 10 modes. On an axisymmetric torus n is conserved, so only poloidal sidebands appear.

**Carry-over and accuracy (criteria 1 and 2).** 64 lines, k = 12. Each level's sums were computed from scratch and carried over from the previous level:

| points per line | N | carried over vs from scratch | time, carried / scratch | waves: Φ | points: Φ | energy (both) |
|---|---|---|---|---|---|---|
| 51 | 3 264 | | | 0.028, 0.10, 0.22 | 0.028, 0.10, 0.22 | 0.028, 0.11, 0.22 |
| 101 | 6 464 | 5.9e−16 | 0.7 s / 0.8 s | 3.0e−3, 2.2e−2, 5.9e−2 | 3.4e−3, 2.3e−2, 6.1e−2 | 2.7e−3, 2.2e−2, 6.0e−2 |
| 201 | 12 864 | 4.7e−16 | 2.4 s / 3.2 s | 3.2e−3, 1.4e−2, 3.4e−2 | 3.6e−3, 1.5e−2, 3.5e−2 | 3.1e−3, 1.4e−2, 3.3e−2 |
| 401 | 25 664 | 6.2e−16 | 9.8 s / 13.4 s | 3.1e−3, 1.4e−2, 3.3e−2 | (not run: 5 GB per matrix) | 2.9e−3, 1.4e−2, 3.2e−2 (waves) |

The wave form is as accurate as the point-value solve, and slightly better for Φ. The energy errors agree to the digits shown.

**Field lines.** 201 points per line, k = 12. The error stops falling with the parallel points from 201 on; the number of lines sets it here (no near-field rule):

| lines | N | waves: energy | points: energy |
|---|---|---|---|
| 32 | 6 432 | 0.025, 0.12, 0.31 | 0.025, 0.12, 0.30 |
| 64 | 12 864 | 3.1e−3, 1.4e−2, 3.3e−2 | 3.1e−3, 1.4e−2, 3.3e−2 |
| 128 | 25 728 | 3.5e−4, 2.9e−3, 8.3e−3 | (not run) |

**Memory (criterion 4).** 64 lines × 801 points, N = 51 264: 50 s, peak memory of the whole process 1.9 GB, against 21 GB for one N × N matrix of doubles. PB3D's 128 × 400 run, of the same size, needed 220 GB. Most of the 1.9 GB is the block of 256 rows that the kernel is evaluated on; smaller blocks need less. Errors at this size: Φ 3.0e−3, 1.4e−2, 3.3e−2; energy 2.9e−3, 1.4e−2, 3.2e−2.

## Not tested

- A near-field rule. PB3D's rule (Gauss–Legendre on a cubic interpolant within four intervals along the line) changes the near pairs when points are added. Those pairs would need their own table, rebuilt at every level; the far pairs carry over.
- A stellarator. The waves then need the toroidal families n ± N_fp as well, so M′ grows; by how much is not known.
- PB3D's own data. In PB3D the normal data are the plasma's modes (`EP` in `calc_vac_res`), so the right-hand side becomes an M′ × M table, carried over in the same way.
- Flux angles. The waves here use the geometric θ; in PB3D they would use the straight-field-line angle.
- Fortran, MPI and the time per level. The pair sums still cost N² per level, of which the carry-over saves a quarter.
