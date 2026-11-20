# Radial schemes on two realistic cases

Measured on 20 November 2026 with PB3D at f9713d7 (one process, Mac) and with `pb3d_coefficients.py`. This note corrects the summary of `../plan_radial_finite_elements.md`: the error proportional to $1/N$ that the plan starts from belongs to the regression decks, whose radial range (0.45 to 0.75) cuts the eigenfunction at both ends. Runs whose range leaves the eigenfunction free do not have it.

$N$ is the number of radial points (`n_r_sol`). "Order" is `norm_disc_prec_sol` of the present finite differences. All errors are relative errors of the lowest eigenvalue $\lambda = \omega^2/\omega_A^2$.

## The two cases

| | cbm18a, interior ballooning | peeling-ballooning, free boundary |
|---|---|---|
| Input | `inputs/input_cbm18a_range_015_085` | `inputs/input_peeling_ballooning_free`, the input of a run of 2017 (`TEST_VAC_PEELING_BALLOONING`) |
| Equilibrium (HELENA) | `cbm18a_extended_vac.12` of the fixtures | the HELENA output (`fort.12`) of that run |
| Radial range | 0.15 to 0.85 | 0.5 to 1.0 |
| Edge | fixed (the eigenfunction is $3 \times 10^{-3}$ there) | free, with the vacuum |
| Fourier modes | $n = 10$, fixed set $m = 7 \ldots 64$ | $n = 5$, fixed set $m = 12 \ldots 37$ |
| Points along the field line | 300 | 500 |
| Where the eigenfunction lives | on the pressure pedestal at 0.59; $10^{-4}$ of its largest size at the inner end | in the pedestal at 0.96 to 1.0, largest at the plasma edge; $10^{-3}$ at the inner end |
| Eigenvalue for many points | $-4.710242 \times 10^{-2}$ | $-1.79753 \times 10^{-2}$ |

The pressure of this cbm18a equilibrium is zero beyond 0.70 and its safety factor reaches 25 at the edge, so a range up to 1.0 needs more than 300 radial points before `setup_modes` accepts it. The range 0.15 to 0.85 avoids that.

## PB3D as it is

| $N$ | cbm18a, order 3 | cbm18a, order 2 | peeling-ballooning, order 3 | peeling-ballooning, order 2 |
|---|---|---|---|---|
| 200 | −3.9e-2 | −2.4e-1 | −1.5e-1 | −1.2e-1 |
| 400 | −2.6e-3 | −6.5e-2 | −2.8e-2 | −5.7e-2 |
| 800 | −1.2e-4 | −1.7e-2 | −3.3e-3 | −1.8e-2 |
| 1600 | +2e-6 | −4.1e-3 | −2.2e-4 | −4.9e-3 |
| 3200 | +3e-6 | −1.0e-3 | −1.1e-5 | −1.2e-3 |

- No term proportional to $1/N$ limits either case down to $10^{-5}$. Order 2 falls as $1/N^2$ and order 3 as about $1/N^4$.
- cbm18a keeps a positive error of about $3 \times 10^{-6}$ at 1600 and 3200 points. It fits the fixed edge at 0.85, where the eigenfunction is small but not zero.
- The archived setting of the peeling-ballooning run (1000 points, order 2) gives $-1.8189 \times 10^{-2}$, which is 1.2 % off.

## The last step of PB3D in Python, with PB3D's own coefficient matrices

`pb3d_coefficients.py` reads the coefficient matrices `PV_0..2` and `KV_0..2` that PB3D stores in `PB3D_out.h5`, interpolates them as `interp_V` does, and builds the two matrices from a table of `fe_model.py`. Everything before that step is PB3D.

**Check.** With the table of the differences, the matrices of the script equal those that PB3D writes with `-eps_view_mat0 binary:A.dat -eps_view_mat1 binary:B.dat`, entry by entry: the largest difference is $7 \times 10^{-15}$ of the largest entry for $\mathsf{A}$ and $1 \times 10^{-13}$ for $\mathsf{B}$ (free-boundary case, 200 and 400 points). The eigenvalues agree with PB3D's to 10 digits at every $N$. The vacuum matrix is the difference in the last diagonal block of $\mathsf{A}$; taken from the runs with 200 and 400 points it is the same to $2 \times 10^{-14}$.

**Result.** The same coefficient matrices with other tables. "Unknowns" is the number of unknowns per Fourier mode, which is $N$ for the differences.

Peeling-ballooning, free boundary:

| unknowns | differences order 3 (PB3D) | elements, degree 2 | elements, degree 3 | cubic Hermite | cubic B-spline |
|---|---|---|---|---|---|
| 100 | −3.5e-1 | +7.4e-1 | +3.6e-1 | | |
| 200 | −1.5e-1 | +1.1e-1 | +2.8e-2 | +8.4e-3 | +3.1e-3 |
| 400 | −2.8e-2 | +9.5e-3 | +1.3e-3 | +1.9e-4 | +2.9e-5 |
| 800 | −3.3e-3 | +6.5e-4 | +3.0e-5 | +4.7e-6 | below 1e-6 |
| 1600 | −2.2e-4 | +4.1e-5 | below 1e-6 | | |

cbm18a, interior ballooning:

| unknowns | differences order 3 (PB3D) | elements, degree 2 | elements, degree 3 | cubic Hermite | cubic B-spline |
|---|---|---|---|---|---|
| 100 | −3.7e-1 | +4.0e-2 | +1.5e-2 | | |
| 200 | −3.9e-2 | +2.9e-3 | +8.4e-4 | +1.4e-4 | +4.4e-5 |
| 400 | −2.6e-3 | +2.0e-4 | +1.9e-5 | +2.8e-6 | +1.4e-7 |
| 800 | −1.2e-4 | +1.3e-5 | +3e-7 | below 1e-7 | below 1e-7 |

- Every error of the elements is positive: the eigenvalue comes from above, as the theory says for a fixed set of Fourier modes. The differences come from below here.
- Elements of degree 3 fall as $1/N^6$, of degree 2 as $1/N^4$. With 6 Gauss points per element the numbers of degree 3 change by less than 7 %.
- For an error of $10^{-3}$ in the peeling-ballooning case, the differences of order 3 need about 1100 radial points, elements of degree 3 about 420, cubic Hermite about 300 and cubic B-splines about 240.
- At 100 unknowns no scheme is usable in the peeling-ballooning case.

## What does not hold any more

- The summary of the plan ("the elements remove an error of $6.3 \times 10^{-3}$ at 200 points"): that error exists only in a deck that cuts the eigenfunction.
- The model problems of `fe_model.py` give the powers of $1/N$ correctly, and the sign of the errors. They do not give the sizes: at 13 radial points per length of change the model has $3.5 \times 10^{-5}$ for the differences of order 3, and PB3D has $2.6 \times 10^{-3}$ on cbm18a.
- "The sliding set of Fourier modes costs $4 \times 10^{-4}$" (plan, section 3.1): the fixed-set runs it was compared with used 100 points along the field line, which is too few for a fixed set of 43 modes (0.3 % error; 8 % for a set of 60 modes; none visible with 300 points).

## Reproducing

```bash
PB3D inputs/input_peeling_ballooning_free <fort.12> -eps_view_mat0 binary:A.dat -eps_view_mat1 binary:B.dat
python -B pb3d_coefficients.py PB3D_out.h5 --compare=A.dat,B.dat,<n_r_sol>,<norm_disc_prec_sol> --vac=vac.npy
python -B pb3d_coefficients.py PB3D_out.h5 -1.797532e-2 200 400 800 1600 --vac=vac.npy
```

The run that provides `PB3D_out.h5` needs `tol_SLEPC` $= 10^{-9}$ or larger, so that PB3D stores it; its number of radial points does not matter, because the coefficient matrices are stored on the perturbation grid.
