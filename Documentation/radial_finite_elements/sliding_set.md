# The sliding set of Fourier modes with PB3D's coefficient matrices

Measured on 20 November 2026 with PB3D at f9713d7 (one process, Mac) and with `sliding_set.py`, on the two realistic cases of `realistic_cases.md`. This note answers the open point of issue #51: "the sliding set (`X_style` 2) is not tested with real coefficient matrices".

$N$ is the number of unknowns per Fourier mode (`n_r_sol` for the differences), $W$ is `n_mod_X`, the number of Fourier modes in the sliding set. All errors are relative errors of the lowest eigenvalue.

## What PB3D does with the sliding set

Read in the code and reproduced by the script (next section).

- **The set of a position.** `init_modes` takes $m_\text{min} = \max(5, \operatorname{nint}(n q - (W-1)/2))$ at the equilibrium points. `setup_modes` interpolates it linearly to the positions of a grid and rounds it; the set is $m_\text{min} \ldots m_\text{min} + W - 1$. Every Fourier mode keeps one slot, $\operatorname{mod}(m - m_\text{first}, W)$, for as long as it is in the set.
- **Coefficient matrices on the solution grid.** `interp_V` interpolates every pair of modes $(m_k, m_l)$ over the perturbation points where both modes are in the set, to the solution points where both are in the set. The spline is cubic (not-a-knot) for 4 points or more, a straight line for 2 or 3, a copy for 1. Outside the input points it extrapolates with the value, the first derivative and half the second derivative at the end point (`calc_extrap` of `spline_utilities`); it does not use the cubic end polynomial.
- **Assembly.** `fill_mat` builds the block of an energy point with the set of that point. `insert_block_mat` drops every entry whose mode is not in the set of the unknown it lands on. So an entry between unknown $(a, m_k)$ and unknown $(b, m_l)$ at energy point $g$ is used when $m_k$ and $m_l$ are both in the set of $g$, $m_k$ is in the set of point $a$ and $m_l$ in the set of point $b$.

**Check.** With the table of the present differences and PB3D's sets, the matrices of the script equal those that PB3D writes with `-eps_view_mat0 binary:A.dat -eps_view_mat1 binary:B.dat`, entry by entry, with the same pattern of non-zero entries:

| run | largest difference in $\mathsf{A}$ | in $\mathsf{B}$ | eigenvalue PB3D | script |
|---|---|---|---|---|
| peeling-ballooning, free, $W = 12$, $N = 200$ | $3.8 \times 10^{-14}$ | $1.5 \times 10^{-13}$ | −1.807758251297e-2 | −1.807758251302e-2 |
| cbm18a, fixed, $W = 20$, $N = 200$ | $6.1 \times 10^{-16}$ | $8.6 \times 10^{-16}$ | −4.896212639671e-2 | −4.896212639681e-2 |

(relative to the largest entry). The vacuum matrix that PB3D adds at the free edge is `vac/res` of `PB3D_out.h5` to $1.7 \times 10^{-14}$; the earlier statement that it is not in the output file is wrong.

**Two facts that make the test simple.**

- The coefficient matrices of a pair of modes do not depend on which other modes are in the set. A run with a fixed set and the sliding-set run of the same case give the same numbers for every pair they share, at every perturbation point, to $1.5 \times 10^{-16}$ (peeling-ballooning, 580 points, 62 321 pairs) and $1.8 \times 10^{-16}$ (cbm18a, 430 points, 152 625 pairs).
- The vacuum matrix behaves the same way: the entries of the sliding-set run equal those of the fixed-set run for the same modes exactly.

So one run with a wide fixed set gives the coefficients of every pair at every position, with no extrapolation, and the vacuum matrix for whatever set the edge unknown has. The series below use such runs: $m = 5 \ldots 48$ for the peeling-ballooning case and $m = 5 \ldots 75$ for cbm18a; every sliding set of the series lies inside them. The same series with the coefficients of the sliding-set runs themselves, interpolated and extrapolated as `interp_V` does, are in the section on extrapolation.

The Python version also reproduces PB3D's own series with the sliding set: with the coefficients of the sliding-set run and the table of the differences, its eigenvalues equal those of PB3D runs at 200 to 1600 points to about $10^{-10}$ (for example −1.8067979777e-2 by PB3D and −1.8067979781e-2 by the script, peeling-ballooning, $W = 16$, $N = 1600$).

## The two rules for elements

The set of an element is the set at the middle of the element; its energy points use that set. A shared unknown needs a rule (plan, section 3.3):

- **Rule A.** An unknown has the set of the innermost element that uses it. For elements with values at the points, a point shared by two elements has the set of the element on its inner side. A mode that leaves the set going outward ends with a free value at the last point of its last element, and its energy beyond that point is left out. A mode that enters starts at zero.
- **Rule B.** An unknown has only the modes that every element using it has. Every mode is then a continuous function that is zero outside its range, and the matrices hold its complete energy. The eigenvalue lies above the one with all modes.

The same two definitions are used for cubic Hermite functions (two unknowns per point) and cubic B-splines (each function spans four elements).

## Reference values

- **Peeling-ballooning.** The fixed set $m = 5 \ldots 48$ gives $-1.842558 \times 10^{-2}$ (B-splines with 1600 and 2400 unknowns and elements of degree 3 with 3199 unknowns agree to $10^{-8}$). The fixed set of the deck, $m = 12 \ldots 37$, gives $-1.79753 \times 10^{-2}$, which is $2.45 \times 10^{-2}$ above. That set is not complete: the eigenvector involves $m = 9 \ldots 38$, and at $\psi = 0.975$ (where $nq = 24.9$) the modes above $10^{-2}$ of the largest one there are $m = 20 \ldots 32$.
- **cbm18a.** The fixed set $m = 5 \ldots 75$ gives $-4.710420 \times 10^{-2}$ (B-splines, 800 unknowns, to $10^{-7}$). The fixed set of the deck, $m = 7 \ldots 64$, gives $-4.710242 \times 10^{-2}$, $3.8 \times 10^{-5}$ above.

All errors below are relative to these two values (the complete fixed sets), not to the decks' fixed sets.

## The error of the sliding set itself

The relative error that is left when the radial resolution is high: elements of degree 3 at 3199 unknowns per mode (rule A and rule B), and PB3D's differences of order 3 at 3200 points.

Peeling-ballooning (free boundary, $n = 5$):

| $W$ | rule A | rule B | PB3D's differences, 3200 points |
|---|---|---|---|
| 8 | | | $-2.8 \times 10^{-4}$ at 200 points, $-4.4 \times 10^{-6}$ at 800 (absolute values, not errors): the instability is almost lost |
| 12 | about $+0.21$ | about $+0.23$ | $+0.204$ |
| 16 | $+2.07 \times 10^{-2}$ | $+2.65 \times 10^{-2}$ | $+2.02 \times 10^{-2}$ |
| 20 | $+1.51 \times 10^{-3}$ | $+1.62 \times 10^{-3}$ | $+1.49 \times 10^{-3}$ |
| 24 | $+1.91 \times 10^{-4}$ | $+2.68 \times 10^{-4}$ | $+1.64 \times 10^{-4}$ |
| 28 | $+8.6 \times 10^{-5}$ | $+9.1 \times 10^{-5}$ | $+6.0 \times 10^{-5}$ |

cbm18a (fixed edge, $n = 10$):

| $W$ | rule A | rule B | PB3D's differences, 3200 points |
|---|---|---|---|
| 10 | $+0.158$ | $+0.186$ | $+0.157$ |
| 15 | $+1.55 \times 10^{-2}$ | $+2.60 \times 10^{-2}$ | $+1.54 \times 10^{-2}$ |
| 20 | between $-2 \times 10^{-5}$ and $+5 \times 10^{-5}$, changing with $N$ | $+4.88 \times 10^{-3}$ | $+4.4 \times 10^{-5}$ |
| 30 | $-9.5 \times 10^{-4}$ | $+4.45 \times 10^{-4}$ | $-9.7 \times 10^{-4}$ |

What the tables say.

- **PB3D's present sliding set is rule A.** In every row its differences approach the number of rule A, not that of rule B; where they differ from it, the differences have not yet converged radially (peeling-ballooning, $W = 24$ and $28$).
- **Rule B lies above the complete value in all nine rows and falls every time $W$ grows**, in both cases.
- **Rule A has no fixed sign and does not always improve when $W$ grows.** In the peeling-ballooning case it is positive and a little closer than rule B for every $W$. In cbm18a it is $5 \times 10^{-5}$ or less at $W = 20$ and $-9.5 \times 10^{-4}$ at $W = 30$: with 30 modes the eigenvalue lies $10^{-3}$ below the complete value, about twenty times farther than with 20 modes. PB3D as it is does the same (−4.71500e-2 at 3200 points with `n_mod_X` = 30, against −4.71021e-2 with 20). A likely reason, not checked: rule A leaves out the energy of a leaving mode beyond its last point, and with $W = 30$ the leaving modes are 15 away from resonance, where that energy is large and stabilizing.
- **The sliding set needs many modes in the peeling-ballooning case.** An error below $10^{-3}$ needs $W = 24$. With $W = 12$ the eigenvalue is 20 % too high and with $W = 8$ the instability is almost lost; the default `n_mod_X` = 10 lies between these two and was not run. With $W = 20$ the sliding set ($1.5 \times 10^{-3}$) is already closer than the deck's fixed set of 26 modes ($2.45 \times 10^{-2}$), because it follows $nq$ from 10.7 at $\psi = 0.5$ to 30.4 at the edge.
- **When the error of the set is large, the eigenvalue changes irregularly with $N$** by several percent, for every scheme (peeling-ballooning, $W = 12$: elements of degree 3, rule A, $+0.32$, $+0.29$, $+0.17$, $+0.21$ at 199, 399, 799, 1600 unknowns). The positions where a mode enters or leaves move with the points. Where the error of the set is below $10^{-3}$ these changes are below $10^{-4}$ (cbm18a, $W = 20$, rule A: between $-2 \times 10^{-5}$ and $+5.3 \times 10^{-5}$ from 799 unknowns on).

## Radial convergence with the sliding set

For the narrowest sets whose error is below $10^{-3}$. Relative errors against the complete values above; "unknowns" is $N$ for the differences.

cbm18a, $W = 20$:

| unknowns | 200 | 400 | 800 | 1600 | 3200 |
|---|---|---|---|---|---|
| differences, order 3 (PB3D) | $-3.94 \times 10^{-2}$ | $-2.55 \times 10^{-3}$ | $-8.4 \times 10^{-5}$ | $+3.8 \times 10^{-5}$ | $+4.4 \times 10^{-5}$ |
| elements, degree 3, rule A | $+6.8 \times 10^{-4}$ | $+9.4 \times 10^{-5}$ | $-2.0 \times 10^{-5}$ | $+2.9 \times 10^{-5}$ | $+5.3 \times 10^{-5}$ |
| elements, degree 2, rule A | $+2.9 \times 10^{-3}$ | $+2.6 \times 10^{-4}$ | $+9.3 \times 10^{-5}$ | $+2.1 \times 10^{-5}$ | |
| cubic Hermite, rule A | $+3.8 \times 10^{-4}$ | $+1.8 \times 10^{-4}$ | $+1.4 \times 10^{-4}$ | $+5.0 \times 10^{-5}$ | |
| cubic B-spline, rule A | $+5.6 \times 10^{-4}$ | $+3.3 \times 10^{-4}$ | $+1.6 \times 10^{-4}$ | $+1.1 \times 10^{-4}$ | |
| elements, degree 3, rule B | $+6.3 \times 10^{-3}$ | $+5.21 \times 10^{-3}$ | $+4.86 \times 10^{-3}$ | $+4.89 \times 10^{-3}$ | $+4.88 \times 10^{-3}$ |
| cubic B-spline, rule B | $+5.8 \times 10^{-3}$ | $+5.4 \times 10^{-3}$ | $+5.1 \times 10^{-3}$ | $+5.0 \times 10^{-3}$ | |

Peeling-ballooning, $W = 24$:

| unknowns | 200 | 400 | 800 | 1600 | 3200 |
|---|---|---|---|---|---|
| differences, order 3 (PB3D) | $-0.154$ | $-2.88 \times 10^{-2}$ | $-3.52 \times 10^{-3}$ | $-1.26 \times 10^{-4}$ | $+1.64 \times 10^{-4}$ |
| elements, degree 3, rule A | $+2.51 \times 10^{-2}$ | $+1.40 \times 10^{-3}$ | $+2.15 \times 10^{-4}$ | $+1.91 \times 10^{-4}$ | $+1.91 \times 10^{-4}$ |
| elements, degree 2, rule A | $+0.101$ | $+8.8 \times 10^{-3}$ | $+7.7 \times 10^{-4}$ | $+2.2 \times 10^{-4}$ | |
| cubic Hermite, rule A | $+7.9 \times 10^{-3}$ | $+3.9 \times 10^{-4}$ | $+1.9 \times 10^{-4}$ | $+1.9 \times 10^{-4}$ | |
| cubic B-spline, rule A | $+3.2 \times 10^{-3}$ | $+2.4 \times 10^{-4}$ | $+2.0 \times 10^{-4}$ | $+2.0 \times 10^{-4}$ | |
| elements, degree 3, rule B | $+2.51 \times 10^{-2}$ | $+1.46 \times 10^{-3}$ | $+2.93 \times 10^{-4}$ | $+2.67 \times 10^{-4}$ | $+2.68 \times 10^{-4}$ |

(The element columns hold 199, 399 or 400, 799 or 800, 1600, 3199 unknowns.)

What the tables say.

- **The radial gain of the elements stays with the sliding set.** Peeling-ballooning, $W = 24$: elements of degree 3 are $2.4 \times 10^{-5}$ from their own limit at 799 unknowns; the differences are $3.7 \times 10^{-3}$ from theirs at 800 points and $3 \times 10^{-4}$ at 1600.
- **Lagrange elements fit the sliding set; Hermite functions and B-splines do so less well.** On cbm18a with rule A, cubic B-splines fall towards their limit only about as $1/N$ ($5.6$, $3.3$, $1.6$, $1.1 \times 10^{-4}$) and cubic Hermite functions irregularly ($3.8$, $1.8$, $1.4$, $0.5 \times 10^{-4}$), while elements of degree 3 are within $10^{-4}$ from 400 unknowns on. With a fixed set the order was the reverse (issue #51). A likely reason, not checked: a B-spline or Hermite function that carries an entering mode starts with zero value and zero derivative, which a Lagrange element does not impose. On the peeling-ballooning case the difference does not show.
- With rule B every scheme converges smoothly, to the larger error of rule B.

## Extrapolation by `interp_V` (risk 2 of the plan)

The same series with the coefficients of the sliding-set runs, which `interp_V` extrapolates where a Gauss point lies outside the perturbation points of a pair, against the series with the wide fixed set (no extrapolation). Largest relative difference of the eigenvalue over the schemes of the tables and $N$ = 200 to 1600:

| case, $W$ | at about 200 unknowns | from 400 unknowns on |
|---|---|---|
| cbm18a, 20 | $2.7 \times 10^{-7}$ | $1.2 \times 10^{-6}$ |
| cbm18a, 30 | $6 \times 10^{-8}$ | $1.2 \times 10^{-8}$ |
| peeling-ballooning, 20 | $7.8 \times 10^{-5}$ | $2.2 \times 10^{-6}$ |
| peeling-ballooning, 24 | $1.2 \times 10^{-5}$ | $1.6 \times 10^{-6}$ |
| peeling-ballooning, 12 and 16; cbm18a, 10 | up to $3 \times 10^{-4}$ | up to $2.4 \times 10^{-3}$ at single $N$, inside the irregular changes of these sets |

Where the set is wide enough to be used, the extrapolation changes the eigenvalue by less than $10^{-5}$ from 400 unknowns on. The fallback of the plan (limit the set of an element to the modes that the perturbation grid has over the whole element) is not needed.

## What this changes in the plan

Worked into `../plan_radial_finite_elements.md` on 20 November 2026 (sections 3.3, 4.5 to 4.7, 5, 7, 8 and 10).

The method of section 4 of the plan is not changed by these results: one set per element taken at its middle, rule A with `insert_block_mat` as it is, the coefficient matrices interpolated by `interp_V` to the Gauss points. What changes is what the plan says about the sliding set and how it is tested.

- **Section 3.3 and decision 4 (rule A or rule B).** The model numbers are replaced by the tables above. Rule A is what PB3D does now, so it keeps the present eigenvalues of `X_style` 2 runs, but its error has no fixed sign and can grow when `n_mod_X` grows (cbm18a: $5 \times 10^{-5}$ with 20 modes, $-9.5 \times 10^{-4}$ with 30). Rule B always lies above the complete value and fell every time `n_mod_X` grew, but at the same `n_mod_X` its error can be a hundred times larger (cbm18a, 20 modes: $4.9 \times 10^{-3}$ against $5 \times 10^{-5}$). The choice stays Toon's; it is now a choice between keeping PB3D's present numbers and having an error with a known sign that falls with `n_mod_X`.
- **Risk 2 (extrapolation at Gauss points).** Measured: below $10^{-5}$ from 400 unknowns on wherever the set is usable. The fallback is not needed.
- **Risk 3 (irregular changes with $N$).** Measured: below $10^{-4}$ where the error of the set is below $10^{-3}$, several percent where it is larger.
- **Gate 3b.** `sliding_set.py` builds the matrices for any table and either rule from `PB3D_out.h5` alone, including the extrapolation of `interp_V`. The debug switch that writes the coefficient matrices at the Gauss points is not needed for the comparison.
- **Gate 3d.** Its numbers for the sliding set come from the cut regression deck. Numbers from these cases: cbm18a, $W = 20$, rule A, degree 3: within $10^{-4}$ of $-4.710420 \times 10^{-2}$ from 400 unknowns on; peeling-ballooning, $W = 24$, rule A, degree 3: within $5 \times 10^{-5}$ of $-1.842206 \times 10^{-2}$ (its own limit) from 799 unknowns on.
- **Decision 1 (functions).** With the sliding set and the two rules as defined here, cubic B-splines and Hermite functions converge more slowly than elements of degree 3 on cbm18a; with a fixed set they were faster. Choosing them would need a rule for the sliding set of their own.
- **The reference of the peeling-ballooning case.** The fixed set of the deck ($m = 12 \ldots 37$) is $2.45 \times 10^{-2}$ above the complete value. The radial series of issue #51 stay valid as series for that set.

Not part of the trial, but found on the way: with PB3D as it is, a larger `n_mod_X` does not always give a better eigenvalue (cbm18a, 20 against 30 modes), so a user cannot check the sliding set by raising `n_mod_X` until the eigenvalue stops changing in one direction.

## Reproducing

```bash
# PB3D runs (tol_SLEPC 1e-10; n_mod_X = W in place of min_sec_X/max_sec_X for the sliding set)
PB3D <input> <equilibrium> -eps_view_mat0 binary:A.dat -eps_view_mat1 binary:B.dat
python -B sliding_set.py compare <run folder>
python -B sliding_set.py scan <PB3D_out.h5 of the wide fixed set> <W> <reference> d3,e3A,e3B,e2A,hA,bA,bB 200 400 800 1600
python -B sliding_set.py scan <PB3D_out.h5 of the sliding run> <W> <reference> d3,e3A,e3B 200 400 800 1600 --vac=<PB3D_out.h5 of the wide fixed set>
```

Wide fixed sets: `min_sec_X = 5`, `max_sec_X = 48` (peeling-ballooning) and `75` (cbm18a); everything else as in `inputs/`. The output files of these runs are 0.4 and 1.4 GB. A Python series of one set takes 1 to 5 minutes on the Mac.
