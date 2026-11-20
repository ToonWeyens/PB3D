# Plan: finite elements in the radial direction

Branch `dev/radial-finite-elements`, started from `dev/vacuum-fixes` at f9713d7 (version 2.50). Written on 20 November 2026. GitHub issue: [#51](https://github.com/ToonWeyens/PB3D/issues/51).

**State.** This is a plan with the measurements it rests on. Nothing in `Modules/` is changed. The branch holds this file and the scripts of `Documentation/radial_finite_elements/`. It is a trial: nothing here is decided or required, and whether the finite elements replace the finite differences in `dev/vacuum-fixes` is Toon's decision after the gates of section 7.

Statements are marked where it matters: **read** (read in the code at f9713d7), **measured** (a run of PB3D or of a script, reproducible with section 12), **standard** (a known result of the finite-element method, not derived here), **judgement** (mine). Two independent reviews, one of the numbers and one of the plan against the code, were worked into this version.

Symbols used throughout: $N$ is `n_r_sol`, the number of radial points; $\Delta$ is the distance between two neighbouring radial points; $p$ is `norm_disc_prec_sol` of the present code, the order of its finite difference; $k$ is the degree of the polynomials of the proposed scheme; $n_\text{el}$ is the number of elements; $n$ is the number of unknowns per Fourier mode in the model script. Eigenvalues are $\lambda = \omega^2/\omega_A^2$. "Relative error" always means the difference to a reference eigenvalue divided by the size of that reference.

## 1. Summary

**What is proposed.** The radial derivative and the radial integral of the energy are no longer a left finite difference and a weighted sum over the points. The displacement of every Fourier mode becomes a continuous function that is a polynomial of degree $k$ on each element, where an element is a group of $k + 1$ neighbouring radial points and neighbouring elements share their end point. The energy of that function is integrated with $k + 1$ Gauss–Legendre points inside every element. The unknowns stay what they are now: the values of $X$ at $N$ equally spaced radial points.

**The present differences have an error proportional to $1/N$; the elements are built so that they do not have it.** With the present differences and $p = 3$, the eigenvalue of the cbm18a fixed-boundary deck is off by $0.0554/N$ (**measured** with PB3D, section 3.1), which is $6.3 \times 10^{-3}$ of the eigenvalue at 200 points and $1.6 \times 10^{-3}$ at 800 points. The cause is measured too: the differences displace the two ends of the plasma. 94 % of the error comes from the inner end and 6 % from the fixed edge:

- PB3D holds the first $p$ points at zero, so the plasma it computes starts $(p - \tfrac{1}{2})$ steps outside `min_r_sol`;
- at a fixed edge all integration weights are 1, which places the edge 0.28 steps inside `max_r_sol` for $p = 3$.

These two displacements, multiplied by the measured dependence of the eigenvalue on the position of each end, give $0.0552/N$.

The elements do not displace the ends: one point is held at each end and the integral runs from `min_r_sol` to `max_r_sol`. In the model the error of the elements has no term proportional to $1/N$: it falls as $1/N^{2k}$ at a free edge, at a fixed edge and when the eigenfunction reaches the inner end with a slope (section 3.2). A PB3D run with elements has not been made; gate 3d is that test.

**What depends on the case: the gain from the higher power.** Away from the ends the error of the elements falls as $1/N^{2k}$, once the eigenfunction is resolved. Call the length of change of an eigenfunction $\ell = \sqrt{\int |X|^2 \mathrm{d}\psi / \int |X'|^2 \mathrm{d}\psi}$. In the model (section 3.2), relative error of the lowest eigenvalue:

| radial points per length $\ell$ | differences $p = 3$ | elements $k = 2$ | elements $k = 3$ |
|---|---|---|---|
| 30 (six Fourier modes, 193 unknowns) | $4.4 \times 10^{-6}$ | $2.7 \times 10^{-7}$ | $1.5 \times 10^{-9}$ |
| 13 (twenty Fourier modes, 193 unknowns) | $3.5 \times 10^{-5}$ | $2.3 \times 10^{-5}$ | $2.5 \times 10^{-6}$ |
| 3.2 (twenty Fourier modes, 49 unknowns) | $2.3 \times 10^{-3}$ | $3.7 \times 10^{-3}$ | $4.9 \times 10^{-3}$ |

With 13 points per length, $k = 2$ is no better than the present differences and $k = 3$ is 14 times better. With 3.2 points per length the elements are 1.6 to 2.1 times worse. The eigenfunction of the cbm18a deck has $\ell = 0.075$ of its radial range (**measured**, section 3.1), which is 15 points at 200 points and 60 at 800, and its derivatives change over lengths that are 1.2 to 1.4 times those of the twenty-mode problem. Its error away from the ends is therefore expected near the middle row of the table at 200 points, which is more than a hundred times smaller than the error of the ends that the elements remove. This expectation is tested by gate 3d before anything is adopted.

**What it costs.**

- The matrices keep their shape: one block of `n_mod_X` unknowns per radial point. For $k = 3$ a row has at most 7 blocks, as now, and 5 on average (section 3.5).
- The coefficient matrices are needed at the Gauss points, which are $\tfrac{4}{3}$ times as many as the radial points for $k = 3$.
- $N$ has to be $k \, n_\text{el} + 1$; a deck with 200 points gets 202 for $k = 3$.
- Every regression value of an eigenvalue moves. POST has to build the same tables of Fourier modes as PB3D, so `min_PB3D_version` has to be raised.
- Chapter 7 of the theory document describes the differences and would have to follow.

**What the elements do not change.**

- The sliding set of Fourier modes (`X_style` 2) keeps an error of its own, which depends on `n_mod_X`. On the cbm18a deck with 20 modes it is $4 \times 10^{-4}$ of the eigenvalue (section 3.1), which is the size of the $1/N$ error of the differences at 3200 points.
- The coefficient matrices are known on the equilibrium grid only. Their interpolation moves the eigenvalue by an amount that no radial scheme removes (section 3.4).

**A cheaper alternative exists.** The two end errors could also be removed inside the present scheme, by a rule for the first points and a rule for a fixed edge of the kind that the sum rule is for a free edge. Neither rule exists yet; section 7.16.5 of the theory document lists both as not settled. Toon asked for the elements. I name the alternative because the certain gain is the same and the work is smaller.

**Decisions that are Toon's** are listed in section 10. The plan assumes my recommendation for each and says where another choice changes the work.

## 2. PB3D as it is

All of this is **read** in the code at f9713d7.

| What | Where | How |
|---|---|---|
| Radial points | `grid_ops.f90:434` | $N$ equally spaced points from `min_r_sol` to `max_r_sol` |
| Coefficient matrices at these points | `driver_sol.f90:175-203`, `sol_utilities.f90:805` | `interp_V` interpolates `PV_0..2`, `KV_0..2` from the perturbation grid, one cubic spline per pair of Fourier modes over the points where both modes exist; outside these points it extrapolates |
| Derivative | `SLEPC_ops.f90:700-704` | left difference of order $p$, coefficients from `calc_coeff_fin_diff` |
| Integral | `SLEPC_ops.f90:685-713` | `step_size` times a weight per point. Free edge: the weights of `calc_end_weights` (the sum rule). Fixed edge: all weights 1 |
| Matrices | `fill_mat`, `SLEPC_ops.f90:605` | per point one block with `V_0`, $p + 1$ blocks with `V_1` and their conjugates, $(p + 1)^2$ blocks with `V_2` |
| Number of blocks per row | `set_nonzeros`, `SLEPC_ops.f90:884-895` | already counted on both sides of the diagonal, $p$ blocks each |
| Inner end | `set_BC`, `SLEPC_ops.f90:989` | the first $p$ points are held at zero with the artificial eigenvalue `EV_BC` |
| Edge | `set_BC_1`, `set_BC_2` | fixed: last point held at zero; free: `vac%res` added to the last diagonal block |
| Sliding set of modes | `X_ops.f90:1059`, `SLEPC_utilities.f90:212-237` | every point has `n_mod_X` slots. For a block that belongs to one point, `insert_block_mat` keeps every entry in the slot of its Fourier mode at the row point and at the column point, and drops the entry if the mode has no slot there. No entry is moved to another slot |
| Tables of modes in POST | `driver_POST.f90:202`, `sol_utilities.f90:344` | POST does not read the table from the output file; it calls `setup_modes` on the solution grid again and labels the slots of the eigenvector with the result |
| Derivative in POST | `sol_utilities.f90:663` | a spline of order `norm_disc_prec_sol` through the values at the points; not the difference of the matrices |
| Stored | `input_ops.f90:1701`, `sol_ops.f90:1574` | `norm_disc_prec_sol` in `misc_sol(3)`; the eigenvectors as values per slot and point |

Three consequences for this plan.

- A scheme whose unknowns are the values at the same points leaves the solution vector and the output file as they are. POST depends on the scheme in two places: it builds the table of modes itself, so PB3D and POST must use the same rule, and it takes spline orders and ghost widths from `norm_disc_prec_sol`.
- `insert_block_mat` already places a block that belongs to one point at rows and columns of other points. An element needs the same operation.
- The coefficient matrices exist only at the radial points. Gauss points lie between them, so the interpolation has to give the coefficient matrices at another set of positions.

## 3. Measurements

### 3.1 The $1/N$ of cbm18a comes from the two ends (PB3D, **measured**)

Deck `tests/regression/input_cbm18a_fixed_bnd` (fixed edge, toroidal mode number 10, range 0.45 to 0.75 in the normalized poloidal flux), `tol_SLEPC` $10^{-12}$, one process, the build of this branch without changes. In this range $n q$ goes from 14.9 to 37.1, so 23 rational surfaces lie inside it, and the deck has 158 equilibrium points in it. Three sets of Fourier modes: the sliding set of 20 modes of the deck, which uses $m = 5 \ldots 47$ over the range; a fixed set $m = 5 \ldots 40$; the complete fixed set $m = 5 \ldots 47$. The eigenvalue is given in units of $10^{-2}$.

| $N$ | sliding, $p = 3$ | $5 \ldots 40$, $p = 3$ | $5 \ldots 47$, $p = 3$ | sliding, $p = 2$ | $5 \ldots 40$, $p = 2$ |
|---|---|---|---|---|---|
| 200 | −4.357290 | −4.350099 | | −4.588865 | −4.582066 |
| 400 | −4.365940 | −4.358696 | | −4.428302 | −4.421165 |
| 800 | −4.372961 | −4.365580 | −4.374661 | −4.389857 | −4.382503 |
| 1600 | −4.376516 | −4.369078 | −4.378223 | −4.381351 | −4.373920 |
| 3200 | −4.378267 | −4.370808 | −4.379985 | −4.379781 | −4.372324 |

A fit $\lambda_\infty + a/N$ to the last two values ($p = 3$), or $\lambda_\infty + a/N + b/N^2$ to the last three ($p = 2$), gives (script `cbm18a_ends.py`; $a$ in absolute units, $\lambda_\infty$ in $10^{-2}$):

| set | $p$ | $a$ | $\lambda_\infty$ |
|---|---|---|---|
| sliding set of 20 | 3 | 0.0560 | −4.38002 |
| fixed set $5 \ldots 40$ | 3 | 0.0554 | −4.37254 |
| complete fixed set $5 \ldots 47$ | 3 | 0.0564 | −4.38175 |
| sliding set of 20 | 2 | 0.0356 | −4.38000 |
| fixed set $5 \ldots 40$ | 2 | 0.0352 | −4.37252 |

The two orders agree on $\lambda_\infty$ to $2 \times 10^{-5}$ in these units. The sliding set has the same $a$ as the fixed set to 1 %, so the sliding set is not the cause of the $1/N$.

**Where the differences place the two ends** (model, **measured**, section `ends` of `fe_model.py`). On a one-mode problem whose eigenvalue depends on the position of an end, the error of the differences equals the error of a problem whose end is displaced:

- inner end: the plasma starts 1.5000 steps ($p = 2$) and 2.5000 steps ($p = 3$) further out, which is $p - \tfrac{1}{2}$, the value of section 7.12 of the theory document;
- fixed edge with all weights 1: the plasma ends 0.0000 ($p = 1$), 0.3750 ($p = 2$), 0.2805 ($p = 3$) steps further in.

**How much the eigenvalue depends on each end** (PB3D, fixed set $5 \ldots 40$, $p = 3$, 800 points). With `min_r_sol` = 0.452, 0.454, 0.458 the eigenvalue is −4.350905, −4.335130, −4.299949. With `max_r_sol` = 0.748, 0.746, 0.742 it is −4.357958, −4.349602, −4.330580. Both are compared with −4.365580 for the range 0.45 to 0.75. A cubic through each set of four runs gives the slope at the nominal end: 6.86 per unit of normalized flux for the inner end moved outwards and 3.64 for the edge moved inwards (in $10^{-2}$).

**Prediction against measurement.** The coefficient $a$ is the slope times the number of steps times the range 0.3:

| $p$ | inner end | fixed edge | their sum, with the same fit as the measurement | measured |
|---|---|---|---|---|
| 3 | 0.0515 (94 %) | 0.0031 (6 %) | 0.0552 | 0.0554 |
| 2 | 0.0309 (88 %) | 0.0041 (12 %) | 0.0350 | 0.0352 |

The column "their sum" applies the fit of the measurement to the predicted eigenvalues, which include the curvature of the two cubics; the plain sums are 0.0546 and 0.0350.

**Why the ends matter in this deck: the shape of the eigenfunction** (PB3D, complete fixed set, 800 points, `tol_SLEPC` $10^{-9}$, script `eigenfunction_profile.py`). Take $|X|$ as the root of the sum of $|X_m|^2$ over the Fourier modes, scaled to a largest value of 1. The largest value lies at 45 % of the range from the inner end, at the normalized flux 0.586. The eigenfunction does not fall to zero towards the ends by itself: from the held points at the inner end it rises with a slope of 1.55 per range, against a largest slope anywhere of 14.7, and reaches 0.026, 0.067 and 0.123 at 0.02, 0.05 and 0.10 of the range from the inner end. From the edge the same three values are 0.014, 0.034 and 0.068. The boundary condition at each end holds at zero a function that would not be zero there, and the eigenvalue then depends on where the end is.

**How fast the eigenfunction changes** (same run). Three lengths, as fractions of the radial range, compared with the model problems of section 3.2 (section `lengths` of `fe_model.py`):

| | $\sqrt{\int \lvert X \rvert^2 / \int \lvert X' \rvert^2}$ | $\sqrt{\int \lvert X' \rvert^2 / \int \lvert X'' \rvert^2}$ | $\sqrt{\int \lvert X'' \rvert^2 / \int \lvert X''' \rvert^2}$ |
|---|---|---|---|
| cbm18a, complete fixed set | 0.0754 | 0.0327 | 0.0152 |
| model, twenty Fourier modes | 0.0675 | 0.0240 | 0.0129 |
| model, six Fourier modes | 0.1549 | 0.0705 | 0.0506 |

The first column is the length of change $\ell$ of section 1. The eigenfunction of cbm18a changes a little more slowly than that of the twenty-mode problem in all three columns, by factors 1.1, 1.4 and 1.2. So the twenty-mode problem at $n$ unknowns is a fair stand-in for cbm18a at about the same number of radial points, as far as the error away from the ends goes.

**The size of the error of the sliding set in this deck.** The sliding set of 20 modes has the limit $-4.38002 \times 10^{-2}$ and the complete fixed set $-4.38175 \times 10^{-2}$: the sliding set lies $4 \times 10^{-4}$ (relative) above the complete set. The fixed set $5 \ldots 40$ lacks $m = 41 \ldots 47$, which is why its limit is the highest of the three.

### 3.2 Orders in the model (**measured**)

`fe_model.py` writes every scheme as one table (section 4.3) and builds the matrices from the table in one loop. The schemes therefore differ in their table only. The differences are those of PB3D: the weights of the sum rule at a free edge, all weights 1 at a fixed edge, the first $p$ points held. Three problems, on $0 \le \psi \le 1$:

- **one Fourier mode**: the test problem of section 7.10.6 of the theory document, $-X'' + V X = \lambda X$ with $V = 4000 (1 - \psi)^6$, whose eigenfunction is small at the inner end, and the same with $V = 0$, whose eigenfunction reaches the inner end with a slope;
- **six Fourier modes**: coefficient matrices that are chosen, not computed from an equilibrium. $\mathsf{P}^2$ and $\mathsf{P}^0$ contain $(n q - m)^2$, $\mathsf{K}^2$ is not zero, $\mathsf{P}^1$ and $\mathsf{K}^1$ are complex and not Hermitian, and neighbouring modes are coupled. Five rational surfaces lie inside the plasma. The eigenfunction changes over $\ell = 0.155$ (free edge) and $0.129$ (fixed edge); it has no thin layers;
- **twenty Fourier modes**: the same construction with 16 rational surfaces inside and $\mathsf{K}^2$ ten times smaller, which gives $\ell = 0.0675$.

The numbers are the relative error of the lowest eigenvalue. The reference comes from polynomials of degree 8 and is known to between $10^{-11}$ and $2 \times 10^{-10}$; values at that size are left out.

Six Fourier modes, free edge:

| scheme | $n = 97$ | 193 | 385 | 769 | falls as |
|---|---|---|---|---|---|
| differences $p = 2$ (PB3D) | 2.3e-3 | 5.9e-4 | 1.5e-4 | 3.7e-5 | $1/n^{2}$ |
| differences $p = 3$ (PB3D) | 7.5e-5 | 4.4e-6 | 2.5e-7 | 2.0e-8 | $1/n^{3.7}$ |
| elements $k = 1$ | 4.1e-4 | 1.0e-4 | 2.5e-5 | 6.4e-6 | $1/n^{2}$ |
| elements $k = 2$ | 4.4e-6 | 2.7e-7 | 1.7e-8 | 1.1e-9 | $1/n^{4}$ |
| elements $k = 3$ | 9.4e-8 | 1.5e-9 | 2.5e-11 | | $1/n^{6}$ |
| elements $k = 2$, Gauss–Lobatto integration (issue #35) | 2.5e-6 | 1.6e-7 | 9.8e-9 | 6.1e-10 | $1/n^{4}$ |
| elements $k = 3$, Gauss–Lobatto integration | 7.4e-8 | 1.2e-9 | 2.0e-11 | | $1/n^{6}$ |
| cubic Hermite | 2.8e-8 | 4.6e-10 | | | $1/n^{6}$ |
| cubic B-spline | 5.3e-10 | | | | |

Six Fourier modes, fixed edge:

| scheme | $n = 97$ | 193 | 385 | 769 | falls as |
|---|---|---|---|---|---|
| differences $p = 2$ (PB3D) | 5.3e-3 | 1.9e-3 | 7.6e-4 | 3.3e-4 | $1/n^{1.2}$ |
| differences $p = 3$ (PB3D) | 1.7e-3 | 8.5e-4 | 4.3e-4 | 2.1e-4 | $1/n$ |
| elements $k = 1$ | 4.7e-4 | 1.2e-4 | 3.0e-5 | 7.4e-6 | $1/n^{2}$ |
| elements $k = 2$ | 5.3e-6 | 3.3e-7 | 2.1e-8 | 1.3e-9 | $1/n^{4}$ |
| elements $k = 3$ | 1.2e-7 | 1.9e-9 | 3.3e-11 | | $1/n^{6}$ |

One Fourier mode with $V = 0$, free edge (the eigenfunction reaches the inner end with a slope):

| scheme | $n = 49$ | 97 | 193 | falls as |
|---|---|---|---|---|
| differences $p = 2$ (PB3D) | 5.6e-2 | 2.7e-2 | 1.3e-2 | $1/n$ |
| differences $p = 3$ (PB3D) | 9.4e-2 | 4.5e-2 | 2.2e-2 | $1/n$ |
| elements $k = 1$ | 1.5e-4 | 3.7e-5 | 9.3e-6 | $1/n^{2}$ |
| elements $k = 2$ | 4.8e-8 | 3.0e-9 | 2.1e-10 | $1/n^{4}$ |

Twenty Fourier modes, free edge; the line "points per length" gives the number of radial points per length $\ell$:

| scheme | $n = 49$ | 97 | 193 | 385 | 769 |
|---|---|---|---|---|---|
| points per length | 3.2 | 6.5 | 13 | 26 | 52 |
| differences $p = 2$ (PB3D) | +1.5e-2 | +3.9e-3 | +9.4e-4 | +2.3e-4 | +5.8e-5 |
| differences $p = 3$ (PB3D) | −2.3e-3 | −5.7e-4 | −3.5e-5 | −7.1e-7 | +1.7e-7 |
| elements $k = 1$ | +1.2e-2 | +3.1e-3 | +7.7e-4 | +1.9e-4 | +4.8e-5 |
| elements $k = 2$ | +3.7e-3 | +3.1e-4 | +2.3e-5 | +1.4e-6 | +9.1e-8 |
| elements $k = 3$ | +4.9e-3 | +4.9e-5 | +2.5e-6 | +4.0e-8 | +6.5e-10 |

What the tables say.

- The elements fall as $1/n^{2k}$ in all four problems, without any rule for the first or the last points.
- The differences of PB3D fall as $1/n$ whenever the eigenfunction reaches a fixed edge or the inner end with a slope. They reach $1/n^3$ to $1/n^{3.7}$ only at a free edge with an eigenfunction that is small at the inner end.
- With Gauss–Legendre integration every error of the elements is positive: the discrete eigenvalue lies above the true one. The differences have errors of both signs. With Gauss–Lobatto integration the sign is not fixed either (it is negative in the one-mode problem); its errors are 0.3 to 0.8 times those of Gauss–Legendre for the lowest eigenvalue and 1.6 times for the lowest ten (section 3.5).
- At equal $n$, the error of $k = 3$ is 2900 times smaller than that of the differences $p = 3$ with six modes and 193 unknowns, 14 times smaller with twenty modes and 193 unknowns, and 2.1 times larger with twenty modes and 49 unknowns. Between 3.2 and 6.5 points per length the elements become the better scheme.

### 3.3 The sliding set of Fourier modes (model, **measured**)

Seventeen Fourier modes, $m = 7 \ldots 23$, of which every position uses the 8, 10 or 12 closest to $n q$. The numbers are the relative distance of the lowest eigenvalue to the one with all 17 modes everywhere. They do not go to zero with $n$: what is left is the error of using fewer modes. Two rules for the elements:

- **Rule A.** One set of modes per element, taken at the middle of the element. A point shared by two elements carries the unknowns of the element on its inner side. A mode is then zero at the first point of its first element and has a free value at the last point of its last element, with no energy beyond it. This is what the present differences do in the limit of many points.
- **Rule B.** One set per element, taken at its middle. A shared point carries only the modes that both elements have. Every mode is then a continuous function that is zero at both ends of its range, and the matrices hold the complete energy of that function. The discretized problem is the problem with all modes, restricted to fewer functions, so its lowest eigenvalue lies above the one with all modes, and it can only fall when `n_mod_X` grows.

$q$ rising from 1.21 to 1.74, largest drive at $\psi = 0.8$:

| set | scheme | $n = 97$ | 193 | 385 | 769 |
|---|---|---|---|---|---|
| 8 | differences $p = 3$, set per point (PB3D now) | +3.6e-5 | −7.9e-5 | −6.7e-5 | −7.3e-5 |
| 8 | elements $k = 3$, rule A | −8.0e-6 | −2.0e-5 | −6.4e-5 | −8.9e-5 |
| 8 | elements $k = 3$, rule A with the set at the last point | −4.2e-4 | −2.0e-4 | −1.1e-4 | −8.8e-5 |
| 8 | elements $k = 3$, rule B | +1.4e-3 | +1.4e-3 | +1.5e-3 | +1.5e-3 |
| 10 | differences $p = 3$, set per point (PB3D now) | +8.4e-5 | +5.7e-6 | +8.6e-7 | +3.2e-7 |
| 10 | elements $k = 3$, rule A | +1.3e-6 | +1.1e-6 | +6.0e-7 | −1.6e-7 |
| 10 | elements $k = 3$, rule A with the set at the last point | −7.0e-6 | −2.9e-6 | −8.5e-7 | −1.5e-7 |
| 10 | elements $k = 3$, rule B | +3.3e-5 | +3.3e-5 | +3.6e-5 | +3.5e-5 |
| 12 | differences $p = 3$, set per point (PB3D now) | +8.3e-5 | +5.5e-6 | +4.1e-7 | +6.1e-8 |
| 12 | elements $k = 3$, rule A | +1.2e-7 | +2.7e-8 | +2.7e-8 | +1.9e-8 |
| 12 | elements $k = 3$, rule A with the set at the last point | +6.1e-8 | −2.1e-9 | +1.1e-8 | +1.9e-8 |
| 12 | elements $k = 3$, rule B | +3.3e-7 | +2.4e-7 | +2.7e-7 | +2.6e-7 |

Three variants of the problem at $n = 769$; the last column is the size of the number for rule A divided by the size of the number for rule B:

| variant | set | PB3D now | rule A | rule B | A / B |
|---|---|---|---|---|---|
| $q$ rising, drive at 0.8 | 8 | −7.3e-5 | −8.9e-5 | +1.5e-3 | 0.06 |
| | 10 | +3.2e-7 | −1.6e-7 | +3.5e-5 | 0.005 |
| | 12 | +6.1e-8 | +1.9e-8 | +2.6e-7 | 0.07 |
| $q$ falling from 1.74 to 1.21, drive at 0.8 | 8 | −9.2e-4 | −9.3e-4 | +4.7e-4 | 2.0 |
| | 10 | −1.2e-5 | −1.2e-5 | +8.2e-6 | 1.5 |
| | 12 | +1.3e-7 | −6.6e-8 | +5.9e-8 | 1.1 |
| $q$ rising, drive at 0.4 | 8 | −2.5e-4 | −2.6e-4 | +3.6e-4 | 0.73 |
| | 10 | −2.7e-5 | −3.5e-5 | +5.7e-5 | 0.62 |
| | 12 | +4.3e-6 | −3.9e-6 | +9.6e-6 | 0.41 |

What the tables say.

- Neither rule is closer to the all-modes eigenvalue in general. In the variant with $q$ rising and the drive at 0.8, rule A is 14 to 220 times closer; with $q$ falling it is 1.1 to 2.0 times farther; with the drive at 0.4 it is 1.4 to 2.4 times closer.
- Rule A has no fixed sign: it is negative in 8 of the 9 rows, which means an eigenvalue below the one with all modes. Rule B is positive in all 9.
- Where the differences of PB3D have stopped changing with $n$ (the set of 8 in all variants, the set of 10 with $q$ falling), rule A agrees with them to between 1 % and 21 % of the number. This fits the argument that both have the same limit; it does not prove it.
- The set of an element has to be taken at its middle. Taken at its last point, the same rule has an extra error proportional to $1/n$ (set of 8: −4.2e-4, −2.0e-4, −1.1e-4 against −8.0e-6, −2.0e-5, −6.4e-5).
- With either rule the eigenvalue changes irregularly with $n$ by an amount of the size of this error, because the positions where a mode enters or leaves move with the points. The differences and the elements have this in common.

### 3.4 Coefficients known on the equilibrium grid only (model, **measured**)

The six-mode problem with its coefficient matrices tabulated on $n_\text{eq}$ equally spaced points and interpolated by a cubic spline, as `interp_V` does.

**The interpolation by itself** moves the eigenvalue by $1.3 \times 10^{-5}$ for $n_\text{eq} = 41$ and by $4.4 \times 10^{-8}$ for $n_\text{eq} = 161$ (relative). No radial scheme gets closer to the true eigenvalue than this.

**The error of each scheme against the problem with the interpolated coefficients:**

| $n_\text{eq}$ | scheme | $n = 97$ | 193 | 385 |
|---|---|---|---|---|
| 41 | differences $p = 3$ (PB3D) | 7.5e-5 | 4.4e-6 | 2.5e-7 |
| 41 | elements $k = 2$ | 4.3e-6 | 2.8e-7 | 1.7e-8 |
| 41 | elements $k = 3$ | 1.2e-7 | 2.7e-9 | 1.2e-10 |
| 41 | elements $k = 3$, 6 Gauss points | 4.7e-8 | 2.5e-10 | |
| 41 | cubic B-spline | 3.8e-9 | | |
| 161 | differences $p = 3$ (PB3D) | 7.5e-5 | 4.4e-6 | 2.5e-7 |
| 161 | elements $k = 2$ | 4.4e-6 | 2.7e-7 | 1.7e-8 |
| 161 | elements $k = 3$ | 1.2e-7 | 3.0e-9 | 1.1e-10 |
| 161 | elements $k = 3$, 6 Gauss points | 6.7e-8 | | |
| 161 | cubic B-spline | 5.6e-10 | | |

- $k = 2$ keeps its $1/n^4$. $k = 3$ with 4 Gauss points falls as $1/n^5$ here, not $1/n^6$: the third derivative of a cubic spline jumps at its points, and the Gauss points do not know where.
- Two things restore the power in the model. With the element boundaries on the points of the spline ($n_\text{eq} = 41$; $n = 121, 241, 481$) the errors of $k = 3$ are 2.5e-8, 3.9e-10, 5.6e-12, which is $1/n^{6}$. With 6 Gauss points per element the error at 193 unknowns is 11 to 16 times smaller than with 4.
- In this model with $n_\text{eq} = 161$, which is close to the 158 equilibrium points of cbm18a, the error of the scheme falls below the $4.4 \times 10^{-8}$ of the interpolation near 120 unknowns for $k = 3$, near 300 for $k = 2$ and near 600 for the differences $p = 3$.

### 3.5 Gauss points, the lowest ten eigenvalues, size of the matrices (model, **measured**)

**Number of Gauss points per element** (six modes, free edge). $k$ points give the same power as $k + 1$ and errors that are 1.5, 1.3 and 1.2 times larger for $k = 1, 2, 3$; $k + 2$ points change nothing in the digits shown. The factor depends on the problem: in the lowest ten eigenvalues below, 3 points give a smaller largest error than 4. The reason to take $k + 1$ is in section 4.7: with $k$ points the term with $\mathsf{K}^0$ is not integrated exactly even for constant coefficients, and the property that the eigenvalue comes from above is no longer guaranteed.

**The lowest ten eigenvalues** (six modes, 193 unknowns; all ten are negative):

| scheme | free edge: largest relative error | number below the reference | fixed edge: largest relative error | number below the reference |
|---|---|---|---|---|
| differences $p = 3$ (PB3D) | 7.7e-4 | 2 | 5.5e-3 | 0 |
| elements $k = 1$ | 3.5e-3 | 0 | 7.0e-3 | 0 |
| elements $k = 2$ | 7.6e-5 | 0 | 2.0e-4 | 0 |
| elements $k = 3$ | 3.3e-6 | 0 | 1.1e-5 | 0 |
| elements $k = 3$, 3 Gauss points | 1.7e-6 | 0 | 5.8e-6 | 0 |
| elements $k = 3$, Gauss–Lobatto integration | 5.4e-6 | 0 | 1.8e-5 | 0 |

No scheme has an eigenvalue among its lowest ten that the reference does not have.

**Two tables that must give the same matrices** (section `same`). Elements with $k = 1$ and one energy point at the outer end of each element give the eigenvalue of the differences with $p = 1$ and all weights 1, with a difference of exactly 0. Equally spaced points inside an element and points at the Gauss–Lobatto positions give the same eigenvalue to $1.6 \times 10^{-13}$: the functions are the same, only their description differs.

**Size of the matrices** for 193 unknowns per Fourier mode:

| scheme | energy points | blocks per row, mean | at most |
|---|---|---|---|
| differences $p = 3$ | 190 | 6.9 | 7 |
| elements $k = 1$ | 384 | 3.0 | 3 |
| elements $k = 2$ | 288 | 4.0 | 5 |
| elements $k = 3$ | 256 | 5.0 | 7 |
| cubic Hermite | 380 | 6.0 | 6 |
| cubic B-spline | 760 | 6.9 | 7 |

**The derivative at the radial points, as POST needs it** (section `post`; one mode, $V = 0$, exact eigenfunction known). Largest error of $X'$ at the points, with $X$ scaled to the exact eigenfunction:

| scheme | unknowns | $X'$ of the element polynomials (mean of the two sides at a shared point) | $X'$ of a cubic spline through the point values (POST now) |
|---|---|---|---|
| elements $k = 2$ | 65 | 6.8e-4 | 5.8e-7 |
| elements $k = 2$ | 129 | 1.7e-4 | 7.3e-8 |
| elements $k = 3$ | 97 | 1.8e-6 | 3.2e-7 |
| elements $k = 3$ | 193 | 2.2e-7 | 4.0e-8 |

The values at the points are more accurate than the polynomial between them, so the spline through the point values gives the better derivative. POST keeps it.

### 3.6 Reference values of the circular cases (PB3D, **measured**)

Circular HELENA equilibrium with inverse aspect ratio 0.33 (deck `tests/regression/input_circ_eps033`), f9713d7, `tol_SLEPC` $10^{-11}$, at 400, 800, 1600, 3200 points.

- Free edge: −0.4378827959, −0.4378826875, −0.4378827902, −0.4378828160. The four values differ by $3 \times 10^{-7}$ at most (relative); no $1/N$ is visible, because this eigenfunction is small at the inner end and the sum rule covers the free edge.
- Fixed edge: 7.839707e-2, 7.838917e-2, 7.838813e-2, 7.838858e-2. This eigenvalue is positive, so it does not lie below the continuous spectrum (section 4.7). It is reported in the gates and is never a condition.

## 4. The method

### 4.1 Picture

Take one Fourier mode. Now, its radial function is known only as values at the points, and the code needs two rules that do not follow from one another: a difference for the derivative and weights for the integral. The errors at a fixed edge and at the inner end, and those at a free edge before the sum rule, have one cause: near an end the two rules do not fit together, and the result is the energy of a plasma whose end lies somewhere else (section 3.1).

With elements the values define a function: on every element, the polynomial of degree $k$ through the $k + 1$ values. The derivative is the derivative of that polynomial and the integral is its integral. The discrete energy is then the energy of an allowed displacement, and an allowed displacement can only have a higher or equal lowest eigenvalue than the true eigenfunction. Nothing special happens at an end: the first element starts at `min_r_sol`, the last one stops at `max_r_sol`.

### 4.2 Formulas

The energy of chapter 7 of the theory document, for the vector $\mathbf{X}(\psi)$ of the Fourier modes and with $\mathsf{L}^j = \mathsf{P}^j - \lambda \mathsf{K}^j$:

$$\mathcal{L}[\mathbf{X}] = \int_{\psi_\text{in}}^{\psi_\text{e}} \left( \mathbf{X}^\dagger \mathsf{L}^0 \mathbf{X} + \mathbf{X}^\dagger \mathsf{L}^1 \mathbf{X}' + \mathbf{X}'^\dagger \mathsf{L}^{1\dagger} \mathbf{X} + \mathbf{X}'^\dagger \mathsf{L}^2 \mathbf{X}' \right) \mathrm{d}\psi + \mathbf{X}_\text{e}^\dagger \, \delta^\text{vac} \, \mathbf{X}_\text{e} .$$

Here $\psi_\text{in}$ is `min_r_sol`, $\psi_\text{e}$ is `max_r_sol` and $\mathbf{X}_\text{e}$ is $\mathbf{X}$ at the edge. All elements have the length $h = (\psi_\text{e} - \psi_\text{in})/n_\text{el}$, computed from the range and not as a difference of two positions, whose rounding error grows with $N$. On element $e$, with the points $a = 0 \ldots k$ of the element and the polynomial $\ell_a(\psi)$ of degree $k$ that is 1 at point $a$ and 0 at the other points of the element,

$$\mathbf{X}(\psi) = \sum_{a=0}^{k} \ell_a(\psi) \, \mathbf{X}_a , \qquad \mathbf{X}'(\psi) = \sum_{a=0}^{k} \ell_a'(\psi) \, \mathbf{X}_a .$$

The integral over the element is a sum over its Gauss points $g$ with positions $\psi_g$ and weights $w_g$ (the Gauss–Legendre weight times $h/2$). Setting the derivative of $\mathcal{L}$ with respect to every $\mathbf{X}_a^\dagger$ to zero gives $\mathsf{A} \mathbf{X} = \lambda \mathsf{B} \mathbf{X}$ with, for every pair of points $a$, $b$ of an element,

$$\mathsf{A}_{ab} \mathrel{+}= \sum_g w_g \left[ \ell_a \ell_b \, \mathsf{P}^0 + \ell_a \ell_b' \, \mathsf{P}^1 + \ell_a' \ell_b \, \mathsf{P}^{1\dagger} + \ell_a' \ell_b' \, \mathsf{P}^2 \right]_{\psi_g} ,$$

and the same for $\mathsf{B}$ with $\mathsf{K}^j$. A point shared by two elements receives terms from both.

### 4.3 The same loop as now

Call an **energy point** a position where the energy density is evaluated. The present code and the elements are two cases of one description:

| | differences (now) | elements |
|---|---|---|
| energy points | the radial points $p + 1 \ldots N$ | the Gauss points of all elements |
| weight | $\Delta$ times the weight of the point | $w_g$ |
| unknowns that the energy point involves | the points $i - p \ldots i$ | the $k + 1$ points of its element |
| $\mathbf{X}$ at the energy point | $\mathbf{X}_i$ | $\sum_a \ell_a(\psi_g) \mathbf{X}_a$ |
| $\mathbf{X}'$ at the energy point | $\sum_s c_s \mathbf{X}_{i-s} / \Delta$ | $\sum_a \ell_a'(\psi_g) \mathbf{X}_a$ |
| coefficient matrices needed | at the radial points | at the Gauss points |

`fill_mat` is a loop over energy points. For each it forms the three blocks from `V_0`, `V_1`, `V_2` and places multiples of them at the rows and columns of the unknowns involved. With elements the loop stays; the multiples change, and the `V_0` block is also placed at pairs of different points. The model script is written this way and reproduces the numbers of chapter 7 for the differences.

### 4.4 The two ends and the vacuum

- **Inner end.** The value at the first point is held at zero, with `EV_BC` as now. One point, not $p$.
- **Fixed edge.** The value at the last point is held at zero.
- **Free edge.** `vac%res` is added to the diagonal block of the last point, as now. The polynomial of the last element takes the value $\mathbf{X}_N$ at the edge, so the vacuum term is exact for the trial function. There are no weights of the last points.

### 4.5 The sliding set of Fourier modes

The plan uses rule A of section 3.3; decision 4 of section 10 asks whether it should be rule B. Rule A needs no change of `insert_block_mat`.

**One table per element.**

- `setup_modes` computes the set of every element at the middle of the element. It needs only positions and the safety factor (**read**: `X_ops.f90:1059-1240` uses `grid%n(3)` and `grid%r_F`), so it can be given the list of the middles.
- Every radial point takes the row of the element on its inner side, which is the element that the point ends or lies inside. The first point takes the row of element 1.
- Every Gauss point takes the row of its element.

The table of the radial points and the table of the Gauss points are therefore copies of rows of the same table, with the same slots and the same order of the rows. The table of the Gauss points holds `n` and `m` per Gauss point, and in `sec` the first and the last Gauss point of every mode. This matters because `interp_V` pairs the rows of its input table and of its output table by their position in the list and checks the mode numbers only in the debug build (**read**: `driver_sol.f90:366-372`). A second, independent call of `setup_modes` on the Gauss positions would give other slots without a message.

**What `fill_mat` hands to `insert_block_mat`.** The last point of the element as the point that the block belongs to. All points of the element except its first have the table of that point, so their entries stay in their slots. At the first point of the element, which has the set of the element before it, `insert_block_mat` drops the entries of a mode that has no slot there: the mode that has just entered. The script `replay_slots.py` replays `setup_local_block` for such tables (the set shifting by $+1$, $+2$, $-1$ and by `n_mod_X` between two elements; $k = 2$ and $3$): no entry lands in a slot other than its own, the pattern of the matrix is symmetric, and every unknown has a diagonal entry. If the first point of all kept a set of its own, one unknown would have no entry, which is why it takes the row of element 1.

**POST.** POST calls `setup_modes` itself (section 2), so the routine gets the degree as an argument and all four calls that build the table of the solution grid pass it: `driver_sol.f90:115` and `128`, `rich_ops.f90:140`, `driver_POST.f90:202`.

Rule B would be a mask on top of this: at a shared point the slots of the modes that only one of the two elements has are held at zero with `EV_BC`.

### 4.6 The coefficient matrices at the Gauss points

- `interp_V` already interpolates to any list of positions. The solution driver gives it the Gauss points in place of the radial points: a grid with $(k + 1) \, n_\text{el}$ positions and the table of section 4.5.
- `setup_mats` and `fill_mat` cut the coefficient arrays with indices of the solution grid (**read**: `SLEPC_ops.f90:517-519`, `709`). `solve_EV_system_SLEPC`, `setup_mats` and `fill_mat` get the grid of Gauss points as an argument and use its indices.
- With several processes, the process that owns the last point of an element assembles the whole element. The division of the Gauss points over the processes follows from that, with the radial points of the element as ghost points (width $k$); the perturbation grid that `redistribute_output_grid` hands to `interp_V` (`driver_sol.f90:192`) has to cover the Gauss points of these elements. A process with a rank at or above `sol_n_procs` gets a defined, empty case.
- `X_grid_style` 1 and 3 work this way. `X_grid_style` 2, in which the perturbation grid is the solution grid, would need the perturbation driver to compute at the Gauss points; it is stage 7 and is refused with a message until then.
- For a mode that enters, the Gauss points in the first half of its first element lie before the first perturbation point at which the mode exists, by up to half an element ($k/2$ steps). `interp_V_spline` extrapolates there, as it does now at the ends of the range of a pair. The same holds behind the last point of a mode that leaves. This is risk 2 of section 8.

### 4.7 What is known about the error

**Standard** (Strang and Fix, *An Analysis of the Finite Element Method*, 1973). The lowest discrete eigenvalue is never below the true one, the same holds for the second, the third and so on counted from the lowest, and the error is proportional to $h^{2k}$ for a smooth eigenfunction. The conditions:

1. *A fixed set of Fourier modes.* With the sliding set there is no such statement for rule A; rule B has it with respect to the problem with all modes.
2. *The kinetic energy is positive.* It is enough that the matrix $\begin{pmatrix} \mathsf{K}^0 & \mathsf{K}^1 \cr \mathsf{K}^{1\dagger} & \mathsf{K}^2 \end{pmatrix}$ is positive semi-definite at every Gauss point and $\mathsf{B}$ is positive definite. The spline interpolation works on every entry separately and does not guarantee this. Stage 3 tests it with a Cholesky factorisation of $\mathsf{B}$, not with its diagonal.
3. *The eigenvalue lies below the continuous spectrum.* The continuous spectrum consists of the $\lambda$ for which $\det(\mathsf{P}^2 - \lambda \mathsf{K}^2) = 0$ somewhere in the plasma. $\mathsf{P}^2$ vanishes for a Fourier mode at its rational surface, so with rational surfaces inside the plasma the statement covers the negative eigenvalues and no others.
4. *Exact integrals.* $k + 1$ Gauss points integrate a polynomial of degree $2k + 1$ exactly, so all four terms are exact when the coefficient matrices are constant over an element. With coefficient matrices that vary, the statement holds up to the error of the integration.

The vacuum term is part of the potential energy and is exact for the trial function, so it does not change these conditions.

**Measured in the model.** The powers 2, 4, 6 for $k = 1, 2, 3$ at both kinds of edge, with the terms $\mathsf{L}^1$, with rational surfaces, and with an eigenfunction that reaches the inner end with a slope. Every error of the elements with Gauss–Legendre integration is positive (section 3.2), for the lowest ten eigenvalues too (section 3.5).

**Spurious eigenvalues.** Under conditions 1 to 4 no discrete eigenvalue can lie below the true eigenvalue of the same number, which excludes eigenvalues that the continuous problem does not have below its continuous spectrum. This follows from the statement above. Gate 3d checks the lowest ten eigenvalues of every case.

**Not known before stage 3.** Whether PB3D's cases follow the power at 200 to 800 points. The measured 15 points per length of change of cbm18a at 200 points (section 3.1) suggests that they do for $k = 3$.

## 5. Choices, with the alternatives

| Choice | Proposed | Alternatives |
|---|---|---|
| Scheme | finite elements | **Repair the two ends inside the differences**: a rule for the first points and a rule for a fixed edge. Removes the same $1/N$; keeps the powers 2 and 3 of the differences and the errors of both signs; needs two rules that do not exist yet (theory document, section 7.16.5). |
| Functions | polynomials of degree $k$ per element, continuous (Lagrange); unknowns are values at points | **Cubic Hermite** ($X$ and $X'$ continuous; unknowns $X$ and $X'$ per point; error 3.3 times smaller than $k = 3$ at equal unknowns in the model): two unknowns per point change the solution vector, the output, POST and the tables of modes. **Cubic B-splines** (error 180 times smaller than $k = 3$ at 97 unknowns, and 210 times with splined coefficients, sections 3.2 and 3.4): the unknowns are no longer values at points, so the solution vector, POST and the rule for the sliding set change, a B-spline belonging to four elements. |
| Points inside an element | equally spaced | Gauss–Lobatto positions give the same eigenvalues (section 3.5) and a grid that is not equally spaced. |
| Integration | Gauss–Legendre, $k + 1$ points per element | **Gauss–Lobatto at the points themselves** (the panels of issue #35): needs no new positions for the coefficient matrices; same power, errors 0.3 to 1.6 times those of Gauss–Legendre in the model; the sign of the error is not fixed. Possible as a later switch; not in this plan. **More Gauss points** ($k + 3$) help when the coefficients are a spline (section 3.4); the number of points is an argument of the new routine, so this stays open at no cost. |
| Degree | input, 1 to 3, default 3 | A degree above 3 widens the band of the matrices beyond the present 7 blocks. Whether it would pay is not known; section 3.4 shows that for $k = 3$ the interpolation of the coefficients already limits the accuracy near 120 unknowns in the model. |
| Sliding set | rule A | rule B; section 3.3 and decision 4 |
| Unknowns | values at $N = k \, n_\text{el} + 1$ equally spaced points | — |

## 6. Changes in the code

**Read** in the code where a line is given; the rest is the plan.

| File | Change |
|---|---|
| `num_utilities.f90` | New: Gauss–Legendre points and weights for any number of points. New: the values and derivatives of the $k + 1$ polynomials at given positions of the reference element. Removed in stage 6: `calc_end_weights` and `calc_coeff_fin_diff`, which are used only by `fill_mat`, by `calc_end_weights` itself and by the unit tests. |
| `SLEPC_ops.f90` `fill_mat` (605) | Stage 1: the loop in the form of section 4.3, with a table. Stage 3: the table of the elements; the element length from the range; no `w_end`, no `ndc`. |
| `SLEPC_ops.f90` `solve_EV_system_SLEPC`, `setup_mats` (410), `fill_mat` | The grid of Gauss points as an argument; the cuts of the coefficient arrays (517-519, 709) use its indices. |
| `SLEPC_ops.f90` `set_nonzeros` (864) | Already counts the blocks on both sides (884-895). Only the width changes from $p$ to $k$. |
| `SLEPC_ops.f90` `set_BC` (989), `set_bulk_lims` (564) | One held point at the inner end. `set_BC_1` already clears as many neighbouring blocks as the width of the difference, which becomes $k$. |
| `SLEPC_ops.f90:97` `step_size` | Replaced by the element length; the line that prints `E_pot` with it (1755) follows. |
| `SLEPC_ops.f90`, new checks | The norm of a matrix minus its Hermitian transpose, divided by the norm of the matrix, with a threshold; a Cholesky factorisation of $\mathsf{B}$. The existing `test_mat_hermiticity` cannot serve: the switches `debug_setup_mats` and `debug_set_BC` that reach it are never set, and it gives no verdict. |
| `driver_sol.f90` (175-203) | The grid handed to `interp_V` is the grid of Gauss points; `X_sol` lives on it. |
| `grid_ops.f90` | New: positions of the Gauss points and their division over the processes. `calc_norm_range_PB3D_sol` (426-427): ghost width $k$. `calc_norm_range_POST` (496-498): ghost widths no longer taken from `norm_disc_prec_sol` (stage 5). |
| `X_ops.f90` `setup_modes` (1059) | An argument for the degree; one set per element at its middle; the rows copied to the radial points and to the Gauss points (section 4.5). Called for the solution grid at `driver_sol.f90:115`, `128`, `rich_ops.f90:140`, `driver_POST.f90:202`. |
| `X_utilities.f90` `trim_modes` (328) | The test checks `id_lim_o`, which is set two lines above, where it means `id_lim_i`: a mode of the output table that the input table lacks is not caught. Corrected in stage 3, because the Gauss points are a new output table. |
| `input_ops.f90` (453, 978-988) | `norm_disc_prec_sol` is the degree, 1 to 3. The two minimum sizes are applied first; then $N$ is raised to the next $k \, n_\text{el} + 1$ with a message. |
| `sol_utilities.f90` (184-199), `sol_ops.f90` (494, 725-730) | Spline orders of POST are taken from `norm_disc_prec_sol` now; a degree of 2 would silently select the Akima spline (`num_utilities.f90:259`). They get an order of their own (stage 5). |
| `sol_utilities.f90` `calc_tot_sol_vec` (663) | Unchanged: the spline through the point values stays (section 3.5). |
| `rich_ops.f90` | 257-258: the printed expected power. 376, and `driver_POST.f90:433`: the memory estimates count the Gauss points. 472: not touched; decision 7. |
| `num_vars.f90:60-61` | `prog_version`, `min_PB3D_version`. |
| tests | Unit tests for the two new routines; in stage 6 they replace those of `calc_end_weights` and `calc_coeff_fin_diff`. `unit_setup_modes_sim` gets a case with elements. Anchors re-recorded in stage 3. A convergence guard per deck. |
| docs | `testing.md`, `README.md`, Doxygen headers of `fill_mat`, `set_BC` and `setup_modes`. |

Not changed: `insert_block_mat`, the solution vector, `print_output_sol`, `reconstruct_PB3D_sol`, `setup_guess`, `store_results` (apart from the printed energy), the vacuum, everything before the solution driver.

## 7. Stages and gates

A stage is finished when its gate passes. Each stage is one or more commits on this branch, pushed, with the gate result as a comment on the issue. No stage starts before the gate of the one before it has passed.

**How two builds are compared.** The matrices are written with the solver options `-eps_view_mat0 binary:A.dat -eps_view_mat1 binary:B.dat` on the command line of PB3D and compared entry by entry with `scripts/snellius/compare_petsc_mats.py`. **Measured** for f9713d7 against itself (cbm18a fixed, 100 points, `tol_SLEPC` $10^{-10}$, two launches on the Mac): none of the 275 200 entries of $\mathsf{A}$ or of $\mathsf{B}$ differs by more than $10^{-12}$ relative, and the eigenvalue is the same in all 17 digits. All gates use `tol_SLEPC` $= 10^{-10}$: the decks have $10^{-6}$, which is too loose for these comparisons, and below about $10^{-11}$ PB3D marks the eigenvalue as complex (`SLEPC_ops.f90:1645`), writes its line with a `#` and does not store the eigenvector.

**Stage 0. Diagnosis and model.** Done: sections 3 and 12.

**Stage 1. The loop of `fill_mat` in the form of section 4.3, with the table of the differences.** The numbers must not change.

*Gate 1.* For the regression decks (cbm18a fixed and free, circular with inverse aspect ratio 0.33 and its VMEC twin, circular with 0.05 free), each with `norm_disc_prec_sol` = 1, 2, 3: no entry of $\mathsf{A}$ or $\mathsf{B}$ differs by more than $10^{-12}$ relative from f9713d7, and the eigenvalues agree to $10^{-10}$ relative. On the Mac with one process; on Snellius with 1, 2 and 4 processes. CI green without re-recording an anchor.

**Stage 2. The two routines of `num_utilities`, with unit tests.**

*Gate 2.* Gauss–Legendre with $n$ points integrates $\psi^{2n-1}$ exactly and $\psi^{2n}$ not. The polynomials for equally spaced points add up to 1 and their derivatives to 0 at every Gauss point. The derivative of $\psi^k$ is exact. The tables for $k = 1, 2, 3$ equal those of `fe_model.table_elements` with equally spaced points to $10^{-14}$; the numbers are written in the Fortran test.

**Stage 3. Elements.** One table of modes per element and its copies; the grid of Gauss points; `interp_V` to it; the table of the elements in `fill_mat`; one held point; the input checks; the degree as an argument of `setup_modes` at all four calls; `trim_modes`; the version; the two new checks; the anchors re-recorded with a table "before, after, reason" in `testing.md`. From this stage on the elements are the scheme of the branch. No switch keeps the differences, because no gate needs one: every comparison is with the separate build of f9713d7.

*Gate 3a, the new parts with numbers that must not change.* With $k = 1$, a fixed set of modes and one energy point at the outer end of each element, the table of the elements is the table of the differences with $p = 1$ and all weights 1 (section 3.5). In that setting the matrices of cbm18a fixed, of the circular fixed case and of its VMEC twin equal those of f9713d7 with `norm_disc_prec_sol` = 1 to $10^{-12}$ per entry. This tests the grid of Gauss points, its table of modes, the interpolation to it, the loop and the preallocation in the real code. The outer-end rule exists in the code for this gate only and is deleted in stage 6. A fixed set is needed because rule A takes the set at the middle of an element, which is not the present table.

*Gate 3b, the matrices for $k = 2, 3$.* For cbm18a fixed with about 50 points, once with a fixed set and once with the sliding set, $\mathsf{A}$ and $\mathsf{B}$ equal to $10^{-12}$ the matrices that the loop of `fe_model.assemble` builds from PB3D's own coefficient matrices at the Gauss points, which a debug switch writes to a file. This tests the polynomials, the weights, the placement and, with the sliding set, the dropped entries.

*Gate 3c, Hermitian and positive.* For both edges and both values of `X_style`: the norm of $\mathsf{A} - \mathsf{A}^\dagger$ over the norm of $\mathsf{A}$ below $10^{-12}$, the same for $\mathsf{B}$; the Cholesky factorisation of $\mathsf{B}$ succeeds.

*Gate 3d, the convergence series.* Series of 100 to 3200 points (the nearest $k \, n_\text{el} + 1$) for $k = 1, 2, 3$ on: cbm18a fixed with the complete fixed set $m = 5 \ldots 47$, cbm18a fixed with the sliding set, cbm18a free, circular 0.33 free and fixed, circular 0.05 free. To pass:

- cbm18a fixed, complete fixed set, $k = 3$: every value from 400 points on within $2 \times 10^{-5}$ relative of $-4.38175 \times 10^{-2}$, the limit of the differences (section 3.1). The differences are $1.6 \times 10^{-3}$ away from it at 800 points. For $k = 1$ and $2$ the limit of a fit agrees with it to $5 \times 10^{-5}$.
- cbm18a fixed, sliding set, $k = 3$: every value from 400 points on within $1 \times 10^{-4}$ relative of $-4.3800 \times 10^{-2}$. The margin is wider because the error of the sliding set itself changes with $N$: in the present code the distance between the sliding set and the fixed set $5 \ldots 40$ at 400 points is $5 \times 10^{-5}$ of the eigenvalue away from its limit.
- circular 0.33 free: the limit equals $-0.4378828$ to $2 \times 10^{-6}$ relative (section 3.6).
- With a fixed set of modes and $k = 2, 3$ the values fall towards their limit from above. Only differences above $10^{-9}$ relative count.
- With `n_sol_requested = 10`: every negative eigenvalue among the lowest ten at 400 points lies within $10^{-3}$ relative of one of the differences at 3200 points, and there is no extra one.
- Positive eigenvalues, such as the circular fixed case, are reported and are never a condition.
- Reported with the gate: the measured power per $k$ on the circular free case with a fixed set, and the largest difference between the eigenvector at the radial points and that of f9713d7 at 3200 points, both scaled to a largest value of 1.

If a case does not follow its power, the gate asks for the reason (equilibrium points, tolerance of the solver, sliding set), found by a run that changes only that.

*Gate 3e, cost.* Run time and memory of the 800-point cbm18a deck within 1.5 times those of f9713d7.

*Gate 3f, the debug build.* cbm18a fixed and free run in the debug build, by hand, without a message from its checks.

**Stage 4. Several processes.**

*Gate 4.* On Snellius, with 1, 2 and 4 processes on cbm18a fixed and on the circular free case: the matrices agree to $10^{-12}$ per entry (`compare_petsc_mats.py`, which was written for this comparison) and the eigenvalues to $10^{-7}$ relative. The present code shifts by $1.5 \times 10^{-8}$ between 1, 2 and 4 processes, so a tighter number would test the solver and not the assembly. CI does not cover this: its rows with the label `mpi` do not assemble $\mathsf{A}$ or $\mathsf{B}$ with several processes unless the build has `-DPB3D_TEST_MPI_REGRESSION=ON`, which CI does not set, and the Mac cannot run `mpirun`.

**Stage 5. POST.** The spline orders and the ghost widths of POST no longer taken from `norm_disc_prec_sol`. The derivative of the solution stays the spline through the point values.

*Gate 5.* POST in the debug build, where `calc_XUQ_arr` checks the mode numbers (`sol_utilities.f90:346-351`), runs without a message on the outputs of stage 3 for cbm18a fixed with the sliding set and for the circular free case. The plotted $X$ at the radial points equals the stored eigenvector. Against POST of f9713d7 on its own output of the same case, the largest difference of $|X|$, both scaled to a largest value of 1, is below $10^{-2}$ and is reported. The ratio of the reconstructed energies against the eigenvalue is reported next to that of f9713d7 (issue #26; not a condition, because an energy consistent with the elements would need $X$ and $X'$ at the Gauss points, which is not part of this plan).

**Stage 6. The finite-difference routines are deleted.** Deleted: `calc_end_weights`, `calc_coeff_fin_diff`, their unit tests, the outer-end rule and the debug output of gate 3b, the comments that describe weights. Added: each free-boundary deck gets the 200-against-400-points guard that `circ_eps005_free_bnd.sh` has. Updated: `testing.md`, `README.md`.

*Gate 6.* CTest and CI green. No occurrence of `calc_end_weights`, `calc_coeff_fin_diff`, `w_end` or `ndc` left in `Modules/` and `tests/unit/`. The Doxygen warning count not above its baseline.

**Stage 7 (optional). `X_grid_style` 2.** The perturbation driver computes at the Gauss points.

*Gate 7.* cbm18a fixed with `X_grid_style` 2 within $10^{-5}$ of `X_grid_style` 3 at 800 points.

**After the branch.** The comparisons with TERPSICHORE and MISHKA-1 of `vacuum_work_overview.md` section 4, rerun; chapter 7 of the theory document.

Sizes, as guesses: stages 1 and 2 a day together; stage 3 three to four days including the series; stage 4 half a day once Snellius is reachable; stage 5 a day; stage 6 half a day.

## 8. Risks and open points

| | Risk or open point | What the plan does |
|---|---|---|
| 1 | A pair of Fourier modes exists on few perturbation points; `interp_V_spline` then copies or interpolates linearly. | Gate 3a compares at identical positions. Gate 3d with the fixed set has no such pairs. The debug statistics of `interp_V` are compared before and after. |
| 2 | For a mode that enters or leaves, Gauss points lie up to half an element outside the perturbation points that have the mode, and the coefficients are extrapolated there (section 4.6). | Gate 3b with the sliding set compares the matrices; gate 3d compares the sliding set with the complete fixed set. If the extrapolation shows in gate 3d, the set of an element can be limited to the modes that the perturbation grid has over the whole element; that is decided then. |
| 3 | The irregular changes of the sliding set (section 3.3) hide the power of the error. | The series are run with the complete fixed set too. |
| 4 | PB3D's cases need more points than the model before the error follows its power. cbm18a has 15 points per length of change at 200 points; below about 5 the elements are worse than the differences (section 3.2). | Gate 3d reports the range. The default degree is decided after it. |
| 5 | An unknown that receives no energy makes both matrices singular. | With rule A every unknown has a diagonal entry (`replay_slots.py`). Gate 3c factorises $\mathsf{B}$. Rule B needs its unused slots held with `EV_BC`. |
| 6 | The interpolated $\mathsf{K}^j$ are not positive semi-definite at some Gauss point (section 4.7, condition 2). | Gate 3c. If it fails, the Gauss point and the pair of modes are reported before anything else is done. |
| 7 | The solver: shift-and-invert on matrices with blocks on both sides of the diagonal. The band is as wide as now for $k = 3$. | Gate 3e. |
| 8 | Several processes: an element whose points belong to two processes. PETSc accepts values for rows of another process and moves them at assembly, as the present code relies on for its differences. | Gate 4, with the matrices compared. |
| 9 | With `sol_n_procs` smaller than `n_procs`, the last energy point and `vac%res` look to be added once more per extra process in the present code (**read**, not run). | Checked in stage 4 before the elements are compared; reported in the issue if confirmed. |
| 10 | Snellius needs Toon's IP registration and the 1Password agent. | Stages 1 to 3 run on the Mac and CI, apart from the Snellius part of gate 1; gate 4 waits for access. |
| 11 | A safety factor that is not monotonic, or a set that shifts by two between elements. | `replay_slots.py` covers shifts of $+2$, $-1$ and `n_mod_X`; the tables of `setup_modes` handle both (**read**: `X_ops.f90:1150-1212`); `unit_setup_modes_sim` gets a case with elements. |

## 9. What this does to the open items of chapter 7

Section 7.16.5 of the theory document lists what is not settled for the differences. With elements: the fixed edge and the inner end need no rule; an error that falls faster than $1/N^3$ is the normal case for $k \ge 2$; the exponent of the Richardson extrapolation along the field line is a separate matter (issues #12, #43).

## 10. Decisions for Toon

| | Decision | Recommended | If otherwise |
|---|---|---|---|
| 1 | Functions | Lagrange, values at the points. **Judgement:** on cbm18a the sliding set of 20 modes costs $4 \times 10^{-4}$, and the model at the same resolution gives $2.5 \times 10^{-6}$ for $k = 3$ at 200 points; the 3 to 200 times of Hermite or B-splines would not be seen. | Hermite or B-splines: stages 3 and 5 and the output format grow; say so before stage 2. |
| 2 | Name of the input | keep `norm_disc_prec_sol`, now the degree | a new name: one more entry in `misc_sol` |
| 3 | Default degree | 3 (same band of the matrices as now), confirmed or changed after gate 3d | 2 if gate 3d shows no gain of 3 |
| 4 | Sliding set: rule A or rule B | **Open.** Rule A keeps the limit of the present code ($-4.3800 \times 10^{-2}$ on cbm18a), needs no change of `insert_block_mat`, and makes gate 3d a comparison with a known number; its error has no fixed sign. Rule B gives an eigenvalue above the one with all modes that can only fall when `n_mod_X` grows, which makes `n_mod_X` a parameter with a convergence test; its error is of the same size, larger or smaller depending on the case (section 3.3). I would start with rule A for the comparison and add rule B as an option after gate 3d. | Rule B from the start: the mask and the held slots come into stage 3, and all `X_style` 2 values move by the amounts of section 3.3. |
| 5 | `X_grid_style` 2 | stage 7, refused with a message until then | remove the style |
| 6 | Where the series run | Mac for stages 1 to 3, Snellius for the several-process parts of gates 1 and 4 | all on Snellius: needs access first |
| 7 | `calc_rich_ex` uses `norm_disc_prec_sol` as the exponent for the points along the field line (`rich_ops.f90:472`) | leave to issue #43; the decks have `max_it_rich` = 1 | decide #43 first |
| 8 | Version | 2.51, `min_PB3D_version` with it | — |
| 9 | Chapter 7 of the theory document | after the branch is accepted | — |

## 11. Seen while reading, not part of this plan

Reported by the review of the plan against the code, **read** and not run. I confirmed the items marked (checked) by reading the lines myself.

- `test_mat_hermiticity` cannot be reached: `debug_setup_mats` and `debug_set_BC` are never set (checked).
- `trim_modes` tests `id_lim_o` where it means `id_lim_i` (`X_utilities.f90:328`; checked). This one is corrected in stage 3.
- The comment in `set_BC_2` about a small non-Hermiticity of the vacuum term is out of date: `vac%res` is made Hermitian at `vac_ops.f90:2021` (checked).
- `Examples/input_cdxu` sets `min_n_r_sol`, which is not an input any more (checked).
- `store_results` (`SLEPC_ops.f90:1755`) multiplies $\mathbf{X}^\dagger \mathsf{A} \mathbf{X}$ by `step_size`, which $\mathsf{A}$ already contains; the number is printed only.
- `calc_tot_sol_vec`: line 663 stops POST with `X_style` 1 on several processes or with a plot range that is not the full one; lines 691-704 use indices of the local grid on positions of the total grid, which is harmless while the grid is equally spaced.
- With `sol_n_procs` smaller than `n_procs`: risk 9 of section 8.
- No regression test runs POST.

## 12. Reproducing the numbers

All scripts are in `Documentation/radial_finite_elements/` and need numpy and scipy (`eigenfunction_profile.py` also h5py).

```bash
cd Documentation/radial_finite_elements && ~/Code/desc-vs-vmec/.venv/bin/python -B fe_model.py
```

runs every section of the model, which takes about half an hour; `sets` is the longest. Sections can be named on the command line:

| section | what it gives | used in |
|---|---|---|
| `one` | one Fourier mode: free edge, fixed edge, $V = 0$ | 3.2 |
| `six` | six Fourier modes: free and fixed edge | 3.2 |
| `hard` | twenty Fourier modes | 1, 3.2 |
| `ends` | where the differences place the two ends | 3.1 |
| `sets` | the sliding set, rules A and B, three variants | 3.3 |
| `spline` | coefficients interpolated from an equilibrium grid | 3.4 |
| `gauss` | number of Gauss points | 3.5 |
| `ten` | the lowest ten eigenvalues | 3.5 |
| `cost` | size of the matrices | 3.5 |
| `same` | two tables with the same matrices | 3.5 |
| `post` | the derivative at the radial points | 3.5 |
| `lengths` | how fast the eigenfunctions of the model change | 3.1 |

`cbm18a_ends.py` holds the PB3D runs of section 3.1 and prints the fits, the two slopes and the prediction. `replay_slots.py` replays `setup_local_block` for the tables of rule A (section 4.5). `eigenfunction_profile.py <PB3D_out.h5>` prints the three lengths of section 3.1 and the shape near the ends of the eigenfunction of a PB3D run.

PB3D series of section 3.1: build this branch (Release, without STRUMPACK, recipe in `Documentation/testing.md`), then run the cbm18a fixed-boundary deck with `n_r_sol` = 200 to 3200, `norm_disc_prec_sol` = 3 or 2 and `tol_SLEPC` = 1.0E-12, once as it is, and once with `n_mod_X` removed and `min_sec_X = 5`, `max_sec_X = 40` or `47`. With this tolerance PB3D writes the line of the eigenvalue with a `#` in `PB3D_out_EV_R_1.txt` and does not store the eigenvector; the run for `eigenfunction_profile.py` used `tol_SLEPC` = 1.0E-9. The decks have leftover lines after the closing `/` of the namelist, which the read ignores and an automatic replacement of a key would hit.
