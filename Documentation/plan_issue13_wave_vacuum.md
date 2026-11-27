# Plan: wave amplitudes as unknowns in the field-line vacuum

Issue #13. The method is section B.12 of the theory document, "A form with wave amplitudes as unknowns". Draft of 27 November 2026 for Toon's review, written against master 31e0354 (version 2.58). Nothing is implemented, and no stage starts before Toon says so. The numbers of sections, equations and tables are those of the theory document that has the section on sums over rows as B.13.

## 1. What is asked

The field-line vacuum (`vac%style` 1) solves $\mathsf{H}\boldsymbol\phi = \mathsf{G}\boldsymbol\phi'$ for the potential $\phi$ at each of the $N$ points of the edge, where $\phi'$ is the normal derivative that the displacement imposes. Section B.12 solves for the amplitudes of $M'$ waves instead, and builds the $M' \times M'$ system from sums over pairs of points, so that no $N \times N$ matrix is stored. This plan puts that form into PB3D as a second way to compute `vac%res`, selected by an input, next to the form with point values, which stays the default.

The plan has five stages. Each stage has a gate with numbers that are fixed here, before any run. Stage 0 measures what nobody has measured and what the other stages depend on; no change to the code is kept from it.

## 2. Decisions for Toon

1. **Start stage 0?** It needs six one-process jobs on Snellius, the largest with 10 GB, and half a day.
2. **Names and defaults of the inputs.** Proposed: `vac_sol_style` (1: point values, the default; 2: wave amplitudes), `n_pol_vac` and `n_tor_vac` for the size of the set of waves (section 4). The defaults of the two numbers come from stage 0.
3. **Stage 2 (reuse of sums from one Richardson level to the next) at all?** It saves at most a quarter of the sums. My recommendation: decide after stage 1, with measured times.
4. **The copy of G and H between levels in the form with point values.** `interlaced_vac_copy` copies the previous G and H into the new ones with four matrix products, and `calc_GH` then sets both to zero. By operation count the four products cost $3N^3$, against $\tfrac{2}{3}N^3$ for the LU. Measured on 27 November with the serial routines of OpenBLAS on one thread, at $N = 6432$ (32 lines, 101 points per line at the previous level and 201 at the new one): the four products take 14.0 s and the LU takes 3.3 s, a ratio of 4.3 (4.1 at $N = 3216$). The matrices that select the points hold a single 1 per column, but the products take the same time as with full matrices. Inside PB3D, with `pdgemm` on several processes, this is not measured. The estimate for `max_tot_mem` counts the copy since version 2.58: 1.96 GB at the second level of the tutorial case. `calc_vac_res` also computes the product `GEP`, which nothing reads. My recommendation: remove the copy of G and H, its part of the estimate, and `GEP`, in a small pull request of its own before stage 1. No number changes.
5. **Rows of one field period only.** This is not in B.12. When a turn by one field period maps the field lines of the grid onto each other, every row of the sums has $N_\text{fp} - 1$ copies, and the rows of one field period are enough (derived, not run). For an axisymmetric edge every line is a copy of the first one, so the rows of one line are enough: the sums then cost $n_\alpha$ times less. For W7-X it needs a number of lines that is a multiple of 5, while PB3D warns to choose `n_alpha` coprime to the number of field periods. Stage 0 checks the statement on the matrices of the twin at no cost. Not planned unless Toon wants it.
6. **Stage 3 (POST) now, or after the work of issue #81.**
7. **The issue.** Its label is `exploration` and its title says "test … on the W7-X boundary". I changed neither. Whether this file stays in the repository after the review, or goes into the issue, is also for Toon to decide.

## 3. What the form with wave amplitudes gives, and what it does not give

| | point values as unknowns (the code now) | wave amplitudes as unknowns (B.12) |
|---|---|---|
| memory in PB3D | Three real $N \times N$ matrices at the solve: `vac%G`, `vac%H`, and the copy of H that `solve_Phi_BEM` factorizes. They are 63 GB at $N = 51\,200$. After the first level the peak is higher, $3.5\,N^2$ numbers, inside `interlaced_vac_copy`. | Two blocks of 64 rows of G and H per process (52 MB together at $N = 51\,200$, from the array sizes) and matrices of the size $M' \times M'$. The Python prototype needed 1.2 GB for the sums at $N = 51\,264$. |
| memory in POST | `reconstruct_PB3D_vac` allocates G and H for every free-boundary run, also when no vacuum plot is asked. | Stage 1 allocates G and H where they are computed, so POST allocates them only for the vacuum plot. |
| work that grows as $N^3$ | The LU factorization of H. | None. The system that is solved has the size $M' \times M'$. |
| work that grows as $N^2$ | The entries of G and H. | The same entries, and the product of every row with the $M'$ waves. By operation count the products cost $4N^2M'$, a fifth of the LU for $M' = 500$ and $N = 15\,488$. |
| a new Richardson level | Everything is computed again. | The sums over pairs of two old points are reused: three quarters of the pairs are computed. The prototype measured 0.69 to 0.74 of the time. |
| accuracy | The error of the sums over the points. On the twin at 64 × 100 points the diagonal elements of the vacuum matrix are 0.98 to 0.99 of those of the axisymmetric version. | The same entries, so the same error of the sums. The result differs from the form with point values by the part of the potential that the $M'$ waves do not hold: $4\cdot10^{-7}$ to $4\cdot10^{-6}$ of the vacuum energy on the torus of the prototype. |

The numbers of the prototype are those of section B.12 of the theory document. The run of the twin on 128 lines × 400 points took 220 GB in total; what the 157 GB beyond the three matrices hold is not recorded.

What it does not give:

- The error of the vacuum does not start to fall as a power of the step along the field lines (issue #11). The entries are the same as now, so the Richardson extrapolation of a free-boundary eigenvalue is as valid or as invalid as it is now.
- Nothing changes for HELENA (`vac%style` 2), where the number of points is hundreds.

What is not known, and what stage 0 measures:

- **The number of waves.** On the circular torus of the prototype 3 extra poloidal mode numbers on each side are enough when $\phi'$ is one wave: the energy then differs from that of the form with point values by $2\cdot10^{-7}$ to $1.2\cdot10^{-6}$. On the cut cbm18a the comment of 24 November on the issue indicates about 10. For an edge that is not axisymmetric there is one prototype (section B.13): an ellipse that turns three times in one toroidal turn, with 15 toroidal mode numbers times 11 poloidal ones, 165 waves, for one wave as the source. For an edge from an equilibrium nothing is measured.
- **The difference between the two forms when both use PB3D's own entries.** The prototype of section B.12 had no rule for intervals that are close, and it multiplied the error at the point $i$ by the area of its piece of the edge, $W_i\,\mathcal{J}|\nabla\psi|$, where equation B.91 has the weight $W_i$ alone.
- **How the time of the vacuum divides between the entries and the LU now.** W7-X on 128 × 121 points takes 13 minutes on one core. Which part of that is the LU is not recorded, and the form with wave amplitudes removes only that part.

## 4. The method in the terms of the code

**The points.** The edge has $n_\alpha$ field lines (`n_alpha`) with $n_\text{par}$ points each, $N = n_\alpha n_\text{par}$ points in all. The point $j$ has the field-line label $\alpha_j$, the parallel coordinate $\theta_j$, the position $\mathbf{r}_j$ (`vac%x_vec`), the vector $\mathbf{n}_j = -\mathcal{J}\nabla\psi$ (`vac%norm`) and the weight $W_j = \Delta\alpha\,\Delta\theta$, half of that for the first and the last point of a line (equation B.22). $q$ is the safety factor at the edge (`vac%jq`), and $N_\text{fp}$ the number of field periods of the equilibrium.

**The waves.** A wave has a toroidal mode number $n'$ and a poloidal mode number $m'$, and the value

$$[\mathbf{e}]_j = e^{i\left(n'\alpha_j + (n'q - m')\theta_j\right)}$$

at the point $j$, which is the expression of `set_up_data` in `calc_vac_res`. The plasma has $M$ waves (`n_mod_X`): all have $n' = n$ (`prim_X`), and their $m'$ are the poloidal mode numbers at the edge, $m_\text{min} \ldots m_\text{max}$ (`vac%sec_X`). $\mathsf{E}$ is the $N \times M$ matrix with these waves as its columns. The potential is written as a sum of $M'$ waves, the columns of the $N \times M'$ matrix $\mathsf{E}'$. The set is proposed as

$$n' = n + pN_\text{fp},\quad p = -P \ldots P, \qquad m' = m_\text{min} - K \ldots m_\text{max} + K,$$

so $M' = (2P+1)(M+2K)$, with the two inputs $K$ (`n_pol_vac`) and $P$ (`n_tor_vac`). The toroidal mode numbers follow from a symmetry: the edge is the same after a turn by one field period, so the exact map from $\phi'$ to $\phi$ turns a wave with $n$ into waves with $n + pN_\text{fp}$ only. For an axisymmetric edge $P = 0$. For the poloidal mode numbers there is no such rule, and the range above, the same for every $p$, is a first proposal only. In the prototype with the turning ellipse the range moves by 2 for every step of $p$, because the shape of that edge depends on $2\theta - N_\text{fp}\zeta$: its waves are $n' = n + pN_\text{fp}$ with $m'$ around $m + 2p$. Stage 0 measures which waves the potential holds on W7-X, and the shape of the set, with such a shift per $p$ if needed, follows from that measurement.

**The sums.** With the entries $H_{ij}$ and $G_{ij}$ of the code,

$$A_{lk} = \sum_{i=1}^N W_i\,[\mathbf{e}_l]_i^*\left(\sum_{j\ne i} H_{ij}\,[\mathbf{e}_k]_j + H_{ii}\,[\mathbf{e}_k]_i\right),$$

$$B_{lm} = \sum_{i=1}^N W_i\,[\mathbf{e}_l]_i^*\left(\sum_{j\ne i} G_{ij}\,[\mathbf{e}_m]_j + G_{ii}\,[\mathbf{e}_m]_i\right),$$

where $k$ and $l$ run over the $M'$ waves and $m$ over the $M$ waves of the plasma. In matrix form $\mathsf{A} = \mathsf{E}'^\dagger\mathsf{W}\mathsf{H}\mathsf{E}'$ and $\mathsf{B} = \mathsf{E}'^\dagger\mathsf{W}\mathsf{G}\mathsf{E}$, with $\mathsf{W}$ the diagonal matrix of the weights. The diagonal entries are those of section B.5 (`vac_diag_style` 1):

$$H_{ii} = -4\pi - \sum_{j\ne i} H_{ij}, \qquad G_{ii} = \frac{\mathbf{n}_i}{|\mathbf{n}_i|^2}\cdot\left[\sum_{j\ne i} H_{ij}(\mathbf{r}_j - \mathbf{r}_i) - \sum_{j\ne i} G_{ij}\,\mathbf{n}_j\right].$$

They need seven numbers per point, which are sums over the same row.

**The vacuum matrix.**

$$\boldsymbol\delta^\text{vac} = \frac{1}{\mu_0}\,\mathsf{F}\left(\mathsf{E}^\dagger\mathsf{W}\mathsf{E}'\right)\mathsf{A}^{-1}\mathsf{B}\,\mathsf{F},$$

with $\mathsf{F}$ the diagonal matrix of the numbers $nq - m$, and $\mathsf{E}^\dagger\mathsf{W}\mathsf{E}'$ a single sum over the points. The last steps are those of `finish_response`: the division by $\mu_0$, the factor $2\pi$ that belongs to the code's weight $\Delta\theta/n_\alpha$, and the Hermitian part.

**What one process does.** The rows $i$ are divided evenly over all MPI processes. A process takes 64 of its rows at a time and computes their entries for all $j$ with `calc_GH_1` itself. That routine already fills the rows it is given into the arrays it is given; only its columns are tied to the ScaLAPACK layout, and one new argument for the column limits removes that. The process then collects the seven sums of every row, sets the two diagonal entries, multiplies the block with the waves, adds the result to its $\mathsf{A}$ and $\mathsf{B}$, and discards the block. A wave is the product of a factor per line and a factor per point along the line, so two small tables hold all wave values. At the end the matrices of the processes are added, and the last process solves the $M' \times M'$ system with LAPACK (`zgesv`). ScaLAPACK is not used in this form.

## 5. Stages and gates

| stage | what | what it gives | new Fortran (my guess) |
|---|---|---|---|
| 0 | measure on PB3D's own G and H, in Python | number of waves, difference of the two forms, share of the LU | none kept |
| 1 | one level without an $N \times N$ matrix | memory $\propto N$, no LU | 300 lines, 300 of tests |
| 2 | reuse of sums from level to level | three quarters of the sums per level | 300 lines |
| 3 | POST from the amplitudes | vacuum plot without an $N \times N$ matrix | 200 lines |
| 4 | documentation, default | | |

One pull request per stage, each into master with its README entry and version number, each merged by Toon.

**Quantities used in the gates.**

- $\boldsymbol\delta^\text{vac}$ is `vac%res`, the $M \times M$ vacuum matrix. Its diagonal elements are positive.
- The distance between two vacuum matrices $\mathsf{X}$ and $\mathsf{Y}$ is $d(\mathsf{X},\mathsf{Y}) = \max_{k,m}\,|X_{km} - Y_{km}| \,/ \sqrt{Y_{kk}Y_{mm}}$. The square root makes the modes near resonance, where $nq - m$ is small, count as much as the others.
- "The eigenvalue" is the lowest $\omega^2/\omega_\text{A}^2$ of the run, in PB3D's normalization. $\Delta_\text{vac}$ is the eigenvalue of a free-boundary run minus the eigenvalue of the same run with a fixed boundary, which is the part that the vacuum causes.
- "The twin" is the circular VMEC case of the regression test `circ_eps033_vmec_twin` ($n = 10$, 10 poloidal mode numbers). "W7-X" is `wout_w7x.03.00llfix_99.nc` with $n = 10$ and $m = -14 \ldots -7$ on 128 lines × 121 points ($N = 15\,488$).
- "The exact potential" is built by the script of stage 0: the potential of point sources on the flux surface with the normalized toroidal flux 0.64 (0.8 of the minor radius), with the strengths $e^{i(n\zeta - m_0\theta)}$ in the angles of VMEC, where $m_0$ is the middle poloidal mode number of the plasma. It is harmonic outside the plasma, and its value and its normal derivative at the points of the edge follow from the formula of a point source. Both forms get the exact $\phi'$ and solve for $\phi$; the error is that of the vacuum energy $\sum_j W_j\,\phi'^*_j\phi_j$ against the same sum with the exact $\phi$. For this potential the set of waves is chosen from its own amplitudes.

### Stage 0. Measure on the matrices that PB3D has now

*Purpose.* Find the number of waves for a three-dimensional edge, the difference between the two forms with PB3D's own entries, and the share of the LU in the time, before any Fortran is written.

*What is done.* A temporary write in `calc_vac_res`, not committed, puts `vac%H`, `vac%G` and the times of `calc_GH` and `solve_Phi_BEM` into files. On one process the two arrays are the whole matrices. The positions, the vectors $\mathbf{n}_j$, the angles, the mode numbers and `vac%res` are in `PB3D_out.h5` already. A Python script forms $\mathsf{A}$ and $\mathsf{B}$ as matrix products for any set of waves, and compares the vacuum matrix with `vac%res` of the same run. The entries are PB3D's, with the rule for intervals that are close and the diagonal from known potentials, so nothing of PB3D is ported to Python. This takes the place of points 1 and 2 of the list "Next, on W7-X" of the issue.

Three cases, each with a free and with a fixed boundary, as six one-process jobs on Snellius at the same time: the twin on 64 × 100 points (0.33 GB per matrix); the tutorial case, the cut cbm18a from VMEC, on 64 × 100 points, as a shaped axisymmetric edge; W7-X on 128 × 121 points (1.9 GB per matrix; the run took 9.7 GB). The runs with a fixed boundary give $\Delta_\text{vac}$ and the memory and time without the vacuum.

The script reports:

1. For every plasma mode, the amplitude of every wave $(n', m')$ in the potential of the form with point values, as a picture. The shape of the set of waves is chosen from it.
2. The amplitude in waves whose $n'$ is not $n + pN_\text{fp}$. The exact potential has none. The sums over 128 lines do not have the symmetry of the 5 field periods exactly, because a turn by one field period does not map 128 lines onto each other, and this amplitude shows by how much.
3. The vacuum matrix against the size of the set, with the weight $W_i$ of equation B.91 and with the area weight of the prototype.
4. The error of both forms for the exact potential.
5. The share of the entries for which the rule for intervals that are close applies. This share is computed again at every level, so it limits what stage 2 can save.
6. The times of the entries and of the LU, and the memory with and without the vacuum.
7. On the twin: whether the rows of one field line give the same sums as all rows (decision 5).

*Gates.*

- **0a, twin.** (i) Let $K_\text{c}$ be the smallest $K$ for which two more poloidal mode numbers on each side change the vacuum matrix by $d < 10^{-5}$; pass when $K_\text{c} \le 6$. (ii) At $K_\text{c}$ the distance to the vacuum matrix of the form with point values is $d < 2.5\cdot10^{-3}$, a quarter of the 1 % by which the diagonal elements of the form with point values differ at least from those of the axisymmetric version at this number of points (0.98 to 0.99, regression test of the twin). (iii) For the exact potential, the error with wave amplitudes is at most 1.5 times the error with point values.
- **0b, W7-X.** (i) There is a set of at most 500 waves for which a larger set ($P+1$, $K+2$) changes the vacuum matrix by $d < 10^{-3}$. (ii) For that set the distance to the vacuum matrix of the form with point values is $d < 3\cdot10^{-3}$, a quarter of the 1.2 % by which the vacuum matrices of the form with point values on six grids differ (theory document, section B.8). (iii) As 0a (iii).

The cut cbm18a has no gate; its numbers are reported. The weight of equation B.91 is kept unless the area weight makes the distance of (ii) or the error of (iii) smaller by more than a factor 2. The measured distances of (ii) are called $d_\text{twin}$ and $d_\text{W7-X}$ below. If a gate fails, the plan stops there, the issue gets the numbers, and Toon decides.

### Stage 1. One level without an N × N matrix (steps 1 to 3 of B.12)

*What changes.*

1. `calc_GH_1` (`vac_ops`) gets the column limits as an argument. `calc_GH` passes `vac%lims_c`, so the form with point values runs the same statements as now.
2. A new routine in `vac_ops` does what section 4 describes and returns the complex $M \times M$ matrix $\mathsf{F}(\mathsf{E}^\dagger\mathsf{W}\mathsf{E}')\mathsf{A}^{-1}\mathsf{B}\,\mathsf{F}$. Inside it, one routine computes for a block of rows the product of the rows with the waves and the seven sums. The sums over rows of section B.13 multiply with the waves at extra points of a row, so they would take the place of that one routine. When the input asks for wave amplitudes, `calc_vac_res` calls it in the place of `calc_GH`, `solve_potential` and `project_response`. `finish_response` is divided in two: the conversion of the real $2M \times 2M$ array stays with the form with point values, and the last steps (section 4) serve both forms. The phase of a wave comes from one function for both forms.
3. `vac%G` and `vac%H` are allocated in `calc_GH`, where they are computed, not in `init_vac`. `copy_vac` and `interlaced_vac_copy` copy them only when they exist. The test `fit_diag_pair`, which fills them after `init`, allocates them itself. The estimate for `max_tot_mem` counts G and H in every equilibrium job (`calc_memory_vac`, version 2.58); with `vac_sol_style` 2 it counts the two blocks of rows and the matrices of the size $M' \times M'$ in their place.
4. The inputs of decision 2, read in `input_ops` and broadcast in `MPI_ops`. With `vac_sol_style` 2, PB3D stops when `check_vac` is on, which needs the stored matrices, and when `vac_diag_style` is 2, a check of the form with point values that this plan does not carry over. It also stops with the toroidal flux as normal coordinate (section 6). It stops when a wave of the set has fewer than 2 points per wavelength across the lines ($|n'|\,\Delta\alpha \ge \pi$) or along them ($|n'q - m'|\,\Delta\theta \ge \pi$), and it warns below 10 points per wavelength. For HELENA the input has no effect, and a message says so.
5. PB3D prints the number of waves, and the distance $d$ between the vacuum matrix and the one without the outermost waves ($P-1$ and $K-1$, where these are not negative). The second matrix needs no new sums, because its $\mathsf{A}$ and $\mathsf{B}$ are parts of the full ones. A distance above $10^{-3}$ gives a warning.

*Tests.* New cases in the full-stack suite `vac_3d`, also under 2 and 4 processes: (i) the sums against the same products formed from the matrices of `calc_GH`, on the circular torus of the suite and on its rippled torus with 5 ripple periods and $P = 1$; (ii) `vac%res` of the two forms on the circular torus, to $2\,d_\text{twin}$. One more row in the regression test of the twin: the free-boundary run with `vac_sol_style` 2, whose eigenvalue has to lie within $2.5\cdot10^{-3}\,|\Delta_\text{vac}|$ of the eigenvalue with point values, with the check of `check_vac_res.py` that the vacuum matrix is Hermitian and not negative.

*Gates.*

- **1a.** Every element of $\mathsf{A}$ and $\mathsf{B}$ from the sums equals the element from the products with the stored matrices to $10^{-10}$ of the largest element, on both test edges, on 1, 2 and 4 processes. (The two ways add the same numbers in another order, and the diagonal entries are differences of sums, so agreement to the last digits is not expected.)
- **1b.** On the twin and on W7-X, `vac%res` from PB3D equals the matrix of the stage 0 script for the same waves to $d < 10^{-10}$.
- **1c.** The twin on 128 lines × 200 points, one level, both forms with the code of that day: the two eigenvalues differ by less than $2.5\cdot10^{-3}\,|\Delta_\text{vac}|$. (The value $-0.43849$ of the run on 128 × 400 points cannot be used: it is older than the weights of the last radial points of a free edge (`calc_end_weights`, 17 November), which moved the eigenvalue at 64 × 100 points by $3\cdot10^{-4}$.)
- **1d.** W7-X: the two eigenvalues differ by less than $3\cdot10^{-3}\,|\Delta_\text{vac}|$.
- **1e.** The twin on 128 lines × 400 points with wave amplitudes: the peak memory (`MaxRSS` of the Slurm job) is less than 1 GB above that of the same run with a fixed boundary, and POST of that run without the vacuum plot allocates no $N \times N$ array.
- **1f.** W7-X on one core: `calc_vac_res` with wave amplitudes takes at most the time of the entries plus 0.3 of the time of the LU, both from stage 0. (By operation count the products with the waves are at most a fifth of the LU.)
- **1g.** With the default `vac_sol_style` 1, `vac%res` and the eigenvalue of the twin are the same to the last bit as before the change, and `ctest` passes on Snellius and in CI.

### Stage 2. From one level to the next (step 4 of B.12)

*What changes.* An entry of the code is written as the kernel at the point times the weight (equation B.24), plus a correction that is not zero only where the rule for intervals that are close applies. The sums with the first part are kept from level to level: for $\mathsf{A}$ and $\mathsf{B}$ one quarter of the previous sums plus the pairs with a new point (equation B.95), and for the seven sums of an old point one half plus the new points. The sums with the correction, in $\mathsf{A}$, $\mathsf{B}$ and the seven sums, are computed again at every level, as section B.12 says under "What is not tested". Two different points at the same position, which exist when $q\,n_\alpha$ is an integer, have no finite kernel; their entry goes into the correction as a whole. The kept sums go into `vac_type`, into `interlaced_vac_copy`, and into the group `vac_R_<level>` of the output file, so that a restart at a later level finds them. When they are not there, or belong to another set of waves, the level is computed from the start, which is the code of stage 1.

*Gates.*

- **2a.** At levels 2 and 3 the sums that are carried equal the sums computed from the start to $10^{-10}$ of the largest element, on the two test edges and on W7-X (61, 121 and 241 points per line).
- **2b.** A run of three levels gives the same eigenvalues with and without carrying, to $10^{-10}$ relative, and so does a restart at level 2 and at level 3.
- **2c.** A level that is carried takes at most 0.85 of the time of the same level computed from the start.

*What this stage is worth.* It saves at most one quarter of the sums per level, less the share of point 5 of stage 0, and about one quarter of the total over all levels. It needs about as much new code as stage 1. To find the intervals that are close to an old point, the code has to compute the distance to every other point again, which costs a part of the time that is saved, unless the points are sorted into cells first. The stage helps only runs with more than one level; issue #11 found that the error of the vacuum does not follow the powers that the extrapolation assumes, and the regression deck of the twin runs one level. With the sums over rows of section B.13 the rule for intervals that are close is not needed, so the part of this stage that computes the correction would not be used with them; the rule of one quarter stays. My recommendation is to decide on stage 2 after stage 1, with the measured times.

### Stage 3. POST without an N × N matrix

`vac_pot_plot` computes G and H again and solves for the potential at the points of the edge. For a run that was made with wave amplitudes it would use the amplitudes instead: PB3D writes $\mathsf{A}^{-1}\mathsf{B}$ and the list of waves into the output file, POST sums the waves at the points of the edge, and `eval_pot_ext` computes rows for the points off the edge only. Whether an output file of an older version is refused (`min_PB3D_version`) or read without these data is decided in this stage.

Gate: on the twin, at the plot points that lie more than 10 cm from the edge, the potential of the two forms differs by less than $2\,d_\text{twin}$ of its largest value, and POST's peak memory for the plot at 128 × 400 points is below 4 GB. Closer to the edge nothing is compared: there the sum over the points of the edge does not resolve the kernels (on the cut cbm18a the vacuum field of POST is off by 2.4 to 15 % within 10 cm, comment of 24 November on the issue), and the nearest points decide the value, among them the ends of the lines, where the two forms differ most. Until this stage, the vacuum plot works on runs of either form, with the matrices and the memory of the form with point values. This stage touches the routines that the open work of issue #81 changes, so it comes after that work.

### Stage 4. Documentation and default

The header of `vac_ops`, the page of the inputs on the documentation site (its own repository), and section B.12 with table B.8 of the theory document, which says now that nothing of the section is in the code. Whether `vac_sol_style` 2 becomes the default is decided here, from the results of stages 1 to 3.

## 6. What can go wrong

- **W7-X needs thousands of waves.** The products with the waves then cost as much as the LU, and only the memory is gained. Gate 0b (i) decides this before any Fortran is written.
- **The waves with $p \ne 0$ change fast along a field line.** Their number $n'q - m'$ is about $pN_\text{fp}q$, which is 5.7 per $p$ on W7-X: with 121 points per line a wave with $p = 2$ has 11 points per wavelength. On the circular torus the theory document gives a diagonal element of the vacuum matrix that is 1.9 % too small with 9 points per wavelength (section B.8). The form with point values holds the same waves on the same points, so this is not a difference between the two forms, but the warning of point 4 of stage 1 will show it.
- **The entries take most of the time now.** The form with wave amplitudes then gains memory and little time. Point 6 of stage 0 measures this.
- **The two forms differ by more than on the prototype.** In the form with point values the potential at the two end points of every line is off by 3.5 to 7 % of its largest value on the cut cbm18a (comment of 24 November on the issue). A sum of waves has one value at each point. Gates 0a and 0b limit the difference, and the exact potential says which of the two forms has the smaller error.
- **The toroidal flux as normal coordinate.** Since version 2.56 a free-boundary run can use it (issue #22). The phase of a wave is then $m'\alpha_j + (n' - \iota m')\zeta_j$, with $\iota$ the rotational transform at the edge and $\zeta_j$ the parallel coordinate, which is the same function of position as for the poloidal flux. Section 4 is written for the poloidal flux, and the three cases of stage 0 use it. The set of waves is a set of pairs $(n', m')$ for both fluxes, and the phase of a wave comes from one function for both forms, so the toroidal flux needs no sums of its own; but stage 0 measures nothing with it. Stage 1 therefore stops for `vac_sol_style` 2 with the toroidal flux. To remove that stop, the two forms have to be compared once on the regression case `circ_eps033_flux_labels`.
- **Results on different numbers of processes** differ by rounding ($10^{-13}$, estimate), because the sums are added in another order. The form with point values has the same property through the LU.
- **No three-dimensional equilibrium is in CI** (issue #20). The rippled torus tests the sums with $P > 0$; W7-X itself is run by hand on Snellius.

## 7. Not in this plan

- The sums over rows of section B.13 of the theory document, which give an error in powers of the step along the field lines. They need wave amplitudes as unknowns, so this plan comes before them, and point 2 of stage 1 names the routine that they would replace.
- The entries themselves: the number $n_\text{near}$ (issue #9), the edge of the zone of the rule (issue #8), the direct diagonal (issue #14).
- The coupling of $n$ and $n \pm N_\text{fp}$ in the plasma (issue #21). The waves with other toroidal mode numbers are unknowns of the vacuum only.
- `check_vac` for the form with wave amplitudes.
