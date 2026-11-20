---
type: "Concept"
title: "Vacuum Response (Boundary Element Method)"
description: "How PB3D computes the free-boundary vacuum energy: the field-line 3-D style for VMEC and the axisymmetric style for HELENA, the sign convention along the stored normal -J grad psi and the Green identities the discrete operators satisfy, the diagonals of H and G (vac_diag_style), the exterior solve, the response matrix added to the edge block of A, the check_vac diagnostics, POST's vacuum potential plot, the tests, and the known limits and cost."
tags: [vacuum, bem, free-boundary, greens-function, field-line, axisymmetric, check_vac]
verified:
  - by: openwiki/0.9.0
    at: 2026-11-20T18:32:33.626Z
sources:
  - id: openwiki-source-4361ace621d24e5f0bcbb1a2
    resource: repo://Modules/driver_sol.f90
  - id: openwiki-source-99e71edd5b71e76c5f03dd84
    resource: repo://Modules/input_ops.f90
  - id: openwiki-source-20a51f9fa5591dbcc1483985
    resource: repo://Modules/SLEPC_ops.f90
  - id: openwiki-source-2cf1be2deef2573fd19af4cf
    resource: repo://Modules/vac_ops.f90
  - id: openwiki-source-259f5b46b16aa5894e6436c2
    resource: repo://tests/fullstack/test_vac_3d.f90
generated: { by: "claude-code", at: "2026-11-20T18:32:33.626Z" }
---

# Vacuum Response (Boundary Element Method)

With a free boundary (`BC_style(2) = 2`), the perturbed magnetic field outside the plasma is the gradient of a potential φ that is harmonic in the vacuum and decays at infinity. PB3D computes φ on the plasma boundary with a boundary element method and turns it into a response matrix: the vacuum energy as a quadratic form in the edge displacement of the modes. That matrix is added to the edge block of the potential-energy matrix A. With a fixed boundary (`BC_style(2) = 1`) the vacuum is skipped: `run_driver_sol` does not call the vacuum operations.

All of it is in `vac_ops.f90`, with the kernel integrals and distributed-matrix helpers in `vac_utilities.f90`, the type `vac_type` in `vac_vars.f90`, and the toroidal functions in `dtorh.f90` (see [Module Map](../architecture/module-map.md)).

## Two styles

The equilibrium code fixes the style (`vac%style`):

| | Style 1: field-line 3-D | Style 2: axisymmetric |
|---|---|---|
| Used for | VMEC equilibria (`eq_style = 1`) | HELENA equilibria (`eq_style = 2`) |
| Boundary points | `n_par_X` points on each of `n_alpha` field lines, stored line by line: point l of line k is column (k − 1) n_p + l of G and H | the poloidal points of the HELENA boundary, with their angles in `vac%ang` |
| Surface integral | trapezoidal rule along each field line, equal weights Δα across the lines, which together cover the surface once | collocation in the poloidal angle; the toroidal integral of the Green's function is analytical, through the toroidal functions P and Q of half-integer degree (`dtorh.f90`) |
| Near the influence point | within four interval lengths, Gauss–Legendre points on the cubic through four points of the same line (`calc_GH_int_1`); at the point itself, the analytical integral over the half-cell, Δα/2 to either side of the line | analytical near-singular integrals and a corner term (`calc_GH_int_2`, `calc_GH_2`) |
| Input requirement | `alpha_style = 2`, so that there are at least two field lines; the input check refuses a free-boundary VMEC run with `alpha_style = 1` | none |

In style 1 the stored normal N = −J∇ψ has the length of the surface element per unit parallel coordinate and per unit field-line label, so the weights need no other Jacobian. The neighbouring field lines of a stellarator can lie much closer together than the points along them; this is why the near-field rule integrates along the line instead of sampling the kernel at the points.

## Sign convention and Green's identities

Both styles take the normal derivative along the stored normal −J∇ψ, which points out of the vacuum, in the kernels, in the corner term and in the data. With the diagonals that `calc_GH` sets, the discrete operators satisfy, on the boundary,

- H φ − G dφ = 0 for a potential that is harmonic in the vacuum and decays at infinity,
- H φ − G dφ = −4π φ for a potential that is harmonic inside the boundary,

where dφ is the normal derivative along −J∇ψ. At a point off the boundary, the same sum without a diagonal gives 4π φ(x) in the vacuum and 0 in the plasma for a vacuum potential (`eval_pot_ext`). The full-stack tests check these identities.

## Diagonals of H and G

The input `vac_diag_style` selects how the diagonal entries are set:

- **1, indirect (default):** from the identities above, for potentials whose boundary values and normal data are known exactly (`diag_indirect`, with the least-squares fit `fit_diag`).
  - Field-line style, in two steps. First H_ii from the constant potential, which has zero normal data; this gives the row-sum rule Σ_j H_ij = −4π. Then G_ii from the potentials x, y and z, with that H_ii, as a least-squares solution of the three equations.
  - Axisymmetric style, in one step: H_ii and G_ii together from R^|n|, R^|n| Z and a ring of sources, each potential scaled to unit maximum. These diagonals converge at third order in the poloidal points; the direct ones converge at first order.
- **2, direct:** a diagnostic, not for production runs; the input check prints a warning for a free-boundary run. Field-line style: H_ii is −2π plus the integral of the kernel over the node's cell (`calc_diag_direct`) where the node has a full cell, and the indirect value at the ends of the lines; G_ii keeps the half-cell integral. Axisymmetric style: the values of the assembly, with the corner term and the analytical near-singular integrals.

<!-- openwiki: broken internal link [/Documentation/w7x_check_vac_diag.md] link "/Documentation/w7x_check_vac_diag.md" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
The difference between the assembled diagonal and the other choice measures the quadrature error of the off-diagonal part of the row; `check_vac` prints it. On W7-X the choice matters far more than on circular cases ([`/Documentation/w7x_check_vac_diag.md`](/Documentation/w7x_check_vac_diag.md)).

## Exterior solve and response matrix

`calc_vac_res` computes the response in four steps:

1. **Data.** The normal data of each mode along the stored normal: −∂X/∂(parallel coordinate) for the field-line style, with the toroidal factor split off for the axisymmetric style. The modes are the secondary mode numbers at the plasma edge.
2. **Solve.** `solve_Phi_BEM` solves H Φ = G R, with R the data, for all modes at once, with the same equation for both styles and no diagonal shift. It uses STRUMPACK-Dense when PB3D was built with it, checks the residual against the uncompressed operator and falls back to ScaLAPACK LU above a relative residual of 1e-6; without STRUMPACK it uses ScaLAPACK LU (`pdgesv`).
3. **Projection.** Φ is projected back on the data with the integration rule: trapezoidal weights along the lines, and for the field-line style the average over the lines.
4. **Finish.** The result is divided by the permeability, multiplied by 2π (the toroidal integral), and replaced by its Hermitian part, which removes a discretisation asymmetry of about 1e-5 that would make SLEPc reject tightly converged eigenvalues as complex.

With this convention the response is the positive Neumann-to-Dirichlet map, and the vacuum adds 2 δW_vac = 2π X_e^H res X_e to A, where X_e is the vector of edge displacements of the modes. `set_BC` inserts `vac%res` in the edge diagonal block of A. On Richardson levels above the first, the matrices of the previous level are copied into the new layout (`interlaced_vac_copy`) and completed; with `jump_to_sol` the vacuum is restored instead of recomputed.

## Diagnostics: `check_vac`

With the input `check_vac`, `check_vac_exact` tests the assembled operators on the real equilibrium surface after `calc_GH`:

- for potentials harmonic inside the boundary (field-line style: 1, x, y, z, x² − y²; axisymmetric style: the potentials of `set_pot_2`) and potentials harmonic outside (a point source on the axis, and rings of sources along the axis with unit weights and with weights cos(nζ) for the primary mode number), it prints the residual of H φ − G dφ − c φ for c = −4π, −2π, 0, 2π, 4π;
- it solves the exterior problem for the exterior potentials and compares the boundary trace with the exact one;
- it compares the assembled diagonals with the alternative ones;
- it evaluates the exact exterior potentials at points off the boundary (`check_vac_ext`), with and without the correction of `sub_nearest`.

Every number is on a line that starts with `CHECKVAC`. `scripts/snellius/checkvac_table.py` and `diag_table.py` collect them from the log.

## POST vacuum potential plot

With `plot_vac_pot`, POST's `vac_pot_plot` plots the potential of every eigenmode outside the plasma, for both styles. The normal data are those of `calc_vac_res` (dφ = −i(nq − m)X, which is −J Q_n, with Q_n the normal magnetic perturbation that POST plots). The boundary potential comes from the exterior solve, and the potential at a point off the boundary from Green's representation (`eval_pot_ext`). The points lie in the planes of constant ζ_F of the plasma plots, on a grid of `n_vac_plot` points in R and Z between `min_Rvac_plot`…`max_Rvac_plot` and `min_Zvac_plot`…`max_Zvac_plot` in metres; by default the box is the boundary's extent widened on both sides by half its larger width. The input check switches the plot off for a fixed boundary.

Outputs per mode: `[id]_vac_pot_RE` and `[id]_vac_pot_IM` (HDF5) and `Data/[id]_vac_bnd.dat`; once: `vac_in_plasma`, 1 in the plasma and 0 in the vacuum, where the potential is set to zero.

Within a few spacings of the boundary points the quadrature does not resolve the kernels. `sub_nearest` subtracts the value at the nearest boundary point through an interior harmonic, which removes most of that error. Two cases stay inaccurate in the axisymmetric style, as the comment in `calc_GH_2` explains: a point that coincides with a boundary point, for a potential harmonic inside the boundary; and a point closer to the boundary than about one point spacing but farther than the tolerance from every boundary point.

## Tests

- **Full-stack** (`tests/fullstack/test_vac_kernels.f90`, `test_vac_greens.f90`, `test_vac_3d.f90`): the kernels and their singular and edge corrections (`G_regular`, `H_regular_numderiv`, `GH_singular_quadrature`, `GH_edge_correction`, `GH_int_1_regular`, `GH_int_1_singular`, `GH_int_1_near`, `H_sing_node`, `solid_angle_node`); the Green identities in both styles (`greens_identity`, `greens_identity_convergence`, `greens_identity_indirect`, `greens_identity_3d`); exterior solves against exact solutions (`solve_Phi_BEM_roundtrip`, `point_source_3d`, `point_source_direct_3d`, `exact_exterior_2d`, `exact_exterior_2d_options`); the diagonals (`diag_direct_3d`); the response against the large-aspect-ratio cylinder (`vac_response_cylinder`); and the field-line style against the axisymmetric one on the same boundary (`response_1_vs_2`, with five recorded ratios checked to 0.5 %).
- **Regression** (`tests/regression/`): `cbm18a_free_bnd` (with the response-matrix check `check_vac_res.py`), `circ_eps005_free_bnd` (the circular external kink, also at 400 radial points), and `circ_eps033_vmec_twin` (the VMEC twin of a HELENA case, free through the 3-D vacuum and compared with the axisymmetric vacuum of its twin).

See [Testing & Extension](../testing-and-extension.md) for how to run them.

## Accuracy, cost and known limits

- The field-line style converges at first order in the number of field lines. On the circular case its external kink agrees with TERPSICHORE to 2e-5 at 128 lines × 400 parallel points; at 64 lines × 100 points expect about 1 %. The working rule is at least 128 field lines and about 36 parallel points per poloidal wavelength.
- That accuracy is expensive: 128 × 400 is 51 000 boundary points, which took 220 GB and 39 minutes on eight cores. Which part of the method sets this cost is not settled.
- The axisymmetric style with the indirect diagonal converges towards a limit 0.13 % from TERPSICHORE and from the 3-D vacuum on the circular case; the cause is not known.
- On W7-X the response converges only when lines and points are refined together; 61 points per line are too few.

<!-- openwiki: broken internal link [/Documentation/vacuum_work_overview.md] link "/Documentation/vacuum_work_overview.md" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/Documentation/vacuum_work_report.md] link "/Documentation/vacuum_work_report.md" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/Documentation/plan_vacuum_bugfix.md] link "/Documentation/plan_vacuum_bugfix.md" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
The numbers come from [`/Documentation/vacuum_work_overview.md`](/Documentation/vacuum_work_overview.md), which also lists the comparisons with TERPSICHORE and MISHKA-1. The detailed report is [`/Documentation/vacuum_work_report.md`](/Documentation/vacuum_work_report.md), and the plan with the night-by-night results is [`/Documentation/plan_vacuum_bugfix.md`](/Documentation/plan_vacuum_bugfix.md). Extend the tests when you change the vacuum.
