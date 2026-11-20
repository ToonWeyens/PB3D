---
type: "Concept"
title: "Physics Domains"
description: "The physics behind PB3D's module groups: the imported VMEC or HELENA equilibrium and its normalisation, the perturbation modes, the radial discretisation with left finite differences and the free edge, the boundary conditions, the generalised eigenvalue problem for omega squared (negative is unstable), Richardson extrapolation, grids and storage, and POST; the vacuum response has its own page."
tags: [physics, mhd, equilibrium, vmec, helena, eigenvalue, richardson, discretisation]
verified:
  - by: openwiki/0.9.0
    at: 2026-11-20T18:32:33.626Z
sources:
  - id: openwiki-source-5a1758b4a37a21c66fce921a
    resource: repo://Modules/driver_POST.f90
  - id: openwiki-source-6ed31598a06aa8e34c662ac6
    resource: repo://Modules/eq_ops.f90
  - id: openwiki-source-99e71edd5b71e76c5f03dd84
    resource: repo://Modules/input_ops.f90
  - id: openwiki-source-ee58f84f32e5ba524ba9d291
    resource: repo://Modules/num_utilities.f90
  - id: openwiki-source-0e5296407709fa9254158dd8
    resource: repo://Modules/num_vars.f90
  - id: openwiki-source-20a51f9fa5591dbcc1483985
    resource: repo://Modules/SLEPC_ops.f90
  - id: openwiki-source-bbe7541e4ee07a1320233e76
    resource: repo://Modules/sol_ops.f90
  - id: openwiki-source-059d81dada4bf2d9219c1c03
    resource: repo://Modules/VMEC_ops.f90
  - id: openwiki-source-59322ed006cffaf717f920e4
    resource: repo://Modules/X_ops.f90
  - id: openwiki-source-23775c3de52f3ab95a13cb8b
    resource: repo://README.md
generated: { by: "claude-code", at: "2026-11-20T18:32:33.626Z" }
---

# Physics Domains

This page explains the physics behind the module groups in [Module Map](../architecture/module-map.md): why the code is structured as it is, not only which files exist.

## The problem PB3D solves

PB3D computes the linear ideal-MHD stability of peeling-ballooning modes. Given a toroidal plasma equilibrium, it finds whether small perturbations grow or oscillate, by discretising the linearised MHD energy along field lines and solving the resulting generalised eigenvalue problem. It extends the usual axisymmetric (tokamak) analysis to 3-D equilibria: tokamaks with ripple and stellarators. The header of `PB3D.f90` cites the theory papers (`weyens2014theory`, `Weyens2017PB3D`).

## Equilibrium (`eq_*`, `VMEC_*`, `HELENA_*`)

PB3D does not compute an equilibrium; it reads one, chosen with `eq_style`:

- **VMEC** (`eq_style = 1`), the 3-D equilibrium code, through LIBSTELL's `read_wout_mod` in `VMEC_ops.f90`. VMEC stores the Jacobian, |B| and the covariant field components on its own Nyquist mode list, which differs from the main mode list when `ntor` > 0; `read_VMEC` keeps that list (`mn_V_nyq`) and matches the modes by (m, n). VMEC stores λ, the Jacobian and the covariant components on its half mesh, at s = (j − 3/2)/(ns − 1), and `read_VMEC` splines them from those positions to the full mesh.
- **HELENA** (`eq_style = 2`), the axisymmetric equilibrium code, through `read_HEL` in `HELENA_ops.f90` from its mapping file.

Both are converted to PB3D's types (`eq_1_type` for flux quantities, `eq_2_type` for metric quantities, in `eq_vars.f90`) and made dimensionless with the constants `R_0`, `B_0`, `psi_0`, `pres_0`, `rho_0`, `T_0` (`calc_normalization_const` in `eq_ops.f90`, selected by `norm_style`). For VMEC with the MISHKA normalisation, the default `B_0` is the (0, 0) Fourier component of |B| on the axis, which is the average over both angles; it does not depend on the choice of toroidal origin. Values given in the input file override the computed constants.

`eq_ops.f90` then derives everything the stability calculation needs: metric coefficients, magnetic shear, curvature, parallel current, safety factor or rotational transform (`calc_derived_q`). Several of these derivations had sign or factor errors in the past (see the README changelog); check the changelog for the quantity before you change `calc_derived_q`.

## Perturbation (`X_*`)

The perturbation is expanded in Fourier harmonics in field-aligned coordinates. There is one primary mode number `prim_X` (the toroidal number n with the poloidal flux as radial coordinate, or the poloidal number m with the toroidal flux) and a set of secondary numbers per flux surface:

- `X_style = 1` prescribes the secondary range `min_sec_X` … `max_sec_X` on all surfaces;
- `X_style = 2` takes a window of `n_mod_X` secondary numbers around the resonant value on each surface (`setup_modes` in `X_ops.f90`).

`setup_modes` puts each mode in the slot computed from its number and checks the two mode tables against each other on every set-up. `calc_X` computes the vectorial and tensorial perturbation quantities, and `calc_magn_ints` integrates them along the field lines into the matrix elements. The window of `X_style = 2` drops the coupling between a mode that leaves the window and one that enters it, so its trial space is not nested; doubling `n_mod_X` once checks the effect. Mode set-up has caused several past bugs; run the regression layer after any change here.

## Radial discretisation and boundary conditions (`SLEPC_*`)

The radial derivative of the perturbation is a left finite difference: the point itself and the `norm_disc_prec_sol` points to its left (an order from 1 to 3). `fill_mat` sums the energy over the radial points with a weight per point. At a free edge with order 2 or 3, the last points carry the end weights of `calc_end_weights` (in `num_utilities.f90`), for which the sum over the points of the finite difference of X equals X at the edge, as the integral of X′ does. With these weights the free-boundary eigenvalue converges at the order of the finite difference.

`set_BC` imposes the boundary conditions, chosen with `BC_style(1:2)`:

- the inner end always takes style 1: the perturbation is set to zero, by giving those rows an artificial eigenvalue `EV_BC`;
- the edge takes style 1 (fixed boundary, the same treatment) or style 2 (free boundary). With style 2 the natural boundary condition is not written into the matrix: the edge row is assembled like any other row and the vacuum term is added to the edge diagonal block, which keeps the matrices A and B Hermitian.

The input check refuses any other value. [Vacuum Response](vacuum-response.md) describes how the vacuum term is computed.

## Eigenvalue problem (`sol_*`, `SLEPC_*`)

The assembled matrices, with the vacuum term and the boundary conditions, define the generalised eigenvalue problem A x = λ B x, solved with SLEPc (`SLEPC_ops.f90`). Here A comes from the potential energy and B from the kinetic energy, and the eigenvalue λ is the normalised square of the frequency, ω². `setup_mats` builds the matrices, `set_BC` applies the boundary conditions, `setup_solver` and `setup_guess` configure the solver (`EV_style`, with an optional first guess `EV_guess`), and `get_solution`, `summarize_solution` and `store_results` extract and store the results.

**A negative ω² means an unstable mode:** the perturbation grows. POST counts the eigenvalues with a negative real part as the unstable ones. `sol_ops.f90` decomposes the energy of each eigenvector into its potential and kinetic parts (`decompose_energy`) and plots the eigenvalues, eigenvectors and harmonics.

## Vacuum response (`vac_*`)

With a free boundary (`BC_style(2) = 2`), the perturbed field outside the plasma must match a vacuum solution. PB3D computes the vacuum with a boundary element method: a field-line 3-D style for VMEC equilibria and an axisymmetric style for HELENA equilibria. The result is a response matrix that adds the vacuum energy to the edge block of A. [Vacuum Response](vacuum-response.md) covers the formulation, the input options, the diagnostics and the tests.

## Richardson extrapolation (`rich_*`)

The perturbation and solution grids have a finite resolution, so the eigenvalues approach the continuum value only as the resolution goes to infinity. `rich_ops.f90` runs PB3D on successively finer grids, up to `max_it_rich` levels, and extrapolates the eigenvalue between levels (`calc_rich_ex`). It stops early when the change is below `tol_rich` (`do_rich`). This is the outer loop of the program; see [Simulation Pipeline](../workflows/simulation-pipeline.md).

## Grids and storage (`grid_*`, `HDF5_*`)

A run uses several related grids: the equilibrium grid, the perturbation grid (often field-aligned) and the solution grid. Each may be distributed over the MPI ranks, and the ghost regions are removed before the final use (`grid_utilities.f90`, `grid_ops.f90`). Every overlap point between two ranks takes its value from the rank that owns it (`redistribute_var`). All results are written as parallel HDF5 with XDMF metadata (`HDF5_ops.f90`), so that they open in tools such as ParaView without a custom reader.

## Post-processing (`driver_POST`)

POST does not recompute the physics. It reconstructs PB3D's state from the HDF5 output (`reconstruct_PB3D_*` in `PB3D_ops.f90`) and produces real-space plots of the field, current and curvature, resonant surfaces, the energy decomposition of every solution (`POST_out_EN_R<level>.txt`), the vacuum potential outside the plasma, and, for HELENA cases, comparison tables against HELENA's own output. The `plot_*` switches of its namelist select the outputs.
