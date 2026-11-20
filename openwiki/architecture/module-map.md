---
type: "Reference"
title: "Module Map (`Modules/`)"
description: "What each of the 50 Fortran modules in Modules/ does, grouped by layer in the order of PB3D_MODULE_SOURCES in CMakeLists.txt, which is the Fortran module dependency order: foundations, messages and numerics, MPI, HDF5, grids, equilibrium, VMEC and HELENA, perturbation, solution, vacuum, Richardson, SLEPc, input/output, PB3D output, drivers and the debug test module."
tags: [architecture, modules, fortran, layering]
verified:
  - by: openwiki/0.9.0
    at: 2026-11-20T18:32:33.626Z
sources:
  - id: openwiki-source-d44494ef3e497fea81240ef8
    resource: repo://CMakeLists.txt
  - id: openwiki-source-d274e5b11d2119c5e0cd6d0d
    resource: repo://Modules/dtorh.f90
  - id: openwiki-source-3877069554ef230e9ffc979c
    resource: repo://Modules/flat_vars.f90
  - id: openwiki-source-2afe37558c38bb26050c81de
    resource: repo://Modules/messages.f90
  - id: openwiki-source-2b951757a1f7314b5a33d0bf
    resource: repo://Modules/MPI_ops.f90
  - id: openwiki-source-ee58f84f32e5ba524ba9d291
    resource: repo://Modules/num_utilities.f90
  - id: openwiki-source-0e5296407709fa9254158dd8
    resource: repo://Modules/num_vars.f90
  - id: openwiki-source-20a51f9fa5591dbcc1483985
    resource: repo://Modules/SLEPC_ops.f90
  - id: openwiki-source-7313dfbfb11179d6dca8938e
    resource: repo://Modules/sol_utilities.f90
  - id: openwiki-source-fcf968e61512cc7d1e7c3774
    resource: repo://Modules/spline_utilities.f90
  - id: openwiki-source-2cf1be2deef2573fd19af4cf
    resource: repo://Modules/vac_ops.f90
  - id: openwiki-source-08dfb6eb1e1e603945fe5a99
    resource: repo://Modules/vac_vars.f90
  - id: openwiki-source-059d81dada4bf2d9219c1c03
    resource: repo://Modules/VMEC_ops.f90
  - id: openwiki-source-59322ed006cffaf717f920e4
    resource: repo://Modules/X_ops.f90
generated: { by: "claude-code", at: "2026-11-20T18:32:33.626Z" }
---

# Module Map (`Modules/`)

All physics, numerics, input/output and driver code lives in the 50 files `Modules/*.f90`. CMake compiles them once into the object library `pb3d_modules`, which both executables link (see [Build & Dependencies](build-and-dependencies.md)). The list `PB3D_MODULE_SOURCES` in `CMakeLists.txt` is in Fortran module dependency order: a module comes before every module that `use`s it. That order is also the clearest map of the layers, and this page follows it. A module named `*_vars` holds types and variables, `*_utilities` holds helpers, and `*_ops` holds the operations that the drivers call.

## 1. Foundation

- `num_vars.f90`: the program-wide variables. It holds the numeric kinds (`dp`), constants (`pi`), the program version (`prog_version`) and the oldest PB3D output version that POST and a Richardson restart accept (`min_PB3D_version`), MPI rank and size, tolerances, memory limits, and the runtime switches read from the input file (`eq_style`, `X_style`, `BC_style`, `EV_style`, `alpha_style`, `norm_style`, `vac_diag_style`, `check_vac`, the plot switches, ...). Read this file first when you need to know what a runtime option controls.
- `str_utilities.f90`: conversions between numbers and strings (`i2str`, `r2str`, ...).

## 2. Messages and numerics

- `messages.f90`: output to the terminal: `writo` (at most 120 characters per line, including the time stamp and indentation), `lvl_ud` for indentation, `print_err_msg`, `print_hello`/`print_goodbye`, and timing with `start_time`/`passed_time`.
- `num_utilities.f90`: numerics without external libraries apart from LAPACK: finite-difference weights, the end weights of the free edge (`calc_end_weights`), integration rules, index helpers, zero finders and small dense determinants and inverses.
- `spline_utilities.f90`: splines for real and complex data (`spline_real`, `spline_complex`), wrapping the PSPLINE/EZspline library. It is a module of its own so that `num_utilities` compiles and is unit-tested without PSPLINE.
- `num_ops.f90`: numerical operations that need more than `num_utilities`; under `ldebug` it uses `output_ops` (HDF5).
- `files_utilities.f90`: low-level file helpers.

## 3. MPI

- `MPI_vars.f90`, `MPI_utilities.f90`, `MPI_ops.f90`: `start_MPI` and `stop_MPI`; `broadcast_input_opts`, which sends the input options to all ranks; the generic `broadcast_var`; `redistribute_var`, which takes every point of an overlap from the rank that owns it; and the error path `sudden_stop` → `abort_MPI`, which ends the run with an exit status that is not zero.

## 4. HDF5 output

- `flat_vars.f90`: the flattened variable type used for HDF5 storage. It has no HDF5 dependency, so that `PB3D_utilities` and other modules whose only HDF5 link is this type compile and are unit-tested without HDF5.
- `HDF5_vars.f90`, `HDF5_utilities.f90`, `HDF5_ops.f90`: all output files, written as HDF5 with XDMF metadata so that they open in tools such as ParaView.

## 5. Grids and equilibrium types

- `grid_vars.f90`: the grid type and grid sizes (`n_r_eq`, `n_r_X`, `n_r_sol`, the number of field lines `n_alpha`, the field-line label `alpha`).
- `eq_vars.f90`: the equilibrium types (`eq_1_type` for flux quantities, `eq_2_type` for metric quantities) and the normalisation constants (`R_0`, `B_0`, `psi_0`, `pres_0`, `rho_0`, `T_0`).
- `grid_utilities.f90`, `grid_ops.f90`: construction of the grids, the field-aligned grid, redistribution over the ranks, and removal of ghost regions (`trim_grid`).

## 6. Equilibrium

- `eq_utilities.f90`, `eq_ops.f90`: `eq_ops.f90` is the largest module (about 7700 lines). It computes the derived equilibrium quantities (`calc_derived_q`), the normalisation, flux quantities and their plots (`flux_q_plot`, `B_plot`, `J_plot`, `kappa_plot`), and divides the equilibrium work into jobs that fit in memory (`divide_eq_jobs`, `calc_eq_jobs_lims`).

## 7. Equilibrium codes

- `VMEC_vars.f90`, `VMEC_utilities.f90`, `VMEC_ops.f90`: `read_VMEC` reads a VMEC `wout` file through LIBSTELL's `read_wout_mod`, keeps the Jacobian, |B| and the covariant field components on VMEC's own Nyquist mode list (matched by (m, n)), and interpolates the half-mesh quantities from VMEC's half-mesh positions. `calc_trigon_factors` stores the one-angle cosine and sine tables, and the generic `fourier2real` forms real-space values for any mode list.
- `HELENA_vars.f90`, `HELENA_ops.f90`: `read_HEL` reads the HELENA mapping file, and `interp_HEL_on_grid` interpolates it onto PB3D's grids.
- The input variable `eq_style` (1 = VMEC, 2 = HELENA) selects between them at run time.

## 8. Perturbation

- `X_vars.f90`, `X_utilities.f90`, `X_ops.f90`: the mode numbers per flux surface (`init_modes`, `setup_modes`, which puts each mode in the slot computed from its number and checks the two mode tables against each other), the vectorial and tensorial perturbation quantities (the generic `calc_X`), and the magnetic integrals (`calc_magn_ints`) that become the matrices of the eigenvalue problem. Mode setup has been the source of several past bugs; treat changes here as risky and run the regression layer.

## 9. Solution

- `sol_vars.f90`, `sol_utilities.f90`, `sol_ops.f90`: post-processing of the eigenvectors: the energy decomposition (`decompose_energy`), the plots of eigenvalues and eigenvectors (`plot_sol_vals`, `plot_sol_vec`, `plot_harmonics`), and the spline interpolation of a solution vector in the radial direction (`interp_V_spline` in `sol_utilities`).

## 10. Vacuum

The free-boundary vacuum response is computed with a boundary element method; [Vacuum Response](../domain/vacuum-response.md) explains it.

- `vac_vars.f90`: the vacuum type `vac_type`. Its `style` is 1 (field-line 3-D, used for VMEC equilibria) or 2 (axisymmetric, used for HELENA equilibria).
- `vac_utilities.f90`: the kernel integrals of both styles (`calc_GH_int_1`, `calc_GH_int_2`), the helpers for distributed ScaLAPACK vectors and matrices, the direct diagonal of the field-line style (`calc_diag_direct`, with the solid angle `solid_angle_node`), the least-squares diagonal fit (`fit_diag`), and the copy of the vacuum to the next Richardson level (`interlaced_vac_copy`).
- `vac_ops.f90`: storing the boundary geometry (`store_vac`, `store_vac_VMEC`, `store_vac_HEL`), assembling the matrices G and H (`calc_GH`, `calc_GH_1`, `calc_GH_2`) with their diagonals (`diag_indirect`), the response matrix that enters the energy (`calc_vac_res`), the exterior solve (`solve_Phi_BEM`), the diagnostics behind the input `check_vac` (`check_vac_exact`, `check_vac_ext`), the potential away from the boundary (`eval_pot_ext`), POST's vacuum potential plot (`vac_pot_plot`) and the vacuum output (`print_output_vac`).
- `dtorh.f90`: the toroidal functions P and Q of half-integer degree, adapted from the DTORH1 routine of Segura and Gil, used by the axisymmetric kernels.

## 11. Richardson extrapolation

- `rich_vars.f90`, `rich_ops.f90`: the outer loop over Richardson levels of increasing resolution (`init_rich`, `do_rich`, `start_rich_lvl`/`stop_rich_lvl`, `calc_rich_ex`, `term_rich`); see [Simulation Pipeline](../workflows/simulation-pipeline.md).

## 12. SLEPc interface

- `SLEPC_utilities.f90`, `SLEPC_ops.f90`: assembly of the generalised eigenvalue problem (`setup_mats`, with `fill_mat` summing the energy over the radial points with their weights), the boundary conditions (`set_BC`), the solver set-up (`setup_solver`), and the results (`get_solution`, `summarize_solution`, `store_results`). The routines here need a trimmed solution grid.

## 13. Input and output

- `input_utilities.f90`, `input_ops.f90`: reading and checking the namelists `inputdata_PB3D` and `inputdata_POST`; see [Input & Configuration](../workflows/input-and-configuration.md).
- `output_ops.f90`: output helpers built on `HDF5_ops`.
- `files_ops.f90`: the command line (`parse_args`, with the list of `--` options per program), and opening the input and output files (`open_input`, `open_output`).

## 14. PB3D output

- `PB3D_utilities.f90`, `PB3D_ops.f90`: reading PB3D output back from HDF5 (`reconstruct_PB3D_in`, `reconstruct_PB3D_grid`, `reconstruct_PB3D_eq_1`/`_2`, `reconstruct_PB3D_sol`, ...). PB3D uses it for Richardson restarts and POST to read PB3D's results.

## 15. Drivers

- `driver_eq.f90` → `run_driver_eq`: the equilibrium grid and quantities for one Richardson level or equilibrium job.
- `driver_X.f90` → `run_driver_X`: the perturbation grid and quantities.
- `driver_sol.f90` → `run_driver_sol`: the solution grid, the eigenvalue problem through `SLEPC_ops`, and Richardson extrapolation.
- `driver_POST.f90` → `init_POST`, `run_driver_POST`, `stop_POST`: all of POST. It reconstructs the PB3D state, builds the output grids, and writes the plots and tables, including the vacuum potential plot.

## 16. Debug test module

- `test.f90` → `generic_tests`: compiled only with `ldebug` and run only with the command-line flag `--test` (`-t`). It is an interactive developer tool, not an automated test; see [Testing & Extension](../testing-and-extension.md).

## Vendored libraries (`Libraries/`, not `Modules/`)

- `dfft.f`: a double-precision FFT package, built as the static library `dfftpack`.
- `foul.f90`: a formatted terminal output library, built as the static library `foul`.
