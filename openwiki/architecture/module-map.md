# Module Map (`Modules/`)

All physics, numerics, I/O, and driver logic lives in `Modules/*.f90` and is compiled into a
single object library shared by both executables (see
[Build & Dependencies](build-and-dependencies.md)). The list in `CMakeLists.txt`
(`PB3D_MODULE_SOURCES`) is ordered by **Fortran module dependency** (a module must appear before
anything that `use`s it) — that ordering is also the clearest map of the codebase's internal
layering, reproduced and annotated below.

## 1. Base / foundation (no PB3D dependencies)
- `num_vars.f90` — the central "god module" of runtime configuration: numeric kinds (`dp`,
  `dpi`), physical constants (`pi`, `mu_0_original`), and **every** runtime "style" switch
  (`eq_style`, `EV_style`, `X_style`, `BC_style`, `norm_style`, plotting flags, tolerances, MPI
  rank/size, memory-budget variables for job splitting). Read this file first when you need to
  know what a runtime option controls.
- `str_utilities.f90` — string/number conversion helpers (e.g. `r2str`, `i2str`).

## 2. Messages & general utilities
- `messages.f90` — leveled logging (`writo`, `lvl_ud` for indentation, `print_err_msg`,
  `print_hello`/`print_goodbye`, timing via `start_time`/`passed_time`).
- `num_utilities.f90`, `num_ops.f90` — numerical helpers (derivatives, splines via PSPLINE/EZspline,
  interpolation — see `interp_V_spline`, fixed in v2.46 for the 2-point case).
- `files_utilities.f90` — low-level file helpers.

## 3. MPI
- `MPI_vars.f90`, `MPI_utilities.f90`, `MPI_ops.f90` — `start_MPI`/`stop_MPI`, broadcasting input
  options to all ranks (`broadcast_input_opts`), and `sudden_stop` (the error-abort path invoked
  by the `CHCKERR` macro).

## 4. HDF5 / output storage
- `HDF5_vars.f90`, `HDF5_utilities.f90`, `HDF5_ops.f90` — all persistent output goes through
  here, written as HDF5 with accompanying XDMF metadata so results can be opened directly in
  tools like ParaView (see the module's header notes on XDMF collections and HDF5 chunking).

## 5. Grid & equilibrium variable types
- `grid_vars.f90` — grid-sizing variables (`n_r_eq`, `n_r_X`, `n_r_sol`, field-line count
  `n_alpha`, field-line label `alpha`, angular ranges).
- `eq_vars.f90` — equilibrium variable *types* (`eq_1_type`, `eq_2_type`) and normalization
  constants (`R_0`, `B_0`, `psi_0`, `pres_0`, `rho_0`, `T_0`).
- `grid_utilities.f90`, `grid_ops.f90` — grid construction, redistribution across MPI ranks,
  field-aligned grid setup, trimming ghost regions (`trim_grid`).

## 6. Equilibrium physics
- `eq_utilities.f90`, `eq_ops.f90` — the largest module in the codebase (`eq_ops.f90`, ~360 KB).
  Computes derived equilibrium quantities (`calc_derived_q`), normalization, flux-surface
  quantities, and various diagnostic plots (`B_plot`, `J_plot`, `kappa_plot`, `flux_q_plot`).
  Also handles splitting equilibrium work into "eq jobs" for memory-constrained parallel runs
  (`divide_eq_jobs`, `calc_eq_jobs_lims`).

## 7. Equilibrium code interfaces (VMEC / HELENA)
- `VMEC_vars.f90`, `VMEC_utilities.f90`, `VMEC_ops.f90` — reads VMEC NetCDF output via LIBSTELL's
  `read_wout_mod` (Fourier-mode equilibrium quantities: `R`, `Z`, `B`, currents, rotational
  transform, etc.).
- `HELENA_vars.f90`, `HELENA_ops.f90` — reads HELENA plain-text equilibrium output
  (`read_HEL`) and interpolates it onto PB3D's internal grids (`interp_HEL_on_grid`).
- Which one is used at runtime is controlled by `eq_style` in `num_vars.f90` (1 = VMEC, 2 =
  HELENA), not by a compile-time switch.

## 8. Perturbation ("X") physics
- `X_vars.f90`, `X_utilities.f90`, `X_ops.f90` — sets up poloidal/toroidal perturbation mode
  numbers (`setup_modes`, `init_modes`), computes vectorial/tensorial perturbation quantities
  (`calc_X`) and magnetic integrals (`calc_magn_ints`) that feed the eigenvalue problem. Mode
  setup has a long, bug-fix-heavy history in the changelog (v2.27, v2.34, v2.35, v2.36) — treat
  changes here as high-risk to correctness.

## 9. Solution physics
- `sol_vars.f90`, `sol_utilities.f90`, `sol_ops.f90` — solution vector post-processing: energy
  decomposition (`decompose_energy`), plotting eigenvectors/eigenvalues (`plot_sol_vec`,
  `plot_sol_vals`, `plot_harmonics`).

## 10. Vacuum response
- `vac_vars.f90`, `vac_utilities.f90`, `vac_ops.f90`, `dtorh.f90` — Boundary Element Method vacuum
  response outside the plasma, solved via STRUMPACK-Dense. **`vac_ops.f90`'s own header marks
  this module as "still under construction and not usable yet."** Treat any related bug report or
  change with extra scrutiny and verify against that caveat first.

## 11. Richardson extrapolation
- `rich_vars.f90`, `rich_ops.f90` — drives convergence-order extrapolation across multiple
  "Richardson levels" of increasing perturbation-grid resolution (`init_rich`, `do_rich`,
  `start_rich_lvl`/`stop_rich_lvl`, `calc_rich_ex`, `term_rich`). This is the outer loop of the
  whole simulation (see [Simulation Pipeline](../workflows/simulation-pipeline.md)).

## 12. SLEPc interface
- `SLEPC_utilities.f90`, `SLEPC_ops.f90` — wraps PETSc/SLEPc to assemble the generalized
  eigenvalue problem matrices (`setup_mats`), apply boundary conditions (`set_BC`), configure the
  solver (`setup_solver`), and extract results (`get_solution`, `summarize_solution`,
  `store_results`). Note the module docstring: routines here require a **trimmed** solution grid.

## 13. Input/output orchestration
- `input_utilities.f90`, `input_ops.f90` — namelist parsing (`inputdata_PB3D`,
  `inputdata_POST`) and printing input-derived output. See
  [Input & Configuration](../workflows/input-and-configuration.md).
- `output_ops.f90` — higher-level output helpers built on `HDF5_ops`.
- `files_ops.f90` — command-line argument parsing (`parse_args`), opening input/output files
  (`open_input`, `open_output`), and the full list of supported `--option` flags per program style.

## 14. PB3D-output specific operations
- `PB3D_utilities.f90`, `PB3D_ops.f90` — reconstructing previously-written PB3D state from HDF5
  (`reconstruct_PB3D_in`, `reconstruct_PB3D_grid`, `reconstruct_PB3D_eq_1/2`,
  `reconstruct_PB3D_sol`), used both for Richardson restarts within PB3D and for POST reading
  PB3D's output.

## 15. Drivers (orchestration layer)
- `driver_eq.f90` → `run_driver_eq`: sets up the equilibrium grid and quantities for one
  Richardson level / equilibrium job.
- `driver_X.f90` → `run_driver_X`: sets up perturbation grid and quantities.
- `driver_sol.f90` → `run_driver_sol`: sets up the solution grid, calls into `SLEPC_ops` to solve
  the eigenvalue problem, and applies Richardson extrapolation.
- `driver_POST.f90` → `init_POST` / `run_driver_POST` / `stop_POST`: the entire postprocessing
  pipeline (by far the largest driver file, ~74 KB): reconstructs PB3D state, builds
  extended/field-aligned output grids, and produces the various 1-D/plot outputs.

See [Simulation Pipeline](../workflows/simulation-pipeline.md) for how these drivers are called in
sequence.

## 16. Test module
- `test.f90` — `generic_tests()`, only compiled/callable when built with `ldebug` and run with
  the `--test`/`-t` CLI flag. See [Testing & Extension](../testing-and-extension.md).

## Vendored libraries (`Libraries/`, not `Modules/`)
- `dfft.f` — double-precision FFT package, built as the `dfftpack` static library.
- `foul.f90` — formatted terminal output library, built as the `foul` static library.
