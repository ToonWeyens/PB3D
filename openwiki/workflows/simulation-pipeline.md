# Simulation Pipeline

This page traces what actually happens at runtime for both executables. Source of truth:
[`/PB3D.f90`](/PB3D.f90) and [`/POST.f90`](/POST.f90).

## Running PB3D

```
PB3D USER_INPUT EQUILIBRIUM_INPUT [OPTIONS]
```

- `USER_INPUT` — a namelist file with group `inputdata_PB3D` (see
  [Input & Configuration](input-and-configuration.md)); if it can't be opened, PB3D falls back to
  defaults ("No input file found. Default used").
- `EQUILIBRIUM_INPUT` — VMEC NetCDF output or HELENA plain-text output, in `.txt` format per the
  usage string built in `files_ops.f90`'s `parse_args`.
- `[OPTIONS]` — flags such as `--no_guess`, `--jump_to_sol`, `--export_HEL`,
  `--plot_VMEC_modes`, `--invert_top_bottom_H`, `--no_plots`, `--no_output`,
  `--do_execute_command_line`, `--mem_info`, `-t`/`--test`, plus pass-through PETSc/SLEPc flags
  (`-eps_type`, `-st_pc_type`, `-log_view`, etc.). Full list and argument counts are in
  `files_ops.f90::init_files`.
- `-h` / `--help` prints usage and exits.

### Startup sequence (`PB3D.f90`)

1. `start_MPI()` — must succeed before anything else.
2. `print_hello()`, `init_output()`, `init_files()`, `init_time()`, `init_HDF5()`.
3. Master rank (`rank == 0`) only:
   - `parse_args()` → `open_input()` → `read_input_opts()`.
   - If `rich_restart_lvl == 1` (starting fresh, not resuming a Richardson level): `read_input_eq()`
     reads the equilibrium file, `calc_normalization_const()` + `normalize_input()` normalize
     physical units, `open_output()` creates the HDF5 output, `print_output_in('in')` records the
     input, then `dealloc_in()` frees the raw equilibrium-code input.
   - Otherwise (resuming a run): just `open_output()`.
4. `broadcast_input_opts()` — sends all input options from master to every MPI rank.
5. **(debug builds only, `ldebug`)** if `ltest` (`-t`/`--test` was passed): run
   `test::generic_tests()` and stop after — see [Testing & Extension](../testing-and-extension.md).
6. `init_rich()` — initialize the Richardson-extrapolation state machine.

### Main loop

```
RICH: do while (do_rich())                     ! Richardson levels (grid refinement)
    start_rich_lvl()
    PAR: do while (do_eq())                    ! "equilibrium jobs" (memory-limited chunks)
        run_driver_eq(...)                     ! equilibrium quantities on this grid/job
        run_driver_X(...)                      ! perturbation modes + magnetic integrals
    end do PAR
    run_driver_sol(...)                        ! assemble + solve eigenvalue problem (SLEPc),
                                                ! Richardson-extrapolate results
    stop_rich_lvl()
end do RICH
term_rich()
```

- The **outer `RICH` loop** implements Richardson extrapolation: PB3D can be run at increasing
  perturbation-grid resolutions (`max_it_rich` levels, `tol_rich` convergence tolerance) to
  extrapolate to the continuum limit. See `rich_ops.f90`.
- The **inner `PAR` loop** iterates over "equilibrium jobs" — the equilibrium/perturbation grid can
  be split into chunks (`eq_jobs_lims`, `divide_eq_jobs` in `eq_ops.f90`) to bound peak memory use
  (`max_tot_mem`, `max_X_mem` in `num_vars.f90`), driven by `do_eq()`.
- `run_driver_sol` is where the vacuum response (`vac_ops::calc_vac_res`) and the SLEPc eigenvalue
  solve (`SLEPC_ops::solve_EV_system_SLEPC`) are actually invoked, followed by Richardson
  extrapolation of the solution (`rich_ops::calc_rich_ex`).

### Shutdown

`stop_MPI(...)` deallocates all major state (grids, equilibrium, perturbation, vacuum, solution
types) and finalizes MPI; `close_output` closes the HDF5 file; `print_goodbye` prints final
timing.

## Running POST

```
POST USER_INPUT PB3D_OUTPUT [OPTIONS]
```

- `USER_INPUT` — namelist file with group `inputdata_POST`.
- `PB3D_OUTPUT` — the `.h5` file produced by a prior PB3D run (must be at least
  `min_PB3D_version` as checked against `prog_version` in `num_vars.f90`).
- Options include `--swap_angles` and `--compare_tor_pos` (see `files_ops.f90::init_files`, case 2).

### Sequence (`POST.f90`)

1. `start_MPI()`, `print_hello()`, `init_output()`, `init_files()`, `init_time()`, `init_HDF5()`.
2. Master rank: `parse_args()` → `open_input()` → `reconstruct_PB3D_in('in')` (rehydrate PB3D's
   recorded input state from HDF5) → `read_input_opts()` → `open_output()` → `dealloc_in()`.
3. `broadcast_input_opts()`.
4. `init_POST()` → `run_driver_POST(...)` → `stop_POST()` (see `driver_POST.f90`): reconstructs
   full grids/equilibrium/solution state, builds either an "extended" or "field-aligned" output
   grid (`POST_style`), and produces the requested plots (magnetic field, current, curvature,
   flux quantities, solution eigenvectors, energy reconstruction, vacuum potential — all gated by
   the `plot_*` flags in `num_vars.f90`).

## Practical notes for future changes

- Both programs share the exact same startup boilerplate (MPI/output/files/time/HDF5 init) — if
  you need to change initialization order, change it in both `PB3D.f90` and `POST.f90` or you will
  create drift between the two.
- `rich_restart_lvl` is the switch between "fresh run" and "restart at a higher Richardson level";
  logic that assumes a fresh run (e.g. reading/normalizing the equilibrium) is guarded by
  `if (rich_restart_lvl.eq.1)` and must stay guarded that way for restarts to work.
- All `CHCKERR` calls immediately return/abort on nonzero `ierr` (`ierr == 66` is the "silent
  stop" convention used for `--help`/usage errors, handled specially — see
  [`/include/PB3D_macros.h`](/include/PB3D_macros.h)).
