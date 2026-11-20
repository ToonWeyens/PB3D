---
type: "Workflow"
title: "Simulation Pipeline"
description: "What happens when PB3D and POST run: the command lines, the start-up sequence, the Richardson loop with its equilibrium jobs and the drivers, where the vacuum geometry is stored and the vacuum response computed, the version check on restart and in POST, the POST sequence, and how errors end the run with a non-zero exit status."
tags: [workflow, pipeline, richardson, drivers, post, error-handling]
verified:
  - by: openwiki/0.9.0
    at: 2026-11-20T18:32:33.626Z
sources:
  - id: openwiki-source-f71e15dc61ec9541c2ea2c88
    resource: repo://include/PB3D_macros.h
  - id: openwiki-source-9676454832f9c783053e9eed
    resource: repo://Modules/driver_eq.f90
  - id: openwiki-source-4361ace621d24e5f0bcbb1a2
    resource: repo://Modules/driver_sol.f90
  - id: openwiki-source-1de3db34d54cdd02011c9204
    resource: repo://Modules/files_ops.f90
  - id: openwiki-source-2b951757a1f7314b5a33d0bf
    resource: repo://Modules/MPI_ops.f90
  - id: openwiki-source-ceb46d56be8dccc0770a071b
    resource: repo://Modules/PB3D_ops.f90
  - id: openwiki-source-ea0f6caa297b3dce09db8548
    resource: repo://Modules/rich_ops.f90
  - id: openwiki-source-5ebda9f0933cb8606c609e8f
    resource: repo://PB3D.f90
  - id: openwiki-source-8257526ded3898f9506f023a
    resource: repo://POST.f90
generated: { by: "claude-code", at: "2026-11-20T18:32:33.626Z" }
---

# Simulation Pipeline

<!-- openwiki: broken internal link [/PB3D.f90] link "/PB3D.f90" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/POST.f90] link "/POST.f90" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
This page traces what happens at run time in both programs. The sources are [`/PB3D.f90`](/PB3D.f90) and [`/POST.f90`](/POST.f90).

## Running PB3D

```
PB3D USER_INPUT EQUILIBRIUM_INPUT [OPTIONS]
```

- `USER_INPUT`: a file with the namelist `inputdata_PB3D` (see [Input & Configuration](input-and-configuration.md)). If the file cannot be opened, PB3D says "No input file found. Default used" and runs with the defaults.
- `EQUILIBRIUM_INPUT`: a VMEC `wout` file (`.nc`) or a HELENA mapping file; PB3D recognises which from the file itself.
- `[OPTIONS]`: PB3D's own options, such as `--no_guess`, `--jump_to_sol`, `--export_HEL`, `--plot_VMEC_modes`, `--invert_top_bottom_H`, `--no_plots`, `--no_output`, `--do_execute_command_line`, `--mem_info` and `-t`/`--test`, and PETSc/SLEPc options that are passed on (`-eps_type`, `-st_pc_type`, `-st_pc_factor_mat_solver_type`, `-log_view`, ...). `parse_args` in `files_ops.f90` holds the full list.
- `-h`/`--help` prints the usage and exits with status 0.

### Start-up (`PB3D.f90`)

1. `start_MPI()`, then `print_hello`, `init_output`, `init_files`, `init_time`, `init_HDF5`.
2. On rank 0 only: `parse_args` → `open_input` → `read_input_opts`. For a fresh run (`rich_restart_lvl = 1`): `read_input_eq` reads the equilibrium, `calc_normalization_const` and `normalize_input` make it dimensionless, `open_output` creates the HDF5 output, `print_output_in('in')` writes the input, and `dealloc_in` frees the equilibrium code's data. For a restart, only `open_output`.
3. `broadcast_input_opts()` sends the input options to every rank.
4. In debug builds with `--test`: `generic_tests()` runs and the program stops (see [Testing & Extension](../testing-and-extension.md)).
5. `init_rich()` sets up the Richardson loop. On a restart it reads the earlier output with `reconstruct_PB3D_in`, which stops with a message if that output is older than `min_PB3D_version`.

### Main loop

```
RICH: do while (do_rich())               ! Richardson levels
    start_rich_lvl()                     ! sets n_r_sol and the other per-level sizes
    PAR: do while (do_eq())              ! equilibrium jobs, sized to fit in memory
        run_driver_eq(...)               ! equilibrium on this job; vacuum geometry (store_vac) for a free boundary
        run_driver_X(...)                ! perturbation quantities and magnetic integrals
    end do PAR
    run_driver_sol(...)                  ! vacuum response, eigenvalue problem, Richardson extrapolation
    stop_rich_lvl()
end do RICH
term_rich()
```

- **The `RICH` loop** is Richardson extrapolation: up to `max_it_rich` levels of increasing resolution, stopping when the extrapolated eigenvalue changes by less than `tol_rich` (`rich_ops.f90`).
- **The `PAR` loop** runs over equilibrium jobs: the grid is divided into pieces (`divide_eq_jobs`, `calc_eq_jobs_lims` in `eq_ops.f90`) so that the memory stays below `max_tot_mem`.
- **The vacuum** enters in two places, only for a free boundary (`BC_style(2) = 2`): `run_driver_eq` stores the boundary geometry (`store_vac`), and `run_driver_sol` computes the response matrix (`calc_vac_res`) before it solves the eigenvalue problem (`solve_EV_system_SLEPC`). See [Vacuum Response](../domain/vacuum-response.md).
- **`run_driver_sol`** also extrapolates the eigenvalues over the levels (`calc_rich_ex`).

### Shutdown

`stop_MPI(...)` frees the grids, equilibrium, perturbation, vacuum and solution and finalises MPI; the HDF5 output is closed and `print_goodbye` prints the timing.

## Running POST

```
POST USER_INPUT PB3D_OUTPUT [OPTIONS]
```

- `USER_INPUT`: a file with the namelist `inputdata_POST`.
- `PB3D_OUTPUT`: the HDF5 output of a PB3D run, of version `min_PB3D_version` or newer.
- Options include `--swap_angles` and `--compare_tor_pos` (`parse_args`, program style 2).

### Sequence (`POST.f90`)

1. `start_MPI()`, `print_hello`, `init_output`, `init_files`, `init_time`, `init_HDF5`.
2. On rank 0: `parse_args` → `open_input` → `reconstruct_PB3D_in('in')`, which reads PB3D's recorded input and checks the version → `read_input_opts` → `open_output` → `dealloc_in`.
3. `broadcast_input_opts()`.
4. `init_POST` → `run_driver_POST` → `stop_POST` (`driver_POST.f90`). POST reconstructs the grids, equilibrium and solutions of the Richardson level `PB3D_rich_lvl`, builds the extended or field-aligned output grid (`POST_style`), and writes the plots and tables the `plot_*` switches ask for: the field, current and curvature, flux quantities, eigenvectors, the energy decomposition of every solution, and the vacuum potential. With `PB3D_rich_lvl` below 1 it does only the equilibrium output; for a level above the highest one in the file it stops and names the highest level found.

## Errors and exit status

<!-- openwiki: broken internal link [/include/PB3D_macros.h] link "/include/PB3D_macros.h" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
Every routine returns an error code `ierr`, and the `CHCKERR` macro returns at once when it is not zero (see [`/include/PB3D_macros.h`](/include/PB3D_macros.h)). In the main programs `CHCKERR` calls `sudden_stop`, which writes out the buffered output and ends the run through `abort_MPI` with a status that is not zero: the error code itself between 1 and 125, and 1 for every other code. The code 66 is the silent stop of the usage and `--help` paths: no error message is printed.

## Notes for changes

- Both programs have the same start-up code (MPI, output, files, time, HDF5). If you change the order, change it in `PB3D.f90` and `POST.f90` together.
- `rich_restart_lvl` separates a fresh run from a restart at a higher Richardson level. Code that assumes a fresh run, such as reading and normalising the equilibrium, sits behind `if (rich_restart_lvl.eq.1)` and has to stay there for restarts to work.
- A change to what PB3D writes can break POST and restarts on older output; raise `min_PB3D_version` when the layout changes (see [Testing & Extension](../testing-and-extension.md#versioning)).
