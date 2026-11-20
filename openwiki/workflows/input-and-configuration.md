---
type: "Guide"
title: "Input & Configuration"
description: "How PB3D and POST are configured: the namelists inputdata_PB3D and inputdata_POST, the defaults and input checks of the main switches (eq_style detection, X_style and the mode-range rules, alpha_style, BC_style, norm_disc_prec_sol, vac_diag_style, check_vac), the POST vacuum plot inputs, the equilibrium and perturbation files, the example decks, and how to add a new option."
tags: [input, namelist, configuration, options, examples]
verified:
  - by: openwiki/0.9.0
    at: 2026-11-20T18:32:33.626Z
sources:
  - id: openwiki-source-1de3db34d54cdd02011c9204
    resource: repo://Modules/files_ops.f90
  - id: openwiki-source-99e71edd5b71e76c5f03dd84
    resource: repo://Modules/input_ops.f90
  - id: openwiki-source-0e5296407709fa9254158dd8
    resource: repo://Modules/num_vars.f90
generated: { by: "claude-code", at: "2026-11-20T18:32:33.626Z" }
---

# Input & Configuration

<!-- openwiki: broken internal link [/Modules/num_vars.f90] link "/Modules/num_vars.f90" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
<!-- openwiki: broken internal link [/Modules/input_ops.f90] link "/Modules/input_ops.f90" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
PB3D and POST each read one Fortran namelist from the input file given on the command line, plus the equilibrium file (PB3D) or the PB3D output (POST). The variables behind the namelists are mostly declared in [`/Modules/num_vars.f90`](/Modules/num_vars.f90), with the perturbation and grid variables in `X_vars.f90` and `grid_vars.f90`. [`/Modules/input_ops.f90`](/Modules/input_ops.f90) sets the defaults (`default_input_PB3D`, `default_input_POST`), reads the namelist, and checks and adapts the values (`adapt_run_PB3D`, `adapt_X_modes`, `adapt_alpha`, `adapt_plot`, ...).

## The two namelists

<!-- openwiki: broken internal link [/Examples/input_cbm18a] link "/Examples/input_cbm18a" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- **`inputdata_PB3D`**, read by `PB3D`. Example: [`/Examples/input_cbm18a`](/Examples/input_cbm18a).
<!-- openwiki: broken internal link [/Examples/input_POST] link "/Examples/input_POST" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
- **`inputdata_POST`**, read by `POST`. Example: [`/Examples/input_POST`](/Examples/input_POST).

Each program reads its namelist once. Fortran stops reading a namelist at the closing `/`, so the lines after it in the example decks are not read; they keep alternative values (for example a prescribed mode range) that you can move into the namelist.

### Main `inputdata_PB3D` variables

| Variable | Meaning | Default and checks |
|---|---|---|
| `min_par_X`, `max_par_X` | range of the parallel coordinate, in units of π | for HELENA a warning unless the range is 2 |
| `min_r_sol`, `max_r_sol`, `n_r_sol` | radial range and number of points of the solution grid | |
| `prim_X` | the primary mode number (toroidal n with the poloidal flux) | must exceed a minimum |
| `n_mod_X`, `min_sec_X`, `max_sec_X` | the secondary mode numbers | see "Mode range" below |
| `alpha`, `alpha_style`, `n_alpha` | field-line label; one field line with many turns (1) or `n_alpha` field lines with one turn each (2) | 2 with `n_alpha = 10` for VMEC, 1 for HELENA |
| `norm_disc_prec_sol` | order of the left finite difference in the radial direction | 3; must be 1 to 3 |
| `BC_style(1:2)` | boundary conditions at the inner end and the edge | `[1, 2]`; the inner end must be 1, the edge 1 (fixed) or 2 (free, with the vacuum) |
| `vac_diag_style` | diagonals of the vacuum operators: 1 indirect, from the identities; 2 direct, a diagnostic | 1; a warning for 2 with a free boundary |
| `check_vac` | test the vacuum operators with exact potentials and print `CHECKVAC` lines | `.false.` |
| `EV_style`, `EV_guess`, `EV_BC` | eigenvalue method, first guess, and the artificial eigenvalue of the zeroed rows | `EV_style = 1` (SLEPc) |
| `n_sol_requested`, `retain_all_sol` | number of solutions to compute, and whether to keep faulty solutions too | |
| `tol_SLEPC`, `max_it_slepc` | SLEPc tolerance (one value per Richardson level) and iteration limit | |
| `max_it_rich`, `tol_rich`, `rich_restart_lvl` | number of Richardson levels, tolerance, and the level to restart from (1 = fresh run) | 1 level, 1e-4 |
| `use_pol_flux_F`, `use_normalization`, `norm_style` | flux coordinate and normalisation | |
| `R_0`, `B_0`, `psi_0`, `pres_0`, `T_0`, `rho_0` | override the computed normalisation constants | computed when not given |
| `max_tot_mem` | memory limit in MB over all processes; drives the division into equilibrium jobs | |
| `plot_*` | diagnostic plots of the equilibrium, such as `plot_B`, `plot_flux_q`, `plot_resonance` | |

The option `norm_disc_style_sol` and the edge boundary style 4 no longer exist: the radial derivative is always a left finite difference, and a deck that asks for `BC_style(2) = 4` is refused with a message and an exit status that is not zero.

### Mode range

`setup_modes` uses one of two styles for the secondary mode numbers, and `input_ops.f90` chooses the style from what the deck sets:

- `min_sec_X` and `max_sec_X` set, `n_mod_X` not set: `X_style = 1`, the prescribed range on all flux surfaces. Both limits must be set, and `|max_sec_X| ≥ |min_sec_X|`.
- Neither limit set: `X_style = 2`, a window of `n_mod_X` numbers around the resonant value on each surface (10 if `n_mod_X` is not set either).
- A limit together with `n_mod_X`, or only one of the two limits: the deck is refused with a message.

### Main `inputdata_POST` variables

| Variable | Meaning |
|---|---|
| `n_sol_plotted` | which solutions to post-process |
| `PB3D_rich_lvl` | which Richardson level of the PB3D output to read; below 1, only the equilibrium output is done |
| `POST_style` | 1 = extended grid (default), 2 = field-aligned grid |
| `plot_sol_xi`, `plot_sol_Q`, `plot_E_rec` | plasma displacement, magnetic perturbation, and energy reconstruction |
| `plot_vac_pot`, `n_vac_plot`, `min_Rvac_plot` … `max_Zvac_plot` | the vacuum potential plot: on by default, 101 × 101 points, and a box in metres that by default follows the boundary (see [Vacuum Response](../domain/vacuum-response.md#post-vacuum-potential-plot)); off for a fixed boundary |
| `pert_mult_factor_POST` | factor applied to the perturbation when a perturbed equilibrium is plotted |

## Command line and PETSc options

`PB3D <input> <equilibrium> [options]` and `POST <input> <PB3D output> [options]`. `files_ops.f90` (`parse_args`) lists the accepted options, including `--test` (debug builds), `--do_execute_command_line`, and PETSc and SLEPc options such as `-st_pc_factor_mat_solver_type`, which are passed on to PETSc.

## Equilibrium input

PB3D recognises the equilibrium file by its content: a name ending in `.nc` is a VMEC `wout` file (`eq_style = 1`, read through LIBSTELL's `read_wout_mod`), and a file whose first item reads as an integer is a HELENA mapping file (`eq_style = 2`, read by `read_HEL`). `open_input` in `files_ops.f90` finds the extension from the last `.` in the name. When PB3D continues from earlier output, it takes `eq_style` from that output.

## Perturbation files (`Examples/pert_*.dat`)

<!-- openwiki: broken internal link [/Examples/pert_topbot.dat] link "/Examples/pert_topbot.dat" is root-absolute, which no real consumer resolves against the repository root (not a coding agent reading the page, not GitHub's Markdown renderer, not a local viewer); use a path relative to this file instead. Fix the href or restore the target, then delete this comment. -->
These files describe a non-axisymmetric boundary perturbation as a sum of Fourier harmonics, for example [`/Examples/pert_topbot.dat`](/Examples/pert_topbot.dat):

```
#   N   M   delta_c     delta_s
    4   5   0.05        0.0
    4   7   0.0        -0.025
    4   3   0.0         0.025
```

Each row is one harmonic (N, M) with a cosine amplitude `delta_c` and a sine amplitude `delta_s`. They feed the equilibrium export routines such as `create_VMEC_input` in `eq_ops.f90`; the README changelog records the earlier problems with mode shifts and symmetry in that code.

## Example and regression decks

- `Examples/` holds sample decks: `input_cbm18a` and `input_cbm18a_ripple`, `input_Hmode` and `input_Hmode_ripple`, `input_cdxu`, `input_qps`, and `input_POST`.
- `tests/regression/input_*` are the decks of the regression tests, with known eigenvalues: `input_cbm18a_fixed_bnd`, `input_cbm18a_free_bnd`, `input_circ_eps005`, `input_circ_eps033`. They are the best-tested starting points.
- `Examples/array_input` is a scan format of the separate `PB3D_tools` repository (one line of `variable = value` overrides per scan point), not read by PB3D or POST.

## Where to look when adding a new option

1. Declare the variable in `num_vars.f90` (or the relevant `*_vars` module), public, with a doc comment.
2. Add it to the namelist `inputdata_PB3D` or `inputdata_POST` in `input_ops.f90`, to the `use` list there, and give it a default in `default_input_PB3D` or `default_input_POST`. Check its values in the `adapt_*` routines of `input_ops.f90`, and stop with a message for a value that is not allowed.
3. If other MPI ranks need it, add it to `broadcast_input_opts` in `MPI_ops.f90`.
4. If the option changes what PB3D writes, check whether POST and a Richardson restart still read the output, and change `min_PB3D_version` if not (see [Testing & Extension](../testing-and-extension.md#versioning)).
5. Add it to an example or regression deck if it is commonly used.
