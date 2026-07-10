# Input & Configuration

PB3D and POST are both configured through **Fortran namelist** files plus a fixed equilibrium
input file and (for PB3D) optional perturbation-shape files. Runtime behavior throughout the
codebase is controlled by a large set of "style" integers and flags defined centrally in
[`/Modules/num_vars.f90`](/Modules/num_vars.f90).

## The two namelists

Defined and parsed in [`/Modules/input_ops.f90`](/Modules/input_ops.f90)
(`read_input_opts`):

- **`inputdata_PB3D`** — used by the `PB3D` executable. Example:
  [`/Examples/input_cbm18a`](/Examples/input_cbm18a).
- **`inputdata_POST`** — used by the `POST` executable. Example:
  [`/Examples/input_POST`](/Examples/input_POST).

Both namelist blocks can be followed by a **second, unlabeled block** in the same file (see the
tail of `input_cbm18a` and `input_POST`) that overrides a subset of variables — this is used to
respecify perturbation mode ranges (`n_mod_X`, `min_sec_X`, `max_sec_X`) or Richardson-restart
(`PB3D_rich_lvl`) without duplicating the whole file. Check `input_ops.f90` for exactly how the
second read is triggered before assuming its semantics.

### Commonly-tuned `inputdata_PB3D` variables (see `/Examples/input_cbm18a`)

| Variable | Meaning |
|---|---|
| `min_par_X` / `max_par_X` | Extent of the field-aligned parallel coordinate (in units of π) |
| `min_r_sol` / `max_r_sol` / `n_r_sol` | Radial (normal-coordinate) range and resolution of the solution grid |
| `prim_X` | Primary (dominant) perturbation mode number |
| `n_mod_X`, `min_sec_X`, `max_sec_X` | Number/range of secondary (coupled) perturbation modes |
| `alpha` / `alpha_style` | Field-line label, and whether one field line with many turns (style 1, legacy) or many field lines with one turn each (style 2) is used. **Style 2 is recommended for 3-D runs**: for the surface-averaged coefficients PB3D computes, its explicit field-line-label grid converges spectrally, whereas the single-line (Weyl) sampling of style 1 converges only at the equidistribution rate, is frozen during Richardson extrapolation (the line length is fixed), and is not supported by the free-boundary vacuum. Style 1 is kept as a cross-check and for axisymmetric cases (where the label is ignorable and one line is exact) |
| `EV_style`, `EV_guess`, `EV_BC` | Eigenvalue-solver method selection and boundary artificial eigenvalue |
| `BC_style(2)` | Boundary-condition style at the two radial boundaries: 1 = eigenvector zeroed (fixed boundary), 4 = explicit surface-energy minimization with the vacuum response (the verified free-boundary style, and the default at the edge). Style 2 is an unimplemented placeholder (reserved for a Hermitian variational variant), style 3 was removed |
| `tol_rich`, `max_it_rich`, `rich_restart_lvl` | Richardson extrapolation tolerance, max levels, and restart level (`1` = fresh start) |
| `tol_SLEPC`, `max_it_slepc` | SLEPc eigenvalue solver tolerance/iteration cap (per Richardson level, hence the array) |
| `use_pol_flux_F`, `use_normalization` | Coordinate/normalization convention switches |
| `max_tot_mem` | Memory budget (MB) across all MPI processes — drives automatic job-splitting (`divide_eq_jobs`, `eq_jobs_lims`) |
| `plot_*` (e.g. `plot_B`, `plot_flux_q`, `plot_resonance`) | Toggle diagnostic plots written alongside normal output |

### Commonly-tuned `inputdata_POST` variables (see `/Examples/input_POST`)

| Variable | Meaning |
|---|---|
| `n_sol_plotted` | Which solved eigenmodes to post-process/plot (`-1` values are sentinel/"skip") |
| `POST_style` | `1` = extended grid, `2` = field-aligned grid, for output |
| `plot_sol_xi`, `plot_sol_Q`, `plot_E_rec` | Plasma displacement, perturbed field, and energy-reconstruction plots |
| `pert_mult_factor_POST` | Scale factor applied when visualizing a perturbed equilibrium |
| `PB3D_rich_lvl` (second block) | Which Richardson level of a PB3D run to post-process |

## Equilibrium input

The second CLI argument to `PB3D` (`EQUILIBRIUM_INPUT`) is:
- a **VMEC** NetCDF `wout_*.nc`-style output (read via LIBSTELL's `read_wout_mod`,
  `VMEC_ops.f90`), or
- a **HELENA** plain-text output file (read via `HELENA_ops::read_HEL`).

Which reader is used is controlled by `eq_style` (set in the namelist or inferred — check
`files_ops.f90::open_input` / `input_ops.f90` for exact detection logic before changing it).
`files_ops.f90::open_input` also explicitly checks for and rejects equilibrium filenames
containing a stray `.` (a known past bug source, see README changelog v2.41).

## Perturbation shape files (`Examples/pert_*.dat`)

Optional files describing an explicit non-axisymmetric boundary perturbation as a sum of
poloidal/toroidal Fourier harmonics, e.g. [`/Examples/pert_topbot.dat`](/Examples/pert_topbot.dat):

```
#   N   M   delta_c     delta_s
    4   5   0.05        0.0
    4   7   0.0        -0.025
    4   3   0.0         0.025
```

Each row is one `(N, M)` harmonic with cosine (`delta_c`) and sine (`delta_s`) amplitudes. These
feed equilibrium export/perturbation routines such as `create_VMEC_input` in `eq_ops.f90` (see
README changelog v2.31/v2.32/v2.42 for the non-trivial history of getting mode-shifting and
symmetry right in this code path — a good example of "why" context before touching it again).

## Parameter-scan input (`Examples/array_input`)

`Examples/array_input` is a lightweight scan-specification format (not a Fortran namelist):
comment lines start with `#`; active lines list comma-separated `variable = value` overrides
(e.g. `prim_X = 25, n_mod_X = 25`) — one line per scan point. This format is consumed by external
run-automation tooling in the separate `PB3D_tools` repository, not by `PB3D`/`POST` directly (see
[Quickstart](../quickstart.md) for that repo split).

## Where to look when adding a new option

1. Add the variable to `num_vars.f90` (with a `public` export and doc comment).
2. Add it to the relevant `namelist /inputdata_PB3D/` or `/inputdata_POST/` list and to the
   `use num_vars, only: ...` clause in `input_ops.f90::read_input_opts`.
3. If it must reach worker MPI ranks, add it to the broadcast list in `MPI_ops.f90`
   (`broadcast_input_opts`).
4. Update an `Examples/` input file if it's a commonly-used option, so it's discoverable.
