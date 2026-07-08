# Physics Domains

This page explains the scientific/business logic behind the module groups in
[Module Map](../architecture/module-map.md) — *why* the code is structured this way, not just
which files exist.

## The problem being solved

PB3D computes **linear ideal-MHD peeling-ballooning stability**: given a 3-D toroidal plasma
equilibrium, it determines whether small perturbations grow (unstable) or oscillate/decay
(stable), by discretizing the linearized MHD energy functional along field lines and solving the
resulting generalized eigenvalue problem. This is standard fusion-plasma stability analysis,
extended here from the traditional 2-D axisymmetric (tokamak) case to fully 3-D equilibria
(tokamaks with ripple, or stellarators). See the citations in the `PB3D.f90` header
(`weyens2014theory`, `Weyens2017PB3D`) for the underlying theory papers — the code implements
that theory.

## Equilibrium (`eq_*`, `VMEC_*`, `HELENA_*`)

PB3D does not compute an equilibrium itself — it **imports** one from an external equilibrium
code:
- **VMEC** (3-D ideal-MHD equilibrium code, widely used for stellarators and non-axisymmetric
  tokamaks): read via LIBSTELL's `read_wout_mod` NetCDF reader in `VMEC_ops.f90`.
- **HELENA** (2-D axisymmetric fixed-boundary equilibrium code): read via `HELENA_ops::read_HEL`
  from its plain-text output.

Both are converted into PB3D's internal representation (`eq_1_type`, `eq_2_type` in
`eq_vars.f90`) and normalized to dimensionless units (`calc_normalization_const`,
`normalize_input` in `eq_ops.f90`, using reference scales `R_0`, `B_0`, `psi_0`, `pres_0`,
`rho_0`, `T_0`). `eq_ops.f90` (the largest module, ~360 KB) then derives all quantities needed
for stability analysis: metric coefficients, magnetic shear, curvature, parallel current, safety
factor/rotational transform, etc. (`calc_derived_q`). Several of these derivations have a long
bug-fix history in the changelog (parallel current in 3-D, v2.32; shear sign errors, v2.37,
v2.40) — when touching `calc_derived_q`, check the changelog entries for the relevant quantity
first.

## Perturbation ("X") (`X_*`)

The perturbation is expanded in poloidal/toroidal Fourier harmonics along field-aligned
coordinates. `X_ops.f90::setup_modes` determines, for each flux surface, which secondary mode
numbers couple to the primary mode (`prim_X`) — this coupling structure is what makes the
peeling-ballooning problem "3-D" (in a purely axisymmetric equilibrium, modes wouldn't couple).
`calc_X` computes the vectorial (`X_1_type`) and tensorial (`X_2_type`) perturbation quantities;
`calc_magn_ints` integrates these over the magnetic geometry to build the matrices later handed to
SLEPc. Mode setup correctness has historically been fragile (see changelog v2.27, v2.34–v2.36) —
this is one of the higher-risk areas to modify without strong regression checks.

## Vacuum response (`vac_*`)

Outside the plasma boundary, the perturbed magnetic field must match a vacuum solution. PB3D
computes this via a **Boundary Element Method**, either using a fully 3-D field-aligned
collocation approach or (for axisymmetric cases) an analytical toroidal Green's-function
integration (see the `vac_ops.f90` module header for the cited references). The resulting dense
linear system — for the exterior (vacuum) side, `(H + 4 pi I) Phi = G dPhi` — is solved with
STRUMPACK-Dense when available, or with ScaLAPACK LU otherwise.

**Verification status**: the axisymmetric building blocks (Green's function kernels, singular
integrals, assembled `G`/`H`, boundary potential solve) are covered by the full-stack test suite
(`tests/fullstack`, see `Documentation/testing.md`); the free-boundary chain is wired for
`BC_style(2) = 4` but not yet benchmarked end-to-end, and the 3-D field-line style is untested
beyond its far-field kernels. Treat vacuum-related work as experimental; verify
current status in `vac_ops.f90` before relying on or extending it.

## Solution (`sol_*`, `SLEPC_*`)

The assembled perturbation matrices (plus vacuum response, plus boundary conditions) define a
generalized eigenvalue problem `A x = λ B x`, solved with **SLEPc** (`SLEPC_ops.f90`):
`setup_mats` builds the matrices, `set_BC` imposes boundary conditions, `setup_solver`/
`setup_guess` configure the solver (with an optional eigenvalue guess `EV_guess`/`EV_style`),
`get_solution`/`summarize_solution`/`store_results` extract and persist results. `sol_ops.f90`
then post-processes solved eigenvectors: decomposing the perturbed energy into potential/kinetic
contributions (`decompose_energy`) and producing harmonic/eigenvector plots.

A positive-growth-rate (unstable) eigenvalue indicates the equilibrium is peeling-ballooning
unstable at that mode number and radial location — this is the ultimate physical answer PB3D is
built to produce.

## Richardson extrapolation (`rich_*`)

Because the perturbation/solution grids are finite-resolution, results converge to the true
continuum answer only in the limit of infinite resolution. `rich_ops.f90` runs PB3D at
successively refined grids (up to `max_it_rich` levels) and extrapolates (`calc_rich_ex`) to
improve the estimate and its error bound, stopping early once `tol_rich` is met (`do_rich`). This
is the outer loop of the whole program — see
[Simulation Pipeline](../workflows/simulation-pipeline.md).

## Grids & storage (`grid_*`, `HDF5_*`)

PB3D operates on multiple related grids per run: the equilibrium grid, the (often field-aligned)
perturbation grid, and the solution grid, each possibly redistributed across MPI ranks and
"trimmed" of ghost regions before final use (`grid_utilities.f90`, `grid_ops.f90`). All persistent
state is written as parallel HDF5 with XDMF sidecar metadata (`HDF5_ops.f90`) so results are
directly viewable in tools like ParaView without a custom reader — a deliberate choice documented
in that module's header notes on hyperslab selection and chunking.

## Postprocessing (`driver_POST`)

POST does not recompute physics — it **reconstructs** PB3D's internal state from HDF5
(`PB3D_ops::reconstruct_PB3D_*`) and derives visualization-ready output: real-space field/current/
curvature plots, resonance surfaces, energy reconstruction, and (for HELENA-derived cases)
comparison tables against HELENA's native output. See `driver_POST.f90`, by far the largest driver
(~74 KB), for the full set of outputs gated by the `plot_*` namelist flags.
