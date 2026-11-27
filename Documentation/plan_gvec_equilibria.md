# Plan: GVEC equilibria in PB3D (issue #104)

Status: plan, 27 November 2026. Nothing is implemented. The decisions in section 7 are Toon's.

## 1. What PB3D needs from an equilibrium

PB3D reads the equilibrium once, in `read_VMEC` (`Modules/VMEC_ops.f90`), and keeps it in the arrays of `Modules/VMEC_vars.f90`. Everything after that works on these arrays.

| array in `VMEC_vars` | what it holds | wout variables it comes from |
|---|---|---|
| `R_V_c`, `R_V_s`, `Z_V_c`, `Z_V_s` | Fourier coefficients of $R$ and $Z$ on VMEC's full mesh, with radial derivatives of order 0 to 3 | `rmnc`, `rmns`, `zmnc`, `zmns` |
| `L_V_c`, `L_V_s` | Fourier coefficients of $\lambda$, moved from the half mesh to the full mesh by a cubic spline that extrapolates at both ends, with radial derivatives of order 0 to 3 | `lmnc`, `lmns` |
| `jac_V_c`, `jac_V_s` | Fourier coefficients of the Jacobian on the Nyquist mode list, moved from the half mesh in the same way, with radial derivatives of order 0 to 2 | `gmnc`, `gmns` |
| `B_V_sub_c`, `B_V_sub_s` | covariant components $B_s$, $B_u$, $B_v$ on the Nyquist list; only the debug tests use them | `bsubsmns`, `bsubumnc`, `bsubvmnc` and the asymmetric ones |
| `flux_t_V`, `flux_p_V`, `rot_t_V`, `q_saf_V`, `pres_V` | toroidal flux, poloidal flux, $\iota$, $q$, pressure, with derivatives of order 0 to 2 | `phi`, `phipf`, `iotaf`, `presf` |
| `B_0_V` | the (0,0) Fourier mode of $\lvert B \rvert$, extrapolated to the axis from the half mesh | `bmnc` |
| debug only | $\lvert B \rvert$ on the Nyquist list, the integrated currents (`J_plot`) | `bmnc`, `bmns`, `jcuru`, `jcurv` |
| scalars and lists | sizes, symmetry, mode lists; `rmin_surf`, `rmax_surf` and `pres_V(1,0)` set the normalization of `norm_style` 2, `betaxis` is printed | `ns`, `mnmax`, `mnmax_nyq`, `mpol`, `ntor`, `nfp`, `lasym`, `lfreeb`, `lrfp`, `version_`, `aspect`, `betaxis`, `gamma`, `rmax_surf`, `rmin_surf`, `zmax_surf`, `xm`, `xn`, `xm_nyq`, `xn_nyq` |

Facts about these arrays that a GVEC route has to respect:

- The radial coordinate is the normalized toroidal flux $s$ on VMEC's `ns` surfaces, equidistant in $s$. The equilibrium grid of PB3D is a subset of these surfaces: `calc_RZL_ind` (`eq_ops.f90`) indexes `R_V_c(:, i_min:i_max, :)` with the grid's radial indices. So PB3D's values of $R$ and $Z$ at its own grid points are the file's values; $\lambda$ and the Jacobian reach the full mesh through a spline from the half mesh, with extrapolation at both ends, and all radial derivatives come from cubic splines (`norm_disc_prec_eq = 3`, the only value `input_ops.f90` accepts). The third radial derivative of a cubic spline is constant between surfaces and jumps at them.
- The Jacobian comes from `gmnc`, not from $R$ and $Z$ (`calc_jac_V_ind`). The vacuum normal uses `jac_FD`, which comes from `gmnc` extrapolated half a step to the boundary (`vac_ops.f90`, around line 430).
- The poloidal flux is the trapezoidal integral of `iotaf` times `phipf` (`VMEC_ops.f90:86-87`), with an error proportional to $h^2$, $h = 1/(\mathrm{ns}-1)$. With `use_pol_flux_F` it is the radial coordinate. No converter changes this.
- PB3D never reads `signgs`. It assumes VMEC's left-handed coordinates $(s, \theta, \varphi)$ with a negative Jacobian, for example in `calc_flux_q_VMEC` ("conversion VMEC LH -> PB3D RH", `eq_ops.f90:368`), in `calc_eq_2` (`pmone = -1`, line 522), in `grid_ops.f90:176` and 889, `eq_utilities.f90:786` and 866, `grid_utilities.f90:211` and 440.
- $R_0$ = `R_V_c(1,1,0)` (`calc_normalization_const_VMEC`) assumes that the mode (0,0) comes first in the main mode list, and `read_VMEC` stops if (0,0) is not the first Nyquist mode. `read_VMEC` checks `mnmax` $= (\mathrm{ntor}+1) + (2\,\mathrm{ntor}+1)(\mathrm{mpol}-1)$, VMEC's full rectangular mode set. The Nyquist list may have any size, in VMEC's order (`VMEC_ops.f90:161-166`).
- `PB3D_out.h5` stores the `VMEC_vars` arrays, and POST reads them back (`PB3D_ops.f90:162-309`); POST never opens the equilibrium file.

`files_ops.f90` treats a file whose name ends in `.nc` as VMEC (`eq_style` 1) and reads it with LIBSTELL's `read_wout_file`. What LIBSTELL's netCDF reader (`read_wout_nc` in `read_wout_mod.f90`, with `Ezcdf/ezcdf_inqvar.f90`) requires of a file:

- A missing variable is skipped silently (`cdfgv` and `cdfInqVar` return without an error), so the file needs only what PB3D uses; the VMEC++ fixture already lacks `lrecon`, `lmove_axis` and `mnmaxpot`.
- `ier_flag` must be 0 (or the "more iterations" or "jac75" flags), otherwise `ierr = -2`.
- `raxis_cc` and `zaxis_cs` must have length `ntor`+1, otherwise LIBSTELL stops.
- Logicals are int32 variables with the suffix `__logical__` (`lasym__logical__`).
- Integers must be NC_INT and reals NC_DOUBLE. Another type prints a message, sets no error and leaves the value unset: an int64 `ns` from numpy would be lost.
- With `version_` of 9 or more LIBSTELL reads `currumnc`; below 9 it computes the currents from the $B$ components.
- `nextcur` = 0.

## 2. What GVEC provides and its conventions

GVEC's unknowns are $X^1$, $X^2$ and $\lambda$ as functions of $(\rho, \vartheta, \zeta)$: B-splines of degree `X1X2_deg` and `LA_deg` in $\rho$ (continuity one less than the degree, the only choice GVEC allows) on `sgrid_nElems` elements, Fourier series in $\vartheta$ and $\zeta$. The position is $x = h(X^1, X^2, \zeta)$ with the map chosen by `which_hmap`. GVEC computes fixed-boundary equilibria only; PB3D's free-boundary stability uses a fixed-boundary equilibrium with a vacuum perturbation outside it, so this is no restriction.

Conventions, read from GVEC 1.5.1's source:

- GVEC's $\Phi$ is the toroidal flux divided by $2\pi$. The input `PhiEdge` is the total toroidal flux; `mhd3d.F90:137-138` divides it by $2\pi$ ("normalization like in VMEC"), and for a run from parameters $\Phi = \Phi_\mathrm{edge}\, s$ with $s = \rho^2$ (`mhd3d.F90:182`). GVEC's CAS3D writer multiplies $\Phi$ and $\chi$ by $2\pi$ (`python/gvec/scripts/cas3d.py:260-267`).
- GVEC's $\mu_0 = 4\pi \cdot 10^{-7}$ (`mhd3d.F90:117`), and the energy has $\lvert B \rvert^2/(2\mu_0)$ and $p$, so $p$ is in Pa, as `presf` in a wout file.
- For map 1, $(x, y, z) = (X^1\cos\zeta, -X^1\sin\zeta, X^2)$: $X^1 = R$, $X^2 = Z$, and $\zeta$ runs clockwise seen from above. $(\rho, \vartheta, \zeta)$ is right-handed.
- GVEC's VMEC reader (`InitVMEC`, `src/vmec/vmec.F90:74-124`) always sets $\zeta = -\varphi$. For VMEC's usual `signgs` = -1 it negates the sign of `gmnc` and of `iotaf` and maps the modes; for `signgs` = +1 it also reverses $\theta$. `pygvec convert-wout` passes `PhiEdge` = $-$`phi`(edge) (`scripts/convert_wout.py`). The CAS3D writer writes $-\iota$.

The transformations a writer from a GVEC state (map 1) to a wout file with `signgs` = -1 needs, with $s = \rho^2$, $\theta = \vartheta$, $\varphi = -\zeta$:

| wout quantity | from GVEC |
|---|---|
| $R(s, \theta, \varphi)$, $Z$, $\lambda$ | $X^1(\sqrt s, \theta, -\varphi)$, $X^2$, $\lambda$ at the same point |
| `phi` | $-2\pi\Phi$ |
| `phipf` | $d\,$`phi`$/ds$ |
| `iotaf` | $-\iota$ |
| `chi` | $-2\pi\chi$ (VMEC writes `chi` $=$ `signgs` $\int$ `chipf` $ds$, while `chipf` $=$ `iotaf` $\cdot$ `phipf` $= 2\pi\chi'$ in $s$) |
| `presf` | $p$ |
| `gmnc` ($\sqrt g$ in $(s, \theta, \varphi)$) | $-\mathcal J/(2\rho)$, negative on every surface |
| $B_s$, $B_u$, $B_v$ | $B_\rho/(2\rho)$, $B_\vartheta$, $-B_\zeta$ |
| $B^u$, $B^v$ | $B^\vartheta$, $-B^\zeta$ |

The signs of `phi`, `iotaf` and `chi` were checked on the circular VMEC twin (`wout_circ_eps033_q104.nc`, VMEC++, `signgs` = -1, 319 surfaces): `phi`(edge) = -0.32945, `iotaf` from -0.9605 to -0.9632, $d\,$`phi`$/ds$ = `phipf` and `chipf` = `iotaf` $\cdot$ `phipf` to rounding, `chi`(edge) = -0.3168 with $d\,$`chi`$/ds = -$`chipf`, and `bsupvmnc`(0,0) = -1.08 where GVEC's $B^\zeta = \Phi'(1+\partial_\vartheta\lambda)/\mathcal J$ is positive. GVEC's own VMEC reader divides `chi` by $2\pi$ without the factor `signgs` (`vmec_readin.F90:638`), so the inverse of that reader would give the wrong sign of `chi`. PB3D does not read `chi`, but other tools do. Gate A2a checks all signs again on the converted files.

pyGVEC (`gvec.load_state`, `state.evaluate`) evaluates $X^1$, $X^2$, $\lambda$, the metric, the Jacobian, $\lvert B \rvert$, $B_\vartheta$, $B_\zeta$, $\Phi$, $\chi$, $\iota$, $p$ at given $(\rho, \vartheta, \zeta)$, with radial derivatives of $X^1$, $X^2$, $\lambda$ up to second order only (`_evaluate.py:73`, the pattern `r{0,2}`). GVEC's CAS3D writer shows the pattern a wout writer can follow: evaluate on a grid, change the angle convention, Fourier transform with GVEC's own `ev2ft`.

GVEC does not write wout files, and I found no other tool that does.

## 3. Three routes

### Route A: a converter from a GVEC state to a wout file, PB3D unchanged

A Python script `gvec_to_wout.py` loads a GVEC state with pyGVEC and writes a wout file on `ns` surfaces equidistant in $s$, with the transformations of section 2:

- $R$, $Z$ on the full mesh and $\lambda$ on the half mesh at $s_{j-1/2} = (j - 3/2)/(\mathrm{ns} - 1)$, from a 2-D FFT of pyGVEC's values on a $(\theta, \varphi)$ grid into VMEC's basis $\cos(m\theta - n\varphi)$, $\sin(m\theta - n\varphi)$, in VMEC's mode order with (0,0) first. The grid is fine enough that the FFT returns the finite Fourier series exactly, so GVEC's internal mode order never enters. Modes outside GVEC's set are zero.
- `gmnc` on the half mesh on a Nyquist list large enough to hold the Jacobian exactly. For map 1, $\sqrt g = R\,(R_\theta Z_s - R_s Z_\theta)$ is a product of three finite Fourier series, so it is itself a finite series with $m \le 3\,m_\mathrm{max}$, $\lvert n \rvert \le 3\,n_\mathrm{max}$. VMEC truncates it (the circular fixture's Nyquist list stops at $m = 15$ for `mpol` = 12, where the exact Jacobian reaches $m = 33$); the converter need not. When $X^1$ and $X^2$ have different limits $M_1$, $M_2$, the bound is $2M_1 + M_2$ (and the same for $n$). For W7-X-sized resolutions the full list has about 2,400 modes, and `fourier2real` costs time in proportion to the number of modes at every grid point, so the converter drops the modes, in VMEC order, whose coefficients stay below a tolerance (for example $10^{-12}$ of the (0,0) coefficient) on all surfaces.
- `bmnc` and $B_s$, $B_u$, $B_v$ on the same list. These are not finite series, so the list truncates them; PB3D uses `bmnc` only for $B_0$ and the $B$ components only in debug tests.
- the profiles and scalars of section 1, `raxis_cc`, `zaxis_cs` of length `ntor`+1, `ier_flag` = 0, `lfreeb` = F, `lrfp` = F, `lasym` from GVEC's sine and cosine choice, `signgs` = -1, `nextcur` = 0, and the netCDF types LIBSTELL needs. Writing the file through `vmecpp.VmecWOut` (already in PB3D's `vmec` dependency group) gives a layout LIBSTELL is known to read; the converter then checks the types.
- before writing, the converter checks that the (0,0) coefficient of `gmnc` is negative on every surface.

PB3D then reads the file as a VMEC file. No Fortran changes, and the file also works in other VMEC-based tools (BOOZ_XFORM, COBRAVMEC, TERPSICHORE).

Limits: only map 1 (map 3 is a straight cylinder). PB3D takes radial derivatives by cubic splines through the surfaces; $\lambda$ and the Jacobian go through VMEC's half mesh and are extrapolated over half a step at both ends; the poloidal flux carries the trapezoidal error of section 1. `ns` is free in the converter, so each of these errors falls when `ns` grows, at the cost of PB3D's memory and time, which grow with the number of surfaces in the equilibrium grid.

Effort (guess): 3 to 5 days for the converter and its checks.

### Route B: a GVEC reader in PB3D for map 1

PB3D fills the same `VMEC_vars` arrays from GVEC data directly: Fourier coefficients of $R$, $Z$, $\lambda$ with radial derivatives of order 0 to 3 taken from GVEC's B-splines, all on the full mesh, without half mesh or extrapolation. `eq_style` stays 1. Only `files_ops.f90` (which file is which) and `read_input_eq` (`input_ops.f90:1175-1178`) change, and one more entry in `misc_in_V` records where the arrays came from. A new `eq_style` value would be harmful: `select case (eq_style)` blocks without a default branch, such as the one in `calc_eq_1` (`eq_ops.f90:310-315`), would skip their work without a message. Two ways to get the data:

- B1: link GVEC's Fortran library into PB3D, read the state with `MODgvec_ReadState`, evaluate with GVEC's basis routines. This adds a compiled dependency to PB3D's CMake build and ties PB3D to GVEC's internal modules, which are not a published interface.
- B2: a Python script writes a small HDF5 file that holds exactly the `VMEC_vars` arrays on a grid equidistant in $s$; PB3D gets a reader of perhaps 200 lines. Because pyGVEC stops at second radial derivatives, the script must evaluate the third derivatives from GVEC's B-spline coefficients itself. The coefficients of odd $m$ behave like $\sqrt s$ near the axis, so their $s$-derivatives are unbounded at $s = 0$: B2 has to decide what to store on the axis surface (GVEC's own writers start at $s = 10^{-8}$, `gvec_to_castor3d`, or avoid the axis, `cas3d.py:119`).

Effort (guess): B2 one to two weeks, B1 more.

### Route C: GVEC's other maps (knot, Frenet frame, G-Frame)

For maps 10, 20 and 21, $X^1$ and $X^2$ are not $R$ and $Z$, and the position needs $h$. PB3D would need:

- the metric of the reference coordinates $q = (X^1, X^2, \zeta)$, $g_Q = (\partial x/\partial q)^T (\partial x/\partial q)$, which is not diagonal in general, with derivatives up to total order 2, so derivatives of $h$ up to order 3. For the G-Frame, $x = X_0(\zeta) + q^1 N(\zeta) + q^2 B(\zeta)$ is linear in $q^1$, $q^2$, and the $\zeta$ derivatives follow from the frame data in the netCDF file of `hmap_ncfile`. For the Frenet frame, $N$ and $B$ already contain up to third derivatives of the axis curve, so the metric needs up to its fifth derivative;
- a transformation from GVEC's logical coordinates to $q$ that replaces $T_{VC}$, and the Jacobian computed from it;
- the 3-D vacuum boundary and its normal in Cartesian coordinates from $h$ (`vac_ops.f90`);
- POST's plots in Cartesian coordinates from $h$;
- a normalization length for configurations without a major radius, for example the axis length divided by $2\pi$ (GVEC's `r_major`);
- the number of field periods of the frame, and a statement of what the toroidal mode number $n$ means when $\zeta$ is the curve parameter of the axis;
- the handedness: GVEC's frames are right-handed while PB3D's VMEC path assumes left-handed coordinates (section 1), so route C needs its own sign conventions in the places listed there.

This is the only route to configurations whose cross-sections cannot be written as $R(\theta, \varphi)$, $Z(\theta, \varphi)$, such as the figure-8 stellarator. A variant, A+, writes a G-Frame equilibrium into a wout file when every plane $\varphi$ = constant cuts each flux surface in one closed curve; it needs many Fourier modes and fails for the figure-8 shape.

Effort (guess): several weeks.

## 4. Recommendation

Route A first. It changes no Fortran, it gives the test cases of section 5 that any later route needs too, and it measures how large the spline, half-mesh and flux-integral errors of section 1 are at a given `ns` (decision D1). Route B only if D1 says so. Route C only if Toon wants configurations that need a flexible frame (decision D2); it is a separate project.

## 5. Stages of route A

Each stage ends with a gate. After each stage: update this plan with what was learned, comment the gate result on #104, commit and push, never merge. Heavy runs go to Snellius as parallel Slurm jobs.

### A0. GVEC installed on the Mac and on Snellius

`pip install gvec==1.5.1` in PB3D's uv environment as an optional dependency group `gvec` (as `vmec = ["vmecpp"]`). The package is a source distribution: it compiles GVEC and needs a Fortran compiler, CMake, LAPACK and netCDF-Fortran (Homebrew on the Mac, modules on Snellius). Gate: GVEC's elliptic-tokamak tutorial runs and `pygvec info` prints its state on both machines.

### A1. A GVEC twin of the circular fixture

`pygvec convert-params --vmec` cannot be used: GVEC 1.5.1 accepts only `power_series` profiles (`util/params.py:455-485`), while `--export_HEL` writes `cubic_spline` profiles (`eq_ops.f90:1949-1985`), and the RBS and ZBC columns of these inputs make it choose a stellarator-asymmetric run. Instead, `pygvec convert-wout tests/fixtures/wout_circ_eps033_q104.nc` gives a GVEC parameter file with the twin's boundary, the profiles as `interpolation` tables, and `_cos_`/`_sin_` symmetry. The parameter file has no `whichInitEquilibrium`, so GVEC starts from its default initial state and solves its own equilibrium, in stages that raise `sgrid_nElems` and the Fourier resolution (`--param` values are written into the file).

`convert-wout` writes the profile tables at 101 values of $s$ (`np.linspace(0, 1, 101)` in `convert_wout.py`). Resampling the twin's profiles there changes $q$ by $7 \cdot 10^{-7}$ but $\iota'$ by up to $2.2 \cdot 10^{-4}$ near $s = 0.59$, where the shear is only between $-0.0035$ and $-0.0018$: 6 to 12 % of the local shear. With 319 values the change in $\iota'$ is $1.9 \cdot 10^{-5}$, with 639 values $4.9 \cdot 10^{-6}$; for cbm18a cut, 101 values change $p'$ by $5 \cdot 10^{-4}$ and $\iota'$ by $1.8 \cdot 10^{-3}$ of their maxima. So the tables in the parameter file are replaced by `iotaf` and `presf` on all the twin's surfaces (a small script after `convert-wout`), or by the AI_AUX and AM_AUX knots with VMEC's end conditions. How GVEC interpolates such a table is checked by evaluating $\iota$, $\iota'$, $p$, $p'$ of the solved state against the twin's.

The twin's VMEC input prescribes $\iota$ (`NCURR = 0`) and the GVEC run prescribes the same table, so $q$ agrees by construction and is no check. Gate: the force residual decreases with resolution, and at equal $s$ the solved GVEC surfaces agree with the VMEC twin (axis position, $R$ and $Z$ on the surfaces $s$ = 0.25, 0.5, 1, $B_{00}$ on the axis) to within the VMEC twin's own resolution error, estimated from VMEC++ runs with 159 and 319 surfaces.

### A2. The converter

Write `gvec_to_wout.py` (route A). Three gates.

Gate A2a, identities inside the written file, independent of GVEC's VMEC reader: $d\,$`phi`$/ds$ = `phipf`; `chipf` = `iotaf` $\cdot$ `phipf`; $d\,$`chi`$/ds$ = `signgs` $\cdot$ `chipf`; `gmnc` negative and equal to $R\,(R_\theta Z_s - R_s Z_\theta)$ built from the file's own `rmnc` and `zmns`; `bmnc` equals $\lvert B \rvert$ built from the file's metric, `gmnc`, `lmns`, `iotaf` and `phipf`. The check script is first run on `wout_circ_eps033_q104.nc` from VMEC++, which fixes the factors $2\pi$ and the signs in the formulas the check uses. On the VMEC twin `gmnc` agrees with $R\,(R_\theta Z_s - R_s Z_\theta)$ to $3 \cdot 10^{-9}$ relative in the interior and $5 \cdot 10^{-7}$ at $s = 0.017$. On the converted file the radial derivatives in these identities are differences between surfaces, so outside the axis region (where the coefficients of odd $m$ behave like $\sqrt s$) the residuals must fall as $h^2$ when `ns` grows; a wrong sign or factor gives a residual of order one that does not fall.

Gate A2b, round trip: `wout_circ_eps033_q104.nc` $\to$ `pygvec convert-wout --param sgrid_nElems=...` $\to$ GVEC state, without a solve $\to$ `gvec_to_wout.py` $\to$ wout'. Only `rmnc`, `zmns` away from the axis, `iotaf`, `presf` and `phi` can converge to the original when `sgrid_nElems` grows. `lmns` keeps a floor (`InitVMEC` fits $\lambda$ on the half mesh, moves it to the full mesh and fits again, `vmec.F90:365-398`), every mode's axis value is a parabolic extrapolation (`vmec.F90:335`), and `gmnc`, `bmnc` and the $B$ components in wout' are GVEC's evaluations, which differ from VMEC's own finite-difference values by VMEC's discretisation error.

Gate A2c: LIBSTELL reads wout', `tests/regression/check_vmec_read.py` passes (mode lists and PB3D's copy of the file), and PB3D's fixed-boundary eigenvalue on wout' approaches that on the original wout (n = 10, m = 11, +7.84303e-2) as `sgrid_nElems` grows: three values of `sgrid_nElems`, the difference falling, and below $10^{-5}$ relative at the largest.

### A3. PB3D on the GVEC twin of the circular case

Convert the state of A1 and run the deck of `regression_circ_eps033_vmec_twin`. The VMEC twin agrees with HELENA to $1.1 \cdot 10^{-6}$ relative (+7.84303e-2 against +7.84304e-2) at this deck and resolution. Gate: the GVEC twin's fixed-boundary eigenvalue approaches the VMEC twin's with the same deck as GVEC's resolution grows, and lies within $10^{-5}$ relative of it at the highest resolution; the difference to HELENA is reported next to it. This is tighter than the twin test's tolerance of $10^{-4}$, and reachable only with the denser profile tables of A1. The table has two values of the converter's `ns` and three GVEC resolutions, so that it shows which of the two limits the agreement. The free-boundary eigenvalue is compared with the VMEC twin's -0.43901; this also tests the half-step extrapolation of the Jacobian to the boundary that the vacuum normal uses.

### A4. Finite pressure and 3-D

Three cases, run side by side on Snellius:

- cbm18a cut at the pedestal foot (axisymmetric, finite pressure; PB3D on HELENA gives -0.04012 at n = 10 fixed boundary, `tests/reference/cbm18a_cut/README.md`; VMEC twin from `input.cbm18a_cut` with 639 surfaces);
- the rippled twin, N = 5, $\delta$ = 0.02 (3-D; PB3D free boundary -0.4188 against TERPSICHORE's family -0.41903);
- W7-X (PB3D fixed-boundary reference 6.62833e-2 on Snellius).

For each: GVEC from `pygvec convert-wout` of the VMEC file, with the profile tables replaced as in A1, and its own solve; PB3D on the converted wout against PB3D on VMEC's wout, at two resolutions of each code. `convert-wout` always prescribes $\iota$. For a VMEC run with a prescribed current (W7-X, if its input has `NCURR` = 1), GVEC prescribes VMEC's converged `iotaf`; both codes then solve the same problem up to discretisation. GVEC's current constraint would need `I_tor` and `picard_current` added by hand and is not needed for this comparison. Gate: the difference between the two eigenvalues decreases when the resolutions grow; the numbers go on #104.

### Decision D1, after A3 and A4

Fix a target accuracy for the eigenvalue first (for example $10^{-5}$ relative). Then, at a fixed GVEC state:

- compare what PB3D stores in `PB3D_out.h5` with GVEC's exact values at the same points, for each radial derivative order of $R$, $Z$ and $\lambda$ separately (the values of $R$ and $Z$ agree by construction, those of $\lambda$ and the Jacobian carry the half-mesh spline error; the first and second derivatives carry the spline error; the third falls only in proportion to $h$), outside the axis region (the deck starts at `min_r_sol` = 0.2), with the exact third derivatives taken from the B-spline coefficients;
- run PB3D at several `ns` of the converted file and find the `ns` at which the eigenvalue stays within the target, and PB3D's memory and time there.

If an `ns` with acceptable cost reaches the target, route A is enough. If not, route B.

### A5. Regression test, documentation

- `regression_circ_eps033_gvec_twin`: the converted wout of A3 as a compressed fixture in `tests/fixtures`, so that CI needs no GVEC; the eigenvalues anchored as in the VMEC twin test.
- The identity check of A2a as a script in `tests/regression`, run on the fixture; the round trip of A2b only when the `gvec` group is installed, skipped with a note otherwise, as `check_vmec_read.py` is without h5py.
- README entry; a page on the site (PB3D.github.io) on how to run PB3D on a GVEC equilibrium: `pygvec run`, `gvec_to_wout.py`, PB3D. No version change, since no Fortran changes.

## 6. Stages of routes B and C, if chosen

Route B2: B2.1 the HDF5 writer in Python, with third radial derivatives from the B-spline coefficients and a decision on the axis surface; B2.2 the reader in PB3D and the provenance entry in `misc_in_V`; B2.3 the cases of A3 and A4 again. Gate: the eigenvalues of routes B2 and A agree to within the spline error measured in D1, and route B2 reaches the target accuracy with fewer surfaces.

Route C: C.1 the metric of the G-Frame and its derivatives in PB3D, tested on a G-Frame description of a configuration that map 1 also describes (GVEC's G-Frame tutorial builds one); C.2 the vacuum boundary, the normal and the plots from $h$; C.3 gate: PB3D gives the same eigenvalue for the same equilibrium computed with map 1 (route A) and with map 21 (route C), to within the difference of the two GVEC solutions, which falls with resolution; C.4 the figure-8 configuration of Plunk et al. (2025).

## 7. Decisions for Toon

- D0: where `gvec_to_wout.py` lives: (a) in PB3D's `scripts/`, (b) in GVEC as a `pygvec` subcommand next to `to-cas3d`, which Toon would propose upstream (MPCDF GitLab or the GitHub mirror), (c) first (a), later (b).
- D1: after A4, whether route B is needed (section 5).
- D2: whether route C is wanted.
- D3: whether the GVEC twin becomes a regression test (A5) or stays a separate study outside this repository.

## 8. Sources

- GVEC JOSS paper, doi:10.21105/joss.09670; flexible frame, doi:10.1088/1361-6587/adba11.
- GVEC documentation, https://gvec.readthedocs.io: parameter list, theory, coordinate conventions, quantities, pygvec CLI, CAS3D and other interfaces.
- GVEC 1.5.1 source, public mirror https://github.com/gvec-group/gvec: `src/functionals/mhd3d/mhd3d.F90`, `src/vmec/vmec.F90`, `python/gvec/scripts/convert_wout.py`, `python/gvec/scripts/cas3d.py`, `python/gvec/util/params.py`, pyGVEC's `_evaluate.py`.
- LIBSTELL (STELLOPT repository): `LIBSTELL/Sources/Modules/read_wout_mod.f90` (`read_wout_nc`) and `LIBSTELL/Sources/Ezcdf/ezcdf_inqvar.f90`.
