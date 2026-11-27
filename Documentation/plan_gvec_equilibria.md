# Plan: GVEC equilibria in PB3D (issue #104)

Status: plan, 27 November 2026. Nothing is implemented. The decisions in section 7 are Toon's.

In short: a Python converter writes a GVEC equilibrium as a VMEC wout file, and PB3D reads it without changes (route A). This covers GVEC equilibria in cylindrical coordinates. All of GVEC's other toroidal coordinate maps have one common form, a closed curve with two frame vectors, and PB3D's VMEC path assumes cylindrical coordinates in six routines only; replacing the cylindrical formulas there by those of this form opens PB3D to every GVEC map, with the same converter (route C).

## 1. What PB3D needs from an equilibrium

### 1.1 What PB3D reads

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
- Row $j$ of a half-mesh array ($j$ = 2 to `ns`) lies at $s_{j-1/2} = (j - 3/2)/(\mathrm{ns} - 1)$; `read_VMEC` does not use row 1.
- The Jacobian comes from `gmnc`, not from $R$ and $Z$ (`calc_jac_V_ind`). The vacuum normal uses `jac_FD`, which comes from `gmnc` extrapolated half a step to the boundary (`vac_ops.f90`, around line 430).
- The poloidal flux is the trapezoidal integral of `iotaf` times `phipf` (`VMEC_ops.f90:86-87`), with an error proportional to $h^2$, $h = 1/(\mathrm{ns}-1)$. With `use_pol_flux_F` it is the radial coordinate. No converter changes this.
- PB3D never reads `signgs`. It assumes VMEC's left-handed coordinates $(s, \theta, \varphi)$ with a negative Jacobian, for example in `calc_flux_q_VMEC` ("conversion VMEC LH -> PB3D RH", `eq_ops.f90:368`), in `calc_eq_2` (`pmone = -1`, line 522), in `grid_ops.f90:176` and 889, `eq_utilities.f90:786` and 866, `grid_utilities.f90:211` and 440.
- $R_0$ = `R_V_c(1,1,0)` (`calc_normalization_const_VMEC`) assumes that the mode (0,0) comes first in the main mode list, and `read_VMEC` stops if (0,0) is not the first Nyquist mode. `read_VMEC` checks `mnmax` $= (\mathrm{ntor}+1) + (2\,\mathrm{ntor}+1)(\mathrm{mpol}-1)$, VMEC's full rectangular mode set. The Nyquist list may have any size, in VMEC's order (`VMEC_ops.f90:161-166`).
- `PB3D_out.h5` stores the `VMEC_vars` arrays, and POST reads them back (`PB3D_ops.f90:162-309`); POST never opens the equilibrium file.
- For every point of the equilibrium grid, `PB3D_out.h5` also holds the point's VMEC coordinates (`r_E`, `theta_E`, `zeta_E`) and the quantities that enter the energy: `g_FD`, `h_FD`, `jac_FD`, `S`, `kappa_n`, `kappa_g`, `sigma` (`print_output_grid`, `print_output_eq_2`). A script can therefore compare them with an independent evaluation at the same points (gate A3a).

### 1.2 Where the VMEC path assumes cylindrical coordinates

`calc_eq_2` builds the metric in three steps: the metric $g_C$ of the cylindrical coordinates $(R, \varphi, Z)$; the matrix $T_{VC}$ of the derivatives of $(R, \varphi, Z)$ with respect to VMEC's $(s, \theta, \varphi)$, which gives $g_V = T_{VC}\, g_C\, T_{VC}^T$ (`calc_g`); and the matrix $T_{VF}$ to flux coordinates. `g_C` is stored as a full symmetric matrix and `calc_g` treats it as one, so nothing after `calc_g_C` assumes that it is diagonal. $T_{VC}$ holds only derivatives of $R$ and $Z$ and the entry $\partial\varphi/\partial\varphi = 1$ (`calc_T_VC_ind`).

The formulas of cylindrical coordinates themselves appear in these places of the VMEC path:

| routine | what it does with cylindrical coordinates |
|---|---|
| `calc_g_C_ind` (`eq_ops.f90:2772`) | $g_C = \mathrm{diag}(1, R^2, 1)$ and its derivatives |
| `calc_derived_DE_epar_VMEC` (`eq_ops.f90:4512`) | derivatives of the basis vectors $e_\theta$, $e_\varphi$ in the components $(e_R, e_\varphi, e_Z)$, with the terms $R_{\varphi\varphi} - R$ and $2R_\varphi/R$; `kappa_n`, `kappa_g` and `sigma` come from them |
| the VMEC block of the 3-D vacuum (`vac_ops.f90:420-460`) | position $(R\cos\varphi, R\sin\varphi, Z)$ and the Cartesian components of the normal on the boundary |
| `calc_XYZ_grid_VMEC` (`grid_utilities.f90:871`) | Cartesian positions for all plots |
| `calc_vec_comp`, step 5 (`grid_utilities.f90:2676-2698`) | Cartesian components of vectors for POST's plots |
| `calc_normalization_const_VMEC` (`eq_ops.f90:5358`) | $R_0$ from the (0,0) mode of $R$ on the axis, or from `rmin_surf` and `rmax_surf` |

Outside the release path: the debug tests `test_g_V` and `test_jac_V`, and the diagnostic `check_vac_exact_1` (sources on the magnetic axis). POST's vacuum plot (`vac_pot_plot`) uses a cylindrical box of points in space, which is valid for any equilibrium.

### 1.3 What LIBSTELL requires of a wout file

`files_ops.f90` treats a file whose name ends in `.nc` as VMEC (`eq_style` 1) and reads it with LIBSTELL's `read_wout_file`. Its netCDF reader (`read_wout_nc` in `read_wout_mod.f90`, with `Ezcdf/ezcdf_inqvar.f90`) reads variables by name:

- A missing variable is skipped silently (`cdfgv` and `cdfInqVar` return without an error), and a variable LIBSTELL does not know is never asked for. So the file needs only what PB3D uses and may hold more; the VMEC++ fixture already lacks `lrecon`, `lmove_axis` and `mnmaxpot`.
- `ier_flag` must be 0 (or the "more iterations" or "jac75" flags), otherwise `ierr = -2`.
- `raxis_cc` and `zaxis_cs` must have length `ntor`+1, otherwise LIBSTELL stops.
- Logicals are int32 variables with the suffix `__logical__` (`lasym__logical__`).
- Integers must be NC_INT and reals NC_DOUBLE. Another type prints a message, sets no error and leaves the value unset: an int64 `ns` from numpy would be lost.
- With `version_` of 9 or more LIBSTELL reads `currumnc`; below 9 it computes the currents from the $B$ components.
- `nextcur` = 0.

## 2. What GVEC provides and its conventions

### 2.1 Unknowns and coordinate maps

GVEC's unknowns are $X^1$, $X^2$ and $\lambda$ as functions of $(\rho, \vartheta, \zeta)$: B-splines of degree `X1X2_deg` and `LA_deg` in $\rho$ (continuity one less than the degree, the only choice GVEC allows) on `sgrid_nElems` elements, Fourier series in $\vartheta$ and $\zeta$. The position is $x = h(X^1, X^2, \zeta)$ with the map chosen by `which_hmap`. GVEC computes fixed-boundary equilibria only; PB3D's free-boundary stability uses a fixed-boundary equilibrium with a vacuum perturbation outside it, so this is no restriction.

Every toroidal map of GVEC has the form

$$x = X_0(\zeta) + X^1\, N(\zeta) + X^2\, B(\zeta),$$

a closed curve $X_0$ and two vectors $N$ and $B$ that are periodic in $\zeta$ (read from the headers of `src/hmap/hmap_*.F90` and the theory page of the documentation). $N$ and $B$ need not be perpendicular to each other or of unit length.

| `which_hmap` | $X_0(\zeta)$ | $N(\zeta)$ | $B(\zeta)$ |
|---|---|---|---|
| 1, cylindrical | 0 | $(\cos\zeta, -\sin\zeta, 0)$ | $(0, 0, 1)$ |
| 10, knot $(k, \ell)$ | $\big((R_0 + \delta\cos\ell\zeta)\cos k\zeta,\ -(R_0 + \delta\cos\ell\zeta)\sin k\zeta,\ \delta\sin\ell\zeta\big)$ | $(\cos k\zeta, -\sin k\zeta, 0)$ | $(0, 0, 1)$ |
| 20, Frenet | the axis curve, from Fourier coefficients of its $R$ and $Z$ | the Frenet normal times a sign $\sigma(\zeta)$ | the Frenet binormal times $\sigma(\zeta)$ |
| 21, G-Frame | `axis/xyz` of the netCDF file `hmap_ncfile` | `axis/Nxyz` | `axis/Bxyz` |

Map 3 is a straight cylinder and not toroidal. For map 1, $X^1 = R$, $X^2 = Z$, and $\zeta$ runs clockwise seen from above. $(\rho, \vartheta, \zeta)$ is right-handed for every map.

pyGVEC returns the position `pos` and the vectors `e_q1` $= \partial x/\partial X^1 = N$ and `e_q2` $= B$ for any map, so a script can obtain $X_0 = x - X^1 N - X^2 B$, $N$ and $B$ on a $\zeta$ grid without knowing which map the equilibrium uses, and can check that the three do not depend on $\rho$ and $\vartheta$.

### 2.2 Conventions, read from GVEC 1.5.1's source

- GVEC's $\Phi$ is the toroidal flux divided by $2\pi$. The input `PhiEdge` is the total toroidal flux; `mhd3d.F90:137-138` divides it by $2\pi$ ("normalization like in VMEC"), and for a run from parameters $\Phi = \Phi_\mathrm{edge}\, s$ with $s = \rho^2$ (`mhd3d.F90:182`). GVEC's CAS3D writer multiplies $\Phi$ and $\chi$ by $2\pi$ (`python/gvec/scripts/cas3d.py:260-267`).
- GVEC's $\mu_0 = 4\pi \cdot 10^{-7}$ (`mhd3d.F90:117`), and the energy has $\lvert B \rvert^2/(2\mu_0)$ and $p$, so $p$ is in Pa, as `presf` in a wout file.
- The field is $B = \Phi'\,\nabla\rho\times\nabla(\vartheta + \lambda) - \chi'\,\nabla\rho\times\nabla\zeta$ with $\iota = \chi'/\Phi'$, so $B^\zeta = \Phi'(1 + \partial_\vartheta\lambda)/\mathcal J$.
- GVEC's VMEC reader (`InitVMEC`, `src/vmec/vmec.F90:74-124`) always sets $\zeta = -\varphi$. For VMEC's usual `signgs` = -1 it negates the sign of `gmnc` and of `iotaf` and maps the modes; for `signgs` = +1 it also reverses $\theta$. `pygvec convert-wout` passes `PhiEdge` = $-$`phi`(edge) (`scripts/convert_wout.py`). The CAS3D writer writes $-\iota$.

### 2.3 From a GVEC state to VMEC's conventions

With $s = \rho^2$, $\theta = \vartheta$ and $\varphi = -\zeta$, the coordinates $(s, \theta, \varphi)$ are left-handed with a negative Jacobian, as PB3D assumes (section 1.1), for every map of section 2.1. The reference coordinates $(X^1, \varphi, X^2)$ are then right-handed, as $(R, \varphi, Z)$ is.

| wout quantity (`signgs` = -1) | from GVEC |
|---|---|
| `rmnc`, `zmns`, `lmns`: series of $R$, $Z$, $\lambda$ in $(\theta, \varphi)$ | $X^1$, $X^2$, $\lambda$ at $(\sqrt s, \theta, -\varphi)$ |
| `phi` | $-2\pi\Phi$ |
| `phipf` | $d\,$`phi`$/ds$ |
| `iotaf` | $-\iota$ |
| `chi` | $-2\pi\chi$ (VMEC writes `chi` $=$ `signgs` $\int$ `chipf` $ds$, while `chipf` $=$ `iotaf` $\cdot$ `phipf` $= 2\pi\, d\chi/ds$) |
| `presf` | $p$ |
| `gmnc` ($\sqrt g$ in $(s, \theta, \varphi)$) | $-\mathcal J/(2\rho)$ with GVEC's full Jacobian $\mathcal J$, negative on every surface |
| $B_s$, $B_u$, $B_v$ | $B_\rho/(2\rho)$, $B_\vartheta$, $-B_\zeta$ |
| $B^u$, $B^v$ | $B^\vartheta$, $-B^\zeta$ |

The signs of `phi`, `iotaf` and `chi` were checked on the circular VMEC twin (`wout_circ_eps033_q104.nc`, VMEC++, `signgs` = -1, 319 surfaces): `phi`(edge) = -0.32945, `iotaf` from -0.9605 to -0.9632, $d\,$`phi`$/ds$ = `phipf` and `chipf` = `iotaf` $\cdot$ `phipf` to rounding, `chi`(edge) = -0.3168 with $d\,$`chi`$/ds = -$`chipf`, and `bsupvmnc`(0,0) = -1.08 where GVEC's $B^\zeta$ is positive. GVEC's own VMEC reader divides `chi` by $2\pi$ without the factor `signgs` (`vmec_readin.F90:638`), so the inverse of that reader would give the wrong sign of `chi`. PB3D does not read `chi`, but other tools do. Gates A2a and A2b check all of this again on the converted files.

### 2.4 What pyGVEC evaluates

pyGVEC (`gvec.load_state`, `state.evaluate`) evaluates at given $(\rho, \vartheta, \zeta)$: $X^1$, $X^2$, $\lambda$ with radial derivatives of order 0 to 2 (`_evaluate.py:73`, the pattern `r{0,2}`); the position, the basis vectors, the metric and the Jacobian; the vectors $B$ and $J$ in Cartesian components, $\lvert B \rvert$, the curvature vector of the field lines (`kappa_B`) and its normal and geodesic parts; $\Phi$, $\chi$, $\iota$, $p$ and their radial derivatives. GVEC's CAS3D writer shows the pattern a wout writer can follow: evaluate on a grid, change the angle convention, Fourier transform.

GVEC does not write wout files: not in version 1.5.1 and not in the development version (1.5.2.dev88 on 27 November 2026, same list of `pygvec` subcommands). A search of GVEC's public issue list for "wout" returns nothing; the nearest entry is issue 22 (open), "Frenet → RZφ conversion", about bringing a solution in the Frenet frame to cylindrical coordinates. I found no other tool that writes a wout file from a GVEC state.

## 3. Three routes

### Route A: a converter from a GVEC state to a wout file, PB3D unchanged

A Python script `gvec_to_wout.py` loads a GVEC state with pyGVEC and writes a wout file on `ns` surfaces equidistant in $s$, with the transformations of section 2.3:

- $R$, $Z$ on the full mesh and $\lambda$ on the half mesh (rows 2 to `ns`, section 1.1), from a 2-D FFT of pyGVEC's values on a $(\theta, \varphi)$ grid into VMEC's basis $\cos(m\theta - n\varphi)$, $\sin(m\theta - n\varphi)$, in VMEC's mode order with (0,0) first. The grid is fine enough that the FFT returns the finite Fourier series exactly, so GVEC's internal mode order never enters. Modes outside GVEC's set are zero.
- `gmnc` on the half mesh on a Nyquist list large enough to hold the Jacobian. For map 1, $\sqrt g = R\,(R_\theta Z_s - R_s Z_\theta)$ is a product of three finite Fourier series, so it is itself a finite series with $m \le 3\,m_\mathrm{max}$, $\lvert n \rvert \le 3\,n_\mathrm{max}$ (with different limits $M_1$, $M_2$ for $X^1$ and $X^2$: $2M_1 + M_2$). VMEC truncates it (the circular fixture's Nyquist list stops at $m = 15$ for `mpol` = 12, where the exact Jacobian reaches $m = 33$); the converter need not. For `mpol` = `ntor` = 12 the full list has about 2,400 modes, and in PB3D `fourier2real` costs time in proportion to the number of modes at every grid point and the trigonometric table grows with the largest $m$ and $n$. So the converter writes the smallest full rectangular Nyquist set outside which all coefficients stay below a tolerance (for example $10^{-12}$ of the (0,0) coefficient) on all surfaces. A rectangular set keeps the file valid for other tools.
- `bmnc`, $B_s$, $B_u$, $B_v$, $B^u$, $B^v$ on the same list. These are not finite series, so the list truncates them; PB3D uses `bmnc` only for $B_0$ and the $B$ components only in debug tests.
- the profiles and scalars of section 1.1, `raxis_cc`, `zaxis_cs` of length `ntor`+1, `ier_flag` = 0, `lfreeb` = F, `lrfp` = F, `lasym` from GVEC's sine and cosine choice, `signgs` = -1, `nextcur` = 0. The file is written with netCDF4 (in PB3D's `scripts` dependency group, as `tests/regression/reverse_field.py` uses it), with the names, dimensions and types of the VMEC++ fixture, which LIBSTELL is known to read (section 1.3).
- before writing, the converter checks that the (0,0) coefficient of `gmnc` is negative on every surface, and that the equilibrium uses the cylindrical map: it obtains $X_0$, $N$, $B$ as in section 2.1 and stops with a message unless $X_0 = 0$, $N = (\cos\zeta, -\sin\zeta, 0)$ and $B = (0, 0, 1)$.

PB3D then reads the file as a VMEC file. No Fortran changes, and the file also works in other VMEC-based tools (BOOZ_XFORM, COBRAVMEC, TERPSICHORE).

Limits: only map 1. PB3D takes radial derivatives by cubic splines through the surfaces; $\lambda$ and the Jacobian go through VMEC's half mesh and are extrapolated over half a step at both ends; the poloidal flux carries the trapezoidal error of section 1.1. `ns` is free in the converter, so each of these errors falls when `ns` grows, at the cost of PB3D's memory and time, which grow with the number of surfaces in the equilibrium grid.

Effort (guess): 3 to 5 days for the converter and its checks, 1 to 2 days for the comparison script of gate A3a.

### Route B: a GVEC reader in PB3D for map 1

PB3D fills the same `VMEC_vars` arrays from GVEC data directly: Fourier coefficients of $R$, $Z$, $\lambda$ with radial derivatives of order 0 to 3 taken from GVEC's B-splines, all on the full mesh, without half mesh or extrapolation. `eq_style` stays 1. Only `files_ops.f90` (which file is which) and `read_input_eq` (`input_ops.f90:1175-1178`) change, and one more entry in `misc_in_V` records where the arrays came from. A new `eq_style` value would be harmful: `select case (eq_style)` blocks without a default branch, such as the one in `calc_eq_1` (`eq_ops.f90:310-315`), would skip their work without a message. Two ways to get the data:

- B1: link GVEC's Fortran library into PB3D, read the state with `MODgvec_ReadState`, evaluate with GVEC's basis routines. This adds a compiled dependency to PB3D's CMake build and ties PB3D to GVEC's internal modules, which are not a published interface.
- B2: the converter of route A writes, next to the wout variables, the same coefficients on the full mesh with their radial derivatives of order 1 to 3 as extra variables (LIBSTELL skips them, section 1.3), and `read_VMEC` uses them instead of its splines when they are present, reading them with LIBSTELL's `ezcdf` routines, which PB3D's LIBSTELL build already holds. One file then serves PB3D and the other tools. Because pyGVEC stops at radial derivatives of order 2, the script must evaluate the third derivatives from GVEC's B-spline coefficients itself. The coefficients of odd $m$ behave like $\sqrt s$ near the axis, so their $s$-derivatives are unbounded at $s = 0$: B2 has to decide what to store on the axis surface (GVEC's own writers start at $s = 10^{-8}$, `gvec_to_castor3d`, or avoid the axis, `cas3d.py:119`).

Effort (guess): B2 one to two weeks, B1 more.

### Route C: GVEC's other maps (knot, Frenet frame, G-Frame)

Because every map has the form of section 2.1, PB3D needs one generalisation, from the cylindrical map to

$$x = X_0(\varphi) + R\, N(\varphi) + Z\, B(\varphi),$$

where $R$ and $Z$ now stand for $X^1$ and $X^2$, and the cylindrical map is the case $X_0 = 0$, $N = e_R(\varphi)$, $B = e_Z$. With $\varphi = -\zeta$ the handedness that PB3D assumes holds for every map (section 2.3), so no sign convention changes.

What stays as it is: the file's `rmnc`, `zmns`, `lmns` (holding $X^1$, $X^2$, $\lambda$), the profiles, `gmnc` (the converter writes GVEC's full Jacobian, which includes the Jacobian of the map), $T_{VC}$, `calc_g` and everything after it.

What changes in PB3D is the list of section 1.2:

- `calc_g_C_ind`: the metric of $(R, \varphi, Z)$ becomes, with a prime for $d/d\varphi$ and $t = \partial x/\partial\varphi$ at fixed $R$ and $Z$, $t = X_0' + R\,N' + Z\,B'$,

$$g_C = \begin{pmatrix} N\cdot N & N\cdot t & N\cdot B \\ \cdot & t\cdot t & t\cdot B \\ \cdot & \cdot & B\cdot B \end{pmatrix},$$

  a polynomial of degree 2 in $R$ and $Z$ with coefficients that depend on $\varphi$. Its derivatives up to total order 2 need $X_0$, $N$, $B$ up to their third derivative. For the cylindrical map $N' = e_\varphi/R$ has unit length, $B' = 0$, and $g_C = \mathrm{diag}(1, R^2, 1)$.
- `calc_derived_DE_epar_VMEC`: the derivatives of the basis vectors of VMEC's coordinates, $e_\theta = R_\theta N + Z_\theta B$ and $e_\varphi = t + R_\varphi N + Z_\varphi B$, in Cartesian components. For example $\partial_\varphi e_\varphi = X_0'' + R\,N'' + Z\,B'' + 2R_\varphi N' + 2Z_\varphi B' + R_{\varphi\varphi} N + Z_{\varphi\varphi} B$, which for the cylindrical map ($N'' = -N$) gives the present terms with $R_{\varphi\varphi} - R$ and $2R_\varphi$.
- the vacuum block, `calc_XYZ_grid_VMEC` and step 5 of `calc_vec_comp`: positions and vector components from $X_0$, $N$, $B$ instead of $\cos\varphi$ and $\sin\varphi$.
- the normalization length: the length of the magnetic axis divided by $2\pi$ (GVEC's `r_major`), which the converter writes.
- the frame itself: the converter adds the nine Cartesian components of $X_0$, $N$, $B$ as Fourier series in $\varphi$ to the file (obtained as in section 2.1, for any map), PB3D reads them with LIBSTELL's `ezcdf` routines, evaluates them and their first three derivatives at the grid's $\varphi$ values, stores them in `PB3D_out.h5`, and POST reads them back. A file without them is a VMEC file and gets the cylindrical frame.
- the debug tests `test_g_V`, `test_jac_V` and the diagnostic `check_vac_exact_1`.

The converter of route A needs one change: it writes the frame instead of stopping when the map is not cylindrical. Such a file is not a VMEC file for other tools.

Two points that are specific to this route:

- The toroidal mode number. PB3D keeps one toroidal mode number $n$ in its flux coordinates, whose toroidal angle is the file's $\varphi$ up to the sign. For a G-Frame, $\zeta$ is the parameter of the curve, which GVEC's documentation recommends to be a Boozer angle, not the geometric toroidal angle. The same 3-D equilibrium described with two different toroidal angles gives PB3D two different sets of displacements for the same $n$, so the eigenvalue can depend on the parametrisation. This follows from the form of PB3D's Fourier series; how large the effect is has not been computed (stage C4).
- A reference. For shapes that cylindrical coordinates cannot describe there is no VMEC equilibrium to compare with. The check that remains valid is the pointwise comparison of PB3D's stored quantities with GVEC's own values (gate A3a), which does not depend on the map.

Effort (guess): two to three weeks, most of it tests.

## 4. Recommendation

Route A first. It changes no Fortran, it gives the checks and test cases of section 5 that any later route needs too, and it measures how large the spline, half-mesh and flux-integral errors of section 1.1 are at a given `ns` (decision D1). Route B only if D1 says so. Route C if Toon wants configurations that need a flexible frame (decision D2): it builds on the converter and the checks of route A, changes six release routines, and can be tested without GVEC (stage C2).

## 5. Stages of route A

Each stage ends with a gate. After each stage: update this plan with what was learned, comment the gate result on #104, commit and push, never merge. Heavy runs go to Snellius as parallel Slurm jobs.

### A0. GVEC installed on the Mac and on Snellius

`pip install gvec==1.5.1` in PB3D's uv environment as an optional dependency group `gvec` (as `vmec = ["vmecpp"]`). The package is a source distribution: it compiles GVEC and needs a Fortran compiler, CMake, LAPACK and netCDF-Fortran (Homebrew on the Mac, modules on Snellius). Gate: GVEC's elliptic-tokamak tutorial runs and `pygvec info` prints its state on both machines.

### A1. A GVEC twin of the circular fixture

`pygvec convert-params --vmec` cannot be used: GVEC 1.5.1 accepts only `power_series` profiles (`util/params.py:455-485`), while `--export_HEL` writes `cubic_spline` profiles (`eq_ops.f90:1949-1985`), and the RBS and ZBC columns of these inputs make it choose a stellarator-asymmetric run. Instead, `pygvec convert-wout tests/fixtures/wout_circ_eps033_q104.nc` gives a GVEC parameter file with the twin's boundary, the profiles as `interpolation` tables, and `_cos_`/`_sin_` symmetry. The parameter file has no `whichInitEquilibrium`, so GVEC starts from its default initial state and solves its own equilibrium, in stages that raise `sgrid_nElems` and the Fourier resolution (`--param` values are written into the file).

`convert-wout` writes the profile tables at 101 values of $s$ (`np.linspace(0, 1, 101)` in `convert_wout.py`). Resampling the twin's profiles there changes $q$ by $7 \cdot 10^{-7}$ but $\iota'$ by up to $2.2 \cdot 10^{-4}$ near $s = 0.59$, where the shear is only between $-0.0035$ and $-0.0018$: 6 to 12 % of the local shear. With 319 values the change in $\iota'$ is $1.9 \cdot 10^{-5}$, with 639 values $4.9 \cdot 10^{-6}$; for cbm18a cut, 101 values change $p'$ by $5 \cdot 10^{-4}$ and $\iota'$ by $1.8 \cdot 10^{-3}$ of their maxima. So the tables in the parameter file are replaced by `iotaf` and `presf` on all the twin's surfaces (a small script after `convert-wout`). How GVEC interpolates such a table is checked by evaluating $\iota$, $\iota'$, $p$, $p'$ of the solved state against the twin's.

The twin's VMEC input prescribes $\iota$ (`NCURR = 0`) and the GVEC run prescribes the same table, so $q$ agrees by construction and is no check. Gate: the force residual decreases with resolution, and at equal $s$ the solved GVEC surfaces agree with the VMEC twin (axis position, $R$ and $Z$ on the surfaces $s$ = 0.25, 0.5, 1, $B_{00}$ on the axis) to within the VMEC twin's own resolution error, estimated from VMEC++ runs with 159 and 319 surfaces.

### A2. The converter

Write `gvec_to_wout.py` (route A). Four gates.

Gate A2a, the file against GVEC in space, with no convention in between. At a few thousand points $(s, \theta, \varphi)$, on full-mesh surfaces and at random angles, a check script computes from the file's variables that PB3D reads (`rmnc`, `zmns`, `lmns`, `gmnc`, `phipf`, `iotaf`, `presf`):

- the position $(R\cos\varphi, R\sin\varphi, Z)$;
- the Cartesian vector $B = B^u\, \partial_\theta x + B^v\, \partial_\varphi x$, with $B^u$ and $B^v$ from VMEC's formulas in `phipf`, `iotaf`, $\lambda$ and $\sqrt g$, where $\lambda$ and $\sqrt g$ come to the full mesh by the cubic spline that `read_VMEC` uses;
- the pressure.

It compares them with pyGVEC's `pos`, `B` and `p` at $(\sqrt s, \theta, -\varphi)$. The positions and the pressure must agree to rounding, and $B$ to the tolerance of the Nyquist cut plus the error of the half-mesh spline, which falls when `ns` grows. A vector in space has no sign convention, so this catches a mirrored equilibrium, a wrong sign of a flux or of $\lambda$, and a missing $2\pi$, none of which an eigenvalue comparison is certain to show (the mirror image of an equilibrium has the same eigenvalues). The script's formulas for $B^u$ and $B^v$ are first calibrated on the VMEC++ fixture, whose `bsupumnc` and `bsupvmnc` they must reproduce.

Gate A2b, identities inside the written file: $d\,$`phi`$/ds$ = `phipf`; `chipf` = `iotaf` $\cdot$ `phipf`; $d\,$`chi`$/ds$ = `signgs` $\cdot$ `chipf`; `gmnc` negative and equal to $R\,(R_\theta Z_s - R_s Z_\theta)$ built from the file's own `rmnc` and `zmns`; `bmnc`, `bsup*` and `bsub*` equal to what the file's other variables give. On the VMEC twin `gmnc` agrees with $R\,(R_\theta Z_s - R_s Z_\theta)$ to $3 \cdot 10^{-9}$ relative in the interior and $5 \cdot 10^{-7}$ at $s = 0.017$. On the converted file the radial derivatives in these identities are differences between surfaces, so outside the axis region (where the coefficients of odd $m$ behave like $\sqrt s$) the residuals must fall as $h^2$ when `ns` grows. This gate covers the variables that PB3D does not read and other tools do.

Gate A2c, round trip: `wout_circ_eps033_q104.nc` $\to$ `pygvec convert-wout --param sgrid_nElems=...` $\to$ GVEC state, without a solve $\to$ `gvec_to_wout.py` $\to$ wout'. Only `rmnc`, `zmns` away from the axis, `iotaf`, `presf` and `phi` can converge to the original when `sgrid_nElems` grows. `lmns` keeps a floor (`InitVMEC` fits $\lambda$ on the half mesh, moves it to the full mesh and fits again, `vmec.F90:365-398`), every mode's axis value is a parabolic extrapolation (`vmec.F90:335`), and `gmnc`, `bmnc` and the $B$ components in wout' are GVEC's evaluations, which differ from VMEC's own finite-difference values by VMEC's discretisation error.

Gate A2d: LIBSTELL reads wout', `tests/regression/check_vmec_read.py` passes (mode lists and PB3D's copy of the file), and PB3D's fixed-boundary eigenvalue on wout' approaches that on the original wout (n = 10, m = 11, +7.84303e-2) as `sgrid_nElems` grows: three values of `sgrid_nElems`, the difference falling, and below $10^{-5}$ relative at the largest.

### A3. PB3D on the GVEC twin of the circular case

Convert the state of A1 and run the deck of `regression_circ_eps033_vmec_twin`.

Gate A3a, PB3D's stored equilibrium against GVEC, point by point. A script reads from `PB3D_out.h5` the VMEC coordinates of every point of the equilibrium grid and PB3D's values there (section 1.1), and evaluates the same quantities with pyGVEC at $(\sqrt{s}, \theta, -\varphi)$: $\lvert B \rvert$ (in PB3D $\lvert B \rvert^2 = g_{\parallel\parallel}/\mathcal J^2$ from `g_FD` and `jac_FD`, because $B = e_\parallel/\mathcal J$ in its flux coordinates), the normal and the geodesic curvature from GVEC's `kappa_B` and $\nabla\rho$ in PB3D's definitions, the parallel current `sigma` from GVEC's $J \cdot B/\lvert B \rvert^2$, and the shear `S` from pyGVEC's first derivatives of the metric. This tests the whole chain from the file to the quantities in the energy (`read_VMEC`, the splines, `calc_eq_2`, `calc_derived_q`) against values that no VMEC file can give, since VMEC has no exact evaluation between its surfaces. `sigma` changes sign under a mirror reflection, so the gate also sees a mirrored reading. Gate: the differences fall when `ns` grows, at the rate of the spline errors of section 1.1, with no difference that stays.

Gate A3b, eigenvalues. The VMEC twin agrees with HELENA to $1.1 \cdot 10^{-6}$ relative (+7.84303e-2 against +7.84304e-2) at this deck and resolution. The GVEC twin's fixed-boundary eigenvalue approaches the VMEC twin's with the same deck as GVEC's resolution grows, and lies within $10^{-5}$ relative of it at the highest resolution; the difference to HELENA is reported next to it. This is tighter than the twin test's tolerance of $10^{-4}$, and reachable only with the denser profile tables of A1. The table has two values of the converter's `ns` and three GVEC resolutions, so that it shows which of the two limits the agreement. The free-boundary eigenvalue lies within $10^{-4}$ relative of the VMEC twin's -0.43901; this also tests the half-step extrapolation of the Jacobian to the boundary that the vacuum normal uses.

### A4. Finite pressure and 3-D

Three cases, run side by side on Snellius, each with gate A3a first and then the eigenvalues:

- cbm18a cut at the pedestal foot (axisymmetric, finite pressure; PB3D on HELENA gives -0.04012 at n = 10 fixed boundary, `tests/reference/cbm18a_cut/README.md`; VMEC twin from `input.cbm18a_cut` with 639 surfaces);
- the rippled twin, N = 5, $\delta$ = 0.02 (3-D; PB3D free boundary -0.4188 against TERPSICHORE's family -0.41903);
- W7-X (PB3D fixed-boundary reference 6.62833e-2 on Snellius). GVEC's repository also holds a W7-X wout file with the GVEC parameters its developers use on it (`test-CI/examples/w7x_from_vmec_*`, `wout_d23p4_tm.nc`), a starting point for the GVEC resolution.

For each: GVEC from `pygvec convert-wout` of the VMEC file, with the profile tables replaced as in A1, and its own solve; PB3D on the converted wout against PB3D on VMEC's wout, at two resolutions of each code. `convert-wout` always prescribes $\iota$. For a VMEC run with a prescribed current (W7-X, if its input has `NCURR` = 1), GVEC prescribes VMEC's converged `iotaf`; both codes then solve the same problem up to discretisation. GVEC's current constraint would need `I_tor` and `picard_current` added by hand and is not needed for this comparison. Gate: the difference between the two eigenvalues decreases when the resolutions grow; the numbers go on #104.

### Decision D1, after A3 and A4

Fix a target accuracy for the eigenvalue first (for example $10^{-5}$ relative). Then, at a fixed GVEC state:

- take the pointwise differences of gate A3a at several `ns`, outside the axis region (the deck starts at `min_r_sol` = 0.2): they show which of PB3D's stored quantities carries the largest error at a given `ns` and how fast it falls;
- run PB3D at the same `ns` values and find the `ns` at which the eigenvalue stays within the target, and PB3D's memory and time there.

If an `ns` with acceptable cost reaches the target, route A is enough. If not, route B.

### A5. Regression test, documentation

- `regression_circ_eps033_gvec_twin`: the converted wout of A3 as a compressed fixture in `tests/fixtures`, so that CI needs no GVEC; the eigenvalues anchored as in the VMEC twin test.
- The identity check of A2b as a script in `tests/regression`, run on the fixture; the checks against GVEC (A2a, A2c, A3a) only when the `gvec` group is installed, skipped with a note otherwise, as `check_vmec_read.py` is without h5py.
- README entry; a page on the site (PB3D.github.io) on how to run PB3D on a GVEC equilibrium: `pygvec run`, `gvec_to_wout.py`, PB3D. No version change, since no Fortran changes.

## 6. Stages of routes B and C, if chosen

### Route B2

B2.1 the converter writes the extra full-mesh variables, with third radial derivatives from the B-spline coefficients and a decision on the axis surface. B2.2 `read_VMEC` uses them when present, and `misc_in_V` records it. B2.3 the cases of A3 and A4 again. Gate: the pointwise differences of A3a fall to the level of rounding and Fourier truncation at any `ns`, and the eigenvalues of routes B2 and A agree to within the error that D1 measured for route A.

### Route C

C1. The frame in PB3D: reading, evaluation with derivatives, storage in `PB3D_out.h5`, and the six routines of section 1.2 written for a general frame, with the cylindrical frame for files that hold none. Gate: every regression and full-stack test gives its present result. Whether VMEC files then go through the general formulas or keep the present cylindrical ones is decision D4: with the general formulas there is one code path, and the results may change at the level of rounding.

C2. Frame invariance, without GVEC. A script rewrites a VMEC wout file in another frame that keeps $\varphi$: a shifted one, $X_0 = (a_0 + a_1\cos N_\mathrm{fp}\varphi)\, e_R(\varphi)$, and a tilted one, $N$ and $B$ turned in the plane $\varphi$ = constant by an angle $\beta_1 \sin N_\mathrm{fp}\varphi$. The surfaces, the angles $(s, \theta, \varphi)$, $\lambda$ and the Jacobian are the same; only $X^1$, $X^2$ and the frame differ. Gate: on the circular twin, the rippled twin and W7-X, PB3D gives the same fixed- and free-boundary eigenvalues, the same stored quantities of A3a, and the same Cartesian positions as with the plain file, to rounding for the shifted frame and to the Fourier truncation of $X^1$, $X^2$ for the tilted one. This tests all six routines on existing cases.

C3. GVEC, the same boundary with map 1 and with a G-Frame whose $\zeta$ is the cylindrical angle. `pygvec load-quasr -f` builds a G-Frame from a boundary given as points in space, and keeps the toroidal angle of the data (GVEC's QUASR page: for surface data computed from cylinder coordinates "the construction will not change this"); the boundary of the rippled twin, sampled at cylindrical angles, is such data. Gate: A3a on both, and equal eigenvalues to within the difference of the two GVEC solutions, which falls with resolution.

C4. A G-Frame whose $\zeta$ is not the cylindrical angle: the G-Frame examples of GVEC's repository (`test-CI/examples/axisNB_*`, among them a knot and a quasi-helical case), a QUASR configuration (`pygvec load-quasr`), then the figure-8 configuration of Plunk et al. (2025). Gate: A3a. The eigenvalues are reported with a study of how they change with $n$ and with the parametrisation of the frame, the open point of route C.

## 7. Decisions for Toon

- D0: where `gvec_to_wout.py` lives: (a) in PB3D's `scripts/`, (b) in GVEC as a `pygvec` subcommand next to `to-cas3d`, which Toon would propose upstream (MPCDF GitLab or the GitHub mirror), (c) first (a), later (b).
- D1: after A4, whether route B is needed (section 5).
- D2: whether route C is wanted.
- D3: whether the GVEC twin becomes a regression test (A5) or stays a separate study outside this repository.
- D4, only with route C: whether VMEC files use the general frame formulas (one code path, results may change at rounding level) or keep the cylindrical formulas.

## 8. Sources

- GVEC JOSS paper, doi:10.21105/joss.09670; flexible frame, doi:10.1088/1361-6587/adba11.
- GVEC documentation, https://gvec.readthedocs.io: parameter list, theory, coordinate conventions, quantities, pygvec CLI (also the `develop` version), QUASR and G-Frame, CAS3D and other interfaces.
- GVEC 1.5.1 source, public mirror https://github.com/gvec-group/gvec: `src/functionals/mhd3d/mhd3d.F90`, `src/vmec/vmec.F90`, `src/vmec/vmec_readin.F90`, `src/hmap/hmap_RZ.F90`, `hmap_knot.F90`, `hmap_frenet.F90`, `hmap_axisNB.F90`, `python/gvec/quantities.py`, `python/gvec/scripts/convert_wout.py`, `python/gvec/scripts/cas3d.py`, `python/gvec/util/params.py`, pyGVEC's `_evaluate.py`; test data under `test-CI/`.
- LIBSTELL (STELLOPT repository): `LIBSTELL/Sources/Modules/read_wout_mod.f90` (`read_wout_nc`) and `LIBSTELL/Sources/Ezcdf/ezcdf_inqvar.f90`.
