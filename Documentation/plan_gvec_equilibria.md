# Plan: GVEC equilibria in PB3D (issue #104)

Status: plan, 27 November 2026. Nothing is implemented. The decisions in section 9 are Toon's.

In short:

- A GVEC run leaves a parameter file and a state file. The state file holds, for every Fourier mode of $X^1$, $X^2$ and $\lambda$, the coefficients of a B-spline in the radius (section 2.5). PB3D needs Fourier coefficients on flux surfaces with their radial derivatives; its energy reads them up to the second order (sections 1.1 and 1.5).
- Section 3 describes five routes from the one to the other. For each route it gives the path a GVEC equilibrium takes until PB3D has its eigenvalues, the place and size of every numerical error on that path, what changes in PB3D, and what speaks for and against the route. Section 4 compares the routes.
- The routes differ most in the radial derivatives. A wout file (routes A and C) hands PB3D values on surfaces, and PB3D takes the derivatives with cubic splines, as it does for VMEC files today. The routes B1, B2 and G take them from GVEC's B-splines, with errors five to seven orders of magnitude smaller (section 4.1).
- What is known of PB3D's sensitivity to errors in the equilibrium rules out no route (section 4.3). How accurate the second radial derivatives must be has not been measured. Stage S3 measures it on an existing VMEC case, without GVEC, and its result weighs on the choice of the route.

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
- The poloidal flux is the trapezoidal integral of `iotaf` times `phipf` (`VMEC_ops.f90:86-87`), with an error proportional to $h^2$, $h = 1/(\mathrm{ns}-1)$. With `use_pol_flux_F` it is the radial coordinate. A wout file cannot change this; the routes B1, B2 and G pass GVEC's own poloidal flux.
- PB3D never reads `signgs`. It assumes VMEC's left-handed coordinates $(s, \theta, \varphi)$ with a negative Jacobian, for example in `calc_flux_q_VMEC` ("conversion VMEC LH -> PB3D RH", `eq_ops.f90:368`), in `calc_eq_2` (`pmone = -1`, line 522), in `grid_ops.f90:176` and 889, `eq_utilities.f90:786` and 866, `grid_utilities.f90:211` and 440.
- $R_0$ = `R_V_c(1,1,0)` (`calc_normalization_const_VMEC`) assumes that the mode (0,0) comes first in the main mode list, and `read_VMEC` stops if (0,0) is not the first Nyquist mode. `read_VMEC` checks `mnmax` $= (\mathrm{ntor}+1) + (2\,\mathrm{ntor}+1)(\mathrm{mpol}-1)$, VMEC's full rectangular mode set. The Nyquist list may have any size, in VMEC's order (`VMEC_ops.f90:161-166`).
- `PB3D_out.h5` stores the `VMEC_vars` arrays, and POST reads them back (`PB3D_ops.f90:162-309`); POST never opens the equilibrium file.
- For every point of the equilibrium grid, `PB3D_out.h5` also holds the point's VMEC coordinates (`r_E`, `theta_E`, `zeta_E`) and the quantities that enter the energy: `g_FD`, `h_FD`, `jac_FD`, `S`, `kappa_n`, `kappa_g`, `sigma` (`print_output_grid`, `print_output_eq_2`). A script can therefore compare them with an independent evaluation at the same points (test T2).

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

### 1.4 Where PB3D branches on the equilibrium style

`files_ops.f90:249-268` sets `eq_style` from the equilibrium file: 1 (VMEC) when the name ends in `.nc`, 2 (HELENA) when the file starts with an integer. PB3D then branches on `eq_style` in 65 `select case` blocks in 17 files and in 38 `if` tests. Every one of the 65 blocks has the cases 1 and 2 and no default, so for any other value a block does nothing and reports nothing. For its other styles PB3D has a default that stops with a message ("No BC style associated with ...", `SLEPC_ops.f90:563`); stage S0 gives the 65 blocks the same.

| what the block decides | blocks | code lines in the VMEC branches |
|---|---|---|
| how the run is organised: which grids exist, what is computed at each Richardson level, the names of variables in the output, default input values, messages | 39 | 158 |
| the equilibrium itself: reading the file, storing the arrays of section 1.1, the fluxes, $R$, $Z$, $\lambda$, the metric, the transformation matrices, the derived quantities, positions in space, the normalization, the debug tests | 26 | 507 |

Of the 507 lines, 275 write the arrays of section 1.1 to `PB3D_in.h5` (`print_output_in`, 156 lines) and read them back (`reconstruct_PB3D_in`, 119 lines). Most of the other VMEC branches are calls of routines. Of the 38 `if` tests, 16 test for HELENA and 22 for VMEC.

The HELENA branches differ from the VMEC branches in nearly all 65 blocks. HELENA has its own arrays (`HELENA_vars`) and its own metric routines, and PB3D computes a HELENA equilibrium once on HELENA's grid and interpolates it to the field-aligned grid, while it computes a VMEC equilibrium on the field-aligned grid at every Richardson level.

Outside `read_VMEC`, these places name the arrays of `VMEC_vars` directly:

| place | code lines |
|---|---|
| VMEC branch of `print_output_in` (`input_ops.f90:1372`) | 156 |
| VMEC branch of `reconstruct_PB3D_in` (`PB3D_ops.f90:162`) | 119 |
| `calc_XYZ_grid_VMEC` (`grid_utilities.f90:871`) | 74 |
| `calc_normalization_const_VMEC` (`eq_ops.f90:5358`) | 65 |
| `coord_F2E_VMEC` (`grid_utilities.f90:186`) | 52 |
| `coord_E2F_VMEC` (`grid_utilities.f90:415`) | 34 |
| `calc_RZL_ind` (`eq_ops.f90:2713`) | 26 |
| `calc_flux_q_VMEC` (`eq_ops.f90:345`) | 19 |
| `normalize_VMEC` (`VMEC_ops.f90:380`) | 19 |
| `calc_jac_V_ind` (`eq_ops.f90:3048`) | 17 |
| VMEC branch of `calc_norm_range_PB3D_in` (`grid_ops.f90:170`) | 10 |
| total | 591 |

The debug tests `test_calc_RZL_VMEC`, `test_sigma_VMEC`, `test_B_F` and `check_vac_exact_1` name them too. The other routines of the VMEC path (`calc_g_C`, `calc_T_VC`, `calc_g_V`, `calc_T_VF`, the `calc_derived_` routines, `fourier2real`, `calc_trigon_factors`) get their data through arguments.

### 1.5 How accurate the equilibrium data must be

Two cases from PB3D's history show how its eigenvalues react to errors in the equilibrium data.

- Noise, an error that changes from one grid point to the next. A HELENA mapping file with 8 significant digits has relative errors up to $5 \cdot 10^{-9}$. In the test of master of 26 and 27 November 2026, runs on such files were compared with runs of the same cases on files with 16 digits. On mappings with 1025 radial points the eigenvalues differed by a median of $10^{-3}$ to $10^{-2}$ relative and by several per cent in the worst cases; the largest move of a benchmark value was 4.2 %. On mappings with 513 points the differences were smaller. PB3D warns about such a file (`read_HEL`, issue #65), and names the rounded poloidal grid as the main cause.
- A smooth error, with the same sign over many surfaces. Until 17 November 2026 `read_VMEC` placed VMEC's half mesh at a slightly wrong radius. This put a relative error of $5.4 \cdot 10^{-6}$ into $\lambda$ and of $1.1 \cdot 10^{-6}$ into the Jacobian. The correction (PR #34) moved the eigenvalues by 2 to $4 \cdot 10^{-6}$ relative on the circular twin, the rippled twin and W7-X.

So a smooth error moves the eigenvalue by about its own relative size, and noise moves it by five to six orders of magnitude more than its own size. The reason is that PB3D interpolates and differentiates tabulated data with cubic splines (for a VMEC file along the radius, section 1.1), and a spline derivative of order $k$ multiplies noise by about the number of points to the power $k$.

The table shows this for the test function $\cos 3s + 0.2\,s^2$, which has values of order 1 and a third derivative up to 27. The function is sampled on `ns` points equidistant in $s$, relative noise is added to the samples, and a cubic spline with not-a-knot ends gives the derivatives. It stands in for PB3D's spline, which is PSPLINE's cubic spline with its default end condition. Each entry is the largest absolute error of the derivative for $0.05 \le s \le 0.95$.

| relative noise in the values | `ns` | derivative of order 1 | order 2 | order 3 |
|---|---|---|---|---|
| none | 319 | $2.0 \cdot 10^{-8}$ | $6 \cdot 10^{-5}$ | 0.12 |
| $10^{-12}$ | 319 | $2.0 \cdot 10^{-8}$ | $6 \cdot 10^{-5}$ | 0.12 |
| $10^{-10}$ | 319 | $8.6 \cdot 10^{-8}$ | $1.1 \cdot 10^{-4}$ | 0.16 |
| $5 \cdot 10^{-9}$ (8 digits) | 319 | $3.5 \cdot 10^{-6}$ | $3.2 \cdot 10^{-3}$ | 2.2 |
| none | 1025 | $6.0 \cdot 10^{-10}$ | $6 \cdot 10^{-6}$ | 0.04 |
| $10^{-12}$ | 1025 | $2.5 \cdot 10^{-9}$ | $1.2 \cdot 10^{-5}$ | 0.05 |
| $10^{-10}$ | 1025 | $2.1 \cdot 10^{-7}$ | $7.1 \cdot 10^{-4}$ | 1.6 |
| $5 \cdot 10^{-9}$ (8 digits) | 1025 | $1.2 \cdot 10^{-5}$ | $3.9 \cdot 10^{-2}$ | 83 |

Which derivative orders the energy reads. PB3D's energy terms read the metric and the Jacobian in flux coordinates with at most one derivative in the radial coordinate: `g_FD` and `jac_FD` with the derivative indices (0,1,0). The second derivatives that `calc_U` reads are in the two angular coordinates (`h_FD` with the indices (1,0,1) and (0,0,2), `X_ops.f90:2286-2292`). The metric holds first radial derivatives of $R$ and $Z$, and the transformation to flux coordinates holds the first radial derivative of $\lambda$. So the energy needs $R$, $Z$ and $\lambda$ up to the second radial derivative. PB3D also computes and stores the third (`max_deriv` + 1 = 3), but no energy term reads a quantity that depends on it. This follows from a search of the source for the derivative indices with which `g_FD`, `h_FD` and `jac_FD` are read; stage S3 confirms it with a run.

The requirements, as far as they are known:

- R1, noise. The values on the flux surfaces carry noise below $10^{-12}$ relative for 319 surfaces and below $10^{-13}$ for 1025, so that the noise does not add to the spline's own error. Rounding in double precision ($10^{-16}$) meets this. A text file with 8 digits, a variable in single precision, or a step that stops at an iteration tolerance does not. Where PB3D does not differentiate the data itself (routes B1, B2 and G), noise in the values is not multiplied.
- R2, smooth errors. PB3D on a HELENA equilibrium and on a VMEC equilibrium of the same case agrees to $1.1 \cdot 10^{-6}$ (circular twin, fixed boundary). A route should add smooth errors well below that, below $10^{-7}$, in the values and in the first radial derivatives.
- R3, the second radial derivatives. What PB3D's eigenvalue needs of them is not measured. With a VMEC file of 319 surfaces PB3D works today with errors of $6 \cdot 10^{-5}$ in the interior and of about $10^{-3}$ at the plasma edge (section 4.1). Stage S3 measures the requirement.

For the comparison of eigenvalues itself: with `tol_SLEPC` = $10^{-6}$ two runs of the same build can differ by more than the differences this plan looks for. Comparisons tighter than $10^{-5}$ use `tol_SLEPC` = $10^{-12}$; at that setting the runs of the test of 26 and 27 November 2026 agreed to $4 \cdot 10^{-7}$ of the scale of the spectrum or better, also between different numbers of processes.

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

The signs of `phi`, `iotaf` and `chi` were checked on the circular VMEC twin (`wout_circ_eps033_q104.nc`, VMEC++, `signgs` = -1, 319 surfaces): `phi`(edge) = -0.32945, `iotaf` from -0.9605 to -0.9632, $d\,$`phi`$/ds$ = `phipf` and `chipf` = `iotaf` $\cdot$ `phipf` to rounding, `chi`(edge) = -0.3168 with $d\,$`chi`$/ds = -$`chipf`, and `bsupvmnc`(0,0) = -1.08 where GVEC's $B^\zeta$ is positive. GVEC's own VMEC reader divides `chi` by $2\pi$ without the factor `signgs` (`vmec_readin.F90:638`), so the inverse of that reader would give the wrong sign of `chi`. PB3D does not read `chi`, but other tools do. Test T1 and gate A1b check all of this again on the converted files.

### 2.4 What pyGVEC evaluates

pyGVEC (`gvec.load_state`, `state.evaluate`) evaluates at given $(\rho, \vartheta, \zeta)$: $X^1$, $X^2$, $\lambda$ with radial derivatives of order 0 to 2 (`_evaluate.py:73`, the pattern `r{0,2}`); the position, the basis vectors, the metric and the Jacobian; the vectors $B$ and $J$ in Cartesian components, $\lvert B \rvert$, the curvature vector of the field lines (`kappa_B`) and its normal and geodesic parts; $\Phi$, $\chi$, $\iota$, $p$ and their radial derivatives. These are direct evaluations of GVEC's B-splines and Fourier series, with rounding as their only error.

pyGVEC also evaluates in other angles (PEST, Boozer). This needs a root search or a transform, and the result is no longer a finite Fourier series: GVEC's CAS3D writer works with a sampling factor and cuts the series (`--sampling`, `--MN_booz` in `cas3d.py`). No route of this plan changes the angles. PB3D keeps GVEC's $\vartheta$ and $\zeta$ up to the sign of $\zeta$ (section 2.3) and makes its own transformation to field-aligned coordinates afterwards, as it does for VMEC.

GVEC does not write wout files: not in version 1.5.1 and not in the development version (1.5.2.dev88 on 27 November 2026, same list of `pygvec` subcommands). A search of GVEC's public issue list for "wout" returns nothing; the nearest entry is issue 22 (open), "Frenet → RZφ conversion", about bringing a solution in the Frenet frame to cylindrical coordinates. I found no other tool that writes a wout file from a GVEC state.

### 2.5 What a GVEC run leaves on disk

A GVEC run leaves two files that together define the equilibrium, read here from the 1.5.1 source; GVEC has not been run for this plan.

- The parameter file holds the map and its parameters, the profiles and the resolution.
- The state file (for example `GVEC_State_0000_00000100.dat`) is a text file, written by `WriteState` in `src/restart/restart.F90`. It holds the radial grid, the number of field periods, the number of the map and, for each of $X^1$, $X^2$ and $\lambda$, the degree and continuity of the B-splines, the list of Fourier modes and the B-spline coefficients of every mode. It also holds a table of $\Phi$, $\chi$, $\iota$ and $p$ at the interpolation points of the splines, and the minor radius, the major radius and the volume.

What matters for accuracy:

- The coefficients are written with 15 significant digits (format `E23.15`), a relative rounding of up to $5 \cdot 10^{-16}$. The same coefficients serve all radii, so this rounding is a smooth error in the sense of section 1.5, not noise.
- GVEC's Fourier basis is the plain $\sin$ and $\cos$ of $m\vartheta - n N_\mathrm{fp}\zeta$ (`src/base/fbase.F90:13`), with $n \ge 0$ for $m = 0$ and $-n_\mathrm{max} \le n \le n_\mathrm{max}$ for $m \ge 1$. With $\varphi = -\zeta$ the mode $(m, n)$ becomes VMEC's mode $(m, -n)$, and the sine coefficients with $m = 0$ change sign; no arithmetic is involved.
- So everything PB3D needs of $X^1$, $X^2$ and $\lambda$ is a finite Fourier series in the angles whose coefficients are piecewise polynomials in $\rho$. Values on a grid followed by an FFT, or the evaluation of each mode's B-spline, both return the coefficients on a surface with rounding as the only error. Measured for the grid and the FFT, on a finite series with $m, n \le 12$ and coefficients between 1 and $10^{-12}$: every coefficient is returned to $5 \cdot 10^{-16}$ absolute.
- The frame of the G-Frame map is a finite Fourier series too: GVEC transforms the points of the netCDF file and works with the modes up to $(n_\zeta - 1)/2$ per field period (`src/hmap/hmap_axisNB.F90:20-21`, 183).
- The state file holds no Jacobian, no metric and no magnetic field, and it is not complete without the parameter file: the parameters of the map are there, and GVEC evaluates the profiles from their definition there, not from the table in the state file.
- pyGVEC's `State` class shows no access to the B-spline coefficients (`python/gvec/core/state/_class.py`, `_binding.py`). A script that wants them reads the state file.

A run also writes plot files (VTK, netCDF, CSV: values on a plot grid) and a log. Files for other codes come from separate tools that read the two files above: `pygvec to-cas3d` (netCDF in Boozer angles), `pygvec to-gist`, a script of SIMPLE, and Fortran converters for CASTOR3D, JOREK, GENE and HOPR, the last four for the cylindrical map only.

## 3. Routes

Every route starts from the same GVEC run and ends in the same part of PB3D.

- The start. GVEC has computed an equilibrium and written its parameter file and its state file (section 2.5). How accurate the equilibrium itself is depends on GVEC's resolution and convergence, which are the same for every route; stage S2 lists what a GVEC run for PB3D needs.
- The end. Once the arrays of section 1.1 are filled, PB3D does what it does for a VMEC file. It normalizes the arrays and stores them in `PB3D_in.h5`. At every Richardson level it sets up the field-aligned grid, sums the Fourier series of $R$, $Z$ and $\lambda$ and of their derivatives at the grid points (`fourier2real`), builds the metric and the transformation to flux coordinates, computes the curvature, the shear and the parallel current, integrates the energy terms of the perturbation along the field lines, and solves the eigenvalue problem, with the vacuum for a free boundary. `PB3D_out.h5` holds the results, and POST reads them.

The routes differ in what lies between: what goes from GVEC to PB3D, and who takes the radial derivatives.

| route | what goes from GVEC to PB3D | the radial derivatives come from | `eq_style` | GVEC's maps | Fortran that changes |
|---|---|---|---|---|---|
| A | a wout file, written by a pyGVEC script | PB3D's cubic splines | 1 | cylindrical | none |
| B1 | GVEC's parameter file and state file | GVEC's library, called by PB3D | 1 | cylindrical | a reader, the build |
| B2 | a wout file with extra variables, written by a pyGVEC script | GVEC's B-splines, in the script | 1 | cylindrical | a second branch in `read_VMEC` |
| C | added to A, B1 or B2: the frame of the map | as in the route it is added to | 1 | all | six routines, for all files or behind a switch |
| G | a file of PB3D's own design, written by a pyGVEC script | GVEC's B-splines, in the script | 3 | all | a reader, 64 blocks, six routines for style 3 |

The errors quoted in this section are for the test function of section 1.5 and 319 surfaces unless stated otherwise; section 4.1 has them side by side.

### Route A: a wout file, PB3D unchanged

The path of an equilibrium:

1. GVEC writes the parameter file and the state file.
2. The user runs `gvec_to_wout.py` on the two files, with the number of surfaces `ns` as an argument. The script loads the state with pyGVEC, evaluates it on `ns` surfaces equidistant in $s$ and writes a wout file with the transformations of section 2.3:
   - $R$, $Z$ on the full mesh and $\lambda$ on the half mesh (rows 2 to `ns`, section 1.1), from a 2-D FFT of pyGVEC's values on a $(\theta, \varphi)$ grid into VMEC's basis $\cos(m\theta - n\varphi)$, $\sin(m\theta - n\varphi)$, in VMEC's mode order with (0,0) first. The grid has more than twice as many points as the largest mode number in each angle, so the FFT returns the finite Fourier series to rounding (section 2.5), and GVEC's internal mode order never enters. Modes outside GVEC's set are zero.
   - `gmnc` on the half mesh on a Nyquist list large enough to hold the Jacobian. For map 1, $\sqrt g = R\,(R_\theta Z_s - R_s Z_\theta)$ is a product of three finite Fourier series, so it is itself a finite series with $m \le 3\,m_\mathrm{max}$, $\lvert n \rvert \le 3\,n_\mathrm{max}$ (with different limits $M_1$, $M_2$ for $X^1$ and $X^2$: $2M_1 + M_2$). VMEC truncates it (the circular fixture's Nyquist list stops at $m = 15$ for `mpol` = 12, where the exact Jacobian reaches $m = 33$); the converter need not. For `mpol` = `ntor` = 12 the full list has about 2,400 modes, and in PB3D `fourier2real` costs time in proportion to the number of modes at every grid point. So the converter writes the smallest full rectangular Nyquist set outside which all coefficients stay below a tolerance on all surfaces. The same set on all surfaces makes the cut a smooth error; the tolerance is $10^{-10}$ of the (0,0) coefficient or smaller (requirement R2).
   - `bmnc`, $B_s$, $B_u$, $B_v$, $B^u$, $B^v$ on the same list. These are not finite series, so the list truncates them; PB3D uses `bmnc` only for $B_0$ and the $B$ components only in debug tests.
   - the profiles and scalars of section 1.1, `raxis_cc`, `zaxis_cs` of length `ntor`+1, `ier_flag` = 0, `lfreeb` = F, `lrfp` = F, `lasym` from GVEC's sine and cosine choice, `signgs` = -1, `nextcur` = 0. The file is written with netCDF4 in double precision, with the names, dimensions and types of the VMEC++ fixture, which LIBSTELL is known to read (section 1.3).
   - before writing, the converter checks that the (0,0) coefficient of `gmnc` is negative on every surface, and that the equilibrium uses the cylindrical map: it obtains $X_0$, $N$, $B$ as in section 2.1 and stops with a message unless $X_0 = 0$, $N = (\cos\zeta, -\sin\zeta, 0)$ and $B = (0, 0, 1)$.
3. The user starts PB3D with the wout file. LIBSTELL reads it and `read_VMEC` fills the arrays of section 1.1. Cubic splines through the `ns` surfaces give the radial derivatives of $R$ and $Z$. A cubic spline moves $\lambda$ and the Jacobian from the half mesh to the full mesh, with an extrapolation over half a step at the axis and at the edge, and gives their radial derivatives. The poloidal flux is the trapezoidal integral of `iotaf` times `phipf`.
4. From the arrays on, PB3D runs as for any VMEC file.

Numerical error on this path:

- Step 2 adds rounding to $R$, $Z$ and $\lambda$, and the cut of the Jacobian's mode list at its tolerance.
- Step 3 adds the errors of PB3D's VMEC reader. They are the errors PB3D has today for every VMEC file.
  - $R$ and $Z$: first radial derivative $2 \cdot 10^{-8}$ in the interior and $4.5 \cdot 10^{-7}$ at the edge; second radial derivative $6 \cdot 10^{-5}$ in the interior and $5.6 \cdot 10^{-4}$ at the edge.
  - $\lambda$ and the Jacobian, through the half mesh: values $2 \cdot 10^{-11}$ in the interior and $1.6 \cdot 10^{-9}$ at the edge; first radial derivative $2 \cdot 10^{-8}$ and $1.8 \cdot 10^{-6}$; second radial derivative of $\lambda$ $6 \cdot 10^{-5}$ and $1.2 \cdot 10^{-3}$.
  - the poloidal flux: the trapezoidal integral differs from the integral of a cubic spline through the same values by $1.2 \cdot 10^{-9}$ relative on the circular VMEC twin, whose $\iota$ is nearly constant; the error is proportional to $h^2$, $h = 1/(\mathrm{ns} - 1)$, and to the curvature of `iotaf` times `phipf`.
- These errors are largest at the plasma edge, where the spline has data on one side only. The eigenfunctions of peeling-ballooning instabilities sit at the edge, and the vacuum normal uses the Jacobian extrapolated to the boundary (section 1.1).
- With twice as many surfaces the first derivative improves by a factor of 8 and the second by a factor of 4 (measured at 639 surfaces: $5.5 \cdot 10^{-8}$ and $1.4 \cdot 10^{-4}$ at the edge for $R$ and $Z$, $2.2 \cdot 10^{-7}$ and $2.9 \cdot 10^{-4}$ for $\lambda$ and the Jacobian). PB3D's memory and time grow with `ns`, because the equilibrium grid is a subset of the file's surfaces.
- Noise: the values carry rounding only, which the splines do not make visible (first rows of the table in section 1.5). Requirement R1 forbids any step in the script that stops at a tolerance or evaluates in other angles.

What changes in PB3D: nothing.

What speaks for route A:

- No Fortran changes, so no result of PB3D on VMEC or HELENA files can change.
- A GVEC equilibrium and a VMEC equilibrium go through the same statements of PB3D. A difference between their eigenvalues then comes from the two equilibria alone.
- The file also works in other tools that read VMEC files (BOOZ_XFORM, COBRAVMEC, TERPSICHORE), and GVEC has no wout writer (section 2.4).
- It gives the first eigenvalue of a GVEC equilibrium soonest.

What speaks against route A:

- It covers the cylindrical map only.
- GVEC has the radial derivatives exactly, and this route discards them: PB3D computes them again from values on surfaces, with errors five to seven orders of magnitude larger than those of routes B1, B2 and G, and largest at the edge. Whether errors of this size matter for the eigenvalue is requirement R3, which is not measured.
- More accuracy costs more surfaces, and with them memory and time in PB3D.
- The script has to reproduce VMEC's file layout and LIBSTELL's rules (section 1.3). LIBSTELL skips a variable with a wrong type or name without a message, so a mistake there shows up only in the checks.

Effort (guess): 3 to 5 days for the converter and its checks, 1 to 2 days for the comparison script of test T2.

### Route B: a GVEC reader in PB3D, `eq_style` stays 1

PB3D fills the arrays of section 1.1 from GVEC data directly: Fourier coefficients of $R$, $Z$, $\lambda$ with radial derivatives taken from GVEC's B-splines, all on the full mesh, without half mesh or extrapolation. `eq_style` stays 1. Only `files_ops.f90` (which file is which) and `read_input_eq` (`input_ops.f90:1175-1178`) change, and one more entry in `misc_in_V` records where the arrays came from. There are two ways to get the data.

#### B1: PB3D links GVEC's library

The path of an equilibrium:

1. GVEC writes the parameter file and the state file.
2. The user starts PB3D with the two files and the number of surfaces. There is no script and no file in between.
3. On PB3D's root process a new reader calls GVEC's library `libgveclib`. The library holds the module that pyGVEC is built on (`MODgvec_py_state`, `src/pygvec/state.F90`): `Init(parameterfile)` and `ReadState(statefile)` load the equilibrium, `evaluate_base_tens` evaluates $X^1$, $X^2$, $\lambda$ and their derivatives, `evaluate_profile` the profiles. The reader evaluates on a grid in the angles on every surface, takes the FFT, and fills the arrays of section 1.1.
4. From the arrays on, PB3D runs as for a VMEC file.

Numerical error on this path: rounding only. The first and second radial derivatives come from GVEC's own B-spline routines (measured with a B-spline of degree 5 in $\rho$ and the chain rule to $s$: below $3 \cdot 10^{-13}$ and $8 \cdot 10^{-11}$). The module's routines offer radial derivatives up to the second order, which is what the energy reads (section 1.5); the third, which PB3D stores, would have to come from GVEC's spline base directly.

What changes in PB3D: the reader, the recognition of the two files, and the build, which needs GVEC as an optional library, as LIBSTELL is needed for VMEC.

What speaks for B1:

- The user runs no script, and no file format lies between the two codes.
- GVEC's own code evaluates the equilibrium, so PB3D cannot differ from GVEC in a convention of the B-splines, the profiles or the map.
- The radial derivatives are exact, and PB3D chooses the number of surfaces in its own input.

What speaks against B1:

- GVEC becomes a compiled dependency: built with PB3D's compiler on the Mac, on Snellius and in CI, for every user who wants GVEC equilibria.
- `MODgvec_py_state` is made for pyGVEC and is not a published interface; it can change between GVEC versions.
- The module keeps its state in module variables, prints to standard output and stops the program itself on an error. It stops when GVEC is compiled with MPI and runs on more than one process, so GVEC has to be built without MPI (its default) and called on the root process only.
- For the other maps PB3D needs the frame with three derivatives at any $\varphi$ (route C). Whether the library gives this for every map has not been checked.

Effort (guess): two to three weeks, a large part of it the build on three systems.

#### B2: the wout file of route A with extra variables

The path of an equilibrium is that of route A, with two differences.

- In step 2 the script writes, next to the wout variables, the Fourier coefficients of $R$, $Z$, $\lambda$ and the Jacobian on the full mesh with their radial derivatives, and GVEC's poloidal flux, as extra variables. LIBSTELL does not ask for variables it does not know (section 1.3), so the file stays readable as a wout file.
- In step 3 `read_VMEC` uses the extra variables when they are present, reading them with LIBSTELL's `ezcdf` routines, which PB3D's LIBSTELL build already holds, and does not use its splines.

Numerical error on this path: rounding only, as in B1. pyGVEC gives the radial derivatives up to the second order in $\rho$; the chain rule gives those in $s$. The third, which the energy does not read, comes from a finite difference of pyGVEC's second derivative or from the B-spline coefficients in the state file; stage S3 decides whether it needs care at all. The coefficients of odd $m$ behave like $\sqrt s$ near the axis, so their derivatives in $s$ are unbounded at $s = 0$: the script has to decide what it stores on the axis surface (GVEC's own writers start at $s = 10^{-8}$, `gvec_to_castor3d`, or avoid the axis, `cas3d.py:119`).

What changes in PB3D: a second branch in `read_VMEC`, and one entry in `misc_in_V`.

What speaks for B2:

- The radial derivatives are exact, with the Fortran change limited to one routine and no edits of the 65 blocks.
- One file serves PB3D and the other tools that read VMEC files.

What speaks against B2:

- The file is a wout file with variables that VMEC does not have, and it has to satisfy LIBSTELL's rules all the same.
- `read_VMEC` gets two ways of filling the same arrays. If an extra variable is missing or misnamed, PB3D falls back to its splines; the reader has to say which way it took.
- PB3D's output does not say that the equilibrium is GVEC's.

Effort (guess): one to two weeks on top of route A.

### Route C: GVEC's other maps (knot, Frenet frame, G-Frame), added to route A or B

Because every map has the form of section 2.1, PB3D needs one generalisation, from the cylindrical map to

$$x = X_0(\varphi) + R\, N(\varphi) + Z\, B(\varphi),$$

where $R$ and $Z$ now stand for $X^1$ and $X^2$, and the cylindrical map is the case $X_0 = 0$, $N = e_R(\varphi)$, $B = e_Z$. With $\varphi = -\zeta$ the handedness that PB3D assumes holds for every map (section 2.3), so no sign convention changes.

The path of an equilibrium is that of the route it is added to, with two additions. The script (or the reader of B1) also passes the frame: the nine Cartesian components of $X_0$, $N$, $B$ as Fourier series in $\varphi$, obtained as in section 2.1 for any map. PB3D reads them, evaluates them and their first three derivatives at the grid's $\varphi$ values, stores them in `PB3D_out.h5`, and POST reads them back. A file without them is a VMEC file and gets the cylindrical frame.

What stays as it is: the file's `rmnc`, `zmns`, `lmns` (holding $X^1$, $X^2$, $\lambda$), the profiles, `gmnc` (the script writes GVEC's full Jacobian, which includes the Jacobian of the map), $T_{VC}$, `calc_g` and everything after it.

What changes in PB3D is the list of section 1.2:

- `calc_g_C_ind`: the metric of $(R, \varphi, Z)$ becomes, with a prime for $d/d\varphi$ and $t = \partial x/\partial\varphi$ at fixed $R$ and $Z$, $t = X_0' + R\,N' + Z\,B'$,

$$g_C = \begin{pmatrix} N\cdot N & N\cdot t & N\cdot B \\ \cdot & t\cdot t & t\cdot B \\ \cdot & \cdot & B\cdot B \end{pmatrix},$$

  a polynomial of degree 2 in $R$ and $Z$ with coefficients that depend on $\varphi$. Its derivatives up to total order 2 need $X_0$, $N$, $B$ up to their third derivative. For the cylindrical map $N' = e_\varphi/R$ has unit length, $B' = 0$, and $g_C = \mathrm{diag}(1, R^2, 1)$.
- `calc_derived_DE_epar_VMEC`: the derivatives of the basis vectors of VMEC's coordinates, $e_\theta = R_\theta N + Z_\theta B$ and $e_\varphi = t + R_\varphi N + Z_\varphi B$, in Cartesian components. For example $\partial_\varphi e_\varphi = X_0'' + R\,N'' + Z\,B'' + 2R_\varphi N' + 2Z_\varphi B' + R_{\varphi\varphi} N + Z_{\varphi\varphi} B$, which for the cylindrical map ($N'' = -N$) gives the present terms with $R_{\varphi\varphi} - R$ and $2R_\varphi$.
- the vacuum block, `calc_XYZ_grid_VMEC` and step 5 of `calc_vec_comp`: positions and vector components from $X_0$, $N$, $B$ instead of $\cos\varphi$ and $\sin\varphi$.
- the normalization length: the length of the magnetic axis divided by $2\pi$ (GVEC's `r_major`).
- the debug tests `test_g_V`, `test_jac_V` and the diagnostic `check_vac_exact_1`.

Numerical error that this route adds: for the knot and the G-Frame the frame is a finite Fourier series (section 2.5), so the script passes it with rounding only. The Frenet frame of map 20 is built from the derivatives of the axis curve and is not a finite series in general; its Fourier list is cut at a tolerance, a smooth error, as for the Jacobian in route A. With a frame that is a finite series the Jacobian is a finite series too, with a longer mode list; with the Frenet frame it is not.

Two points that are specific to the other maps, in every route that includes them:

- The toroidal mode number. PB3D keeps one toroidal mode number $n$ in its flux coordinates, whose toroidal angle is the file's $\varphi$ up to the sign. For a G-Frame, $\zeta$ is the parameter of the curve, which GVEC's documentation recommends to be a Boozer angle, not the geometric toroidal angle. The same 3-D equilibrium described with two different toroidal angles gives PB3D two different sets of displacements for the same $n$, so the eigenvalue can depend on the parametrisation. This follows from the form of PB3D's Fourier series; how large the effect is has not been computed (stage C4).
- A reference. For shapes that cylindrical coordinates cannot describe there is no VMEC equilibrium to compare with. The check that remains valid is the pointwise comparison of PB3D's stored quantities with GVEC's own values (test T2), which does not depend on the map.

What speaks for route C: it opens every GVEC map to PB3D with one change of six routines, and it can be tested without GVEC (stage C2).

What speaks against route C:

- With the wout file as the carrier, the file is no longer a VMEC file for other tools once it holds a frame.
- The six routines are on the path of every VMEC run. Either VMEC files go through the general formulas too, and their results can change at the level of rounding, or a switch keeps the cylindrical statements for them (decision D4).

Effort (guess): two to three weeks, most of it tests.

### Route G: a file of PB3D's own design, read as a third equilibrium style

The path of an equilibrium:

1. GVEC writes the parameter file and the state file.
2. The user runs `gvec_to_pb3d.py` on the two files, with the number of surfaces `ns` as an argument. The script loads the state with pyGVEC and writes a netCDF file in double precision that holds what PB3D keeps of an equilibrium (section 1.1), in SI units, in the coordinates $(s, \theta, \varphi)$ of section 2.3, on `ns` surfaces equidistant in $s$:
   - the Fourier coefficients of $X^1$, $X^2$ and $\lambda$, all three on the same surfaces, with their derivatives in $s$. The values and the derivatives of order 1 and 2 are pyGVEC's, from a grid in the angles and an FFT (section 2.5) and the chain rule from $\rho$ to $s$. For the order 3 the same holds as in B2;
   - the toroidal flux, GVEC's poloidal flux $\chi$, $\iota$ and the pressure, with derivatives of order 0 to 2;
   - $\lvert B \rvert$ on the axis, the number of field periods, the symmetry, the mode list and the number of GVEC's map;
   - the covariant components of $B$, which only the debug tests use;
   - from stage G3 on, the frame $X_0$, $N$, $B$ of section 2.1.
3. The user starts PB3D with this file. PB3D recognises it, sets `eq_style` = 3, and a new routine `read_GVEC` copies the file's arrays into the arrays of `VMEC_vars`. No spline and no integral is involved.
4. From the arrays on, PB3D runs the statements of the VMEC path, with two differences for style 3. It computes the Jacobian at its grid points from the derivatives of $R$ and $Z$ that it has there,

$$\mathcal J = (R_s Z_\theta - R_\theta Z_s)\; N \cdot (B \times t),$$

   which for the cylindrical map is $R\,(R_\theta Z_s - R_s Z_\theta)$; the file holds no Jacobian. And from stage G3 on it uses the frame formulas of route C in the six routines of section 1.2, also for a GVEC equilibrium with the cylindrical map.

Numerical error on this path: rounding only, as in B1 and B2, and no mode list is cut: the Jacobian is computed where it is needed, and it agrees with the metric by construction. On the axis surface the same decision is needed as in B2. For the other maps the remarks of route C on the frame hold.

A VMEC file and this file both end in `.nc`, so PB3D tells them apart by a global attribute of the file. netCDF is the format of GVEC's own files for other codes (the CAS3D file, the file of the G-Frame); xarray writes it, and PB3D reads it with netCDF-Fortran, which it links already.

What changes in PB3D. How much Fortran a third style repeats depends on where `read_GVEC` puts the data (the counts are in section 1.4):

- With arrays of its own, as HELENA has, the 11 places that name the arrays of `VMEC_vars` would exist twice with other array names: 591 code lines, plus the debug tests.
- With the arrays of `VMEC_vars`, nothing exists twice. 64 of the 65 blocks get `case (1,3)`, and one block (`read_input_eq`) gets a case of its own that calls `read_GVEC`. The 22 `if` tests for VMEC are extended to style 3 where they concern the Fourier series.

This plan takes the arrays of `VMEC_vars`. Inside PB3D a GVEC equilibrium then differs from a VMEC equilibrium in four things:

- the reader and the recognition of the file (`read_GVEC` in a new module, `files_ops.f90`, `read_input_eq`). My estimate for `read_GVEC` is 100 to 150 code lines; `read_VMEC` has 263, a large part of them for the half mesh and the splines;
- the Jacobian, computed at the grid points for style 3;
- the frame, from stage G3 on: the six routines of section 1.2 get the formulas of route C for style 3. The statements for VMEC files stay as they are, so decision D4 does not arise;
- messages and the debug tests that name VMEC.

The arrays and routines keep their names with `V`. `VMEC_vars` then holds an equilibrium as Fourier series in VMEC's coordinates, whichever code computed it (decision D6).

What speaks for route G:

- The radial derivatives are exact, the Jacobian needs no mode list, and the poloidal flux is GVEC's.
- The statements that run for VMEC and HELENA files are not edited, so their results cannot change; the edits of the 64 blocks add a case and change none.
- All of GVEC's maps fit without a file that pretends to be VMEC's, and the frame formulas are selected by the style.
- The script uses pyGVEC's public interface only, and CI needs no GVEC: a file is the fixture.
- PB3D's messages, input checks and output say that the equilibrium is GVEC's.

What speaks against route G:

- It needs the most Fortran before the first eigenvalue of a GVEC equilibrium: the reader, the 64 blocks, the 22 tests and the Jacobian.
- A third style is one more thing that every later change of PB3D has to consider in 65 blocks, even where the case is shared with VMEC.
- The file is PB3D's own: other tools cannot read it, and PB3D has to define and keep the format.
- A GVEC equilibrium and a VMEC equilibrium no longer go through the same statements (derivatives, Jacobian), so a comparison of their eigenvalues also contains the difference between PB3D's two ways of reading.
- The arrays and routines named after VMEC then also serve GVEC (decision D6).

Effort (guess): the script 3 to 5 days; the reader, the edits of the blocks, the Jacobian and their tests about one week; the comparison script of test T2 1 to 2 days; the frame of stage G3 two to three weeks, most of it tests.

A variant of this route, PB3D reading GVEC's state file with Fortran of its own and no script, is not followed. The state file lacks the parameters of the map and the definition of the profiles (section 2.5), so PB3D would have to read GVEC's parameter file too and repeat GVEC's B-splines, its profile types and its treatment of the G-Frame in Fortran, following file formats that are not a published interface. Taking the profiles from the table in the state file would break requirement R2: the table has about as many points as the radial grid has elements, and stage S2 shows that even 101 points change $\iota'$ by 6 to 12 % of the local shear on the circular twin.

## 4. Comparison

### 4.1 Accuracy

The table gives the error that each route adds between GVEC's solution and PB3D's arrays. GVEC's own discretisation error is the same for every route and is not in the table. The Jacobian itself holds first radial derivatives of $R$ and $Z$, so its first radial derivative stands at the level of their second. The numbers are measured on test functions with values of order 1, not on GVEC data: for the wout routes the function of section 1.5 on 319 surfaces through a cubic spline with not-a-knot ends, and for the half mesh the same spline on the half-mesh points with extrapolation over half a step; for the other routes a B-spline of degree 5 in $\rho$ on 20 and on 100 elements, with coefficients rounded to 15 digits and the chain rule from $\rho$ to $s$, for $0.05 \le \rho \le 1$.

| quantity | route A, and C on A | routes B1, B2, G |
|---|---|---|
| values of $R$, $Z$ on the surfaces | rounding | rounding |
| values of $\lambda$ and the Jacobian | $2 \cdot 10^{-11}$ interior, $1.6 \cdot 10^{-9}$ at the edge | rounding |
| first radial derivative of $R$, $Z$ | $2 \cdot 10^{-8}$ interior, $4.5 \cdot 10^{-7}$ at the edge | below $3 \cdot 10^{-13}$ |
| first radial derivative of $\lambda$ and of the Jacobian | $2 \cdot 10^{-8}$ interior, $1.8 \cdot 10^{-6}$ at the edge | below $3 \cdot 10^{-13}$ |
| second radial derivative of $R$, $Z$ | $6 \cdot 10^{-5}$ interior, $5.6 \cdot 10^{-4}$ at the edge | below $8 \cdot 10^{-11}$ |
| second radial derivative of $\lambda$ | $6 \cdot 10^{-5}$ interior, $1.2 \cdot 10^{-3}$ at the edge | below $8 \cdot 10^{-11}$ |
| third radial derivative of $R$, $Z$, $\lambda$, which the energy does not read (section 1.5) | 0.12 interior, 0.32 at the edge, of values up to 27 | below $2 \cdot 10^{-11}$ (20 elements), $1.5 \cdot 10^{-8}$ (100 elements) |
| poloidal flux | trapezoid, $1.2 \cdot 10^{-9}$ on the circular twin, proportional to $h^2$ | GVEC's $\chi$ |
| Jacobian's mode list | cut at a tolerance in the script | B1, B2: cut at a tolerance; G: none, computed at the grid points |
| noise in the values, multiplied by PB3D's splines | rounding, not visible (section 1.5) | PB3D does not differentiate the data |
| more accuracy by | more surfaces: factor 8 (first derivative) and 4 (second) for twice the surfaces, with more memory and time in PB3D | not needed |

Read against the requirements of section 1.5:

- R1 (noise) is met by every route, as long as the scripts only evaluate GVEC's B-splines and Fourier series directly and write double precision.
- R2 (smooth errors below $10^{-7}$ in values and first derivatives) is met by the routes B1, B2 and G everywhere. Route A meets it in the interior. At the edge its first derivatives have errors of $4.5 \cdot 10^{-7}$ ($R$, $Z$) and $1.8 \cdot 10^{-6}$ ($\lambda$, Jacobian) at 319 surfaces. From these values, those at 639 surfaces and the factor 8 for twice the surfaces, route A meets R2 at the edge from about 530 surfaces on for $R$ and $Z$ and from about 830 surfaces on for $\lambda$ and the Jacobian.
- R3 (second derivatives) is the open one. The routes B1, B2 and G have errors of $10^{-10}$ there. Route A has $6 \cdot 10^{-5}$ to $10^{-3}$ at 319 surfaces, the same as PB3D has for VMEC files now.

### 4.2 Other criteria

| | A | B1 | B2 | C, added to A or B | G |
|---|---|---|---|---|---|
| steps for the user after the GVEC run | script, PB3D | PB3D | script, PB3D | as A or B | script, PB3D |
| Fortran that changes | none | new reader, build | `read_VMEC` | six routines on the path of every VMEC run | new reader, 64 blocks, 22 tests, Jacobian; six routines for style 3 |
| can results on VMEC or HELENA files change | no | no | no, if the branch without extra variables is untouched | yes, at rounding level, unless a switch keeps the cylindrical statements (D4) | no |
| Fortran that exists twice | none | none | two ways to fill the arrays in `read_VMEC` | none | none with the arrays of `VMEC_vars`; 591 lines with arrays of its own |
| new dependency of PB3D's build | none | GVEC's library | none | none | none |
| needed to run the scripts | pyGVEC | no script | pyGVEC | as A or B | pyGVEC |
| GVEC in CI | no, a file is the fixture | yes, or the tests are skipped | no | no | no |
| GVEC's maps | cylindrical | cylindrical | cylindrical | all | all |
| the file outside PB3D | other VMEC-based tools read it | no file | other VMEC-based tools read it | no longer a VMEC file | PB3D only |
| GVEC and VMEC equilibria through the same statements of PB3D | yes | no: derivatives | no: derivatives | as A or B, plus the frame | no: derivatives, Jacobian, frame |
| depends on GVEC interfaces that are not published | no | yes, `MODgvec_py_state` | only if the third derivative is taken from the state file | no | only if the third derivative is taken from the state file |
| first eigenvalue of a GVEC equilibrium after (guess) | 1 week | 2 to 3 weeks | 2 to 3 weeks | 2 to 3 weeks more | 2 to 3 weeks |

### 4.3 What is feasible

- No route fails a requirement that is known. Route A is PB3D's present way of reading an equilibrium, applied to data that are smoother than a VMEC file's, because they are evaluations of one B-spline and not the result of a radial finite-difference scheme.
- Route A's accuracy at the edge is the weakest point of the comparison. If stage S3 shows that the eigenvalue feels errors of $10^{-3}$ in the second radial derivatives at the edge at the level of the target accuracy, route A needs many more surfaces than 319, or it is not good enough; the same then holds for PB3D's results on VMEC files.
- Not feasible, in any route:
  - profile tables with about 100 points in the GVEC run, as `pygvec convert-wout` writes them: they change $\iota'$ by 6 to 12 % of the local shear on the circular twin (stage S2);
  - the profiles from the table in the state file, for the same reason;
  - a path through other angles (PEST, Boozer) or through GVEC's plot files: both hold resampled values with errors far above R1;
  - text or single precision between the script and PB3D.
- Feasible but not followed: PB3D reading GVEC's files with Fortran of its own (end of section 3).

### 4.4 Recommendation

This is my reading of the comparison; the decision is D5 in section 9.

- Stage S0 (the error for an unknown `eq_style`) and stage S3 (how accurate the second radial derivatives must be) come first. Neither needs GVEC, and S3 turns the open requirement R3 into a number that weighs on the choice between route A and the routes with exact derivatives.
- If only GVEC equilibria in cylindrical coordinates are wanted and S3 shows that errors of $10^{-3}$ in the second derivatives at the edge do not move the eigenvalue at the target accuracy: route A. It changes no Fortran, and it compares GVEC and VMEC through the same statements.
- If S3 shows that such errors matter, or if GVEC's other maps are wanted (decision D2): route G, or route B2 with route C. Both need about the same amount of new Fortran. Route G leaves the statements for VMEC files untouched and needs no wout file with variables that VMEC does not have; it costs the edits of 64 blocks and a third style to maintain. Route B2 with C keeps `eq_style` 1 and edits `read_VMEC` and the six routines, behind a switch or for all files.
- Route B1 only if running PB3D on GVEC's files without a script is worth a compiled dependency on a GVEC module that is not a published interface. I would not choose it.
- Without the result of S3, and because the maps other than the cylindrical one are the reason GVEC exists, I expect the choice to fall on route G. The first week of work (S0, S1, S2, S3 and the script's evaluation of the state) is the same for the routes A, B2 and G, so the choice can wait for S3.

## 5. Stages and tests shared by the routes

Each stage ends with a gate. After each stage: update this plan with what was learned, comment the gate result on #104, commit and push, never merge. Heavy runs go to Snellius as parallel Slurm jobs.

### S0. An error for an unknown `eq_style`

Every one of the 65 `select case (eq_style)` blocks of section 1.4 gets a `case default` that stops with the message "No equilibrium style associated with ...", in the form PB3D uses for its other styles. This is independent of GVEC and of the route: its own issue and pull request. Gate: every regression and full-stack test gives its present result.

### S1. GVEC installed on the Mac and on Snellius

`pip install gvec==1.5.1` in PB3D's uv environment as an optional dependency group `gvec` (as `vmec = ["vmecpp"]`). The package is a source distribution: it compiles GVEC and needs a Fortran compiler, CMake, LAPACK and netCDF-Fortran (Homebrew on the Mac, modules on Snellius). Gate: GVEC's elliptic-tokamak tutorial runs and `pygvec info` prints its state on both machines.

### S2. A GVEC twin of the circular fixture, and what a GVEC run for PB3D needs

`pygvec convert-params --vmec` cannot be used: GVEC 1.5.1 accepts only `power_series` profiles (`util/params.py:455-485`), while `--export_HEL` writes `cubic_spline` profiles (`eq_ops.f90:1949-1985`), and the RBS and ZBC columns of these inputs make it choose a stellarator-asymmetric run. Instead, `pygvec convert-wout tests/fixtures/wout_circ_eps033_q104.nc` gives a GVEC parameter file with the twin's boundary, the profiles as `interpolation` tables, and `_cos_`/`_sin_` symmetry. The parameter file has no `whichInitEquilibrium`, so GVEC starts from its default initial state and solves its own equilibrium, in stages that raise `sgrid_nElems` and the Fourier resolution (`--param` values are written into the file).

What the GVEC run needs, for every route:

- Profile tables that are dense enough. `convert-wout` writes the profile tables at 101 values of $s$ (`np.linspace(0, 1, 101)` in `convert_wout.py`). Resampling the twin's profiles there changes $q$ by $7 \cdot 10^{-7}$ but $\iota'$ by up to $2.2 \cdot 10^{-4}$ near $s = 0.59$, where the shear is only between $-0.0035$ and $-0.0018$: 6 to 12 % of the local shear. With 319 values the change in $\iota'$ is $1.9 \cdot 10^{-5}$, with 639 values $4.9 \cdot 10^{-6}$; for cbm18a cut, 101 values change $p'$ by $5 \cdot 10^{-4}$ and $\iota'$ by $1.8 \cdot 10^{-3}$ of their maxima. So the tables in the parameter file are replaced by `iotaf` and `presf` on all the twin's surfaces (a small script after `convert-wout`). How GVEC interpolates such a table is checked by evaluating $\iota$, $\iota'$, $p$, $p'$ of the solved state against the twin's.
- B-splines of degree 5 for $X^1$, $X^2$ and $\lambda$, GVEC's default. With GVEC's continuity of one less than the degree, the second radial derivative, the highest that PB3D's energy reads, is then smooth, and so is the function that the cubic splines of route A interpolate. The scripts read the degrees from the header of the state file and stop below degree 4.
- An equilibrium that is converged in GVEC's resolution: PB3D assumes force balance. The state is computed at three resolutions (radial elements and Fourier modes), and the tests below are made at all three.

The twin's VMEC input prescribes $\iota$ (`NCURR = 0`) and the GVEC run prescribes the same table, so $q$ agrees by construction and is no check. Gate: the force residual decreases with resolution, and at equal $s$ the solved GVEC surfaces agree with the VMEC twin (axis position, $R$ and $Z$ on the surfaces $s$ = 0.25, 0.5, 1, $B_{00}$ on the axis) to within the VMEC twin's own resolution error, estimated from VMEC++ runs with 159 and 319 surfaces.

### S3. How accurate the data must be: measured on VMEC cases, without GVEC

This stage turns the requirements of section 1.5 into numbers. It uses the circular VMEC twin (fixed and free boundary) and the cbm18a cut, whose eigenfunction sits at the edge, with `tol_SLEPC` = $10^{-12}$. The Fortran changes are switches in a test build and are not merged.

- S3a, the third radial derivative. A switch sets the third radial derivatives of $R$, $Z$ and $\lambda$ to zero after `read_VMEC`. If the eigenvalues do not change, the energy does not read them (section 1.5), and no route has to compute them with care.
- S3b, the first and the second radial derivative. A switch multiplies the radial derivatives of one order by $1 + \delta\, w(s)$: the first derivatives of $R$, $Z$, $\lambda$ and the Jacobian for $\delta$ from $10^{-8}$ to $10^{-4}$, and the second derivatives of $R$, $Z$ and $\lambda$ for $\delta$ from $10^{-6}$ to $10^{-2}$. Each is done once with $w = 1$ and once with a weight $w$ that is 1 in the outermost 5 % of the radius and 0 elsewhere. The eigenvalue's change per $\delta$ checks requirement R2 for the first derivatives and gives requirement R3 for the second. With the errors of section 4.1 it says how many surfaces route A needs for a given accuracy of the eigenvalue, and what the same errors cost PB3D's results on VMEC files today.
- S3c, noise. A script adds random relative noise of $10^{-14}$ to $10^{-8}$ to `rmnc`, `zmns`, `lmns` and `gmnc` of the wout file; no Fortran is involved. The eigenvalue's change against the size of the noise is the measured form of requirement R1.

Gate: the three results are on #104, with the number of surfaces that route A needs for an eigenvalue accurate to $10^{-5}$ and to $10^{-6}$.

### S4. The cases with finite pressure and in 3-D

Three cases, run side by side on Snellius, each with test T2 first and then the eigenvalues:

- cbm18a cut at the pedestal foot (axisymmetric, finite pressure; PB3D on HELENA gives -0.04012 at n = 10 fixed boundary, `tests/reference/cbm18a_cut/README.md`; VMEC twin from `input.cbm18a_cut` with 639 surfaces);
- the rippled twin, N = 5, $\delta$ = 0.02 (3-D; PB3D free boundary -0.4188 against TERPSICHORE's family -0.41903);
- W7-X (PB3D fixed-boundary reference 6.62833e-2 on Snellius). GVEC's repository also holds a W7-X wout file with the GVEC parameters its developers use on it (`test-CI/examples/w7x_from_vmec_*`, `wout_d23p4_tm.nc`), a starting point for the GVEC resolution.

For each: GVEC from `pygvec convert-wout` of the VMEC file, with the profile tables replaced as in S2, and its own solve; PB3D on the GVEC equilibrium against PB3D on VMEC's wout, at two resolutions of each code. `convert-wout` always prescribes $\iota$. For a VMEC run with a prescribed current (W7-X, if its input has `NCURR` = 1), GVEC prescribes VMEC's converged `iotaf`; both codes then solve the same problem up to discretisation. GVEC's current constraint would need `I_tor` and `picard_current` added by hand and is not needed for this comparison. Gate: the difference between the two eigenvalues decreases when the resolutions grow; the numbers go on #104.

### Tests used by every route

T1, the file against GVEC in space, with no convention in between (routes A, B2 and G, which have a file). At a few thousand points $(s, \theta, \varphi)$, on the file's surfaces and at random angles, a check script computes from the file's arrays the position, the Cartesian vector $B = B^\theta\, \partial_\theta x + B^\varphi\, \partial_\varphi x$ (with $B^\theta$ and $B^\varphi$ from the fluxes, $\lambda$ and the Jacobian) and the pressure, and compares them with pyGVEC's `pos`, `B` and `p` at $(\sqrt s, \theta, -\varphi)$. A vector in space has no sign convention, so this catches a mirrored equilibrium, a wrong sign of a flux or of $\lambda$, and a missing $2\pi$, none of which an eigenvalue comparison is certain to show (the mirror image of an equilibrium has the same eigenvalues). Gate: the positions agree to $10^{-13}$ of the major radius and the pressure to $10^{-13}$ relative. $B$ agrees to $10^{-11}$ relative for the files of B2 and G. For the wout file of route A the check takes $\lambda$ and the Jacobian from the half mesh with the spline that `read_VMEC` uses, and $B$ agrees to $10^{-8}$ at 319 surfaces, 16 times better at 639; the check's formulas are first calibrated on the VMEC++ fixture, whose `bsupumnc` and `bsupvmnc` they must reproduce.

T2, PB3D's stored equilibrium against GVEC, point by point (every route), on the GVEC twin of S2 with the deck of `regression_circ_eps033_vmec_twin`. A script reads from `PB3D_out.h5` the coordinates of every point of the equilibrium grid and PB3D's values there (section 1.1), and evaluates the same quantities with pyGVEC at $(\sqrt{s}, \theta, -\varphi)$: $\lvert B \rvert$ (in PB3D $\lvert B \rvert^2 = g_{\parallel\parallel}/\mathcal J^2$ from `g_FD` and `jac_FD`, because $B = e_\parallel/\mathcal J$ in its flux coordinates), the normal and the geodesic curvature from GVEC's `kappa_B` and $\nabla\rho$ in PB3D's definitions, the parallel current `sigma` from GVEC's $J \cdot B/\lvert B \rvert^2$, and the shear `S` from pyGVEC's first derivatives of the metric. The curvature, the current and the shear hold second radial derivatives, so this tests the whole chain from GVEC to the quantities in the energy against values that no VMEC file can give, since VMEC has no exact evaluation between its surfaces. `sigma` changes sign under a mirror reflection, so the test also sees a mirrored reading. Gate for the routes B1, B2 and G: agreement to $10^{-9}$ relative at any number of surfaces. Gate for route A: the differences are of the size of section 4.1 ($6 \cdot 10^{-5}$ in the interior, $10^{-3}$ at the edge at 319 surfaces) and fall by a factor of 4 for twice the surfaces, with no difference that stays.

T3, eigenvalues (every route), with the deck of `regression_circ_eps033_vmec_twin` and `tol_SLEPC` = $10^{-12}$. The eigenfunction is plotted, and the radial range of the deck does not cut it. The VMEC twin agrees with HELENA to $1.1 \cdot 10^{-6}$ relative (+7.84303e-2 against +7.84304e-2) at this deck and resolution. Gate: the GVEC twin's fixed-boundary eigenvalue approaches the VMEC twin's as GVEC's resolution grows, and lies within $10^{-5}$ relative of it at the highest of the three resolutions of S2; the difference to HELENA is reported next to it. The free-boundary eigenvalue lies within $10^{-4}$ relative of the VMEC twin's -0.43901. For routes with a number of surfaces `ns` the table has two values of `ns`, so that it shows whether GVEC's resolution or `ns` limits the agreement.

T4, from VMEC through GVEC without a solve (every route): `wout_circ_eps033_q104.nc` $\to$ `pygvec convert-wout --param sgrid_nElems=...` $\to$ GVEC state $\to$ the route's path into PB3D. Gate: PB3D's fixed-boundary eigenvalue approaches that of PB3D on the original file (n = 10, m = 11, +7.84303e-2) as `sgrid_nElems` grows: three values of `sgrid_nElems`, the difference falling, and below $10^{-5}$ relative at the largest. GVEC's fit of a wout file keeps a floor in $\lambda$ (`InitVMEC` fits $\lambda$ on the half mesh, moves it to the full mesh and fits again, `vmec.F90:365-398`) and extrapolates every mode to the axis by a parabola (`vmec.F90:335`); if the eigenvalue stops approaching, these two are the first places to look.

## 6. Stages of route A

A1. The converter `gvec_to_wout.py`. Gates:

- A1a: test T1.
- A1b, identities inside the written file: $d\,$`phi`$/ds$ = `phipf`; `chipf` = `iotaf` $\cdot$ `phipf`; $d\,$`chi`$/ds$ = `signgs` $\cdot$ `chipf`; `gmnc` negative and equal to $R\,(R_\theta Z_s - R_s Z_\theta)$ built from the file's own `rmnc` and `zmns`; `bmnc`, `bsup*` and `bsub*` equal to what the file's other variables give. On the VMEC twin `gmnc` agrees with $R\,(R_\theta Z_s - R_s Z_\theta)$ to $3 \cdot 10^{-9}$ relative in the interior and $5 \cdot 10^{-7}$ at $s = 0.017$. On the converted file the radial derivatives in these identities are differences between surfaces, so outside the axis region (where the coefficients of odd $m$ behave like $\sqrt s$) the residuals must fall as $h^2$ when `ns` grows. This gate covers the variables that PB3D does not read and other tools do.
- A1c: the file is read back and equals what the script held in memory, bit for bit; LIBSTELL reads it, and `tests/regression/check_vmec_read.py` passes (mode lists and PB3D's copy of the file).

A2. PB3D on the GVEC twin: tests T2, T3 and T4. In T4 only `rmnc`, `zmns` away from the axis, `iotaf`, `presf` and `phi` of the file after the round trip can converge to the original file when `sgrid_nElems` grows; `gmnc`, `bmnc` and the $B$ components are GVEC's evaluations, which differ from VMEC's own finite-difference values by VMEC's discretisation error.

A3. The cases of stage S4.

A4. The number of surfaces. At a fixed GVEC state, the differences of test T2 and PB3D's eigenvalue at several `ns`, with PB3D's memory and time. With the result of S3 this answers decision D1: whether an `ns` with acceptable cost reaches the target accuracy.

A5. Regression test and documentation:

- `regression_circ_eps033_gvec_twin`: the converted wout of A2 as a compressed fixture in `tests/fixtures`, so that CI needs no GVEC; the eigenvalues anchored as in the VMEC twin test.
- The identity check of A1b as a script in `tests/regression`, run on the fixture; the checks against GVEC (T1, T2, T4) only when the `gvec` group is installed, skipped with a note otherwise, as `check_vmec_read.py` is without h5py.
- README entry; a page on the site (PB3D.github.io) on how to run PB3D on a GVEC equilibrium: `pygvec run`, `gvec_to_wout.py`, PB3D. No version change, since no Fortran changes.

## 7. Stages of routes B and C

### Route B1

B1.1. GVEC's library built with PB3D's compilers on the Mac, on Snellius and in CI, without MPI; PB3D's CMake finds it as an optional package. Gate: a small Fortran program linked against it loads the state of S2 and returns pyGVEC's values at the same points to rounding.

B1.2. The reader, on the root process, and the recognition of the two files. Gate: every regression and full-stack test gives its present result in a build without GVEC and in a build with it; tests T2, T3 and T4.

B1.3. The cases of S4, the regression test (skipped in a build without GVEC) and the documentation, with a new version number.

### Route B2

B2.1. The converter of route A writes the extra full-mesh variables, with the decision on the axis surface. Gate: test T1 with the extra variables, and each stored derivative of order $k+1$ agrees with a finite difference of the stored order $k$ at the rate of that difference.

B2.2. `read_VMEC` uses them when present, says in its output which way it took, and `misc_in_V` records it. Gate: every regression and full-stack test gives its present result; tests T2, T3 and T4.

B2.3. The cases of S4 again. Gate: the eigenvalues of routes B2 and A agree to within the error that stage A4 measured for route A.

### Route C

C1. The frame in PB3D: reading, evaluation with derivatives, storage in `PB3D_out.h5`, and the six routines of section 1.2 written for a general frame, with the cylindrical frame for files that hold none. Gate: every regression and full-stack test gives its present result. Whether VMEC files then go through the general formulas or keep the present cylindrical ones is decision D4: with the general formulas there is one code path, and the results may change at the level of rounding.

C2. Frame invariance, without GVEC. A script rewrites a VMEC wout file in another frame that keeps $\varphi$: a shifted one, $X_0 = (a_0 + a_1\cos N_\mathrm{fp}\varphi)\, e_R(\varphi)$, and a tilted one, $N$ and $B$ turned in the plane $\varphi$ = constant by an angle $\beta_1 \sin N_\mathrm{fp}\varphi$. The surfaces, the angles $(s, \theta, \varphi)$, $\lambda$ and the Jacobian are the same; only $X^1$, $X^2$ and the frame differ. Gate: on the circular twin, the rippled twin and W7-X, PB3D gives the same fixed- and free-boundary eigenvalues, the same stored quantities of test T2, and the same Cartesian positions as with the plain file, to rounding for the shifted frame and to the Fourier truncation of $X^1$, $X^2$ for the tilted one. This tests all six routines on existing cases.

C3. GVEC, the same boundary with map 1 and with a G-Frame whose $\zeta$ is the cylindrical angle. `pygvec load-quasr -f` builds a G-Frame from a boundary given as points in space, and keeps the toroidal angle of the data (GVEC's QUASR page: for surface data computed from cylinder coordinates "the construction will not change this"); the boundary of the rippled twin, sampled at cylindrical angles, is such data. Gate: test T2 on both, and equal eigenvalues to within the difference of the two GVEC solutions, which falls with resolution.

C4. A G-Frame whose $\zeta$ is not the cylindrical angle: the G-Frame examples of GVEC's repository (`test-CI/examples/axisNB_*`, among them a knot and a quasi-helical case), a QUASR configuration (`pygvec load-quasr`), then the figure-8 configuration of Plunk et al. (2025). Gate: test T2. The eigenvalues are reported with a study of how they change with $n$ and with the parametrisation of the frame, the open point of the other maps.

## 8. Stages of route G

G1. The script `gvec_to_pb3d.py`. Gates:

- G1a: test T1.
- G1b, identities inside the file: the derivative of the poloidal flux equals $\iota$ times the derivative of the toroidal flux, to rounding; each stored derivative of order $k+1$ agrees with a finite difference of the stored order $k$ at the rate of that difference.
- G1c: the file is read back and equals what the script held in memory, bit for bit.

G2. Style 3 in PB3D, cylindrical map: `read_GVEC`, the recognition of the file, the edits of the 64 blocks and of the tests for style 1 (section 1.4), and the Jacobian at the grid points. Until stage G3, `read_GVEC` stops with a message when the file's map is not the cylindrical one. Gates:

- G2a, the existing tests. Every regression and full-stack test gives its present result, bit for bit where the test compares that way: for a VMEC or HELENA file PB3D runs the statements it ran before.
- G2b, a VMEC equilibrium through style 3, without GVEC. A script writes the arrays that `read_VMEC` made of the circular VMEC twin into a file of route G. PB3D stores these arrays in `PB3D_in.h5`, normalized and for the radial range of the run; the script undoes the normalization, and the stage settles how the radial range is handled. The two runs then differ in one thing: style 1 reads VMEC's Jacobian from the file, and style 3 computes the Jacobian from $R$ and $Z$. On the VMEC twin the two Jacobians agree to $3 \cdot 10^{-9}$ in the interior and to $5 \cdot 10^{-7}$ near the axis (gate A1b). So the stored metric of the two runs must agree to rounding, and the stored Jacobian and the eigenvalues to $10^{-6}$ relative. This tests the reader and every edited block with no change of the equilibrium.
- G2c, G2d, G2e: tests T2, T3 and T4.

G3. The frame: GVEC's other maps.

- G3.1. The frame in the file and in PB3D: the script writes $X_0$, $N$, $B$ as Fourier series in $\varphi$ for any map (obtained as in section 2.1); PB3D reads them, evaluates them and their first three derivatives at the grid's $\varphi$ values, stores them in `PB3D_out.h5`, and uses the formulas of route C in the six routines of section 1.2 for style 3. The normalization length is the length of the magnetic axis divided by $2\pi$ (GVEC's `r_major`). Gate: the tests with VMEC and HELENA files give their present results bit for bit, and the GVEC cases of G2 and S4, now with the cylinder as their frame, give their earlier results to rounding.
- G3.2. Frame invariance, without GVEC: stage C2, with the script of gate G2b writing the VMEC twin's arrays in the shifted and in the tilted frame.
- G3.3 and G3.4: stages C3 and C4, with the file of route G.

G4. The cases of stage S4, the regression test and the documentation:

- `regression_circ_eps033_gvec_twin`: the GVEC file of G2 as a compressed fixture in `tests/fixtures`, so that CI needs no GVEC; the eigenvalues anchored as in the VMEC twin test.
- The identity check of G1b as a script in `tests/regression`, run on the fixture; the checks against GVEC (T1, T2, T4) only when the `gvec` group is installed, skipped with a note otherwise, as `check_vmec_read.py` is without h5py.
- README entry with a new version number, since the Fortran changes; a page on the site (PB3D.github.io) on how to run PB3D on a GVEC equilibrium: `pygvec run`, `gvec_to_pb3d.py`, PB3D.

## 9. Decisions for Toon

- D0: where the script lives (`gvec_to_wout.py` or `gvec_to_pb3d.py`): (a) in PB3D's `scripts/`, (b) in GVEC as a `pygvec` subcommand next to `to-cas3d`, which Toon would propose upstream (MPCDF GitLab or the GitHub mirror), (c) first (a), later (b).
- D1: after stage S3, and with route A after stage A4: whether the accuracy of radial derivatives taken by splines is enough, and at which number of surfaces.
- D2: whether GVEC's maps other than the cylindrical one are wanted: route C, or stage G3 of route G.
- D3: whether the GVEC twin becomes a regression test or stays a separate study outside this repository.
- D4, only with route C: whether VMEC files use the general frame formulas (one code path, results may change at rounding level) or keep the cylindrical formulas behind a switch.
- D5: the route (section 4). On 27 November 2026 Toon said that he leans towards a GVEC reader as a third equilibrium style, route G.
- D6, only with route G: whether the arrays and routines of the VMEC path keep their names with `V` when they also hold a GVEC equilibrium, or get other names. This plan keeps the names: other names change about 40 `use` statements and the names of the variables in `PB3D_in.h5` and `PB3D_out.h5`, for no change in what the code does.

## 10. Sources

- GVEC JOSS paper, doi:10.21105/joss.09670; flexible frame, doi:10.1088/1361-6587/adba11.
- GVEC documentation, https://gvec.readthedocs.io: parameter list, theory, coordinate conventions, quantities, pygvec CLI (also the `develop` version), QUASR and G-Frame, CAS3D and other interfaces.
- GVEC 1.5.1 source, public mirror https://github.com/gvec-group/gvec: `src/functionals/mhd3d/mhd3d.F90`, `src/vmec/vmec.F90`, `src/vmec/vmec_readin.F90`, `src/hmap/hmap_RZ.F90`, `hmap_knot.F90`, `hmap_frenet.F90`, `hmap_axisNB.F90`, `src/base/fbase.F90`, `src/restart/restart.F90` (the state file), `src/pygvec/state.F90` (`MODgvec_py_state`), `CMakeLists.txt` (`libgveclib`), `python/gvec/quantities.py`, `python/gvec/scripts/convert_wout.py`, `python/gvec/scripts/cas3d.py`, `python/gvec/util/params.py`, `python/gvec/core/state/` (`_evaluate.py`, `_class.py`, `_binding.py`); test data under `test-CI/`.
- LIBSTELL (STELLOPT repository): `LIBSTELL/Sources/Modules/read_wout_mod.f90` (`read_wout_nc`) and `LIBSTELL/Sources/Ezcdf/ezcdf_inqvar.f90`.
- PB3D: issue #65 and `read_HEL` (mapping files with 8 digits), PR #34 (the half mesh in `read_VMEC`), the source files named in sections 1.1 to 1.5.
- The error measurements of sections 1.5, 3 and 4.1 are made with numpy and scipy on the test functions named there (cubic splines: `scipy.interpolate.CubicSpline`; B-splines: `scipy.interpolate.BSpline` of degree 5 with clamped ends, holding $s^2 - 0.3\,s$ exactly) and on the circular VMEC twin of `tests/fixtures`.
