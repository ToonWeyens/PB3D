# Survey B: 3-D linear ideal-MHD stability codes and published 3-D benchmarks

Prepared 31 October 2026 for the PB3D fixed-boundary 3-D benchmark. Scope: stellarator and 3-D tokamak
linear stability codes, their intercomparisons, reference cases, the theory of mode families, and the
infinite-n ballooning limit.

How each entry was checked:

- **[read]** I read the web page, abstract page or PDF text myself.
- **[library]** I read the text of the paper in Toon's own library (`~/Code/papers`, the `texts` table).
- **[repo]** I read it in the PB3D worktree.
- **[snippet]** I only saw a search-engine summary. Treat numbers and titles from these as leads to check.
- **[inference]** My own reasoning from the cited material, not a statement in any paper.

All web requests were anonymous (no name, address or token sent). OSTI (osti.gov) refused connections
during this session, several IOP, Elsevier, MPG-PuRe and EPS pages returned 403 or a bot wall, and the
W7-X web service host does not resolve from outside IPP.

---

## 1. Three-dimensional linear stability codes

### Summary table

| Code | Physics | Boundary | Mode families | Equilibrium input | Availability |
|---|---|---|---|---|---|
| CAS3D | Global linear ideal MHD, full compressible δW | Fixed, or free with Green's-function vacuum (no wall needed) | Couples the toroidal harmonics of one family | VMEC, mapped to Boozer | IPP (Greifswald/Garching); no public release found |
| TERPSICHORE | Global linear ideal MHD; VVBAL module for infinite-n ballooning | Fixed, or vacuum as pseudo-plasma inside a conducting wall | One family per run | VMEC wout, own Boozer mapping | **Open source since 2025**, EPFL GitLab, Apache 2.0 |
| CASTOR3D | Linear full MHD, ideal and resistive, viscosity, flow, diamagnetic drift; STARWALL vacuum and resistive walls | Free boundary with 3-D walls | One family per run; reduced ansatz at high n* | NEMEC/VMEC or GVEC | IPP Garching; no public release found |
| COBRA / COBRAVMEC | Infinite-n ideal ballooning equation, one field line at a time | Not applicable (local) | Not applicable | VMEC wout | Open source in STELLOPT (GitHub, MIT) |
| DESC ballooning objective | Infinite-n ideal ballooning; max over several field lines | Local | Not applicable | DESC equilibrium | Open source (GitHub PlasmaControl/DESC) |
| ballopt (Gaur) | Infinite-n ballooning with adjoint gradient, uses SIMSOPT + VMEC2000 | Local | Not applicable | VMEC | GitHub, Zenodo DOI |
| JOREK3D | Nonlinear reduced MHD (stellarator model) | Not checked | Simulates several n of a family together | GVEC | JOREK consortium, by contact |
| M3D-C1 (stellarator) | Nonlinear extended MHD | Fixed in the published W7-X runs | All n in the run | Initialised from free-boundary VMEC | Public GitHub, BSD-style licence |
| PB3D | High-n (peeling-ballooning ordering) ideal MHD, one n, field-line based | Fixed or free (BEM vacuum) | One n only | VMEC, HELENA | This code |

### 1.1 CAS3D (Code for the Analysis of the MHD Stability of 3-D Equilibria)

- C. Schwab, "Ideal magnetohydrodynamics: Global mode analysis of three-dimensional plasma configurations",
  *Phys. Fluids B* **5**, 3195–3206 (1993). [snippet: citation confirmed by two independent reference lists;
  also cited in Lazerson et al. 2026, read in library]
- C. Nührenberg, "Global ideal magnetohydrodynamic stability analysis for the configurational space of
  Wendelstein 7-X", *Phys. Plasmas* **3**, 2401–2410 (1996). [library: cited as ref. 16 of the PB3D JCP
  paper and ref. 3 of the NCSX review slides; content snippet only] https://www.osti.gov/biblio/284302
- Note: CAS3D computes the potential and kinetic energy by finite elements radially and Fourier in both
  angles, and "uses Green's function method to solve the vacuum problem and thus can evaluate stability
  without a conducting wall" (G.Y. Fu, NCSX Physics Validation Review, 2001). [read]
  https://ncsx.pppl.gov/pvr/PVR_Presentations/Fu_PVR.pdf
- A variant with non-equidistant radial grid, CAS3DN, was used for the CASTOR3D benchmark (section 2.3). [read]
- CAS3D-K is a kinetic (energetic-particle, Alfvén-eigenmode) extension listed by IPP Greifswald. [snippet]
- What it gives the benchmark: the reference global ideal code for W7-X; all W7-X β limits in the
  literature come from it. It computes one mode family at a time and resolves the coupling of n with
  n ± N_fp, so a PB3D-versus-CAS3D comparison measures what PB3D leaves out. Access only through IPP.
- Correction for the PB3D thesis (§1.4, [library]): it calls CAS3D "an extension of CASTOR". CAS3D is
  Schwab's independent code; the CASTOR extension is CASTOR3D (section 1.3). [read: CASTOR3D preprint]

### 1.2 TERPSICHORE (with VVBAL)

- D.V. Anderson, W.A. Cooper, R. Gruber, S. Merazzi, U. Schwenn, "TERPSICHORE: A Three-Dimensional Ideal
  Magnetohydrodynamic Stability Program", in *Scientific Computing on Supercomputers II*, Springer (1990)
  159–174, doi:10.1007/978-1-4613-0659-7_8. Also *Int. J. Supercomput. Appl.* **4**, 34 (1990). [read: EPFL page; library]
- Code page: https://www.epfl.ch/research/domains/swiss-plasma-center/research/theory/codes/research_theory_codes_terpsichore/ [read]
- Repository: https://gitlab.epfl.ch/spc/public/terpsichore [read via the GitLab API]: public, created
  May 2025, last activity September 2026, Apache 2.0 (tag 1.1 of 1 June 2025 is "the first open-source
  TERPSICHORE version, based on T. Cooper's version"). Build: `make` in `srcs/`. Contains a fixed-boundary
  benchmark `QAS3_fx_benchmark/` with `input.QAS3_fx_cur75_bench`, a text wout file and expected output,
  and a regression test `tests/QAS3_fx_cur_bench_n1/`. The README names a contact at SPC-EPFL.
- What it solves: global linear ideal MHD; hybrid finite elements radially, Fourier in both angles; maps
  the VMEC solution to Boozer coordinates; outputs normal and binormal displacement (EPFL page) [read].
  "Treats vacuum as a pseudo-plasma" so a conducting wall is needed (Fu 2001) [read]. Lazerson et al.
  2026 put a conformal wall 30 % further out to avoid wall stabilisation [library].
- VVBAL module: solves a modified ballooning equation with the multiplier (1 − λ) on the pressure term, so
  only its marginal points match COBRA, not its eigenvalues (Fu 2001; Sanchez et al. 2001) [read; library].
- What it gives the benchmark: the only global 3-D ideal code that can now be downloaded and run without
  an agreement. It reads VMEC wout files, runs one mode family per call, and has a fixed-boundary
  benchmark case. This is the most practical reference for PB3D's fixed-boundary 3-D runs.

### 1.3 CASTOR3D (with STARWALL)

- E. Strumberger, S. Günter, "CASTOR3D: linear stability studies for 2D and 3D tokamak equilibria",
  *Nucl. Fusion* **57**, 016032 (2017), doi:10.1088/0029-5515/57/1/016032. Preprint EUROfusion
  WP15ER-PR(16)16343: https://scipub.euro-fusion.org/wp-content/uploads/eurofusion/WP15ERPR16_16343_submitted-2.pdf [read]
- E. Strumberger, S. Günter, *Nucl. Fusion* **59**, 106008 (2019), doi:10.1088/1741-4326/ab314b
  (QA stellarator, parallel viscosity, flow, resistive walls). [read: reference list]
- E. Strumberger, S. Günter, W7-X team, "Linear, resistive stability studies for Wendelstein 7-X-type
  equilibria with external current drive", *Nucl. Fusion* **60**, 106013 (2020), doi:10.1088/1741-4326/aba9ea. [read: reference list]
- E. Strumberger, S. Günter, K. Lackner, J. Puchmayr, "CASTOR3D: linear magnetohydrodynamics and
  diamagnetic drift effects", *J. Plasma Phys.* **89**, 905890309 (2023), doi:10.1017/S0022377823000508. [read: reference list]
- E. Strumberger, J. Puchmayr, F. Hindenlang, "Linear stability studies for a quasi-axisymmetric
  equilibrium including plasma flow, drift and viscosity effects", *J. Plasma Phys.* **91**(6), E157
  (2025), doi:10.1017/S0022377825101013, open access.
  https://www.cambridge.org/core/product/3EBD9D62179F8287E3219E82CDA2F567/core-reader [read]
- What it solves: CASTOR3D is a hybrid of CASTOR_3DW (an extension of CASTOR and CASTOR_FLOW) and STARWALL;
  ideal and resistive; plasma inertia and wall resistivity in one extended eigenvalue problem; general 3-D
  flux coordinates (NEMEC, Boozer, straight-field-line). [read: preprint] The 2025 version adds drift and
  E×B effects, uses GVEC equilibria, labels a mode by its dominant harmonic n*, keeps only harmonics of one
  family (for a 2-period QA: odd n or even n), and switches to a reduced ansatz at high n* that cuts the
  matrix size by about four. [read]
- What it gives the benchmark: the global code that has published **peeling-ballooning modes in 3-D
  geometry** — 3-D tokamaks with magnetic perturbations (Puchmayr et al. 2023, 2024; section 2.8) and a
  2-period QA stellarator with n* = 19 to 47 (Strumberger et al. 2025; section 3). Those are the cases
  closest to what PB3D computes. Access through IPP Garching only.

### 1.4 COBRA / COBRAVMEC (infinite-n ballooning)

- R. Sanchez, S.P. Hirshman, J.C. Whitson, A.S. Ware, "COBRA: An Optimized Code for Fast Analysis of Ideal
  Ballooning Stability of Three-Dimensional Magnetic Equilibria", *J. Comput. Phys.* **161**, 576–588
  (2000), doi:10.1006/jcph.2000.6514. [library]
- R. Sanchez, S.P. Hirshman, H.V. Wong, "Improved magnetic coordinate representation for ideal ballooning
  stability calculations with the COBRA code", *Comput. Phys. Commun.* **135**, 82–92 (2001),
  doi:10.1016/S0010-4655(00)00225-3 — this is COBRAVMEC, working directly in VMEC coordinates. [library]
- Code: folder `COBRAVMEC` in https://github.com/PrincetonUniversity/STELLOPT (MIT licence, active). [read]
  A Julia rewrite exists in the Wisconsin WISTELL group; a 2022 report mentions a bug passed to the
  STELLOPT maintainers that affected the Fortran-versus-Julia comparison. [snippet]
  https://wistell.engr.wisc.edu/wp-content/uploads/sites/741/2022/05/wistell_20220318.pdf
- What it solves: the ballooning ODE along a field line labelled by (s, α, θk); a matrix estimate with
  error proportional to h² (h the mesh step), refined by a variational principle to error proportional to
  h⁴, then Richardson extrapolation. Growth rate normalised to τ_A = (μ0 ρ R0²/B0²)^½ with B0 on axis (COBRA paper Fig. 4). [library]
- PB3D already has a "COBRA normalisation" (`norm_style = 2`, `Modules/eq_ops.f90` around line 5495):
  it uses a modified B_0 and ψ_0 and states EV_COBRA = EV_PB3D · β/2. [repo]
- What it gives the benchmark: the n → ∞ limit per field line, cheap, already used in both PB3D papers.
  It is the natural first 3-D check on W7-X (scan α and compare with PB3D's α-dependence).

### 1.5 Other infinite-n ballooning solvers for VMEC or DESC equilibria

- **DESC `BallooningStability` objective and "ideal ballooning lambda" quantity.** Tutorial:
  https://desc-docs.readthedocs.io/en/latest/notebooks/tutorials/ideal_ballooning_stability.html [read].
  Solves the infinite-n equation with Dirichlet ends; reports max over α of λ = γ² per surface; example
  uses 8 field lines on [0, π), 3 toroidal turns, 200 points per turn; also a cheaper "Newcomb ballooning
  metric" (sign only). Method paper: R. Gaur et al., "Omnigenous stellarators with improved ideal and
  kinetic ballooning stability", *Plasma Phys. Control. Fusion* **67**, 125015 (2025),
  doi:10.1088/1361-6587/ae1e9e, arXiv:2410.04576 [read: arXiv HTML]: finite differences in ζ giving a
  symmetric tridiagonal generalised eigenproblem; 14 field lines, 15 values of ζ0, ζ ∈ [−4π, 4π] in the
  convergence study. Data and scripts on Zenodo, doi:10.5281/zenodo.13887566. [read]
- **ballopt** (R. Gaur): https://github.com/rahulgaur104/ideal-ballooning-solver [read]. Python + SIMSOPT +
  VMEC2000, adjoint gradient of the ballooning eigenvalue. Paper: R. Gaur, S. Buller, M.E. Ruth,
  M. Landreman, I.G. Abel, W.D. Dorland, "An adjoint-based method for optimising MHD equilibria against the
  infinite-n, ideal ballooning mode", *J. Plasma Phys.* **89**, 905890518 (2023),
  doi:10.1017/S0022377823000995, arXiv:2302.07673. [read: arXiv abstract; journal data snippet] The paper
  states its numerical method is nearly identical to COBRAVMEC's. [snippet]
- **STESA** and its successor (Hudson & Hegna; Bird & Hegna; Cote et al.): COBRA-based ballooning solver;
  "STESA having been benchmarked directly against the COBRA solver" (Cote et al. 2019). S.R. Hudson,
  C.C. Hegna, *Phys. Plasmas* **10**, 4716 (2003) ("Marginal stability boundaries for infinite-n
  ballooning modes in a quasiaxisymmetric stellarator", title from search) [library: reference; snippet: title]
  https://ncsx.pppl.gov/StelPub/HH_03.pdf
- **3DBALLOON** (R. Sanchez, J.A. Jimenez, L. Garcia, A. Varias, *Nucl. Fusion* **37** (1997)):
  shooting code with fixed step and error proportional to h², which COBRA replaced. [library: COBRA reference list]
- No infinite-n ballooning solver inside SIMSOPT itself was found; ballopt sits on top of SIMSOPT. [snippet]

### 1.6 Newer global codes with 3-D equilibria

- **JOREK3D.** N. Nikulsin, R. Ramasamy, M. Hoelzl, F. Hindenlang, E. Strumberger, K. Lackner, S. Günter,
  *Phys. Plasmas* **29**, 063901 (2022), doi:10.1063/5.0087104, arXiv:2201.12033. [read] Reduced MHD,
  GVEC equilibria, benchmarked linearly against CASTOR3D (section 2.4). Follow-ups: R. Ramasamy et al.,
  "Nonlinear MHD modeling of soft β limits in W7-AS", *Nucl. Fusion* **64**, 086030 (2024),
  doi:10.1088/1741-4326/ad56a1 [read: reference list]; R. Ramasamy et al., "A comparison of low-n
  Mercier unstable Wendelstein stellarators and quasi-interchange modes in tokamaks", *J. Plasma Phys.*
  **91**, E119 (2025), doi:10.1017/S0022377825100512, arXiv:2505.03987 [read: abstract]. Access: the
  website https://www.jorek.eu/ offers only an e-mail contact; the wiki is restricted. [read]
- **M3D-C1, stellarator version.** Y. Zhou, N.M. Ferraro, S.C. Jardin, H.R. Strauss, "Approach to nonlinear
  magnetohydrodynamic simulations in stellarator geometry", *Nucl. Fusion* **61**, 086015 (2021),
  arXiv:2105.01186 [snippet]. Y. Zhou, K. Aleynikova, C. Liu, N.M. Ferraro, "Benign saturation of ideal
  ballooning instability in a high-performance stellarator", *Phys. Rev. Lett.* **133**, 135102 (2024),
  doi:10.1103/PhysRevLett.133.135102, arXiv:2404.01966 [read; details in section 3]. Source code:
  https://github.com/PrincetonUniversity/M3DC1, public, BSD-style PPPL licence. [read] No linear growth-rate
  benchmark against another code is reported in the 2024 paper; verification against MIPS, JOREK and NIMROD
  is listed as future work. [read]
- **NIMSTELL** (NIMROD for stellarators, Sovinec group, Wisconsin): tearing benchmark (Patil & Sovinec,
  APS DPP 2024) [read: reference list]; a 2026 Sherwood abstract reports that ideal ballooning growth rates
  rise with resolution and do not converge without damping. [snippet]
- **SPEC linear stability** (Hessian of multi-region relaxed MHD; A. Kumar et al. 2021 cylinder, 2022
  toroidal). The 2022 paper compares marginal stability with CAS3D and MISHKA-1 in simple tokamak and
  stellarator cases. [snippet] https://researchportalplus.anu.edu.au/en/publications/nature-of-ideal-mhd-instabilities-as-described-by-multi-region-re/
- **FAR3d** (Varela, Spong): gyro-Landau-fluid code with VMEC input, used for Alfvén and pressure-driven
  modes in LHD and other devices; no ballooning benchmark against an ideal code found. [snippet]
- **MISHKA or NOVA-K with 3-D equilibria:** none found. The 3-D CASTOR line is CASTOR_3DW → CASTOR3D.
- **Gyrokinetic codes** (GENE, GS2, stella, GX) compute kinetic ballooning modes on flux tubes; they are
  the physics that sets a lower β threshold than ideal MHD in W7-X (section 3), not a check on ideal codes.

### 1.7 Perturbative 3-D methods built on axisymmetric codes

- C.C. Hegna, "Effects of a weakly 3-D equilibrium on ideal magnetohydrodynamic instabilities",
  *Phys. Plasmas* **21**, 072502 (2014). Preprint: https://cptc.wisc.edu/wp-content/uploads/sites/327/2017/09/UW-CPTC_14-1.pdf [read]
- M.S. Anastopoulos-Tzanis, B.D. Dudson, C.J. Ham, C.C. Hegna, P.B. Snyder, H.R. Wilson, "Non-axisymmetric
  equilibrium and stability using the ELITE stability code", *Nucl. Fusion* **59**, 126028 (2019),
  doi:10.1088/1741-4326/ab40ef [read]; follow-up "Peeling-ballooning stability of tokamak plasmas with
  applied 3D magnetic fields", *Nucl. Fusion* **60**, 106003 (2020) [read: citation only].
- What these give the benchmark: an independent calculation of the same quantity PB3D computes for the
  ripple tokamak, built from axisymmetric ELITE eigenmodes and explicit sideband coupling (section 4).

---

## 2. Published 3-D benchmarks and intercomparisons

### 2.1 CAS3D versus TERPSICHORE

- **LHD test case, several codes.** Y. Nakamura, T. Matsumoto, M. Wakatani, S.A. Galkin, V.V. Drozdov,
  A.A. Martynov, Yu.Yu. Poshekhonov, K. Ichiguchi, L. Garcia, B.A. Carreras, C. Nührenberg, W.A. Cooper,
  J.L. Johnson, "Comparison of the calculations of the stability properties of a specific stellarator
  equilibrium with different MHD stability codes", *J. Comput. Phys.* **128**(1), 43–57 (1996),
  doi:10.1006/jcph.1996.0195. https://infoscience.epfl.ch/record/119493/ [read: abstract] An LHD
  configuration with an unusually flat pressure profile; the fastest growing n = 2 and n = 3 modes; 3-D
  codes and stellarator-averaged codes; the abstract reports good agreement, no numbers on the page.
  Note: the same equilibrium was the 10-period "testbed" in Cooper, Singleton & Dewar 1996 (section 5). [library]
- **Large-aspect-ratio configurations.** Sanchez et al. 2001 (CPC) state that TERPSICHORE "has also been
  benchmarked against the well-known CAS3D stability code for several large aspect ratio configurations",
  citing J. Nührenberg et al., *Plasma Phys. Control. Fusion* **35**, B115 (1993). [library: statement; the
  1993 title could not be confirmed]
- **NCSX / QAS, external kink.** M. Redi, E. Fredrickson, G. Fu, J. Canik, W.A. Cooper, C. Nührenberg,
  R.L. Dewar, APS-DPP 2000 talk, https://ncsx.pppl.gov/Scientificconf/2000_10_23_APS/redi.pdf [read]:
  N = 1 kink of QAS3_C82, CAS3D with wall at infinity stabilises at β ≈ 3.9 %, TERPSICHORE with wall at
  2.5 a is marginal at β = 4.03 %; the slides call the agreement "excellent". Fu 2001 (PVR slides) [read]:
  TERPSICHORE agrees with CAS3D on the critical β of C82, and with PEST and CAS3D on the n = 1 kink β limit
  of an ARIES-RS tokamak.
- A 1999 EPS paper says the CAS3D and TERPSICHORE packages "are found in agreement in numerous
  benchmarking studies". [snippet] https://lac913.epfl.ch/epsppd3/1999/web/pdf/p4085.pdf
- What it gives: the two global codes agree on marginal β to a few percent when walls are handled the same
  way; differences come from wall position (TERPSICHORE needs a wall) and Boozer-mapping resolution.

### 2.2 COBRA versus TERPSICHORE-VVBAL (ballooning limit)

- R. Sanchez, M.Yu. Isaev, S.P. Hirshman, W.A. Cooper, G.Y. Fu, J.A. Jimenez, L.P. Ku, M.I. Mikhailov,
  D.A. Monticello, A.H. Reiman, A.A. Subbotin, "Ideal MHD stability calculations for compact
  stellarators", *Comput. Phys. Commun.* **141**, 55–65 (2001), doi:10.1016/S0010-4655(01)00396-4. [library]
  - Cases: QAS C82 (3 periods; β = 3–5 %, ι from 0.26 to 0.47) with interchange-type modes at ι = 3/7 … 3/11;
    a QOS case (3 periods, A ≈ 3.5, ι 0.55–0.65, β = 0–3 %) with localised ballooning modes.
  - Agreement: both codes find the same unstable radial regions for every β and "similar critical β";
    eigenvalues are not comparable because the equations differ away from marginality (VVBAL has the
    (1 − λ) multiplier and no inertia term).
  - Causes of residual differences named by the authors: different normalisation; COBRA uses the VMEC
    solution while TERPSICHORE uses its own Boozer reconstruction, which differs most at rational surfaces;
    number of Boozer modes; box length (COBRA's K_w helical wells must match TERPSICHORE's S0 ≈ K_w ι/M
    poloidal transits; C82 needed K_w ≈ 65–70, i.e. S0 = 12).
- Fu 2001 PVR slides [read]: "Terpsichore-VVBAL's marginal stability points agree well with those of
  Cobra" for a quasi-helical and a QAS case; γ depends on (s, θk, α).

### 2.3 CASTOR3D versus CAS3D(N)

From the 2016 preprint of Strumberger & Günter 2017 (section 3.3 there) [read]:
- Axisymmetric circular equilibrium, internal kinks n = 1 … 7: CASTOR3D (NEMEC, Boozer, straight-field-line
  coordinates), original CASTOR and CAS3DN agree within 0.5 %.
- AUG-like shaped equilibrium, n = 3 and n = 5: coordinates agree within 0.25 %; CASTOR3D versus CAS3DN
  within 3.5 %, largest for n = 5; the authors attribute this to the narrower mode at higher n and to the
  different radial coordinate (CASTOR3D s = √S, CAS3DN S).
- **3-D helical-core tokamak** (NEMEC, axis shift ΔR): eigenvalues of the dominant-n = 3 mode agree within
  1 %; growth rates of dominant n = 2, 3, 4 modes agree within 1 % for ΔR ≤ 3.16 cm and within 4 % up to
  ΔR = 5.6 cm, where 95 (NEMEC coordinates) to 147 (Boozer) harmonics were needed and accuracy drops.
- What it gives: a published 3-D global-versus-global agreement level (1–4 %) to set expectations for
  what "agreement" means between global codes before the single-n question enters.

### 2.4 JOREK3D versus CASTOR3D (W7-A)

Nikulsin et al. 2022 [read: full text]: W7-A (5 periods), low β.
- Tearing (2,1): maximum deviation 13.2 % over a β scan and 26.2 % over a resistivity scan; attributed to
  the neglected parallel velocity in reduced MHD.
- **Ballooning, mode family n = 0, 5, 10:** β = 0.11 % and 0.21 %; maximum deviation 6.4 %, the other three
  cases below 3 %. Both codes truncate the family at the same number of toroidal harmonics, because
  ballooning growth rates diverge as n → ∞ without flow.

### 2.5 DESC versus COBRAVMEC

Gaur et al. 2025 (PPCF, Fig. 3a) [read: arXiv HTML]: radial profiles compared after scaling COBRAVMEC by
λ_max,DESC/λ_max,COBRAVMEC, because of normalisation differences; no agreement numbers. So the published
check covers the profile shape, not absolute growth rates.

### 2.6 Other code checks in the ballooning limit

- STESA versus COBRA (Hudson & Hegna 2003), as stated by Cote et al. 2019. [library]
- SPEC versus CAS3D and MISHKA-1, marginal stability. [snippet]

### 2.7 3-D tokamaks: RMP and ripple

- **CASTOR3D, AUG with magnetic perturbations.** J. Puchmayr, M. Dunne, E. Strumberger, M. Willensdorfer,
  H. Zohm, F. Hindenlang, ASDEX Upgrade Team, "Helical mode localization and mode locking of ideal MHD
  instabilities in magnetically perturbed tokamak plasmas", *Nucl. Fusion* **63**, 086008 (2023),
  doi:10.1088/1741-4326/acdd12; and "Effect of symmetry-breaking on the MHD edge stability limit of tokamak
  plasmas", *Nucl. Fusion* **64**, 086013 (2024), doi:10.1088/1741-4326/ad54d6 [read: citations]. Reported:
  symmetry breaking lowers the critical edge pressure gradient by up to 30 %, 3-D geometry raises the
  edge growth rates, consistent with experiment. [snippet] PhD thesis: J. Puchmayr, LMU Munich (2023),
  https://edoc.ub.uni-muenchen.de/33197/1/Puchmayr_Jonas.pdf (over 10 MB; not read).
  What it gives: a published global, family-coupled peeling-ballooning calculation on 3-D tokamak
  equilibria — the closest existing counterpart to PB3D's ripple/RMP runs.
- **Infinite-n ballooning on 3-D VMEC tokamak equilibria.** C.J. Ham, I.T. Chapman, A. Kirk, S. Saarelma,
  "Modelling of three dimensional equilibrium and stability of MAST plasmas with magnetic perturbations
  using VMEC and COBRA", *Phys. Plasmas* **21**, 102501 (2014), doi:10.1063/1.4895602, arXiv:1312.3104
  [read: abstract]. M. Willensdorfer et al., *Phys. Rev. Lett.* **119**, 085002 (2017),
  doi:10.1103/PhysRevLett.119.085002 [library]; T.B. Cote et al., "Helically localized ballooning
  instabilities in three-dimensional tokamak pedestals", *Nucl. Fusion* **59**, 016015 (2019),
  doi:10.1088/1741-4326/aaf01d [library]: 27 free-boundary VMEC equilibria of AUG with n = 2 MP,
  STESA-type solver, growth rate depends on α with the N = 2 periodicity of the perturbation.
  M. Willensdorfer et al., *Plasma Phys. Control. Fusion* **61**, 014019 (2019) compared measured boundary
  displacement with VMEC and MARS-F. [snippet]
- **ELITE perturbative**, Anastopoulos-Tzanis et al. 2019/2020 (sections 1.7, 4).
- **MARS-F** solves the linear response of an axisymmetric equilibrium; it is not a stability code for 3-D
  equilibria, and no 3-D stability benchmark with it was found. No M3D-C1 linear stability benchmark on
  3-D tokamak equilibria was found either.
- **Ripple:** the only 3-D peeling-ballooning ripple study found is PB3D's own (section 2.8). JET has
  about 0.08 % ripple, JT-60U about 1 %. [snippet]

### 2.8 PB3D's own published benchmarks

- T. Weyens, R. Sánchez, G. Huijsmans, A. Loarte, L. García, "PB3D: A new code for edge 3-D ideal linear
  peeling-ballooning stability", *J. Comput. Phys.* **330**, 997–1009 (2017), doi:10.1016/j.jcp.2016.10.054. [library]
  - Axisymmetric CBM18 (circular, ballooning unstable; HELENA and a VMEC port): PB3D follows MISHKA's
    n-dependence and marginal n but is about 10 % more unstable, the same kind of offset ELITE shows
    against MISHKA; COBRA marks the n → ∞ limit (Fig. 3).
  - 3-D test: the VMEC CBM18 with minor radius varied by 10 % over ζ = 0 … 2π. PB3D and COBRA both find the
    3-D case more unstable than either axisymmetric end case (Fig. 6). Agreement is qualitative.
  - α-dependence: for small parallel box θ_max the result depends on α; the difference falls linearly with
    box size, and the paper states that for an integration over −∞ … ∞ "the field line label α
    [should] have no influence on the final results" (Fig. 8).
  - Peeling verification was postponed until the vacuum term was complete.
- T. Weyens, PhD thesis "3D Ideal Linear Peeling Ballooning Theory in Magnetic Fusion Devices", UC3M/TU/e
  (2016–17), doi:10.13140/RG.2.2.26457.11361. [library] Chapter 5 (toroidal field ripple on H-mode
  equilibria, Table 5.3): COBRA results "similar but slightly more unstable"; some cases stable in PB3D but
  unstable in COBRA's n → ∞ limit; parallel range −10π … 10π.
- T. Weyens, A. Loarte, G. Huysmans, R. Sánchez, L. García, "Influence of toroidal field ripple on 3-D
  ballooning stability of H-mode plasma configurations", submitted to PPCF (2017); publication status not
  found. [library] EPS 2016 P4.067, EPS 2017 O4.128, EPS 2018 P1.1044. [snippet]
- Current repo (README, issue #4) [repo]: after the vacuum sign and 2π fix, the circular external kink
  agrees with MISHKA to 1 % (ε = 0.05) and 2 % (ε = 0.33), and "the 2017 MISHKA comparison of the thesis is
  resolved in MISHKA's favour". The published 3-D benchmarks therefore predate the vacuum correction.
- No comparison of PB3D 3-D results with CAS3D, TERPSICHORE or CASTOR3D exists.

---

## 3. Stellarator ideal-MHD reference cases

### 3.1 W7-X

- **Configuration labels.** Three-letter codes: first letter the mirror ratio in 1 % steps (A = 0 %,
  K ≈ 10 %), second letter the ι class (B low-ι with ι = 5/6 at the boundary, I/J standard with 5/5
  islands, T… high-ι), third the horizontal shift. EIM is the standard configuration; Aleynikova et al.
  use AIM (zero mirror), KIM (high mirror, 10 %), EBM and ETM (low and high ι(0)). K. Aleynikova, A. Zocco, J. Geiger, "Influence of magnetic configuration properties on
  kinetic ballooning modes in W7-X", *J. Plasma Phys.* **88**, 905880411 (2022),
  doi:10.1017/S0022377822000745, open access. [read]
- **Ideal-MHD β limits.**
  - EIM "is MHD-stable up to a volume-averaged β (≈5 %)" (Aleynikova et al. 2022). [read]
  - Zhou et al. 2024 PRL [read: arXiv HTML]: free-boundary VMEC equilibria of EIM at ⟨β⟩ = 4.9, 5.4 and
    6 %, B ≈ 2.3 T, zero toroidal current, near-parabolic pressure flattened in the core; COBRA finds the
    high-β equilibria locally ballooning unstable but Mercier stable. Earlier design numbers quoted there:
    ideal ballooning at β = 4.7 % and resistive interchange at 7.5 % (C. Beidler et al., *Fusion Technol.*
    **17**, 148 (1990)), β limit 4.3 % (G. Grieger et al., *Phys. Fluids B* **4**, 2081 (1992)). The
    alternative "TEH" configuration (higher mirror, ι ≈ 0.8, nearly zero shear) is Mercier unstable at β = 1 %.
  - Nührenberg 1996: CAS3D survey of the W7-X configuration class; "scenarios of safely stable operation".
    [snippet]
  - CAS3D global-mode β limit of about 1.7 % for low-shear configurations (J. Geiger, C. Nührenberg et al.,
    EPS 2023 poster P4.064). [snippet] https://pure.mpg.de/pubman/item/item_3587739_1/component/file_3587740/P4_064.pdf
  - Follow-up M3D-C1 study: Y. Zhou et al., *Phys. Plasmas* **33**, 032505 (2026), arXiv:2603.00869. [snippet]
- **Kinetic ballooning thresholds** (GENE), well below the ideal limits: local β_crit ≈ 0.82 % (AIM),
  0.98 % (EIM), 1.05 % (KIM), 0.85 % (EBM), 1.02 % (ETM) (Aleynikova et al. 2022). [read] Also K. Aleynikova,
  A. Zocco, P. Xanthopoulos, P. Helander, C. Nührenberg, "Kinetic ballooning modes in tokamaks and
  stellarators", *J. Plasma Phys.* **84**, 745840602 (2018), doi:10.1017/S0022377818001186 [read: abstract],
  which relates the KBM thresholds to CAS3D ideal results (numbers not on the open page).
- **Equilibrium files you can get today.**
  - SIMSOPT test files, https://github.com/hiddenSymmetries/simsopt/tree/master/tests/test_files [read]:
    `input.W7-X_without_coil_ripple_beta0p05_d23p4_tm` and
    `wout_W7-X_without_coil_ripple_beta0p05_d23p4_tm_reference.nc` (NFP = 5, PHIEDGE = 2.15, fixed boundary,
    MPOL = NTOR = 10, zero current, ⟨β⟩ = 5 % by the file name), plus `input.W7-X_standard_configuration`.
    These run directly in PB3D, COBRAVMEC and (after conversion to the text wout format of the test case)
    TERPSICHORE.
  - DESC examples, https://github.com/PlasmaControl/DESC/tree/master/desc/examples [read]: `W7-X`, and
    also NCSX, HSX, ATF, ARIES-CS, WISTELL-A, HELIOTRON, ESTELL, precise QA/QH.
  - VMEC++ ships a W7-X input in `examples/data`. [snippet] https://pypi.org/project/vmecpp/
  - IPP W7-X reference equilibria, https://www.ipp.mpg.de/4631651/Equilibria [read]: verified collection
    "available for internal and external researchers", served by the W7-X VMEC web service, which did not
    resolve from outside IPP in this session; no DOI or bulk download.
  - QUASR: about 370 000 vacuum QA/QH stellarators with coils, SIMSOPT and VMEC formats, Zenodo
    doi:10.5281/zenodo.13717741 (vacuum only, so a pressure profile must be added). [snippet]
  - ConStellaration (Proxima Fusion, arXiv:2506.19583): about 158 000 QI-like boundaries with VMEC++
    equilibria at five β levels, on Hugging Face. [snippet]

### 3.2 Other stellarator cases with published numbers

- **NCSX QAS3_LI383** (Redi et al. 2000 slides [read]): COBRA stable at β = 4.0 %, unstable at 4.25 %;
  TERPSICHORE-VVBAL marginal at 4 %. DESC has an NCSX example equilibrium.
- **QAS3_C82**: kink marginal β 3.9 % (CAS3D, wall at infinity) and 4.03 % (TERPSICHORE, wall at 2.5 a);
  high-n global ballooning (TERPSICHORE, n = 19, 22, 25) more stable than the infinite-n limit (section 5).
- **TERPSICHORE QAS3 fixed-boundary benchmark** in the public repository, with VMEC input, text wout,
  `fort.18` and expected output (n = 1 family, current-carrying). [read: file listing]
- **CASTOR3D 2-period QA** (Strumberger & Günter 2019; Strumberger et al. 2025 [read]): ⟨β⟩ = 2.9 %, edge
  q = 1.79 with flat edge q, R0 = 12.29 m, GVEC fixed boundary; unstable to external kinks and to
  peeling-ballooning modes whose growth rate rises with n* (n* = 19, 22, 44, 47 shown). This is a
  stellarator peeling-ballooning case from a global code; no data-availability statement was found.
- **LHD flat-pressure case** (Nakamura et al. 1996) and the 10-period β0 = 4 % testbed of Cooper et al. 1996.
- **W7-A** (Nikulsin et al. 2022): β = 0.11 % and 0.21 %, ballooning family n = 0, 5, 10.
- **GIGA** (Gauss Fusion, 4 periods): S.A. Lazerson et al., "The fixed boundary plasma equilibrium basis for
  a one Gigawatt electric stellarator power plant", arXiv:2607.09346 (2026) [library]: COBRAVMEC (outboard
  midplane flux tube at φ = 0 found most unstable by a scan) and TERPSICHORE for the n = 0, 1, 2 families,
  Boozer m = 0…32, n = −16…16, conformal wall at +30 %. No public equilibrium file stated.

---

## 4. Mode families and what a single-n result means

### 4.1 Definitions

- W.A. Cooper, D.B. Singleton, R.L. Dewar, "Spectrum of ballooning instabilities in a stellarator",
  *Phys. Plasmas* **3**, 275 (1996), doi:10.1063/1.871853 [library]: with L field periods and stellarator
  symmetry, a family with index 1 ≤ K ≤ L/2 contains the toroidal numbers {L·N_e ± K}; families decouple
  exactly, like n in axisymmetry. Example for L = 10: modes dominated by n = 2, 8 and 12 are one family.
  To match WKB theory they had to keep 30 harmonics of the dominant n and 64 of its sidebands; leaving out
  the sidebands spoils even the qualitative agreement.
- Lazerson et al. 2026 [library]: for N_fp = 4 the families are n = 0, 1, 2.
- Strumberger et al. 2025 [read]: for a 2-period QA, odd and even n form the two families.
- Coupling is not limited to stellarators: a 3-D wall or a helical core couples toroidal harmonics in a
  tokamak (CASTOR3D preprint [read]).

### 4.2 Perturbation theory in the 3-D amplitude δ

- **Hegna 2014** [read: preprint]: expanding in the axisymmetric eigenmodes, the eigenvalue is
  ω_n² = ω_n0² + V_nn + Σ_k |V_nk|²/(ω_n² − ω_k0²) to order δ². For a purely 3-D correction V_nn = 0, so
  the change is quadratic in δ and always destabilising for the most unstable mode. Example: an N = 3
  distortion couples n = 10 to n = 7 and 13; n = 10 is destabilised, its neighbours stabilised, the
  spectrum becomes more peaked.
- **Anastopoulos-Tzanis et al. 2019** [read]: same result with ELITE. If the δ-order force operator F(1)
  is non-axisymmetric, the eigenfunction changes but the eigenvalue does not; at order δ² there are two
  terms, the axisymmetric part F(2) (dropped by the authors) and the sideband sum; the mode becomes a
  triplet {n − N, n, n + N} and localises poloidally "at sufficiently large applied field".
- **EPS 2019** (Anastopoulos-Tzanis, Ham, Wilson, UKAEA-CCFE-CP(19)51) [read: abstract page]:
  perturbation theory is accurate for δB_N/B0 ~ 10⁻⁴; for a growth rate that rises with n, the lower
  sideband destabilises and the upper one stabilises.
  https://scientific-publications.ukaea.uk/papers/the-impact-of-symmetry-breaking-on-pedestal-stability-of-tokamak-plasmas/
- **Near rational surfaces** small 3-D fields change the infinite-n boundary strongly through
  Pfirsch–Schlüter-like currents (T.M. Bird, C.C. Hegna, *Nucl. Fusion* **53**, 013004 (2013),
  doi:10.1088/0029-5515/53/1/013004, arXiv:1208.4313 [library: abstract]); a δ expansion fails there.

### 4.3 High n in 3-D: ballooning WKB, localisation and chaos

- R.L. Dewar, A.H. Glasser, "Ballooning mode spectrum in general toroidal systems", *Phys. Fluids* **26**,
  3038 (1983), doi:10.1063/1.864028 [library]: in axisymmetry the discrete growth-rate spectrum is bounded
  above by an accumulation point (the n = ∞ value); without axisymmetry the rays are unbounded, semiclassical
  quantisation breaks down, and the accumulation point broadens into a continuum band.
- R.L. Dewar, P. Cuthbert, R. Ball, "Strong 'quantum' chaos in the global ballooning mode spectrum of
  three-dimensional plasmas", *Phys. Rev. Lett.* **86**, 2321 (2001), arXiv:physics/0102065. [snippet]
- M.H. Redi, J. Canik, R.L. Dewar, J.L. Johnson, S. Klasky, W.A. Cooper, W. Kernbichler, "Ballooning
  stability of the compact quasiaxially symmetric stellarator", PPPL-3612 (2001),
  https://bp-pub.pppl.gov/pub_report/2001/PPPL-3612-abs.html [read: abstract]: broken symmetry localises the
  ballooning eigenfunction; unstable isosurfaces in (s, α, θk) are topologically spherical and rays run away
  in k-space; the global spectrum "must be described using the language of quantum chaos theory".
- N. Nakajima, S.R. Hudson, C.C. Hegna, *Fusion Sci. Technol.* **51**(1), 79–91 (2007): local quasi-modes
  in the covering space cannot be superposed into global modes when the unstable domains are spherical. [snippet]
- In 3-D tokamaks, field-line (helical) localisation is measured and modelled: Willensdorfer 2017,
  Cote 2019 [library]; Puchmayr 2023 [read: title].

### 4.4 What this means for a single-n PB3D result [inference]

- PB3D integrates each field line over −∞ … ∞ and becomes independent of α in that limit (JCP 2017,
  Fig. 8). If that α-independent limit equals a toroidal average of the 3-D coefficients, which the JCP
  statement suggests but I have not checked in the code, then in Hegna's language PB3D keeps the averaged
  equilibrium terms and leaves out the sideband sum.
  Both missing and kept corrections are of order δ². So for small δ, PB3D and a family-coupled code should
  agree at order δ (both equal the axisymmetric value plus nothing linear in δ), and their difference
  should scale as δ². A ripple-amplitude scan on the circular tokamak would test this directly.
- The sideband sum has denominators ω_n0² − ω_{n±N}². For high-n ballooning modes neighbouring n have
  nearly equal growth rates, so these denominators shrink and the expansion fails once the coupling V
  exceeds the n-spacing of the eigenvalues. Then the global mode localises on the worst field lines (the
  field-line localisation above), and the coupled growth rate approaches the maximum over α of the local
  eigenvalue, which can differ from the α-independent value already at order δ. The conditions under
  which a single-n result is accurate are therefore: small δ (δB/B ~ 10⁻⁴ per the 2019 EPS paper), no
  strong resonant shielding currents at rational surfaces in the region of interest (Bird & Hegna), and moderate n where |V_{n,n±N}| is small compared
  with |ω_n0² − ω_{n±N,0}²|.
- For W7-X the 3-D shaping is not small, so no δ ordering exists; a comparison there should be stated
  per family (CAS3D, TERPSICHORE) and per field line (COBRA), and differences are expected.
- No paper was found that quantifies the sideband effect specifically for a single-n high-n code against
  a family-coupled code in a stellarator.

---

## 5. The infinite-n ballooning limit as a benchmark

- **What COBRA returns.** One eigenvalue per (s, α, θk) — per flux surface, field line and ballooning
  angle (COBRA 2000 [library]; Fu 2001 [read]). COBRAVMEC users usually scan a few lines; Lazerson et al.
  2026 found the outboard midplane line at φ = 0 most unstable. DESC reports the maximum over 8–14 lines
  and several ζ0. [read]
- **Upper bound.** In axisymmetry the n = ∞ value is the upper bound of the discrete spectrum (Dewar &
  Glasser 1983). In 3-D, Cooper, Singleton & Dewar 1996 [library] compared WKB ballooning (ray tracing plus
  semiclassical quantisation) with TERPSICHORE for n = 1 … 14 and radial numbers N = 0, 1, 2 on a 10-period
  equilibrium at β0 = 4 %: "excellent qualitative agreement, and good quantitative agreement"; the n = ∞
  continuum band (between −λ = 0.004176 and 0.004187) lies above all finite-n eigenvalues found.
  [inference] So the maximum over field lines of COBRA's growth rate is the number PB3D's growth rates
  should stay below and approach as n grows.
- **Finite-n margin in 3-D.** G.Y. Fu, IAEA FEC 2000 talk,
  https://ncsx.pppl.gov/Scientificconf/2000_10_02_IAEA/fu_talk.pdf [read]: TERPSICHORE global ballooning
  modes with n = 19, 22, 25 in QAS C82 are "significantly more stable than the infinite-n limit"; the
  slides plot eigenvalue versus β with the infinite-n β limit marked lower. Redi et al. (APS 2000) show the
  same plot. [read]
- **Tokamak 1/n behaviour.** R.L. Dewar, J. Manickam, R.C. Grimm, M.S. Chance, "n-dependence of ballooning
  instabilities", *Nucl. Fusion* **21**, 493 (1981), doi:10.1088/0029-5515/21/4/008 [library]: PEST-II
  (n < 20) and WKB (n > 5) agree for 5 < n < 20; the WKB formula reduces to the 1/n correction only at very
  high n and is much better at moderate n; critical β versus n oscillates at low n, so extrapolation from
  low n is unreliable; the lower bound on β_crit comes from n = ∞. P.B. Snyder, H.R. Wilson, *Plasma Phys.
  Control. Fusion* **45**, 1671 (2003), doi:10.1088/0741-3335/45/9/308: finite-n corrections are
  stabilising for narrow pedestals. [library: abstract]
- **PB3D versus COBRA so far.** JCP 2017 Figs. 3 and 6: PB3D growth rates rise with n towards the COBRA
  line; thesis Table 5.3: COBRA "often give[s] some overestimation of the instability". [library]
- **Normalisation traps.** COBRA: γ in units of τ_A = (μ0 ρ R0²/B0²)^½; PB3D's COBRA mode states
  EV_COBRA = EV_PB3D · β/2 and uses a modified B_0 and ψ_0 [repo]; VVBAL eigenvalues differ by definition;
  DESC versus COBRAVMEC needed a rescaling. Comparing marginal points (β or pressure gradient where γ = 0)
  avoids all of these.
- **Gap.** No paper was found that extrapolates finite-n 3-D growth rates in 1/n and compares the limit
  with COBRA. The nearest are Cooper et al. 1996 (WKB quantisation against TERPSICHORE), Fu 2000
  (qualitative), and PB3D's own figures.

---

## 6. Practical routes to a reference run

| Code | Route | What you need |
|---|---|---|
| TERPSICHORE | Clone https://gitlab.epfl.ch/spc/public/terpsichore (Apache 2.0); `make` in `srcs/`; run the QAS3 test | VMEC wout in the text format of the test case; a wall (put it on the boundary for fixed boundary); choose the family and Boozer spectrum. Contributions need the CLA in the repo. Contact: SPC, EPFL Lausanne (named in the README). |
| COBRAVMEC | STELLOPT on GitHub (MIT) | VMEC wout; scan α and θk; match normalisation |
| DESC ballooning | GitHub PlasmaControl/DESC; tutorial above | DESC equilibrium (DESC has a W7-X example); import of VMEC files not checked here |
| ballopt | GitHub rahulgaur104/ideal-ballooning-solver | SIMSOPT, VMEC2000 |
| CAS3D | Max Planck Institute for Plasma Physics, Stellarator Theory division (Greifswald) [snippet]; no public distribution or request procedure found | Collaboration with IPP; CAS3D is the code behind the published W7-X β limits and has a no-wall free-boundary option |
| CASTOR3D | IPP Garching; no public distribution found | Collaboration; this group has published 3-D tokamak and QA-stellarator peeling-ballooning results |
| JOREK3D | JOREK consortium, contact via https://www.jorek.eu/ | Not a linear eigenvalue code; growth rates from initial-value runs |
| M3D-C1 | Public GitHub (PrincetonUniversity/M3DC1) | Nonlinear initial-value; stellarator capability per Zhou et al. 2021 |
| GVEC | Zenodo doi:10.5281/zenodo.15026781 (Hindenlang et al. 2025) | Equilibrium only (used by CASTOR3D and JOREK3D) |
| W7-X equilibria | SIMSOPT test files (β = 5 %), DESC examples; IPP reference equilibria through IPP | — |

Note: R. Sanchez, first author of COBRA, is a co-author of both PB3D papers (UC3M Madrid at the time).
[library: author lists]

---

## Gaps I could not fill

- No published comparison of a single-n high-n code with a family-coupled global code in a stellarator,
  and no 1/n extrapolation of 3-D finite-n results to COBRA.
- CAS3D access terms: none published. CASTOR3D: same.
- Numbers behind several W7-X statements are snippet-only: Nührenberg 1996 content, the 1.7 % low-shear
  CAS3D limit, the CAS3D critical β values in Aleynikova et al. 2018, Puchmayr 2024's 30 %.
- Title of J. Nührenberg et al., PPCF 35 B115 (1993), cited as the TERPSICHORE–CAS3D benchmark.
- Whether "C. Nührenberg, Nucl. Fusion 56 (2016) 076010, free-boundary stability of W7-X divertor
  equilibria" exists: not found; not cited above.
- W7-X reference equilibria via the IPP web service could not be reached from outside IPP.
- OSTI was unreachable; several publisher pages blocked automated reading.
