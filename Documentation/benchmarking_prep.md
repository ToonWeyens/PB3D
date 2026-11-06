# PB3D benchmarking: preparation and priorities

Status on 31 October 2026, branch `dev/vacuum-fixes` (issue #4). This document decides what to compare PB3D with, what each comparison needs, and in which order to do them. It condenses three literature surveys written the same afternoon; they hold the full notes and links, and the decisions are here:

- `benchmarking/lit_A_axisymmetric.md`: axisymmetric peeling-ballooning benchmarks, the cbm18 family, ELITE, MISHKA, GATO, M3D-C1, JOREK;
- `benchmarking/lit_B_3d_codes.md`: 3-D stability codes (CAS3D, TERPSICHORE, CASTOR3D, COBRA, DESC), 3-D benchmarks, mode families;
- `benchmarking/lit_C_vacuum_bem.md`: 3-D vacuum and boundary-element literature and tools.

**Tags.** (this week): a result of this week's validation work on PB3D. [read]: a survey author read the source text; [read, fig] marks a value read off a figure by eye, accurate to about ±0.005 in γ/ω_A. [meta]: authors, journal, volume, pages and DOI checked against Crossref, content from the abstract or from a citing paper. [snippet]: search snippet, a fetch tool's summary or a secondary source only; treat it as a lead. [inference]: reasoning by a survey author or by this document, not stated in any source. (local): a file in this repository or on Toon's Mac. Effort estimates in days are [inference] unless a survey gave them.

---

## 1. Purpose and present state

**What PB3D computes.** PB3D solves the linear ideal-MHD eigenvalue problem for one toroidal mode number n in the high-n (peeling-ballooning) ordering. It reads HELENA (2-D) and VMEC (3-D) equilibria and prints ω²/ω_A² in the MISHKA normalisation: ω_A uses B_0, the (0,0) Fourier mode of |B| on the magnetic axis, R_0, the major radius of the axis, and a flat density. PB3D has no wall (the vacuum extends to infinity), does not couple different n, keeps only the perpendicular kinetic energy (no parallel inertia) and uses the incompressible ordering. On a free boundary a boundary-element method gives the vacuum response: style 2 (axisymmetric, toroidal harmonics) for HELENA, style 1 (3-D collocation on field lines, kernel 1/|x−x'|) for VMEC.

**Validation so far.** Every number in this table comes from this week's work (this week).

| Quantity | Reference | Agreement | Status |
|---|---|---|---|
| Fixed boundary, circular case: HELENA twin against VMEC twin | the same equilibrium through PB3D's two input paths | 2e-6 | done |
| Radial convergence, fixed boundary | PB3D grid sequence | order 2 on the circular case; order 1 on cbm18a (22 rational surfaces in the window), where the 200-point regression anchors sit 0.5 % from the converged value | done |
| Free boundary, circular ε = 0.33, q_a ≈ 1.04 external kink, n = 10, 15, 20, 200 radial points | MISHKA-1 (2017 source, flat density, conformal wall at 10 a, which equals no wall) | +0.005 %, +0.63 %, +0.89 % | done |
| Same at ε = 0.05 | MISHKA-1 | −1.2 %, −0.55 %, −0.26 %; the gap falls as 1/n² | done |
| cbm18a with the plasma taken to ψ_N = 1 (pseudo-vacuum, flat density), n = 10, 15, 20 | MISHKA-1, adiabatic index 5/3 | PB3D more unstable by 6.6 %, 15 %, 18 % in ω²; MISHKA with γ = 0 identical to all digits (MISHKA-1 is incompressible) | open: cause not found; the pseudo-vacuum's rational surfaces are suspected (section 2a) |
| cbm18a free boundary, plasma to ψ_N = 0.75, true vacuum beyond, n = 10 | published ELITE, GATO, M3D-C1 curves | PB3D −4.5039e−2 at 200 points, −4.526e−2 converged, γ/ω_A = 0.212; ω_A of the curves unknown | not compared |
| 3-D vacuum: Green identities for 1, x, y, z, x²−y² on the real W7-X surface | exact | 1e-3 | done |
| 3-D vacuum: exterior point and ring sources | exact potentials | 0.6 to 1.6 % | done |
| 3-D vacuum: W7-X response matrix on six grids, 64 to 256 field lines × 121 points per line | grid sequence | 1.2 % | done; the trapezoidal sum along the open lines is still order 1 (end corrections) |
| Circular VMEC twin, free boundary through the 3-D vacuum, 64 lines | the same twin through the axisymmetric vacuum | 1.2 % below; radial limits −0.43901 and −0.43366 (radial order 2) | done |
| 3-D vacuum against another code | none yet | — | open (section 4) |
| 3-D stability: W7-X fixed boundary; circular tokamak with an n = 5 helical ripple of amplitude 0.01 to 0.04 | none | never compared | open (section 3) |
| Rotation invariance; scaling of R_0 and B_0 | exact | round-off | done |

**Two corrections to the surveys.** Survey A (section 5) quotes PB3D at 2.1 to 2.3 % from MISHKA, and survey B (section 2.8) quotes the README at 1 % and 2 %. Both predate today's change. Before the free edge got its trapezoidal quadrature weight, the circular eigenvalue converged only at order 1 in the radial grid and sat 2.3 % (this week) from MISHKA at 200 points; the 2.3 % had been misread as the toroidal correction that the high-n ordering drops. Survey A (section 6) also looks for the cause of MISHKA's NaN on the extended cbm18a in the density options. The causes were two: a square root of 2P/ρ with P a few 1e−20 below zero in the extension (clamped), and MISHKA's FFT, which needs a mapping with 2^k + 1 poloidal points; a 513-point (this week) HELENA mapping fixed that.

**What is not validated.**

- The pedestal case. The one code comparison, MISHKA on cbm18a with a pseudo-vacuum, has PB3D 6.6 to 18 % more unstable.
- The free-boundary cbm18a value. No other code has computed it at the same plasma edge and normalisation.
- The 3-D vacuum operator. Every check is internal: identities, exact sources, grid sequences, the axisymmetric twin.
- 3-D stability. No comparison with any other code. VMEC's Mercier criterion does not help: DMerc says the rippled twin is Mercier-unstable everywhere and W7-X Mercier-stable outside s = 0.07 (this week), while PB3D at n = 20 and 40 (this week) finds no instability on the rippled twin. That is no contradiction at finite n, because Mercier modes are localised and grow very slowly, but it means the Mercier test cannot check PB3D.

**Toon's tools.** MISHKA-1 and HELENA (hel12a, PB3D-patched) build and run on his Mac. VMEC++ runs in a Python venv. PB3D writes a VMEC input from a HELENA equilibrium (`--export_HEL`). Snellius (SURF) runs PB3D under Slurm. The code is on GitHub (ToonWeyens/PB3D) with CI on every push. The plan of issue #4 asks for local checks before requests to other groups (rule 8 of `plan_vacuum_bugfix.md`, local); section 5 keeps that order.

---

## 2. Axisymmetric benchmarks

### 2a. MISHKA-1: done, and what remains

Setup (this week): MISHKA-1 rebuilt from the 2017 source, flat density, conformal wall at 10 a; PB3D style 2 at 200 radial points. Values are −ω²/ω_A² in the MISHKA normalisation.

| n | ε = 0.33: PB3D | MISHKA | PB3D − MISHKA | ε = 0.05: PB3D | MISHKA | PB3D − MISHKA |
|---|---|---|---|---|---|---|
| 10 | 0.43380 | 0.43378 | +0.005 % | 0.42168 | 0.42679 | −1.2 % |
| 15 | 0.44562 | 0.44283 | +0.63 % | 0.44599 | 0.44844 | −0.55 % |
| 20 | 0.32548 | 0.32262 | +0.89 % | 0.33704 | 0.33792 | −0.26 % |

cbm18a (this week): the HELENA equilibrium extended with a zero-pressure, zero-current region from ψ_N = 0.70 to 1, q_a = 25, with the plasma taken to ψ_N = 1 in both codes. MISHKA's result does not depend on the wall radius beyond ψ_N = 1 and barely on the density in the extension. The γ column is the square root of the ω² ratio [inference].

| n | MISHKA, adiabatic index 5/3 | PB3D | PB3D/MISHKA − 1 in ω² | same in γ | MISHKA with γ = 0 |
|---|---|---|---|---|---|
| 10 | 0.04213 | 0.0449 | +6.6 % | +3.2 % | 0.04213 |
| 15 | 0.07284 | 0.0835 | +15 % | +7.1 % | 0.07284 |
| 20 | 0.09403 | 0.1112 | +18 % | +8.7 % | 0.09403 |

With γ = 0 MISHKA gives the same values to all printed digits (this week): MISHKA-1 is incompressible, as the literature says, and the adiabatic index in its namelist is not used by the equations. PB3D stays 6.6, 15 and 18 % more unstable.

What remains:

1. **Not compressibility, not resolution.** The γ = 0 run changed nothing (this week), so MISHKA-1 is incompressible (Pamela 2020; Snyder 2002 [read]) and the gap is not the γp|∇·ξ|² term. MISHKA is converged in its radial grid (NG 201 → 401: 4e−4 at n = 10, 1e−5 at n = 15; at n = 20 NG 201 did not converge and 401 did) and in its harmonics (41 instead of 21: +0.45 % at n = 10, +0.16 % at n = 20, with the in-core workspace enlarged) (this week); PB3D in its modes (20 → 30 → 40: 0.16 %) and radial points (this week). What is left [inference]: PB3D's pseudo-vacuum run equals its true-vacuum run (−0.11125 against −0.11124 at n = 20), so PB3D does not see the rational surfaces of the extension (about 230 of them for n = 10, q from 1.6 to 25), while an ideal pseudo-plasma cannot carry a radial field at a rational surface, which makes a pseudo-vacuum more stable than a vacuum, and MISHKA's grid (accumulated at ψ_N = 0.70) resolves them. The gap should then grow with n (more and closer resonances), as observed, and be absent on the circular cases, which have no extension, as observed. The test is item 4 below: a comparison with no pseudo-vacuum on either side. On cbm18 the literature's compressibility effect (Γ = 5/3 lowers γ by 13 to 18 % at n = 10 and about 6 % at n = 40 against Γ = 0 [read, fig]) is a separate matter for the comparison with ELITE and M3D-C1, which do have it.
2. **MISHKA's radial grid and harmonics.** Done (item 1); the table's values are the NG = 401, 21-harmonic ones.
3. **The ε = 0.33 gap grows with n** (0.005 % to 0.89 %), while the ε = 0.05 gap falls as 1/n². Running both codes at 400 radial points at n = 15 and 20 shows whether the growth is radial resolution [inference].
4. **A free-boundary MISHKA run on cbm18a.** A HELENA run whose boundary is the pedestal-foot surface, with MISHKA's own vacuum and a large `RWALL`, is the direct counterpart of PB3D's free-boundary run (survey A, section 6) [inference]. It would check −4.526e−2 without the pseudo-vacuum. Cut at ψ_N = 0.75 to match the regression anchor, or at 0.701 to match the literature (section 2b).
5. **History.** The 2017 JCP paper found PB3D about 10 % more unstable than MISHKA on CBM18, "the same kind of offset ELITE shows against MISHKA" [read]; that comparison predates the vacuum correction [read]. ELITE was benchmarked against MISHKA for 5 < n < 50 without a published percentage (Snyder 2002) [read].

### 2b. The published cbm18 curves

**The equilibrium.**

| Item | cbm18_dens8 in the literature | PB3D's cbm18a (local) |
|---|---|---|
| Origin | made by P. B. Snyder (GA) with TOQ; no journal paper defines it [read] | HELENA mapping `tests/fixtures/cbm18a_extended_vac.12.xz` |
| Shape | shifted circle, no X-point [read] | same |
| Plasma edge | ψ_N = 0.701 of the TOQ boundary flux, the "standard cutoff" [read]; Ferraro writes Ψ_w = 10/7 in pedestal-foot units [read] | pressure gradient ends at ψ_N = 0.70; zero pressure and current out to ψ_N = 1, q = 25 there |
| q at the edge | q_a = 3.03 [read] | 3.02 at ψ_N = 0.70; 1.05 on the axis |
| β_N; pedestal width | 1.51; L_ped/a = 0.049 [read] | — |
| Outer TOQ boundary | circle of radius 2 m about R = 3 m [read] | HELENA boundary at ψ_N = 1; R_geo/R_mag = 0.8319, B_geo/B_mag = 1.1635 |
| "dens8" | orders the series by pedestal height: dens6 has lower pedestal pressure and edge current [read]; scaling of density at fixed temperature is suggested by the name, not verified [inference] | which member cbm18a is has not been established; q_a and the cutoff point to dens8 [inference] |
| Density | not part of the equilibrium; flat in all ideal benchmarks [read]; M3D-C1 uses 100 × lower outside the plasma [read] | flat |
| Γ | ELITE and GATO with Γ = 0 and 5/3; M3D-C1 with 0, 1 and 5/3 [read] | incompressible ordering |
| First unstable n | n = 4 [read] | — |
| Sensitivity to the cutoff | ELITE scan over ψ_N = 0.675 to 0.727: curves differ visibly only below n ≈ 10 [read, fig] | PB3D's anchor uses 0.75, outside that scan |

Survey A found no public EQDSK or TOQ file of cbm18; the only public files are BOUT++ grid files on GitHub (BOUT-dev, `examples/elm-pb`) [snippet], whose domain covers 0.4 < ψ_N < 1.2 in pedestal-foot units (Dudson 2011) [read], so the core is probably missing [inference]. Survey A's geometry estimate puts the magnetic axis at R_mag ≈ 3.61 m [inference]; with R_geo = 3 m the ratio above gives R_mag = 3.606 m, and with B_geo = 2 T, if HELENA's B_geo is the vacuum field at R = 3 m, it gives B_mag = 1.719 T [inference]. The ratios come from the local golden values in `tests/fullstack/test_read_HEL.f90`; `Modules/HELENA_vars.f90` defines `RMtoG_H` = R_geo/R_mag and `BMtoG_H` = B_geo/B_mag (local).

**Which ω_A each curve uses.**

| Source | ω_A as stated | R/B | Density |
|---|---|---|---|
| Ferraro 2010 | B0/(R0√(4πρ_p)), B0 and R0 at the magnetic axis, ρ_p at the pedestal top [read] | not given | uniform inside [read] |
| Squires 2008, Dudson 2011 (ELITE) | label "R/B = 1.5" [read] | probably R = 3 m, B = 2 T [inference] | flat [inference] |
| Dudson 2011 (BOUT++) | V_A/R; τ_A = 0.37 μs at 10¹⁹ m⁻³ [read] | 1.80 m/T if deuterium [inference] | flat [read] |
| ELITE default (Lunniss 2016) | B/(R√(μ0ρ_a)), R the mean of the innermost and outermost R of the boundary, B = f/R there [read] | R ≈ 3.49 m for cbm18 [inference] | matters only if ρ varies [read] |
| MISHKA and PB3D | B and R on the axis | R_mag ≈ 3.61 m [inference] | flat |
| Pamela 2020 (JOREK) | physical units, s⁻¹, at ρ = 10¹⁹ m⁻³ [read] | needs R, B and ion mass | flat |

Ferraro's ELITE Γ = 0 curve and the Squires "R/B = 1.5" points agree within about 0.005 at n = 10 and 20 [read, fig], so the two figures share one normalisation. The papers do not say which [inference].

**Conversion.** For the same mass density,

    (γ/ω_A)_published = (γ/ω_A)_PB3D × (R_pub/B_pub) / (R_mag/B_mag),
    R_mag/B_mag = (R_geo/B_geo) × BMtoG/RMtoG = 1.3986 × R_geo/B_geo   (cbm18a, local values).

**Worked example.** PB3D's converged free-boundary value at n = 10 is −ω²/ω_A² = 0.04526, so γ/ω_A = 0.2127, the 0.212 of section 1. ELITE Γ = 0 at n = 10 is 0.160 (Squires) to 0.165 (Ferraro) [read, fig]. All rows below are [inference].

| Assumed ω_A of the published curve | Factor (R_pub/B_pub)/(R_mag/B_mag) | PB3D converted | PB3D against ELITE |
|---|---|---|---|
| Axis values (Ferraro's stated formula) | 1 | 0.213 | +29 to +33 % |
| R_geo = 3 m and B_geo = 2 T (the "R/B = 1.5" label) | RMtoG/BMtoG = 0.7150 | 0.152 | −5 to −8 % |
| ELITE default, R ≈ 3.49 m, B = f/R with f = R_geo·B_geo | (3.49/3)² × 0.7150 = 0.968 | 0.206 | +25 to +29 % |

The spread, −8 % to +33 %, comes from the choice of R and B alone, so for now it says nothing about PB3D. Two further differences remain in this example: PB3D cut the plasma at ψ_N = 0.75, not 0.701, and its radial window starts at ψ_N = 0.45 (local, `tests/regression/cbm18a_free_bnd.sh`). Two comparisons need no ω_A at all [inference]: the shape γ(n)/γ(n = 20), and the first unstable n (n = 4 [read]).

**What PB3D must run.** A free boundary with the plasma–vacuum interface at ψ_N = 0.701, no wall, flat density, the Γ = 0 curves as reference (PB3D's ordering is incompressible), n = 5, 7, 10, 15, 20, 30, 40 (the n with markers), each n converged radially (the 200-point anchor is 0.5 % off at n = 10). Keep the 0.75 run for continuity. ELITE carries two orders in 1/n and is accurate for n ≳ 5; arbitrary-n ELITE and the original agree for n ≳ 7 [read]. Compare at n ≥ 10 first.

**Figures to digitise.** All three curve figures are vector graphics (no raster image on those pages) [read], so the marker coordinates can be read from the PDF drawing commands, for instance with PyMuPDF's drawing list, and mapped through the axis ticks [inference]. Downloading the PDFs needs Toon's go-ahead. Contents as listed in survey A [read]:

| Source | Where | Codes and n | PDF |
|---|---|---|---|
| Ferraro, Jardin, Snyder 2010 | Fig. 6 | ELITE Γ = 0 and 5/3 (n ≈ 5 to 40), GATO Γ = 0 and 5/3 (to n = 10), M3D-C1 Γ = 0, 1, 5/3 (n = 4 to 40) | https://m3dc1.pppl.gov/Papers/ferraro_PP_10.pdf |
| Squires et al. 2008 | slide 9 | ELITE "R/B = 1.5", GATO, NIMROD, n = 1 to 20; ELITE cutoff scan | https://w3.pppl.gov/cemm/APS2008/Squires.pdf |
| Dudson et al. 2011 | Fig. 2 | ELITE "R/B = 1.5" to n = 60, GATO to n = 10, BOUT++ | arXiv:1008.4554 |
| Dudson et al. 2009 | text | printed: ELITE 0.239, BOUT++ 0.245 at n = 20 | arXiv:0810.5757 |
| Squires et al. 2008 | text | printed: ELITE and NIMROD ≈ 0.132 at n = 8 | as above |

Survey A's readings by eye, for orientation only [read, fig]:

| n | ELITE "R/B = 1.5" | ELITE Γ = 0 (Ferraro) | GATO Γ = 0 | M3D-C1 Γ = 0 | ELITE Γ = 5/3 | M3D-C1 Γ = 5/3 |
|---|---|---|---|---|---|---|
| 5 | 0.059 | – | 0.06 | 0.055 | 0.049 | ≈0.02 |
| 10 | 0.160 | 0.165 | 0.157 | 0.160 | 0.135 | 0.139 |
| 15 | 0.208 | 0.207 | – | 0.207 | ≈0.183 | 0.189 |
| 20 | 0.237 | 0.239 | – | 0.236 | 0.215 | 0.221 |
| 30 | – | 0.272 | – | 0.269 | 0.250 | ≈0.253 |
| 40 | – | 0.290 | – | 0.286 | 0.271 | 0.270 |

**Agreement the literature reaches between codes.**

| Pair | Case | Agreement | Source |
|---|---|---|---|
| ELITE and GATO | cbm18_dens8, n = 5 to 10 | within ≈0.005 in γ/ω_A | Squires 2008; Ferraro 2010 [read, fig] |
| ELITE, GATO, NIMROD | cbm18_dens8, n ≥ 5 | within ≈0.01 | Squires 2008 [read, fig] |
| BOUT++ and ELITE | shifted circle, n = 20 | 2.5 % | Dudson 2009 [read] |
| NIMROD, mapping only | CBM18 | ≈10 % from small mapping differences | Burke 2010, via Ferraro 2010 [read] |
| MARG2D and ELITE | Meudas-1, cutoff 99.6 and 99.8 % | within 5 %; ILSA ≈10 % lower | Ferraro 2010 [read] |
| NIMROD and ELITE | Meudas-1 | ≈5 %, 27 % at n = 4 | King 2016 [snippet] |
| JOREK and BOUT++ | circular, q95 = 3.0 | 4 % | Hoelzl 2021 [read] |

Target for PB3D: a few percent at n ≥ 10, once ω_A, density, Γ and the cutoff match [inference].

### 2c. The ELITE request

ELITE is not public; its developers are at General Atomics and the University of York [inference]. ELITE reads HELENA output (Lunniss 2016 ran it on HELENA files) [read], so it can run on PB3D's own equilibrium.

| Item | What to send or ask |
|---|---|
| Equilibrium | the cbm18a HELENA output used by PB3D (`tests/fixtures/cbm18a_extended_vac.12.xz`, or the 513-point mapping made for MISHKA); ask which HELENA file ELITE reads |
| Plasma–vacuum interface | ψ_N = 0.701 (standard cutoff); a second run at ψ_N = 0.75 (PB3D's anchor). If ELITE cannot cut a HELENA file at an interior surface, Toon makes a HELENA file whose boundary is the 0.701 surface |
| n | 5, 10, 15, 20, 30, 40 |
| Wall, density | no wall (infinite vacuum); flat density |
| Γ | 0 and 5/3 (ELITE's Γ = 5/3 is accurate only to lower order in 1/n [read]) |
| Ask also | the R, B and ρ that define ω_A in the published "R/B = 1.5" curves and in Ferraro's Fig. 6; the numbers behind those curves; the TOQ files of cbm18_dens8; the unpublished "A Practical User's Guide for ELITE" (Snyder and Wilson 2006, cited in Lunniss 2016 [read]) |
| Send in return | PB3D's values for the same runs, with the normalisation spelled out |

Text Toon can paste:

> I maintain PB3D, a high-n ideal-MHD code for peeling-ballooning modes in 2-D and 3-D equilibria (Weyens et al., J. Comput. Phys. 330 (2017) 997). This autumn I corrected its vacuum term; it now agrees with MISHKA-1 to within 1.2 % on circular external kinks at n = 10 to 20, and I would like to compare it with ELITE on the cbm18_dens8 case. Could you tell me which R, B and mass density define ω_A in the ELITE curves labelled "R/B = 1.5" (Squires et al., APS 2008; Dudson et al., PPCF 2011) and in Fig. 6 of Ferraro et al. (Phys. Plasmas 2010), and send the numbers behind those curves if you still have them? I would also be grateful for an ELITE run on the attached HELENA file of the same equilibrium, with the plasma-vacuum interface at ψ_N = 0.701 and, separately, at 0.75, no wall, flat density, Γ = 0 and Γ = 5/3, for n = 5, 10, 15, 20, 30 and 40, giving γ/ω_A with ω_A stated. If the TOQ files of cbm18_dens8 and the ELITE user's guide that describes the cbm and dbm cases can be shared, they would settle the remaining questions about the equilibrium. I will send PB3D's results for the same cases in return.

Timing is Toon's decision. Rule 8 of the plan puts requests after the local checks; the question about ω_A concerns published data only and blocks the absolute comparison of section 2b, so it could go earlier.

### 2d. Pamela 2020: the Zenodo CBM18 dataset

Pamela et al. 2020 (Fig. 4) compare JOREK reduced and full MHD with MISHKA (incompressible ideal) and CASTOR (ideal and resistive, η0 = 6 × 10⁻⁶ Ω m) on CBM18 for n = 1 to 20, γ in 10⁶ s⁻¹ for ρ = 10¹⁹ m⁻³, log axis; first unstable n = 4; "reasonable" agreement [read]. The figure is a raster image, but the numbers are in the dataset doi:10.5281/zenodo.3971940 (CC BY 4.0): `CBM18.zip`, 136.6 kB, with `growth_rates_castor.txt` (742 B), `growth_rates_cbm18_new.txt` (746 B) and `plot_growth_rates.py` [snippet]. Nobody has downloaded it; that needs Toon's approval.

How to use it:

1. Unpack into its own empty directory. Read `plot_growth_rates.py` for the columns and units rather than running it.
2. Convert: γ/ω_A = γ R √(μ0ρ)/B with ρ = 10¹⁹ m⁻³ × ion mass. R, B and the ion mass must come from the paper or the script; the surveys do not have them. For R/B = 1.5 m/T and deuterium, 1/ω_A = 0.31 μs, and the n = 20 value of 0.75 to 0.8 × 10⁶ s⁻¹ becomes 0.23 to 0.25 [read, fig; inference].
3. Compare its MISHKA curve with Toon's MISHKA on cbm18a at n = 10, 15, 20. The paper does not say where the plasma boundary was [read]. The same n-dependence would confirm Toon's setup; a different one points to the boundary or the density.
4. Then compare PB3D with the MISHKA and ideal CASTOR columns. JOREK reduced MHD is a different model; use it only as a bracket.

---

## 3. 3-D fixed-boundary benchmarks

Agreement between global codes in 3-D, to set expectations before the single-n question enters:

| Pair | Case | Agreement | Source |
|---|---|---|---|
| CASTOR3D, CASTOR, CAS3DN | circular, internal kinks n = 1 to 7 | 0.5 % | Strumberger 2017 [read] |
| CASTOR3D and CAS3DN | AUG-like, n = 3 and 5 | 3.5 % | same [read] |
| CASTOR3D and CAS3DN | 3-D helical core, dominant n = 2, 3, 4 | 1 % for ΔR ≤ 3.16 cm, 4 % up to 5.6 cm | same [read] |
| JOREK3D and CASTOR3D | W7-A ballooning, family n = 0, 5, 10 | 6.4 % at most, three cases below 3 % | Nikulsin 2022 [read] |
| CAS3D and TERPSICHORE | QAS3_C82 n = 1 kink, marginal β | 3.9 % (no wall) against 4.03 % (wall at 2.5 a) | Redi 2000 [read] |
| COBRA and TERPSICHORE-VVBAL | QAS C82, QOS | same unstable regions, "similar critical β"; eigenvalues not comparable | Sanchez 2001 [read] |
| 3-D and stellarator-averaged codes | LHD flat-pressure case, n = 2 and 3 | "good agreement"; no numbers on the abstract page | Nakamura 1996 [read, abstract] |

### 3a. TERPSICHORE (open source since 2025)

**Status, night of 31 October 2026 (this week):** built on Snellius without changes, benchmark reproduced; `scripts/terpsichore/` holds the netCDF-to-fort.18 converter and the card generator; the shipped 1.2 release has a parity error in its MODELK 2 kinetic energy (`tprgl5_block_ap.f` 776 and 779, sine and cosine tables swapped for the η–η and ξ–η terms; the poloidal displacement then carries no mass), which made the circular twin come out at λ = +0.151 fixed and −1.725 with the wall at 2a; with the two lines swapped TERPSICHORE gives +0.07851 and −0.43847 against PB3D's +0.07843 and −0.43914 (0.10 % and 0.15 %), MISHKA 0.43378 (1.1 %). The cbm18a twin shows core resonant modes at q = 1.1 that PB3D's window excludes; with the corrected kinetic energy it finds no eigenvalue below −0.045 (PB3D −0.0446, MISHKA −0.0421) and a spurious −0.011 inside the zero-pressure extension, so the pedestal comparison needs a cbm18a without the extension, as for MISHKA. Steps 1 to 3 below remain as written.

**Later that night, the first 3-D case:** the rippled circular twin (δ = 0.04, five-fold ripple) through `wout2terps.py --nfp 5`, n = 10 family n = 0..20. The external kink with the wall at 2a: TERPSICHORE −0.41903 against PB3D's −0.42163 at 64 field lines and −0.42532 at 128 (0.6 % and 1.5 %, PB3D more unstable, growing with the line count, so not the line discretisation); the single n = 10 row gives −0.41246, so the sideband coupling that PB3D carries along field lines is 1.6 % here. The fixed-boundary family spectrum is dense and gives nothing to compare. W7-X at n = 10 is dense and stable (ι = 10/11 inside the plasma); a W7-X comparison needs an unstable equilibrium. Item C7 (CASTOR3D on the rippled case) now has a TERPSICHORE number to stand next to.

**1 November 2026, the convergence study:** PB3D's 3-D vacuum at 128 field lines × 400 parallel points gives −0.43849 for the circular kink against TERPSICHORE's −0.43848 (converged over 159 to 639 surfaces and 16 to 32 vacuum intervals), and its response is within 0.4 % of the axisymmetric style on a 1025-point HELENA mapping. The axisymmetric style is first order in HELENA's poloidal points (−0.43380 at the fixture's 257, −0.43695 at 1025), so the earlier 1e−4 agreement with MISHKA-1 at 257 points was a coincidence; MISHKA-1 is converged at 0.43377 and 1.1 % less unstable than the other two codes, probably a model difference in the inertia. The rippled twin follows the same pattern against TERPSICHORE's family value. Rule for the 3-D vacuum: at least 128 lines and about 36 parallel points per poloidal wavelength.

**1 November 2026, morning, the pedestal case settled:** the cbm18a plasma cut at the pedestal foot (`tests/reference/cbm18a_cut/`) gives, at n = 10, 15, 20, PB3D −0.04012, −0.08128, −0.11120 fixed and −0.05002, −0.08577, −0.11209 free against TERPSICHORE −0.03991, −0.07951, −0.10701 and −0.04979, −0.08441, −0.10753 (0.5 % to 4 %), while MISHKA-1 gives 0.03637, 0.07009, 0.09290 and 0.04492, 0.07407, 0.09453, 9 to 16 % less unstable than both. Items 1 and 3a are therefore done at the few-percent level; what remains open is MISHKA-1's offset, and ELITE (item 2) would be the natural third opinion on the pedestal mode.

TERPSICHORE solves global linear ideal MHD with hybrid finite elements radially and Fourier in both angles, maps the VMEC solution to Boozer coordinates and runs one mode family per call [read]. It is the only global 3-D ideal code that can be downloaded and run without an agreement [inference, survey B]. It treats the vacuum as a pseudo-plasma inside a wall [read], so only fixed-boundary results compare with PB3D. Fu (IAEA 2000) ran it at n = 19, 22 and 25 [read], so n = 20 is within its range.

Install plan:

1. Clone https://gitlab.epfl.ch/spc/public/terpsichore anonymously (Apache 2.0; created May 2025; tag 1.1 of 1 June 2025, "the first open-source TERPSICHORE version") [read]. The clone needs Toon's go-ahead.
2. `make` in `srcs/` [read].
3. Run `QAS3_fx_benchmark/` (`input.QAS3_fx_cur75_bench`, a text wout file, expected output) and the regression test `tests/QAS3_fx_cur_bench_n1/`, and compare with the expected output [read]. This is an n = 1, current-carrying case: it checks the installation, not PB3D, whose ordering needs high n [inference].
4. Convert PB3D's netCDF wout files to the text wout format of the test case [read]. How to produce that format from VMEC++ is not in the surveys; the test case's wout and TERPSICHORE's reader show it [inference].
5. Choose the wall on the plasma boundary (fixed boundary), the Boozer spectrum (Lazerson 2026 used m = 0 to 32, n = −16 to 16 for a four-period stellarator [read]) and the family. Whether TERPSICHORE can set Γ = 0 is not in the surveys; check the input, since PB3D is incompressible [inference].

What to compare, in this order:

1. **Axisymmetric twins** (circular `wout_circ_eps033_q104.nc`, local; a cbm18a VMEC twin from `--export_HEL`): n = 10, 15, 20. No coupling exists here, so TERPSICHORE, PB3D and MISHKA should agree to the level global codes reach in the table above (0.5 to 3.5 % [read]) [inference]. This step also fixes the normalisation conversion.
2. **Rippled twin**, ripple amplitudes from 0.01 to 0.04 (this week's cases), same n. The ripple repeats five times around the torus, so the families are {5N}, {5N ± 1} and {5N ± 2} (Cooper, Singleton, Dewar 1996 [read]). n = 10, 15 and 20 share the family that also holds n = 0 and 5; n = 11 or 12 lie in families without n = 0. Run both kinds [inference].
3. **W7-X**, Toon's wout and the SIMSOPT ⟨β⟩ = 5 % file (section 3e), per family.

Quantities: ω² at matching n and family after conversion, and marginal points (β or pressure-gradient scale at which ω² = 0). Marginal points do not depend on ω_A, the density or the kinetic-energy form, so they avoid the normalisation problem of section 2b [inference].

**How PB3D's single-n result relates to a family-coupled result.** Hegna 2014 and Anastopoulos-Tzanis et al. 2019 expand in the 3-D amplitude δ [read]: a purely 3-D perturbation changes the eigenvalue only at order δ², through an averaged term and a sum over the sidebands n ± N; for the most unstable mode the sideband sum destabilises. If PB3D's α-independent limit equals a toroidal average of the 3-D coefficients (the JCP 2017 statement suggests it; nobody has checked it in the code), PB3D keeps the averaged term and leaves out the sideband sum; then PB3D and a family-coupled code agree to order δ and differ at order δ² [inference]. The ripple scan of step 2 tests exactly this.

Expected size of the coupling, from the literature:

- perturbation theory is accurate for δB_N/B0 ~ 10⁻⁴ (Anastopoulos-Tzanis, Ham, Wilson, EPS 2019) [read, abstract];
- the sideband sum has denominators ω_n0² − ω_{n±N}²; neighbouring high-n modes have nearly equal growth rates, so the expansion fails once the coupling exceeds their spacing, and the mode then localises on the worst field lines [inference]; near rational surfaces the δ expansion fails (Bird and Hegna 2013) [read, abstract];
- without the sidebands even the qualitative agreement between WKB and TERPSICHORE was lost in Cooper et al. 1996 [read];
- symmetry breaking lowers the critical edge pressure gradient by up to 30 % in CASTOR3D on AUG (Puchmayr 2024) [snippet];
- survey B found no paper that quantifies the difference between a single-n high-n code and a family-coupled code in a stellarator.

For W7-X the 3-D shaping is not small, no δ ordering exists, and differences are expected; their size is unknown. Effort: 1 day install and QAS3, 1 to 2 days for the twins and the conversion, 2 days for the ripple scan, 2 days for W7-X.

### 3b. COBRAVMEC for the n → ∞ limit (plan item L11)

Plan item L11 (local): PB3D's n → ∞ growth rate on a 3-D equilibrium is at most COBRA's worst-field-line value and equals it in axisymmetry; item T8 is the axisymmetric version. The plan lists a STELLOPT clone at `~/Code/stellopt` (local, item C2); whether COBRAVMEC is built there was not checked for this document.

COBRAVMEC (folder `COBRAVMEC` in https://github.com/PrincetonUniversity/STELLOPT, MIT [read]) solves the ballooning equation on one field line (s, α, θ_k) in VMEC coordinates, with error ∝ h⁴ after a variational step and Richardson extrapolation [read]. Both PB3D papers used it [read].

What to run:

1. COBRAVMEC on every surface of PB3D's radial window, for the circular and cbm18a VMEC twins, the rippled twin and W7-X; scan α over one period (8 to 16 lines) and θ_k; take the maximum over α per surface [inference]. Make the line long enough: QAS C82 needed K_w ≈ 65 to 70 helical wells [read].
2. PB3D at n = 10, 20, 40, 80 and, if the resolution allows, 160, same window.
3. Fit ω²(n) = c1 + c2/n, the form Snyder 2009 used to reach n = ∞ with ELITE [read], and repeat with a c3/n² term to see how stable c1 is. Compare c1 with COBRA: equal in axisymmetry, at most COBRA's worst line in 3-D.

Caveats:

- Normalisation: COBRA uses τ_A = (μ0ρR0²/B0²)^½ with B0 on the axis [read]; PB3D's `norm_style = 2` uses a modified B_0 and ψ_0 and states EV_COBRA = EV_PB3D · β/2 (local, `Modules/eq_ops.f90`). Marginal surfaces avoid both.
- Extrapolation: the 1/n form holds only at very high n, critical β against n oscillates at low n, and extrapolation from low n is unreliable (Dewar, Manickam, Grimm, Chance 1981) [read].
- In 3-D the n = ∞ value broadens into a band (Dewar and Glasser 1983) [read]; that band lies above all finite-n eigenvalues found (Cooper 1996) [read], and TERPSICHORE at n = 19 to 25 was "significantly more stable" than the infinite-n limit in QAS C82 (Fu 2000) [read]. In 3-D only the inequality of L11 is a test.
- PB3D is α-independent only for an integration over −∞ to ∞; on a finite parallel range its α-dependence falls linearly with the range (JCP 2017, Fig. 8) [read].

Effort: 0.5 to 1 day to build, 1 to 2 days for scans and fits.

### 3c. DESC ballooning

Usable. DESC (GitHub PlasmaControl/DESC, open source) has a `BallooningStability` objective and an "ideal ballooning lambda" quantity: the infinite-n equation with Dirichlet ends, λ = γ² per surface, maximum over α; the tutorial uses 8 lines on [0, π), 3 toroidal turns and 200 points per turn [read]. The method (Gaur et al. 2025) uses finite differences along ζ and a symmetric tridiagonal generalised eigenproblem [read]. DESC reads VMEC wout files (survey C) [read]; whether the ballooning quantity runs on an equilibrium loaded that way was not checked [inference]. The published DESC–COBRAVMEC comparison matched only the profile shape, after a rescaling by the ratio of maxima [read]. Use DESC as a second n → ∞ code for the unstable surfaces and the worst α, not for absolute growth rates. Effort: 0.5 to 1 day. The ballopt code (Gaur) uses a method nearly identical to COBRAVMEC's [snippet].

### 3d. CAS3D and CASTOR3D through IPP

CAS3D (Greifswald) is the reference global ideal code for W7-X (Nührenberg 1996 [meta]); it uses the full compressible δW, one family at a time, and a Green's-function vacuum that needs no wall [read]. CASTOR3D (Garching) has published peeling-ballooning modes in 3-D tokamaks with magnetic perturbations and in a two-period QA stellarator [read]. Survey B found no public release and no request procedure for either code.

A request would contain:

- the VMEC wout files: circular twin, rippled twins with amplitudes 0.01 to 0.04, the W7-X file PB3D uses;
- for each, n = 10, 15, 20 and one n outside the n = 0 family, PB3D's ω²/ω_A² with the normalisation spelled out (B_0 the (0,0) mode of |B| on the axis, R_0 the axis, flat density), the radial window, fixed boundary;
- the questions: growth rate of the most unstable mode of the family that contains n, its dominant harmonics, the radial and Fourier resolution, the compressibility setting (Γ = 0 if possible) and the normalisation;
- a second stage for CAS3D: free boundary without a wall, so that free minus fixed checks PB3D's 3-D vacuum (plan items C7 and C8, local).

Effort for Toon: 1 day to prepare files and text. By rule 8 this goes out after the local checks.

### 3e. W7-X and other published reference points

| Case | Published result | Code | Equilibrium available | Use for PB3D |
|---|---|---|---|---|
| W7-X standard (EIM) | MHD-stable up to ⟨β⟩ ≈ 5 % (Aleynikova 2022) [read] | not stated in the survey | SIMSOPT test files `input.W7-X_without_coil_ripple_beta0p05_d23p4_tm` and its reference wout (NFP 5, fixed boundary, ⟨β⟩ = 5 % by the name) [read]; DESC example W7-X [read]; VMEC++ example input [snippet] | β scan with VMEC++, PB3D against COBRA marginal β |
| W7-X EIM at ⟨β⟩ = 4.9, 5.4, 6 %, B ≈ 2.3 T | COBRA: locally ballooning-unstable, Mercier-stable (Zhou 2024) [read] | COBRA, M3D-C1 | free-boundary VMEC; availability not stated | PB3D should find high-n instability there |
| W7-X design | ideal ballooning at 4.7 % (Beidler 1990), β limit 4.3 % (Grieger 1992) [read, quoted by Zhou 2024] | — | no | context |
| W7-X low-shear configurations | CAS3D global-mode β limit ≈ 1.7 % (EPS 2023 poster) [snippet] | CAS3D | no | ask IPP |
| W7-X reference equilibria | verified set, "available for internal and external researchers" [read] | VMEC | only through the IPP web service, unreachable from outside IPP [read] | ask IPP |
| NCSX QAS3_LI383 | COBRA stable at β = 4.0 %, unstable at 4.25 %; VVBAL marginal at 4 % [read] | COBRA, TERPSICHORE | DESC example NCSX [read] (same version not checked) | marginal β at high n |
| QAS3_C82 | n = 19, 22, 25 more stable than n = ∞ [read] | TERPSICHORE | not found | context for the n → ∞ gap |
| W7-A, β = 0.11 and 0.21 % | JOREK3D and CASTOR3D ballooning, 6.4 % at most [read] | JOREK3D, CASTOR3D | not stated | expected agreement |
| CASTOR3D two-period QA | peeling-ballooning at n* = 19 to 47, ⟨β⟩ = 2.9 %, edge q = 1.79, R0 = 12.29 m [read] | CASTOR3D | no data statement found [read] | closest stellarator case; ask IPP |
| AUG with magnetic perturbations | 3-D geometry raises edge growth rates (Puchmayr 2023, 2024) [snippet] | CASTOR3D | not public | 3-D tokamak counterpart of PB3D's ripple runs |
| GIGA (Gauss Fusion, 4 periods) | COBRAVMEC and TERPSICHORE on families n = 0, 1, 2 (Lazerson 2026) [read] | — | no public file | the Gauss Fusion inquiry |

PB3D's own published 3-D checks are qualitative: in JCP 2017 PB3D and COBRA both find a CBM18 with 10 % minor-radius variation more unstable than either axisymmetric end case [read], and in the thesis COBRA is "similar but slightly more unstable" on rippled H-mode cases [read].

---

## 4. 3-D vacuum operator checks

Survey C found no published 3-D Chance-style vacuum response for a stellarator with accuracy numbers, so every check below is a comparison with a code, not with a published number.

### 4a. DESC `singular_integral`: the first independent check

**What it computes.** `singular_integral(eval_data, source_data, kernel, interpolator)` in `desc/integrals/singularities.py` (MIT) implements the partition-of-unity polar quadrature of Malhotra et al. 2020, with the kernels `'1_over_r'` (PB3D's G) and `'nr_over_r3'` = n'·(x−x')/|x−x'|³ (PB3D's H) [read]. `best_params` and `best_ratio` choose a patch that is round in real space, which handles anisotropic grids [read]. DESC's own test reaches |2π + ∫ n'·(x−x')/|x−x'|³ dS'| < 6e-7 on its W7-X equilibrium with `LinearGrid(M=50, N=50)` per field period [read]. Install: `pip install desc-opt` [read].

**How to feed PB3D's W7-X boundary and points.**

1. Load the wout file PB3D used (`wout_w7xflip.nc`, local) into DESC; the boundary is then the same Fourier series. Check that DESC's poloidal angle on the boundary is VMEC's, or map between them, before comparing point values [inference].
2. Export from PB3D: collocation points (x, y, z) and their (θ_V, ζ), normals, the G and H matrices, and the normal data ∂Φ/∂n that a computed eigenmode puts on the boundary.
3. Evaluate in DESC on its own regular grid and carry the smooth results to PB3D's points by Fourier interpolation [inference]. Whether DESC accepts arbitrary targets on the surface directly was not checked.

**What to compare.**

- *Operator actions.* G·σ and H·φ for smooth test densities, for instance cos(mθ − nN_fpζ) over PB3D's mode range and the traces of interior multipoles (exact exterior data on the real surface, survey C section 5.3 [inference]). This shows which of G and H carries an error, and where on the surface.
- *The Neumann-to-Dirichlet map on PB3D's own normal data.* Solve the exterior second-kind equation in DESC (GMRES or dense LU, N ≈ 10⁴ [inference, survey C]) for the ∂Φ/∂n of PB3D's eigenmode, and compare Φ at PB3D's points and the vacuum energy ∫ Φ ∂Φ/∂n dS, the quantity that enters δW. The energy needs no interpolation.
- *Sign and factors.* With the normal out of the plasma and Φ decaying outside, (½I − K)Φ = −V ∂Φ/∂n [inference, survey C, checked on the unit sphere]; match DESC's kernel signs and 2π factors to PB3D's H before comparing.

Expected: DESC below 1e-6; PB3D within its own grid spread, 1.2 % (this week's W7-X sequence). A larger difference points to a PB3D defect [inference]. Effort: 1 to 2 days (survey C's estimate).

### 4b. GPEC's Julia Vacuum module and BIEST

**GPEC Julia `Vacuum`** (https://github.com/OpenFUSIONToolkit/GPEC, MIT; `Kernel3D.jl` Apache-2.0, translated from BIEST) [read]. It builds the 3-D single- and double-layer matrices on a regular (θ, ζ) grid: trapezoidal far field, Malhotra polar quadrature in the near field (patch radius 11 points, 20 radial Gauss points, 5th-order Lagrange interpolation), free term as a unit diagonal, field-period folding, stellarator symmetry, and W_v = (4π²/M) E^H (D⁻¹S) E as in Chance 1997; walls "nowall" and "conformal" only [read]. Its tests put the 3-D W_v within 1e-3 of the 2-D Chance path on a circular torus (R0 = 1.7, a = 0.3, 48 × 96) [read]; the Sherwood 2026 abstract reports "less than 0.1%" for an axisymmetric case discretised in 3-D, and for a stellarator only convergence with toroidal resolution [read]. It computes the quantity closest to PB3D's, a vacuum energy matrix in a Fourier basis. Before matrices can be compared, the Fourier basis, the angle convention (straight-field-line angle with ν) and the (2π)² factors must match [inference, survey C]. Order: circular torus against PB3D's style 2 and style 1 on the circular twin (plan items C5 and C6, local), then W7-X. Effort: 2 to 4 days (survey C's estimate).

**BIEST** (https://github.com/dmalhotra/BIEST, Apache-2.0, header-only C++17) [read]. `ExtVacuumField` solves the exterior problem for given B·n; Φ on the surface follows as S0[σ] through `Laplace3D::FxU`, a use the documentation does not describe [inference]. On-surface error 2.6e-7 at 40 × 40 per period and about 1e-10 at 80 × 80 on a five-period test surface [read]; built-in W7-X boundary and VMEC RBC/ZBS input [read]. Its patch is square in index space and it warns at a spacing ratio above 4 [read]; PB3D's lattice has about 8.6 [inference], which matters only if BIEST were run on PB3D's points. Compare the energy ∫ Φ ∂Φ/∂n on BIEST's own grid. Effort: 1 to 2 days including the build (survey C's estimate).

A low-order third check is bempp-cl on a triangulation of the VMEC grid: error ∝ h², 10⁴ to 10⁵ triangles with the FMM assembler, numba backend on a Mac, about 1 day [read for the library; inference for the accuracy].

### 4c. Quadrature literature for the two remaining approximations

**End correction along the open field lines.** PB3D sums along each open line with the trapezoidal rule, and the line ends leave an error of order 1 in the grid spacing (this week). Gregory end corrections (Fornberg 2021 [meta]; Javed and Trefethen 2016 [meta]) change the weights of a few points at each end, so the rule reaches order p + 1 on a non-periodic interval; the weights grow with p, and stabilised versions use extra end points. In PB3D this changes the weights of the first and last two to four points of each line, with no new points and no new kernel evaluations [inference]. Alpert's hybrid Gauss–trapezoidal rules (Alpert 1999 [meta]) add a few nodes near each end, with positive weights, any order, and variants for x^(−1/2) and log x end singularities; in PB3D the new nodes need the density interpolated from the line's grid values, which the cubic interpolant of the near field could supply [inference]. Kapur and Rokhlin 1997 [meta] correct the weights of existing equispaced points near a singularity and near the ends, with no new points, though at high order some weights are large or negative [inference]. Hao et al. 2014 [meta] compare these rules side by side. A different route removes the ends: choose the number of lines and the parallel length so that the end of each line lands on the start of another, as in the twist-and-shift condition of flux-tube codes (Beer, Cowley, Hammett 1995 [meta]). The field-line grid is then a closed, doubly periodic lattice, the trapezoidal rule needs no end correction, and the smooth far field converges fast (Trefethen and Weideman 2014 [read]) [inference, survey C]. It needs a rational approximation of the rotational transform over the line length on the boundary surface [inference].

**Diagonal on the anisotropic lattice.** Wu and Martinsson 2021 [read] correct the punctured trapezoidal rule for the Laplace single layer by one diagonal weight, τ = −Z_A(1)·h, where Z_A is the Epstein zeta function of the local first fundamental form; the diagonal correction alone gives order 3 for the single and the double layer (the double-layer weight uses the second fundamental form and the mean curvature), and a 9-point stencil gives order 5. Their Remark 4 covers unequal spacings h1, h2 with τ = −Z_Ã(1)·√(h1h2), Ã = [[E h1/h2, F], [F, G h2/h1]] [read]. PB3D now sets the G diagonal from the linear potentials (this week); Sun, Klaseboer, Khoo and Chan 2014 [read] apply the same idea row by row with the local normal [inference, survey C]. Wu and Martinsson give that diagonal in closed form from the local metric at each point, with F ≠ 0 because the along-line and between-line directions are not at right angles and h1/h2 ≈ 8.6 [inference]. PB3D can use it as an independent check of its diagonal, or replace the diagonal by it. The method assumes a locally uniform lattice and a periodic or compactly supported density, so the open ends again need separate treatment [inference, survey C]. ZetaTrap3D (MATLAB, MIT) is the reference implementation [read].

Barnett's "5h rule" (Barnett 2014 [meta, snippet]) explains this week's near-field fix: a smooth rule with spacing h is accurate only for targets at least 5h away. On PB3D's W7-X lattice h along a line is about 0.06 R0 while neighbouring lines are 0.007 R0 apart [inference], so neighbouring points sit far inside 5h, and the near field must be integrated along the line.

### 4d. Chance 1997 and NESTOR, for context

Chance 1997 (VACUUM) is axisymmetric [read]: Green's second identity, collocation on the plasma and wall contours, the toroidal Fourier component of 1/|x−x'| as kernel, the log singularity subtracted and integrated analytically against 5-point Lagrange polynomials, Gauss–Legendre for the remainder. Its double-layer self-term comes from the constant-potential identity, the one PB3D uses for the H diagonal; it reports accuracy within 1e-5 for a simple closed wall and uses the symmetry of the vacuum matrices as a routine accuracy test [read]. PB3D can do the same at no cost: the antisymmetric part of its discrete vacuum energy matrix is an error estimate (GPEC Julia warns when it is large [read]) [inference].

**Correction.** Chance, Turnbull and Snyder 2007 (J. Comput. Phys. 221, 330) is titled "Calculation of the vacuum Green's function valid even for high toroidal mode numbers in tokamaks". It is an axisymmetric paper on the toroidal Fourier components G_n of 1/|x−x'| at high n, accurate to 12 to 16 digits, and does not treat 3-D surfaces [read]. For PB3D it is a cross-check of the style-2 G_n (dtorh) at high n [inference]. No public VACUUM3D code exists; the public 3-D descendant of VACUUM is the 2026 Julia rewrite inside GPEC [read].

NESTOR (Merkel 1986 [meta]) writes Green's identity on the surface and solves for Φ in a truncated double Fourier basis, after subtracting surrogate kernels built from the first fundamental form; Malhotra et al. 2020 show it converges only as h² on the constant-potential test, on a circular torus and on W7-X [read]. VMEC++ contains a C++ NESTOR, and survey C found Schilling 2026 the clearest derivation of its singularity treatment [read]. CAS3D's vacuum is credited to P. Merkel, with no paper on its discretisation or accuracy [read]. TERPSICHORE's pseudo-plasma vacuum failed for walls beyond about 1.7 to 2.7 minor radii [snippet], so it cannot check a no-wall Neumann-to-Dirichlet map [inference].

---

## 5. Priority list

| # | Action | Reference quantity | Code or paper | Equilibrium | Expected agreement | Effort (days) | Who |
|---|---|---|---|---|---|---|---|
| 1 | Settle the cbm18a gap with MISHKA: a HELENA cbm18a that ends at the pedestal foot (no pseudo-vacuum on either side), MISHKA with a far wall against PB3D free; both codes at 400 points at n = 15, 20 on ε = 0.33 | −ω²/ω_A² at n = 10, 15, 20 | MISHKA-1 | cbm18a cut at the pedestal foot; circular ε = 0.33 | within 1.2 %, as on the circular cases, if the pseudo-vacuum's rational surfaces explain the gap [inference] | 1 to 2 | assistant (HELENA input needs Toon's EPS, BETAP, XIAB for the cut plasma) |
| 2 | PB3D cbm18a free boundary at ψ_N = 0.701, n = 5 to 40, each converged radially; digitise the vector figures; compare γ(n)/γ(20) and the first unstable n at once | γ/ω_A against n | Ferraro 2010 Fig. 6; Squires 2008; Dudson 2011 Fig. 2 | cbm18a (HELENA) | a few % at n ≥ 10 (ELITE–GATO 0.005, BOUT++–ELITE 2.5 % [read]); absolute values wait for ω_A | 2 | assistant; PDF downloads after Toon's approval |
| 3 | Download and convert Pamela 2020's CBM18 data | γ in s⁻¹, n = 1 to 20 | JOREK, MISHKA, CASTOR (Zenodo 3971940) | CBM18 | "reasonable" between codes [read]; Toon's MISHKA should match theirs in shape | 0.5 | Toon approves the download; assistant |
| 4 | DESC `singular_integral` check of the 3-D vacuum | G·σ, H·φ on test densities; Φ and ∫ Φ ∂Φ/∂n on PB3D's normal data | DESC (Malhotra 2020 quadrature) | W7-X (PB3D's wout); circular twin | DESC < 1e-6; PB3D within its 1.2 % grid spread | 1 to 2 (survey C) | assistant |
| 5 | COBRAVMEC n → ∞ limit (L11, T8), DESC ballooning as cross-check | ω² at n → ∞; unstable surfaces | COBRAVMEC (STELLOPT); DESC | circular and cbm18a VMEC twins; rippled twin; W7-X | axisymmetric: PB3D's 1/n fit equals COBRA within the fit error; 3-D: PB3D at most COBRA's worst line (L11, local) | 2 to 3 | assistant |
| 6 | TERPSICHORE: install, QAS3 test, axisymmetric twins, ripple scan, W7-X | ω² per family at n = 10 to 20; marginal β | TERPSICHORE (EPFL GitLab, 2025) | QAS3; twins; rippled twin 0.01 to 0.04; W7-X | twins 0.5 to 3.5 % [read]; ripple: difference ∝ δ² [inference]; W7-X unknown | 6 to 7 | assistant; Toon approves the clone and reviews |
| 7 | End corrections, or a closed field-line lattice, in the style-1 BEM; Wu–Martinsson diagonal as a check | measured order along the lines | Fornberg 2021; Alpert 1999; Kapur–Rokhlin 1997; Beer et al. 1995; Wu–Martinsson 2021 | circular twin; W7-X | order ≥ 2 along the lines [inference] | 2 to 4 | Toon picks the method; assistant implements |
| 8 | GPEC Julia `Vacuum`, then BIEST | W_v matrix; Neumann-to-Dirichlet map | GPEC Julia (2026); BIEST | circular torus; W7-X | GPEC 3-D against 2-D < 1e-3; BIEST 1e-7 to 1e-10 [read] | 2 to 4, plus 1 to 2 (survey C) | assistant |
| 9 | ELITE request (section 2c) | γ/ω_A at n = 5 to 40, Γ = 0 and 5/3; definition of ω_A; TOQ files | ELITE | PB3D's HELENA cbm18a cut at 0.701 and 0.75; cbm18_dens8 | a few % at n ≥ 10 [read, intercomparisons] | 0.5 to prepare | Toon with the ELITE group |
| 10 | IPP request (section 3d) | growth rate per family; free minus fixed without a wall | CAS3D (Greifswald); CASTOR3D (Garching) | circular twin; rippled twin; W7-X | 1 to 4 % between global codes [read]; single-n effect unknown | 1 to prepare | Toon with the IPP groups |

Items 1 to 8 run on the Mac, on Snellius or on CI and need nobody outside; items 2, 3 and 6 need Toon's approval of a download. Item 9's question about ω_A could go out before item 2 ends. Items 9 and 10 otherwise follow rule 8 of the plan.

---

## 6. References

### Section 2 (axisymmetric)

- Ferraro, Jardin, Snyder, "Ideal and resistive edge stability calculations with M3D-C1", Phys. Plasmas 17, 102508 (2010). doi:10.1063/1.3492727; PDF https://m3dc1.pppl.gov/Papers/ferraro_PP_10.pdf [read]
- Squires, Kruger, Hegna, Held, Snyder, Sovinec, Zhu, "Dynamic behavior of peeling-ballooning modes in a shifted-circle tokamak equilibrium", APS-DPP poster, Dallas (2008). https://w3.pppl.gov/cemm/APS2008/Squires.pdf [read]
- Dudson, Xu, Umansky, Wilson, Snyder, "Simulation of edge localized modes using BOUT++", Plasma Phys. Control. Fusion 53, 054005 (2011). arXiv:1008.4554 [read] (volume from a repository record [snippet]; DOI not confirmed)
- Dudson, Umansky, Xu, Snyder, Wilson, "BOUT++: a framework for parallel plasma fluid simulations", Comput. Phys. Commun. 180, 1467 (2009). doi:10.1016/j.cpc.2009.03.008; arXiv:0810.5757 [read]
- Burke, Kruger, Hegna, Zhu, Snyder, Sovinec, Howell, "Edge localized linear ideal magnetohydrodynamic instability studies in an extended-magnetohydrodynamic code", Phys. Plasmas 17, 032103 (2010). doi:10.1063/1.3309732 [snippet, abstract only]
- Pamela, Bhole, Huijsmans, Nkonga, Hoelzl, Krebs, Strumberger, "Extended full-MHD simulation of non-linear instabilities in tokamak plasmas", Phys. Plasmas 27, 102510 (2020). doi:10.1063/5.0018208 [read]; data doi:10.5281/zenodo.3971940 [snippet]
- Hoelzl et al., "The JOREK non-linear extended MHD code and applications to large-scale instabilities and their control in magnetically confined fusion plasmas", Nucl. Fusion 61, 065001 (2021). doi:10.1088/1741-4326/abf99f; arXiv:2011.09120 [read]
- Snyder et al., "Edge localized modes and the pedestal: A model based on coupled peeling-ballooning modes", Phys. Plasmas 9, 2037 (2002). doi:10.1063/1.1449463 [read]
- H. R. Wilson, P. B. Snyder, G. T. A. Huysmans, R. L. Miller, "Numerical studies of edge localized instabilities in tokamaks", Phys. Plasmas 9, 1277–1286 (2002). doi:10.1063/1.1459058 [meta]
- Snyder et al., "Pedestal stability comparison and ITER pedestal prediction", Nucl. Fusion 49, 085035 (2009). doi:10.1088/0029-5515/49/8/085035 [read, preprint EFDA-JET-CP(08)05/29]
- Lunniss, "Modelling eruptions and edge stability in tokamak plasmas", PhD thesis, University of York (2016). https://etheses.whiterose.ac.uk/16638/ [read]
- King, Pankin, Kruger, Snyder, "The impact of collisionality, FLR and parallel closure effects on instabilities in the tokamak pedestal: numerical studies with the NIMROD code", Phys. Plasmas 23, 062123 (2016). doi:10.1063/1.4954302; arXiv:1702.00042 [snippet]
- Mikhailovskii, Huysmans, Kerner, Sharapov, MISHKA-1, Plasma Phys. Rep. 23, 844 (1997) [snippet, not checked]
- Kerner et al., CASTOR, J. Comput. Phys. 142, 271 (1998). doi:10.1006/jcph.1998.5910 [snippet, not read]
- R. L. Dewar, J. Manickam, R. C. Grimm, M. S. Chance, "n-dependence of ballooning instabilities", Nucl. Fusion 21, 493 (1981). doi:10.1088/0029-5515/21/4/008 [read]
- BOUT++ grid files for cbm18_dens8, https://github.com/boutproject/BOUT-dev (examples/elm-pb) [snippet]

### Section 3 (3-D codes and cases)

- T. Weyens, R. Sánchez, G. Huijsmans, A. Loarte, L. García, "PB3D: A new code for edge 3-D ideal linear peeling-ballooning stability", J. Comput. Phys. 330, 997–1009 (2017). doi:10.1016/j.jcp.2016.10.054 [read]
- T. Weyens, "3D Ideal Linear Peeling Ballooning Theory in Magnetic Fusion Devices", PhD thesis, UC3M/TU/e (2016–17). doi:10.13140/RG.2.2.26457.11361 [read]
- C. Schwab, "Ideal magnetohydrodynamics: Global mode analysis of three-dimensional plasma configurations", Phys. Fluids B 5, 3195–3206 (1993). doi:10.1063/1.860656 [meta]
- C. Nührenberg, "Global ideal magnetohydrodynamic stability analysis for the configurational space of Wendelstein 7-X", Phys. Plasmas 3, 2401–2410 (1996). doi:10.1063/1.871924 [meta]
- D. V. Anderson, W. A. Cooper, R. Gruber, S. Merazzi, U. Schwenn, "TERPSICHORE: A Three-Dimensional Ideal Magnetohydrodynamic Stability Program", in Scientific Computing on Supercomputers II, Springer, 159–174 (1990). doi:10.1007/978-1-4613-0659-7_8 [read]; and "Methods for the efficient calculation of the (MHD) magnetohydrodynamic stability properties of magnetically confined fusion plasmas", Int. J. Supercomput. Appl. 4(3), 34–47 (1990). doi:10.1177/109434209000400305 [meta]
- TERPSICHORE repository, https://gitlab.epfl.ch/spc/public/terpsichore (Apache 2.0, 2025) [read]
- E. Strumberger, S. Günter, "CASTOR3D: linear stability studies for 2D and 3D tokamak equilibria", Nucl. Fusion 57, 016032 (2017). doi:10.1088/0029-5515/57/1/016032 [read, preprint]
- E. Strumberger, J. Puchmayr, F. Hindenlang, "Linear stability studies for a quasi-axisymmetric equilibrium including plasma flow, drift and viscosity effects", J. Plasma Phys. 91(6), E157 (2025). doi:10.1017/S0022377825101013 [read]
- J. Puchmayr et al., Nucl. Fusion 63, 086008 (2023), doi:10.1088/1741-4326/acdd12; Nucl. Fusion 64, 086013 (2024), doi:10.1088/1741-4326/ad54d6 [read, citations; content snippet]
- N. Nikulsin et al., Phys. Plasmas 29, 063901 (2022). doi:10.1063/5.0087104; arXiv:2201.12033 [read]
- R. Sanchez, S. P. Hirshman, J. C. Whitson, A. S. Ware, "COBRA: An Optimized Code for Fast Analysis of Ideal Ballooning Stability of Three-Dimensional Magnetic Equilibria", J. Comput. Phys. 161, 576–588 (2000). doi:10.1006/jcph.2000.6514 [read]
- R. Sanchez, S. P. Hirshman, H. V. Wong, Comput. Phys. Commun. 135, 82–92 (2001). doi:10.1016/S0010-4655(00)00225-3 [read]
- R. Sanchez et al., "Ideal MHD stability calculations for compact stellarators", Comput. Phys. Commun. 141, 55–65 (2001). doi:10.1016/S0010-4655(01)00396-4 [read]
- STELLOPT (COBRAVMEC), https://github.com/PrincetonUniversity/STELLOPT [read]
- R. Gaur et al., "Omnigenous stellarators with improved ideal and kinetic ballooning stability", Plasma Phys. Control. Fusion 67, 125015 (2025). doi:10.1088/1361-6587/ae1e9e; arXiv:2410.04576 [read]; DESC tutorial https://desc-docs.readthedocs.io/en/latest/notebooks/tutorials/ideal_ballooning_stability.html [read]
- R. Gaur et al., J. Plasma Phys. 89, 905890518 (2023). doi:10.1017/S0022377823000995; arXiv:2302.07673 [read, abstract]
- W. A. Cooper, D. B. Singleton, R. L. Dewar, "Spectrum of ballooning instabilities in a stellarator", Phys. Plasmas 3, 275 (1996). doi:10.1063/1.871853 [read]
- C. C. Hegna, "Effects of a weakly 3-D equilibrium on ideal magnetohydrodynamic instabilities", Phys. Plasmas 21, 072502 (2014). Preprint https://cptc.wisc.edu/wp-content/uploads/sites/327/2017/09/UW-CPTC_14-1.pdf [read]
- M. S. Anastopoulos-Tzanis et al., "Non-axisymmetric equilibrium and stability using the ELITE stability code", Nucl. Fusion 59, 126028 (2019). doi:10.1088/1741-4326/ab40ef [read]; EPS 2019, UKAEA-CCFE-CP(19)51 [read, abstract]
- T. M. Bird, C. C. Hegna, Nucl. Fusion 53, 013004 (2013). doi:10.1088/0029-5515/53/1/013004; arXiv:1208.4313 [read, abstract]
- R. L. Dewar, A. H. Glasser, "Ballooning mode spectrum in general toroidal systems", Phys. Fluids 26, 3038 (1983). doi:10.1063/1.864028 [read]
- M. Redi et al., APS-DPP 2000 slides, https://ncsx.pppl.gov/Scientificconf/2000_10_23_APS/redi.pdf [read]; G. Y. Fu, IAEA FEC 2000 slides, https://ncsx.pppl.gov/Scientificconf/2000_10_02_IAEA/fu_talk.pdf [read]; G. Y. Fu, NCSX PVR slides (2001), https://ncsx.pppl.gov/pvr/PVR_Presentations/Fu_PVR.pdf [read]
- K. Aleynikova, A. Zocco, J. Geiger, J. Plasma Phys. 88, 905880411 (2022). doi:10.1017/S0022377822000745 [read]
- Y. Zhou, K. Aleynikova, C. Liu, N. M. Ferraro, "Benign saturation of ideal ballooning instability in a high-performance stellarator", Phys. Rev. Lett. 133, 135102 (2024). doi:10.1103/PhysRevLett.133.135102; arXiv:2404.01966 [read]
- J. Geiger, C. Nührenberg et al., EPS 2023 poster P4.064 [snippet]
- S. A. Lazerson et al., "The fixed boundary plasma equilibrium basis for a one Gigawatt electric stellarator power plant", arXiv:2607.09346 (2026) [read]
- Y. Nakamura et al., "Comparison of the calculations of the stability properties of a specific stellarator equilibrium with different MHD stability codes", J. Comput. Phys. 128, 43–57 (1996). doi:10.1006/jcph.1996.0195 [read, abstract]
- SIMSOPT test files, https://github.com/hiddenSymmetries/simsopt/tree/master/tests/test_files [read]; DESC examples, https://github.com/PlasmaControl/DESC/tree/master/desc/examples [read]; IPP W7-X equilibria page, https://www.ipp.mpg.de/4631651/Equilibria [read]

### Section 4 (vacuum)

- M. S. Chance, "Vacuum calculations in azimuthally symmetric geometry", Phys. Plasmas 4, 2161–2180 (1997). doi:10.1063/1.872380 [read]
- M. S. Chance, A. D. Turnbull, P. B. Snyder, "Calculation of the vacuum Green's function valid even for high toroidal mode numbers in tokamaks", J. Comput. Phys. 221, 330–348 (2007). doi:10.1016/j.jcp.2006.06.025 [read]
- P. Merkel, "An integral equation technique for the exterior and interior Neumann problem in toroidal regions", J. Comput. Phys. 66, 83–98 (1986). doi:10.1016/0021-9991(86)90055-0 [meta]
- J. Schilling, "The Numerics of VMEC++", arXiv:2502.04374v3 (2026); code https://github.com/proximafusion/vmecpp [read]
- GPEC Julia, https://github.com/OpenFUSIONToolkit/GPEC, docs https://openfusiontoolkit.github.io/GPEC/dev/vacuum/ [read]; J. Halpern et al., Sherwood Theory Conference 2026 abstract, https://sherwoodtheory.org/sw2026/abstract.php?view=247 [read]
- D. Malhotra, A. J. Cerfon, M. O'Neil, E. Toler, "Efficient high-order singular quadrature schemes in magnetic fusion", Plasma Phys. Control. Fusion 62, 024004 (2020). doi:10.1088/1361-6587/ab57f4; arXiv:1909.07417 [read]
- D. Malhotra, A. Cerfon, L.-M. Imbert-Gérard, M. O'Neil, "Taylor states in stellarators: a fast high-order boundary integral solver", J. Comput. Phys. 397, 108791 (2019). doi:10.1016/j.jcp.2019.06.067; arXiv:1902.01205 [read, abstract]; BIEST https://github.com/dmalhotra/BIEST, https://biest.readthedocs.io [read]
- DESC `desc/integrals/singularities.py`, https://github.com/PlasmaControl/DESC [read]; R. Conlin et al., "High order free boundary MHD equilibria in DESC", arXiv:2412.05680 (2024) [read, abstract]; D. Panici et al., J. Plasma Phys. 89 (2023), doi:10.1017/S0022377823000272 [meta]
- B. Wu, P.-G. Martinsson, "Corrected trapezoidal rules for boundary integral equations in three dimensions", Numer. Math. 149, 1025–1071 (2021). doi:10.1007/s00211-021-01244-1; arXiv:2007.02512 [read]; ZetaTrap3D https://github.com/bobbielf2/ZetaTrap3D [read]
- B. K. Alpert, "Hybrid Gauss-trapezoidal quadrature rules", SIAM J. Sci. Comput. 20, 1551–1584 (1999). doi:10.1137/S1064827597325141 [meta]
- S. Kapur, V. Rokhlin, "High-order corrected trapezoidal quadrature rules for singular functions", SIAM J. Numer. Anal. 34, 1331–1356 (1997). doi:10.1137/S0036142995287847 [meta]
- B. Fornberg, "Improving the accuracy of the trapezoidal rule", SIAM Rev. 63, 167–180 (2021). doi:10.1137/18M1229353 [meta]; M. Javed, L. N. Trefethen, Numer. Math. 132, 201–216 (2016), doi:10.1007/s00211-015-0713-x [meta]
- S. Hao, A. H. Barnett, P. G. Martinsson, P. Young, Adv. Comput. Math. 40, 245–272 (2014). doi:10.1007/s10444-013-9306-3 [meta]
- A. H. Barnett, SIAM J. Sci. Comput. 36, A427–A451 (2014). doi:10.1137/120900253; arXiv:1310.5352 [meta]
- Q. Sun, E. Klaseboer, B. C. Khoo, D. Y. C. Chan, "A robust and non-singular formulation of the boundary integral method for the potential problem", Eng. Anal. Bound. Elem. 43, 117–123 (2014). doi:10.1016/j.enganabound.2014.03.010; arXiv:1902.04168 [read]
- M. A. Beer, S. C. Cowley, G. W. Hammett, "Field-aligned coordinates for nonlinear simulations of tokamak turbulence", Phys. Plasmas 2, 2687–2700 (1995). doi:10.1063/1.871232 [meta]
- L. N. Trefethen, J. A. C. Weideman, "The exponentially convergent trapezoidal rule", SIAM Rev. 56, 385–458 (2014). doi:10.1137/130932132 [read]
- T. Betcke, M. W. Scroggs, "Bempp-cl: a fast Python based just-in-time compiling boundary element library", J. Open Source Softw. 6, 2879 (2021). doi:10.21105/joss.02879 [read]

Compiled 31 October 2026 from the three surveys in Documentation/benchmarking/, machine-assisted; verification tags as in the surveys.
