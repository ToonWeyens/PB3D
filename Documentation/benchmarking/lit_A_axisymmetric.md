# Axisymmetric peeling-ballooning benchmarks for PB3D: literature survey A

Survey date: 31 October 2026. Scope: the circular cbm18 benchmark family and the published
tokamak (2-D) edge-stability code comparisons, as references for PB3D's cbm18a runs.

**How each statement was checked** (tag after each claim):

- **[R]** I read the text myself (PDF converted to text, or the page itself).
- **[F]** I read the value off a figure by eye, from a 300 to 500 dpi rendering of the PDF page. Accuracy about ±0.005 in γ/ω_A.
- **[W]** The page was loaded, but I saw it only through the fetch tool's summary, not the full text.
- **[S]** Search-engine snippet or abstract only; the source itself was not opened.
- **[L]** Local file on this machine (PB3D repository or the MISHKA source), not literature.
- **[I]** My own inference or arithmetic; not stated in any source.

---

## 0. What this gives the cbm18a benchmark

1. **The three most useful growth-rate curves for cbm18_dens8** are figures, not tables:
   Ferraro, Jardin and Snyder 2010, Fig. 6 (ELITE with Γ = 0 and Γ = 5/3, GATO, M3D-C1, n up to 40);
   the Squires et al. 2008 APS poster (ELITE, GATO and NIMROD for n = 1 to 20, plus an ELITE scan of the
   plasma-vacuum cutoff); Dudson et al. 2011, Fig. 2 (ELITE, GATO and BOUT++, n up to 60).
   All three figures are vector graphics inside the PDFs (no raster image on those pages, checked with
   `pdfimages -list`), so the marker coordinates can be read from the PDF drawing commands exactly instead of by eye [R].
2. **Printed numbers** exist for two points only: ELITE γ/ω_A = 0.239 at n = 20 (Dudson et al. 2009, CPC; the
   equilibrium is not named there, but the value matches the cbm18_dens8 ELITE curve at n = 20) [R, I], and
   γ/ω_A ≈ 0.132 at n = 8 for both ELITE and NIMROD (Squires 2008 poster) [R].
3. **Numbers as data files**: the JOREK full-MHD paper (Pamela et al. 2020) has a Zenodo dataset
   (doi:10.5281/zenodo.3971940, licence CC BY 4.0). Its archive `CBM18.zip` (136.6 kB) contains
   `growth_rates_castor.txt` (742 B), `growth_rates_cbm18_new.txt` (746 B) and `plot_growth_rates.py`,
   i.e. the CBM18 growth rates of JOREK (reduced and full MHD), MISHKA and CASTOR for n = 1 to 20 [W].
   I did not download the archive; downloading it needs your approval.
4. **The normalisation is the first thing to pin down.** The published curves use at least four
   different choices of R and B in ω_A = B/(R√(μ0ρ)) (table in section 2.2). PB3D's cbm18a free-boundary
   anchor at n = 10, ω²/ω_A² = −4.50e−2, is γ/ω_A = 0.212 in the MISHKA normalisation [L, I]; the published
   ELITE value at n = 10 is 0.160 to 0.165 [F]. The ratio, 1.32, is of the size that the choice of R and B
   alone can produce for this equilibrium (section 2.2) [I]. Until the R, B and density used by each
   curve are known, a disagreement of this size says nothing about PB3D.
5. **Other differences to control**, each worth 5 to 30 %: the density profile (flat in all ideal
   benchmarks; a profile that falls in the pedestal raises γ roughly as 1/√ρ at mid-pedestal) [R];
   the adiabatic index (Γ = 0 and Γ = 5/3 differ by 13 to 18 % at n = 10) [F]; the position of the plasma
   boundary (cbm18 was built to be insensitive to it) [R]; the equilibrium mapping (about 10 % in NIMROD
   runs) [R].
6. **For a request to the ELITE group**: ask for cbm18_dens8 at n = 5, 10, 15, 20, 30, 40 with Γ = 0,
   flat density, infinite vacuum and the standard cutoff ψ_N = 0.701, and ask for the R, B and ρ that
   define their ω_A. ELITE also reads HELENA equilibria (it was run on HELENA files in Lunniss 2016) [R],
   so a second run on PB3D's own HELENA cbm18a file would separate mapping effects from code effects.
7. **MISHKA has run CBM18 before**: Pamela et al. 2020 compare JOREK with MISHKA and CASTOR on CBM18 for
   n = 1 to 20 [R]. So MISHKA can handle this equilibrium in some form; the paper does not say where the
   plasma boundary was placed. A HELENA run with its boundary at the pedestal foot (ψ_N = 0.70) and
   MISHKA's own vacuum (large `RWALL`) is the direct counterpart of a PB3D free-boundary run cut at 0.70 [I].

---

## 1. The cbm18 equilibrium family

### 1.1 Who made it and where it is described

- **Origin.** P. B. Snyder (General Atomics) made the cbm18 series with the TOQ equilibrium code [R:
  Dudson 2011, sec. 3; Squires 2008]. TOQ is cited as R. L. Miller and J. W. Van Dam, Nucl. Fusion 27, 2101 (1987)
  [R: Xia 2010 talk]. No journal paper defines cbm18 by itself. The ELITE test cases of the same origin
  ("circa", "dbm8", "dbm9") are documented in an unpublished "A Practical User's Guide for ELITE"
  (Snyder and Wilson, 2006), cited in Lunniss 2016 [R]. That guide is the likely primary description of
  the cbm/dbm cases; I found no public copy.
- **Name.** Ferraro 2010 calls the circular case CBM18 and the elongated case DBM18 [R]; Lunniss 2016
  calls TOQ's D-shaped test cases dbm8 and dbm9 [R]. "cbm" = circular benchmark and "dbm" = D-shaped
  benchmark is my reading of the names; no source spells it out [I].
- **The "dens" suffix.** The Squires 2008 poster calls cbm18 a TOQ-generated series of equilibria that
  scans across the stability boundary, and describes dens6 as a similar equilibrium with lower pedestal
  pressure and lower edge current than dens8 [R]. So the number orders the series by pedestal height
  (dens6 below dens8). No source I found says whether the series was made by scaling the density at
  fixed temperature; the name suggests it, but that is unverified [I]. Xia et al. 2010 write
  "cbm18_dan8", evidently a typing error [R].
- **First published uses.** Squires et al., APS-DPP 2008 poster (NIMROD) [R]; Burke et al. 2010, PoP
  (NIMROD, ELITE, GATO; the abstract speaks of two shifted-circle equilibria, which the poster suggests
  are dens8 and dens6) [W, I]; Ferraro et al. 2010 (M3D-C1) [R]; Dudson et al. 2011 (BOUT++) [R];
  Pamela et al. 2020 and Hoelzl et al. 2021 (JOREK) [R].

### 1.2 Parameters of cbm18_dens8

| Quantity | Value | Source |
|---|---|---|
| Cross-section | circular, shifted circle, no X-point | Dudson 2011 [R]; Ferraro 2010 [R] |
| Edge safety factor q_a | 3.03 | Dudson 2011 [R]; PB3D cbm18a gives 3.02 at ψ_N = 0.70 [L] |
| q on axis | 1.05 | PB3D cbm18a [L] |
| Normalised beta β_N | 1.51 | Dudson 2011 [R] (Pamela 2020 quotes about 1.3 %, with its own definition of β_N [R]) |
| Pedestal toroidal beta β_t0 | 1.9 × 10⁻² | Dudson 2011 [R] |
| Pedestal width | L_ped/a = 0.049 | Dudson 2011 [R] |
| Pedestal width (FWHM of ⟨α⟩) | ΔΨ ≈ 0.12, with Ψ normalised to the pedestal foot | Ferraro 2010 [R] |
| Plasma edge (pedestal foot) | ψ_N = 0.701 of the TOQ boundary flux ("standard cutoff") | Squires 2008 [R] |
| Same in pedestal-foot normalisation | outer boundary at Ψ_w = 10/7 | Ferraro 2010 [R] |
| p′ and J_∥ at the plasma edge | both zero; zero pressure and current outside | Squires 2008 [R] |
| Outer (TOQ) boundary | circle with centre R = 3 m and radius 2 m; M3D-C1 puts its conducting wall there | Ferraro 2010 [R] |
| Magnetic axis, plasma minor radius | R_mag ≈ 3.61 m; plasma edge ≈ 1.20 m radius about R ≈ 3.49 m | [I] from PB3D's cbm18a ratios (HELENA boundary 0.555 R_mag about 0.832 R_mag; edge 0.333 R_mag about 0.969 R_mag) [L] and the 2 m circle about 3 m |
| Field | not stated for dens8. The poster gives B0 = 2 T, R0 = 3 m, a = 1 m and n = 1.06 × 10²⁰ m⁻³ for a modified TOQ shifted-circle equilibrium with edge currents set to zero, probably a different case; the ELITE label R/B = 1.5 fits R = 3 m and B = 2 T | Squires 2008 [R], [I] |
| Pressure | about 3 × 10⁴ Pa at ψ = 0.4; J_∥ peaks near 5.5 × 10⁵ A/m² at ψ ≈ 0.85 (pedestal-foot normalisation) | Dudson 2011, Fig. 1 [F] |
| First unstable n | n = 4 | Pamela 2020 [R]; GATO in Squires 2008 [F] |
| Density | not part of the equilibrium; each code chose its own (section 2.2) | |

The cbm18 equilibrium was built so that the growth rates change little with the position of the
plasma-vacuum interface: the ELITE cutoff scan in the Squires poster (ψ_N = 0.675, 0.688, 0.701, 0.714,
0.727) shows curves that differ visibly only below n ≈ 10 and merge above [R, F].

**Relation to PB3D's cbm18a [L, I].** PB3D's fixture has the pressure gradient ending at ψ_N = 0.70 and
q = 3.02 there, a pseudo-vacuum with no pressure gradient and no current out to the HELENA boundary at
ψ_N = 1 (q = 25), and an outer boundary that is the 2 m circle about R = 3 m if R_mag ≈ 3.61 m. This
matches the TOQ domain of cbm18_dens8 (cutoff 0.701, q_a = 3.03, outer boundary at Ψ_w = 10/7).
Which member of the series cbm18a is has not been established; q_a and the cutoff point to dens8.

### 1.3 Public files

- **BOUT++ grid files (netCDF), public on GitHub** [W]:
  `examples/elm-pb/cbm18_dens8.grid_nx68ny64.nc`, `examples/elm-pb/cbm18_8_y064_x516_090309.nc`,
  the same two files in `examples/tokamak-grids/`, and `cbm18_dens8.grid_nx68ny64.nc` in
  `examples/6field-simple/`, all in https://github.com/boutproject/BOUT-dev.
  These hold the equilibrium on BOUT++'s field-aligned grid (R, Z, ψ, pressure, J_∥ and metric).
  Dudson 2011 says the BOUT++ domain covers 0.4 < ψ_N < 1.2 in pedestal-foot normalisation [R], so the
  core is probably missing [I]. I did not open the files.
- **No EQDSK, TOQ or ELITE input file for cbm18 found in public.** Not in the M3D-C1, NIMROD, JOREK
  or GPEC material I could reach. The JOREK Zenodo archive (`CBM18.zip`) lists only growth rates,
  plots and a plot script, no equilibrium [W].
- **BOUT++ example growth rates are not ideal-MHD references.** The elm-pb example README quotes
  0.275 (68 radial points) and 0.187 (516 points) "× Alfvén time", but the input deck has diamagnetic
  effects, Lundquist number 10⁸ and hyper-resistivity 10⁻⁴ switched on [W].

### 1.4 Related TOQ benchmark equilibria (same origin)

- **DBM18**: as CBM18 but elongation 1.5, triangularity 0.2, narrower pedestal ΔΨ ≈ 0.08 [R: Ferraro 2010].
- **MEUDAS1 (Meudas-1)**: JT-60U-like diverted case, pedestal ΔΨ ≈ 0.06, wall close to the pedestal foot [R: Ferraro 2010]; a 1025 × 1025 EFIT version "meudas1025" exists [R: Lunniss 2016].
- **circa**: circular, q95 = 1.80, β_N = 1.88, tanh pedestals [R: Lunniss 2016].
- **dbm9**: aspect ratio 3, elongation 1.5, triangularity 0.2, q95 = 2.38, β_N = 2.83; built for the original ELITE versus GATO benchmark [R: Lunniss 2016].
- **dbm8**: aspect ratio 3, elongation 1.5, upper triangularity 0, lower 0.3 [R: Lunniss 2016].

---

## 2. Published growth-rate curves for cbm18

### 2.1 The curves

**(a) Ferraro, Jardin, Snyder, "Ideal and resistive edge stability calculations with M3D-C1",
Phys. Plasmas 17, 102508 (2010), doi:10.1063/1.3492727.**
PDF: https://m3dc1.pppl.gov/Papers/ferraro_PP_10.pdf ; preprint: https://w3.pppl.gov/cemm/Project/ELMbenchmark.pdf [R, F]
Fig. 6 shows γ/ω_A against n for CBM18: M3D-C1 with Γ = 5/3, 0 and 1 at n = 4, 5, 7, 10, 15, 20, 30, 40;
ELITE with Γ = 0 and Γ = 5/3 as curves from n ≈ 5 to 40; GATO with Γ = 0 and Γ = 5/3 as curves up to n = 10.
ω_A = B0/(R0√(4πρ_p)) with B0 and R0 at the magnetic axis and ρ_p the mass density at the pedestal top;
density uniform inside the plasma and 100 times lower outside (in M3D-C1). ELITE has an infinite vacuum;
M3D-C1 has a conducting wall at the 2 m circle; the GATO boundary is not stated. Figure only, vector
graphics, easy to digitise. DBM18 and MEUDAS1 curves are in the same figure.

**(b) Squires, Kruger, Hegna, Held, Snyder, Sovinec, Zhu, "Dynamic behavior of peeling-ballooning
modes in a shifted-circle tokamak equilibrium", APS-DPP poster, Dallas, November 2008.**
https://w3.pppl.gov/cemm/APS2008/Squires.pdf [R, F]
Slide 9: cbm18_dens8 with standard cutoff ψ_N = 0.701; GATO (480 × 960 grid), ELITE labelled "R/B = 1.5",
NIMROD, n = 1 to 20; the three agree within about 0.01 from n = 5 up, ELITE and GATO within about 0.005. A second panel is the ELITE cutoff
scan (n ≈ 5 to 30). A third panel is NIMROD with and without the resistive halo; the halo raises γ.
NIMROD used 1.06 × 10²⁰ m⁻³ inside and a hundred times less in the halo. Slide 12: NIMROD γτ_A for dens6
(n = 11 to 20 converged, 0.05 to 0.143; n ≤ 6 slow oscillating modes). Figures only; slide 9 is vector graphics.

**(c) Dudson, Xu, Umansky, Wilson, Snyder, "Simulation of edge localized modes using BOUT++",
Plasma Phys. Control. Fusion 53, 054005 (2011); arXiv:1008.4554.**
https://arxiv.org/abs/1008.4554 [R, F] (journal volume and number confirmed via a repository record [S]; DOI not confirmed)
Fig. 2: cbm18_dens8, BOUT++ reduced ideal MHD (stars), ELITE "R/B = 1.5" (open circles, n up to 60),
GATO 480 × 960 (filled squares, n up to 10). ω_A = V_A/R; flat density 10¹⁹ m⁻³, so τ_A = 0.37 μs.
BOUT++ treats the outside as a low-pressure plasma. BOUT++ rises above ELITE above n ≈ 30 because the
radial grid is fixed. Figure only, vector graphics. ELITE reaches about 0.30 at n = 60 [F].

**(d) Dudson, Umansky, Xu, Snyder, Wilson, "BOUT++: a framework for parallel plasma fluid simulations",
Comput. Phys. Commun. 180, 1467 (2009), doi:10.1016/j.cpc.2009.03.008; arXiv:0810.5757.**
https://arxiv.org/abs/0810.5757 [R]
At n = 20: BOUT++ γ = 0.245 ω_A, ELITE γ = 0.239 ω_A, with ω_A = V_A/R and V_A = B0/√(μ0ρ0), constant
density. Numbers printed in the text. The test case is a strongly ballooning shifted circle and is not
named; the ELITE value equals the cbm18_dens8 ELITE curve at n = 20 [I].

**(e) Burke, Kruger, Hegna, Zhu, Snyder, Sovinec, Howell, "Edge localized linear ideal
magnetohydrodynamic instability studies in an extended-magnetohydrodynamic code",
Phys. Plasmas 17, 032103 (2010), doi:10.1063/1.3309732.**
https://pubs.aip.org/aip/pop/article/106661 [W, abstract only]
NIMROD against ELITE and GATO on two ballooning-unstable shifted-circle equilibria, n = 1 to 20.
NIMROD reproduces the marginal n and the growth spectrum if the outside is a resistive "halo" with
critical Lundquist numbers of about 10⁸ (core) and 10³ (halo), sharp transitions and a little
dissipation. Ferraro 2010 cites it for two further results: ideal behaviour needs ρ_v/ρ_p ≲ 10⁻²,
S_p ≳ 10⁸, S_v ≲ 10³, and small differences in how CBM18 is mapped onto the grid change γ by about 10 % [R].
Full text not read; whether it tabulates numbers is unknown.

**(f) Pamela, Bhole, Huijsmans, Nkonga, Hoelzl, Krebs, Strumberger, "Extended full-MHD simulation of
non-linear instabilities in tokamak plasmas", Phys. Plasmas 27, 102510 (2020), doi:10.1063/5.0018208.**
Open PDF: https://scientific-publications.ukaea.uk/wp-content/uploads/5-0018208.PDF [R];
data: https://doi.org/10.5281/zenodo.3971940 [W]
Fig. 4: CBM18, JOREK reduced and full MHD, MISHKA (incompressible ideal) and CASTOR (ideal and resistive,
η0 = 6 × 10⁻⁶ Ω m), n = 1 to 20, γ in 10⁶ s⁻¹ for ρ = 10¹⁹ m⁻³ (log axis). First unstable mode n = 4;
"reasonable" agreement. The figure is a raster image, but the numbers are in the Zenodo archive
(section 0, item 3). To convert: γ/ω_A = γ R√(μ0ρ)/B. If R/B = 1.5 m/T and the ions are deuterium,
1/ω_A = 0.31 μs, and the n = 20 value of about 0.75 to 0.8 × 10⁶ s⁻¹ becomes 0.23 to 0.25 [F, I].

**(g) Hoelzl et al., "The JOREK non-linear extended MHD code and applications to large-scale
instabilities and their control in magnetically confined fusion plasmas", Nucl. Fusion 61, 065001 (2021),
doi:10.1088/1741-4326/abf99f; arXiv:2011.09120.** [R]
Uses CBM18 at n = 20 as the ballooning convergence test: the growth-rate error falls as
(√(n_flux n_θ))⁻⁵ in reduced MHD, only third order in full MHD. No cross-code numbers for CBM18.

**(h) Xia, Xu, Dudson, Li, "Five-field peeling-ballooning modes simulation ... with BOUT++",
APS-DPP 2010 talk, LLNL-PRES-512837.** https://bout.llnl.gov/sites/bout/files/APS-xia_web.pdf [R]
Same pressure ("cbm18_dan8"), five density profiles; in ideal MHD one non-uniform density case grows
6.2 % faster than the constant-density case. Shows the size of density effects; no ideal reference curve.

**Not found:** any cbm18 curve from DCON/GPEC, MARG2D, KINX, ILSA, CASTOR3D or MISHKA other than (f).

### 2.2 Normalisations used

| Source | Stated ω_A | R/B used | Density |
|---|---|---|---|
| Ferraro 2010 | B0/(R0√(4πρ_p)), B0 and R0 at the magnetic axis | not given | uniform in plasma, step down ×100 outside [R] |
| Squires 2008, Dudson 2011 (ELITE) | label "R/B = 1.5" | 1.5 as printed, probably m/T, i.e. R = 3 m and B = 2 T [I] | flat [I] |
| Dudson 2011 (BOUT++) | V_A/R; τ_A = 0.37 μs at n = 10¹⁹ m⁻³ | 1.80 m/T if deuterium, 2.55 m/T if hydrogen [I] | flat 10¹⁹ m⁻³ [R] |
| Dudson 2009 | V_A/R with V_A = B0/√(μ0ρ0) | not given | constant [R] |
| ELITE default (Lunniss 2016) | B/(R√(μ0ρ_a)), R = mean of innermost and outermost R of the boundary surface, B = f/R there, ρ_a at the edge | for cbm18: R ≈ 3.49 m [I] | eigenvalue is γ²4πρ, so ρ matters only if it varies [R] |
| MISHKA and PB3D | axis values of B and R, central density | R_mag ≈ 3.61 m [I] | MISHKA: set by `DSURF`, `IDPOW`, `IEQ` [L] |
| Pamela 2020 (JOREK) | physical units, s⁻¹ | needs R, B and ion mass | 10¹⁹ m⁻³ for the plot [R] |
| Squires 2008 (NIMROD) | γ/ω_A (dens8), γτ_A (dens6) | not given | 1.06 × 10²⁰ m⁻³ inside [R] |

Ferraro's ELITE Γ = 0 curve and the Squires ELITE "R/B = 1.5" points agree within about 0.005 at n = 10
and 20 [F], so the two figures share one normalisation. Whether that is Ferraro's axis-based formula or
R/B = 1.5 m/T cannot be decided from the papers [I].

**Size of the effect [I].** With the same mass density, γ/ω_A in PB3D's units equals the published value
times (R_mag/B_mag)/(R/B)_published. If the published curves use R/B = 1.5 m/T and the vacuum field is
2 T at R = 3 m (f = 6 T m), then B at R_mag ≈ 3.61 m is about 1.66 T and the factor is about 1.45. With
B_mag nearer 1.8 T it is about 1.3. PB3D's n = 10 value over ELITE's is 1.32. PB3D's own fixture stores
`RMtoG` = 0.832 and `BMtoG` = 1.163 [L]; these, with the R and B the published curves used, give the
factor exactly.

### 2.3 My readings of the curves (by eye, ±0.005 in γ/ω_A)

| n | ELITE "R/B = 1.5" (Squires) | ELITE Γ = 0 (Ferraro) | M3D-C1 Γ = 0 | GATO Γ = 0 (Ferraro) | ELITE Γ = 5/3 (Ferraro) | M3D-C1 Γ = 5/3 |
|---|---|---|---|---|---|---|
| 4 | (GATO 0.015) | – | 0.013 | 0.015 | – | ≈0.01 |
| 5 | 0.059 | – | 0.055 | 0.06 | 0.049 | ≈0.02 |
| 7 | 0.111 | 0.11 | 0.112 | 0.11 | ≈0.085 | ≈0.085 |
| 8 | 0.131 (poster text: 0.132) | – | – | – | – | – |
| 10 | 0.160 | 0.165 | 0.160 | 0.157 | 0.135 | 0.139 |
| 15 | 0.208 | 0.207 | 0.207 | – | ≈0.183 | 0.189 |
| 20 | 0.237 (Dudson 2009 text: 0.239) | 0.239 | 0.236 | – | 0.215 | 0.221 |
| 30 | – | 0.272 | 0.269 | – | 0.250 | ≈0.253 |
| 40 | – | 0.290 | 0.286 | – | 0.271 | 0.270 |

These are for orientation only; extract the vector paths before using them as a reference.

---

## 3. Published edge-stability code benchmarks

**Snyder et al., "Edge localized modes and the pedestal: A model based on coupled peeling-ballooning
modes", Phys. Plasmas 9, 2037 (2002), doi:10.1063/1.1449463.**
https://scientific-publications.ukaea.uk/wp-content/uploads/Published/POPVOL9p2037.pdf [R]
ELITE carries the expansion through two orders in 1/n for accuracy at n ≳ 5. ELITE was benchmarked
against the compressionless MISHKA code for 5 < n < 50 (in Wilson 2002) and against GATO on a
low-squareness DIII-D equilibrium, with good agreement in threshold and mode structure in their overlap
4 ≤ n ≤ 10. No percentage given.

**Wilson, Snyder, Huysmans, Miller, "Numerical studies of edge localized instabilities in tokamaks",
Phys. Plasmas 9, 1277 (2002), doi:10.1063/1.1459058.**
Record: https://eprints.whiterose.ac.uk/7255/ [S; full text not available there]
The ELITE formalism paper and the ELITE versus MISHKA comparison cited above. A snippet of a related
Wilson EPS paper speaks of good agreement for n > 4 [S]. Not read.

**Snyder et al., "Pedestal stability comparison and ITER pedestal prediction", Nucl. Fusion 49, 085035
(2009), doi:10.1088/0029-5515/49/8/085035.**
Preprint EFDA-JET-CP(08)05/29: https://scipub.euro-fusion.org/wp-content/uploads/2014/11/EFDC080529.pdf [R]
The ITPA pedestal group paper. On a D-shaped limiter equilibrium with pedestal profiles, GATO (n = 1 to 19),
MARG2D (n = 5 to 160) and ELITE (n = 4 to 1000) agree well; ELITE at n = 800 to 1000 extrapolated with
γ² = c1 + c2/n agrees with BAL-MSC at n = ∞. On JT-60U-like Meudas-1 cut at the 99.6 % surface, MARG2D and
ELITE agree, also at 99.8 % with a 1025 × 1025 equilibrium; ILSA gives somewhat lower γ. A cutoff beyond
about 99.4 % changes the stability boundary less than measurement error. Figures only.

**Ferraro, Jardin, Snyder 2010 (see 2.1 a).** [R]
M3D-C1 agrees with ELITE and GATO on CBM18 and DBM18; the largest disagreement is MEUDAS1 at low n. Causes
named: infinite vacuum (ELITE) versus a nearby conducting wall (M3D-C1), but moving the wall out by 8 %
changes n = 10 by only a few percent; the plasma-vacuum interface position (MEUDAS1: +30 % when the cutoff
moves from the separatrix to the 99.7 % surface; moving the density and resistivity steps inward together
by 10 % of the distance from pedestal foot to pedestal centre doubles n = 10, moving only the resistivity
step adds 5 %, so the density position, i.e. the inertia, causes most of it); and the equilibrium mapping
(about 10 %, from Burke 2010). For MEUDAS1, MARG2D is within a few percent of ELITE and ILSA about 10 %
lower; ELITE and MARG2D agree within 5 % for cutoffs at 99.6 % and 99.8 %. A non-uniform density raises γ
roughly by the square root of the density drop at mid-pedestal.

**Burke et al. 2010 (see 2.1 e).** [W] NIMROD versus ELITE and GATO, n = 1 to 20, on two shifted circles.

**Dudson et al. 2009 and 2011 (see 2.1 c, d).** [R] BOUT++ versus ELITE: 0.245 against 0.239 at n = 20
(2.5 %). Named causes: different reductions of the MHD equations at finite n, and the vacuum, which
BOUT++ treats as an ideal low-pressure plasma.

**King, Pankin, Kruger, Snyder, "The impact of collisionality, FLR and parallel closure effects on
instabilities in the tokamak pedestal: numerical studies with the NIMROD code", Phys. Plasmas 23, 062123
(2016), doi:10.1063/1.4954302; arXiv:1702.00042.** https://ar5iv.arxiv.org/html/1702.00042 [W]
NIMROD versus ELITE on Meudas-1 (not cbm18): typical difference about 5 %, all below 8 % except 27 % at
n = 4. τ_A = R0/v_A with axis values, R0 = 2.936 m; flat density inside the boundary for the ideal case.

**Pamela et al. 2020 (see 2.1 f).** [R] JOREK reduced and full MHD versus MISHKA and CASTOR on CBM18; also
a JET-like X-point case (reduced versus full MHD, below 5 %) and MAST.

**Hoelzl et al. 2021 (see 2.1 g).** [R] JOREK versus BOUT++ on a circular case with aspect ratio 3.3,
q95 = 3.0, uniform density 10¹⁹ m⁻³: linear agreement within 4 % with and without diamagnetic effects
(from the York PhD thesis of S. Smith, 2019, https://etheses.whiterose.ac.uk/26275/, not read).

**Lunniss, "Modelling eruptions and edge stability in tokamak plasmas", PhD thesis, University of York
(2016).** https://etheses.whiterose.ac.uk/16638/ [R]
Arbitrary-n ELITE versus original ELITE, GATO and MARG2D on circa, dbm9, dbm8 and meudas1025. On circa
with vacuum, the two ELITE versions agree very well for n ≳ 7 and show a small difference for 5 ≤ n ≤ 7;
on dbm8 GATO (n = 1 to 8) and arbitrary-n ELITE agree, better with rising n; on meudas1025 a small
systematic ELITE–MARG2D difference is seen, as in an earlier benchmark. An ideal wall on circa makes
n ≤ 12 stable. All figures, no tables.

**Kwon and Saarelma, EPS 2005 poster P1.054 (ELITE versus MISHKA on JET 55937).**
https://lac913.epfl.ch/epsppd3/2005/pdf/P1_054.pdf [S; the site blocked the fetch]
Snippet: ELITE results agree qualitatively and quantitatively with earlier MISHKA results. No numbers seen.

**Not found:** a 2015 to 2026 paper that benchmarks several ideal edge-stability codes (ELITE, MISHKA,
GATO, DCON, CASTOR) on one tokamak equilibrium with tabulated growth rates. The multi-code benchmarks
since 2015 that I found are about vertical displacement events (JOREK, M3D-C1, NIMROD) or 3-D cases.

---

## 4. The codes

| Code | Equilibrium input | Vacuum and wall | Model and n range | Availability |
|---|---|---|---|---|
| ELITE (Wilson 2002; Snyder 2002) | TOQ files, EFIT g-files and HELENA output, all used in Lunniss 2016 [R] | infinite vacuum, analytic [R: Ferraro 2010, Lunniss 2016]; the arbitrary-n version has an ideal-wall option [R: Lunniss] | ideal MHD expanded to two orders in 1/n, n ≳ 5; Γ = 0 preferred, Γ = 5/3 only to lower order in 1/n [R]; arbitrary-n version since 2016 [R] | not public; General Atomics and University of York (developers); no public repository found [I] |
| GATO (Bernard, Helton, Moore, CPC 24, 377 (1981)) | GA mapping; not checked | not checked | full ideal MHD, global, n = 1 to about 20 [R: Snyder 2009] | General Atomics; page https://fusion.gat.com/THEORY/gato/ cited in Lunniss 2016 [R], not opened |
| MISHKA-1 (Mikhailovskii, Huysmans, Kerner, Sharapov, Plasma Phys. Rep. 23, 844 (1997)) | HELENA mapping (`fort.12`) [L] | vacuum with a conformal wall at `RWALL` (in units of the plasma radius) [L] | full ideal MHD, incompressible [R: Pamela 2020]; all n | not public; from the developers (FOM, JET, CEA lineage) or EUROfusion code platforms [S, I] |
| CASTOR (Kerner et al., J. Comput. Phys. 142, 271 (1998), doi:10.1006/jcph.1998.5910) | HELENA mapping | vacuum and wall (resistive-wall versions exist); not checked | full resistive MHD, all n | as MISHKA [S] |
| CASTOR3D (Strumberger, Günter, Nucl. Fusion 57, 016032 (2017), doi:10.1088/0029-5515/57/1/016032) | 2-D and 3-D equilibria | not checked here; see survey of 3-D codes | linear resistive full MHD | IPP Garching [I] |
| M3D-C1 (Ferraro, Jardin, J. Comput. Phys. 228, 7742 (2009)) | g-file or its own Grad–Shafranov solve [I] | outer region is a resistive plasma; conducting wall at the domain edge [R: Ferraro 2010]; resistive wall since Ferraro et al., Phys. Plasmas 23, 056114 (2016) [R: title only] | initial-value extended MHD, full or reduced; one n at a time in linear mode | read-only access to the Princeton GitHub repository under a licence; https://m3dc1.pppl.gov [S] |
| NIMROD (Glasser et al., PPCF 41, A747 (1999)) | g-file and others, not checked | resistive "halo" plasma outside, wall at domain edge [R: Burke abstract, Squires] | initial-value extended MHD | NIMROD team, https://nimrodteam.org ("Code access" page, terms not read) [W] |
| JOREK (Hoelzl et al., NF 61, 065001 (2021)) | own Grad–Shafranov solve from profiles and boundary | scrape-off layer as plasma; free boundary via STARWALL coupling [S] | reduced and full MHD, initial value | JOREK team, https://www.jorek.eu; not open source as far as the site shows [W] |
| BOUT++ (Dudson et al., CPC 180, 1467 (2009)) | grid netCDF made from g-files or ELITE/GATO files | outside treated as low-pressure plasma [R] | reduced MHD, high-n field-aligned | open source, https://github.com/boutproject/BOUT-dev [W] |
| DCON / GPEC | EFIT g-file and other formats [W] | VACUUM code (Chance), wall options [I, not checked] | Newcomb-type ideal MHD, one n per run, mainly low n | open source, MIT licence: https://github.com/PrincetonUniversity/GPEC ; docs http://princetonuniversity.github.io/GPEC/ [W] |
| MARG2D (Tokuda, Watanabe, PoP 6, 3012 (1999); Aiba et al., CPC 175, 269 (2006)) | not checked | not checked | 2-D Newcomb eigenvalue problem, n up to 160 used [R: Snyder 2009] | JAEA/QST (Japan) [I] |
| KINX (Medvedev et al., PPCF 48, 927 (2006)) | not checked | not checked | ideal MHD with X-point | Keldysh Institute / EPFL [I] |
| ILSA (Konz et al., EPS 2008, Vol. 32D, P2.071) | not checked | not checked | ideal MHD | IPP Garching [I] |

---

## 5. Pitfalls when comparing a high-n code with a full-MHD code at finite n

1. **Validity of the 1/n expansion.** ELITE keeps two orders in 1/n and is stated accurate for n ≳ 5
   (Snyder 2002) [R]. Arbitrary-n ELITE and original ELITE agree very well for n ≳ 7 and differ slightly for
   5 ≤ n ≤ 7 (Lunniss 2016) [R]; GATO and ELITE agree for 4 ≤ n ≤ 10 (Snyder 2002) [R]; on cbm18_dens8 the
   ELITE and GATO curves coincide from n = 5 to 10 within about 0.005 (Squires; Ferraro) [F]. BOUT++ and
   ELITE differ by 2.5 % at n = 20 (Dudson 2009) [R]. PB3D's own measure: 2.1 to 2.3 % against MISHKA at
   n = 10 to 20 on a circular ε = 0.33 external kink [L]. So for a ballooning-dominated case like cbm18,
   n ≥ 10 is safe for a few-percent comparison; below n ≈ 7 expect larger differences.
2. **Extrapolation in n.** Snyder 2009 fits γ² = c1 + c2/n to reach n = ∞ [R]. The same fit applied to
   PB3D minus MISHKA would test whether their difference falls as 1/n [I].
3. **Compressibility.** On CBM18, Γ = 5/3 lowers γ by 13 % (M3D-C1) to 18 % (ELITE) at n = 10 and by about
   6 % at n = 40, against Γ = 0 [F]. ELITE's Γ = 5/3 is accurate only to lower order in 1/n, so it is less reliable at low n
   (Ferraro 2010) [R]. MISHKA-1 is incompressible [R]; compare PB3D with the Γ = 0 curves unless PB3D's
   model is known to correspond to another Γ [I].
4. **Density.** Growth rates scale as 1/√ρ_eff, with ρ_eff about the density at mid-pedestal; a density
   that falls through the pedestal raises γ noticeably (Ferraro 2010) [R]; 6 % in one BOUT++ case (Xia 2010) [R].
   All ideal cbm18 benchmarks used a flat density.
5. **Normalisation.** Section 2.2. Differences of 20 to 45 % are possible for cbm18 from R and B alone [I].
6. **Plasma boundary.** cbm18 has p′ = J_∥ = 0 at the edge and was built to be insensitive to the cutoff [R];
   Dudson 2009 chose a ballooning case for the same reason [R]. Diverted cases are sensitive: +30 % in
   MEUDAS1 for a cutoff at 99.7 % instead of the separatrix, more sensitive to inward than outward shifts
   (Ferraro 2010) [R]; cutoffs beyond about 99.4 % are adequate for stability boundaries (Snyder 2009) [R].
   Low-n peeling modes driven by edge current are strongly stabilised near an X-point, while
   peeling-ballooning thresholds move a few percent (Ferraro 2010, citing Huysmans 2005, Medvedev 2006,
   Webster and Gimblett 2009) [R].
7. **Wall.** PB3D has no wall, so ELITE (infinite vacuum) is the natural partner. M3D-C1's wall sits at the
   2 m circle for CBM18 [R]. A wall can matter a lot: on the circular "circa" case an ideal wall makes
   n ≤ 12 stable (Lunniss 2016) [R].
8. **Edge current and q_a.** cbm18 has zero current at the edge and q_a = 3.03, so the peeling drive at
   the edge is small [R, I]. I found no cbm18-specific study of nq_a effects.
9. **Mapping and resolution.** About 10 % from small mapping differences in NIMROD (Burke 2010, via
   Ferraro) [R]; BOUT++ drifts high at large n at fixed radial resolution (Dudson 2011) [R]; JOREK needs
   resolution extrapolation for CBM18 at n = 20 (Hoelzl 2021) [R]; MARG2D and ELITE need 1025 × 1025
   near an X-point (Snyder 2009) [R].

---

## 6. MISHKA: the NaN on the extended cbm18a, and density options

- **Literature.** I found no published report of MISHKA failing on an equilibrium with an extended
  zero-pressure (pseudo-vacuum) region. The one published MISHKA run on CBM18 (Pamela 2020) does not say
  where the plasma boundary was placed [R].
- **Density options in MISHKA-1** (read in the local source `mishka1_fast.f`, November 2017 version) [L].
  s = √ψ_N, namelist variables `DSURF`, `DSURF1` to `DSURF5`, `IDPOW`, `IEQ`:
  - default (`IEQ` = 1): ρ = (1 − (1 − `DSURF`) s²)^`IDPOW`; `DSURF` = 1 gives a flat density;
  - `IEQ` = 2: ρ = (1 − `DSURF1` s²)(1 + `DSURF` s²) (JET case);
  - `IEQ` = 3: ρ = 1 + `DSURF` s² + `DSURF1` s⁴ + `DSURF2` s⁶ + `DSURF3` s⁸;
  - `IEQ` = 4: polynomial in s up to s⁶ with `DSURF` to `DSURF5`;
  - `IEQ` = 6: ρ = `DSURF1` + (1 − `DSURF1`)(p/p0)^`DSURF`, a power of the pressure plus a floor;
  - `IEQ` = 7: parabolic up to s = `DSURF`, then linear between `DSURF1` and `DSURF2`;
  - `IEQ` = 8: tanh step in ψ: ρ = `DSURF` + (1 − `DSURF`)(1 − tanh((ψ − `DSURF1`)/`DSURF2`))/2.
  With `IEQ` = 6 and `DSURF1` = 0 the density, and the inertia, vanish wherever p = 0, so this option
  cannot be used on cbm18a's pseudo-vacuum [I]. The profile table MISHKA prints contains
  `RLARM` = `CWW`·√(2p/ρ), which turns into NaN wherever the pressure spline dips below zero; it is only
  printed, so it does not by itself explain a NaN eigenvalue, but the printed table is the first place
  to look for where a NaN starts [L, I].
- **CASTOR.** The CASTOR paper (Kerner et al. 1998) was not read; its density options are not checked here.
- **Practical route [I].** Run HELENA for cbm18 with the boundary at the pedestal foot (ψ_N = 0.70 of the
  TOQ flux), MISHKA with a flat density (`DSURF` = 1, `IEQ` = 1) and a distant wall, and compare with PB3D
  free boundary cut at the same surface. This is the setup that matches ELITE's (cutoff 0.701, infinite
  vacuum) and avoids the pseudo-vacuum altogether.

---

## 7. Gaps I could not fill

1. The exact definition of the cbm18 series and of "dens": the unpublished ELITE user's guide (Snyder and
   Wilson 2006) is the probable source; ask the ELITE group for it with the TOQ files.
2. The R, B and ρ behind the "R/B = 1.5" label and behind Ferraro's γ/ω_A for cbm18.
3. The full text of Burke 2010 (paywalled; OSTI was unreachable): its figures, the second equilibrium,
   and whether it gives numbers.
4. Any public EQDSK or TOQ file for cbm18; only BOUT++ grid files covering ψ_N 0.4 to 1.2 are public.
5. The contents of the JOREK `CBM18.zip` growth-rate files (not downloaded; needs approval).
6. A recent (2015 to 2026) multi-code ideal edge-stability benchmark with tables.
7. Any DCON/GPEC peeling-ballooning result for cbm18.

## 8. Method note

All web requests were anonymous (WebSearch and WebFetch; no name, address or token sent). Two arXiv PDFs
(0810.5757 and 2011.09120) were fetched with curl and a generic "Mozilla/5.0" User-Agent instead of
WebFetch. OSTI, the EPFL EPS archive, the IAEA proceedings and fire.pppl.gov refused or blocked access,
so sources there are marked [S].
