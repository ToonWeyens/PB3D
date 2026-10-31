# Running TERPSICHORE on PB3D's VMEC equilibria

TERPSICHORE (EPFL, open source since 2025, Apache 2.0, https://gitlab.epfl.ch/spc/public/terpsichore)
is a global linear ideal-MHD stability code for 3-D equilibria. These two scripts feed it PB3D's
VMEC files (31 October 2026, issue #4; the notes of the first runs are in
`Documentation/plan_vacuum_bugfix.md`, evening of 31 October, and in `Documentation/benchmarking_prep.md`).

- `wout2terps.py wout.nc fort.18 [--flip]`: writes TERPSICHORE's equilibrium file from a netCDF
  wout (VMEC++ or VMEC 9), in the layout of its own `srcs/vmec2terps.f`, which only reads the
  VMEC 8 text format. Stellarator-symmetric files only. `--flip` reverses the sign of the toroidal
  flux (iota and phips), for files such as PB3D's HELENA twins that carry a negative flux; the
  shipped benchmark has both positive.
- `terps_data.py out.data title MM NMIN NMAX MMS NSMIN NSMAX NJ NK IVAC NOWALL AWALL MODELK NSTA AL0 [m_lo m_hi]`:
  writes the fixed-column control cards (not a namelist) for one toroidal mode number per row
  N = NSMIN..NSMAX (n = N·NFP/NSTA), poloidal modes m_lo..m_hi of the 56 columns (plus MPINIT,
  which the script leaves at 0), the Boozer table for N = NMIN..NMAX, fixed boundary with
  IVAC = 0 or a conformal wall at AWALL minor radii with IVAC > 0 and NOWALL = -1, the kinetic
  energy model MODELK (2 = physical) and the inverse-iteration shift AL0. The other cards are
  those of the shipped QAS3 benchmark. A conformal wall on an axisymmetric equilibrium needs the
  Boozer rows N = -1..1 (the wall routine indexes the n = NFP harmonic), and RPLMIN must then be
  negative so that the empty harmonics are kept.

Run (Snellius, `~/pb3d/ext/terpsichore`, built with `cd srcs && make`, 22 s, no changes):

    mkdir runs/<case> && cd runs/<case> && cp ../fort18_<eq> fort.18
    rm -f tpr16_dat_pvi tpr16_dat_wall; ../../srcs/tpr_ap.x < ../<case>.data > <case>.out

The lines "EIGENVALUE FROM WP / WK", "NUMBER OF NEGATIVE EIGENVALUES" (below the shift) and
"GROWTH RATE" (sqrt(|λ|/WALF), WALF printed as "AFVEN FREQUENCY SQUARED x MU_0 x MASS DENSITY")
carry the result; fort.22 holds the eigenvector. Each run converges to the eigenvalue nearest
the shift, so a shift scan with the count of eigenvalues below each shift brackets the spectrum.

## Rippled twins and W7-X (31 October 2026, late)

- `--nfp K` on `wout2terps.py` declares K field periods for a file written with nfp = 1 whose
  boundary has that periodicity (PB3D's rippled twins, a(1 + δ cos 5ζ) with NTOR = 6): only the
  modes with n a multiple of K are kept, and TERPSICHORE then solves one mode family (rows
  N = n/K in the stability table) instead of every n. Rippled twin, δ = 0.04, n = 10: fixed
  +0.078853 with the n = 10 row alone; free (wall at 2a) −0.412459 with the row alone and
  −0.419032 with the family n = 0, 5, 10, 15, 20, against PB3D's −0.421625 at 64 field lines and −0.425321 at 128
  (λ = ω² (B_0/R_0)², factor 1.006652 for this file).
- The Boozer table must hold at least as many modes as the VMEC input: `fvjac` is allocated
  with the Boozer count and filled with the VMEC count (`tprgl0_block_ap.f:666`), and the
  release build stops with a bare backtrace. W7-X with 288 modes needs n/5 = −12..12, m ≤ 14.
- A stable shift in a family run returns a continuum point, not a mode: the stable spectrum of
  a family is dense (the rippled twin's fixed-boundary family gave 0.07766, 0.07797 and 0.08601
  for the shifts 0.070, 0.078 and 0.086). Compare discrete (unstable) eigenvalues, or a single
  row when the continuum is empty near the shift.
- `srcs_fix/tpr_ap.x` is the corrected build (the two kinetic-energy lines swapped);
  `runs/run_terps_cases.sh case ...` runs a list of case directories (each with fort.18 and
  <case>.data) and prints the result lines.
