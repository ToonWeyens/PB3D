# cbm18a cut at the pedestal foot

The committed cbm18a fixture (`tests/fixtures/cbm18a_extended_vac.12.xz`) extends the plasma beyond the
pedestal foot with a current-free, pressure-free region to ψ_N = 1 (q rises to 25). MISHKA-1 and
TERPSICHORE find spurious modes in that extension, so this directory holds the same plasma ending at the
foot, made on 1 November 2026 (issue #4):

- `cbm18a_cut.fort.10`: HELENA input (hel12a). The p' and FF' tables are the first 217 of the extended
  input's 259 points (HELENA's tables live on a uniform grid in sqrt(ψ_N), so the foot at index 217 is
  ψ_N = (216/258)² = 0.7009 of the extended flux); the boundary is the foot flux surface of the extended
  mapping as 24 cosine coefficients (`ISHAPE = 2`, `FM(1) = 2` is a unit circle); `EPS = 0.344262` from
  that surface; `B` and `XIAB` scaled from the extended values by (R0_new/R0_old)² and
  (a_old R0_new)/(a_new R0_old). No fitting was needed: q on axis and at the foot agree with the
  extended run to 6e−6, the q and p profiles inside to 3e−4 and 7e−6, the axis position to 3e−7.
  Run: copy to `fort.10`, `hel12a < /dev/null` (the plot prompt at the end is harmless), mapping in
  `fort.12`; NR 251, NP 201, NRMAP 401, NPMAP 513, NCHI 513 (power of two plus one, for MISHKA's FFT).
- `foot_surface.py <extended fort.12> <psi_c>`, `make_cut_input.py <extended fort.10> <out> B XIAB coarse|fine`,
  `parse_hel12.py`: the scripts that made it (PB3D's `read_HEL` order for the mapping).
- `input.cbm18a_cut`: the VMEC INDATA from PB3D's `--export_HEL` on the cut mapping (639 surfaces with
  VMEC++ reproduce q to five digits; `B_0 = R_0 = 1`).

Results, n = 10, 15, 20, ω²/ω_A² (PB3D units; TERPSICHORE's λ equals ω² for this twin):

| | PB3D, HELENA path | TERPSICHORE on the twin | MISHKA-1 |
|---|---|---|---|
| fixed boundary | −0.04012, −0.08128, −0.11120 | −0.03991, −0.07951, −0.10701 | 0.03637, 0.07009, 0.09290 (wall on the boundary) |
| free boundary (TERPSICHORE: wall at 2a) | −0.05002, −0.08577, −0.11209 | −0.04979, −0.08441, −0.10753 | 0.04492, 0.07407, 0.09453 |

PB3D and TERPSICHORE agree to 0.5 % at n = 10 and 3 to 4 % at n = 20 (PB3D is first order in its radial
grid on cbm18a, 1 % between 200 and 400 points at n = 20); MISHKA-1 is 9 to 16 % less unstable than both.
TERPSICHORE needs `LCURRF = 9` (parallel current from the geometry): its default reconstruction from
charge conservation loses the pedestal's edge current and gives one eighth of the eigenvalue at n = 10
(about a third of the growth rate).
