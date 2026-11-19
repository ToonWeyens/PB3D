#!/usr/bin/env python3
"""Generate PB3D decks for the scans over the number of points along the field line.
usage: gen.py <family> ; writes decks/<name> and prints list lines "<name> <deck> <eq>"."""
import sys, os
ROOT = os.environ["PB3D_ROOT"]              # the PB3D folder on the cluster, e.g. /home/<user>/pb3d
FX = f"{ROOT}/fixtures"
RD = f"{ROOT}/rich/decks"
BASE = dict(min_n_par_X=100, min_par_X="0.0", max_par_X="2.0", alpha="0.0", min_r_sol="0.2", max_r_sol="1.0",
            n_r_sol=200, max_it_zero=500, use_pol_flux_F=".true.", use_normalization=".true.", prim_X=10,
            n_mod_X=10, tol_rich="1.0E-14", tol_SLEPC="1.0E-12", n_sol_requested=1, BC_style="1 1", K_style=1,
            norm_style=1, U_style=3, max_it_slepc=10000, max_tot_mem=3000, plot_magn_grid=".false.",
            plot_B=".false.", plot_J=".false.", plot_kappa=".false.", plot_flux_q=".false.",
            plot_resonance=".false.", max_it_rich=1, rich_restart_lvl=1, n_zeta_plot=1, magn_int_style=1,
            ex_plot_style=2)
EQ = dict(H257=f"{FX}/circ_eps033_q104.12", H513=f"{FX}/circ_eps033_nchi513.12",
          H1025=f"{FX}/circ_eps033_nchi1025.12", H2049=f"{FX}/circ_eps033_nchi2049.12",
          V=f"{FX}/wout_circ_eps033_q104.nc", R04=f"{FX}/wout_circ_eps033_rip04.nc",
          W7X=f"{FX}/wout_w7xflip.nc", CBM=f"{FX}/cbm18a_extended_vac.12")
CASE = dict(W7X=dict(min_r_sol="0.85", n_r_sol=60, n_mod_X=8),
            CBM=dict(min_r_sol="0.45", max_r_sol="0.75", n_mod_X=20))
out = []
def deck(name, eq, **kw):
    d = dict(BASE); d.update(CASE.get(eq, {})); d.update(kw)
    os.makedirs("decks", exist_ok=True)
    with open(f"decks/{name}", "w") as f:
        f.write("&inputdata_PB3D\n")
        for k, v in d.items():
            f.write(f"    {k:<24}= {v}\n")
        f.write("/\n")
    out.append(f"{name} {RD}/{name} {EQ[eq]}")
P2 = [25, 49, 97, 193, 385, 769, 1537]
P3 = [31, 61, 121, 241, 481, 961]
fam = sys.argv[1]
if fam == "pilot":
    for tol in ("1.0E-6", "1.0E-9", "1.0E-12", "1.0E-14"):
        for eq in ("H257", "V"):
            deck(f"pilot_{eq}_s2_n97_tol{tol[-2:].strip('-')}", eq, min_n_par_X=97, magn_int_style=2, tol_SLEPC=tol)
elif fam == "axi":        # axisymmetric: HELENA at four poloidal resolutions, and the VMEC twin
    for eq in ("H257", "H513", "H1025", "H2049", "V"):
        for s in (1, 2):
            for n in sorted(P2 + P3):
                deck(f"{eq}_s{s}_n{n:04d}", eq, min_n_par_X=n, magn_int_style=s)
elif fam == "cbm":
    for s in (1, 2):
        for n in sorted(P2 + P3):
            deck(f"CBM_s{s}_n{n:04d}", "CBM", min_n_par_X=n, magn_int_style=s)
elif fam == "rip":        # rippled twin, several field lines of one poloidal turn each
    for na, plist in ((1, P2), (2, P2), (4, P2), (8, P2[:6]), (16, P2[:6]), (32, P2[:5])):
        for s in (1, 2):
            for n in plist:
                deck(f"R04_a{na:02d}_s{s}_n{n:04d}", "R04", min_n_par_X=n, magn_int_style=s, alpha_style=2, n_alpha=na)
elif fam == "w7x":
    for na, plist in ((1, P2), (4, P2), (8, P2), (16, P2[:6]), (32, P2[:6]), (64, P2[:5])):
        for s in (1, 2):
            for n in plist:
                deck(f"W7X_a{na:02d}_s{s}_n{n:04d}", "W7X", min_n_par_X=n, magn_int_style=s, alpha_style=2, n_alpha=na)
elif fam == "turns":      # one field line followed over 10 poloidal turns
    for eq in ("R04", "W7X", "V"):
        for s in (1, 2):
            for per in (24, 48, 96, 192, 384):
                n = 10 * per + 1
                deck(f"{eq}_t10_s{s}_n{n:04d}", eq, min_n_par_X=n, magn_int_style=s, alpha_style=1,
                     min_par_X="-10.0", max_par_X="10.0")
elif fam == "richlvl":    # the code's own Richardson loop, 6 levels from 25 points
    for eq, kw in (("H257", {}), ("V", {}), ("CBM", {}), ("R04", dict(alpha_style=2, n_alpha=4)),
                   ("W7X", dict(alpha_style=2, n_alpha=16))):
        for s in (1, 2):
            deck(f"RL_{eq}_s{s}", eq, min_n_par_X=25, magn_int_style=s, max_it_rich=6, **kw)
elif fam == "radial":     # second question: radial points and the order of the difference
    for eq, kw0 in (("H257", {}), ("V", {}), ("CBM", {})):
        for bc, bcs in (("fix", "1 1"), ("free", "1 2")):
            if eq != "H257" and bc == "free":
                continue
            for p in (1, 2, 3):
                for nr in (100, 200, 400, 800, 1600, 3200):
                    deck(f"RAD_{eq}_{bc}_p{p}_r{nr:04d}", eq, min_n_par_X=100, magn_int_style=2, n_r_sol=nr,
                         norm_disc_prec_sol=p, BC_style=bcs, max_tot_mem=6000, **kw0)
elif fam == "extra":
    for eq, kw in (("H257", {}), ("V", {}), ("CBM", {}), ("R04", dict(alpha_style=2, n_alpha=4)),
                   ("W7X", dict(alpha_style=2, n_alpha=16)), ("W7X8", dict(alpha_style=2, n_alpha=8))):
        for s in (1, 2):
            deck(f"RL2_{eq}_s{s}", eq.replace("W7X8", "W7X"), min_n_par_X=25, magn_int_style=s, max_it_rich=6,
                 tol_SLEPC="6*1.0E-12", **kw)
    for na, plist in ((1, (3073, 6145)), (4, (3073,)), (8, (3073,))):
        for s in (1, 2):
            for n in plist:
                deck(f"W7X_a{na:02d}_s{s}_n{n:04d}", "W7X", min_n_par_X=n, magn_int_style=s, alpha_style=2, n_alpha=na)
    for na in (5, 10, 12, 20, 24, 40, 48, 64):
        deck(f"R04_a{na:02d}_s1_n0097", "R04", min_n_par_X=97, magn_int_style=1, alpha_style=2, n_alpha=na)
    for na in (64, 96, 128):
        deck(f"R04_a{na:03d}_s1_n0049", "R04", min_n_par_X=49, magn_int_style=1, alpha_style=2, n_alpha=na)
    for n in (33, 65, 129, 257, 513, 1025, 2049, 4097):
        deck(f"H257_s1_n{n:04d}", "H257", min_n_par_X=n, magn_int_style=1)
    deck("H257_s1_n3073", "H257", min_n_par_X=3073, magn_int_style=1)
    deck("H257_s2_n3073", "H257", min_n_par_X=3073, magn_int_style=2)
elif fam == "ends":      # which end gives the 1/N at a fixed edge: move the inner end, same step
    # same radial step 0.8/(N-1) as the reference scan: range 0.1-1.0 with 9/8 as many intervals, 0.4-1.0 with 6/8
    for p in (2, 3):
        for rmin, fac in (("0.1", 9), ("0.4", 6)):
            for nr0 in (400, 800, 1600):
                nr = (nr0 * fac) // 8 + 1
                deck(f"END_H257_fix_p{p}_in{rmin}_r{nr0:04d}", "H257", min_n_par_X=100, magn_int_style=2,
                     n_r_sol=nr, min_r_sol=rmin, norm_disc_prec_sol=p, BC_style="1 1", max_tot_mem=6000)
        for nr0 in (400, 800, 1600):
            deck(f"END_H257_fix_p{p}_in0.2_r{nr0:04d}", "H257", min_n_par_X=100, magn_int_style=2,
                 n_r_sol=nr0 + 1, norm_disc_prec_sol=p, BC_style="1 1", max_tot_mem=6000)
print("\n".join(out))
