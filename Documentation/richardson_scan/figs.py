#!/usr/bin/env python3
"""figs.py -- figures and tables of Documentation/richardson_exponent.md.

Reads results.txt (one line per PB3D run: "EV <name> rc=.. t=..s <eigenvalue per Richardson level>") and
writes fig/*.png. With the argument "tables" it also prints the numbers that the report quotes.

Run names:
  <case>_s<rule>_n<N>          one run with N points along the field line; rule 1 = trapezoid,
                               2 = Simpson 3/8. Cases: V (axisymmetric VMEC twin), H257 ... H2049 (HELENA,
                               circular, with that many poloidal points), CBM (HELENA cbm18a),
                               R04_a<n> (rippled twin, n field lines), W7X_a<n>, <case>_t10 (one field
                               line over 10 poloidal turns).
  RL2_<case>_s<rule>           PB3D's own Richardson loop, 6 levels from 25 points.
  RAD_<case>_<edge>_p<p>_r<N>  N radial points, order p of the radial difference, edge fix or free.
  END_...                      as RAD, with the inner end of the radial range moved.
"""
import collections
import math
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"          # validated categorical slots 1 to 3 (light surface)
RULE = {1: ("trapezoid", C1, "o"), 2: ("Simpson 3/8", C2, "s")}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 9.5,
                     "axes.titleweight": "bold", "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "savefig.facecolor": SURF, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7, "legend.frameon": False})


# ---------------------------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------------------------
def read():
    res = {}
    for line in open(os.path.join(HERE, "results.txt")):
        p = line.split()
        if len(p) >= 5 and p[0] == "EV" and "NO_EV" not in line:
            res[p[1]] = [float(v) for v in p[4:]]
    return res


RES = read()
SER = collections.defaultdict(dict)                    # (case, rule) -> {N: eigenvalue}
for name, v in RES.items():
    m = re.match(r"(.+)_s([12])_n(\d+)$", name)
    if m and not name.startswith(("pilot", "RL")):
        case = re.sub(r"_a0*(\d+)$", lambda k: f"_a{int(k.group(1)):02d}", m.group(1))
        SER[(case, int(m.group(2)))][int(m.group(3))] = v[0]


def limit(case):
    """Best estimate of the eigenvalue for an infinite number of points along the line: the Simpson value at
    the largest N, or the trapezoid value there if it has changed less in the last doubling."""
    best = None
    for rule in (1, 2):
        d = SER.get((case, rule), {})
        chain = sorted(n for n in d if (n - 1) % 24 == 0 or (n - 1) % 240 == 0)
        if len(chain) < 2:
            continue
        last = abs(d[chain[-1]] - d[chain[-2]])
        if best is None or last < best[0]:
            best = (last, d[chain[-1]])
    return best[1], best[0]


def chain(case, rule, base=24):
    d = SER[(case, rule)]
    ns = sorted(n for n in d if (n - 1) % base == 0 and math.log2((n - 1) / base).is_integer())
    return ns, [d[n] for n in ns]


def guides(ax, x0, y0, powers, x1, ymin=None):
    """Thin lines proportional to N^-q through (x0, y0), with the label at their visible end."""
    for q in powers:
        xe = x1 if ymin is None else min(x1, x0 * (y0 / ymin) ** (1.0 / q))
        ye = y0 * (x0 / xe) ** q
        ax.plot([x0, xe], [y0, ye], color="#a9a8a2", lw=0.8, zorder=1)
        ax.annotate(f"$N^{{-{q}}}$", (xe, ye), xytext=(3, 3 if xe < x1 else 0), textcoords="offset points",
                    fontsize=7.5, color=INK2, va="bottom" if xe < x1 else "center")


def style_line(ax, x, y, color, marker, label):
    ax.plot(x, y, color=color, lw=1.6, marker=marker, ms=6.5, mec=SURF, mew=1.2, label=label, zorder=3,
            solid_capstyle="round", solid_joinstyle="round")


# ---------------------------------------------------------------------------------------------------------
# figure 1: error against the number of points along the line, VMEC cases
# ---------------------------------------------------------------------------------------------------------
def fig1():
    panels = [("V", "Axisymmetric VMEC twin, 1 field line"),
              ("R04_a01", "Rippled twin, 1 field line"),
              ("R04_a32", "Rippled twin, 32 field lines"),
              ("W7X_a01", "W7-X, 1 field line"),
              ("W7X_a08", "W7-X, 8 field lines"),
              ("W7X_a64", "W7-X, 64 field lines")]
    fig, axs = plt.subplots(2, 3, figsize=(10.5, 6.4), sharex=True, sharey=True)
    for ax, (case, title) in zip(axs.flat, panels):
        ref, _ = limit(case)
        ax.axhspan(1e-15, 1e-12, color="#efeeea", zorder=0)
        guides(ax, 25, 3e-3, (2, 4, 6), 6145, ymin=4e-14)
        for rule in (1, 2):
            ns, ev = chain(case, rule)
            pts = [(n, abs(e - ref) / abs(ref)) for n, e in zip(ns, ev) if e != ref]
            style_line(ax, [p[0] for p in pts], [p[1] for p in pts], RULE[rule][1], RULE[rule][2],
                       RULE[rule][0])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(title, loc="left")
        ax.set_ylim(1e-14, 3e-1)
        ax.set_xlim(20, 9000)
        ax.set_xticks([25, 49, 97, 193, 385, 769, 1537, 3073, 6145])
        ax.set_xticklabels(["25", "49", "97", "193", "385", "769", "1537", "3073", "6145"], rotation=45)
        ax.minorticks_off()
    axs[0, 1].text(22, 1.1e-13, "below the solver tolerance", fontsize=7.5, color=INK2, va="center")
    for ax in axs[1]:
        ax.set_xlabel("points along the field line, $N$")
    for ax in axs[:, 0]:
        ax.set_ylabel("relative error of the eigenvalue")
    axs[0, 0].legend(loc="upper right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig", "fig1_parallel_vmec.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------
# figure 2: HELENA
# ---------------------------------------------------------------------------------------------------------
def helena_levels():
    node = SER[("H257", 1)][513]                       # all points of the line are HELENA points
    interp = SER[("H257", 2)][3073]                    # integral of the linear interpolant
    return node, interp


def fig2():
    node, interp = helena_levels()
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.9))
    ax = axs[0]
    for rule in (1, 2):
        d = SER[("H257", rule)]
        ns = sorted(d)
        ax.plot(ns, [(d[n] - node) / node * 1e5 for n in ns], ls="none", marker=RULE[rule][2], ms=6.5,
                color=RULE[rule][1], mec=SURF, mew=1.2, label=RULE[rule][0], zorder=3)
    top = (interp - node) / node * 1e5
    for y in (0.0, top):
        ax.axhline(y, color="#a9a8a2", lw=0.8, zorder=1)
    ax.text(150, top + 0.22, "upper line: integral of the linear interpolant", fontsize=7.5, color=INK2)
    ax.text(700, 0.06, "lower line: every point of the\nfield line is a HELENA point", fontsize=7.5,
            color=INK2, va="bottom")
    ax.set_xscale("log")
    ax.set_ylim(-0.35, 1.75)
    ax.set_xlim(20, 5200)
    ax.set_xticks([25, 49, 97, 193, 385, 769, 1537, 4097])
    ax.set_xticklabels(["25", "49", "97", "193", "385", "769", "1537", "4097"])
    ax.minorticks_off()
    ax.set_title("Circular case, HELENA with 257 poloidal points", loc="left")
    ax.set_xlabel("points along the field line, $N$")
    ax.set_ylabel(r"relative difference to the lower line ($10^{-5}$)")
    ax.legend(loc="center right", fontsize=8.5)
    ax = axs[1]
    ref = SER[("CBM", 1)][1537]
    for rule in (1, 2):
        d = SER[("CBM", rule)]
        ns = sorted(n for n in d if n >= 61)
        ax.plot(ns, [(d[n] - ref) / abs(ref) * 1e5 for n in ns], ls="none", marker=RULE[rule][2], ms=6.5,
                color=RULE[rule][1], mec=SURF, mew=1.2, label=RULE[rule][0], zorder=3)
    ax.axhline(0.0, color="#a9a8a2", lw=0.8, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(52, 1900)
    ax.set_xticks([61, 97, 193, 385, 769, 1537])
    ax.set_xticklabels(["61", "97", "193", "385", "769", "1537"])
    ax.minorticks_off()
    ax.set_title("cbm18a, HELENA with 401 poloidal points", loc="left")
    ax.set_xlabel("points along the field line, $N$")
    ax.set_ylabel(r"relative difference to the value at $N = 1537$ ($10^{-5}$)")
    ax.legend(loc="lower right", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig", "fig3_helena.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------
# figure 3: what the extrapolation does
# ---------------------------------------------------------------------------------------------------------
def richardson(levels, q):
    """The table of calc_rich_ex: s(l,1) = level value, s(l,r) = s(l,r-1) + (s(l,r-1) - s(l-1,r-1)) /
    (2^(q (r-1)) - 1). Returns the best value of every level, s(l,l)."""
    tab = []
    for l, v in enumerate(levels):
        row = [v]
        for r in range(1, l + 1):
            row.append(row[r - 1] + (row[r - 1] - tab[l - 1][r - 1]) / (2 ** (q * r) - 1))
        tab.append(row)
    return [row[-1] for row in tab]


def w7x8_ref():
    d = SER[("W7X_a08", 2)]
    return d[3073] + (d[3073] - d[1537]) / 15.0        # Simpson at 3073 points, one step with exponent 4


def fig3():
    ref = w7x8_ref()
    unc = abs(SER[("W7X_a08", 2)][3073] - SER[("W7X_a08", 2)][1537]) / 15.0 / abs(ref)
    ns = [97, 193, 385, 769, 1537]
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.0), sharey=True)
    for ax, rule, qtrue in ((axs[0], 1, 2), (axs[1], 2, 4)):
        lev = [SER[("W7X_a08", rule)][n] for n in ns]
        curves = [("each level by itself", lev, "#8a8984", "o"),
                  ("extrapolated with exponent 6 (the code)", richardson(lev, 6), C2, "s"),
                  (f"extrapolated with exponent {qtrue} (measured)", richardson(lev, qtrue), C3, "^")]
        ax.axhspan(1e-14, 3 * unc, color="#efeeea", zorder=0)
        for lab, vals, col, mk in curves:
            style_line(ax, ns, [abs(v - ref) / abs(ref) for v in vals], col, mk, lab)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(ns)
        ax.set_xticklabels([str(n) for n in ns])
        ax.minorticks_off()
        ax.set_xlim(85, 1800)
        ax.set_ylim(1e-12, 1e-4)
        ax.set_title(f"W7-X, 8 field lines, {RULE[rule][0]}", loc="left")
        ax.set_xlabel("points along the field line at this Richardson level")
        ax.legend(loc="upper right", fontsize=8.5)
    axs[0].text(90, 2.2e-12, "within the uncertainty of the reference value", fontsize=7.5, color=INK2)
    axs[0].set_ylabel("relative error of the eigenvalue")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig", "fig4_extrapolation.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------
# figure 4: number of field lines
# ---------------------------------------------------------------------------------------------------------
def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.9), sharey=True)
    for ax, pre, title, nas, many in ((axs[0], "W7X", "W7-X", (1, 4, 8, 16, 32, 64), None),
                                      (axs[1], "R04", "Rippled twin", (1, 2, 4, 8, 16, 32), 7.7287903011863945e-02)):
        if many is None:
            many = limit(f"{pre}_a64")[0]
        par, lin, xs = [], [], []
        for na in nas:
            case = f"{pre}_a{na:02d}"
            lim, _ = limit(case)
            xs.append(na)
            e = abs(SER[(case, 1)][193] - lim)
            if e == 0.0:                                # the limit is this value: take the two rules apart
                e = abs(SER[(case, 1)][193] - SER[(case, 2)][193])
            par.append(e / abs(lim))
            lin.append(abs(lim - many) / abs(many))
        style_line(ax, xs, par, C1, "o", "left at 193 points along the line (trapezoid)")
        pts = [(x, y) for x, y in zip(xs, lin) if y > 0]
        style_line(ax, [p[0] for p in pts], [p[1] for p in pts], C2, "s", "from the number of field lines")
        ax.axhspan(1e-15, 1e-12, color="#efeeea", zorder=0)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xticks(list(nas))
        ax.set_xticklabels([str(n) for n in nas])
        ax.set_ylim(1e-14, 1)
        ax.set_title(title, loc="left")
        ax.set_xlabel("number of field lines")
        ax.legend(loc="upper right", fontsize=8.5)
    axs[0].set_ylabel("relative error of the eigenvalue")
    axs[0].text(1.0, 1.1e-13, "below the solver tolerance", fontsize=7.5, color=INK2, va="center")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig", "fig2_field_lines.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------
# figure 5: radial points and the order of the difference
# ---------------------------------------------------------------------------------------------------------
RAD = collections.defaultdict(dict)                    # (case, edge, p) -> {N: eigenvalue}
for name, v in RES.items():
    m = re.match(r"RAD_(.+)_(fix|free)_p(\d)_r(\d+)$", name)
    if m:
        RAD[(m.group(1), m.group(2), int(m.group(3)))][int(m.group(4))] = v[0]


def fig5():
    panels = [(("H257", "free"), "Circular case, free edge"), (("H257", "fix"), "Circular case, fixed edge"),
              (("CBM", "fix"), "cbm18a, fixed edge")]
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.9), sharey=True)
    for ax, ((case, edge), title) in zip(axs, panels):
        scale = abs(RAD[(case, edge, 3)][3200])
        guides(ax, 100, 3e-2, (1, 2, 3), 1600)
        for p, col, mk in ((1, C1, "o"), (2, C2, "s"), (3, C3, "^")):
            d = RAD[(case, edge, p)]
            ns = sorted(d)
            style_line(ax, ns[:-1], [abs(d[ns[i + 1]] - d[ns[i]]) / scale for i in range(len(ns) - 1)], col, mk,
                       f"order {p}")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks([100, 200, 400, 800, 1600])
        ax.set_xticklabels(["100", "200", "400", "800", "1600"])
        ax.minorticks_off()
        ax.set_xlim(85, 2300)
        ax.set_ylim(1e-8, 3e-1)
        ax.set_title(title, loc="left")
        ax.set_xlabel("radial points $N$ (change from $N$ to $2N$)")
    axs[0].set_ylabel("|change of the eigenvalue| / |eigenvalue|")
    axs[0].legend(loc="lower left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig", "fig5_radial.png"), dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------------------------------------
def tables():
    print("## limits and exponents (doubling chain N-1 = 24 * 2^k)")
    for case in ("V", "R04_a01", "R04_a02", "R04_a04", "R04_a08", "R04_a16", "R04_a32", "W7X_a01", "W7X_a04",
                 "W7X_a08", "W7X_a16", "W7X_a32", "W7X_a64"):
        ref, unc = limit(case)
        print(f"\n{case}: limit {ref:.13e} (last change {unc / abs(ref):.1e} relative)")
        for rule in (1, 2):
            ns, ev = chain(case, rule)
            err = [(e - ref) / abs(ref) for e in ev]
            ch = [ev[i + 1] - ev[i] for i in range(len(ev) - 1)]
            expo = ["%.2f" % math.log2(abs(ch[i] / ch[i + 1])) if ch[i + 1] != 0 else "-" for i in range(len(ch) - 1)]
            print(f"  rule {rule}: N = {ns}")
            print("     errors    " + " ".join(f"{e:+.1e}" for e in err))
            print("     exponents " + " ".join(expo))
    for case in ("R04_t10", "W7X_t10", "V_t10"):
        for rule in (1, 2):
            ns, ev = chain(case, rule, base=240)
            ch = [ev[i + 1] - ev[i] for i in range(len(ev) - 1)]
            expo = ["%.2f" % math.log2(abs(ch[i] / ch[i + 1])) if ch[i + 1] != 0 else "-" for i in range(len(ch) - 1)]
            print(f"{case} rule {rule}: N = {ns}  changes " + " ".join(f"{c / abs(ev[-1]):+.1e}" for c in ch) +
                  "  exponents " + " ".join(expo))
    node, interp = helena_levels()
    print(f"\n## HELENA 257: node-aligned {node:.10e}, interpolant {interp:.10e}, difference {(interp - node) / node:.3e}")
    for n, k in ((1025, 2), (1537, 3), (2049, 4), (3073, 6), (4097, 8)):
        pred = node + (1 - 1 / k ** 2) * (interp - node)
        print(f"   N = {n}: k = {k}, measured {SER[('H257', 1)][n]:.10e}, predicted {pred:.10e}")
    for case in ("H257", "H513", "H1025", "H2049", "CBM"):
        for rule in (1, 2):
            d = SER[(case, rule)]
            vals = [d[n] for n in sorted(d) if 49 <= n <= 481]
            print(f"   {case} rule {rule}: spread over N = 49..481: {(max(vals) - min(vals)) / abs(vals[0]):.2e}")
    print("\n## the extrapolation on PB3D's own levels (25 to 769 points)")
    for run, case, refv in (("RL2_W7X8", "W7X_a08", w7x8_ref()), ("RL2_R04", "R04_a04", None)):
        if refv is None:
            d = SER[(case, 2)]
            refv = d[1537] + (d[1537] - d[769]) / 15.0
        for rule, qtrue in ((1, 2), (2, 4)):
            lev = RES[f"{run}_s{rule}"]
            single = [SER[(case, rule)][n] for n in (25, 49, 97, 193, 385, 769)]
            print(f"{run} rule {rule}: largest difference level against single run "
                  f"{max(abs(a - b) / abs(b) for a, b in zip(lev, single)):.1e}")
            for lab, vals in (("levels", lev), ("q=6", richardson(lev, 6)), (f"q={qtrue}", richardson(lev, qtrue))):
                print(f"   {lab:7s} " + " ".join(f"{(v - refv) / abs(refv):+.1e}" for v in vals))
    print("\n## the extrapolation from 97 points (single runs as levels), W7-X with 8 lines")
    refv = w7x8_ref()
    for rule, qtrue in ((1, 2), (2, 4)):
        lev = [SER[("W7X_a08", rule)][n] for n in (97, 193, 385, 769, 1537)]
        for lab, vals in (("levels", lev), ("q=6", richardson(lev, 6)), (f"q={qtrue}", richardson(lev, qtrue))):
            print(f"   rule {rule} {lab:7s} " + " ".join(f"{(v - refv) / abs(refv):+.1e}" for v in vals))
    print("\n## field lines")
    many = {"W7X": limit("W7X_a64")[0], "R04": 7.7287903011863945e-02}
    for pre, nas in (("W7X", (1, 4, 8, 16, 32, 64)), ("R04", (1, 2, 4, 8, 16, 32))):
        for na in nas:
            case = f"{pre}_a{na:02d}"
            lim, _ = limit(case)
            print(f"   {case}: limit {lim:.11e}  parallel error at 97: {abs(SER[(case, 1)][97] - lim) / abs(lim):.1e}"
                  f"  at 193: {abs(SER[(case, 1)][193] - lim) / abs(lim):.1e}  line error {abs(lim - many[pre]) / abs(many[pre]):.1e}")
    for na in (5, 10, 12, 20, 24, 40, 48, 64):
        v = RES[f"R04_a{na:02d}_s1_n0097"][0]
        print(f"   R04 {na} lines at N = 97: {v:.11e}  ({(v - many['R04']) / many['R04']:+.1e})")
    print("\n## radial")
    for key in sorted(RAD):
        d = RAD[key]
        ns = sorted(d)
        ch = [d[ns[i + 1]] - d[ns[i]] for i in range(len(ns) - 1)]
        expo = ["%.2f" % math.log2(abs(ch[i] / ch[i + 1])) for i in range(len(ch) - 1)]
        print(f"   {key}: " + " ".join(f"{d[n]:+.7e}" for n in ns) + "  exponents " + " ".join(expo))
    print("\n## inner end moved")
    for name in sorted(n for n in RES if n.startswith("END_")):
        print(f"   {name} {RES[name][0]:.12e}")


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "fig"), exist_ok=True)
    fig1(); fig2(); fig3(); fig4(); fig5()
    if len(sys.argv) > 1 and sys.argv[1] == "tables":
        tables()
