"""Area 2c: who fills which point in MPI_utilities.redistribute_var (commit 798fe50).

Every process hands in its whole local range, ghost points included, so the
source ranges overlap. The destination ranges do not overlap. The routine
decides, for every destination point, which process sends it.

2.47: for each destination, walk the sources in rank order and take from the
first one whose range still covers the next point. In an overlap the point
comes from the LOWER rank, i.e. from that rank's ghost region, where radial
derivatives and interpolations are one-sided and less accurate.

Now: split every overlap in the middle (the convention of trim_grid) and take
each point from its owner.

The two loops below are copied from the Fortran (old and new) with 1-based
inclusive limits. The ranges are an example with 3 processes, 60 points and
a ghost region of 4 points on each side.

Second panel: grid_ops now widens the range of the spline that interp_V fits
by 20 points (ghost_interp). A cubic not-a-knot spline fitted on a window
differs from the spline on the full range near the window ends; the
difference falls by about 2 - sqrt(3) = 0.27 per interval inward.
"""
import math

import numpy as np
from scipy.interpolate import CubicSpline

from _common import BLUE, GREY, ORANGE, plt, save

n_pts, n_procs, ghost = 60, 3, 4
own = [(1, 20), (21, 40), (41, 60)]
lims_tot = [(max(a - ghost, 1), min(b + ghost, n_pts)) for a, b in own]   # source ranges, overlapping
lims_dis_tot = own                                                         # destination ranges


def old_rule():
    src = np.zeros(n_pts + 1, int)
    for jd in range(n_procs):
        max_loc = lims_dis_tot[jd][0]
        for id_ in range(n_procs):
            if lims_tot[id_][1] >= max_loc:
                n = min(lims_tot[id_][1], lims_dis_tot[jd][1]) - max_loc + 1
                if n > 0:
                    src[max_loc:max_loc + n] = id_
                max_loc += max(n, 0)
    return src[1:]


def new_rule():
    eff = [list(r) for r in lims_tot]
    for id_ in range(1, n_procs):
        if lims_tot[id_][0] <= lims_tot[id_ - 1][1]:
            eff[id_][0] = lims_tot[id_][0] + math.ceil((lims_tot[id_ - 1][1] - lims_tot[id_][0] + 1.0) / 2)
            eff[id_ - 1][1] = eff[id_][0] - 1
    src = np.zeros(n_pts + 1, int)
    for jd in range(n_procs):
        for id_ in range(n_procs):
            start = max(lims_dis_tot[jd][0], eff[id_][0])
            n = max(0, min(lims_dis_tot[jd][1], eff[id_][1]) - start + 1)
            src[start:start + n] = id_
    return src[1:], eff


old = old_rule()
new, eff = new_rule()
in_ghost_old = [i + 1 for i in range(n_pts) if not (own[old[i]][0] <= i + 1 <= own[old[i]][1])]
print("source ranges incl. ghosts:", lims_tot)
print("owned ranges after the split:", [tuple(e) for e in eff])
print("points the 2.47 rule took from another process' ghost region:", in_ghost_old)
print("points the new rule takes from a non-owner:",
      [i + 1 for i in range(n_pts) if not (own[new[i]][0] <= i + 1 <= own[new[i]][1])])

fig, axes = plt.subplots(1, 2, figsize=(11, 3.4), gridspec_kw={"width_ratios": [1.4, 1]})
ax = axes[0]
cols = [BLUE, ORANGE, "#3a9a5b"]
for p in range(n_procs):
    a, b = lims_tot[p]
    ax.barh(3 - p * 0.5, b - a + 1, left=a - 0.5, height=0.35, color=cols[p], alpha=0.25)
    a, b = own[p]
    ax.barh(3 - p * 0.5, b - a + 1, left=a - 0.5, height=0.35, color=cols[p], alpha=0.8)
    ax.text(0, 3 - p * 0.5, f"rank {p}", va="center", ha="right", fontsize=8)
x = np.arange(1, n_pts + 1)
ax.scatter(x, np.full(n_pts, 1.2), c=[cols[s] for s in old], s=14, marker="s")
ax.scatter(x, np.full(n_pts, 0.7), c=[cols[s] for s in new], s=14, marker="s")
for i in in_ghost_old:
    ax.plot(i, 1.2, "kx", ms=6)
ax.text(0, 1.2, "source, 2.47", va="center", ha="right", fontsize=8)
ax.text(0, 0.7, "source, now", va="center", ha="right", fontsize=8)
ax.set_yticks([]); ax.set_xlim(-12, n_pts + 1); ax.set_xlabel("global radial index")
ax.set_title("pale: ghost points handed in; x: taken from a ghost region")
ax.grid(False)

# panel 2: local spline window vs full-range spline
xs = np.linspace(0, 1, 201)
f = np.sin(7 * xs) * np.exp(xs)
full = CubicSpline(xs, f)
i0, i1 = 80, 140
fine_idx = np.arange(i0, i1 + 1)
for margin, c in ((0, ORANGE), (20, BLUE)):
    a, b = i0 - margin, i1 + margin
    loc = CubicSpline(xs[a:b + 1], f[a:b + 1])
    xm = 0.5 * (xs[fine_idx[:-1]] + xs[fine_idx[1:]])               # midpoints, where the splines differ
    d = np.abs(loc(xm) - full(xm))
    lab = f"window widened by {margin}" + (" (round-off everywhere)" if margin else "")
    axes[1].semilogy(fine_idx[:-1] - i0, np.maximum(d, 2e-17), color=c, label=lab)
    if margin == 0:
        d0 = d[0]
axes[1].semilogy(np.arange(0, 30), d0 * (2 - math.sqrt(3)) ** np.arange(0, 30), color=GREY, ls=":", label="0.27 per interval")
axes[1].set_ylim(1e-17, 1e-1)
axes[1].set_xlabel("interval inside the window of one process")
axes[1].set_title("spline on a window vs on the full range")
axes[1].legend(fontsize=8)
save(fig, "redistribute_owner.png")
