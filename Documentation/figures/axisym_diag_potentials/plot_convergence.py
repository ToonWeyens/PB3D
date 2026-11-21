#!/usr/bin/env python3
"""Solve error of the exterior problem for the exact toroidal harmonics against the number of
boundary points: 2.51, direct, option A, option B. Reads exact_harmonics.txt (lines of
study_diag: label shape diag n k n_bnd ierr interior exterior solve). Writes 5_convergence.png.
Usage: python3 plot_convergence.py <folder of this script>"""
import sys, os, collections
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
root = sys.argv[1]
D = collections.defaultdict(list)
for l in open(os.path.join(root, 'exact_harmonics.txt')):
    p = l.split()
    lab, shp, n, k, nb, ierr, sol = p[0], p[1], int(p[3]), int(p[4]), int(p[5]), int(p[6]), float(p[9])
    if ierr == 0 and sol > 0:
        D[(lab, shp, n, k)].append((nb, sol))
S = {'base': ('2.51 (Rⁿ, RⁿZ, Q, unit maximum)', '#2a78d6', 'o', '-'),
     'direct': ('direct diagonal', '#eb6834', 's', '-'),
     'A': ('A: P and Q, exact', '#1baf7a', 'D', '-'),
     'B': ('B: Rⁿ, RⁿZ, P, Q, least squares', '#eda100', '^', '--')}
panels = [('circle', 10, 10), ('D', 10, 10), ('circle', 30, 30), ('D', 30, 30)]
fig, axs = plt.subplots(2, 2, figsize=(9.5, 7.2), sharex=True)
for ax, (shp, n, k) in zip(axs.flat, panels):
    for lab, (name, col, mk, ls) in S.items():
        pts = sorted(D.get((lab, shp, n, k), []))
        if pts:
            ax.loglog([a for a, _ in pts], [b for _, b in pts], ls, color=col, marker=mk, ms=6, lw=2,
                      label=name, mfc=col if lab != 'B' else 'none', mew=1.5)
    for ref, slope, lab2 in [(2e-3, 1, 'slope −1'), (2e-4, 3, 'slope −3')]:
        x = [300, 1200]
        ax.loglog(x, [ref * (x[0] / xx) ** slope for xx in x], ':', color='#888888', lw=1)
        ax.text(1250, ref * (300 / 1200) ** slope, lab2, color='#666666', fontsize=8, va='center')
    ax.set_title('%s, n = %d, k = %d' % ('circle' if shp == 'circle' else 'D-shape (κ 1.7, δ 0.3)', n, k), fontsize=10)
    ax.grid(True, which='major', color='#e6e6e6', lw=0.6)
    ax.set_ylim(1e-7, 2e3)
    ax.set_xticks([151, 301, 601, 1201, 2401])
    ax.set_xticklabels(['151', '301', '601', '1201', '2401'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
for ax in axs[1]:
    ax.set_xlabel('boundary points')
for ax in axs[:, 0]:
    ax.set_ylabel('max error of the solved boundary values')
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='upper center', ncol=2, frameon=False, fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(os.path.join(root, '5_convergence.png'), dpi=150)
print('written')
