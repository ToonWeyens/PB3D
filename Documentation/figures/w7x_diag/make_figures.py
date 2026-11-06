"""Figures of the W7-X check_vac runs with vac_diag_style 1 and 2 (issue #4, 6 November 2026).
Usage: w7x_diag_figs.py <pts-dir> <rows.json> <out-dir>
  maps:        w7x_diag_maps.png  (style 1 at 128 x 121: the diagonal difference, the cell integral,
               the axis-ring solve error; and the axis-ring solve error of style 2), in the index plane
               and in the (toroidal angle within a period, poloidal position) plane
  convergence: w7x_diag_convergence.png (max and median of the diagonal difference, identity residual,
               exterior solve errors, eigenvalue, against the number of boundary points)"""
import sys, os, json, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

pts_dir, rows_json, out_dir = sys.argv[1:4]
os.makedirs(out_dir, exist_ok=True)
NFP = 5

def load(run_glob, fn):
    d = sorted(glob.glob(os.path.join(pts_dir, run_glob)))
    if not d:
        return None
    return np.loadtxt(os.path.join(d[-1], fn))

def panel_index(ax, a, title, mask_zero=False, vmin=None, vmax=None):
    n_par = int(a[:, 1].max()); n_line = int(a[:, 2].max())
    v = a[:, 6].reshape(n_line, n_par)        # id = i_par + (i_line-1)*n_par
    v = np.ma.masked_where(v <= 0, v) if mask_zero else np.ma.masked_where(v <= 0, v)
    im = ax.imshow(v, origin='lower', aspect='auto', norm=LogNorm(vmin=vmin, vmax=vmax),
                   extent=[1, n_par, 1, n_line], cmap='viridis')
    ax.set_xlabel('point along the line, i_par'); ax.set_ylabel('field line, i_line'); ax.set_title(title, fontsize=9)
    return im

def panel_zeta(ax, a, title, vmin=None, vmax=None):
    """mean of the value in bins of the (toroidal angle within a period, position along the line) plane"""
    n_par = int(a[:, 1].max())
    zeta = np.degrees(np.arctan2(a[:, 4], a[:, 3])) % (360.0 / NFP)
    th = (a[:, 1] - 1) / (n_par - 1)
    v = a[:, 6]; sel = v > 0
    nb = (72, 60)
    s, xe, ye = np.histogram2d(zeta[sel], th[sel], bins=nb, range=[[0, 360.0 / NFP], [0, 1]], weights=v[sel])
    c, _, _ = np.histogram2d(zeta[sel], th[sel], bins=nb, range=[[0, 360.0 / NFP], [0, 1]])
    m = np.ma.masked_where(c == 0, s / np.maximum(c, 1))
    im = ax.imshow(m.T, origin='lower', aspect='auto', extent=[0, 360.0 / NFP, 0, 1],
                   norm=LogNorm(vmin=vmin, vmax=vmax), cmap='viridis')
    ax.set_xlabel('toroidal angle in a period (deg)')
    ax.set_ylabel('position along the line (0 to 1)'); ax.set_title(title, fontsize=9)
    return im

# --- maps at 128 x 121
d1 = 'd8s1_w7x_na128_p121-*'; d2 = 'd8s2_w7x_na128_p121-*'
items = [(d1, 'check_vac_diag_diff.dat', 'style 1: |row sum - direct| / 4pi'),
         (d1, 'check_vac_diag_Ising.dat', 'cell integral |I_sing| / 4pi'),
         (d1, 'check_vac_solve_axis_ring.dat', 'style 1: axis-ring solve error / max|phi|'),
         (d2, 'check_vac_solve_axis_ring.dat', 'style 2 (direct): axis-ring solve error / max|phi|')]
fig, axs = plt.subplots(2, 4, figsize=(18, 8))
for k, (rg, fn, title) in enumerate(items):
    a = load(rg, fn)
    if a is None:
        continue
    vmin, vmax = (1e-4, 1e-1) if 'solve' in fn else (1e-4, 1e-1)
    im = panel_index(axs[0, k], a, title, vmin=vmin, vmax=vmax)
    plt.colorbar(im, ax=axs[0, k], fraction=0.046)
    sc = panel_zeta(axs[1, k], a, title, vmin=vmin, vmax=vmax)
    plt.colorbar(sc, ax=axs[1, k], fraction=0.046)
fig.suptitle('W7-X (wout_w7xflip), 128 field lines x 121 points, check_vac per-point files: top row in the index plane, bottom row binned means over (toroidal angle in a period: bean at 0 and 72 deg, triangle at 36; poloidal position along the line, 0.5 = inboard); line ends (no full cell) left out', fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'w7x_diag_maps.png'), dpi=85)
print('wrote maps')

# --- convergence
rows = json.load(open(rows_json))
rows = [r for r in rows if not (r['np'] == 'np1' and any(x['run'] == r['run'] and x['np'] == 'np4' for x in rows))]
series = {'d8s1': ('style 1, indirect (default)', 'C0', 'o'), 'd8s2': ('style 2, direct', 'C3', 's'),
          'cv7': ('1 November series (older binary)', 'C7', '^')}
fig, axs = plt.subplots(2, 2, figsize=(13, 9))
def pick(key, f):
    out = {}
    for r in rows:
        if r['series'] != key:
            continue
        try:
            out[r['n_bnd']] = (f(r), '%dx%d' % (r['n_alpha'], r['n_par']))
        except (KeyError, TypeError, IndexError):
            pass
    ks = sorted(out)
    return np.array(ks), np.array([out[k][0] for k in ks]), [out[k][1] for k in ks]
def plot_series(ax, key, f, label_suffix='', ls='-', annotate=False):
    x, y, lab = pick(key, f)
    if len(x) == 0:
        return
    name, col, mk = series[key]
    ax.loglog(x, y, ls, color=col, marker=mk, ms=5, label=name + label_suffix)
    if annotate:
        for xi, yi, li in zip(x, y, lab):
            ax.annotate(li, (xi, yi), fontsize=6, textcoords='offset points', xytext=(3, 3))
ax = axs[0, 0]
plot_series(ax, 'd8s1', lambda r: r['pt_diag']['max'], ' max', annotate=True)
plot_series(ax, 'd8s1', lambda r: r['pt_diag']['med'], ' median', '--')
plot_series(ax, 'd8s1', lambda r: r['pt_Ising']['max'], ', max |I_sing|/4pi', ':')
ax.set_title('|row-sum diagonal - direct diagonal| / 4pi (style 1)', fontsize=9); ax.legend(fontsize=7)
ax = axs[0, 1]
for key in series:
    plot_series(ax, key, lambda r: r['ident']['x'][0], ', x', '-')
    plot_series(ax, key, lambda r: r['ident']['z'][0], ', z', '--')
ax.set_title('identity residual, interior harmonics x and z (c = -4pi), max norm', fontsize=9); ax.legend(fontsize=7)
ax = axs[1, 0]
for key in series:
    plot_series(ax, key, lambda r: r['solve']['axis_ring'][0], ', axis ring', '-')
    plot_series(ax, key, lambda r: r['solve']['axis_ring_n'][0], ', ring cos(n zeta)', '--')
ax.set_title('exterior solve error (max norm, relative to max|phi|)', fontsize=9); ax.legend(fontsize=7)
ax = axs[1, 1]
for key in series:
    x, y, lab = pick(key, lambda r: r['EV'])
    if len(x):
        ax.semilogx(x, y, '-', color=series[key][1], marker=series[key][2], ms=5, label=series[key][0])
ax.set_title('free-boundary eigenvalue, n = 10, BC_style 1 2, window 0.85-1', fontsize=9); ax.legend(fontsize=7)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=False)
for ax in axs.flat:
    ax.set_xlabel('boundary points n_bnd'); ax.grid(True, which='both', alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'w7x_diag_convergence.png'), dpi=100)
print('wrote convergence')
