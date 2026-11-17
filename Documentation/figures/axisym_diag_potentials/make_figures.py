#!/usr/bin/env python3
"""Test potentials of the indirect diagonal of the axisymmetric vacuum
(diag_indirect_2 in vac_ops), now and as proposed (17 November 2026).

    python make_figures.py [--out DIR]

Now: R^|n| and R^|n| Z (harmonic inside the boundary, c = -4 pi) and the ring
source Q_{|n|-1/2}(gamma)/sqrt(R R_s) (harmonic outside and decaying, c = 0).
Proposed: the same Q and P_{|n|-1/2}(gamma)/sqrt(R R_s), the other solution of
the same Legendre equation (regular at the ring, harmonic inside, c = -4 pi),
with gamma = 1 + |x - x_s|^2/(2 R R_s) and the ring x_s = (R_s, Z_s) at the
centroid of the boundary points, as in diag_indirect_2.

Geometry: a circle R = 1 + a cos t, Z = a sin t and a D-shape
R = 1 + a cos(t + 0.3 sin t), Z = 1.7 a sin t, a = 0.33, 257 points, the normal
into the plasma (the direction of the stored normal -J grad psi). Normal data
of P and Q by a central difference with step 1e-6.

Figures (PNG, default: the directory of this script):
  1_potential_maps_n10.png   log10 |f| / max over the boundary of |f| in the
                             poloidal plane, n = 10;
  2_boundary_values.png      |f| / max over the boundary along the boundary,
                             n = 10 and 20;
  3_row_conditioning.png     per row i, the smallest singular value of the
                             matrix with rows (f_k, -f'_k a/|n|), one row per
                             potential k, each scaled to unit maximum as
                             diag_indirect_2 does: how well the identities of
                             row i fix H_ii and G_ii.

Plotted is the amplitude f(R, Z) of f e^{i n zeta}. Needs numpy, scipy and
matplotlib.
"""
import argparse
import os
import numpy as np
from scipy.special import hyp2f1, gamma as Gam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ap = argparse.ArgumentParser(description='Test potentials of the axisymmetric indirect diagonal.')
ap.add_argument('--out', default=os.path.dirname(os.path.abspath(__file__)), help='where the PNGs go')
out = ap.parse_args().out
os.makedirs(out, exist_ok=True)

# ---- style (reference palette, light mode) ----
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
C = {'Rn': '#2a78d6', 'RnZ': '#eb6834', 'Q': '#1baf7a', 'P': '#eda100'}
SEQ = LinearSegmentedColormap.from_list('blue', ['#fcfcfb', '#cde2fb', '#9ec5f4', '#6da7ec',
                                                 '#3987e5', '#256abf', '#184f95', '#0d366b'])
plt.rcParams.update({'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
                     'text.color': INK, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.edgecolor': GRID, 'font.size': 10, 'axes.titlesize': 10.5,
                     'axes.spines.top': False, 'axes.spines.right': False})

# ---- geometry: N boundary points, counterclockwise, normal into the plasma ----
A, N = 0.33, 257
def boundary(kind):
    th = 2*np.pi*np.arange(N)/N
    k, d = (1.0, 0.0) if kind == 'circle' else (1.7, 0.3)
    R = 1 + A*np.cos(th + d*np.sin(th)); Z = k*A*np.sin(th)
    dR = -A*np.sin(th + d*np.sin(th))*(1 + d*np.cos(th)); dZ = k*A*np.cos(th)
    t = np.hypot(dR, dZ)
    return th, R, Z, -dZ/t, dR/t

# ---- potentials ----
def gam(R, Z, Rs, Zs): return 1 + ((R-Rs)**2 + (Z-Zs)**2)/(2*R*Rs)
def P_ring(n, R, Z, Rs, Zs):                       # P_{|n|-1/2}(gamma)/sqrt(R Rs)
    nu = abs(n) - 0.5
    return hyp2f1(-nu, nu+1, 1, (1-gam(R, Z, Rs, Zs))/2)/np.sqrt(R*Rs)
def Q_ring(n, R, Z, Rs, Zs):                       # Q_{|n|-1/2}(gamma)/sqrt(R Rs)
    nu, g = abs(n) - 0.5, gam(R, Z, Rs, Zs)
    q = np.sqrt(np.pi)*Gam(nu+1)/(Gam(nu+1.5)*(2*g)**(nu+1))*hyp2f1((nu+1)/2, (nu+2)/2, nu+1.5, 1/g**2)
    return q/np.sqrt(R*Rs)
def potentials(n, R, Z, Rs, Zs):
    n = abs(n)
    return {'Rn': R**n, 'RnZ': R**n*Z, 'Q': Q_ring(n, R, Z, Rs, Zs), 'P': P_ring(n, R, Z, Rs, Zs)}
NAME = {'Rn': r'$R^{|n|}$', 'RnZ': r'$R^{|n|}Z$', 'Q': r'$Q_{|n|-1/2}(\gamma)/\sqrt{RR_s}$',
        'P': r'$P_{|n|-1/2}(\gamma)/\sqrt{RR_s}$'}
ROLE = {'Rn': 'now · interior, c = −4π', 'RnZ': 'now · interior, c = −4π',
        'Q': 'now and proposed · exterior, c = 0', 'P': 'proposed · interior, c = −4π'}

def bnd_data(kind, n):
    th, R, Z, nR, nZ = boundary(kind)
    Rs, Zs = R.mean(), Z.mean()                    # centroid of the boundary points, as diag_indirect_2
    e = 1e-6
    f0 = potentials(n, R, Z, Rs, Zs)
    fp = potentials(n, R+e*nR, Z+e*nZ, Rs, Zs); fm = potentials(n, R-e*nR, Z-e*nZ, Rs, Zs)
    d = {k: (f0[k], (fp[k]-fm[k])/(2*e)) for k in f0}
    return th, R, Z, Rs, Zs, d

# sanity: Q from the hypergeometric formula against Heine's integral
th, R, Z, nR, nZ = boundary('D')
phi = np.linspace(0, np.pi, 4001)
g = gam(R, Z, 1.0, 0.0)
heine = np.trapezoid(np.cos(10*phi)[None, :]/np.sqrt(g[:, None]-np.cos(phi)[None, :]), phi, axis=1)/np.sqrt(2)
print('Q check (max rel. diff, hyp2f1 vs Heine):',
      np.abs(Q_ring(10, R, Z, 1.0, 0.0)*np.sqrt(R) / heine - 1).max())

# ============ figure 1: maps in the poloidal plane, n = 10 ============
n = 10
Rg, Zg = np.meshgrid(np.linspace(0.45, 1.55, 441), np.linspace(-0.78, 0.78, 625))
fig, axs = plt.subplots(2, 4, figsize=(13, 7.6), constrained_layout=True)
for row, kind in enumerate(('circle', 'D')):
    th, R, Z, Rs, Zs, d = bnd_data(kind, n)
    F = potentials(n, Rg, Zg, Rs, Zs)
    for col, k in enumerate(('Rn', 'RnZ', 'P', 'Q')):
        ax = axs[row, col]
        L = np.log10(np.abs(F[k])/np.abs(d[k][0]).max())
        if k == 'Q':
            L = np.where(np.hypot(Rg-Rs, Zg-Zs) < 0.015, np.nan, L)
        im = ax.pcolormesh(Rg, Zg, np.clip(L, -4, 1), cmap=SEQ, vmin=-4, vmax=1, shading='auto', rasterized=True)
        ax.contour(Rg, Zg, L, levels=np.arange(-4, 1.01, 0.5), colors=INK2, linewidths=0.5, alpha=0.6)
        ax.contour(Rg, Zg, L, levels=[0], colors=INK2, linewidths=1.2)
        ax.plot(np.r_[R, R[0]], np.r_[Z, Z[0]], color=INK, lw=2)
        ax.plot(Rs, Zs, marker='+', color=INK, ms=9, mew=1.6)
        ax.set_aspect('equal'); ax.set_xlim(0.45, 1.55); ax.set_ylim(-0.78, 0.78)
        span = np.log10(np.abs(d[k][0]).max()/np.abs(d[k][0][np.abs(d[k][0]) > 0]).min())
        if k == 'RnZ':
            txt = 'zero on Z = 0'
        else:
            txt = f'max/min on boundary: {10**span:.3g}'
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', fontsize=9, color=INK,
                bbox=dict(boxstyle='round,pad=0.25', fc=SURF, ec=GRID))
        if row == 0:
            ax.set_title(f'{NAME[k]}\n{ROLE[k]}', color=INK)
        if col == 0:
            ax.set_ylabel(f'{"circle" if kind == "circle" else "D-shape"}\nZ / R₀')
        if row == 1:
            ax.set_xlabel('R / R₀')
cb = fig.colorbar(im, ax=axs, shrink=0.75, pad=0.01)
cb.set_label('log₁₀ |f| / max over the boundary of |f|')
fig.suptitle('Test potentials for n = 10, amplitude $f(R,Z)$ of $f\\,e^{in\\zeta}$, '
             'with $\\gamma = 1 + |\\mathbf{x}-\\mathbf{x}_s|^2/(2RR_s)$\n'
             'black: boundary ($a$ = 0.33 $R_0$)   +: ring source at the centroid of the boundary points   '
             'thick gray contour: the largest value on the boundary   thin contours: every half decade',
             fontsize=10, color=INK2)
fig.savefig(os.path.join(out, '1_potential_maps_n10.png'), dpi=150)
plt.close(fig)

# ============ figure 2: values on the boundary ============
LS = {'Rn': (0, (5, 2)), 'RnZ': (0, (5, 2)), 'Q': '-', 'P': '-'}
fig, axs = plt.subplots(2, 2, figsize=(11, 7.2), sharex=True, constrained_layout=True)
for row, kind in enumerate(('circle', 'D')):
    for col, n in enumerate((10, 20)):
        ax = axs[row, col]
        th, R, Z, Rs, Zs, d = bnd_data(kind, n)
        for k in ('Rn', 'RnZ', 'Q', 'P'):
            y = np.abs(d[k][0])/np.abs(d[k][0]).max()
            y = np.where(y < 1e-12, np.nan, y)
            ax.plot(th, y, color=C[k], ls=LS[k], lw=2, label=f'{NAME[k]}  ({ROLE[k].split(" ·")[0]})')
        ax.set_yscale('log'); ax.set_ylim(1e-7 if n == 20 else 1e-5, 8)
        ax.grid(True, which='major', color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi], ['0\noutboard', 'π/2\ntop', 'π\ninboard', '3π/2\nbottom', '2π'])
        ax.set_title(f'{"circle" if kind == "circle" else "D-shape"}, n = {n}', color=INK)
        if col == 0: ax.set_ylabel('|f| / max over the boundary')
        # direct labels: P above its maximum and Q below its minimum on the upper half,
        # R^|n|Z on the left flank (above R^|n|), R^|n| on the right flank (below R^|n|Z)
        half = slice(0, N//2 + 1)
        yv = {k: np.abs(d[k][0])/np.abs(d[k][0]).max() for k in d}
        i_lab = {'P': int(np.argmax(yv['P'][half])), 'Q': int(np.argmin(yv['Q'][half])),
                 'RnZ': int(0.40*N), 'Rn': int(0.70*N)}
        off = {'P': (0, 7), 'Q': (0, -15), 'RnZ': (14, 8), 'Rn': (14, -14)}
        for k in ('Rn', 'RnZ', 'P', 'Q'):
            i = i_lab[k]
            ax.annotate({'Rn': '$R^{|n|}$', 'RnZ': '$R^{|n|}Z$', 'Q': '$Q$', 'P': '$P$'}[k], (th[i], yv[k][i]),
                        textcoords='offset points', xytext=off[k], fontsize=10.5, color=INK, ha='center')
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='outside upper center', ncol=4, frameon=False, fontsize=9.5)
fig.supxlabel('boundary parameter θ (counterclockwise from the outboard midplane), 257 points', color=INK2, fontsize=10)
fig.savefig(os.path.join(out, '2_boundary_values.png'), dpi=150)
plt.close(fig)

# ============ figure 3: how well each row fixes H_ii and G_ii ============
SETS = [('now: $R^{|n|}$, $R^{|n|}Z$, $Q$', ('Rn', 'RnZ', 'Q'), C['Rn'], (0, (5, 2))),
        ('proposed: $P$, $Q$', ('P', 'Q'), C['RnZ'], '-'),
        ('all four', ('Rn', 'RnZ', 'P', 'Q'), C['Q'], (0, (1, 1.5)))]
fig, axs = plt.subplots(2, 2, figsize=(11, 7.2), sharex=True, constrained_layout=True)
for row, kind in enumerate(('circle', 'D')):
    for col, n in enumerate((10, 20)):
        ax = axs[row, col]
        th, R, Z, Rs, Zs, d = bnd_data(kind, n)
        for lab, ks, c, ls in SETS:
            P = [(d[k][0]/np.abs(d[k][0]).max(), d[k][1]/np.abs(d[k][0]).max()) for k in ks]
            sv = np.array([np.linalg.svd(np.array([[p[i], -dp[i]*A/n] for p, dp in P]), compute_uv=False)[-1]
                           for i in range(N)])
            sv = np.where(sv < 1e-14, np.nan, sv)
            ax.plot(th, sv, color=c, ls=ls, lw=2, label=lab)
            i = np.nanargmin(sv)
            ax.plot(th[i], sv[i], 'o', ms=6, color=c, mec=SURF, mew=1.5)
            ax.annotate(f'{np.nanmin(sv):.1e}', (th[i], sv[i]), textcoords='offset points', xytext=(8, -4),
                        fontsize=8.5, color=INK)
        ax.set_yscale('log'); ax.set_ylim(1e-9 if n == 20 else 1e-6, 2)
        ax.grid(True, which='major', color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi], ['0\noutboard', 'π/2\ntop', 'π\ninboard', '3π/2\nbottom', '2π'])
        ax.set_title(f'{"circle" if kind == "circle" else "D-shape"}, n = {n}', color=INK)
        if col == 0: ax.set_ylabel('smallest singular value of row i')
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='outside upper center', ncol=3, frameon=False, fontsize=10)
fig.supxlabel('boundary parameter θ of row i; matrix rows: (f_k, −f′_k a/|n|) per potential k, each scaled to unit '
              'maximum as diag_indirect_2 does; dots: worst row', color=INK2, fontsize=9.5)
fig.savefig(os.path.join(out, '3_row_conditioning.png'), dpi=150)
plt.close(fig)
print('written to', out)
