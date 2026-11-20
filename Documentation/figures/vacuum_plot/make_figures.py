#!/usr/bin/env python3
"""Figures of the vacuum potential plot of POST (issue #4, 6 November 2026).

    python make_figures.py RUN_HEL RUN_VMEC ZOOM_HEL ZOOM_VMEC [--id 1] [--out DIR]

RUN_HEL and RUN_VMEC are POST run directories (with Data/) of the circular
eps = 0.33 HELENA case (axisymmetric vacuum, vac%style 2) and of its VMEC twin
(field-line vacuum, vac%style 1), run with plot_vac_pot and plot_sol_Q at one
plane, zeta_F = 0; ZOOM_HEL and ZOOM_VMEC the same with a narrow box across
the outboard midplane (min_Rvac_plot 1.25, max_Rvac_plot 1.45, min_Zvac_plot
-0.01, max_Zvac_plot 0.01, n_vac_plot 301 3). --id is the eigenmode, --out
where the PNG goes (default: the directory of this script).

Per vacuum style: a map with the normal component Q_n of the magnetic
perturbation inside the plasma and the vacuum potential phi outside, each
divided by its largest value on the boundary; and the potential along the
outboard midplane against the distance from the boundary, starting from the
boundary potential that the exterior system gives (marker) with the slope that
its normal data give (dashed line), dphi = norm.grad phi = -J Q_n with norm =
-J grad psi. The eigenvector's complex normalisation is arbitrary: all
quantities of a run are multiplied by the same unit factor that makes the
boundary potential at the outboard midplane real and positive.

Needs numpy, matplotlib and h5py.
"""
import argparse, os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.path import Path
from matplotlib.patches import PathPatch

INK, INK2, MUTED, GRID, BASE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
BLUE, ORANGE = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': BASE, 'axes.labelcolor': INK2, 'xtick.color': INK2,
    'ytick.color': INK2, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
    'axes.spines.right': False, 'figure.facecolor': 'white', 'axes.facecolor': 'white', 'legend.frameon': False,
    'text.color': INK, 'axes.titlecolor': INK, 'lines.linewidth': 1.6})
DIVERGING = LinearSegmentedColormap.from_list('div', ['#1c5cab', '#2a78d6', '#86b6ef', '#f0efec', '#f3a07f',
                                                      '#eb6834', '#b8431a'])


def read2(run, name):
    """R, Z and the complex field of a POST plot at one plane (RE and IM files)."""
    out = []
    for part in ['RE', 'IM']:
        with h5py.File(os.path.join(run, 'Data', name + '_' + part + '.h5'), 'r') as f:
            out.append((f['R_1'][:].squeeze(), f['Z_1'][:].squeeze(), f['var_1'][:].squeeze()))
    return out[0][0], out[0][1], out[0][2] + 1j*out[1][2]


def load(run, mode):
    """Plasma Q_n, vacuum potential, plasma indicator and boundary data of a run."""
    Rq, Zq, Qn = read2(run, '%d_sol_Qn' % mode)
    Rv, Zv, phi = read2(run, '%d_vac_pot' % mode)
    with h5py.File(os.path.join(run, 'Data', 'vac_in_plasma.h5'), 'r') as f:
        inpl = f['var_1'][:].squeeze()
    b = np.loadtxt(os.path.join(run, 'Data', '%d_vac_bnd.dat' % mode))
    dim = (b.shape[1] - 6)//2
    x, nrm = b[:, 2:2+dim], b[:, 2+dim:2+2*dim]
    dphi = b[:, -4] + 1j*b[:, -3]
    Phi = b[:, -2] + 1j*b[:, -1]
    R_b = np.hypot(x[:, 0], x[:, 1]) if dim == 3 else x[:, 0]
    Z_b = x[:, -1]
    # length unit of the plots against the normalised boundary points
    L = Rq[-1].max()/R_b.max()
    # the boundary point at the outboard midplane, zeta_F = 0: the first point
    i0 = 0
    rot = np.exp(-1j*np.angle(Phi[i0]))
    # boundary potential at the points of the plasma edge in this plane: of
    # the nearest boundary point (display only)
    if dim == 3:
        pe = np.stack([Rq[-1], np.zeros_like(Rq[-1]), Zq[-1]], 1)/L
    else:
        pe = np.stack([Rq[-1], Zq[-1]], 1)/L
    ke = np.argmin(((pe[:, None, :]-x[None, :, :])**2).sum(-1), 1)
    return dict(Rq=Rq, Zq=Zq, Qn=Qn*rot, Rv=Rv, Zv=Zv, phi=phi*rot, inpl=inpl, R_b=R_b*L, Z_b=Z_b*L,
                Phi=Phi*rot, dphi=dphi*rot, nrm=nrm, i0=i0, L=L, dim=dim, phi_edge=Phi[ke]*rot)


def midplane(d):
    """Potential along the outboard midplane outside the boundary, and the boundary value and slope there."""
    # the vacuum arrays are (Z, R) after h5py's reversal of the Fortran (R, Z) order
    Zrow = d['Zv'][:, 0]
    jz = np.argmin(np.abs(Zrow))
    R = d['Rv'][jz, :]
    ph = d['phi'][jz, :]
    out = d['inpl'][jz, :] < 0.5
    R_b0 = d['R_b'][d['i0']]
    m = out & (R > R_b0)
    Phi0 = d['Phi'][d['i0']]
    # dphi along the inward normal -J grad psi, which points in -R at the outboard midplane
    slope = -d['dphi'][d['i0']]/np.linalg.norm(d['nrm'][d['i0']])/d['L']
    return R[m] - R_b0, ph[m], Phi0, slope


def panel_map(ax, d, title):
    Qs = np.abs(d['Qn'][-1]).max()
    Ps = np.abs(d['Phi']).max()
    # the vacuum layer: in the plasma, where vac_pot_plot writes zero, the
    # potential of the nearest point of the plasma edge (display only, so
    # that the cells cut by the edge do not blend towards zero), clipped at
    # the edge of the plasma plot; imshow on the regular (R, Z) grid clips per
    # pixel, a gouraud-shaded mesh drops whole cells that touch the clip path
    Re, Ze, Pe = d['Rq'][-1], d['Zq'][-1], d['phi_edge']
    vac = d['phi'].real.copy()
    ins = d['inpl'] > 0.5
    k = np.argmin(np.hypot(d['Rv'][ins][:, None]-Re[None, :], d['Zv'][ins][:, None]-Ze[None, :]), 1)
    vac[ins] = Pe[k].real
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    m = ax.imshow(vac/Ps, cmap=DIVERGING, norm=norm, origin='lower', interpolation='bilinear',
                  extent=[d['Rv'].min(), d['Rv'].max(), d['Zv'].min(), d['Zv'].max()])
    box = [(d['Rv'].min(), d['Zv'].min()), (d['Rv'].max(), d['Zv'].min()), (d['Rv'].max(), d['Zv'].max()),
           (d['Rv'].min(), d['Zv'].max()), (d['Rv'].min(), d['Zv'].min())]
    edge = list(zip(Re[::-1], Ze[::-1])) + [(Re[-1], Ze[-1])]
    clip = Path(box + edge, [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY] + [Path.MOVETO] +
                [Path.LINETO]*(len(edge)-2) + [Path.CLOSEPOLY])
    m.set_clip_path(PathPatch(clip, transform=ax.transData, fc='none', ec='none'))
    ax.pcolormesh(d['Rq'], d['Zq'], d['Qn'].real/Qs, cmap=DIVERGING, norm=norm, shading='gouraud')
    ax.plot(np.append(Re, Re[0]), np.append(Ze, Ze[0]), color=INK2, lw=0.6)
    ax.set_aspect('equal')
    ax.set_xlabel('R')
    ax.set_ylabel('Z')
    ax.set_title(title, fontsize=9)
    ax.grid(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_hel')
    ap.add_argument('run_vmec')
    ap.add_argument('zoom_hel')
    ap.add_argument('zoom_vmec')
    ap.add_argument('--id', type=int, default=1)
    ap.add_argument('--out', default=os.path.dirname(os.path.abspath(__file__)))
    a = ap.parse_args()
    runs = [(a.run_hel, a.zoom_hel, 'HELENA, axisymmetric vacuum (vac%style 2)'),
            (a.run_vmec, a.zoom_vmec, 'VMEC twin, field-line vacuum (vac%style 1)')]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[1, 1, 1.1]))
    for k, (run, zoom, title) in enumerate(runs):
        d = load(run, a.id)
        panel_map(axs[k], d, title + '\ninside: Re $Q_n$, outside: Re $\\phi$ (each / max on the boundary)')
        s, ph, Phi0, slope = midplane(load(zoom, a.id))
        Ps = np.abs(d['Phi']).max()
        col = [BLUE, ORANGE][k]
        lab = ['HELENA', 'VMEC twin'][k]
        axs[2].plot(s/0.33, ph.real/Ps, color=col, label=lab + ': $\\phi$ on the plot grid')
        axs[2].plot([0], [Phi0.real/Ps], 'o', color=col, ms=5, label=lab + ': boundary potential')
        ss = np.linspace(0, 0.03, 5)
        axs[2].plot(ss/0.33, (Phi0.real + slope.real*ss)/Ps, '--', color=col, lw=1,
                    label=lab + ': slope from $-\\mathcal{J}Q_n$')
        for j in range(3):
            print('%s: Phi_s(outboard)/max = %.4f; at %.4f a from the boundary: phi/max = %.4f, tangent %.4f' %
                  (lab, Phi0.real/Ps, s[j]/0.33, ph[j].real/Ps, (Phi0.real + slope.real*s[j])/Ps))
        # decay rate of the potential just outside the boundary against the
        # one that the normal data give, -slope/Phi_s, both in units of 1/a
        rate = -np.polyfit(s[:10], np.log(np.abs(ph[:10].real)), 1)[0]
        print('%s: decay rate just outside %.3f / a, from the normal data %.3f / a' %
              (lab, rate*0.33, -slope.real/Phi0.real*0.33))
    axs[2].set_xlabel('distance from the boundary along the outboard midplane / a')
    axs[2].set_ylabel('Re $\\phi$ / max on the boundary')
    axs[2].set_xlim(0, 0.42)
    axs[2].legend(fontsize=7)
    axs[2].set_title('vacuum potential leaving the boundary', fontsize=9)
    fig.tight_layout()
    fn = os.path.join(a.out, 'vac_pot_circular.png')
    fig.savefig(fn, dpi=150)
    print('wrote', fn)


if __name__ == '__main__':
    main()
