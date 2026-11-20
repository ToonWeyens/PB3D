#!/usr/bin/env python3
"""Figures for the step-0 baseline of the vacuum bug-fix plan (issue #4).

    python make_figures.py --fixtures DIR [--runs DIR] [--out DIR]

--fixtures  directory with the decompressed HELENA fixtures cbm18a_extended_vac.12,
            circ_eps033_q104.12 and circ_eps005_q104.12 (default: $PB3D_FIXTURE_DIR)
--runs      directory holding the step-0 run directories fix075 and free075_s2, each
            with its PB3D_out.h5 (regression decks at psi_N 0.45 to 0.75, code 17b712c);
            without it the eigenfunction figure is skipped
--out       where the PNG files go (default: the directory of this script)

Needs numpy, matplotlib and h5py. The HELENA reader follows PB3D's read_HEL
(Modules/HELENA_ops.f90); lengths come out in units of R_mag, as in PB3D. The mode
numbers of the eigenvectors are rebuilt as setup_modes (Modules/X_ops.f90) assigns
them for X_style 2: at every normal point the 20 slots hold m = nint(n q - 9.5) + 0..19,
and slot(m) = (m - m_min(first point)) mod 20.
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

INK, INK2, MUTED, GRID, BASE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
BLUE, ORANGE = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': BASE, 'axes.labelcolor': INK2, 'xtick.color': INK2,
    'ytick.color': INK2, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
    'axes.spines.right': False, 'figure.facecolor': 'white', 'axes.facecolor': 'white', 'legend.frameon': False,
    'text.color': INK, 'axes.titlecolor': INK, 'lines.linewidth': 1.6})
BLUE_RAMP = LinearSegmentedColormap.from_list('blue', ['#86b6ef', '#2a78d6', '#184f95', '#0d366b'])
DIVERGING = LinearSegmentedColormap.from_list('div', ['#1c5cab', '#2a78d6', '#86b6ef', '#f0efec', '#f3a07f', '#eb6834', '#b8431a'])


# ----------------------------------------------------------------------------- HELENA mapping file
def read_helena(path):
    """Read a HELENA mapping file (fort.12, PB3D-patched format) in the order of read_HEL."""
    tok = open(path).read().split()
    pos = [2]

    def take(n):
        out = np.array(tok[pos[0]:pos[0] + n], dtype=float)
        pos[0] += n
        return out

    js0, ias = int(tok[0]), int(tok[1])
    nr = js0 + 1
    s = take(nr)
    q = take(nr)
    take(2)                                   # dq on axis and edge (finite-element coordinate)
    take(nr - 1)                              # dq
    curj = take(nr)
    take(2)                                   # dj on axis and edge
    nchi = int(float(tok[pos[0]])); pos[0] += 1
    chi = take(nchi)
    nmap = (nr - 1) * nchi
    take(nmap); take(nmap)                    # gem11, gem12
    cpsurf, radius = take(2)
    take(nmap)                                # gem33
    raxis, b0 = take(2)
    p = take(nr)
    take(2)                                   # dp on axis and edge
    f = take(nr)
    take(2)                                   # dF on axis and edge (unused)
    vx = take(nchi); vy = take(nchi)
    eps = take(1)[0]
    xout = take(nmap).reshape(nr - 1, nchi)
    yout = take(nmap).reshape(nr - 1, nchi)
    assert pos[0] == len(tok), 'HELENA file not read to the end: %d of %d tokens' % (pos[0], len(tok))

    a = radius * raxis                        # minor radius of the HELENA boundary over R_mag
    R = np.empty((nr, nchi)); Z = np.empty((nr, nchi))
    R[1:] = a * (1.0 / eps + xout); Z[1:] = a * yout
    R[0] = raxis; Z[0] = 0.0
    return dict(nr=nr, ias=ias, nchi=nchi, s=s, psi_n=s**2, q=q, curj=curj, chi=chi, cpsurf=cpsurf,
                radius=radius, raxis=raxis, b0=b0, p=p, f=f, eps=eps, a=a, R=R, Z=Z)


def full_RZ(h):
    """Mirror the top half of a top-bottom symmetric equilibrium (ias = 0) to a closed cross-section."""
    R, Z = h['R'], h['Z']
    if h['ias'] == 0:
        R = np.concatenate([R, R[:, -2::-1]], axis=1); Z = np.concatenate([Z, -Z[:, -2::-1]], axis=1)
    return R, Z


def plasma_edge(h):
    """Normalised flux where the pressure gradient vanishes (start of the pseudo-vacuum)."""
    psi, p = h['psi_n'], h['p']; dp = np.gradient(p, psi)
    idx = np.where(np.abs(dp) > 1e-3 * np.abs(dp).max())[0]
    return psi[idx[-1]]


def surf(h, R, Z, psi_n):
    k = np.argmin(abs(h['psi_n'] - psi_n)); return R[k], Z[k]


# ----------------------------------------------------------------------------- figure 1: cbm18a cross-section
def fig_cross_section(h, out):
    R, Z = full_RZ(h); pe = plasma_edge(h)
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    for pn in np.arange(0.05, 1.0, 0.05):
        r, z = surf(h, R, Z, pn); ax.plot(r, z, color=GRID, lw=0.7, zorder=1)
    r1, z1 = surf(h, R, Z, 0.45); r2, z2 = surf(h, R, Z, 0.75)
    ax.fill(np.concatenate([r2, r1[::-1]]), np.concatenate([z2, z1[::-1]]), color=BLUE, alpha=0.10, lw=0, zorder=0)
    for pn, col, lab, ls in [(0.45, BLUE, r'$\psi_N$ = 0.45, inner edge of the solution window', '-'),
                             (pe, ORANGE, r'$\psi_N$ = %.2f, plasma edge ($p^\prime$ = 0 beyond)' % pe, '-'),
                             (0.75, BLUE, r'$\psi_N$ = 0.75, deck boundary (fixed / free)', '--'),
                             (1.0, INK, r'$\psi_N$ = 1, HELENA boundary (pseudo-vacuum)', '-')]:
        r, z = surf(h, R, Z, pn); ax.plot(r, z, color=col, lw=1.6, ls=ls, label=lab, zorder=3)
    ax.plot([1.0], [0.0], marker='+', color=INK, ms=9, mew=1.4, ls='none', label='magnetic axis', zorder=4)
    ax.set_aspect('equal'); ax.set_xlabel(r'$R\,/\,R_{\rm mag}$'); ax.set_ylabel(r'$Z\,/\,R_{\rm mag}$')
    ax.set_title('cbm18a with extended pseudo-vacuum: HELENA flux surfaces (every 0.05 in $\\psi_N$)', fontsize=9.5, loc='left')
    ax.legend(loc='upper left', fontsize=7.6, handlelength=1.8)
    ax.set_xlim(0.2, 1.45); ax.set_ylim(-0.6, 0.6)
    fig.tight_layout(); fig.savefig(os.path.join(out, 'cbm18a_cross_section.png'), dpi=130); plt.close(fig)
    return pe


# ----------------------------------------------------------------------------- figure 2: cbm18a profiles
def fig_profiles(h, pe, out):
    psi, q, p, F = h['psi_n'], h['q'], h['p'], h['f']
    shear = 2 * psi / q * np.gradient(q, psi)               # s = 2 (psi/q) dq/dpsi, the flux form of (r/q) dq/dr
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.4), sharex=True)

    def mark_windows(ax):
        ax.axvspan(0.45, 0.75, color=BLUE, alpha=0.10, lw=0)
        ax.axvline(pe, color=ORANGE, lw=1.2); ax.axvline(0.75, color=BLUE, lw=1.0, ls='--')

    a = axs[0, 0]; mark_windows(a); a.plot(psi, p / p[0], color=INK); a.set_ylabel(r'$p\,/\,p(0)$')
    a.set_title('pressure', loc='left', fontsize=9)
    a.text(0.77, 0.62, 'plasma edge\n$\\psi_N$ = %.2f' % pe, color=ORANGE, fontsize=7.6, va='center')
    a.text(0.46, 0.50, 'solution window\n0.45 to 0.75', color=BLUE, fontsize=7.6, va='center')
    a = axs[0, 1]; mark_windows(a); a.plot(psi, q, color=INK); a.set_ylabel('$q$'); a.set_ylim(0, 10)
    a.set_title('safety factor (q = %.2f on axis, %.1f at $\\psi_N$ = 1)' % (q[0], q[-1]), loc='left', fontsize=9)
    win = (psi >= 0.45) & (psi <= 0.75)
    ms = np.arange(np.ceil(10 * q[win].min()), np.floor(10 * q[win].max()) + 1)
    a.plot(np.interp(ms / 10, q, psi), ms / 10, marker='|', color=BLUE, ms=7, mew=1.2, ls='none',
           label='rational surfaces $q = m/10$, m = %d to %d' % (ms[0], ms[-1]))
    a.legend(loc='upper left', fontsize=7.4)
    a = axs[1, 0]; mark_windows(a); a.plot(psi, shear, color=INK); a.set_ylabel(r'$s = 2\,(\psi/q)\,{\rm d}q/{\rm d}\psi$')
    a.set_title('magnetic shear', loc='left', fontsize=9); a.set_xlabel(r'$\psi_N$'); a.set_ylim(0, 12)
    a = axs[1, 1]; mark_windows(a); a.plot(psi, F / F[0], color=INK); a.set_ylabel(r'$F = R B_\phi$, over $F(0)$')
    a.set_title('poloidal current function', loc='left', fontsize=9); a.set_xlabel(r'$\psi_N$')
    for a in axs.ravel(): a.set_xlim(0, 1)
    fig.suptitle('cbm18a profiles against normalised poloidal flux, from the HELENA fixture (%d radial points)' % h['nr'],
                 fontsize=9.5, x=0.01, ha='left')
    fig.tight_layout(); fig.savefig(os.path.join(out, 'cbm18a_profiles.png'), dpi=130); plt.close(fig)


# ----------------------------------------------------------------------------- figure 3: circular external-kink cases
def fig_circular(fixtures, out):
    cases = [('circ_eps033_q104.12', r'$\varepsilon$ = 0.33, flat current', BLUE),
             ('circ_eps005_q104.12', r'$\varepsilon$ = 0.05, $\Gamma = 1 - 0.02\,\psi$', ORANGE)]
    fig = plt.figure(figsize=(7.4, 6.2)); gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.95])
    hs = []
    for i, (fn, lab, col) in enumerate(cases):
        hc = read_helena(os.path.join(fixtures, fn)); hs.append(hc); Rc, Zc = full_RZ(hc); ax = fig.add_subplot(gs[0, i])
        for pn in np.arange(0.1, 1.0, 0.1):
            r, z = surf(hc, Rc, Zc, pn); ax.plot(r, z, color=GRID, lw=0.7)
        ax.plot(Rc[-1], Zc[-1], color=col, lw=1.6); ax.plot([1], [0], marker='+', color=INK, ms=8, mew=1.3, ls='none')
        ax.set_aspect('equal'); ax.set_title(lab, loc='left', fontsize=9, color=col)
        ax.set_xlabel(r'$R\,/\,R_{\rm mag}$'); ax.set_ylabel(r'$Z\,/\,R_{\rm mag}$')
        cen = (hc['R'].max() + hc['R'].min()) / 2; half = 1.12 * hc['a']
        ax.set_xlim(cen - half, cen + half); ax.set_ylim(-half, half); ax.tick_params(labelsize=7.5)
    ax = fig.add_subplot(gs[1, :])
    for hc, (fn, lab, col) in zip(hs, cases):
        ax.plot(hc['psi_n'], hc['q'], color=col, label=lab + ': $q_0$ = %.4f, $q_a$ = %.4f' % (hc['q'][0], hc['q'][-1]))
    for n, m in [(10, 11), (15, 16), (20, 21)]:
        ax.axhline(m / n, color=MUTED, lw=0.9, ls=':')
        ax.text(0.01, m / n + 0.001, 'q = %d/%d (lowest rational surface of the n = %d deck)' % (m, n, n), color=INK2, fontsize=7.4, va='bottom')
    ax.set_xlim(0, 1); ax.set_ylim(1.02, 1.115); ax.set_xlabel(r'$\psi_N$'); ax.set_ylabel('$q$')
    ax.set_title('safety factor: below every m/n of the decks, so the modes are external kinks', loc='left', fontsize=9)
    ax.legend(loc='center', bbox_to_anchor=(0.62, 0.66), fontsize=7.4)
    fig.suptitle('R2 circular cases: HELENA flux surfaces (every 0.1 in $\\psi_N$) and safety factor', fontsize=9.5, x=0.01, ha='left')
    fig.tight_layout(); fig.savefig(os.path.join(out, 'circular_cases.png'), dpi=130); plt.close(fig)


# ----------------------------------------------------------------------------- figure 4: baseline eigenfunctions
def nint(x):
    return np.floor(x + 0.5).astype(int)      # Fortran nint for positive arguments


def load_run(path, n=10, n_mod=20):
    import h5py
    with h5py.File(os.path.join(path, 'PB3D_out.h5')) as f:
        rF_eq = f['grid_eq/r_F/var'][()]; q_eq = f['eq_1/q_saf_FD/var'][()][:rF_eq.size]     # derivative 0
        rF = f['grid_sol/r_F/var'][()]
        X = (f['sol_R_1/RE_sol_vec/var'][()] + 1j * f['sol_R_1/IM_sol_vec/var'][()]).reshape(rF.size, n_mod)
        ev = f['sol_R_1/RE_sol_val/var'][()][0]
    min_m_eq = nint(n * q_eq - (n_mod - 1) * 0.5)
    min_m = nint(np.interp(rF, rF_eq, min_m_eq.astype(float)))
    slots = np.arange(n_mod)
    m_of = min_m[:, None] + np.mod(slots[None, :] - (min_m[:, None] - min_m[0]), n_mod)
    k = np.unravel_index(np.argmax(np.abs(X)), X.shape)
    X = X * np.exp(-1j * np.angle(X[k])); X = X / np.abs(X).max()                           # one global phase, max 1
    return dict(rF=rF, X=X, m_of=m_of, ev=ev)


def fig_eigenfunctions(h, runs, out):
    cps = h['cpsurf']; pe = plasma_edge(h)
    specs = [('fix075', 'fixed boundary (BC_style 1 1)'), ('free075_s2', 'free boundary, vacuum style 2 (BC_style 1 2)')]
    data = [load_run(os.path.join(runs, r)) for r, _ in specs]
    all_m = sorted(set(np.unique(np.concatenate([d['m_of'].ravel() for d in data]))))
    mlo, mhi = min(all_m), max(all_m)
    chi = h['chi']; chi_full = np.concatenate([chi, 2 * np.pi - chi[-2::-1]])
    Rh, Zh = full_RZ(h)
    fig, axs = plt.subplots(2, 2, figsize=(8.0, 7.6), gridspec_kw=dict(height_ratios=[0.9, 1.25]))
    for i, (d, (run, lab)) in enumerate(zip(data, specs)):
        psi = d['rF'] / cps; X = d['X']; m_of = d['m_of']
        ax = axs[0, i]
        ax.axvline(pe, color=ORANGE, lw=1.0); ax.axvline(0.75, color=BLUE, lw=1.0, ls='--')
        for m in all_m:
            kd, sl = np.where(m_of == m)
            if kd.size == 0: continue
            y = np.full(psi.size, np.nan); y[kd] = np.abs(X[kd, sl])
            col = BLUE_RAMP((m - mlo) / (mhi - mlo)); ax.plot(psi, y, color=col, lw=1.1)
            if np.nanmax(y) > 0.25:
                j = np.nanargmax(y); ax.text(psi[j], y[j] + 0.02, str(m), color=col, fontsize=6.8, ha='center')
        ax.set_xlim(0.45, 0.78); ax.set_ylim(0, 1.08); ax.set_xlabel(r'$\psi_N$')
        ax.set_title('%s: $\\omega^2/\\omega_A^2$ = %.4e' % (lab, d['ev']), loc='left', fontsize=8.6)
        if i == 0: ax.set_ylabel(r'$|X_m|\,/\,\max|X|$')
        ax.text(pe, 1.0, ' plasma edge', color=ORANGE, fontsize=7, va='top', ha='right', rotation=90)
        ax.text(0.75, 1.0, ' deck boundary', color=BLUE, fontsize=7, va='top', ha='right', rotation=90)
        # poloidal cross-section at zeta = 0: Re sum_m X_m(psi) exp(-i m theta), theta the HELENA angle
        Rs = np.array([np.interp(psi, h['psi_n'], Rh[:, j]) for j in range(chi_full.size)]).T
        Zs = np.array([np.interp(psi, h['psi_n'], Zh[:, j]) for j in range(chi_full.size)]).T
        field = np.zeros((psi.size, chi_full.size))
        for m in all_m:
            kd, sl = np.where(m_of == m)
            field[kd] += np.real(X[kd, sl][:, None] * np.exp(-1j * m * chi_full)[None, :])
        ax = axs[1, i]; ax.grid(False)
        pc = ax.pcolormesh(Rs, Zs, field / np.abs(field).max(), cmap=DIVERGING, norm=TwoSlopeNorm(0, -1, 1),
                           shading='gouraud', rasterized=True)
        k70 = np.argmin(abs(h['psi_n'] - pe)); ax.plot(Rh[k70], Zh[k70], color=ORANGE, lw=1.0)
        ax.plot(Rs[-1], Zs[-1], color=BLUE, lw=1.0, ls='--'); ax.plot(Rs[0], Zs[0], color=BLUE, lw=0.8)
        ax.set_aspect('equal'); ax.set_xlabel(r'$R\,/\,R_{\rm mag}$')
        if i == 0: ax.set_ylabel(r'$Z\,/\,R_{\rm mag}$')
        ax.set_title(r'Re $X(R,Z)$ at $\zeta$ = 0, same run', loc='left', fontsize=8.6)
        ax.set_xlim(0.55, 1.36); ax.set_ylim(-0.42, 0.42)
    cb = fig.colorbar(pc, ax=axs[1, :].tolist(), fraction=0.025, pad=0.02); cb.set_label(r'Re $X$ / max', color=INK2)
    cb.outline.set_visible(False)
    fig.suptitle('Step-0 baseline eigenfunctions, cbm18a, n = 10, 20 poloidal modes following $nq$ (m = %d to %d), code 17b712c unchanged'
                 % (mlo, mhi), fontsize=9.2, x=0.01, ha='left')
    fig.savefig(os.path.join(out, 'cbm18a_eigenfunctions.png'), dpi=130, bbox_inches='tight'); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fixtures', default=os.environ.get('PB3D_FIXTURE_DIR'), help='directory with the HELENA fixtures')
    ap.add_argument('--runs', default=None, help='directory with the step-0 run directories fix075 and free075_s2')
    ap.add_argument('--out', default=os.path.dirname(os.path.abspath(__file__)), help='output directory')
    args = ap.parse_args()
    if not args.fixtures:
        sys.exit('give --fixtures or set PB3D_FIXTURE_DIR')
    os.makedirs(args.out, exist_ok=True)
    h = read_helena(os.path.join(args.fixtures, 'cbm18a_extended_vac.12'))
    pe = fig_cross_section(h, args.out)
    fig_profiles(h, pe, args.out)
    fig_circular(args.fixtures, args.out)
    if args.runs:
        fig_eigenfunctions(h, args.runs, args.out)
    else:
        print('no --runs given: eigenfunction figure skipped')
    print('figures written to', args.out)


if __name__ == '__main__':
    main()
