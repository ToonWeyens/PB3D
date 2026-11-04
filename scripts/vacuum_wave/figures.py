"""Data and figures of Documentation/vacuum_wave_report (issue #4, exploration only).

    .venv/bin/python scripts/vacuum_wave/figures.py data    # about 9 minutes on a Mac
    .venv/bin/python scripts/vacuum_wave/figures.py plot    # PDFs next to report.tex

`data` writes figdata.json next to the report; `plot` reads it. Every number
of the report comes from that file. The parts of `data` can be rerun alone and
rewrite only their own entry: `single` (single-wave data, bands of A and of the
mass matrix), `compare` (wave form against point values directly), `memsplit`
(where the memory goes at 64 x 801).
"""
import json
import pathlib
import resource
import subprocess
import sys
import time

import numpy as np

import prototype as pt

OUT = pathlib.Path(__file__).resolve().parents[2]/'Documentation'/'vacuum_wave_report'
DATA = OUT/'figdata.json'
K_MAIN = 12                                     # sidebands on each side of the window
COLORS = {6: '#0072B2', 11: '#D55E00', 15: '#009E73'}   # Okabe-Ito, one per test harmonic


# --- data ----------------------------------------------------------------------

def spectrum():
    """Fourier content in m of the exact Phi on the torus, at phi = 0."""
    th = np.linspace(-np.pi, np.pi, 1024, endpoint=False)
    ms = np.arange(-20, 46)
    out = {'m': ms.tolist(), 'theta': th[::8].tolist()}
    for m0 in pt.M_TEST:
        f = pt.exact_on_torus(th, 0*th, m0)[0]
        c = np.abs(np.exp(1j*np.outer(ms, th))@f)/len(th)          # coefficient of exp(-i m theta)
        out[f'spec{m0}'] = (c/c.max()).tolist()
        f = f/np.abs(f).max()
        out[f'phi{m0}'] = f.real[::8].tolist()
        out[f'wave{m0}'] = np.cos(m0*th[::8]).tolist()
    return out


def matrices():
    """|H| of the point-value method on a small grid, |A| of the wave form at 64 x 201."""
    g = pt.Grid(16, 21)
    every = np.arange(g.n)
    Kw, _ = pt.kernel_rows(g, every, every)
    H = -Kw
    H[np.diag_indices(g.n)] = 4*np.pi+Kw.sum(1)
    g = pt.Grid(64, 201)
    ws = pt.wave_sums_from_scratch(g, pt.window(K_MAIN), g.data()[1])
    Amat = ws.system()[0]
    return {'H_small': round3(np.abs(H)), 'A': round3(np.abs(Amat)), 'A_m': ws.m_list.tolist()}


def round3(M):
    """Three significant digits: the matrices only feed grey-scale images."""
    return [[float(f'{v:.3g}') for v in row] for row in np.asarray(M)]


def sidebands(e_ex):
    g = pt.Grid(64, 201)
    ph_ex, qd = g.data()
    ws = pt.wave_sums_from_scratch(g, pt.window(20), qd)
    rows = []
    for k in range(21):
        sel = np.isin(ws.m_list, pt.window(k))
        sub = pt.WaveSums(g, ws.m_list[sel], qd)
        sub.A_off, sub.B_off = ws.A_off[np.ix_(sel, sel)], ws.B_off[sel]
        sub.S, sub.Kx, sub.Gn = ws.S, ws.Kx, ws.Gn
        e_phi, e_en = pt.errors(g, sub.solve(), ph_ex, qd, e_ex)
        rows.append({'k': k, 'waves': int(sel.sum()), 'phi': e_phi.tolist(), 'energy': e_en.tolist()})
    return rows


def parallel_points(e_ex):
    rows, prev = [], None
    for n_par in (51, 101, 201, 401, 801):
        g = pt.Grid(64, n_par)
        ph_ex, qd = g.data()
        t0 = time.time()
        ws = pt.wave_sums_from_scratch(g, pt.window(K_MAIN), qd)
        row = {'n_par': n_par, 'N': g.n, 't_scratch': time.time()-t0}
        e_phi, e_en = pt.errors(g, ws.solve(), ph_ex, qd, e_ex)
        row.update(phi=e_phi.tolist(), energy=e_en.tolist())
        if prev is not None:
            t0 = time.time()
            wc = pt.wave_sums_carried(prev, g, pt.window(K_MAIN), qd)
            row['t_carried'] = time.time()-t0
            row['carry_diff'] = max(abs(wc.A_off-ws.A_off).max()/abs(ws.A_off).max(),
                                    abs(wc.B_off-ws.B_off).max()/abs(ws.B_off).max(),
                                    abs(wc.S-ws.S).max()/abs(ws.S).max())
        if g.n <= 13000:
            e_phi, e_en = pt.errors(g, pt.solve_points(g, qd), ph_ex, qd, e_ex)
            row.update(phi_points=e_phi.tolist(), energy_points=e_en.tolist())
        rows.append(row)
        prev = ws
        print(f'  parallel points {n_par} done', flush=True)
    return rows


def field_lines(e_ex):
    rows = []
    for n_lines in (16, 32, 64, 128):
        g = pt.Grid(n_lines, 201)
        ph_ex, qd = g.data()
        ws = pt.wave_sums_from_scratch(g, pt.window(K_MAIN), qd)
        e_phi, e_en = pt.errors(g, ws.solve(), ph_ex, qd, e_ex)
        row = {'lines': n_lines, 'N': g.n, 'phi': e_phi.tolist(), 'energy': e_en.tolist()}
        if g.n <= 13000:
            e_phi, e_en = pt.errors(g, pt.solve_points(g, qd), ph_ex, qd, e_ex)
            row.update(phi_points=e_phi.tolist(), energy_points=e_en.tolist())
        rows.append(row)
        print(f'  field lines {n_lines} done', flush=True)
    return rows


def single_wave():
    """Potential for normal data that are single waves, as PB3D's plasma modes are."""
    g = pt.Grid(64, 201)
    ms = pt.window(20)
    ws = pt.wave_sums_from_scratch(g, ms, g.waves(np.array(pt.M_TEST)))
    a = np.linalg.solve(*ws.system())
    out = {'offset': list(range(-8, 9))}
    for c, m0 in enumerate(pt.M_TEST):
        amp = abs(a[:, c])/abs(a[:, c]).max()
        out[f'amp{m0}'] = [float(amp[list(ms).index(m0+d)]) for d in out['offset']]
    return out


def matrix_bands(A):
    """Median ratio of the k-th off-diagonal of |A| to its diagonal."""
    dia = np.diag(A)
    return [float(np.median(np.concatenate([np.diag(A, k)/dia[:-k], np.diag(A, -k)/dia[k:]])))
            for k in (1, 2, 3)]


def mass_bands():
    """The part of A that the surface element alone gives: M = E^H I E on 64 x 201,
    and what remains of A once 2 pi M is taken out (the double-layer part)."""
    g = pt.Grid(64, 201)
    ws = pt.wave_sums_from_scratch(g, pt.window(K_MAIN), g.data()[1])
    Amat = ws.system()[0]
    M = (g.w[:, None]*np.conj(ws.E)).T@ws.E
    rest = Amat-2*np.pi*M
    return {'mass': matrix_bands(np.abs(M)), 'A': matrix_bands(np.abs(Amat)),
            'rest_over_A_diag': [float(np.median(np.abs(np.diag(rest))/np.abs(np.diag(Amat))))]
            + matrix_bands_against(np.abs(rest), np.abs(np.diag(Amat)))}


def matrix_bands_against(B, dia):
    """Median ratio of the k-th off-diagonal of B to a given diagonal."""
    return [float(np.median(np.concatenate([np.diag(B, k)/dia[:-k], np.diag(B, -k)/dia[k:]])))
            for k in (1, 2, 3)]


def add_single():
    """Add the single-wave results and the bands of |A| to an existing figdata.json."""
    d = json.loads(DATA.read_text())
    d['single_wave'] = single_wave()
    d['A_bands'] = matrix_bands(np.array(d['matrices']['A']))
    d['mass_bands'] = mass_bands()
    DATA.write_text(json.dumps(d))
    print(d['A_bands'], d['mass_bands'], d['single_wave'])


def compare():
    """Direct difference between the wave form and the point-value form on the same rows."""
    d = json.loads(DATA.read_text())
    e_ex = np.array([pt.exact_energy(m) for m in pt.M_TEST])
    out = []
    for n_lines, n_par in ((64, 101), (64, 201)):
        g = pt.Grid(n_lines, n_par)
        ph_ex, qd = g.data()
        ph_p = pt.solve_points(g, qd)
        w = g.w[:, None]
        en_p = -(w*np.conj(qd)*ph_p).sum(0)
        for k in (6, 12, 20):
            ph_w = pt.wave_sums_from_scratch(g, pt.window(k), qd).solve()
            en_w = -(w*np.conj(qd)*ph_w).sum(0)
            out.append({'lines': n_lines, 'n_par': n_par, 'k': k,
                        'energy_diff': (abs(en_w-en_p)/abs(en_p)).tolist(),
                        'phi_diff': np.sqrt((w*abs(ph_w-ph_p)**2).sum(0)/(w*abs(ph_p)**2).sum(0)).tolist(),
                        'energy_err_points': (abs(en_p-e_ex)/abs(e_ex)).tolist()})
            print(out[-1], flush=True)
    d['compare'] = out
    DATA.write_text(json.dumps(d))


def memsplit():
    """Where the memory goes at the largest grid: the data alone, then the pair sums with two block sizes."""
    d = json.loads(DATA.read_text())
    d['memsplit'] = []
    for stage, chunk in (('data', 256), ('all', 256), ('all', 64)):
        res = subprocess.run([sys.executable, __file__, 'mem', '801', stage, str(chunk)],
                             capture_output=True, text=True, check=True)
        d['memsplit'].append(json.loads(res.stdout.strip().splitlines()[-1]))
        print(d['memsplit'][-1], flush=True)
    DATA.write_text(json.dumps(d))


def memory_run(n_par, stage='all', chunk=256):
    """One process per grid, so that the peak memory belongs to that grid alone."""
    g = pt.Grid(64, n_par)
    t0 = time.time()
    qd = g.data()[1]
    if stage == 'all':
        ws = pt.WaveSums(g, pt.window(K_MAIN), qd)
        ws.add_pairs(np.arange(g.n), np.arange(g.n), chunk=chunk)
        ws.solve()
    wall = time.time()-t0
    ru = resource.getrusage(resource.RUSAGE_SELF)
    print(json.dumps({'n_par': n_par, 'N': g.n, 'stage': stage, 'chunk': chunk, 'seconds': wall,
                      'cpu_fraction': (ru.ru_utime+ru.ru_stime)/wall, 'peak_MiB': pt.peak_mb()}))


def data():
    e_ex = np.array([pt.exact_energy(m) for m in pt.M_TEST])
    out = {'harmonic_check': [pt.check_harmonic(m) for m in pt.M_TEST], 'spectrum': spectrum()}
    print('spectrum done', flush=True)
    out['matrices'] = matrices()
    out['sidebands'] = sidebands(e_ex)
    print('sidebands done', flush=True)
    out['parallel'] = parallel_points(e_ex)
    out['lines'] = field_lines(e_ex)
    out['memory'] = []
    for n_par in (51, 101, 201, 401, 801):
        res = subprocess.run([sys.executable, __file__, 'mem', str(n_par)],
                             capture_output=True, text=True, check=True)
        out['memory'].append(json.loads(res.stdout.strip().splitlines()[-1]))
        print(f'  memory {n_par} done', flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.write_text(json.dumps(out))
    add_single()
    compare()
    memsplit()


# --- figures -------------------------------------------------------------------

def setup():
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'text.usetex': True, 'text.latex.preamble': r'\usepackage{amsmath}',
                         'font.family': 'serif', 'font.size': 8,
                         'axes.labelsize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7,
                         'ytick.labelsize': 7, 'lines.linewidth': 1.0, 'axes.linewidth': 0.5,
                         'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                         'legend.frameon': False, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})
    return plt


def fig_grid(plt):
    """Field lines on the unrolled torus; one Richardson step."""
    fig, ax = plt.subplots(figsize=(3.4, 2.3))
    n_lines, n_old = 6, 7
    for k in range(n_lines):
        th = np.linspace(-np.pi, np.pi, 400)
        ph = 2*np.pi*k/n_lines+pt.Q_EDGE*th
        wrap = np.mod(ph, 2*np.pi)
        cut = np.where(np.diff(wrap) < 0)[0]+1
        for seg_th, seg_ph in zip(np.split(th, cut), np.split(wrap, cut)):
            ax.plot(seg_ph, seg_th, color='0.7', lw=0.6, zorder=1)
        th_new = np.linspace(-np.pi, np.pi, 2*n_old-1)
        ph_new = np.mod(2*np.pi*k/n_lines+pt.Q_EDGE*th_new, 2*np.pi)
        ax.scatter(ph_new[1::2], th_new[1::2], s=14, facecolor='white', edgecolor='#D55E00', lw=0.8,
                   zorder=3, label='new points' if k == 0 else None)
        ax.scatter(ph_new[::2], th_new[::2], s=14, color='k', zorder=3, label='old points' if k == 0 else None)
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(-np.pi*1.05, np.pi*1.05)
    ax.set_xticks([0, np.pi, 2*np.pi], ['$0$', r'$\pi$', r'$2\pi$'])
    ax.set_yticks([-np.pi, 0, np.pi], [r'$-\pi$', '$0$', r'$\pi$'])
    ax.set_xlabel(r'toroidal angle $\varphi$')
    ax.set_ylabel(r'poloidal angle $\theta$')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.16), ncol=2, handletextpad=0.2)
    fig.savefig(OUT/'fig_grid.pdf')


def fig_spectrum(plt, d):
    s = d['spectrum']
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.4, 3.6))
    th = np.array(s['theta'])
    a1.plot(th, s['phi11'], color=COLORS[11], label=r'$\Re\,\Phi$, exact, $m=11$')
    a1.plot(th, s['wave11'], color='0.4', ls='--', lw=0.8, label=r'single wave $\cos 11\theta$')
    a1.set_xlim(-np.pi, np.pi)
    a1.set_xticks([-np.pi, 0, np.pi], [r'$-\pi$', '$0$', r'$\pi$'])
    a1.set_xlabel(r'poloidal angle $\theta$ (at $\varphi = 0$)')
    a1.set_ylabel(r'$\Re\,\Phi / \max|\Phi|$')
    a1.legend(loc='lower left', ncol=1, bbox_to_anchor=(0, 1.0))
    ms = np.array(s['m'])
    lo, hi = min(pt.M_TEST)-K_MAIN, max(pt.M_TEST)+K_MAIN
    a2.axvspan(lo-0.5, hi+0.5, color='0.92', zorder=0)
    for m0, mk in zip(pt.M_TEST, 'osD'):
        a2.semilogy(ms, s[f'spec{m0}'], marker=mk, ms=2.5, color=COLORS[m0], label=f'$m={m0}$')
    a2.set_ylim(1e-7, 2)
    a2.set_xlim(ms[0], ms[-1])
    a2.set_xlabel(r'poloidal mode number $m$ of the wave $e^{\imath(n\varphi - m\theta)}$')
    a2.set_ylabel('amplitude / largest')
    a2.legend(loc='upper left', title='exact harmonic', title_fontsize=7)
    fig.tight_layout(h_pad=0.6)
    fig.savefig(OUT/'fig_spectrum.pdf')


def fig_matrices(plt, d):
    from matplotlib.colors import LogNorm
    m = d['matrices']
    H, A = np.array(m['H_small']), np.array(m['A'])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={'width_ratios': [1, 1]})
    im = a1.imshow(H/H.max(), norm=LogNorm(1e-4, 1), cmap='Greys', interpolation='nearest')
    a1.set_xlabel(r'source point $j$')
    a1.set_ylabel(r'collocation point $i$')
    a1.set_title(rf'point values: $|H_{{ij}}|$, $N = {H.shape[0]}$ (16 lines $\times$ 21 points)', fontsize=8)
    fig.colorbar(im, ax=a1, fraction=0.046, pad=0.03, label='relative size')
    ms = np.array(m['A_m'])
    im = a2.imshow(A/A.max(), norm=LogNorm(1e-4, 1), cmap='Greys', interpolation='nearest',
                   extent=[ms[0]-0.5, ms[-1]+0.5, ms[-1]+0.5, ms[0]-0.5])
    a2.set_xlabel(r'wave $m_k$')
    a2.set_ylabel(r'test wave $m_l$')
    a2.set_title(rf'waves: $|A_{{lk}}|$, $M^\prime = {len(ms)}$ (any $N$; here 64 $\times$ 201)', fontsize=8)
    fig.colorbar(im, ax=a2, fraction=0.046, pad=0.03, label='relative size')
    fig.tight_layout(w_pad=2.0)
    fig.savefig(OUT/'fig_matrices.pdf')


def fig_carry(plt):
    """Pair sums of one Richardson step, old points listed first."""
    from matplotlib.patches import Rectangle
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    blue, orange = '#cfe3f1', '#f6d5c0'
    blocks = [(0, 1, blue, (0.3, 1.3), '$\\tfrac14\\times$\nprevious\nlevel'),
              (1, 1, orange, (1.5, 1.5), 'new\npairs'), (0, 0, orange, (0.5, 0.5), 'new\npairs'),
              (1, 0, orange, (1.33, 0.32), 'new\npairs')]
    for x, y, c, xy, t in blocks:
        ax.add_patch(Rectangle((x, y), 1, 1, facecolor=c, edgecolor='k', lw=0.5))
        ax.text(*xy, t, ha='center', va='center', fontsize=7, linespacing=1.1)
    ax.plot([0, 2], [2, 0], color='0.3', lw=2.5, alpha=0.5)
    ax.text(0.5, 2.06, r'old $j$', ha='center', va='bottom', fontsize=7)
    ax.text(1.5, 2.06, r'new $j$', ha='center', va='bottom', fontsize=7)
    ax.text(-0.06, 1.5, r'old $i$', ha='right', va='center', fontsize=7)
    ax.text(-0.06, 0.5, r'new $i$', ha='right', va='center', fontsize=7)
    ax.annotate('', xy=(2.5, 1.0), xytext=(2.1, 1.0), arrowprops={'arrowstyle': '->', 'lw': 0.7})
    ax.add_patch(Rectangle((2.55, 0.75), 0.5, 0.5, facecolor='white', edgecolor='k', lw=0.5))
    ax.text(2.8, 1.0, r'$\underline{\mathrm{A}}$', ha='center', va='center', fontsize=8)
    ax.text(2.8, 0.66, r"$M^\prime\times M^\prime$", ha='center', va='top', fontsize=7)
    ax.text(1.0, -0.08, r'grey diagonal: self terms, from per-point sums ($\tfrac12\times$ old $+$ new $j$)',
            ha='center', va='top', fontsize=6.5)
    ax.set_xlim(-0.45, 3.15)
    ax.set_ylim(-0.3, 2.22)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.savefig(OUT/'fig_carry.pdf')


def fig_accuracy(plt, d):
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.35), sharey=True)
    sb = d['sidebands']
    for i, m0 in enumerate(pt.M_TEST):
        axs[0].semilogy([r['k'] for r in sb], [r['energy'][i] for r in sb], marker='o', ms=2.5,
                        color=COLORS[m0], label=f'$m={m0}$')
    axs[0].set_xlabel(r'sidebands $k$ on each side of $m = 6\ldots15$')
    axs[0].set_ylabel(r'relative error of $-\!\int \bar q\,\Phi\,\mathrm{d}S$')
    axs[0].axvline(K_MAIN, color='0.6', lw=0.5, ls=':')
    axs[0].legend(loc='upper right')
    for ax, rows, key, xl in ((axs[1], d['parallel'], 'n_par', 'parallel points per field line (64 lines)'),
                              (axs[2], d['lines'], 'lines', 'field lines (201 points per line)')):
        x = [r[key] for r in rows]
        for i, m0 in enumerate(pt.M_TEST):
            ax.loglog(x, [r['energy'][i] for r in rows], marker='o', ms=2.5, color=COLORS[m0],
                      label='waves' if i == 0 else None)
            xp = [r[key] for r in rows if 'energy_points' in r]
            ax.loglog(xp, [r['energy_points'][i] for r in rows if 'energy_points' in r], ls='none',
                      marker='s', ms=6, mfc='none', mec=COLORS[m0], mew=0.6, label='points' if i == 0 else None)
        ax.set_xticks(x, [str(v) for v in x])
        ax.minorticks_off()
        ax.set_xlabel(xl)
        ax.legend(loc='upper right')
    for ax in axs[1:]:
        ax.tick_params(labelleft=True)
    fig.tight_layout(w_pad=0.8)
    fig.savefig(OUT/'fig_accuracy.pdf')


def fig_cost(plt, d):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.6))
    mem = d['memory']
    N = np.array([r['N'] for r in mem], float)
    Ns = np.geomspace(2e3, 1.2e5, 50)
    a1.loglog(Ns, 8*Ns**2/1e9, color='k', lw=0.8, label=r'one $N\times N$ matrix of doubles')
    a1.loglog([51200], [220], ls='none', marker='*', ms=7, color='k', label=r'PB3D, 128 lines $\times$ 400 points')
    a1.loglog(N, [r['peak_MiB']*2**20/1e9 for r in mem], marker='o', ms=3, color='#0072B2',
              label='wave form, whole process')
    a1.set_ylabel('peak memory of the process (GB)')
    a1.set_ylim(0.02, 800)
    a1.legend(loc='upper left')
    a1.set_xlabel(r'boundary points $N$')
    par = d['parallel']
    a2.loglog([r['N'] for r in par], [r['t_scratch'] for r in par], marker='o', ms=3, color='#0072B2',
              label='from scratch')
    a2.loglog([r['N'] for r in par if 't_carried' in r], [r['t_carried'] for r in par if 't_carried' in r],
              marker='s', ms=3, color='#D55E00', label='carried over from the previous level')
    a2.loglog(Ns, par[-1]['t_scratch']*(Ns/par[-1]['N'])**2, color='0.6', ls=':', lw=0.8, label=r'$\propto N^2$')
    a2.set_xlabel(r'boundary points $N$')
    a2.set_ylabel('time of the pair sums (s, one core)')
    a2.legend(loc='upper left')
    fig.tight_layout(w_pad=2.0)
    fig.savefig(OUT/'fig_cost.pdf')


def plot():
    plt = setup()
    d = json.loads(DATA.read_text())
    fig_grid(plt)
    fig_spectrum(plt, d)
    fig_matrices(plt, d)
    fig_carry(plt)
    fig_accuracy(plt, d)
    fig_cost(plt, d)


if __name__ == '__main__':
    {'data': data, 'plot': plot, 'single': add_single, 'compare': compare, 'memsplit': memsplit}.get(
        sys.argv[1], lambda: memory_run(int(sys.argv[2]), *sys.argv[3:4], *map(int, sys.argv[4:5])))()
