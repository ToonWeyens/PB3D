"""One row per run from the W7-X check_vac runs with vac_diag_style 1 (indirect: row sum and
identities) and 2 (direct: corner term plus cell integral), and from the earlier cv7 series (the
code before this branch): the CHECKVAC lines of np*/run.log (identity residuals at the expected
free term, exact-solve errors, the condition estimate), the free-boundary eigenvalue, and the
per-point files check_vac_*.dat (max, mean, median, and the position of the largest value: i_par,
i_line, R and Z in units of R_0). The diagonal comparison is taken from the per-point files, not
from the CHECKVAC diagonal lines, which writo cut at 100 characters until 6 November 2026 (the
lines are split since; the identity residuals come on two lines per potential).
Usage: diag_table.py [--json out.json] <run-dir>..."""
import sys, re, os, glob, json
import numpy as np

def num(x):
    """float, or None for a token that writo cut at the line limit"""
    try:
        return float(x)
    except ValueError:
        return None

def stats(path, mask_zero=False, offset=0.0):
    try:
        a = np.loadtxt(path)
    except OSError:
        return None
    if a.ndim != 2 or a.shape[0] == 0:
        return None
    v = np.abs(a[:, 6] - offset)
    sel = (a[:, 6] != 0.0) if mask_zero else np.ones(len(v), bool)
    if not sel.any():
        return None
    k = int(np.flatnonzero(sel)[np.argmax(v[sel])])
    return dict(max=float(v[k]), mean=float(v[sel].mean()), med=float(np.median(v[sel])),
                i_par=int(a[k, 1]), i_line=int(a[k, 2]), R=float(np.hypot(a[k, 3], a[k, 4])),
                Z=float(a[k, 5]), zeta=float(np.degrees(np.arctan2(a[k, 4], a[k, 3])) % 72.0), n=int(sel.sum()))

args = sys.argv[1:]
json_out = None
if args and args[0] == '--json':
    json_out, args = args[1], args[2:]
rows = []
for p in args:
    name = os.path.basename(p.rstrip('/'))
    m = re.match(r'(\w+?)_w7x_na(\d+)_p(\d+)', name)
    series, na, npar = (m.group(1), int(m.group(2)), int(m.group(3))) if m else (name, 0, 0)
    for rdir in sorted(glob.glob(os.path.join(p, 'np*'))):
        log = os.path.join(rdir, 'run.log')
        if not os.path.isfile(log):
            continue
        info = dict(run=name, series=series, n_alpha=na, n_par=npar, np=os.path.basename(rdir), ident={}, solve={})
        for line in open(log, errors='replace'):
            m = re.search(r'CHECKVAC identity (\S+) n_bnd=(\d+) res\(c=([^)]*)\)=\s*(.*)', line)
            if m:  # two lines per potential since 6 November 2026: c = -4pi, -2pi, 0 and c = 2pi, 4pi
                vals = [num(x) for x in m.group(4).split()]
                if m.group(3).startswith('-4pi'):
                    info['ident'][m.group(1)] = vals
                else:
                    info['ident'].setdefault(m.group(1), []).extend(vals)
                info['n_bnd'] = int(m.group(2))
            m = re.search(r'CHECKVAC solve (\S+) n_bnd=\d+ min_dist=\S+ err_max=(\S+) err_2=(\S+)', line)
            if m:
                info['solve'][m.group(1)] = [float(m.group(2)), float(m.group(3))]
            m = re.search(r'CHECKVAC diagonal n_bnd=\d+ n_points=(\d+) vac_diag_style=(\d)', line)
            if m:
                info['n_avail'], info['style'] = int(m.group(1)), int(m.group(2))
            m = re.search(r'CHECKVAC condition .*cond1=\s*(\S+)', line)
            if m:
                info['cond'] = float(m.group(1))
            if 'ERROR' in line and 'error' not in info:
                info['error'] = line.strip()[:80]
        ev = os.path.join(rdir, 'PB3D_out_EV_R_1.txt')
        if os.path.isfile(ev):
            for l in open(ev):
                f = l.split()
                if f and f[0] == '1':
                    info['EV'] = float(f[1])
        for tag, fn, mz, off in [('diag', 'check_vac_diag_diff.dat', True, 0.0),
                                 ('facet', 'check_vac_diag_diff_facet.dat', True, 0.0),
                                 ('Ising', 'check_vac_diag_Ising.dat', True, 0.0),
                                 ('omega', 'check_vac_diag_omega.dat', True, 0.5),
                                 ('point', 'check_vac_solve_point_source.dat', False, 0.0),
                                 ('ring', 'check_vac_solve_axis_ring.dat', False, 0.0),
                                 ('ring_n', 'check_vac_solve_axis_ring_n.dat', False, 0.0),
                                 ('id_x', 'check_vac_identity_x.dat', False, 0.0),
                                 ('id_z', 'check_vac_identity_z.dat', False, 0.0)]:
            e = stats(os.path.join(rdir, fn), mz, off)
            if e:
                info['pt_' + tag] = e
        rows.append(info)

if json_out:
    json.dump(rows, open(json_out, 'w'), indent=1)

def g(d, k, i=None, fmt='%.2e'):
    try:
        v = d[k] if i is None else d[k][i]
        return fmt % v if isinstance(v, float) else str(v)
    except (KeyError, IndexError, TypeError):
        return '-'

print('# Table 1: identity residuals (interior harmonics at c=-4pi, exterior sources at c=0) and exact-solve errors (max norm), eigenvalue')
print('%-28s %-4s %-2s %6s | %-8s %-8s %-8s | %-8s %-8s %-8s | %-8s %-8s %-8s | %-5s %s' % (
    'run', 'np', 'st', 'n_bnd', 'id x', 'id z', 'id x2-y2', 'id point', 'id ring', 'id ring_n',
    'sv point', 'sv ring', 'sv ring_n', 'cond1', 'EV'))
for r in rows:
    print('%-28s %-4s %-2s %6s | %-8s %-8s %-8s | %-8s %-8s %-8s | %-8s %-8s %-8s | %-5s %s %s' % (
        r['run'][:28], r['np'], g(r, 'style'), g(r, 'n_bnd'),
        g(r['ident'], 'x', 0), g(r['ident'], 'z', 0), g(r['ident'], 'x2-y2', 0),
        g(r['ident'], 'point_source', 2), g(r['ident'], 'axis_ring', 2), g(r['ident'], 'axis_ring_n', 2),
        g(r['solve'], 'point_source', 0), g(r['solve'], 'axis_ring', 0), g(r['solve'], 'axis_ring_n', 0),
        g(r, 'cond', fmt='%.2f'), g(r, 'EV', fmt='%.7e'), r.get('error', '')))

print()
print('# Table 2: assembled diagonal of H against the direct one, |assembled - direct|/4pi over the points with a full cell')
print('#   style 1: row sum vs -2pi + cell integral (the quadrature error of the off-diagonal row); style 2: zero by construction')
print('%-28s %-2s %6s %7s | %-8s %-8s %-8s | %-8s %-8s | %-8s %-8s | %s' % (
    'run', 'st', 'n_bnd', 'n_cell', 'max', 'mean', 'median', 'facet max', 'Ising max', 'omega max', 'omega med',
    'largest at (i_par, i_line, R, Z, zeta in period)'))
for r in rows:
    if r['np'] == 'np1' and any(x['run'] == r['run'] and x['np'] == 'np4' for x in rows):
        continue
    e = r.get('pt_diag')
    loc = '(%d, %d, R=%.3f, Z=%+.3f, %.0f deg)' % (e['i_par'], e['i_line'], e['R'], e['Z'], e['zeta']) if e else '-'
    print('%-28s %-2s %6s %7s | %-8s %-8s %-8s | %-8s %-8s | %-8s %-8s | %s' % (
        r['run'][:28], g(r, 'style'), g(r, 'n_bnd'), g(r, 'n_avail'),
        g(e or {}, 'max'), g(e or {}, 'mean'), g(e or {}, 'med'),
        g(r.get('pt_facet', {}), 'max'), g(r.get('pt_Ising', {}), 'max'),
        g(r.get('pt_omega', {}), 'max'), g(r.get('pt_omega', {}), 'med'), loc))

print()
print('# Table 3: where the largest exact-solve error sits (relative to max|phi|): max (i_par, i_line, R, Z, zeta in period) median')
print('%-28s %-2s | %-46s | %-46s | %s' % ('run', 'st', 'axis ring', 'ring with cos n zeta', 'point source'))
for r in rows:
    if r['np'] == 'np1' and any(x['run'] == r['run'] and x['np'] == 'np4' for x in rows):
        continue
    def f(tag):
        e = r.get('pt_' + tag)
        return '%.2e (%d, %d, %.3f, %+.3f, %.0f) %.1e' % (e['max'], e['i_par'], e['i_line'], e['R'], e['Z'], e['zeta'], e['med']) if e else '-'
    print('%-28s %-2s | %-46s | %-46s | %s' % (r['run'][:28], g(r, 'style'), f('ring'), f('ring_n'), f('point')))
