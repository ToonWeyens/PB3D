"""One table row per run from the CHECKVAC lines that check_vac = .true. writes: identity residuals
(c = -4 pi for the interior harmonics x, z, x^2 - y^2; c = 0 for the exterior potentials) and the
exterior solve errors (max norm) for the point source, the axis ring and the ring with cos n zeta
weights, plus the 1-norm condition estimate when present. Each argument is a run directory (the
first np*/run.log, PB3D's standard output, is read) or such a log file. Usage: checkvac_table.py <run>..."""
import sys, re, os, glob
cols_id = [('x', 0), ('z', 0), ('x2-y2', 0), ('point_source', 2), ('axis_ring', 2), ('axis_ring_n', 2)]
cols_sv = ['point_source', 'axis_ring', 'axis_ring_n']
print('%-30s %6s | %-9s %-9s %-9s | %-9s %-9s %-9s | %-9s %-9s %-9s | %s' % ('run', 'n_bnd', 'id x', 'id z', 'id x2-y2',
      'id point', 'id ring', 'id ring_n', 'sv point', 'sv ring', 'sv ring_n', 'cond1'))
for p in sys.argv[1:]:
    f = p if os.path.isfile(p) else sorted(glob.glob(os.path.join(p, 'np*', 'run.log')))[0]
    name = os.path.basename(p.rstrip('/'))
    ident, solve, nb, cond = {}, {}, '', ''
    for line in open(f):
        m = re.search(r'CHECKVAC identity (\S+) n_bnd=(\d+) res\(c=[^)]*\)=\s*(.*)', line)
        if m:
            vals = m.group(3).split(); ident[m.group(1)] = vals; nb = m.group(2)
        m = re.search(r'CHECKVAC solve (\S+) n_bnd=\d+ min_dist=\S+ err_max=(\S+)', line)
        if m: solve[m.group(1)] = m.group(2)
        m = re.search(r'CHECKVAC condition .*cond1=\s*(\S+)', line)
        if m: cond = m.group(1)
    def g(d, k, i=None):
        try: return d[k][i] if i is not None else d[k]
        except (KeyError, IndexError): return '-'
    print('%-30s %6s | %-9s %-9s %-9s | %-9s %-9s %-9s | %-9s %-9s %-9s | %s' % (name[:30], nb,
          *[g(ident, k, i) for k, i in cols_id], *[g(solve, k) for k in cols_sv], cond))
