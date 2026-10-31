"""Write TERPSICHORE's equilibrium file fort.18 from a VMEC netCDF wout file.

Usage: python wout2terps.py wout.nc fort.18 [--flip] [--nfp K]
Follows srcs/vmec2terps.f of TERPSICHORE 1.2 (EPFL, 2025): header, then per surface and
mode "m n*NFP rmnc zmns gmnc" + "sigma tau 0 0", then per surface "iotas mass mu0*pres -phips vp"
with the js = 1 profile values zero. gmnc is taken at the (m, n) of the R, Z spectrum from the
Nyquist set. Stellarator-symmetric files only. --flip changes the sign of phi (and so of iota
and phips) for files written with a negative toroidal flux. --nfp K declares K field periods for a file
written with nfp = 1 whose boundary has that periodicity (PB3D's rippled twins): only the modes with
n a multiple of K are kept, so TERPSICHORE solves one mode family instead of every n.
"""
import sys, numpy as np, netCDF4 as nc
mu0 = 4e-7*np.pi
src, dst = sys.argv[1], sys.argv[2]; flip = '--flip' in sys.argv
f = nc.Dataset(src); g = lambda k: np.array(f.variables[k][:])
ns, nfp, mpol, ntor, mnmax = (int(g(k)) for k in ('ns','nfp','mpol','ntor','mnmax'))
xm, xn = g('xm').astype(int), g('xn').astype(int)
rmnc, zmns, gmnc = g('rmnc'), g('zmns'), g('gmnc')          # (ns, mnmax), gmnc (ns, mnmax_nyq)
if '--nfp' in sys.argv:
    nfp = int(sys.argv[sys.argv.index('--nfp') + 1]); keep = (xn % nfp == 0)
    xm, xn, rmnc, zmns = xm[keep], xn[keep], rmnc[:, keep], zmns[:, keep]
    mnmax, ntor = int(keep.sum()), int(abs(xn).max()) // nfp
xmn, xnn = g('xm_nyq').astype(int), g('xn_nyq').astype(int)
idx = [int(np.where((xmn == m) & (xnn == n))[0][0]) for m, n in zip(xm, xn)]
gm = gmnc[:, idx]
iotas, mass, pres, phips, vp = g('iotas'), g('mass'), g('pres'), g('phips'), g('vp')
wb = float(g('wb')); gam = float(g('gamma')) if 'gamma' in f.variables else 0.0
lasym = int(g('lasym__logical__')) if 'lasym__logical__' in f.variables else 0
assert lasym == 0, 'stellarator-symmetric files only'
sgn = -1.0 if flip else 1.0
with open(dst, 'w') as o:
    o.write(' %22.12E %12.5E%12.5E%12.5E%4d%4d%4d%4d%4d%4d%4d%4d%1d\n' % (wb, gam, 1.0/nfp, 0.0, mnmax, ns, mpol, ntor, 1, 1, 100, 0, 0))
    for js in range(ns):
        for mn in range(mnmax):
            m, n = xm[mn], xn[mn]
            o.write(' %22.14E%22.14E%22.14E%22.14E%22.14E\n' % (m, n, rmnc[js, mn], zmns[js, mn], gm[js, mn] if js > 0 else 0.0))
            st = 1.0 if (m == 0 and n == 0) else 0.0
            o.write(' %22.14E%22.14E%22.14E%22.14E\n' % (st, st, 0.0, 0.0))
    for js in range(ns):
        if js == 0:
            o.write(' %22.14E%22.14E%22.14E%22.14E%22.14E\n' % (0, 0, 0, 0, 0))
        else:
            o.write(' %22.14E%22.14E%22.14E%22.14E%22.14E\n' % (sgn*iotas[js], mass[js], mu0*pres[js], -sgn*phips[js], vp[js]))
print('wrote', dst, ': ns', ns, 'nfp', nfp, 'mpol', mpol, 'ntor', ntor, 'mnmax', mnmax, 'wb %.6g' % wb,
      'iota edge %.5f' % (sgn*iotas[-1]), 'phips edge %.5g' % phips[-1], 'pres axis %.6g Pa' % pres[1], 'vp edge %.5g' % vp[-1],
      '| flip' if flip else '')
