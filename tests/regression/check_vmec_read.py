#!/usr/bin/env python3
"""Check that PB3D stored a VMEC file as the file says (issue #16).

Compares the group 'in' of PB3D_out.h5 with the wout file of the same run:

 1. the file's mode lists against VMEC's conventions: each list in VMEC's
    order (m rising; n = 0 .. N for m = 0, n = -N .. N for m >= 1, n a
    multiple of N_fp), its length equal to the count its limits give, the
    Nyquist list containing every (m, n) of the main list, (0, 0) first, and
    the arrays as long as their list;
 2. PB3D's mode lists mn_V and mn_V_nyq against xm, xn and xm_nyq, xn_nyq,
    and PB3D's sizes (misc_in_V) against the file, exactly;
 3. the full-mesh coefficients rmnc, zmns (and rmns, zmnc for lasym) against
    PB3D's R_V and Z_V times R_0, and bsubsmns (and bsubsmnc) against
    component 1 of PB3D's B_V_sub times R_0 B_0, to rounding;
 4. the half-mesh coefficients lmns, gmnc, bsubumnc, bsubvmnc (and their
    lasym partners) against PB3D's full-mesh values, through an independent
    cubic spline of the half-mesh points placed at
      (a) VMEC's positions s = (j - 1.5)/(ns - 1), what the file means, and
      (b) s_j - 0.5/ns, where read_VMEC put them before the fix of #16
          (0.5/(ns (ns - 1)) too far out), to show what that error was;
 5. the Jacobian in real space, sum over the Nyquist list, on one field
    period of every written surface: PB3D's coefficients with PB3D's labels
    against the file's coefficients (spline (a) and (b)) with the file's
    labels, and the size of the modes that only the Nyquist list has.

Items 1 to 3 decide PASS or FAIL (exact, or relative 1e-13 for 3), and so
does item 4 (a) in the interior (relative 1e-12) when PB3D's radial
interpolation is cubic (norm_disc_prec_eq = 3, the default), the order of
the spline here; item 5 is a measurement. Differences are relative to the
largest coefficient (4) or the largest value on the surface (5). "interior"
leaves out the two surfaces at each end, where PSPLINE and scipy extrapolate
differently.

Usage: check_vmec_read.py PB3D_out.h5 wout.nc [out.json]
Exit status 0 and "PASS" on success, 1 and "FAIL" otherwise.
"""
import json
import sys

import h5py
import numpy as np
from scipy.interpolate import CubicSpline


def read_pb3d(f, name):
    """Read a PB3D array with its lim (lower bounds, then upper bounds)."""
    g = f[name]
    lim = g['lim'][...]
    n = len(lim)//2
    shape = tuple(lim[n:] - lim[:n] + 1)
    return g['var'][...].reshape(shape, order='F')


def read_wout(path):
    """Read the variables used here from a wout file (netCDF 3 or 4)."""
    names = ['ns', 'mpol', 'ntor', 'nfp', 'mnmax', 'mnmax_nyq', 'xm', 'xn',
             'xm_nyq', 'xn_nyq', 'rmnc', 'zmns', 'lmns', 'gmnc', 'bmnc',
             'bsubumnc', 'bsubvmnc', 'rmns', 'zmnc', 'lmnc', 'gmns', 'bmns',
             'bsubumns', 'bsubvmns', 'bsubsmns', 'bsubsmnc',
             'lasym__logical__', 'lasym']
    out = {}
    try:
        import netCDF4
        with netCDF4.Dataset(path) as d:
            for k in names:
                if k in d.variables:
                    out[k] = np.array(d.variables[k][...])
    except ImportError:
        if h5py.is_hdf5(path):                                                  # netCDF 4 (e.g. VMEC++)
            with h5py.File(path, 'r') as d:
                for k in names:
                    if k in d:
                        out[k] = np.array(d[k][...])
        else:                                                                   # netCDF 3 (e.g. VMEC2000)
            from scipy.io import netcdf_file
            with netcdf_file(path, 'r', mmap=False) as d:
                for k in names:
                    if k in d.variables:
                        out[k] = np.array(d.variables[k][...])
    lasym = out.get('lasym__logical__', out.get('lasym', 0))
    out['lasym'] = bool(int(np.ravel(lasym)[0]))
    for k in ['ns', 'mpol', 'ntor', 'nfp', 'mnmax', 'mnmax_nyq']:
        out[k] = int(np.ravel(out[k])[0])
    for k in ['xm', 'xn', 'xm_nyq', 'xn_nyq']:
        out[k] = np.rint(out[k]).astype(int)
    return out


def vmec_order(mpol, ntor, nfp):
    """(m, n) in VMEC's order: m rising, n = 0..N for m = 0, -N..N else."""
    mn = [(m, n*nfp) for m in range(mpol)
          for n in (range(0, ntor+1) if m == 0 else range(-ntor, ntor+1))]
    return np.array(mn, dtype=int).reshape(-1, 2)


def check_lists(w):
    """Item 1: the file's mode lists against VMEC's conventions."""
    fails = []
    nfp = w['nfp']
    main = np.column_stack([w['xm'], w['xn']])
    nyq = np.column_stack([w['xm_nyq'], w['xn_nyq']])
    for name, lst, nmax in [('main', main, w['mnmax']),
                            ('Nyquist', nyq, w['mnmax_nyq'])]:
        if len(lst) != nmax:
            fails.append('%s list has %d modes, mnmax says %d' % (name, len(lst), nmax))
        if np.any(lst[:, 1] % nfp):
            fails.append('%s list has n not a multiple of nfp = %d' % (name, nfp))
        mpol = lst[:, 0].max() + 1
        ntor = np.abs(lst[:, 1]).max()//nfp
        if name == 'main' and (mpol != w['mpol'] or ntor != w['ntor']):
            fails.append('main list spans m < %d, |n| <= %d nfp; file says mpol %d, ntor %d'
                         % (mpol, ntor, w['mpol'], w['ntor']))
        exp = vmec_order(mpol, ntor, nfp)
        if exp.shape != lst.shape or np.any(exp != lst):
            fails.append('%s list is not in VMEC order for m < %d, |n| <= %d nfp'
                         % (name, mpol, ntor))
        if tuple(lst[0]) != (0, 0):
            fails.append('%s list does not start with (0, 0)' % name)
    missing = set(map(tuple, main)) - set(map(tuple, nyq))
    if missing:
        fails.append('%d main modes missing from the Nyquist list, e.g. %s'
                     % (len(missing), sorted(missing)[:3]))
    arrays = [('rmnc', 'mnmax'), ('zmns', 'mnmax'), ('lmns', 'mnmax'),
              ('gmnc', 'mnmax_nyq'), ('bmnc', 'mnmax_nyq'),
              ('bsubumnc', 'mnmax_nyq'), ('bsubvmnc', 'mnmax_nyq')]
    if w['lasym']:
        arrays += [('rmns', 'mnmax'), ('zmnc', 'mnmax'), ('lmnc', 'mnmax'),
                   ('gmns', 'mnmax_nyq'), ('bmns', 'mnmax_nyq'),
                   ('bsubumns', 'mnmax_nyq'), ('bsubvmns', 'mnmax_nyq')]
    for a, n in arrays:
        if a not in w:
            fails.append('%s missing from the file' % a)
        elif w[a].shape != (w['ns'], w[n]):
            fails.append('%s has shape %s, expected (ns, %s) = (%d, %d)'
                         % (a, w[a].shape, n, w['ns'], w[n]))
    return fails


def half_to_full(arr, ns, offset):
    """Cubic spline of the half-mesh rows 2..ns at s_j - offset, on the full mesh."""
    s_full = np.arange(ns)/(ns - 1)
    s_half = s_full[1:] - offset
    return CubicSpline(s_half, arr[1:, :], axis=0, bc_type='not-a-knot',
                       extrapolate=True)(s_full)


def rel_diff(pb, ref, interior):
    """Max |pb - ref| / max |ref|, over all rows and over the interior rows."""
    scale = np.abs(ref).max()
    if scale == 0:
        return 0.0, 0.0
    d = np.abs(pb - ref)/scale
    return float(d.max()), float(d[interior].max())


def main():
    fname, wname = sys.argv[1], sys.argv[2]
    jname = sys.argv[3] if len(sys.argv) > 3 else None
    w = read_wout(wname)
    ns, nfp = w['ns'], w['nfp']
    res = {'file': wname, 'ns': ns, 'mnmax': w['mnmax'],
           'mnmax_nyq': w['mnmax_nyq'], 'lasym': w['lasym']}
    fails = []

    # 1. the file's mode lists
    fl = check_lists(w)
    fails += fl
    print('1. file mode lists: mnmax %d, mnmax_nyq %d (m < %d, |n| <= %d nfp), nfp %d, ns %d, lasym %s: %s'
          % (w['mnmax'], w['mnmax_nyq'], w['xm_nyq'].max() + 1,
             np.abs(w['xn_nyq']).max()//nfp, nfp, ns, w['lasym'],
             'ok' if not fl else '; '.join(fl)))

    with h5py.File(fname, 'r') as f:
        misc = read_pb3d(f, 'in/misc_in')
        misc_V = read_pb3d(f, 'in/misc_in_V')
        mn_V = read_pb3d(f, 'in/mn_V').astype(int)
        mn_nyq = read_pb3d(f, 'in/mn_V_nyq').astype(int)
        RZL = read_pb3d(f, 'in/RZL_V')
        jac = read_pb3d(f, 'in/jac_V')
        Bsub = read_pb3d(f, 'in/B_V_sub')
    R_0 = misc[3] if misc[12] == 1 else 1.0
    B_0 = misc[5] if misc[12] == 1 else 1.0
    n_r_in, n_r_eq = int(misc[14]), int(misc[15])

    # 2. PB3D's sizes and mode lists
    f2 = []
    pb_sizes = dict(lasym=misc_V[0] == 1, mnmax=int(misc_V[2]), mpol=int(misc_V[3]),
                    ntor=int(misc_V[4]), nfp=int(misc_V[5]), mnmax_nyq=int(misc_V[7]))
    for k, v in pb_sizes.items():
        if v != w[k]:
            f2.append('%s: PB3D %s, file %s' % (k, v, w[k]))
    if n_r_in != ns:
        f2.append('n_r_in %d, file ns %d' % (n_r_in, ns))
    if mn_V.shape != (w['mnmax'], 2) or np.any(mn_V != np.column_stack([w['xm'], w['xn']])):
        f2.append('mn_V differs from (xm, xn)')
    nyq_file = np.column_stack([w['xm_nyq'], w['xn_nyq']])
    if mn_nyq.shape != nyq_file.shape or np.any(mn_nyq != nyq_file):
        bad = (int(np.sum(np.any(mn_nyq != nyq_file, axis=1)))
               if mn_nyq.shape == nyq_file.shape else -1)
        f2.append('mn_V_nyq differs from (xm_nyq, xn_nyq) (%d modes)' % bad)
    mpol_tf = max(w['mpol'], w['xm_nyq'].max() + 1)
    ntor_tf = max(w['ntor'], np.abs(w['xn_nyq']).max()//nfp)
    if (int(misc_V[8]), int(misc_V[9])) != (mpol_tf, ntor_tf):
        f2.append('mpol_tf, ntor_tf = %d, %d; expected %d, %d'
                  % (misc_V[8], misc_V[9], mpol_tf, ntor_tf))
    fails += f2
    print('2. PB3D sizes and mode lists against the file: %s' % ('ok' if not f2 else '; '.join(f2)))

    # written radial range: the offset at which R matches the file best
    Rc = RZL[:, :, 0, 0].T*R_0
    errs = [np.abs(Rc - w['rmnc'][i0:i0 + n_r_eq, :]).max()
            for i0 in range(ns - n_r_eq + 1)]
    i0 = int(np.argmin(errs))
    rows = slice(i0, i0 + n_r_eq)
    jj = np.arange(i0, i0 + n_r_eq)
    interior = (jj >= 2) & (jj <= ns - 3)
    res['rows'] = [i0 + 1, i0 + n_r_eq]
    print('   written surfaces %d..%d of %d (1-based)' % (i0 + 1, i0 + n_r_eq, ns))

    # 3. full-mesh coefficients
    full = [('rmnc', RZL[:, :, 0, 0], R_0), ('zmns', RZL[:, :, 0, 3], R_0),
            ('bsubsmns', Bsub[:, :, 0, 1], R_0*B_0)]
    if w['lasym']:
        full += [('rmns', RZL[:, :, 0, 1], R_0), ('zmnc', RZL[:, :, 0, 2], R_0),
                 ('bsubsmnc', Bsub[:, :, 0, 0], R_0*B_0)]
    full = [x for x in full if x[0] in w]                                       # B_s is optional in a wout file
    f3 = []
    res['full'] = {}
    for name, pb, sc in full:
        d, _ = rel_diff(pb.T*sc, w[name][rows, :], interior)
        res['full'][name] = d
        if d > 1e-13:
            f3.append('%s %.2e' % (name, d))
    fails += f3
    print('3. full-mesh coefficients, max relative difference: %s'
          % ', '.join('%s %.1e' % kv for kv in res['full'].items()))

    # 4. half-mesh coefficients
    half = [('lmns', RZL[:, :, 0, 5], 1.0), ('gmnc', jac[:, :, 0, 0], R_0**3),
            ('bsubumnc', Bsub[:, :, 1, 0], R_0*B_0), ('bsubvmnc', Bsub[:, :, 2, 0], R_0*B_0)]
    if w['lasym']:
        half += [('lmnc', RZL[:, :, 0, 4], 1.0), ('gmns', jac[:, :, 0, 1], R_0**3),
                 ('bsubumns', Bsub[:, :, 1, 1], R_0*B_0),
                 ('bsubvmns', Bsub[:, :, 2, 1], R_0*B_0)]
    offsets = {'VMEC': 0.5/(ns - 1), 'old': 0.5/ns}
    cubic = int(misc[13]) == 3
    spl = {}
    res['half'] = {}
    f4 = []
    print('4. half-mesh coefficients, max relative difference (all surfaces / interior):')
    print('   %-9s %-25s %-25s' % ('', "(a) VMEC's s=(j-1.5)/(ns-1)", '(b) old s_j-0.5/ns'))
    for name, pb, sc in half:
        line = []
        for key, off in offsets.items():
            spl[name, key] = half_to_full(w[name], ns, off)
            d_all, d_int = rel_diff(pb.T*sc, spl[name, key][rows, :], interior)
            res['half'][name + '_' + key] = [d_all, d_int]
            line.append('%.1e / %.1e' % (d_all, d_int))
            if key == 'VMEC' and cubic and d_int > 1e-12:
                f4.append('%s %.2e' % (name, d_int))
        print('   %-9s %-25s %-25s' % (name, line[0], line[1]))
    if not cubic:
        print('   (a) not decisive: norm_disc_prec_eq = %d, the spline here is cubic' % misc[13])
    elif f4:
        f4 = ['half-mesh coefficients off VMEC\'s positions (interior): ' + ', '.join(f4)]
    fails += f4

    # 5. the Jacobian in real space on one field period
    th = 2*np.pi*np.arange(32)/32
    ze = 2*np.pi*np.arange(32)/(32*nfp)
    TH, ZE = [a.ravel() for a in np.meshgrid(th, ze, indexing='ij')]

    def trig(mn):
        ang = np.outer(mn[:, 0], TH) - np.outer(mn[:, 1], ZE)
        return np.cos(ang), np.sin(ang)

    c_pb, s_pb = trig(mn_nyq)
    c_f, s_f = trig(nyq_file)
    J_pb = (jac[:, :, 0, 0].T*R_0**3) @ c_pb
    if w['lasym']:
        J_pb += (jac[:, :, 0, 1].T*R_0**3) @ s_pb
    res['jac_real'] = {}
    print('5. Jacobian in real space, max |J_PB3D - J_file| / max |J_file| per surface'
          ' (all surfaces / interior):')
    for key in offsets:
        J_f = spl['gmnc', key][rows, :] @ c_f
        if w['lasym']:
            J_f += spl['gmns', key][rows, :] @ s_f
        d = np.abs(J_pb - J_f).max(axis=1)/np.abs(J_f).max(axis=1)
        res['jac_real'][key] = [float(d.max()), float(d[interior].max())]
        print('   against spline %s: %.1e / %.1e' % ('(a)' if key == 'VMEC' else '(b)',
                                                  d.max(), d[interior].max()))
    extra = ~np.array([tuple(mn) in set(map(tuple, np.column_stack([w['xm'], w['xn']])))
                       for mn in nyq_file])
    if extra.any():
        J_x = spl['gmnc', 'VMEC'][rows][:, extra] @ c_f[extra]
        J_t = spl['gmnc', 'VMEC'][rows] @ c_f
        size = float((np.abs(J_x).max(axis=1)/np.abs(J_t).max(axis=1)).max())
    else:
        size = 0.0
    res['nyq_only'] = {'modes': int(extra.sum()), 'max_rel_size': size}
    print('   modes only in the Nyquist list: %d, largest share of |J| on a surface %.1e'
          % (extra.sum(), size))

    ok = not fails
    res['pass'] = ok
    res['fails'] = fails
    if jname:
        with open(jname, 'w') as fj:
            json.dump(res, fj, indent=1)
    print('PASS' if ok else 'FAIL: ' + '; '.join(fails))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
