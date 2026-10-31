#!/usr/bin/env python3
"""Check the vacuum response matrix written by PB3D (group vac of PB3D_out.h5).

vac%res is the surface-energy Hessian as it enters the potential-energy
matrix A: X_e^dagger res X_e is the (doubled) vacuum magnetic energy of an
edge displacement X_e, which is the integral of |grad Phi|^2 over the vacuum
and therefore non-negative. The exact operator is Hermitian. The test asserts

  - Hermiticity, max |res - res^dagger| / max |res| below tol_herm,
  - positive semi-definiteness, min eigenvalue of (res + res^dagger)/2
    above -tol_psd * max eigenvalue.

Before issue #4 step 3 the matrix entered A with the wrong sign (every
eigenvalue negative) and an asymmetry of 2.7e-5.

Usage: check_vac_res.py PB3D_out.h5 [tol_psd] [tol_herm]
Exit status 0 and "PASS" on success, 1 and "FAIL" otherwise.
"""
import sys

import h5py
import numpy as np


def read(f, name):
    """Read a PB3D array with its lim (lower bounds, then upper bounds)."""
    g = f[name]
    lim = g['lim'][...]
    n = len(lim)//2
    shape = tuple(lim[n:] - lim[:n] + 1)
    return g['var'][...].reshape(shape, order='F')


def main():
    fname = sys.argv[1]
    tol_psd = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-8
    tol_herm = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-10
    with h5py.File(fname, 'r') as f:
        # HELENA runs write the group 'vac', VMEC runs one 'vac_R_<level>' per
        # Richardson level: take the last one
        groups = sorted(k for k in f.keys() if k == 'vac' or k.startswith('vac_R_'))
        if not groups:
            print('FAIL: no vacuum group in', fname)
            return 1
        grp = groups[-1]
        res = read(f, grp + '/RE_res') + 1j*read(f, grp + '/IM_res')
        sec = read(f, grp + '/sec_X').astype(int)
    print('group', grp)
    res = res[:, :, 0, 0] if res.ndim == 4 else res
    scale = np.abs(res).max()
    herm = np.abs(res - res.conj().T).max()/scale
    ev = np.linalg.eigvalsh(0.5*(res + res.conj().T))
    print('modes m =', sec.ravel().tolist())
    print('max|res - res^dagger|/max|res| = %.3e' % herm)
    print('eigenvalues of the Hermitian part: min %.4e, max %.4e, %d of %d negative'
          % (ev.min(), ev.max(), int((ev < 0).sum()), ev.size))
    ok = herm < tol_herm and ev.min() > -tol_psd*ev.max()
    print('PASS' if ok else 'FAIL', '(tolerances: psd %.1e, hermitian %.1e)' % (tol_psd, tol_herm))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
