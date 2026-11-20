#!/usr/bin/env python3
"""Check the HELENA shear S written by PB3D against the direct formula.

PB3D (eq_ops.calc_derived_q, HELENA branch) obtains the shear by splining
Theta^alpha = h^{alpha psi}/h^{psi psi} in theta and differentiating. On a
top-bottom symmetric equilibrium that spline must carry the odd boundary
condition (zero second derivative), since Theta^alpha is odd in theta.
The direct formula uses the metric derivatives that PB3D also writes:

    S = -(1/J) d/dtheta (h^{12}/h^{22})
      = -(D_theta h^{12}/h^{22} - D_theta h^{22} h^{12}/(h^{22})^2)/J.

With the even condition (the bug of June 2018, issue #4 step 2) S is forced
to zero at the outboard and inboard midplanes and deviates from the direct
formula in a boundary layer of a few theta intervals; in the interior the
two agree to about 1e-5. The test asserts agreement on the whole grid and a
non-zero outboard-midplane value.

Usage: check_shear_HEL.py PB3D_out.h5 [tolerance]
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
    tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4
    with h5py.File(fname, 'r') as f:
        h = read(f, 'eq_2/h_FD')           # (theta, zeta, r, component, d1, d2, d3)
        J = read(f, 'eq_2/jac_FD')         # (theta, zeta, r, d1, d2, d3)
        S = read(f, 'eq_2/S')              # (theta, zeta, r)
    h12, h22 = h[..., 1, 0, 0, 0], h[..., 3, 0, 0, 0]
    Dh12, Dh22 = h[..., 1, 0, 0, 1], h[..., 3, 0, 0, 1]
    S_direct = -(Dh12/h22 - Dh22*h12/h22**2)/J[..., 0, 0, 0]

    scale = np.abs(S).max()
    diff = np.abs(S - S_direct)/scale
    per_theta = diff.max(axis=(1, 2))
    edge = np.abs(S[0]).max()/scale
    edge_direct = np.abs(S_direct[0]).max()/scale
    print('max |S - S_direct|/max|S| = %.3e (theta index %d of %d); interior %.3e'
          % (diff.max(), per_theta.argmax(), len(per_theta), per_theta[3:-3].max()))
    print('outboard midplane: max_r |S|/max|S| = %.3e (direct formula: %.3e)'
          % (edge, edge_direct))
    ok = diff.max() < tol and edge > 0.5*edge_direct
    print('PASS' if ok else 'FAIL', '(tolerance %.1e)' % tol)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
