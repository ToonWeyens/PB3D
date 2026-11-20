"""Finite differences against finite elements with the coefficient matrices of a PB3D run.

This is not PB3D in Python. PB3D does everything up to the coefficient
matrices PV_0..2 and KV_0..2 (equilibrium, field-aligned grid, integrals along
the field line, vacuum) and writes them on the perturbation grid to
PB3D_out.h5. This script redoes only the last step, what fill_mat and set_BC
do: it reads the coefficient matrices, interpolates them with a cubic spline
(as interp_V does) to wherever a scheme needs them, builds the two matrices
from a table of fe_model.py (the differences of PB3D, or elements) and finds
the lowest eigenvalue.

The eigenvalue of the differences must equal the one PB3D printed for the same
number of radial points; that checks the reading, the interpolation and the
table. The eigenvalue of the elements is then what an implementation in PB3D
would give for this case.

Only for runs with a fixed set of Fourier modes (X_style 1). The vacuum matrix
of a free-boundary run is not in the output file. It is taken from the
matrices that PB3D writes with the options
`-eps_view_mat0 binary:A.dat -eps_view_mat1 binary:B.dat`:

    python -B pb3d_coefficients.py <PB3D_out.h5> --compare=<A.dat>,<B.dat>,<N>,<p> [--vac=<file.npy to write>]

compares the matrices of this script with those of PB3D entry by entry and
writes the vacuum matrix. Then

    python -B pb3d_coefficients.py <PB3D_out.h5> <reference eigenvalue or 0> <N> [<N> ...] [--vac=<file.npy>] [--only=<part of a label>]
"""
import sys
import time

import h5py
import numpy as np
import scipy.sparse as sp
from scipy.interpolate import CubicSpline
from scipy.sparse.linalg import eigsh

import fe_model as fm


def load(path):
    """The problem of a PB3D run, on the radial range of its solution grid mapped to 0 ... 1."""
    f = h5py.File(path, 'r')

    def arr(name):
        lim = f[name + '/lim'][:]
        n = len(lim)//2
        return f[name + '/var'][:].reshape([lim[n+i] - lim[i] + 1 for i in range(n)], order='F')

    misc_X, misc_sol = arr('in/misc_X'), arr('in/misc_sol')
    M = int(misc_X[1])
    r_X, r_sol = arr('grid_X/r_F'), arr('grid_sol/r_F')
    a, L = r_sol[0], r_sol[-1] - r_sol[0]

    def dense(name, sym):
        v = (arr('X_2_int_R_1/RE_' + name) + 1j*arr('X_2_int_R_1/IM_' + name)).sum(axis=(0, 1))    # (radial point, pair)
        V = np.zeros((v.shape[0], M, M), complex)
        for m in range(1, M + 1):
            for k in range(m if sym else 1, M + 1):
                if sym:                                # stored: the lower triangle, column by column
                    x = v[:, (m - 1)*M + k - (m - 1)*m//2 - 1]
                    V[:, k-1, m-1] = x
                    V[:, m-1, k-1] = x.conj() if k != m else x
                else:
                    V[:, k-1, m-1] = v[:, (m - 1)*M + k - 1]
        return V

    stack = np.stack([dense('PV_0', True), dense('PV_1', False), dense('PV_2', True),
                      dense('KV_0', True), dense('KV_1', False), dense('KV_2', True)], axis=1)
    spl = CubicSpline(r_X, stack, axis=0)             # not-a-knot at both ends, as PB3D's spline
    scale = np.array([L, 1.0, 1/L, L, 1.0, 1/L])[:, None, None]    # d/dr = (1/L) d/dt and dr = L dt

    def coeff(t):
        return tuple(spl(a + L*t)*scale)

    edge = 'free' if int(round(misc_sol[4])) == 2 else 'fixed'
    prob = fm.Problem(M, None, None, None, edge, path)
    prob.coeff = coeff
    info = dict(lam=arr('sol_R_1/RE_sol_val')[0], N=len(r_sol), p=int(misc_sol[2]), n=int(misc_X[0]),
                m=(int(misc_X[2]), int(misc_X[3])), n_X=len(r_X), edge=edge)
    return prob, info


def matrices(prob, tab):
    """A and B of a table for all unknowns, without the vacuum term and without held unknowns removed."""
    M, n = prob.M, tab.n_dof
    w = max(int(np.max(d) - np.min(d)) for d in tab.dof)
    A = np.zeros((n, 2*w + 1, M, M), complex)
    B = np.zeros((n, 2*w + 1, M, M), complex)
    for g in range(len(tab.psi)):
        P0, P1, P2, K0, K1, K2 = prob.coeff(tab.psi[g])
        P1h, K1h = P1.conj().T, K1.conj().T
        f, d, dof = tab.phi[g], tab.dphi[g], tab.dof[g]
        for i in range(len(dof)):
            for j in range(len(dof)):
                c = f[i]*f[j], f[i]*d[j], d[i]*f[j], d[i]*d[j]
                if not any(c):
                    continue
                A[dof[i], dof[j] - dof[i] + w] += tab.w[g]*(c[0]*P0 + c[1]*P1 + c[2]*P1h + c[3]*P2)
                B[dof[i], dof[j] - dof[i] + w] += tab.w[g]*(c[0]*K0 + c[1]*K1 + c[2]*K1h + c[3]*K2)
    indptr, indices, sel = [0], [], []
    for i in range(n):
        for o in range(2*w + 1):
            j = i + o - w
            if 0 <= j < n:
                indices.append(j), sel.append((i, o))
        indptr.append(len(indices))
    ii, oo = np.array(sel).T
    return [sp.bsr_matrix((blocks[ii, oo], np.array(indices), np.array(indptr)), shape=(n*M, n*M)).tocsr() for blocks in (A, B)]


def lowest(prob, tab, sigma, vac=None):
    """Lowest eigenvalue of the matrices of a table."""
    M, n = prob.M, tab.n_dof
    A, B = matrices(prob, tab)
    keep = np.ones(n*M, bool)
    held = list(tab.held) + (list(tab.held_edge) if prob.edge == 'fixed' else [])
    for h in held:
        keep[h*M:(h + 1)*M] = False
    if prob.edge == 'free':
        e = tab.edge*M
        A = A.tolil()
        A[e:e+M, e:e+M] = A[e:e+M, e:e+M].toarray() + vac
        A = A.tocsr()
    A, B = A[keep][:, keep].tocsc(), B[keep][:, keep].tocsc()
    herm = abs(A - A.conj().T).max()/abs(A).max()
    val = eigsh(A, k=1, M=B, sigma=sigma, which='LM', tol=1e-12)[0]
    return val[0], herm


def read_petsc(path):
    """A matrix written by PETSc in binary form (complex scalars)."""
    b = open(path, 'rb').read()
    cid, rows, cols, nnz = [int(x) for x in np.frombuffer(b, dtype='>i4', count=4)]
    off = 16
    row_len = np.frombuffer(b, dtype='>i4', count=rows, offset=off)
    off += 4*rows
    col = np.frombuffer(b, dtype='>i4', count=nnz, offset=off)
    off += 4*nnz
    val = np.frombuffer(b, dtype='>f8', count=2*nnz, offset=off)
    return sp.csr_matrix((val[0::2] + 1j*val[1::2], col, np.r_[0, np.cumsum(row_len)]), shape=(rows, cols))


def compare(prob, files, vac_file):
    """Matrices of this script against those of PB3D, entry by entry; the vacuum matrix is what is left at the edge."""
    fa, fb, N, p = files.split(',')
    N, p, M = int(N), int(p), prob.M
    tab = fm.table_fd(N, p, prob.edge)
    A, B = matrices(prob, tab)
    inner = np.arange(p*M, (N - 1)*M if prob.edge == 'fixed' else N*M)       # rows and columns that PB3D does not overwrite
    for name, mine, theirs in (('A', A, read_petsc(fa)), ('B', B, read_petsc(fb))):
        d = (theirs - mine).tocsr()[inner][:, inner].toarray()
        scale = abs(mine).max()
        if name == 'A' and prob.edge == 'free':
            vac = d[-M:, -M:].copy()
            d[-M:, -M:] = 0
            print('vacuum matrix (the difference in the last diagonal block): largest entry %.3e, minus its Hermitian '
                  'transpose %.1e' % (abs(vac).max(), abs(vac - vac.conj().T).max()))
            if vac_file:
                np.save(vac_file, vac)
        print('%s: %d x %d, largest entry %.3e; largest difference between PB3D and this script, relative to it: %.1e'
              % (name, len(inner), len(inner), scale, abs(d).max()/scale))


if __name__ == '__main__':
    opts = dict(s[2:].split('=', 1) for s in sys.argv[1:] if s.startswith('--'))
    args = [s for s in sys.argv[1:] if not s.startswith('--')]
    prob, info = load(args[0])
    print('PB3D run: %d radial points, order %d, %s edge, n = %d, m = %d ... %d, %d perturbation points; its eigenvalue %.10e'
          % (info['N'], info['p'], info['edge'], info['n'], *info['m'], info['n_X'], info['lam']))
    if 'compare' in opts:
        compare(prob, opts['compare'], opts.get('vac'))
        sys.exit()
    vac = np.load(opts['vac']) if 'vac' in opts else None
    ref = float(args[1]) or info['lam']
    print('reference eigenvalue for the relative errors: %.8e' % ref)
    edge = prob.edge
    schemes = [('differences p = 2 (PB3D)', lambda N: fm.table_fd(N, 2, edge)),
               ('differences p = 3 (PB3D)', lambda N: fm.table_fd(N, 3, edge)),
               ('elements k = 1', lambda N: fm.table_elements(N - 1, 1)),
               ('elements k = 2', lambda N: fm.table_elements((N - 1)//2, 2)),
               ('elements k = 3', lambda N: fm.table_elements((N - 1)//3, 3)),
               ('elements k = 3, 6 Gauss points', lambda N: fm.table_elements((N - 1)//3, 3, n_q=6)),
               ('cubic Hermite', lambda N: fm.table_hermite(N//2 - 1)),
               ('cubic B-spline', lambda N: fm.table_bspline(N - 3, 3))]
    for label, make in schemes:
        if 'only' in opts and opts['only'] not in label:
            continue
        for N in [int(s) for s in args[2:]]:
            t0 = time.time()
            tab = make(N)
            lam, herm = lowest(prob, tab, 3*ref, vac)
            print('%-32s unknowns %5d  lambda = %.10e  relative error %+.2e  (A minus its Hermitian transpose: %.0e; %.0f s)'
                  % (label, tab.n_dof, lam, (lam - ref)/abs(ref), herm, time.time() - t0), flush=True)
