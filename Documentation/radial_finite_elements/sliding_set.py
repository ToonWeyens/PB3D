"""The sliding set of Fourier modes (X_style 2) with the coefficient matrices of PB3D runs.

Companion of pb3d_coefficients.py, which does the same for a fixed set of
Fourier modes. PB3D does everything up to the coefficient matrices PV_0..2 and
KV_0..2 and stores them in PB3D_out.h5; this script redoes only the last step,
what interp_V, fill_mat and set_BC do with them, for any table of fe_model.py.

What PB3D does with the sliding set, and what this script copies:

- The set of a position is the n_mod_X = W modes m = lo ... lo + W - 1, with
  lo = nint(linear interpolation of min_m_X between the equilibrium points),
  min_m_X = max(min_nm_X, nint(n q - (W - 1)/2)) at the equilibrium points
  (init_modes, setup_modes). Every mode keeps one slot:
  slot = mod(m - lo at the first point of the grid, W).
- interp_V interpolates every pair of modes (m_k, m_l) over the perturbation
  points where both modes are in the set, to the solution points where both
  are in the set: a not-a-knot cubic spline for 4 points or more, straight
  lines for 2 or 3, a copy for 1 (interp_V_spline). Outside the input points
  it extrapolates with value, first and half the second derivative at the end
  point (spline_utilities, calc_extrap), not with the cubic end polynomial.
- fill_mat builds the block of an energy point with the set of that point, and
  insert_block_mat drops every entry whose mode is not in the set of the
  unknown it lands on. In this script: an entry between unknown (a, m_k) and
  unknown (b, m_l) at energy point g is used if m_k and m_l are both in the set
  of g, m_k is a mode of unknown a and m_l a mode of unknown b.

The sets of the elements (rules of section 3.3 of the plan): the set of an
element is the set at its middle; its energy points use it.
  rule A: an unknown has the set of the first (innermost) element that uses it,
  rule B: an unknown has the modes that all elements using it have.
For elements of degree k with values at the points, rule A gives a shared point
the set of the element on its inner side; for cubic Hermite and B-splines the
same two definitions are used.

The coefficients can come from the sliding-set run itself (as PB3D would use
them, with the extrapolation of interp_V) or from a run with a fixed set that
contains every mode the sets need (no extrapolation). A free edge needs the
vacuum matrix; vac/res of a run is used, for the modes it was computed for.

    python -B sliding_set.py compare <run folder>      matrices of the differences against A.dat, B.dat of the run
    python -B sliding_set.py scan ...                  see the bottom of this file
"""
import sys
import time

import h5py
import numpy as np
import scipy.sparse as sp
from scipy.interpolate import CubicSpline
from scipy.sparse.linalg import eigsh

import fe_model as fm
from pb3d_coefficients import read_petsc

NAMES = ('PV_0', 'PV_1', 'PV_2', 'KV_0', 'KV_1', 'KV_2')
SYM = (True, False, True, True, False, True)
MIN_NM_X = 5                                    # X_vars.f90: min_nm_X


def fnint(x):
    """Fortran NINT: halves away from zero."""
    x = np.asarray(x, float)
    return (np.sign(x)*np.floor(np.abs(x) + 0.5)).astype(int)


def read(f, name):
    lim = f[name + '/lim'][:]
    n = len(lim)//2
    return f[name + '/var'][:].reshape([lim[n+i] - lim[i] + 1 for i in range(n)], order='F')


class Pair:
    """interp_V_spline for one pair of modes: the 6 coefficients over the input points where both modes exist."""

    def __init__(self, x, y):
        self.x, self.y, self.n = x, y, len(x)
        if self.n >= 4:
            self.s = CubicSpline(x, y, axis=0)    # not-a-knot, as PSPLINE with PB3D's default boundary conditions

    def __call__(self, xn):
        if self.n == 1:
            return np.repeat(self.y[:1], len(xn), axis=0)
        out = np.empty((len(xn), self.y.shape[1]), complex)
        lo, hi = xn < self.x[0], xn > self.x[-1]
        mid = ~(lo | hi)
        if self.n >= 4:
            out[mid] = self.s(xn[mid])
            for m, xb in ((lo, self.x[0]), (hi, self.x[-1])):
                if m.any():
                    d = (xn[m] - xb)[:, None]
                    out[m] = self.s(xb)[None] + self.s(xb, 1)[None]*d + 0.5*self.s(xb, 2)[None]*d**2
        else:                                       # straight lines, extrapolated with the slope of the end segment
            for c in range(self.y.shape[1]):
                out[mid, c] = np.interp(xn[mid], self.x, self.y[:, c].real) + 1j*np.interp(xn[mid], self.x, self.y[:, c].imag)
            for m, (i0, i1) in ((lo, (0, 1)), (hi, (-2, -1))):
                if m.any():
                    slope = (self.y[i1] - self.y[i0])/(self.x[i1] - self.x[i0])
                    xb, yb = (self.x[0], self.y[0]) if i0 == 0 else (self.x[-1], self.y[-1])
                    out[m] = yb[None] + slope[None]*(xn[m] - xb)[:, None]
        return out


class Run:
    """The coefficient matrices of one PB3D run, its sets of modes, its vacuum matrix."""

    def __init__(self, path):
        f = h5py.File(path, 'r')
        misc_X, misc_sol = read(f, 'in/misc_X'), read(f, 'in/misc_sol')
        self.path = path
        self.n, self.W, self.X_style = int(misc_X[0]), int(misc_X[1]), int(misc_X[7])
        self.p, self.edge = int(misc_sol[2]), 'free' if int(round(misc_sol[4])) == 2 else 'fixed'
        self.r_eq, self.q_eq = read(f, 'grid_eq/r_F'), read(f, 'eq_1/q_saf_FD')[:, 0]
        self.r_X, self.r_sol = read(f, 'grid_X/r_F'), read(f, 'grid_sol/r_F')
        self.a, self.L = self.r_sol[0], self.r_sol[-1] - self.r_sol[0]
        if self.X_style == 1:
            self.lo_eq = np.full(len(self.r_eq), int(misc_X[2]))
        else:
            assert self.n*self.q_eq[0] > 0
            self.lo_eq = np.maximum(MIN_NM_X, fnint(self.n*self.q_eq - (self.W - 1)*0.5))
        self.lam = read(f, 'sol_R_1/RE_sol_val')[0] if 'sol_R_1' in f else np.nan
        W = self.W
        V = np.zeros((len(self.r_X), 6, W, W), complex)
        for i, (name, sym) in enumerate(zip(NAMES, SYM)):
            v = (read(f, 'X_2_int_R_1/RE_' + name) + 1j*read(f, 'X_2_int_R_1/IM_' + name)).sum(axis=(0, 1))
            for m in range(1, W + 1):
                for k in range(m if sym else 1, W + 1):
                    if sym:                         # stored: the lower triangle, column by column
                        x = v[:, (m - 1)*W + k - (m - 1)*m//2 - 1]
                        V[:, i, k-1, m-1] = x
                        V[:, i, m-1, k-1] = x.conj() if k != m else x
                    else:
                        V[:, i, k-1, m-1] = v[:, (m - 1)*W + k - 1]
        self.V_X = V                                # (perturbation point, coefficient, slot, slot)
        self.lo_X = self.lo(self.r_X)
        self.vac = self.vac_m = None
        if 'vac' in f:
            self.vac = read(f, 'vac/RE_res') + 1j*read(f, 'vac/IM_res')
            self.vac_m = fnint(read(f, 'vac/sec_X'))
        self.pairs = {}

    def lo(self, r):
        """Lowest mode of the set at the positions r (setup_modes)."""
        return fnint(np.interp(r, self.r_eq, self.lo_eq*1.0))

    def lo_t(self, t):
        """Lowest mode of the set at positions 0 ... 1 of the solution range."""
        return self.lo(self.a + self.L*np.asarray(t))

    def pair(self, mk, ml):
        """interp_V of the pair (mk, ml) of the sliding-set run."""
        key = (mk, ml)
        if key not in self.pairs:
            W, lo = self.W, self.lo_X
            have = (lo <= min(mk, ml)) & (max(mk, ml) <= lo + W - 1)
            idx = np.nonzero(have)[0]
            if len(idx) == 0:
                raise ValueError('modes %d and %d are never together in the set of the perturbation grid' % (mk, ml))
            assert idx[-1] - idx[0] + 1 == len(idx), 'pair over a range that is not connected'
            sk, sl = (mk - lo[0]) % W, (ml - lo[0]) % W
            self.pairs[key] = Pair(self.r_X[idx], self.V_X[idx, :, sk, sl])
        return self.pairs[key]

    def coeff(self, t, modes):
        """The 6 coefficient matrices at positions t (0 ... 1) for the given modes, scaled to t: (len(t), 6, W', W')."""
        r = self.a + self.L*np.asarray(t, float)
        out = np.empty((len(r), 6, len(modes), len(modes)), complex)
        for a, mk in enumerate(modes):
            for b, ml in enumerate(modes):
                out[:, :, a, b] = self.pair(mk, ml)(r)
        scale = np.array([self.L, 1.0, 1/self.L, self.L, 1.0, 1/self.L])
        return out*scale[None, :, None, None]

    def vacuum(self, modes):
        """The vacuum matrix for the given modes, from vac/res of this run."""
        pos = {m: i for i, m in enumerate(self.vac_m)}
        missing = [m for m in modes if m not in pos]
        if missing:
            raise ValueError('no vacuum entries for modes %s in %s' % (missing, self.path))
        i = np.array([pos[m] for m in modes])
        return self.vac[np.ix_(i, i)]


# ----------------------------------------------------------------------------
# sets of modes per energy point and per unknown
# ----------------------------------------------------------------------------


class Sets:
    """The sliding set of W modes of a run's equilibrium, on its solution range (init_modes, setup_modes)."""

    def __init__(self, run, W):
        self.W = W
        self.r_eq, self.a, self.L = run.r_eq, run.a, run.L
        self.lo_eq = np.maximum(MIN_NM_X, fnint(run.n*run.q_eq - (W - 1)*0.5))
        self.m_ref = self.lo_t(0.0)

    def lo_t(self, t):
        return fnint(np.interp(self.a + self.L*np.asarray(t), self.r_eq, self.lo_eq*1.0))


def sets_of(run, W):
    """The sliding set of W modes, or for W = 0 the fixed set of the run."""
    S = Sets(run, run.W if W == 0 else W)
    if W == 0 or run.X_style == 1 and W == run.W:
        S.lo_eq = run.lo_eq.copy()
        S.m_ref = S.lo_t(0.0)
    return S


def window(lo, W, m_ref):
    """The set lo ... lo + W - 1 in slots: mode per slot (slot = mod(m - m_ref, W))."""
    out = np.empty(W, int)
    m = lo + np.arange(W)
    out[(m - m_ref) % W] = m
    return out


def sets_points(S, tab, W=None, m_ref=None):
    """PB3D's differences: every point and every energy point has the set at its own position."""
    W, m_ref = S.W, S.m_ref
    e = np.array([window(l, W, m_ref) for l in S.lo_t(tab.psi)])
    pts = np.linspace(0, 1, tab.n_dof)
    d = np.array([window(l, W, m_ref) for l in S.lo_t(pts)])
    return e, d


def element_of(tab, n_el):
    return np.minimum((np.array(tab.psi)*n_el).astype(int), n_el - 1)


def sets_elements(S, tab, n_el, rule, where=0.5):
    """One set per element, taken at the middle of the element; the sets of the unknowns by rule A or B."""
    W, m_ref = S.W, S.m_ref
    el = element_of(tab, n_el)
    lo_el = S.lo_t((np.arange(n_el) + where)/n_el)
    win_el = np.array([window(l, W, m_ref) for l in lo_el])
    e = win_el[el]
    d = np.full((tab.n_dof, W), -1)
    first = np.full(tab.n_dof, -1)
    for g in range(len(tab.psi)):
        for dof in tab.dof[g]:
            if first[dof] < 0:
                first[dof] = el[g]
                d[dof] = win_el[el[g]]
            elif rule == 'B':
                d[dof] = np.where(d[dof] == win_el[el[g]], d[dof], -1)
    assert (first >= 0).all()
    return e, d


# ----------------------------------------------------------------------------
# matrices and eigenvalue
# ----------------------------------------------------------------------------


def matrices(run, tab, e, d, W):
    """A and B in blocks of W slots per unknown, without the vacuum and without removing anything."""
    n = tab.n_dof
    nd = max(len(x) for x in tab.dof)
    w = max(int(np.max(x) - np.min(x)) for x in tab.dof)
    A = np.zeros((n, 2*w + 1, W, W), complex)
    B = np.zeros((n, 2*w + 1, W, W), complex)
    psi = np.array(tab.psi)
    lo_g = e.min(axis=1)
    chunk = max(64, int(2e8//(6*W*W*16)))           # energy points per chunk of coefficients (about 200 MB)
    for g in range(len(tab.psi)):
        if g % chunk == 0:                          # coefficients at the next energy points, grouped by set
            gc = np.arange(g, min(g + chunk, len(psi)))
            C = np.empty((len(gc), 6, W, W), complex)
            for lo in np.unique(lo_g[gc]):
                sel = np.nonzero(lo_g[gc] == lo)[0]
                C[sel] = run.coeff(psi[gc[sel]], e[gc[sel[0]]])      # modes in slot order
        P0, P1, P2, K0, K1, K2 = C[g - gc[0]]
        f, df, dof = tab.phi[g], tab.dphi[g], tab.dof[g]
        here = [(e[g] == d[x]) for x in dof]        # slots whose mode the energy point and the unknown both have
        for i in range(len(dof)):
            for j in range(len(dof)):
                c = f[i]*f[j], f[i]*df[j], df[i]*f[j], df[i]*df[j]
                if not any(c):
                    continue
                use = np.outer(here[i], here[j])
                A[dof[i], dof[j] - dof[i] + w] += tab.w[g]*(c[0]*P0 + c[1]*P1 + c[2]*P1.conj().T + c[3]*P2)*use
                B[dof[i], dof[j] - dof[i] + w] += tab.w[g]*(c[0]*K0 + c[1]*K1 + c[2]*K1.conj().T + c[3]*K2)*use
    indptr, indices, sel = [0], [], []
    for i in range(n):
        for o in range(2*w + 1):
            j = i + o - w
            if 0 <= j < n:
                indices.append(j), sel.append((i, o))
        indptr.append(len(indices))
    ii, oo = np.array(sel).T
    return [sp.bsr_matrix((blocks[ii, oo], np.array(indices), np.array(indptr)), shape=(n*W, n*W)).tocsr()
            for blocks in (A, B)]


def lowest(run, tab, e, d, W, sigma, vac_run=None, number=1):
    """Lowest eigenvalue: unknowns without a mode, held unknowns removed; vacuum added at a free edge."""
    n = tab.n_dof
    A, B = matrices(run, tab, e, d, W)
    keep = (d >= 0).ravel()
    held = list(tab.held) + (list(tab.held_edge) if run.edge == 'fixed' else [])
    for h in held:
        keep[h*W:(h + 1)*W] = False
    if run.edge == 'free':
        vr = vac_run or run
        ed = d[tab.edge]
        s = np.nonzero(ed >= 0)[0]
        V = vr.vacuum(ed[s])
        idx = tab.edge*W + s
        A = A.tolil()
        A[np.ix_(idx, idx)] = A[np.ix_(idx, idx)].toarray() + V
        A = A.tocsr()
    A, B = A[keep][:, keep].tocsc(), B[keep][:, keep].tocsc()
    herm = abs(A - A.conj().T).max()/abs(A).max()
    val = np.sort(eigsh(A, k=max(number, 1), M=B, sigma=sigma, which='LM', tol=1e-12)[0])
    return (val[0] if number == 1 else val[:number]), herm, int(keep.sum())


def compare(folder):
    """Matrices of the differences with PB3D's sets against A.dat and B.dat of a run, entry by entry."""
    run = Run(folder + '/PB3D_out.h5')
    N, p, W = len(run.r_sol), run.p, run.W
    tab = fm.table_fd(N, p, run.edge)
    e, d = sets_points(sets_of(run, W), tab)
    A, B = matrices(run, tab, e, d, W)
    inner = np.arange(p*W, (N - 1)*W if run.edge == 'fixed' else N*W)
    print('%s: N = %d, p = %d, W = %d, %s edge, sets from %d ... %d to %d ... %d'
          % (folder, N, p, W, run.edge, d[0].min(), d[0].max(), d[-1].min(), d[-1].max()))
    for name, mine, theirs in (('A', A, read_petsc(folder + '/A.dat')), ('B', B, read_petsc(folder + '/B.dat'))):
        diff = (theirs - mine).tocsr()[inner][:, inner].toarray()
        scale = abs(mine).max()
        if name == 'A' and run.edge == 'free':
            vac = diff[-W:, -W:].copy()
            diff[-W:, -W:] = 0
            mine_vac = run.vacuum(d[-1])
            print('  vacuum: difference in the last diagonal block against vac/res of the run: %.1e of its largest entry %.3e'
                  % (abs(vac - mine_vac).max()/abs(vac).max(), abs(vac).max()))
        nz_theirs = (abs(theirs.tocsr()[inner][:, inner]) > 0).sum()
        nz_mine = (abs(mine.tocsr()[inner][:, inner]) > 0).sum()
        print('  %s: largest entry %.3e; largest difference, relative to it: %.1e; non-zero entries PB3D %d, script %d'
              % (name, scale, abs(diff).max()/scale, nz_theirs, nz_mine))
    lam, herm, nk = lowest(run, tab, e, d, W, 3*run.lam if run.lam < 0 else -abs(run.lam))
    print('  eigenvalue PB3D %.12e, script %.12e, relative difference %.1e' % (run.lam, lam, (lam - run.lam)/abs(run.lam)))


SCHEMES = {
    'd3': ('differences p = 3 (PB3D)', lambda N, edge: (fm.table_fd(N, 3, edge), None, None)),
    'd2': ('differences p = 2 (PB3D)', lambda N, edge: (fm.table_fd(N, 2, edge), None, None)),
    'e2A': ('elements k = 2, rule A', lambda N, edge: (fm.table_elements((N - 1)//2, 2), (N - 1)//2, 'A')),
    'e2B': ('elements k = 2, rule B', lambda N, edge: (fm.table_elements((N - 1)//2, 2), (N - 1)//2, 'B')),
    'e3A': ('elements k = 3, rule A', lambda N, edge: (fm.table_elements((N - 1)//3, 3), (N - 1)//3, 'A')),
    'e3B': ('elements k = 3, rule B', lambda N, edge: (fm.table_elements((N - 1)//3, 3), (N - 1)//3, 'B')),
    'hA': ('cubic Hermite, rule A', lambda N, edge: (fm.table_hermite(N//2 - 1), N//2 - 1, 'A')),
    'hB': ('cubic Hermite, rule B', lambda N, edge: (fm.table_hermite(N//2 - 1), N//2 - 1, 'B')),
    'bA': ('cubic B-spline, rule A', lambda N, edge: (fm.table_bspline(N - 3, 3), N - 3, 'A')),
    'bB': ('cubic B-spline, rule B', lambda N, edge: (fm.table_bspline(N - 3, 3), N - 3, 'B')),
}


def scan(coeff_path, W, ref, Ns, keys, vac_path=None):
    """Lowest eigenvalue for the schemes `keys` at N unknowns per mode, sets of W modes (0: the run's fixed set)."""
    run = Run(coeff_path)
    vr = Run(vac_path) if vac_path else run
    S = sets_of(run, W)
    print('coefficients from %s (%s set of %d modes), sets of %s, vacuum from %s'
          % (coeff_path, 'fixed' if run.X_style == 1 else 'sliding', run.W,
             'the run (fixed)' if W == 0 else '%d modes' % W, vr.path if run.edge == 'free' else '-'), flush=True)
    out = {}
    for key in keys:
        label, make = SCHEMES[key]
        for N in Ns:
            t0 = time.time()
            tab, n_el, rule = make(N, run.edge)
            e, d = sets_points(S, tab) if n_el is None else sets_elements(S, tab, n_el, rule)
            lam, herm, nk = lowest(run, tab, e, d, S.W, 3*ref, vr)
            out[key, N] = lam
            print('%-28s N %5d  unknowns %6d  lambda %.10e  relative to the reference %+.3e  (A-A^H %.0e, %.0f s)'
                  % (label, tab.n_dof, nk, lam, (lam - ref)/abs(ref), herm, time.time() - t0), flush=True)
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'compare':
        for folder in sys.argv[2:]:
            compare(folder)
    elif sys.argv[1] == 'scan':
        # scan <PB3D_out.h5 with the coefficients> <W, 0 = fixed set of that run> <reference> <keys,...> <N> [<N> ...] [--vac=<PB3D_out.h5>]
        opts = dict(s[2:].split('=', 1) for s in sys.argv[2:] if s.startswith('--'))
        args = [s for s in sys.argv[2:] if not s.startswith('--')]
        scan(args[0], int(args[1]), float(args[2]), [int(n) for n in args[4:]], args[3].split(','), opts.get('vac'))
