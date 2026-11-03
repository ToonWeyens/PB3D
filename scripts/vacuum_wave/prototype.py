"""Prototype: the 3-D vacuum with wave amplitudes as unknowns (issue #4).

Exploration only, nothing decided. PB3D's field-line vacuum (`vac%style` 1)
solves for the potential at every boundary point, with N x N matrices G and H
(N = field lines x parallel points). Here the same rows of G and H are
projected on waves exp(i(n phi - m theta)) from both sides (Galerkin), so the
unknowns are the amplitudes of M' waves and every matrix entry is a sum over
pairs of points. A Richardson level (parallel points 2 n_par - 1) then carries
the pair sums over as 1/4 of the old ones plus the pairs with a new point, and
no N x N matrix is stored.

Test case: circular torus R0 = 3, a = 1, field lines phi = alpha + q theta,
exact exterior harmonics sqrt(cosh eta - cos u) P^n_{|m|-1/2}(cosh eta)
exp(i(n phi - m u)) in toroidal coordinates (eta, u, phi). Given their normal
derivative on the torus, both methods solve for the potential; the error is
measured against the exact potential and against the exact vacuum energy
-int conj(q) Phi dS.

Success criteria, fixed before the runs (3 November 2026):
1. the carried-over sums equal the ones computed from scratch to 1e-12;
2. with enough waves, the wave form is as accurate as the point-value solve
   on the same grid (error at most 1.5 times larger);
3. the number of waves stays a small multiple of the plasma's mode window;
4. the peak memory of the wave form stays far below one N x N matrix.

Conventions: n is the normal out of the plasma (into the vacuum), q = dPhi/dn,
K(x, y) = (x - y).n_y/|x - y|^3, G(x, y) = 1/|x - y|. Exterior boundary
equation: 2 pi Phi(x) - PV int Phi K dS = - int G q dS. The self terms come
from the interior identities, as in calc_GH_1: the constant for K, the linear
potentials x, y, z for G.

Usage: .venv/bin/python scripts/vacuum_wave/prototype.py
"""
import resource
import sys
import time

import numpy as np

R0, A = 3.0, 1.0
C = np.sqrt(R0**2 - A**2)                   # focal radius of toroidal coordinates
Q_EDGE = 1.04296875                         # q n_alpha = k + 1/2, 3/4, 3/8 for 128, 64, 32 lines: no line end meets a line start
N_TOR = 10
M_TEST = (6, 11, 15)                        # poloidal mode numbers of the exact harmonics


# --- exact exterior harmonics --------------------------------------------------

def legendre_p(nu, n, eta, n_t=512):
    """P^n_nu(cosh eta) up to an eta-independent factor (Laplace's integral)."""
    eta = np.asarray(eta)
    t = np.linspace(0, np.pi, n_t)
    wt = np.full(n_t, np.pi/(n_t-1))
    wt[[0, -1]] *= 0.5                      # trapezoid on an even periodic integrand: spectral
    base = np.cosh(eta)[..., None] + np.sinh(eta)[..., None]*np.cos(t)
    return (base**nu*np.cos(n*t)*wt).sum(-1)


def exact_phi(R, Z, phi, m, n=N_TOR):
    """Decaying harmonic outside the torus, regular on the z axis."""
    eta = 0.5*np.log(((R+C)**2+Z**2)/((R-C)**2+Z**2))
    u = np.arctan2(2*C*Z, R**2+Z**2-C**2)
    return (np.sqrt(np.cosh(eta)-np.cos(u))*legendre_p(abs(m)-0.5, n, eta)
            * np.exp(1j*(n*phi-m*u)))


def exact_on_torus(theta, phi, m, h=1e-3):
    """Phi and q = dPhi/dr on r = A (fourth-order differences in r)."""
    def at(r):
        return exact_phi(R0+r*np.cos(theta), r*np.sin(theta), phi, m)
    q = (-at(A+2*h)+8*at(A+h)-8*at(A-h)+at(A-2*h))/(12*h)
    return at(A), q


def check_harmonic(m, h=1e-2):
    """Relative Laplacian of exact_phi at a few exterior points (should be ~h^4)."""
    out = []
    for R, Z in ((4.5, 0.3), (1.6, -0.8), (3.2, 1.7)):
        f = lambda dR, dZ: exact_phi(R+dR, Z+dZ, 0.0, m)
        d2 = lambda a, b: (-a(2)+16*a(1)-30*a(0)+16*a(-1)-a(-2))/(12*h**2)
        fRR = d2(lambda k: f(k*h, 0), None)
        fZZ = d2(lambda k: f(0, k*h), None)
        fR = (-f(2*h, 0)+8*f(h, 0)-8*f(-h, 0)+f(-2*h, 0))/(12*h)
        lap = fRR+fR/R+fZZ-N_TOR**2/R**2*f(0, 0)
        out.append(abs(lap)/(abs(fRR)+abs(fZZ)+abs(N_TOR**2/R**2*f(0, 0))))
    return max(out)


def exact_energy(m, n_th=1024):
    """-int conj(q) Phi dS on a periodic (theta, phi) grid: spectral."""
    th = np.linspace(-np.pi, np.pi, n_th, endpoint=False)
    ph, q = exact_on_torus(th, 0.0*th, m)
    return -(np.conj(q)*ph*A*(R0+A*np.cos(th))).sum()*(2*np.pi/n_th)*2*np.pi


# --- field-line grid ------------------------------------------------------------

class Grid:
    """n_alpha field lines, n_par points per line over one poloidal turn."""

    def __init__(self, n_alpha, n_par):
        self.n_alpha, self.n_par = n_alpha, n_par
        alpha = 2*np.pi*np.arange(n_alpha)/n_alpha
        theta = np.linspace(-np.pi, np.pi, n_par)
        al, th = np.meshgrid(alpha, theta, indexing='ij')      # line-major order
        self.theta, self.phi = th.ravel(), (al+Q_EDGE*th).ravel()
        R = R0+A*np.cos(self.theta)
        self.x = np.stack([R*np.cos(self.phi), R*np.sin(self.phi), A*np.sin(self.theta)], 1)
        self.nrm = np.stack([np.cos(self.theta)*np.cos(self.phi),
                             np.cos(self.theta)*np.sin(self.phi), np.sin(self.theta)], 1)
        w_th = np.full(n_par, 2*np.pi/(n_par-1))
        w_th[[0, -1]] *= 0.5
        self.w = (np.tile(w_th, n_alpha)*(2*np.pi/n_alpha)*A*R)
        self.n = n_alpha*n_par

    def old_points(self):
        """Indices of the previous Richardson level's points (even along each line)."""
        return (np.arange(self.n_alpha)[:, None]*self.n_par
                + np.arange(0, self.n_par, 2)[None, :]).ravel()

    def waves(self, m_list):
        return np.exp(1j*(N_TOR*self.phi[:, None]-np.outer(self.theta, m_list)))

    def data(self):
        """Exact (Phi, q) of the test harmonics at the points, one column each."""
        cols = [exact_on_torus(self.theta, self.phi, m) for m in M_TEST]
        return np.stack([c[0] for c in cols], 1), np.stack([c[1] for c in cols], 1)


def kernel_rows(g, I, J):
    """w_j K(x_i, x_j) and w_j G(x_i, x_j) for i in I, j in J; zero for i = j."""
    d = [g.x[I, k, None]-g.x[None, J, k] for k in range(3)]
    r2 = d[0]**2+d[1]**2+d[2]**2
    r2[I[:, None] == J[None, :]] = np.inf
    rinv = 1/np.sqrt(r2)
    dot = d[0]*g.nrm[J, 0]+d[1]*g.nrm[J, 1]+d[2]*g.nrm[J, 2]
    return dot*rinv**3*g.w[J], rinv*g.w[J]


# --- the two methods ------------------------------------------------------------

def self_terms(g, S, Kx, Gn):
    """Self-cell corrections c (of K) and gs (of G) from the interior identities."""
    c = -2*np.pi-S                                         # sum_j K_ij = -2 pi for Phi = 1
    gs = (g.nrm*(2*np.pi*g.x+Kx+c[:, None]*g.x-Gn)).sum(1) # Phi = x, y, z; |n| = 1
    return c, gs


def solve_points(g, qd):
    """Point values: N x N matrices and a dense solve, as PB3D does now."""
    Kw, Gw = kernel_rows(g, np.arange(g.n), np.arange(g.n))
    S, Kx, Gn = Kw.sum(1), Kw@g.x, Gw@g.nrm
    c, gs = self_terms(g, S, Kx, Gn)
    M = -Kw
    M[np.diag_indices(g.n)] = 4*np.pi+S
    return np.linalg.solve(M, -(Gw@qd+gs[:, None]*qd))


class WaveSums:
    """The pair sums of the wave form: everything a Richardson level keeps."""

    def __init__(self, g, m_list, qd):
        self.g, self.m_list, self.qd = g, np.asarray(m_list), qd
        self.E = g.waves(self.m_list)
        mw = len(self.m_list)
        self.A_off = np.zeros((mw, mw), complex)                # -sum_{i!=j} w_i e_l*(x_i) w_j K_ij e_k(x_j)
        self.B_off = np.zeros((mw, qd.shape[1]), complex)      # -sum_{i!=j} w_i e_l*(x_i) w_j G_ij q(x_j)
        self.S, self.Kx, self.Gn = np.zeros(g.n), np.zeros((g.n, 3)), np.zeros((g.n, 3))

    def add_pairs(self, I, J, chunk=256):
        """Add the pairs with observation point in I and source point in J."""
        g, E = self.g, self.E
        Er = np.concatenate([E.real, E.imag], 1)
        Qr = np.concatenate([self.qd.real, self.qd.imag], 1)
        mw, nq = E.shape[1], self.qd.shape[1]
        for s in range(0, len(I), chunk):
            Ic = I[s:s+chunk]
            Kw, Gw = kernel_rows(g, Ic, J)
            self.S[Ic] += Kw.sum(1)
            self.Kx[Ic] += Kw@g.x[J]
            self.Gn[Ic] += Gw@g.nrm[J]
            vK, vG = Kw@Er[J], Gw@Qr[J]
            test = g.w[Ic, None]*np.conj(E[Ic])
            self.A_off -= test.T@(vK[:, :mw]+1j*vK[:, mw:])
            self.B_off -= test.T@(vG[:, :nq]+1j*vG[:, nq:])

    def solve(self):
        g, E = self.g, self.E
        c, gs = self_terms(g, self.S, self.Kx, self.Gn)
        test = g.w[:, None]*np.conj(E)
        Amat = self.A_off+(test*(4*np.pi+self.S)[:, None]).T@E
        rhs = self.B_off-test.T@(gs[:, None]*self.qd)
        return E@np.linalg.solve(Amat, rhs)


def wave_sums_from_scratch(g, m_list, qd):
    ws = WaveSums(g, m_list, qd)
    every = np.arange(g.n)
    ws.add_pairs(every, every)
    return ws


def wave_sums_carried(prev, g, m_list, qd):
    """Next Richardson level from the previous level's sums: old pairs x 1/4."""
    ws = WaveSums(g, m_list, qd)
    old = g.old_points()
    new = np.setdiff1d(np.arange(g.n), old)
    ws.A_off, ws.B_off = 0.25*prev.A_off, 0.25*prev.B_off
    ws.S[old], ws.Kx[old], ws.Gn[old] = 0.5*prev.S, 0.5*prev.Kx, 0.5*prev.Gn
    ws.add_pairs(old, new)
    ws.add_pairs(new, np.arange(g.n))
    return ws


# --- measurements ---------------------------------------------------------------

def errors(g, ph, ph_ex, qd, e_ex):
    """Relative L2 error of Phi and relative error of the energy, per harmonic."""
    w = g.w[:, None]
    e_phi = np.sqrt((w*abs(ph-ph_ex)**2).sum(0)/(w*abs(ph_ex)**2).sum(0))
    energy = -(w*np.conj(qd)*ph).sum(0)
    return e_phi, abs(energy-e_ex)/abs(e_ex)


def window(k):
    return np.arange(min(M_TEST)-k, max(M_TEST)+k+1)


def peak_mb():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss/2**20 if sys.platform == 'darwin' else rss/2**10


def main():
    print(f'torus R0 = {R0}, a = {A}, q = {Q_EDGE}, n = {N_TOR}, test m = {M_TEST}')
    print('relative Laplacian of the exact harmonics:',
          ', '.join(f'm={m}: {check_harmonic(m):.1e}' for m in M_TEST))
    e_ex = np.array([exact_energy(m) for m in M_TEST])
    print('exact energy -int conj(q) Phi dS:', ', '.join(f'{e.real:.6e}{e.imag:+.0e}i' for e in e_ex))
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'

    n_alpha = 64
    if mode in ('all', 'waves'):
        print(f'\n1. waves needed: {n_alpha} lines x 201 points, sidebands k around m = '
              f'{min(M_TEST)}..{max(M_TEST)}')
        g = Grid(n_alpha, 201)
        ph_ex, qd = g.data()
        ws = wave_sums_from_scratch(g, window(20), qd)
        for k in (0, 1, 2, 4, 6, 8, 12, 16, 20):
            sel = np.isin(ws.m_list, window(k))
            sub = WaveSums(g, ws.m_list[sel], qd)
            sub.A_off, sub.B_off = ws.A_off[np.ix_(sel, sel)], ws.B_off[sel]
            sub.S, sub.Kx, sub.Gn = ws.S, ws.Kx, ws.Gn
            e_phi, e_en = errors(g, sub.solve(), ph_ex, qd, e_ex)
            print(f'   k = {k:2d}, {sel.sum():2d} waves: Phi error ' + ' '.join(f'{e:.2e}' for e in e_phi)
                  + ' | energy error ' + ' '.join(f'{e:.2e}' for e in e_en))

    if mode in ('all', 'grid'):
        print(f'\n2. grid: {n_alpha} lines, parallel points doubled, waves k = 12')
        prev = None
        for n_par in (51, 101, 201, 401):
            g = Grid(n_alpha, n_par)
            ph_ex, qd = g.data()
            t0 = time.time()
            ws = wave_sums_from_scratch(g, window(12), qd)
            t_scratch = time.time()-t0
            e_phi, e_en = errors(g, ws.solve(), ph_ex, qd, e_ex)
            line = (f'   n_par {n_par:3d}, N {g.n:5d}: waves: Phi ' + ' '.join(f'{e:.2e}' for e in e_phi)
                    + ' | energy ' + ' '.join(f'{e:.2e}' for e in e_en))
            if prev is not None:
                t0 = time.time()
                wc = wave_sums_carried(prev, g, window(12), qd)
                t_carry = time.time()-t0
                diff = max(abs(wc.A_off-ws.A_off).max()/abs(ws.A_off).max(),
                           abs(wc.B_off-ws.B_off).max()/abs(ws.B_off).max(),
                           abs(wc.S-ws.S).max()/abs(ws.S).max())
                line += f'\n      carried over: max rel. diff {diff:.1e}, time {t_carry:.1f} s vs {t_scratch:.1f} s'
            if g.n <= 13000:
                e_phi, e_en = errors(g, solve_points(g, qd), ph_ex, qd, e_ex)
                line += ('\n      points: Phi ' + ' '.join(f'{e:.2e}' for e in e_phi)
                         + ' | energy ' + ' '.join(f'{e:.2e}' for e in e_en))
            print(line)
            prev = ws

    if mode in ('all', 'lines'):
        print('\n3. field lines: 201 points per line, waves k = 12')
        for n_lines in (32, 64, 128):
            g = Grid(n_lines, 201)
            ph_ex, qd = g.data()
            ws = wave_sums_from_scratch(g, window(12), qd)
            e_phi, e_en = errors(g, ws.solve(), ph_ex, qd, e_ex)
            line = (f'   {n_lines:3d} lines, N {g.n:5d}: waves: Phi ' + ' '.join(f'{e:.2e}' for e in e_phi)
                    + ' | energy ' + ' '.join(f'{e:.2e}' for e in e_en))
            if g.n <= 13000:
                e_phi, e_en = errors(g, solve_points(g, qd), ph_ex, qd, e_ex)
                line += ('\n      points: Phi ' + ' '.join(f'{e:.2e}' for e in e_phi)
                         + ' | energy ' + ' '.join(f'{e:.2e}' for e in e_en))
            print(line)

    if mode == 'memory':
        n_par = int(sys.argv[2])
        g = Grid(n_alpha, n_par)
        ph_ex, qd = g.data()
        t0 = time.time()
        ws = wave_sums_from_scratch(g, window(12), qd)
        e_phi, e_en = errors(g, ws.solve(), ph_ex, qd, e_ex)
        print(f'memory run: N = {g.n}, {time.time()-t0:.0f} s, peak {peak_mb():.0f} MB '
              f'(one N x N real matrix: {8*g.n**2/2**20:.0f} MB); Phi '
              + ' '.join(f'{e:.2e}' for e in e_phi) + ' | energy ' + ' '.join(f'{e:.2e}' for e in e_en))


if __name__ == '__main__':
    main()
