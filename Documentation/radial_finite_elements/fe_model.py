"""Model problems for the radial discretisation of PB3D: finite differences against finite elements.

The continuous problem has M Fourier modes X(psi) on 0 <= psi <= 1 and the
energy of chapter 7 of the theory document,

    L[X] = int ( X^H L0 X + X^H L1 X' + X'^H L1^H X + X'^H L2 X' ) dpsi
           + X_e^H dvac X_e,          L^j = P^j - lambda K^j,

with X = 0 at psi = 0 and either a free edge (the term with dvac) or X = 0 at
psi = 1. Stationarity gives A X = lambda B X.

Every discretisation here is written as one table: a list of energy points,
each with a position, a weight, the unknowns it involves, and the numbers that
give X and X' at that position from those unknowns. The matrices follow from
the table in one loop (function assemble). The finite differences of PB3D and
all finite elements differ only in their table.

    python -B fe_model.py                 every section (about half an hour)
    python -B fe_model.py six hard        the named sections

Sections: one six hard ends sets spline gauss ten cost same post lengths.
It needs numpy and scipy.
"""
import sys

import numpy as np
import scipy.sparse as sp
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import BSpline, CubicSpline
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.sparse.linalg import eigsh

# ----------------------------------------------------------------------------
# continuous problems
# ----------------------------------------------------------------------------


class Problem:
    """Coefficient matrices as functions of psi, the vacuum matrix and the kind of edge."""

    def __init__(self, M, P, K, dvac, edge='free', name=''):
        self.M, self.P, self.K, self.dvac, self.edge, self.name = M, P, K, dvac, edge, name

    def coeff(self, psi):
        """P0, P1, P2, K0, K1, K2 at one position, each an (M, M) array."""
        return self.P(psi) + self.K(psi)


def one_mode(edge='free', V=lambda psi: 4000*(1 - psi)**6, kappa=1.0, a=0.0, b=1.0):
    """The one-mode test problem of section 7.10.6: -X'' + V X = lambda X on a <= psi <= b.

    The interval a ... b is mapped to 0 ... 1, so that the ends can be moved.
    """
    L = b - a
    one, zero = np.ones((1, 1)), np.zeros((1, 1))
    return Problem(1, lambda t: (V(a + L*t)*one, zero, one/L**2), lambda t: (one, zero, zero),
                   kappa/L*one, edge, 'one mode')


def coupled(edge='free', M=6, n=10, m0=12, q_in=1.21, q_e=1.74, drive=38.0, drive_at=0.8, k2_scale=1.0):
    """M coupled Fourier modes with rational surfaces inside the plasma.

    The coefficients are chosen, not computed from an equilibrium. They copy
    four properties of PB3D's matrices: P2 and P0 contain (n q - m)^2, which
    vanishes at the rational surface of the mode; K2 is not zero, so that
    L2 = P2 - lambda K2 stays positive for lambda < 0; P1 and K1 are complex
    and not Hermitian; neighbouring modes are coupled, with a drive that is
    largest at drive_at. A smaller k2_scale makes the eigenfunction change
    over a shorter length near the rational surfaces.
    """
    m = m0 + np.arange(M)
    off = np.eye(M, k=1) + np.eye(M, k=-1)
    skew = 1j*(np.eye(M, k=1) - np.eye(M, k=-1))

    def q(psi):
        return q_in + (q_e - q_in)*psi**2

    def P(psi):
        f = n*q(psi) - m
        D = drive*np.exp(-((psi - drive_at)/0.12)**2)
        P0 = np.diag(4.0*f**2 - D) - 0.45*D*off + 0.3*np.sin(3*psi)*off
        P1 = np.diag(0.6*f*(1 + 0.5j)) + 0.25*(1 + psi)*skew + 0.2*np.eye(M, k=1)*f[:, None]
        P2 = np.diag(f**2*(0.05 + 0.02*psi))
        return P0, P1, P2

    def K(psi):
        K0 = np.diag(1.0 + 0.3*psi + 0*m) + 0.1*off
        K1 = np.sqrt(k2_scale)*(0.02*skew + 0.01*np.eye(M, k=1))
        K2 = k2_scale*np.diag(0.004*(1 + psi) + 0*m)
        return K0, K1, K2

    dvac = np.diag(0.02*np.abs(m).astype(float)) + 0.004*off
    prob = Problem(M, P, K, dvac, edge, '%d modes' % M)
    prob.q, prob.n, prob.m = q, n, m
    return prob


def rational_surfaces(prob):
    """Number of modes of the problem whose rational surface lies inside the plasma."""
    f0, f1 = prob.n*prob.q(0.0) - prob.m, prob.n*prob.q(1.0) - prob.m
    return int(np.sum(f0*f1 < 0))


def splined(prob, n_eq):
    """The same problem with coefficients known only on n_eq points and interpolated.

    PB3D has its coefficient matrices on the radial points of the equilibrium
    grid and interpolates them with a cubic spline to where the energy needs them.
    """
    x = np.linspace(0, 1, n_eq)
    spl = CubicSpline(x, np.array([np.array(prob.coeff(xi)) for xi in x]), axis=0)
    return Problem(prob.M, lambda psi: tuple(spl(psi)[:3]), lambda psi: tuple(spl(psi)[3:]),
                   prob.dvac, prob.edge, prob.name + ', %d equilibrium points' % n_eq)


# ----------------------------------------------------------------------------
# tables
# ----------------------------------------------------------------------------


class Table:
    """One discretisation.

    psi[g], w[g]            position and weight of energy point g
    dof[g], phi[g], dphi[g] the unknowns it involves, and the numbers that
                            give X and X' there:  X = sum_a phi[a] X_dof[a]
    n_dof                   number of unknowns per Fourier mode
    held                    unknowns that are held at zero at the inner end
    held_edge               unknowns that are held at zero at a fixed edge
    edge                    the unknown that is X at the edge
    """

    def __init__(self, n_dof, held, held_edge, edge, label):
        self.n_dof, self.held, self.held_edge, self.edge, self.label = n_dof, held, held_edge, edge, label
        self.psi, self.w, self.dof, self.phi, self.dphi = [], [], [], [], []

    def add(self, psi, w, dof, phi, dphi):
        self.psi.append(psi), self.w.append(w), self.dof.append(np.asarray(dof))
        self.phi.append(np.asarray(phi, float)), self.dphi.append(np.asarray(dphi, float))


def fd_coefficients(p):
    """c_0 ... c_p of the left difference (D X)_i = sum_s c_s X_{i-s} / Delta."""
    s = np.arange(p + 1.0)
    rhs = np.zeros(p + 1)
    rhs[1] = -1.0
    return np.linalg.solve(np.vander(s, p + 1, increasing=True).T, rhs)


def fd_end_weights(p, n=60):
    """The weights of the sum rule, from the last point inwards."""
    c = fd_coefficients(p)
    w = [1/c[0]]
    while len(w) < n:
        k = len(w)
        w.append(-sum(c[s]*w[k-s] for s in range(1, min(p, k) + 1))/c[0])
    return np.array(w)


def table_fd(N, p=3, edge='free', weights=None):
    """PB3D as it is: N equally spaced points, left difference of order p.

    The energy points are the points p+1 ... N. The first p points are held at
    zero. At a free edge the last points have the weights of the sum rule
    (p = 2, 3) or the last point has the weight 1/2 (p = 1). At a fixed edge
    all weights are 1. `weights` ('sum rule', 'last 1/2', 'all 1') overrides this.
    """
    delta = 1.0/(N - 1)
    c = fd_coefficients(p)
    w = np.ones(N)
    weights = weights or ('all 1' if edge == 'fixed' else 'sum rule')
    if weights == 'sum rule' and p >= 2:
        ww = fd_end_weights(p)[:N-p-1]
        w[N-1-np.arange(len(ww))] = ww
    elif weights in ('sum rule', 'last 1/2'):
        w[-1] = 0.5
    tab = Table(N, list(range(p)), [N-1], N-1, 'differences p = %d' % p)
    for i in range(p, N):
        phi = np.zeros(p + 1)
        phi[0] = 1.0
        tab.add(i*delta, delta*w[i], i - np.arange(p + 1), phi, c/delta)
    return tab


def lobatto(n):
    """Gauss-Lobatto points and weights on -1 ... 1."""
    if n == 2:
        return np.array([-1.0, 1.0]), np.array([1.0, 1.0])
    Pn = np.polynomial.legendre.Legendre.basis(n - 1)
    x = np.r_[-1.0, np.sort(Pn.deriv().roots().real), 1.0]
    return x, 2/(n*(n - 1)*Pn(x)**2)


def rule(kind, n):
    """Points and weights on -1 ... 1."""
    if kind == 'gauss':
        return leggauss(n)
    if kind == 'lobatto':
        return lobatto(n)
    if kind == 'outer end':
        return np.array([1.0]), np.array([2.0])
    raise ValueError(kind)


def table_elements(n_el, k, quad='gauss', n_q=None, points='equal'):
    """Polynomials of degree k on every element, continuous between elements.

    The unknowns are the values of X at k+1 points per element; neighbouring
    elements share their end point. The points inside an element are equally
    spaced ('equal') or at the Gauss-Lobatto positions ('lobatto'). The energy
    points are n_q Gauss-Legendre points per element (default k + 1).
    """
    h = 1.0/n_el
    xn = lobatto(k + 1)[0] if points == 'lobatto' else np.linspace(-1, 1, k + 1)
    C = np.linalg.inv(np.vander(xn, k + 1, increasing=True))       # column a: coefficients of polynomial a
    xq, wq = rule(quad, n_q or k + 1)
    pw = np.vander(xq, k + 1, increasing=True)
    phi = pw@C
    dpw = np.zeros_like(pw)
    dpw[:, 1:] = pw[:, :-1]*np.arange(1, k + 1)
    dphi = dpw@C
    tab = Table(n_el*k + 1, [0], [n_el*k], n_el*k, 'elements k = %d, %s %d' % (k, quad, len(xq)))
    for e in range(n_el):
        for g in range(len(xq)):
            tab.add((e + (xq[g] + 1)/2)*h, wq[g]*h/2, e*k + np.arange(k + 1), phi[g], dphi[g]*2/h)
    return tab


def table_hermite(n_el, n_q=4):
    """Cubic polynomials on every element, with X and X' continuous.

    Two unknowns per element boundary: X, and X' times the length of an element.
    """
    h = 1.0/n_el
    xq, wq = leggauss(n_q)
    t = (xq + 1)/2
    phi = np.array([1 - 3*t**2 + 2*t**3, t - 2*t**2 + t**3, 3*t**2 - 2*t**3, -t**2 + t**3]).T
    dphi = np.array([-6*t + 6*t**2, 1 - 4*t + 3*t**2, 6*t - 6*t**2, -2*t + 3*t**2]).T
    tab = Table(2*(n_el + 1), [0], [2*n_el], 2*n_el, 'cubic Hermite')
    for e in range(n_el):
        for g in range(len(xq)):
            tab.add((e + t[g])*h, wq[g]*h/2, 2*e + np.arange(4), phi[g], dphi[g]/h)
    return tab


def table_bspline(n_el, k=3, n_q=None):
    """Splines of degree k with k-1 continuous derivatives; one unknown per element plus k.

    The first and the last unknown are the values of X at the two ends.
    """
    xb = np.linspace(0, 1, n_el + 1)
    knots = np.r_[[0.0]*k, xb, [1.0]*k]
    n = n_el + k
    xq, wq = leggauss(n_q or k + 1)
    tab = Table(n, [0], [n-1], n-1, 'B-spline k = %d' % k)
    spl = BSpline(knots, np.eye(n), k)
    dspl = spl.derivative()
    h = 1.0/n_el
    for e in range(n_el):
        x = (e + (xq + 1)/2)*h
        val, der = spl(x), dspl(x)
        for g in range(len(xq)):
            tab.add(x[g], wq[g]*h/2, e + np.arange(k + 1), val[g, e:e+k+1], der[g, e:e+k+1])
    return tab


# ----------------------------------------------------------------------------
# from the table to the matrices
# ----------------------------------------------------------------------------


def assemble(prob, tab, masks=None):
    """A and B from the table: one loop over the energy points.

    masks = (at_point, at_dof) says which Fourier modes exist where: boolean
    arrays of shape (number of energy points, M) and (n_dof, M). A coefficient
    is used only if both of its modes exist at the energy point and at the two
    unknowns that it couples. An unknown whose mode does not exist is left out.
    Without masks all modes exist everywhere.
    """
    M, n = prob.M, tab.n_dof
    rows, cols, va, vb = [], [], [], []
    for g in range(len(tab.psi)):
        P0, P1, P2, K0, K1, K2 = prob.coeff(tab.psi[g])
        f, d, dof = tab.phi[g], tab.dphi[g], tab.dof[g]
        for a in range(len(dof)):
            for b in range(len(dof)):
                ca = f[a]*f[b], f[a]*d[b], d[a]*f[b], d[a]*d[b]
                if not any(ca):
                    continue
                blkA = tab.w[g]*(ca[0]*P0 + ca[1]*P1 + ca[2]*P1.conj().T + ca[3]*P2)
                blkB = tab.w[g]*(ca[0]*K0 + ca[1]*K1 + ca[2]*K1.conj().T + ca[3]*K2)
                if masks is not None:
                    use = np.outer(masks[0][g] & masks[1][dof[a]], masks[0][g] & masks[1][dof[b]])
                    blkA, blkB = blkA*use, blkB*use
                r = dof[a]*M + np.arange(M)
                c = dof[b]*M + np.arange(M)
                rows.append(np.repeat(r, M)), cols.append(np.tile(c, M))
                va.append(blkA.ravel()), vb.append(blkB.ravel())
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    A = sp.csr_matrix((np.concatenate(va), (rows, cols)), shape=(n*M, n*M), dtype=complex)
    B = sp.csr_matrix((np.concatenate(vb), (rows, cols)), shape=(n*M, n*M), dtype=complex)
    A = A.tolil()
    if prob.edge == 'free':
        e = tab.edge*M
        A[e:e+M, e:e+M] = A[e:e+M, e:e+M].toarray() + prob.dvac
    held = list(tab.held) + (list(tab.held_edge) if prob.edge == 'fixed' else [])
    keep = np.ones(n*M, bool) if masks is None else masks[1].ravel().copy()
    for a in held:
        keep[a*M:(a+1)*M] = False
    A, B = A.tocsr()[keep][:, keep], B[keep][:, keep]
    return A, B, keep


def lowest(prob, tab, guess=None, masks=None, number=1):
    """Lowest eigenvalue of the matrices of a table (the lowest `number` as an array if number > 1)."""
    A, B, keep = assemble(prob, tab, masks)
    if A.shape[0] <= 1200 or guess is None:
        val = eigh(A.toarray(), B.toarray(), eigvals_only=True)
    else:
        val = np.sort(eigsh(A.tocsc(), k=max(2, number), M=B.tocsc(), sigma=guess, which='LM')[0])
    return val[0] if number == 1 else val[:number]


def reference(prob, n_el=None):
    """True lowest eigenvalue: polynomials of degree 8 on many elements, checked by doubling."""
    n_el = n_el or (40 if prob.M == 1 else 160)
    high = lambda n: table_elements(n, 8, n_q=11, points='lobatto')
    coarse = lowest(prob, high(max(8, 1000//(8*prob.M))))
    guess = coarse - 0.05*abs(coarse)
    a = lowest(prob, high(n_el), guess=guess)
    b = lowest(prob, high(2*n_el), guess=guess)
    return b, abs(a - b)


def length_of_change(prob):
    """sqrt(int |X|^2 / int |X'|^2) of the lowest eigenfunction: the length over which it changes."""
    tab = table_elements(max(8, 2400//(8*prob.M)), 8, n_q=11, points='lobatto')
    A, B, keep = assemble(prob, tab)
    vec = eigh(A.toarray(), B.toarray())[1][:, 0]
    full = np.zeros(tab.n_dof*prob.M, complex)
    full[keep] = vec
    X = full.reshape(tab.n_dof, prob.M)
    num = sum(tab.w[g]*np.sum(abs(tab.phi[g]@X[tab.dof[g]])**2) for g in range(len(tab.psi)))
    den = sum(tab.w[g]*np.sum(abs(tab.dphi[g]@X[tab.dof[g]])**2) for g in range(len(tab.psi)))
    return np.sqrt(num/den)


# ----------------------------------------------------------------------------
# the sliding set of Fourier modes
# ----------------------------------------------------------------------------


def set_at(prob, width):
    """The sliding set of Fourier modes: the `width` modes closest to n q(psi)."""
    def w(psi):
        lo = int(np.rint(prob.n*prob.q(psi) - (width - 1)/2))
        return (prob.m >= lo) & (prob.m < lo + width)
    return w


def sets_per_point(prob, tab, width):
    """PB3D as it is: every point has its own set of modes."""
    w = set_at(prob, width)
    return np.array([w(x) for x in tab.psi]), np.array([w(x) for x in np.linspace(0, 1, tab.n_dof)])


def sets_rule_A(prob, tab, width, n_el, k, where='middle'):
    """One set per element; a point shared by two elements has the unknowns of the element on its inner side.

    The set of an element is the one at its middle (or at its last point,
    where='last'). A mode is zero at the first point of its first element, and
    has a free value at the last point of its last element, with no energy
    beyond it. For many points this is what the differences of PB3D do with
    one set per point.
    """
    w = set_at(prob, width)
    h = 1.0/n_el
    at_el = np.array([w((e + (0.5 if where == 'middle' else 1.0))*h) for e in range(n_el)])
    at_point = np.repeat(at_el, len(tab.psi)//n_el, axis=0)
    at_dof = np.ones((tab.n_dof, prob.M), bool)
    for e in range(n_el):
        at_dof[e*k+1:e*k+k+1] = at_el[e]
    at_dof[0] = at_el[0]
    return at_point, at_dof


def sets_rule_B(prob, tab, width, n_el, k):
    """One set per element; a shared point has only the modes that both elements have.

    Every mode is then a continuous function that is zero outside the elements
    whose set contains it, and the matrices hold the complete energy of that
    function.
    """
    w = set_at(prob, width)
    h = 1.0/n_el
    at_el = np.array([w((e + 0.5)*h) for e in range(n_el)])
    at_point = np.repeat(at_el, len(tab.psi)//n_el, axis=0)
    at_dof = np.ones((tab.n_dof, prob.M), bool)
    for e in range(n_el):
        at_dof[e*k:e*k+k+1] &= at_el[e]
    return at_point, at_dof


# ----------------------------------------------------------------------------
# sections
# ----------------------------------------------------------------------------


def power(errs, ns):
    e, n = np.abs(np.array(errs[-2:], float)), np.array(ns[-2:], float)
    return np.log(e[0]/e[1])/np.log(n[1]/n[0]) if e.min() > 0 else np.nan


def scan(prob, ref, makers, sizes, guess=None):
    """Relative error of the lowest eigenvalue against the number of unknowns per mode."""
    print('  %-36s' % 'scheme' + ''.join('%11s' % ('n = %d' % n) for n in sizes) + '   falls as')
    for label, make in makers:
        errs, ns = [], []
        for n in sizes:
            tab = make(n)
            errs.append((lowest(prob, tab, guess=guess) - ref)/abs(ref)), ns.append(tab.n_dof)
        print('  %-36s' % label + ''.join('%+11.1e' % e for e in errs) + '   1/n^%.1f' % power(errs, ns))


def makers(edge, short=False):
    """Schemes at (about) the same number n of unknowns per mode."""
    out = [
        ('differences p = 2 (PB3D)', lambda n: table_fd(n, 2, edge)),
        ('differences p = 3 (PB3D)', lambda n: table_fd(n, 3, edge)),
        ('elements k = 1', lambda n: table_elements(n - 1, 1)),
        ('elements k = 2', lambda n: table_elements((n - 1)//2, 2)),
        ('elements k = 3', lambda n: table_elements((n - 1)//3, 3)),
    ]
    if not short:
        out += [
            ('elements k = 2, Lobatto 3 (panels)', lambda n: table_elements((n - 1)//2, 2, quad='lobatto')),
            ('elements k = 3, Lobatto 4 (panels)', lambda n: table_elements((n - 1)//3, 3, quad='lobatto')),
            ('cubic Hermite', lambda n: table_hermite(n//2 - 1)),
            ('B-spline k = 3', lambda n: table_bspline(n - 3, 3)),
        ]
    return out


def section_one():
    for edge in ('free', 'fixed'):
        prob = one_mode(edge)
        ref, acc = reference(prob)
        print('\nOne mode, %s edge: lowest eigenvalue %.10f (known to %.0e)' % (edge, ref, acc))
        scan(prob, ref, makers(edge), (49, 97, 193))
    prob = one_mode('free', V=lambda psi: 0*psi)
    ref, acc = reference(prob)
    print('\nOne mode, free edge, V = 0 (eigenfunction not small at the inner end): %.10f' % ref)
    scan(prob, ref, makers('free', short=True), (49, 97, 193))


def section_six():
    for edge in ('free', 'fixed'):
        prob = coupled(edge)
        ref, acc = reference(prob)
        print('\nSix modes, %d rational surfaces inside, %s edge: lowest eigenvalue %.10f (known to %.0e)'
              % (rational_surfaces(prob), edge, ref, acc))
        print('  the eigenfunction changes over the length %.3f' % length_of_change(prob))
        scan(prob, ref, makers(edge), (97, 193, 385, 769), guess=ref - 0.3*abs(ref))


def section_hard():
    prob = coupled('free', M=20, n=30, m0=37, k2_scale=0.1)
    ref, acc = reference(prob, n_el=120)
    L = length_of_change(prob)
    print('\nTwenty modes, %d rational surfaces inside, K2 ten times smaller, free edge: %.10f (known to %.0e)'
          % (rational_surfaces(prob), ref, acc))
    print('  the eigenfunction changes over the length %.4f: %s points per length for n = 49, 193, 769'
          % (L, ', '.join('%.1f' % (L*(n - 1)) for n in (49, 193, 769))))
    scan(prob, ref, makers('free', short=True), (49, 97, 193, 385, 769), guess=ref - 0.3*abs(ref))


def section_ends():
    """By how many steps the differences of PB3D displace the two ends."""
    h = 1e-3
    V0 = lambda psi: 0*psi
    Ve = lambda psi: 4000*(1 - psi)**6
    print('\nThe ends of the differences (one mode). Error of the eigenvalue = Delta x (steps) x d lambda / d (position)')
    ref = reference(one_mode('free', V=V0))[0]
    slope = (reference(one_mode('free', V=V0, a=h))[0] - reference(one_mode('free', V=V0, a=-h))[0])/(2*h)
    for p in (2, 3):
        x = [(lowest(one_mode('free', V=V0), table_fd(N, p), guess=0.7*ref) - ref)*(N - 1)/slope for N in (769, 1537)]
        print('  inner end, p = %d: the plasma starts %.4f steps further out (p - 1/2 = %.1f)' % (p, 2*x[1] - x[0], p - 0.5))
    ref = reference(one_mode('fixed', V=Ve))[0]
    slope = (reference(one_mode('fixed', V=Ve, b=1 - h))[0] - reference(one_mode('fixed', V=Ve, b=1 + h))[0])/(2*h)
    for p in (1, 2, 3):
        x = [(lowest(one_mode('fixed', V=Ve), table_fd(N, p, 'fixed'), guess=0.7*ref) - ref)*(N - 1)/slope
             for N in (769, 1537)]
        print('  fixed edge with all weights 1, p = %d: the plasma ends %.4f steps further in' % (p, 2*x[1] - x[0]))


def section_sets():
    M_all, m0 = 17, 7
    variants = [('q rising, drive at 0.8', dict(), (97, 193, 385, 769)),
                ('q falling, drive at 0.8', dict(q_in=1.74, q_e=1.21), (385, 769)),
                ('q rising, drive at 0.4', dict(drive_at=0.4), (385, 769))]
    print('\nSliding set of Fourier modes out of %d, free edge' % M_all)
    print('  relative distance of the lowest eigenvalue to the one with all %d modes everywhere' % M_all)
    for name, kw, sizes in variants:
        prob = coupled('free', M=M_all, m0=m0, **kw)
        full, acc = reference(prob, n_el=120)
        g = full - 0.3*abs(full)
        print('  %s: all modes %.10f (known to %.0e)' % (name, full, acc))
        for width in (8, 10, 12):
            print('    set of %d' % width)
            print('      %-46s' % 'scheme' + ''.join('%11s' % ('n = %d' % n) for n in sizes))
            k = 3
            rows = [('differences p = 3, set per point (PB3D now)',
                     lambda n: (table_fd(n, 3), lambda t: sets_per_point(prob, t, width))),
                    ('elements k = 3, rule A',
                     lambda n: (table_elements((n - 1)//k, k), lambda t: sets_rule_A(prob, t, width, (n - 1)//k, k))),
                    ('elements k = 3, rule A, set at the last point',
                     lambda n: (table_elements((n - 1)//k, k),
                                lambda t: sets_rule_A(prob, t, width, (n - 1)//k, k, 'last'))),
                    ('elements k = 3, rule B',
                     lambda n: (table_elements((n - 1)//k, k), lambda t: sets_rule_B(prob, t, width, (n - 1)//k, k)))]
            for label, make in rows:
                out = []
                for n in sizes:
                    tab, mk = make(n)
                    out.append('%+11.2e' % ((lowest(prob, tab, guess=g, masks=mk(tab)) - full)/abs(full)))
                print('      %-46s' % label + ''.join(out))


def section_spline():
    base = coupled('free')
    exact = reference(base)[0]
    print('\nCoefficients known on n_eq points and interpolated by a cubic spline (six modes, free edge)')
    for n_eq in (41, 161):
        prob = splined(base, n_eq)
        ref = lowest(prob, table_elements(8*(n_eq - 1), 6, n_q=8, points='lobatto'), guess=exact - 0.3*abs(exact))
        print('  n_eq = %d: lowest eigenvalue %.10f; with the true coefficients %.10f (distance %.1e)'
              % (n_eq, ref, exact, (ref - exact)/abs(exact)))
        print('  error against the problem with the interpolated coefficients:')
        scan(prob, ref, [
            ('differences p = 3 (PB3D)', lambda n: table_fd(n, 3)),
            ('elements k = 2', lambda n: table_elements((n - 1)//2, 2)),
            ('elements k = 3', lambda n: table_elements((n - 1)//3, 3)),
            ('elements k = 3, 6 Gauss points', lambda n: table_elements((n - 1)//3, 3, n_q=6)),
            ('B-spline k = 3', lambda n: table_bspline(n - 3, 3)),
        ], (97, 193, 385), guess=ref - 0.3*abs(ref))
        if n_eq == 41:
            print('  element boundaries on the points of the spline:')
            scan(prob, ref, [('elements k = 3', lambda n: table_elements((n - 1)//3, 3))],
                 (121, 241, 481), guess=ref - 0.3*abs(ref))


def section_gauss():
    prob = coupled('free')
    ref = reference(prob)[0]
    print('\nNumber of Gauss points per element (six modes, free edge)')
    rows = []
    for k in (1, 2, 3):
        for n_q in (k, k + 1, k + 2):
            rows.append(('elements k = %d, Gauss %d' % (k, n_q),
                         lambda n, k=k, n_q=n_q: table_elements((n - 1)//k, k, n_q=n_q)))
    scan(prob, ref, rows, (97, 193, 385), guess=ref - 0.3*abs(ref))


def section_ten():
    """The lowest ten eigenvalues, not only the first."""
    print('\nThe lowest ten eigenvalues (six modes), 193 unknowns per mode, against polynomials of degree 8')
    print('  %-30s %-8s %22s %26s' % ('scheme', 'edge', 'largest relative error', 'number below the reference'))
    for edge in ('free', 'fixed'):
        prob = coupled(edge)
        ref = lowest(prob, table_elements(60, 8, n_q=11, points='lobatto'), number=10)
        for label, tab in [('differences p = 3 (PB3D)', table_fd(193, 3, edge)),
                           ('elements k = 1', table_elements(192, 1)),
                           ('elements k = 2', table_elements(96, 2)),
                           ('elements k = 3', table_elements(64, 3)),
                           ('elements k = 3, Gauss 3', table_elements(64, 3, n_q=3)),
                           ('elements k = 3, Lobatto 4', table_elements(64, 3, quad='lobatto'))]:
            val = lowest(prob, tab, number=10)
            rel = (val - ref)/abs(ref)
            print('  %-30s %-8s %22.1e %26d' % (label, edge, abs(rel).max(), int(np.sum(rel < -1e-9))))
        print('  (reference: %.6f ... %.6f, all negative)' % (ref[0], ref[-1]))


def section_cost():
    print('\nSize of the matrices for 193 unknowns per mode (blocks of M x M)')
    print('  %-30s %14s %22s' % ('scheme', 'energy points', 'blocks per row'))
    n = 193
    for label, tab in [('differences p = 3', table_fd(n, 3)),
                       ('elements k = 1', table_elements(n - 1, 1)),
                       ('elements k = 2', table_elements((n - 1)//2, 2)),
                       ('elements k = 3', table_elements((n - 1)//3, 3)),
                       ('cubic Hermite', table_hermite(n//2 - 1)),
                       ('B-spline k = 3', table_bspline(n - 3, 3))]:
        pat = sp.lil_matrix((tab.n_dof, tab.n_dof))
        for dof in tab.dof:
            for a in dof:
                pat[a, dof] = 1
        nnz = pat.tocsr().getnnz(axis=1)
        print('  %-30s %14d %22s' % (label, len(tab.psi), '%.1f mean, %d at most' % (nnz.mean(), nnz.max())))


def section_same():
    """Two tables that must give the same matrices."""
    prob = coupled('fixed')
    print('\nSame matrices from two tables (six modes, fixed edge, 97 points)')
    a = lowest(prob, table_fd(97, 1, 'fixed'))
    b = lowest(prob, table_elements(96, 1, quad='outer end'))
    print('  differences p = 1 with all weights 1                    %.14f' % a)
    print('  elements k = 1, one energy point at the outer end       %.14f   (difference %.1e)' % (b, abs(a - b)))
    a = lowest(prob, table_elements(32, 3, points='equal'))
    b = lowest(prob, table_elements(32, 3, points='lobatto'))
    print('  elements k = 3, equally spaced points                   %.14f' % a)
    print('  elements k = 3, points at the Gauss-Lobatto positions   %.14f   (difference %.1e)' % (b, abs(a - b)))


def section_post():
    """X' at the radial points, from the point values of the computed eigenfunction, in two ways."""
    s = brentq(lambda s: s*np.cos(s) + np.sin(s), 1.6, 3.0, xtol=1e-15)
    prob = one_mode('free', V=lambda psi: 0*psi)
    print('\nX\' at the radial points from the values at the points (one mode, V = 0, free edge; exact X = sin(%.6f psi))' % s)
    print('  largest error at the points that are at least k points from an end, X scaled to the exact one')
    print('  %-16s %10s %24s %26s %32s' % ('scheme', 'unknowns', 'error of the eigenvalue', 'X\' of the element polynomials',
                                         'X\' of a cubic spline through X'))
    for k in (1, 2, 3):
        xn = np.linspace(-1, 1, k + 1)
        C = np.linalg.inv(np.vander(xn, k + 1, increasing=True))
        dpw = np.zeros((k + 1, k + 1))
        dpw[:, 1:] = np.vander(xn, k + 1, increasing=True)[:, :-1]*np.arange(1, k + 1)
        for n_el in (16, 32, 64):
            tab = table_elements(n_el, k)
            A, B, keep = assemble(prob, tab)
            val, vec = eigh(A.toarray(), B.toarray())
            X = np.zeros(tab.n_dof, complex)
            X[keep] = vec[:, 0]
            x = np.linspace(0, 1, tab.n_dof)
            ex, dex = np.sin(s*x), s*np.cos(s*x)
            X = (X*(X.conj()@ex)/(X.conj()@X)).real
            dphi = dpw@C*2*n_el                       # d/dpsi of polynomial b at point a of an element
            d_el, count = np.zeros(tab.n_dof), np.zeros(tab.n_dof)
            for e in range(n_el):
                idx = e*k + np.arange(k + 1)
                d_el[idx] += dphi@X[idx]              # a point shared by two elements gets the mean of the two sides
                count[idx] += 1
            d_el /= count
            d_sp = CubicSpline(x, X)(x, 1)
            inner = slice(k, -k)
            print('  %-16s %10d %24.1e %26.1e %32.1e' % ('elements k = %d' % k, tab.n_dof, (val[0] - s*s)/s**2,
                                                       abs(d_el - dex)[inner].max(), abs(d_sp - dex)[inner].max()))


def lengths(X):
    """Three lengths of X[point, mode], given on equally spaced points over a range of 1:
    sqrt(int |X|^2 / int |X'|^2), sqrt(int |X'|^2 / int |X''|^2), sqrt(int |X''|^2 / int |X'''|^2)."""
    d = 1.0/(X.shape[0] - 1)
    s = [np.sum(abs(np.diff(X, j, axis=0)/d**j)**2) for j in range(4)]
    return [np.sqrt(s[j]/s[j+1]) for j in range(3)]


def section_lengths():
    """How fast the lowest eigenfunction and its derivatives change (compare with eigenfunction_profile.py)."""
    print('\nLengths of the lowest eigenfunction, free edge, from its values at 769 points (elements k = 3)')
    print('  %-16s %22s %24s %26s' % ('problem', 'sqrt(|X|^2/|X\'|^2)', 'sqrt(|X\'|^2/|X\'\'|^2)', 'sqrt(|X\'\'|^2/|X\'\'\'|^2)'))
    for name, prob in [('six modes', coupled('free')), ('twenty modes', coupled('free', M=20, n=30, m0=37, k2_scale=0.1))]:
        tab = table_elements(256, 3)
        A, B, keep = assemble(prob, tab)
        guess = lowest(prob, table_elements(40, 3))
        vec = eigsh(A.tocsc(), k=1, M=B.tocsc(), sigma=guess - 0.05*abs(guess), which='LM')[1][:, 0]
        full = np.zeros(tab.n_dof*prob.M, complex)
        full[keep] = vec
        print('  %-16s %22.4f %24.4f %26.4f' % (name, *lengths(full.reshape(tab.n_dof, prob.M))))


SECTIONS = dict(one=section_one, six=section_six, hard=section_hard, ends=section_ends, sets=section_sets,
                spline=section_spline, gauss=section_gauss, ten=section_ten, cost=section_cost, same=section_same,
                post=section_post, lengths=section_lengths)

if __name__ == '__main__':
    np.set_printoptions(linewidth=150)
    for name in sys.argv[1:] or SECTIONS:
        SECTIONS[name]()
