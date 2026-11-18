"""Area D: the weights of the last points of a free edge (num_utilities.calc_end_weights,
used by SLEPC_ops.fill_mat since cba0762).

PB3D sums the energy over the radial points i = 0 ... N with weights h w_i, and
takes the radial derivative at point i as a left finite difference over p + 1
points (p = norm_disc_prec_sol):

    (D X)_i = (1/h) sum_{m=0}^{p} c_m X_{i-m}.

Part 1 computes the end weights with the recurrence of calc_end_weights,

    c_0 w_N = 1,   sum_{m=0}^{min(p,k)} c_m w_{N-k+m} = 0  for k >= 1,

and checks the property they are made for: h sum_i w_i (D X)_i = X_N for any
grid function X that vanishes at the first p points (the discrete form of the
integral of X' over the plasma).

Part 2 is a one-mode copy of the assembly of fill_mat, set_BC_1 and set_BC_2
for the model problem

    -X'' + V X = lambda X  on 0 < r < 1,   X(0) = 0,   X'(1) + kappa X(1) = 0,

which is PB3D's Euler equation and natural boundary condition with P2 = 1,
P0 = V, K0 = 1, P1 = K1 = K2 = 0 and vac%res = kappa. V = 4000 (1 - r)^6 puts
the mode against the edge, so that the edge alone sets the error. It compares
three choices for the weights of the last points, all others being 1:

    'h'    w_N = 1                       (2.47)
    'h/2'  w_N = 1/2                     (4915d4b, F5)
    'end'  the weights of part 1         (cba0762, for p = 2 and 3)

The reference eigenvalue comes from Chebyshev collocation on 200 points.

Part 3 repeats the end weights with V = 0, so that the mode also reaches the
inner end, where set_BC_1 sets the first p points to zero.
"""
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eig
from scipy.sparse.linalg import eigsh

from _common import BLUE, GREEN, ORANGE, plt, save

KAPPA = 1.0
EV_BC = 1.0e6                                  # artificial eigenvalue of the inner boundary rows


def V(r):
    return 4000 * (1 - r) ** 6


def fd(offsets):
    """Coefficients of the first derivative at offset 0 from the given offsets (step 1)."""
    o = np.asarray(offsets, float)
    rhs = np.zeros(len(o))
    rhs[1] = 1.0
    return np.linalg.solve(np.vander(o, len(o), increasing=True).T, rhs)


def end_weights(p, tol=1e-12):
    """The recurrence of calc_end_weights: w[0] is the last point, w[1] the one before it, ..."""
    c = fd(np.arange(0, -p - 1, -1))           # c[m] multiplies X_{i-m}
    w = [1 / c[0]]
    while True:
        k = len(w)
        w.append(-sum(c[m] * w[k - m] for m in range(1, min(p, k) + 1)) / c[0])
        if k + 1 > p and max(abs(np.array(w[-p:]) - 1)) < tol:
            return np.array(w[:-p])


def weights(N, p, kind):
    w = np.ones(N + 1)
    tail = {"h": [1.0], "h/2": [0.5], "end": end_weights(p)}[kind]
    for m, wm in enumerate(tail):
        w[N - m] = wm
    return w


def assemble(N, p, kind, V=V):
    h = 1.0 / N
    r = np.arange(N + 1) * h
    w = weights(N, p, kind) * h
    c = fd(np.arange(0, -p - 1, -1)) / h
    A = sp.lil_matrix((N + 1, N + 1))
    B = sp.lil_matrix((N + 1, N + 1))
    for i in range(p, N + 1):                  # fill_mat: rows from ndps+1 (1-based) to the edge
        st = i - np.arange(p + 1)
        A[i, i] += w[i] * V(r[i])
        B[i, i] += w[i]
        for a, ca in zip(st, c):
            for b, cb in zip(st, c):
                A[a, b] += w[i] * ca * cb
    for k in range(p):                         # set_BC_1 at the inner end
        for d in range(1, p + 1):
            A[k, k + d] = A[k + d, k] = B[k, k + d] = B[k + d, k] = 0.0
        A[k, k] = EV_BC
        B[k, k] = 1.0
    A[N, N] += KAPPA                           # set_BC_2: vac%res on the edge diagonal
    return A.tocsc(), B.tocsc()


def reference(n=200, V=V):
    x = np.cos(np.pi * np.arange(n + 1) / n)
    cc = np.r_[2, np.ones(n - 1), 2] * (-1) ** np.arange(n + 1)
    xx = np.tile(x, (n + 1, 1)).T
    D = np.outer(cc, 1 / cc) / (xx - xx.T + np.eye(n + 1))
    D -= np.diag(D.sum(axis=1))
    r = (x + 1) / 2
    D = 2 * D                                  # r = 1 is index 0, r = 0 is index n
    L = -D @ D + np.diag(V(r))
    M = np.eye(n + 1)
    L[n] = 0; L[n, n] = 1; M[n] = 0            # X(0) = 0
    L[0] = D[0]; L[0, 0] += KAPPA; M[0] = 0    # X'(1) + kappa X(1) = 0
    ev = eig(L, M, right=False)
    ev = ev[np.isfinite(ev) & (ev.real > 0)].real
    return ev.min()


def lowest(N, p, kind, near, V=V):
    A, B = assemble(N, p, kind, V)
    return eigsh(A, k=1, M=B, sigma=near, which="LM", return_eigenvectors=False)[0]


# part 1: the weights and the property they are made for
print("part 1: end weights of calc_end_weights, last point first")
rng = np.random.default_rng(1)
for p in (1, 2, 3):
    w = end_weights(p)
    if np.allclose(w, 1.0):
        print(f"  p = {p}: all weights are 1 (PB3D does not use them: it keeps 1/2 at the edge point)")
    else:
        print(f"  p = {p}: {len(w)} weights differ from 1 by more than 1e-12: "
              + ", ".join(f"{x:.4f}" for x in w[:6]) + " ...")
    print(f"         sum of (w - 1) = {np.sum(w - 1):+.6f}   (the trapezoidal rule: -0.5)")
    N = 50
    h = 1.0 / N
    X = rng.standard_normal(N + 1)
    X[:p] = 0.0
    c = fd(np.arange(0, -p - 1, -1)) / h
    DX = np.array([sum(c[m] * X[i - m] for m in range(p + 1)) for i in range(p, N + 1)])
    for kind in ("h/2", "end"):
        lhs = h * np.sum(weights(N, p, kind)[p:] * DX)
        print(f"         {kind:>3} weights: h sum w (DX) - X_N = {lhs - X[N]:+.2e}")
print("  (p = 3: 6/11 = %.4f, 108/121 = %.4f)" % (6 / 11, 108 / 121))

# part 2: the lowest eigenvalue of the model
lam = reference()
Ns = np.array([100, 200, 400, 800, 1600, 3200])
print(f"\npart 2: model eigenvalue, reference {lam:.9f} (Chebyshev, 200 points)")
print("  error of the lowest eigenvalue; error*N constant means first order")
res = {}
for p, kind in ((3, "h"), (3, "h/2"), (3, "end"), (2, "h/2"), (2, "end")):
    err = np.array([lowest(N, p, kind, lam) - lam for N in Ns])
    res[(p, kind)] = err
    if len(set(np.sign(err[-3:]))) == 1:
        order = -np.polyfit(np.log(Ns[-3:]), np.log(np.abs(err[-3:])), 1)[0]
        tail = f"   order from the last three: {order:.1f}"
    else:
        tail = "   (changes sign)"
    print(f"  p = {p}, weights {kind:>3}: error " + " ".join(f"{e:+9.2e}" for e in err) + tail)
    print(f"  {'':22} error*N " + " ".join(f"{e * N:+9.3f}" for e, N in zip(err, Ns)))

# part 3: the inner end
lam0 = reference(V=lambda r: 0 * r)
err0 = np.array([lowest(N, 3, "end", lam0, V=lambda r: 0 * r) - lam0 for N in Ns])
print(f"\npart 3: V = 0, reference {lam0:.9f}: the mode reaches the inner end")
print(f"  p = 3, weights end: error " + " ".join(f"{e:+9.2e}" for e in err0))
print(f"  {'':22} error*N " + " ".join(f"{e * N:+9.3f}" for e, N in zip(err0, Ns)))

fig, ax = plt.subplots(figsize=(5.6, 3.8))
ax.loglog(Ns, np.abs(res[(3, "h")]), "o-", color=ORANGE, label="p = 3, w_N = 1 (2.47)")
ax.loglog(Ns, np.abs(res[(3, "h/2")]), "s-", color=BLUE, label="p = 3, w_N = 1/2 (F5)")
ax.loglog(Ns, np.abs(res[(3, "end")]), "^-", color=GREEN, label="p = 3, end weights (cba0762)")
ax.set_xlabel("number of radial intervals N")
ax.set_ylabel("|error| of the lowest eigenvalue")
ax.set_title("one-mode copy of fill_mat at a free edge")
ax.legend(fontsize=8)
save(fig, "end_weights_1d.png")
