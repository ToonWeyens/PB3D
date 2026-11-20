"""Area 2e: the quadrature weight of the free edge point (SLEPC_ops.fill_mat, w_q).

PB3D assembles the energy integrals over the radial points as sum_k h f_k.
At 2.47 every point had the weight h, also the last one. With a fixed boundary
the last row is replaced by the boundary condition, so its weight never
matters. With a free boundary (BC_style(2) = 2 or 4) the last point stays in
the system and the trapezoidal rule gives it h/2. Commit 4915d4b sets
w_q = 1/2 there.

This is a MODEL of the principle, not PB3D's assembly: the 1-D eigenproblem

    -u'' = lambda u on [0, 1],  u(0) = 0,  u'(1) + kappa u(1) = 0,

which is the stationary point of W/K with
    W = int_0^1 u'^2 dx + kappa u(1)^2   (kappa u(1)^2 plays the vacuum term)
    K = int_0^1 u^2 dx.
Discretized with N intervals: W with differences on the intervals (second
order), K with point values u_k and weights w_k = h, except w_N, which is
either h (2.47) or h/2 (now). The exact lambda solves sqrt(l) cos(sqrt(l)) +
kappa sin(sqrt(l)) = 0.

With u(1) = 0 instead (fixed boundary) the weight of the last point drops out:
both versions give the same result, converging like 1/N^2.
"""
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

from _common import BLUE, GREY, ORANGE, plt, save

KAPPA = 1.0


def exact_free():
    return brentq(lambda s: s * np.cos(s) + KAPPA * np.sin(s), 1.0, 3.0) ** 2


def discrete(N, w_end, free=True):
    h = 1.0 / N
    n = N if free else N - 1                       # unknowns u_1..u_N (free) or u_1..u_{N-1}
    A = np.zeros((n, n)); B = np.zeros((n, n))
    # stiffness: sum over intervals of ((u_k - u_{k-1})/h)^2 h, u_0 = 0
    for k in range(1, N + 1):
        i, j = k - 1, k - 2                        # indices of u_k and u_{k-1}
        for a, ca in ((i, 1.0), (j, -1.0)):
            for b, cb in ((i, 1.0), (j, -1.0)):
                if 0 <= a < n and 0 <= b < n:
                    A[a, b] += ca * cb / h
    if free:
        A[n - 1, n - 1] += KAPPA                   # the surface ("vacuum") term
    w = np.full(n, h)
    if free:
        w[-1] = w_end
    B[np.arange(n), np.arange(n)] = w
    return eigh(A, B, eigvals_only=True)[0]


lam = exact_free()
lam_fix = np.pi ** 2
Ns = np.array([25, 50, 100, 200, 400, 800])
err_old = np.array([abs(discrete(N, 1.0 / N) - lam) / lam for N in Ns])
err_new = np.array([abs(discrete(N, 0.5 / N) - lam) / lam for N in Ns])
err_fix = np.array([abs(discrete(N, 1.0 / N, free=False) - lam_fix) / lam_fix for N in Ns])

print(f"exact lambda (free, kappa = {KAPPA}) = {lam:.10f}")
print(f"{'N':>5} {'rel err, w_N = h (2.47)':>24} {'rel err, w_N = h/2 (now)':>25} {'fixed bnd':>11}")
for N, a, b, c in zip(Ns, err_old, err_new, err_fix):
    print(f"{N:5d} {a:24.3e} {b:25.3e} {c:11.3e}")
p_old = np.polyfit(np.log(Ns), np.log(err_old), 1)[0]
p_new = np.polyfit(np.log(Ns), np.log(err_new), 1)[0]
print(f"observed order: 2.47 {-p_old:.2f}, now {-p_new:.2f}")

fig, ax = plt.subplots(figsize=(5.6, 3.8))
ax.loglog(Ns, err_old, "o-", color=ORANGE, label="free end, weight h (2.47)")
ax.loglog(Ns, err_new, "s-", color=BLUE, label="free end, weight h/2 (now)")
ax.loglog(Ns, err_fix, "^--", color=GREY, label="fixed end (either weight)")
ax.set_xlabel("number of radial intervals N")
ax.set_ylabel("relative error of the lowest eigenvalue")
ax.set_title("model problem: the end weight sets the order")
ax.legend(fontsize=8)
save(fig, "edge_weight_1d.png")
