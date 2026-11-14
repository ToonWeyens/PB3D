"""Area 4c: why the field-line kernels need a near-field rule (commit 0bbda47).

On W7-X the field lines at the inboard bean tips run almost toroidally, so
neighbouring lines lie much closer together (median 0.007 R_0 at 128 lines)
than the points along each line (0.06 R_0). Seen from a point on one line,
the kernel along the neighbouring line is a peak of width ~d (the line
distance) sampled at spacing h >> d.

Model: a straight neighbouring line at distance d from the influence point,
points x_k = (k + f) h, k = -K..K, where f in [0, 1) is the offset of the
points relative to the closest approach (on W7-X, set by frac(q dtheta /
dalpha) of the grid). Integrals over |x| <= L:

    G-like:  int dx / sqrt(d^2 + x^2)          = 2 asinh(L/d)
    H-like:  int d dx / (d^2 + x^2)^(3/2)      = 2 L / (d sqrt(d^2 + L^2))

Rules: (1) the point rule (trapezoidal along the line, as calc_GH_int_1 at
2.47 and now far from the influence point); (2) the near-field rule of
calc_GH_int_1 now: within n_near = 4 interval lengths, 8-point Gauss-Legendre
on each half interval (here on the exact straight line; PB3D uses the cubic
interpolant through four points, exact for a straight line).
The error of the point rule is of order one and its sign depends on f: this
is how every W7-X grid got sorted into too large or plausible responses.
"""
import numpy as np

from _common import BLUE, GREEN, ORANGE, RED, plt, save

xg, wg = np.polynomial.legendre.leggauss(8)
h, K = 1.0, 40
L = K * h


def rules(d, f):
    xk = (np.arange(-K, K + 1) + f) * h
    xk = xk[np.abs(xk) <= L]
    exact_G = 2 * np.arcsinh(L / d)
    exact_H = 2 * L / (d * np.sqrt(d**2 + L**2))
    gk = 1 / np.sqrt(d**2 + xk**2)
    hk = d / (d**2 + xk**2) ** 1.5
    # point rule (trapezoid) on the sampled part + exact tails to |x| = L
    def tails(fun_int):
        return fun_int(L) - fun_int(xk[-1]) + fun_int(-xk[0]) - fun_int(L)
    G_int = lambda x: np.arcsinh(x / d)                     # noqa: E731
    H_int = lambda x: x / (d * np.sqrt(d**2 + x**2))        # noqa: E731
    w = np.full(xk.size, h); w[0] = w[-1] = h / 2
    pt_G = np.sum(w * gk) + (G_int(L) - G_int(xk[-1])) + (G_int(xk[0]) - G_int(-L))
    pt_H = np.sum(w * hk) + (H_int(L) - H_int(xk[-1])) + (H_int(xk[0]) - H_int(-L))
    # near-field rule: Gauss on the intervals within 4 h of the closest approach
    nf_G, nf_H = pt_G, pt_H
    for i in range(xk.size - 1):
        a, b = xk[i], xk[i + 1]
        if min(abs(a), abs(b)) < 4 * h:
            # remove the trapezoid of this interval, add Gauss on its two halves
            nf_G -= h / 2 * (gk[i] + gk[i + 1]); nf_H -= h / 2 * (hk[i] + hk[i + 1])
            for lo, hi in ((a, 0.5 * (a + b)), (0.5 * (a + b), b)):
                x = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
                nf_G += 0.5 * (hi - lo) * np.sum(wg / np.sqrt(d**2 + x**2))
                nf_H += 0.5 * (hi - lo) * np.sum(wg * d / (d**2 + x**2) ** 1.5)
    return [(pt_G - exact_G) / exact_G, (pt_H - exact_H) / exact_H,
            (nf_G - exact_G) / exact_G, (nf_H - exact_H) / exact_H]


fs = np.linspace(0, 1, 101)[:-1]
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6), sharex=True)
print(f"{'d/h':>5} {'point rule G: min..max':>24} {'point rule H: min..max':>24} {'near rule H: max':>17}")
for d, c in ((0.12, RED), (0.25, ORANGE), (0.5, GREEN), (1.0, BLUE)):
    e = np.array([rules(d, f) for f in fs])
    axes[0].plot(fs, 100 * e[:, 1], color=c, label=f"d/h = {d}")
    axes[1].semilogy(fs, np.abs(e[:, 3]) + 1e-17, color=c, label=f"d/h = {d}")
    print(f"{d:5.2f} {100*e[:,0].min():10.2f} .. {100*e[:,0].max():6.2f} %   "
          f"{100*e[:,1].min():8.1f} .. {100*e[:,1].max():8.1f} %   {np.abs(e[:,3]).max():14.1e}")
axes[0].set_title("point rule, H-like kernel: error in %")
axes[0].set_xlabel("offset f of the points (fraction of h)")
axes[0].set_ylim(-110, 400)
axes[1].set_title("near-field rule (8-pt Gauss within 4 h): |relative error|")
axes[1].set_xlabel("offset f of the points (fraction of h)")
axes[0].legend(fontsize=8); axes[1].legend(fontsize=8)
save(fig, "near_field_lines.png")
print("W7-X bean tips: d/h ~ 0.007/0.06 = 0.12 (median at 128 lines).")
