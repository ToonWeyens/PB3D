"""Area 4c: the singular half-cell integral of the field-line G (vac_utilities.calc_GH_int_1).

Near the influence point the squared distance on the surface is the quadratic
form g_tt dp^2 + 2 g_at dp da + g_aa da^2 (p along the field line, a the line
label). With the metric-normalised half-cell sizes
    U = sqrt(g_tt) dp/2,   V = sqrt(g_aa) da/2,   s = g_at / sqrt(g_aa g_tt),
the code uses the primitive
    F = int_0^U int_-V^V du dv / sqrt(u^2 + 2 s u v + v^2)
      = U ln((V + sU + D+)/(-V + sU + D-)) + V ln((U + sV + D+)/(-U + sV + D-)),
    D+- = sqrt(U^2 +- 2 s U V + V^2).
This script checks the closed form against brute-force quadrature in polar
coordinates around the singular corner (which removes the 1/r singularity),
for perpendicular and strongly skewed cells and very different aspect ratios.

It also evaluates the bracket of the 2.47 expression,
    ln|(d2 + E + s)/(d1 - E + s)| - E/s^2 ln|(d2 + E + s)/(d1 - E + s)|,
with E = V/U (the same logarithm twice), to show that it diverges like 1/s^2
for a cell with perpendicular sides (s -> 0), where the integral is finite.
"""
import numpy as np
from scipy.integrate import quad

from _common import BLUE, ORANGE, plt, save


def F_closed(U, V, s):
    Dp = np.sqrt(U**2 + 2 * s * U * V + V**2)
    Dm = np.sqrt(U**2 - 2 * s * U * V + V**2)
    return U * np.log(abs((V + s * U + Dp) / (-V + s * U + Dm))) + \
        V * np.log(abs((U + s * V + Dp) / (-U + s * V + Dm)))


def F_brute(U, V, s):
    """Rectangle [0,U] x [-V,V] split into two triangles-fans from the origin (polar)."""
    def over(phi_a, phi_b, rmax):
        def inner(phi):
            c, d = np.cos(phi), np.sin(phi)
            q = np.sqrt(c * c + 2 * s * c * d + d * d)        # r^2 (c,d) -> distance factor
            return rmax(phi) / q                              # int_0^rmax r dr / (r q) = rmax / q
        return quad(inner, phi_a, phi_b, epsabs=1e-13, epsrel=1e-13, limit=200)[0]
    pc = np.arctan2(V, U)                                     # corner angles
    tot = 0.0
    tot += over(-pc, pc, lambda p: U / np.cos(p))             # right edge u = U
    tot += over(pc, np.pi / 2, lambda p: V / np.sin(p))       # top edge v = V
    tot += over(-np.pi / 2, -pc, lambda p: V / -np.sin(p))    # bottom edge v = -V
    return tot


def old_bracket(U, V, s):
    E = V / U
    d1 = np.sqrt(1 + 2 * E * s + E**2)
    d2 = np.sqrt(1 - 2 * E * s + E**2)
    L = np.log(abs((d2 + E + s) / (d1 - E + s)))
    return L - E / s**2 * L


print(f"{'U':>6} {'V':>6} {'s':>6} {'closed form':>14} {'brute force':>14} {'rel diff':>9}")
for U, V, s in [(1, 1, 0.0), (1, 1, 0.5), (1, 1, -0.9), (0.1, 1, 0.3), (1, 0.05, 0.95), (2.5, 0.4, -0.6)]:
    a, b = F_closed(U, V, s), F_brute(U, V, s)
    print(f"{U:6.2f} {V:6.2f} {s:6.2f} {a:14.10f} {b:14.10f} {abs(a - b) / b:9.1e}")

ss = np.linspace(-0.95, 0.95, 400)
ss = ss[np.abs(ss) > 0.01]
fig, ax = plt.subplots(figsize=(5.6, 3.6))
ax.plot(ss, [F_closed(1, 1, s) for s in ss], color=BLUE, label="F, closed form (now)")
ax.plot(ss, [old_bracket(1, 1, s) for s in ss], color=ORANGE, label="2.47 bracket")
ax.set_ylim(-20, 20)
ax.set_xlabel("skew s = g_at / sqrt(g_aa g_tt)")
ax.set_title("half-cell integral for U = V = 1")
ax.legend(fontsize=8)
save(fig, "halfcell_integral.png")
