"""Areas 4b and 4e: the vacuum-side sign convention and the indirect diagonal, on a toy.

A boundary-element toy on a circular torus (R0 = 3, a = 1) with the SAME
conventions as PB3D's vac_ops now (both styles):

  * stored normal  n = x_theta x x_phi = -a (R0 + a cos theta) e_r:
    it points INTO the plasma, i.e. OUT OF the vacuum, and carries the area
    element (PB3D: -J grad psi, which carries |J grad psi| for dS = |J grad psi|
    dtheta dalpha);
  * normal data    dphi_j = n_j . grad phi(x_j);
  * G_ij = -w / |x_i - x_j|,   H_ij = w n_j . (x_j - x_i) / |x_j - x_i|^3,
    w = dtheta dphi (plain trapezoidal rule on a doubly periodic grid);
  * diagonal of G: the integral of -1/r over the singular cell (flat-cell
    closed form, as calc_GH_int_1's half cells);
  * diagonals, three variants:
      "row sum + G lin": H_ii from the row sum H 1 = -4 pi and a correction of
          G_ii from the identities of the linear harmonics x, y, z projected on
          the normal (PB3D's default for the field-line style, vac_diag_style 1);
      "row sum": the row sum for H_ii, G_ii from the flat cell only;
      "-2 pi": the corner term of a smooth surface alone for H_ii (a direct
          diagonal WITHOUT the cell integral of the H kernel that PB3D's
          vac_diag_style 2 adds), G_ii from the flat cell.

Identities that the assembled operators should satisfy on the boundary:
    interior harmonic phi (1, x, z):          H phi - G dphi = -4 pi phi
    exterior decaying phi (point source on
    the magnetic axis):                       H phi - G dphi = 0
and the exterior problem H Phi = G dphi returns phi on the boundary (no shift
of the diagonal: PB3D's solve_Phi_BEM since 5 November 2026).

Then the same with the normal flipped (n -> -n, outward), as store_vac_VMEC
had it until commit e3b4847 while the kernels assumed the inward one.

This is a toy for understanding: plain trapezoidal rule on a (theta, phi)
grid, not PB3D's field-line grid with its near-field rule.
"""
import numpy as np

from _common import BLUE, GREY, ORANGE, RED, plt, save

R0, A = 3.0, 1.0
X0 = np.array([R0, 0.0, 0.0])            # point source on the magnetic axis (inside the plasma)


def F_closed(U, V, s):
    Dp = np.sqrt(U**2 + 2 * s * U * V + V**2)
    Dm = np.sqrt(U**2 - 2 * s * U * V + V**2)
    return U * np.log(abs((V + s * U + Dp) / (-V + s * U + Dm))) + \
        V * np.log(abs((U + s * V + Dp) / (-U + s * V + Dm)))


def surface(nt):
    nphi = 3 * nt                                             # roughly square cells
    th, ph = np.meshgrid(2 * np.pi * np.arange(nt) / nt, 2 * np.pi * np.arange(nphi) / nphi, indexing="ij")
    th, ph = th.ravel(), ph.ravel()
    rho = R0 + A * np.cos(th)
    x = np.stack([rho * np.cos(ph), rho * np.sin(ph), A * np.sin(th)], 1)
    x_t = np.stack([-A * np.sin(th) * np.cos(ph), -A * np.sin(th) * np.sin(ph), A * np.cos(th)], 1)
    x_p = np.stack([-rho * np.sin(ph), rho * np.cos(ph), 0 * th], 1)
    n = np.cross(x_t, x_p)                                    # = -A rho e_r: into the plasma
    dt, dp = 2 * np.pi / nt, 2 * np.pi / nphi
    g_tt, g_pp = A**2 * np.ones_like(th), rho**2
    G_diag = -2 * F_closed(np.sqrt(g_tt) * dt / 2, np.sqrt(g_pp) * dp / 2, 0.0) / np.sqrt(g_tt * g_pp)
    return x, n, dt * dp, G_diag


def operators(x, n, w, G_diag, diag):
    d = x[None, :, :] - x[:, None, :]                         # x_j - x_i
    r = np.sqrt(np.sum(d**2, 2)); np.fill_diagonal(r, 1.0)
    G = -w / r
    H = w * np.einsum("jk,ijk->ij", n, d) / r**3
    np.fill_diagonal(G, G_diag); np.fill_diagonal(H, 0.0)
    if diag.startswith("row sum"):
        np.fill_diagonal(H, -4 * np.pi - H.sum(1))
    else:
        np.fill_diagonal(H, -2 * np.pi)
    if diag == "row sum + G lin":
        # calc_GH_1 (vac_diag_style 1): add to G_ii what makes the identity
        # H x - G n_x = -4 pi x hold for x, y, z, projected on the normal
        lin = np.zeros(len(x))
        for k in range(3):
            lin += n[:, k] * (H @ x[:, k] + 4 * np.pi * x[:, k] - G @ n[:, k])
        G[np.arange(len(x)), np.arange(len(x))] += lin / np.sum(n**2, 1)
    return G, H


def potentials(x, n):
    """(phi, dphi, free term c) for interior 1, x, z and the exterior point source."""
    out = [(np.ones(len(x)), np.zeros(len(x)), -4 * np.pi, "1 (interior)"),
           (x[:, 0], n[:, 0], -4 * np.pi, "x (interior)"),
           (x[:, 2], n[:, 2], -4 * np.pi, "z (interior)")]
    dd = x - X0; r = np.linalg.norm(dd, axis=1)
    out.append((1 / r, np.sum(n * (-dd / r[:, None]**3), 1), 0.0, "1/|x-x0| (exterior)"))
    return out


rows = []
nts = (8, 12, 16, 24, 32)
for nt in nts:
    x, n_in, w, Gd = surface(nt)
    for normal, n in (("inward (now)", n_in), ("outward (VMEC before e3b4847)", -n_in)):
        for diag in ("row sum + G lin", "row sum", "-2 pi"):
            if normal.startswith("outward") and diag != "row sum + G lin":
                continue
            G, H = operators(x, n, w, Gd, diag)
            res = {}
            for phi, dphi, c, name in potentials(x, n):
                res[name] = np.max(np.abs(H @ phi - G @ dphi - c * phi)) / np.max(np.abs(4 * np.pi * phi))
            phi, dphi, _, _ = potentials(x, n)[-1]
            Phi = np.linalg.solve(H, G @ dphi)
            res["exterior solve"] = np.max(np.abs(Phi - phi)) / np.max(np.abs(phi))
            rows.append((nt, len(x), normal, diag, res))

names = ["1 (interior)", "x (interior)", "z (interior)", "1/|x-x0| (exterior)", "exterior solve"]
print("max residual of the identity / max|4 pi phi|, and the error of the exterior solve")
print(f"{'n_theta':>7} {'N':>5} {'normal':30} {'diagonals':16} " + " ".join(f"{s:>20}" for s in names))
for nt, N, normal, diag, res in rows:
    print(f"{nt:7d} {N:5d} {normal:30} {diag:16} " + " ".join(f"{res[s]:20.2e}" for s in names))

fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.7))
for (normal, diag, c, mk) in (("inward (now)", "row sum + G lin", BLUE, "o"),
                              ("inward (now)", "row sum", GREY, "d"),
                              ("inward (now)", "-2 pi", ORANGE, "s"),
                              ("outward (VMEC before e3b4847)", "row sum + G lin", RED, "^")):
    sel = [r for r in rows if r[2] == normal and r[3] == diag]
    N = [np.sqrt(r[1] / 3) for r in sel]                       # = n_theta
    axes[0].loglog(N, [r[4]["1/|x-x0| (exterior)"] for r in sel], mk + "-", color=c, label=f"{normal}, H diag {diag}")
    axes[1].loglog(N, [r[4]["exterior solve"] for r in sel], mk + "-", color=c, label=f"{normal}, H diag {diag}")
for ax in axes:
    ax.loglog(nts, 0.5 * (nts[0] / np.array(nts)), ":", color=GREY, label="first order")
    ax.loglog(nts, 0.5 * (nts[0] / np.array(nts))**2, "--", color=GREY, label="second order")
    ax.set_xlabel("points per poloidal turn (toroidal: 3x)")
axes[0].set_title("exterior identity H phi - G dphi = 0: residual")
axes[1].set_title("exterior solve H Phi = G dphi: error of Phi")
axes[1].legend(fontsize=7)
save(fig, "torus_bem_toy.png")
