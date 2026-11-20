"""Area 2a: the HELENA shear and the end condition of an odd spline.

PB3D computes the magnetic shear of a HELENA equilibrium in
eq_ops.calc_derived_S_from_deriv_HEL as

    S = -(1/J) d/dchi (h^{alpha psi} / h^{psi psi}),

by fitting a cubic spline to h12/h22 over the poloidal grid chi_H and
differentiating it. For an up-down symmetric equilibrium (ias = 0) the grid
covers only the half period [0, pi], so the spline needs end conditions at
chi = 0 and chi = pi. h12/h22 is odd in chi. An odd function has a non-zero
slope at 0 and pi and a zero second derivative there.

At 2.47 the call passed bcs(:,1) = [1,1] (PSPLINE code 1: prescribed first
derivative, value 0), the condition for EVEN quantities. Commit 28f8dab passes
bcs(:,2) = [2,2] (code 2: prescribed second derivative, value 0).

This script repeats that with scipy's CubicSpline on a model odd function
and shows the error of the spline derivative with both end conditions.
"""
import numpy as np
from scipy.interpolate import CubicSpline

from _common import BLUE, ORANGE, plt, save


def model(chi):
    """An odd model of h12/h22 with a few harmonics (shape is not important)."""
    return 0.6 * np.sin(chi) + 0.25 * np.sin(2 * chi) - 0.08 * np.sin(3 * chi)


def dmodel(chi):
    return 0.6 * np.cos(chi) + 0.5 * np.cos(2 * chi) - 0.24 * np.cos(3 * chi)


fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
print(f"{'nchi':>6} {'max err, even BC (2.47)':>26} {'max err, odd BC (now)':>24}")
for nchi in (33, 65, 129, 257):
    chi = np.linspace(0, np.pi, nchi)
    f = model(chi)
    wrong = CubicSpline(chi, f, bc_type=((1, 0.0), (1, 0.0)))   # PSPLINE bcs code 1, value 0
    right = CubicSpline(chi, f, bc_type=((2, 0.0), (2, 0.0)))   # PSPLINE bcs code 2, value 0
    fine = np.linspace(0, np.pi, 2001)
    ew = np.max(np.abs(wrong(fine, 1) - dmodel(fine)))
    er = np.max(np.abs(right(fine, 1) - dmodel(fine)))
    print(f"{nchi:6d} {ew:26.3e} {er:24.3e}")
    if nchi == 65:
        axes[0].plot(fine, dmodel(fine), color="k", lw=1, label="exact d/dchi")
        axes[0].plot(fine, wrong(fine, 1), color=ORANGE, label="spline, zero 1st derivative (2.47)")
        axes[0].plot(fine, right(fine, 1), color=BLUE, ls="--", label="spline, zero 2nd derivative (now)")
        axes[1].semilogy(chi, np.abs(wrong(chi, 1) - dmodel(chi)) + 1e-16, color=ORANGE, label="2.47")
        axes[1].semilogy(chi, np.abs(right(chi, 1) - dmodel(chi)) + 1e-16, color=BLUE, label="now")

axes[0].set_xlabel("chi"); axes[0].set_title("derivative of an odd function on [0, pi], 65 points")
axes[0].legend(fontsize=8)
axes[1].set_xlabel("chi"); axes[1].set_title("|error| of the derivative at the grid points")
axes[1].legend(fontsize=8)
for ax in axes:
    ax.set_xticks([0, np.pi / 2, np.pi]); ax.set_xticklabels(["0", "pi/2", "pi"])
save(fig, "shear_odd_spline.png")
print("The wrong end condition forces the slope to 0 at chi = 0 and pi, where the")
print("odd function has its largest slope: an O(1) error that does not shrink with")
print("nchi at the ends and decays over a few intervals inward.")
