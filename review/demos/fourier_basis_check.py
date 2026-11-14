"""Area 3: the trigonometric basis of fourier2real (VMEC_utilities, commits 3b99109, 1a53f50).

At 2.47 calc_trigon_factors stored cos(m theta - n zeta) and sin(...) for
every mode of the MAIN list at every grid point, and fourier2real summed
coefficient * derivative factor * stored factor.

Now calc_trigon_factors stores the one-angle basis
    index m+1                       : cos(m theta), sin(m theta),   m = 0..mpol_tf-1
    index mpol_tf+ntor_tf+1+n       : cos(n nfp zeta), sin(n nfp zeta), n = -ntor_tf..ntor_tf
(grid points first) and fourier2real_1 forms, per mode,
    cos(m theta - n zeta) = cos(m theta) cos(n zeta) + sin(m theta) sin(n zeta)
    sin(m theta - n zeta) = sin(m theta) cos(n zeta) - cos(m theta) sin(n zeta)
for any mode list (main or Nyquist) within the basis.

This script copies the index arithmetic of the new fourier2real_1 line by line
(1-based indices kept as in the Fortran, shifted by -1 when used) and checks it
against direct evaluation of the derivatives of
    f = sum_k c_k cos(m_k theta - n_k zeta) + s_k sin(m_k theta - n_k zeta)
for random coefficients, the W7-X-like mode extent, and every combination of
up to two derivatives in each angle and of the sym flags.
It checks the algorithm, not the compiled Fortran: that is what the unchanged
eigenvalues on the three VMEC fixtures check (README, 2.49).
"""
import itertools

import numpy as np

rng = np.random.default_rng(1)
nfp, mpol_V, ntor_V = 5, 12, 12
mn_V = [(0, n * nfp) for n in range(0, ntor_V + 1)] + \
       [(m, n * nfp) for m in range(1, mpol_V) for n in range(-ntor_V, ntor_V + 1)]
mn_nyq = [(0, n * nfp) for n in range(0, 16 + 1)] + \
         [(m, n * nfp) for m in range(1, 19) for n in range(-16, 16 + 1)]
mpol_tf = max(mpol_V, max(m for m, _ in mn_nyq) + 1)                 # as read_VMEC
ntor_tf = max(ntor_V, max(abs(n) for _, n in mn_nyq) // nfp)

theta = rng.uniform(0, 2 * np.pi, (7, 3))
zeta = rng.uniform(0, 2 * np.pi, (7, 3))

# calc_trigon_factors: (n1, n2, mpol_tf+2*ntor_tf+1, 2), one normal point
tf = np.zeros(theta.shape + (mpol_tf + 2 * ntor_tf + 1, 2))
for m in range(mpol_tf):
    tf[..., m, 0], tf[..., m, 1] = np.cos(m * theta), np.sin(m * theta)
for n in range(-ntor_tf, ntor_tf + 1):
    tf[..., mpol_tf + ntor_tf + n, 0] = np.cos(n * nfp * zeta)
    tf[..., mpol_tf + ntor_tf + n, 1] = np.sin(n * nfp * zeta)


def fourier2real_1(c, s, mn, deriv, sym):
    """Line-by-line copy of the new Fortran (one normal point)."""
    i_c, i_s = (1, 2) if sum(deriv) % 2 == 0 else (2, 1)
    sym = [sym[0] and np.any(np.abs(c) > 0), sym[1] and np.any(np.abs(s) > 0)]
    need = [(sym[0] and i_c == 1) or (sym[1] and i_s == 1),
            (sym[0] and i_c == 2) or (sym[1] and i_s == 2)]
    varr = np.zeros(theta.shape)
    for k, (m, n) in enumerate(mn):
        i_m = m + 1
        i_n = mpol_tf + ntor_tf + 1 + n // nfp
        fac_c = m ** deriv[0] * (-n) ** deriv[1] * (-1) ** ((sum(deriv) + 1) // 2)
        fac_s = m ** deriv[0] * (-n) ** deriv[1] * (-1) ** (sum(deriv) // 2)
        if abs(c[k]) <= 0 and abs(s[k]) <= 0:
            continue
        cs = [None, None]
        if need[0]:
            cs[0] = tf[..., i_m - 1, 0] * tf[..., i_n - 1, 0] + tf[..., i_m - 1, 1] * tf[..., i_n - 1, 1]
        if need[1]:
            cs[1] = tf[..., i_m - 1, 1] * tf[..., i_n - 1, 0] - tf[..., i_m - 1, 0] * tf[..., i_n - 1, 1]
        if sym[0]:
            varr += c[k] * fac_c * cs[i_c - 1]
        if sym[1]:
            varr += s[k] * fac_s * cs[i_s - 1]
    return varr


def direct(c, s, mn, deriv, sym):
    """d^a/dtheta^a d^b/dzeta^b of the Fourier sum, written out."""
    a, b = deriv
    out = np.zeros(theta.shape)
    for k, (m, n) in enumerate(mn):
        ph = m * theta - n * zeta
        amp = m ** a * (-n) ** b                     # chain rule factor of (m theta - n zeta)
        p = a + b                                    # d^p/dph^p of cos and sin
        dcos = [np.cos(ph), -np.sin(ph), -np.cos(ph), np.sin(ph)][p % 4]
        dsin = [np.sin(ph), np.cos(ph), -np.sin(ph), -np.cos(ph)][p % 4]
        if sym[0]:
            out += c[k] * amp * dcos
        if sym[1]:
            out += s[k] * amp * dsin
    return out


worst = 0.0
for mn, label in ((mn_V, "main list"), (mn_nyq, "Nyquist list")):
    c = rng.normal(size=len(mn)) / (1 + np.arange(len(mn)))
    s = rng.normal(size=len(mn)) / (1 + np.arange(len(mn)))
    for deriv in itertools.product(range(3), range(3)):
        for sym in ((True, True), (True, False), (False, True)):
            new = fourier2real_1(c, s, mn, deriv, sym)
            ref = direct(c, s, mn, deriv, sym)
            err = np.max(np.abs(new - ref)) / max(np.max(np.abs(ref)), 1e-300)
            worst = max(worst, err)
    print(f"{label:13}: {len(mn)} modes, max relative difference over 27 cases so far {worst:.1e}")
print(f"mpol_tf = {mpol_tf}, ntor_tf = {ntor_tf}: basis of {mpol_tf + 2 * ntor_tf + 1} functions "
      f"instead of {len(mn_V)} (main) or {len(mn_nyq)} (Nyquist) per grid point")
assert worst < 1e-10, "basis form disagrees with direct evaluation"
print("PASS")
