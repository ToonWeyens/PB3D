"""Where the 1/N error of the cbm18a fixed-boundary eigenvalue comes from.

The numbers are runs of PB3D at f9713d7 (deck tests/regression/input_cbm18a_fixed_bnd,
tol_SLEPC 1e-12, one process); see section 12 of plan_radial_finite_elements.md.
Eigenvalues are omega^2/omega_A^2 in units of 1e-2; N is n_r_sol; p is
norm_disc_prec_sol.

The differences displace the two ends of the plasma (fe_model.py, section 'ends'):
the first p points are held at zero, so the plasma starts p - 1/2 steps further
out; at a fixed edge with all weights 1 it ends 0.2805 (p = 3) or 0.375 (p = 2)
steps further in. This script multiplies these displacements by the measured
sensitivity of the eigenvalue to the position of each end, and compares with the
measured dependence on N.

    python -B cbm18a_ends.py
"""
import numpy as np

N = np.array([200, 400, 800, 1600, 3200])
SERIES = {  # (set of Fourier modes, p): eigenvalue at the five N
    ('sliding set of 20', 3): [-4.357290, -4.365940, -4.372961, -4.376516, -4.378267],
    ('fixed set m = 5..40', 3): [-4.350099, -4.358696, -4.365580, -4.369078, -4.370808],
    ('sliding set of 20', 2): [-4.588865, -4.428302, -4.389857, -4.381351, -4.379781],
    ('fixed set m = 5..40', 2): [-4.582066, -4.421165, -4.382503, -4.373920, -4.372324],
    ('fixed set m = 5..47', 3): [np.nan, np.nan, -4.374661, -4.378223, -4.379985],
}
# fixed set m = 5..40, p = 3, 800 points: one end moved
PSI_IN, LAM_IN = np.array([0.450, 0.452, 0.454, 0.458]), np.array([-4.365580, -4.350905, -4.335130, -4.299949])
PSI_E, LAM_E = np.array([0.750, 0.748, 0.746, 0.742]), np.array([-4.365580, -4.357958, -4.349602, -4.330580])
STEPS_IN = {3: 2.5, 2: 1.5}
STEPS_E = {3: 0.2805, 2: 0.375}


def fit(lam, p):
    """Coefficient of 1/N and limit: two points for p = 3, three points with a 1/N^2 term for p = 2."""
    lam = np.array(lam)
    if p == 3:
        a = (lam[3] - lam[4])/(1/N[3] - 1/N[4])
        return a, lam[4] - a/N[4]
    c = np.linalg.solve(np.array([[1, 1/n, 1/n**2] for n in N[2:]]), lam[2:])
    return c[1], c[0]


print('Measured: lambda(N) = limit + a/N (+ b/N^2 for p = 2)')
for (name, p), lam in SERIES.items():
    a, lim = fit(lam, p)
    print('  %-20s p = %d:  a = %.4f   limit = %.5f' % (name, p, a/100, lim))

# sensitivities at the nominal ends; every 800-point run has its ends displaced by the steps above
d_in = (0.75 - PSI_IN)/799
f_in = np.poly1d(np.polyfit(PSI_IN + 2.5*d_in - 0.45, LAM_IN, 3))
d_e = (PSI_E - 0.45)/799
lam_e = LAM_E - (f_in(2.5*d_e) - f_in(2.5*d_e[0]))        # without the small move of the inner end
f_e = np.poly1d(np.polyfit(0.75 - (PSI_E - 0.2805*d_e), lam_e, 3))
s_in, s_e = f_in.deriv()(0.0), f_e.deriv()(0.0)
print('\nSensitivity of the eigenvalue (in 1e-2) per unit of normalized flux, fixed set m = 5..40:')
print('  inner end moved outwards: %.2f     fixed edge moved inwards: %.2f' % (s_in, s_e))

print('\nPredicted coefficient of 1/N = sensitivity x steps x (range 0.3):')
for p in (3, 2):
    a_in, a_e = s_in*STEPS_IN[p]*0.3/100, s_e*STEPS_E[p]*0.3/100
    lam = lambda n: f_in(STEPS_IN[p]*0.3/(n - 1)) + f_e(STEPS_E[p]*0.3/(n - 1))
    if p == 3:
        a_fit = (lam(1600) - lam(3200))/(1/1600 - 1/3200)/100
    else:
        a_fit = np.linalg.solve(np.array([[1, 1/n, 1/n**2] for n in N[2:]]), [lam(n) for n in N[2:]])[1]/100
    print('  p = %d: inner end %.4f (%.0f %%), fixed edge %.4f (%.0f %%); the same fit applied to this prediction '
          'gives %.4f; measured %.4f'
          % (p, a_in, 100*a_in/(a_in + a_e), a_e, 100*a_e/(a_in + a_e), a_fit, fit(SERIES[('fixed set m = 5..40', p)], p)[0]/100))
