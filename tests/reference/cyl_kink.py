"""T1/T2: external kink of a uniform-current, zero-pressure, large-aspect-ratio
cylinder (Shafranov), derived from the reduced-MHD energy principle with the
exterior Laplace solution written out, for a conformal wall at radius b.

Conventions: xi_r = xi(r) exp(i(m theta - n z/R0)); F = k.B = (B_theta/r)(m - n q);
q = r B_z/(R0 B_theta) = q_a (uniform current, B_z uniform); rho uniform;
omega_A^2 = B_z^2/(mu0 rho R0^2); omega^2 = -gamma^2 (PB3D/MISHKA sign).
"""
import sympy as sp

r, a, b, R0, Bz, Bta, mu0, rho, xa = sp.symbols('r a b R_0 B_z B_theta_a mu_0 rho xi_a', positive=True)
m = sp.symbols('m', positive=True, integer=True)
n, x = sp.symbols('n x', real=True)             # x = n q_a
qa = a*Bz/(R0*Bta)

# --- plasma: reduced-MHD fluid energy for the straight tokamak (Newcomb's f, g in the
#     limit k^2 r^2 << m^2, beta -> 0). F is constant for uniform q.
F = Bta/a*(m - x)                                # k.B, uniform q_a = a Bz/(R0 Bta)
xi = xa*(r/a)**(m-1)                             # exact interior eigenfunction (harmonic)
integrand = (r**3*sp.diff(xi, r)**2 + (m**2-1)*r*xi**2)*F**2/m**2
dW_F = sp.pi**2*R0/mu0*sp.integrate(integrand, (r, 0, a))
dW_F = sp.simplify(dW_F)

# --- vacuum: chi harmonic for r > a, d chi/dr = 0 at r = b; B_1r = i F xi at r = a.
C = sp.symbols('C')
chi = C*((r/a)**(-m) + (a/b)**(2*m)*(r/a)**m)
Csol = sp.solve(sp.Eq(sp.diff(chi, r).subs(r, a), F*xa), C)[0]  # |B_1r| = F xi_a
chi = chi.subs(C, Csol)
# W_V = (1/2mu0) int |grad chi|^2 dV = (1/2mu0) oint chi (-d chi/dr) dS, with the theta average 1/2
dW_V = sp.simplify(sp.Rational(1, 2)/mu0*2*sp.pi*R0*2*sp.pi*a*sp.Rational(1, 2)*(-chi*sp.diff(chi, r)).subs(r, a))
lam = (1 + (a/b)**(2*m))/(1 - (a/b)**(2*m))
assert sp.simplify(dW_V - sp.pi**2*R0*a**2*lam/(mu0*m)*F**2*xa**2) == 0

# --- surface term of the current discontinuity at r = a (Wesson 6.9.3, Freidberg 11.5):
#     (2/q_a)(n/m - 1/q_a) + (n/m - 1/q_a)^2, in units pi^2 Bz^2 a^2 xi_a^2/(mu0 R0)
nm_1q = (x - m)/(m*qa)                           # n/m - 1/q_a
dW_S = sp.pi**2*Bz**2*a**2*xa**2/(mu0*R0)*(2/qa*nm_1q + nm_1q**2)

dW = sp.simplify(dW_F + dW_V + dW_S)
# kinetic energy: incompressible, |xi_theta| = |xi_r| for the r^(m-1) eigenfunction
K = sp.simplify(sp.Rational(1, 2)*rho*2*sp.pi*R0*2*sp.pi*sp.Rational(1, 2)*sp.integrate(2*xi**2*r, (r, 0, a)))
gamma2 = sp.simplify(-dW/K)
omegaA2 = Bz**2/(mu0*rho*R0**2)
w2 = sp.simplify(-gamma2/omegaA2)                # omega^2/omega_A^2
w2 = sp.simplify(w2.subs(Bta, a*Bz/(R0*sp.Symbol('q_a', positive=True))))
print('omega^2/omega_A^2 (wall at b) =', sp.factor(w2))
w2_nowall = sp.simplify(sp.limit(w2, b, sp.oo))
print('omega^2/omega_A^2 (no wall)   =', sp.factor(w2_nowall))

# --- numbers for the HELENA cases
qa_sym = sp.Symbol('q_a', positive=True)
f_nowall = sp.lambdify((x, m, qa_sym), w2_nowall)
f_wall = sp.lambdify((x, m, qa_sym, a, b), w2)
for label, q in [('eps=0.05 (q_a=1.0378)', 1.0378), ('eps=0.33 (q_a=1.0382)', 1.0382)]:
    print(f'\n{label}: no wall, mode m = ceil(n q_a); wall factor for b/a = 3, 2, 1.5, 1.2')
    for nn in [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]:
        xx = nn*q
        mm = int(xx) + 1
        row = [f_nowall(xx, mm, q)] + [f_wall(xx, mm, q, 1.0, bb) for bb in (3, 2, 1.5, 1.2)]
        print(f'  n={nn:2d} m={mm:2d} nq_a={xx:7.4f}  ' + '  '.join(f'{v:+.5f}' for v in row))
