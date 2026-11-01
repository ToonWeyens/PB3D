"""The pedestal-foot flux surface of the extended cbm18a: geometry, Fourier radius, targets for the cut equilibrium."""
import sys, numpy as np
from scipy.interpolate import CubicSpline
sys.path.insert(0, '.'); from parse_hel12 import parse
d = parse(sys.argv[1]); psi_c = float(sys.argv[2]); s_c = np.sqrt(psi_c)
cs, qs, p0, F = d['cs'], d['qs'], d['p0'], d['rbphi']
a_old = d['radius'] * d['raxis']; eps_old = d['eps']; R0_old = a_old / eps_old
# surfaces 2..nr hold xout/yout; interpolate in s
x = np.array([CubicSpline(cs[1:], d['xout'][:, j])(s_c) for j in range(d['nchi'])])
y = np.array([CubicSpline(cs[1:], d['yout'][:, j])(s_c) for j in range(d['nchi'])])
R = a_old * (1 / eps_old + x); Z = a_old * y
# close the curve by symmetry (ias = 0: chi from 0 to pi)
Rf = np.concatenate([R, R[-2:0:-1]]); Zf = np.concatenate([Z, -Z[-2:0:-1]])
R0 = 0.5 * (R.max() + R.min()); a = 0.5 * (R.max() - R.min()); eps = a / R0
th = np.arctan2(Zf / a, (Rf - R0) / a); rho = np.hypot((Rf - R0) / a, Zf / a)
o = np.argsort(np.mod(th, 2 * np.pi)); th, rho = np.mod(th, 2 * np.pi)[o], rho[o]
thu = np.linspace(0, 2 * np.pi, 4096, endpoint=False); rhou = np.interp(thu, np.concatenate([th, th[:1] + 2 * np.pi]), np.concatenate([rho, rho[:1]]))
c = np.fft.rfft(rhou) / len(thu); am = 2 * c.real  # rho = am[0]/2 + sum am[m] cos(m th) (sine parts vanish by symmetry)
print('psi_c %.6f s_c %.6f | old: R0 %.6f a %.6f eps %.6f | new: R0 %.6f a %.6f eps %.6f Zmax/a %.5f' % (psi_c, s_c, R0_old, a_old, eps_old, R0, a, eps, Z.max() / a))
print('rho Fourier (am[0]/2 + sum am[m] cos m theta):', ' '.join('%.3e' % v for v in am[:16]))
print('max |sin parts| %.2e' % np.abs(c.imag).max(), '| reconstruction error with 24 modes: %.2e' % np.abs(rhou - (am[0] / 2 + sum(am[m] * np.cos(m * thu) for m in range(1, 24)))).max())
q_c = CubicSpline(cs, qs)(s_c); F_c = CubicSpline(cs, F)(s_c); p_c = CubicSpline(cs, p0)(s_c)
print('targets: q0 %.5f  q(psi_c) %.5f | F(psi_c) %.6f (axis %.6f) | p0 axis %.6g p(psi_c) %.3g | cpsurf*psi_c %.6f' % (qs[0], q_c, F_c, F[0], p0[0], p_c, d['cpsurf'] * psi_c))
np.save('foot_am.npy', am); np.save('foot_targets.npy', np.array([qs[0], q_c, eps, R0, a, F_c, p0[0]]))
# the q profile inside psi_c against the new normalised flux, for the later comparison
import json; json.dump({'psi_new': (cs**2 / psi_c)[cs**2 <= psi_c].tolist(), 'q': qs[cs**2 <= psi_c].tolist(), 'p_over_p0': (p0 / p0[0])[cs**2 <= psi_c].tolist()}, open('extended_inside.json', 'w'))
