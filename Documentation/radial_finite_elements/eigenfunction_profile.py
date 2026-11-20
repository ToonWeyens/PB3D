"""Shape of the lowest eigenfunction of a PB3D run along the radius, from PB3D_out.h5.

Prints the length over which the eigenfunction changes,

    sqrt( int sum_m |X_m|^2 dpsi / int sum_m |X_m'|^2 dpsi ),

as a fraction of the radial range of the run, the number of radial points per
such length, and the size and slope of the eigenfunction near the two ends.
Use it on a run with a fixed set of Fourier modes (X_style 1), where a slot of
the solution vector holds the same mode at every point. PB3D does not store
the eigenvector of an eigenvalue that it marks as complex, so tol_SLEPC must
be 1e-9 or larger.

    python -B eigenfunction_profile.py <path to PB3D_out.h5>
"""
import sys

import h5py
import numpy as np

f = h5py.File(sys.argv[1], 'r')
g = f['sol_R_1']
lim = g['RE_sol_vec/lim'][:]
shape = [lim[len(lim)//2 + i] - lim[i] + 1 for i in range(len(lim)//2)]
X = (g['RE_sol_vec/var'][:] + 1j*g['IM_sol_vec/var'][:]).reshape(shape, order='F')[:, :, 0]
M, N = X.shape
amp = np.sqrt(np.sum(abs(X)**2, axis=0))
X, amp = X/amp.max(), amp/amp.max()
t = np.arange(N)/(N - 1)                                    # fraction of the radial range
dX = np.diff(X, axis=1)*(N - 1)
slope = np.sqrt(np.sum(abs(dX)**2, axis=0))
length = np.sqrt(np.sum(abs(X)**2)/np.sum(abs(dX)**2))
print('%d Fourier modes, %d radial points; |X| is the root of the sum over the modes, scaled to a largest value of 1'
      % (M, N))
print('length of change: %.4f of the radial range, which is %.1f points here' % (length, length*(N - 1)))
held = int(np.sum(amp[:10] < 1e-12))
s = [np.sum(abs(np.diff(X[:, held:], j, axis=1)*(N - 1)**j)**2) for j in range(4)]
print('the same without the held points, and for the derivatives (section lengths of fe_model.py): '
      'sqrt(|X|^2/|X\'|^2) = %.4f, sqrt(|X\'|^2/|X\'\'|^2) = %.4f, sqrt(|X\'\'|^2/|X\'\'\'|^2) = %.4f'
      % tuple(np.sqrt(s[j]/s[j+1]) for j in range(3)))
print('largest |X| at %.3f of the range; largest slope %.1f per range' % (t[amp.argmax()], slope.max()))
print('inner end: %d points held at zero; slope over the next five steps %.2f per range' % (held, slope[held:held+5].mean()))
print('edge:      |X| at the last point %.1e; slope over the last five steps %.2f per range' % (amp[-1], slope[-5:].mean()))
print('|X| at 0.02, 0.05, 0.10 of the range from the inner end: %s; from the edge: %s'
      % (', '.join('%.3f' % np.interp(x, t, amp) for x in (0.02, 0.05, 0.10)),
         ', '.join('%.3f' % np.interp(1 - x, t, amp) for x in (0.02, 0.05, 0.10))))
