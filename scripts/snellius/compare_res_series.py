"""Convergence of the stored 3-D vacuum response: for every run dir given, print the eigenvalues of
the Hermitian part of res and its diagonal; then the Frobenius-norm distance of each res to the
last (finest) one. Usage: compare_res_series.py <run dir>..."""
import sys, h5py, numpy as np
def load(path):
    f = h5py.File(path + '/np0/PB3D_out.h5'); g = f['vac_R_1'] if 'vac_R_1' in f else f['vac']
    re, im = g['RE_res/var'][()], g['IM_res/var'][()]; k = int(round(np.sqrt(re.size)))
    return (re + 1j*im).reshape((k, k), order='F'), g['sec_X/var'][()].astype(int), g['misc_vac/var'][()]
runs = [(p, *load(p)) for p in sys.argv[1:]]
ref = runs[-1][1]
for p, r, sec, misc in runs:
    ev = np.linalg.eigvalsh(0.5*(r + r.conj().T)); d = np.real(np.diag(r))
    name = p.rstrip('/').split('/')[-1]
    print(f"{name:28s} lines {int(misc[4]):4d} par {int(misc[3]):4d} | diag {np.array2string(d, precision=3, max_line_width=200)} | eig min {ev.min():+.3g} max {ev.max():.3g} | dist to finest {np.linalg.norm(r-ref)/np.linalg.norm(ref):.3e}")
print("m =", runs[0][2].tolist())
