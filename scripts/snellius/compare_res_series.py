"""Convergence of the stored 3-D vacuum response: for every PB3D output given (a PB3D_out.h5 file,
or a run directory whose first np*/PB3D_out.h5 is taken), print the diagonal of res, the extreme
eigenvalues of its Hermitian part and the Frobenius-norm distance to the last (finest) one.
Usage: compare_res_series.py <PB3D_out.h5 or run dir>..."""
import sys, glob, os, h5py, numpy as np
def load(path):
    if os.path.isdir(path): path = sorted(glob.glob(os.path.join(path, 'np*', 'PB3D_out.h5')))[0]
    f = h5py.File(path); g = f['vac_R_1'] if 'vac_R_1' in f else f['vac']
    re, im = g['RE_res/var'][()], g['IM_res/var'][()]; k = int(round(np.sqrt(re.size)))
    return (re + 1j*im).reshape((k, k), order='F'), g['sec_X/var'][()].astype(int), g['misc_vac/var'][()]
runs = [(p, *load(p)) for p in sys.argv[1:]]
ref = runs[-1][1]
for p, r, sec, misc in runs:
    ev = np.linalg.eigvalsh(0.5*(r + r.conj().T)); d = np.real(np.diag(r))
    name = os.path.basename(os.path.dirname(os.path.dirname(p))) if p.endswith('.h5') else os.path.basename(p.rstrip('/'))
    print(f"{name:32s} lines {int(misc[4]):4d} par {int(misc[3]):4d} | diag {np.array2string(d, precision=2, max_line_width=200)} | eig {ev.min():+.3g} {ev.max():.3g} | dist to last {np.linalg.norm(r-ref)/np.linalg.norm(ref):.3e}")
print("m =", runs[0][2].tolist())
