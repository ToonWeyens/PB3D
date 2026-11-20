"""Diagonal of the stored vacuum response of one PB3D output against another, matched by |m|.
Usage: compare_res.py <reference PB3D_out.h5> <other PB3D_out.h5> [label]"""
import sys, h5py, numpy as np
def load(path):
    f = h5py.File(path); g = f['vac_R_1'] if 'vac_R_1' in f else f['vac']
    re, im = g['RE_res/var'][()], g['IM_res/var'][()]; k = int(round(np.sqrt(re.size)))
    return (re + 1j*im).reshape((k, k), order='F'), g['sec_X/var'][()].astype(int)
r0, s0 = load(sys.argv[1]); r1, s1 = load(sys.argv[2]); d0, d1 = np.real(np.diag(r0)), np.real(np.diag(r1))
ratio = np.array([d1[i]/d0[list(np.abs(s0)).index(abs(m))] for i, m in enumerate(s1)])
ev = np.linalg.eigvalsh(0.5*(r1 + r1.conj().T))
print(f"{sys.argv[3] if len(sys.argv) > 3 else ''} m = {np.abs(s1).tolist()}")
print(f"   diag ratio: {np.array2string(ratio, precision=3, max_line_width=200)}   eig min {ev.min():+.3g} max {ev.max():.3g}")
