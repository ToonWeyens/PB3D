"""Where do the stored integrated perturbation matrices differ between process counts?
Usage: compare_np.py <np1 dir> <npK dir> ...  (each holding PB3D_out.h5). For every dataset
under X_2_int_R_1, reshape by its lim and print the radial indices where |a-b| > 1e-9 max|a|."""
import sys, h5py, numpy as np
def load(path):
    f = h5py.File(path + '/PB3D_out.h5'); g = f['X_2_int_R_1']; out = {}
    for k in g:
        var = g[k]['var'][()]; lim = g[k]['lim'][()].astype(int)
        shape = tuple(lim[1::2] - lim[0::2] + 1)   # lim = [min1,max1,min2,max2,...]
        out[k] = var.reshape(shape[::-1]).T if var.size == np.prod(shape) else var
    return out
ref = load(sys.argv[1])
for other in sys.argv[2:]:
    o = load(other); print(f"=== {other.split('/')[-1]} against {sys.argv[1].split('/')[-1]} ===")
    for k in sorted(ref):
        a, b = ref[k], o[k]
        if a.shape != b.shape: print(f"{k}: shape {a.shape} vs {b.shape}"); continue
        mx = np.max(np.abs(a)); d = np.abs(a - b)
        if mx == 0 or d.max() <= 1e-9 * mx: print(f"{k:10s} {a.shape} identical (max rel diff {d.max()/mx if mx else 0:.1e})"); continue
        # radial axis: the last axis
        rad = tuple(range(a.ndim - 1)); prof = d.max(axis=rad) / mx
        idx = np.where(prof > 1e-9)[0]
        print(f"{k:10s} {a.shape} max rel diff {prof.max():.2e} at radial index {prof.argmax()}; {len(idx)} of {prof.size} radial points differ: {idx[:12].tolist()}{' ...' if len(idx) > 12 else ''}")
