"""Compare the equilibrium (and grid) datasets of two PB3D output files.
Usage: compare_eq.py <a.h5> <b.h5>  -- prints, per dataset under grid_eq, eq_1, eq_2, X_1 and grid_X,
max|a-b|/max|a| with shapes, largest differences first."""
import sys, h5py, numpy as np
a, b = h5py.File(sys.argv[1]), h5py.File(sys.argv[2])
rows = []
def walk(ga, gb, p):
    for k in ga:
        if k not in gb: rows.append((np.inf, p + '/' + k, 'missing in b', '')); continue
        oa, ob = ga[k], gb[k]
        if isinstance(oa, h5py.Dataset):
            if oa.dtype.kind in 'fi' and oa.shape == ob.shape:
                x, y = oa[()], ob[()]; mx = np.max(np.abs(x)) if x.size else 0
                d = np.max(np.abs(x - y)) if x.size else 0
                rows.append((d / mx if mx > 0 else d, p + '/' + k, str(oa.shape), f"max|a| {mx:.3e}"))
            elif oa.shape != ob.shape: rows.append((np.inf, p + '/' + k, f'shape {oa.shape} vs {ob.shape}', ''))
        else: walk(oa, ob, p + '/' + k)
for g in ('grid_eq', 'eq_1', 'eq_2', 'grid_X_1', 'X_1', 'grid_sol', 'sol_R_1'):
    if g in a and g in b: walk(a[g], b[g], g)
    elif g in a or g in b: rows.append((np.inf, g, 'only in one file', ''))
rows.sort(key=lambda r: -r[0])
for r in rows[:28]: print(f"{r[0]:10.3e}  {r[1]:45s} {r[2]:22s} {r[3]}")
print(f"... {len(rows)} datasets compared; identical ones (rel diff 0) : {sum(1 for r in rows if r[0] == 0)}")
