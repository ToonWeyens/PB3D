"""Compare the matrices PB3D assembled with 1 and 2 processes (PETSc binary dumps, complex scalars).
Usage: compare_petsc_mats.py mats_np1.dat mats_np2.dat [block size]"""
import sys, numpy as np
def read_mats(path):
    raw = np.fromfile(path, dtype='>i4'); b = open(path, 'rb').read(); off = 0; mats = []
    while off < len(b):
        hdr = np.frombuffer(b, dtype='>i4', count=4, offset=off); off += 16
        cid, rows, cols, nnz = [int(x) for x in hdr]; assert cid == 1211216, (cid, off)
        rl = np.frombuffer(b, dtype='>i4', count=rows, offset=off); off += 4*rows
        ci = np.frombuffer(b, dtype='>i4', count=nnz, offset=off); off += 4*nnz
        vals = np.frombuffer(b, dtype='>f8', count=2*nnz, offset=off); off += 16*nnz
        v = vals[0::2] + 1j*vals[1::2]
        ri = np.repeat(np.arange(rows), rl); mats.append((rows, cols, ri, ci, v))
    return mats
a, b2 = read_mats(sys.argv[1]), read_mats(sys.argv[2]); bs = int(sys.argv[3]) if len(sys.argv) > 3 else 20
print(f"{len(a)} matrices in file 1, {len(b2)} in file 2")
for k, ((r, c, ri, ci, v), (r2, c2, ri2, ci2, v2)) in enumerate(zip(a, b2)):
    print(f"--- matrix {k+1}: {r}x{c}, nnz {v.size} vs {v2.size} ---")
    d1 = {}; d2 = {}
    for i, j, x in zip(ri, ci, v): d1[(int(i), int(j))] = d1.get((int(i), int(j)), 0) + x
    for i, j, x in zip(ri2, ci2, v2): d2[(int(i), int(j))] = d2.get((int(i), int(j)), 0) + x
    keys = set(d1) | set(d2); mx = max(abs(x) for x in d1.values())
    diffs = [(abs(d1.get(kk, 0) - d2.get(kk, 0)) / mx, kk) for kk in keys]
    diffs = [t for t in diffs if t[0] > 1e-12]
    diffs.sort(reverse=True)
    print(f"   max |A| {mx:.3e}; entries differing by more than 1e-12 relative: {len(diffs)} of {len(keys)}")
    rows_d = sorted({kk[0] // bs for _, kk in diffs}); print(f"   block rows with differences: {rows_d[:40]}{' ...' if len(rows_d) > 40 else ''}")
    for dd, kk in diffs[:8]: print(f"   row {kk[0]:5d} (block {kk[0]//bs:3d}) col {kk[1]:5d} (block {kk[1]//bs:3d}): {d1.get(kk,0):.6e} vs {d2.get(kk,0):.6e}  rel {dd:.2e}")
