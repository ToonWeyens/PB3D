"""Area 2d: VMEC's two mode lists (VMEC_ops.read_VMEC, commit 0d49f19).

A VMEC wout file stores R, Z and lambda on the main list (xm, xn; mnmax
modes, set by mpol and ntor) and the Jacobian gmn, |B| bmn and the covariant
components bsub*mn on the Nyquist list (xm_nyq, xn_nyq; mnmax_nyq modes, set
by the Nyquist limits of the real-space grid). Both lists are ordered
m = 0: n = 0..ntor; m >= 1: n = -ntor..ntor (xn = n * nfp), each with its
own limits.

At 2.47 read_VMEC took gmnc(1:mnmax,:), bmnc(1:mnmax,:) and bsub*(1:mnmax,:)
and labelled them with the MAIN list. This script builds both lists and
counts how many of the first mnmax Nyquist entries carry a different (m, n)
than the main-list label they were given.

If the fixture wout_circ_eps033_q104.nc.xz can be found, it also checks the
claim on a real file (ntor = 0).
"""
import lzma
import os
import tempfile


def vmec_list(m_max, n_max, nfp):
    """VMEC ordering: m = 0 with n = 0..n_max, then m = 1..m_max with n = -n_max..n_max."""
    out = [(0, n * nfp) for n in range(0, n_max + 1)]
    for m in range(1, m_max + 1):
        out += [(m, n * nfp) for n in range(-n_max, n_max + 1)]
    return out


cases = [
    ("circular twin (ntor 0)", 12, 0, 1, 15, 0),
    ("rippled twin (mpol 12, ntor 6, nfp 5)", 12, 6, 5, 15, 8),
    ("W7-X STELLOPT file (mpol 12, ntor 12, nfp 5)", 12, 12, 5, 11, 12),
]
print(f"{'case':46} {'mnmax':>6} {'mnmax_nyq':>10} {'mislabelled at 2.47':>20}")
for name, mpol, ntor, nfp, mnyq, nnyq in cases:
    main = vmec_list(mpol - 1, ntor, nfp)
    nyq = vmec_list(mnyq, nnyq, nfp)
    k = min(len(main), len(nyq))
    bad = sum(1 for a, b in zip(main[:k], nyq[:k]) if a != b)
    print(f"{name:46} {len(main):6d} {len(nyq):10d} {bad:20d}")
    if ntor == 6:
        print("   first entries:  main", main[7:11], " nyquist", nyq[7:11])
print("(The W7-X row assumes a Nyquist list of the same size as the main list,")
print(" as in that STELLOPT file, mnmax_nyq = mnmax = 288: then nothing is mislabelled.)")

here = os.path.dirname(os.path.abspath(__file__))
fx = os.path.join(here, "..", "..", "tests", "fixtures", "wout_circ_eps033_q104.nc.xz")
if os.path.exists(fx):
    try:
        import netCDF4
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "wout.nc")
            with open(p, "wb") as fo:
                fo.write(lzma.open(fx).read())
            ds = netCDF4.Dataset(p)
            xm, xn = ds["xm"][:], ds["xn"][:]
            xmn, xnn = ds["xm_nyq"][:], ds["xn_nyq"][:]
            same = all(int(xm[i]) == int(xmn[i]) and int(xn[i]) == int(xnn[i]) for i in range(len(xm)))
            print(f"fixture: mnmax {len(xm)}, mnmax_nyq {len(xmn)}; first mnmax Nyquist entries equal the main list: {same}")
            ds.close()
    except ImportError:
        print("netCDF4 not installed: fixture check skipped")
