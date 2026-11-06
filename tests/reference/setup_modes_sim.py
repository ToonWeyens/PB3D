"""Simulate the mode bookkeeping of PB3D's `setup_modes` (Modules/X_ops.f90).

With X_style = 2 ("fast", the default) the window of n_mod_X poloidal modes
follows nq across the radius. At each radial point kd the lowest mode number of
the window is nint(min_nm_X(kd)). When it changes by delta, |delta| modes end
and |delta| new modes start. PB3D stores the modes in n_mod_X table slots:

* the nm_X table puts the mode of rank ld (1 = lowest) in slot
  modulo(ld + ld_shift - 1, n_mod_X) + 1, where ld_shift is the running sum of
  delta (X_ops.f90, "set n or m at local normal position");
* the ind_tot / mds%sec table records, for every mode, its number, its radial
  range and the slot it is stored in (column 4).

Both tables must agree. This script copies the index formulas of

* the CURRENT code (since 63ce0f4, 2018-03-27; X_ops.f90:1169-1200 on
  dev/ci-quality): the slot of a new mode is computed from the slot of the
  mode that was created last;
* the 2020 rewrite afc9d17 (branch origin/issue_vac_not_used_VMEC, never
  merged): at every sub-step it searches the active modes for the lowest and
  highest mode number, ends one of them, gives the new mode the freed slot and
  the number (extreme +/- jd);
* a CORRECTED version: as afc9d17, but the new number is extreme +/- 1,
  because the extreme is searched again at every sub-step;
* the FIXED code (issue #4, step 5): as the current code, but the slot of a
  new mode is modulo(m_new - nm_X(1,1), n_mod_X) + 1, the formula the nm_X
  table uses, so the two tables agree by construction.

It reports, for a set of min_nm_X sequences, how many (radial point, mode)
pairs end up with a wrong slot or a wrong set of mode numbers, and exits
with status 1 if the fixed code has any (so that CTest can run it).

Run: python3 setup_modes_sim.py
"""
import sys

N = 5  # n_mod_X


def check(kd, active, tot, minseq, ld_shift, errors):
    """Compare the active modes with what the nm_X table holds at kd."""
    want = list(range(minseq[kd], minseq[kd] + N))
    got = sorted(tot[i][0] for i in active)
    if got != want:
        errors.append((kd, "mode numbers", got, want))
        return
    for i in active:
        rank = tot[i][0] - minseq[kd] + 1
        true_slot = (rank + ld_shift - 1) % N + 1
        if tot[i][3] != true_slot:
            errors.append((kd, "slot", tot[i][0], tot[i][3], true_slot))


def current(minseq):
    ld_shift = 0
    tot = [[minseq[0] + ld - 1, 1, 0, ld] for ld in range(1, N + 1)]  # [m, start, end, slot]
    cur = list(range(N))   # ind_cur: indices into tot, kept sorted from low to high m
    ind_id = N             # number of entries created so far
    errors = []
    for kd in range(1, len(minseq)):
        d = minseq[kd] - minseq[kd - 1]
        ld_shift += d
        for ld in range(1, abs(d) + 1):
            last_slot = tot[ind_id - 1][3]   # ind_tot(ind_id,4): slot of the last created mode
            if d > 0:
                tot[cur[0]][2] = kd
                tot.append([tot[cur[-1]][0] + 1, kd + 1, 0, (last_slot + ld - 1) % N + 1])
                cur = cur[1:] + [len(tot) - 1]
            else:
                tot[cur[-1]][2] = kd
                tot.append([tot[cur[0]][0] - 1, kd + 1, 0, (last_slot - ld) % N + 1])
                cur = [len(tot) - 1] + cur[:-1]
        ind_id = len(tot)
        check(kd, cur, tot, minseq, ld_shift, errors)
    return errors


def fixed(minseq):
    """The current code with the slot computed from the mode number."""
    ld_shift = 0
    m0 = minseq[0]        # nm_X(1,1)
    tot = [[minseq[0] + ld - 1, 1, 0, ld] for ld in range(1, N + 1)]
    cur = list(range(N))
    errors = []
    for kd in range(1, len(minseq)):
        d = minseq[kd] - minseq[kd - 1]
        ld_shift += d
        for ld in range(1, abs(d) + 1):
            if d > 0:
                tot[cur[0]][2] = kd
                m_new = tot[cur[-1]][0] + 1
                tot.append([m_new, kd + 1, 0, (m_new - m0) % N + 1])
                cur = cur[1:] + [len(tot) - 1]
            else:
                tot[cur[-1]][2] = kd
                m_new = tot[cur[0]][0] - 1
                tot.append([m_new, kd + 1, 0, (m_new - m0) % N + 1])
                cur = [len(tot) - 1] + cur[:-1]
        check(kd, cur, tot, minseq, ld_shift, errors)
    return errors


def searched(minseq, step_is_jd):
    ld_shift = 0
    tot = [[minseq[0] + ld - 1, 1, 0, ld] for ld in range(1, N + 1)]
    active = list(range(N))
    errors = []
    for kd in range(1, len(minseq)):
        d = minseq[kd] - minseq[kd - 1]
        ld_shift += d
        for jd in range(1, abs(d) + 1):
            step = jd if step_is_jd else 1
            imin = min(active, key=lambda i: tot[i][0])
            imax = max(active, key=lambda i: tot[i][0])
            if d > 0:
                tot[imin][2] = kd
                tot.append([tot[imax][0] + step, kd + 1, 0, tot[imin][3]])
                active[active.index(imin)] = len(tot) - 1
            else:
                tot[imax][2] = kd
                tot.append([tot[imin][0] - step, kd + 1, 0, tot[imax][3]])
                active[active.index(imax)] = len(tot) - 1
        check(kd, active, tot, minseq, ld_shift, errors)
    return errors


CASES = {
    "monotonic up (q rising, m > 0)": [10, 10, 11, 11, 12, 13, 13, 14],
    "monotonic down (nq < 0, e.g. VMEC)": [14, 13, 13, 12, 11, 11, 10],
    "up then down (q has a maximum)": [10, 11, 12, 12, 11, 10],
    "down then up (reversed shear)": [12, 11, 10, 10, 11, 12, 13],
    "jitter in flat q": [10, 10, 11, 10, 11, 11, 12],
    "jumps of 2 up": [10, 12, 14, 16],
    "jumps of 2 down": [16, 14, 12, 10],
    "jumps of 3 up": [10, 13, 16],
}

if __name__ == "__main__":
    print(f"{'case':38s} {'current':>8s} {'afc9d17':>8s} {'corrected':>10s} {'fixed':>6s}")
    n_bad = 0
    for name, seq in CASES.items():
        a = len(current(seq))
        b = len(searched(seq, step_is_jd=True))
        c = len(searched(seq, step_is_jd=False))
        d = len(fixed(seq))
        n_bad += d
        print(f"{name:38s} {a:8d} {b:8d} {c:10d} {d:6d}")
    print("\n(numbers are counts of wrong (radial point, mode) entries; 0 = correct)")
    print("PASS" if n_bad == 0 else "FAIL")
    sys.exit(0 if n_bad == 0 else 1)
