"""Does rule A of the plan work with insert_block_mat as it is?

PB3D with a sliding set of Fourier modes (X_style 2) gives every radial point
n_mod_X slots. setup_modes (X_ops.f90) decides which Fourier mode a slot holds
at a point; setup_local_block (SLEPC_utilities.f90:212-237) takes a block that
belongs to point r_id and decides, entry by entry, in which slot it goes at
the row point and at the column point, or that it is dropped.

Rule A of the plan: the set of modes is the same at all points of an element
except its first point, which has the set of the element on its inner side;
the first point of all has the set of the first element. fill_mat hands
insert_block_mat the last point of the element as r_id.

This script replays the two Fortran routines for such a table and prints, for
every case: which entries are dropped, whether any entry lands in a slot other
than its own, whether the pattern of the matrix is symmetric, and whether an
unknown is left without any entry.

    python -B replay_slots.py
"""
import numpy as np


def tables(lowest_mode, n):
    """nm[point, slot] = Fourier mode, as setup_modes fills n and m for a set that shifts by `lowest_mode`."""
    N = len(lowest_mode)
    nm = np.zeros((N, n), int)
    nm[0] = lowest_mode[0] + np.arange(n)
    shift = 0
    for kd in range(1, N):
        shift += lowest_mode[kd] - lowest_mode[kd - 1]
        for ld in range(1, n + 1):
            nm[kd, (ld + shift - 1) % n] = lowest_mode[kd] + ld - 1
    return nm


def local_block(nm, r_id, ind, n):
    """Entries (k, m) of a block of point r_id that are placed at row point r_id + ind[0] and column point
    r_id + ind[1], with the slots they get there; and the entries that are dropped. Points count from 0."""
    kept, dropped = [], []
    for m in range(1, n + 1):
        for k in range(1, n + 1):
            k_loc = k + nm[r_id, k - 1] - nm[min(r_id + ind[0], len(nm) - 1), k - 1]
            m_loc = m + nm[r_id, m - 1] - nm[min(r_id + ind[1], len(nm) - 1), m - 1]
            if 1 <= k_loc <= n and 1 <= m_loc <= n:
                kept.append(((k, m), (k_loc, m_loc)))
            else:
                dropped.append((k, m))
    return kept, dropped


def run(label, lowest_per_element, k, n=3, first='element'):
    """lowest_per_element: lowest Fourier mode of the set of every element. first = 'element': the first point
    has the set of element 1 (the plan); a number: the first point has a set of its own with that lowest mode."""
    n_el = len(lowest_per_element)
    N = n_el*k + 1
    lowest = [lowest_per_element[0] if first == 'element' else first]
    for e in range(n_el):
        lowest += [lowest_per_element[e]]*k
    nm = tables(lowest, n)
    print('\n%s; degree %d, n_mod_X = %d' % (label, k, n))
    print('  lowest mode of the set at every point:', lowest)
    pattern = np.zeros((N*n, N*n), int)
    moved, n_dropped = 0, 0
    for e in range(n_el):
        points = list(range(e*k, e*k + k + 1))
        r_id = points[-1]
        for a in points:
            for b in points:
                kept, dropped = local_block(nm, r_id, (a - r_id, b - r_id), n)
                n_dropped += len(dropped)
                for (kk, mm), (kl, ml) in kept:
                    pattern[a*n + kl - 1, b*n + ml - 1] += 1
                    moved += (kl, ml) != (kk, mm)
                if dropped and a == b:
                    modes = sorted({int(nm[r_id, kk - 1]) for kk, mm in dropped if kk == mm})
                    print('  element %d, diagonal block of point %d: the entries of mode(s) %s are dropped'
                          % (e + 1, a + 1, modes))
    filled = pattern.astype(bool)
    print('  entries dropped: %d; entries placed in a slot other than their own: %d; pattern symmetric: %s; '
          'unknowns without any entry: %d' % (n_dropped, moved, bool((filled == filled.T).all()),
                                              int((pattern.sum(axis=1) == 0).sum())))


if __name__ == '__main__':
    run('set shifts by +1 between elements 1 and 2', [10, 11], 2)
    run('set shifts by +1 between elements 1 and 2', [10, 11], 3)
    run('set shifts by +2 between elements 1 and 2', [10, 12], 2)
    run('set shifts by -1 between elements 1 and 2 (q falling)', [10, 9], 2)
    run('set shifts by +3 = n_mod_X between elements 1 and 2', [10, 13], 2)
    run('three elements, shifts +1 and +1', [10, 11, 12], 3)
    run('NOT the plan: the first point keeps a set of its own', [11, 11], 2, first=10)
