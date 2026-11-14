"""Area M4: prove that the two module splits moved code without changing it.

Commit 8ac37d7 moved the PSPLINE wrapper `spline` out of num_utilities into
the new spline_utilities; commit 17b712c moved var_1D_type and
dealloc_var_1D out of HDF5_vars into the new var_1D_vars. For each split this
compares, after normalising whitespace, every non-blank line removed from the
old module with the lines of the new module, and lists what is NOT a pure
move in either direction (module boilerplate is expected there).

Run from the repository root:  python3 review/tools/pure_move_check.py
"""
import re
import os
import subprocess
from collections import Counter


os.chdir(subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                      check=True, cwd=os.path.dirname(os.path.abspath(__file__))).stdout.strip())


def git(*a):
    return subprocess.run(["git", *a], capture_output=True, text=True, check=True).stdout


def lines(diff, sign):
    return [l[1:] for l in diff.splitlines() if l.startswith(sign) and not l.startswith(sign * 3)]


def norm(L):
    return [re.sub(r"\s+", " ", x).strip() for x in L if x.strip()]


for sha, old, new in (("8ac37d7", "Modules/num_utilities.f90", "Modules/spline_utilities.f90"),
                      ("17b712c", "Modules/HDF5_vars.f90", "Modules/var_1D_vars.f90")):
    d_old = git("show", "--format=", sha, "--", old)
    d_new = git("show", "--format=", sha, "--", new)
    rem, add_old = Counter(norm(lines(d_old, "-"))), Counter(norm(lines(d_old, "+")))
    add = Counter(norm(lines(d_new, "+")))
    not_moved = rem - add                       # removed but not found in the new module
    only_new = add - rem                        # in the new module but not removed from the old one
    print(f"== {sha}: {old} -> {new}")
    print(f"   {sum(rem.values())} non-blank lines removed from the old module, "
          f"{sum(add.values())} lines in the new module, {sum((rem & add).values())} identical")
    print("   removed but not in the new module:")
    for k, v in not_moved.items():
        print(f"      {k}" + (f"   (x{v})" if v > 1 else ""))
    print("   new lines (not moved):")
    for k, v in only_new.items():
        print(f"      {k}")
    print("   lines added to the old module:")
    for k in add_old:
        print(f"      {k}")
