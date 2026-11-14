"""What a commit changed in the Fortran statements, ignoring comments,
blank lines, whitespace, and the lines that only restructure code.

For a refactoring commit (code moved into contained procedures, a repeated
expression replaced by a function) the full diff is long while the change of
statements is short. This prints, for the given commit and files, the code
lines that were removed and not added back, and those added and not removed
before (multisets, comments and trailing comments stripped, whitespace
normalised). What remains is what to read.

Run from the repository root:
    python3 review/tools/statement_diff.py 225aa72 Modules/vac_ops.f90
    python3 review/tools/statement_diff.py 028735c Modules/vac_ops.f90 Modules/vac_utilities.f90 Modules/vac_vars.f90
"""
import re
import os
import subprocess
import sys
from collections import Counter


os.chdir(subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                      check=True, cwd=os.path.dirname(os.path.abspath(__file__))).stdout.strip())


def git(*a):
    return subprocess.run(["git", *a], capture_output=True, text=True, check=True).stdout


def code_part(line):
    """The statement part of a Fortran line: trailing comment (outside quotes) removed."""
    out, quote = [], None
    for ch in line.strip():
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return re.sub(r"\s+", " ", "".join(out)).strip()


def statements(diff, sign):
    lines = [l[1:] for l in diff.splitlines() if l.startswith(sign) and not l.startswith(sign * 3)]
    return Counter(c for c in (code_part(l) for l in lines) if c)


def main(sha, files):
    diff = git("show", "--format=", "-U0", sha, "--", *files)
    rem, add = statements(diff, "-"), statements(diff, "+")
    subj = git("show", "-s", "--format=%h %ad %s", "--date=short", sha).strip()
    print(f"== {subj}")
    print(f"   statements removed {sum(rem.values())}, added {sum(add.values())}, "
          f"identical on both sides {sum((rem & add).values())}")
    print("   removed and not added back:")
    for k, v in (rem - add).items():
        print(f"     - {k}" + (f"   (x{v})" if v > 1 else ""))
    print("   added and not there before:")
    for k, v in (add - rem).items():
        print(f"     + {k}" + (f"   (x{v})" if v > 1 else ""))


if __name__ == "__main__":
    args = sys.argv[1:] or ["225aa72", "Modules/vac_ops.f90"]
    main(args[0], args[1:] or ["Modules/"])
