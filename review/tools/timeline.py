"""Orientation figure: the 92 commits from tag 2.47 to 9e9aba5 by date and kind.

Each commit is a dot in the row of its kind (the areas of the review guide);
commits that touch no Fortran module go into "tests, CI, build" or "docs,
plans, scripts" by the files they touch. Two time windows, because nothing
happened between 10 July and 27 October 2026.

Run from the repository root:  python3 review/tools/timeline.py
"""
import datetime as dt
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "demos"))
from _common import BLUE, GREEN, GREY, ORANGE, RED, plt, save  # noqa: E402

AREA = {  # commit -> row (Fortran-touching commits; from the review guide)
    "a8c07b6": "mechanical", "775219a": "mechanical", "d4812f7": "mechanical", "1862df0": "mechanical",
    "8ac37d7": "mechanical", "6ce5029": "mechanical", "17716c4": "mechanical", "17b712c": "mechanical",
    "28f8dab": "fixes outside the vacuum", "e4e3fdd": "fixes outside the vacuum",
    "798fe50": "fixes outside the vacuum", "0d49f19": "fixes outside the vacuum",
    "4915d4b": "fixes outside the vacuum", "41b5274": "fixes outside the vacuum",
    "1d74b8f": "fixes outside the vacuum", "83c02df": "fixes outside the vacuum",
    "3b99109": "Fourier basis (VMEC)", "329b7a4": "Fourier basis (VMEC)", "1a53f50": "Fourier basis (VMEC)",
    "ad2313b": "vacuum", "429f003": "vacuum", "3679ee7": "vacuum", "d6d49e6": "vacuum", "cc4c96d": "vacuum",
    "212e3ec": "vacuum", "50c9511": "vacuum", "e3b4847": "vacuum", "6894199": "vacuum", "612fc1f": "vacuum",
    "fe305fd": "vacuum", "13bdef4": "vacuum", "83724cc": "vacuum", "0bbda47": "vacuum", "d8a1c6c": "vacuum",
    "f407afa": "vacuum", "f3ebf11": "vacuum", "e07f0fc": "vacuum", "e027949": "vacuum", "8d64b87": "vacuum",
    "49e2441": "vacuum", "028735c": "vacuum", "225aa72": "vacuum", "9af0dcc": "vacuum", "a8804ac": "vacuum",
    "3247911": "vacuum", "9e9aba5": "vacuum",
}
ROWS = ["vacuum", "Fourier basis (VMEC)", "fixes outside the vacuum", "mechanical",
        "tests, CI, build", "docs, plans, scripts"]
COL = {"vacuum": RED, "Fourier basis (VMEC)": ORANGE, "fixes outside the vacuum": BLUE,
       "mechanical": GREY, "tests, CI, build": GREEN, "docs, plans, scripts": "#8a6fb0"}

out = subprocess.run(["git", "log", "--reverse", "--format=%h|%aI|%an", "2.47..9e9aba5"],
                     capture_output=True, text=True, check=True).stdout.splitlines()
pts = []
for line in out:
    h, d, who = line.split("|", 2)
    if h in AREA:
        row = AREA[h]
    else:
        files = subprocess.run(["git", "show", "--format=", "--name-only", h],
                               capture_output=True, text=True, check=True).stdout.splitlines()
        test_like = any(f.startswith(("tests/", ".github/", "CMakeLists", "Libraries/", "spack", "fortitude"))
                        for f in files)
        row = "tests, CI, build" if test_like else "docs, plans, scripts"
    pts.append((dt.datetime.fromisoformat(d), row, who))

wins = [(dt.datetime(2026, 7, 5), dt.datetime(2026, 7, 11)), (dt.datetime(2026, 10, 26), dt.datetime(2026, 11, 7))]
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True, gridspec_kw={"width_ratios": [1, 2]})
for ax, (a, b) in zip(axes, wins):
    for t, row, who in pts:
        if a <= t.replace(tzinfo=None) <= b:
            y = len(ROWS) - 1 - ROWS.index(row)
            ax.plot(t.replace(tzinfo=None), y, "o", ms=6, color=COL[row],
                    mfc=COL[row] if who != "Claude" else "white", mew=1.5)
    ax.set_xlim(a, b)
    ax.set_yticks(range(len(ROWS)))
    ax.set_yticklabels(ROWS[::-1])
    ax.tick_params(axis="x", labelrotation=0, labelsize=8)
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter("%d %b"))
axes[0].set_title("July 2026 (open dots: author 'Claude', cloud session)")
axes[1].set_title("27 October to 6 November 2026 (issue #4)")
save(fig, "timeline.png")
n = {r: sum(1 for p in pts if p[1] == r) for r in ROWS}
print("commits per row:", n, "total", len(pts))
