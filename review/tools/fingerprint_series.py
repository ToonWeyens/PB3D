"""Appendix A as a figure: the first eigenvalue of each of the nine decks at every build.

Reads review/fingerprint/commits.txt and results.txt (the Snellius run of 14 November 2026)
the same way build_report.py does, and draws one panel per deck: the real part of the first
eigenvalue against the builds in commit order. A label marks each build where a value first
appears or moves by more than FP_ROUND (relative), with the relative change.

Run from the repository root:  python3 review/tools/fingerprint_series.py
"""
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "demos"))
from _common import BLUE, plt, save  # noqa: E402

BASE = os.path.join("review", "fingerprint")
DECKS = [  # the decks of build_report.FP_DECKS, in the same order
    ("hel_fixed", "cbm18a (HELENA), fixed boundary"),
    ("hel_free4", "cbm18a (HELENA), free boundary, BC style 4"),
    ("hel_free2", "cbm18a (HELENA), free boundary, BC style 2"),
    ("c005_fixed", r"circular $\varepsilon = 0.05$ (HELENA), fixed boundary"),
    ("c005_free2", r"circular $\varepsilon = 0.05$ (HELENA), free boundary, BC style 2"),
    ("twin_fixed", r"VMEC twin $\varepsilon = 0.33$, fixed boundary"),
    ("twin_free2", r"VMEC twin $\varepsilon = 0.33$, free boundary, BC style 2 (64 field lines)"),
    ("w7x_fixed", "W7-X, fixed boundary"),
    ("w7x_free2", r"W7-X, free boundary, BC style 2 (16 field lines $\times$ 61 points)"),
]
FP_ROUND = 1e-10  # as in build_report.py
INK, INK2, MUTED, BAND = "#0b0b0b", "#52514e", "#b9b8b3", "#e6e5e1"

rows = [l.split() for l in open(os.path.join(BASE, "commits.txt")) if l.strip()]
res, failed = {}, set()
for line in open(os.path.join(BASE, "results.txt")):
    line = re.sub(r"\x1b?\[[0-9;]*m", "", line)
    line = re.sub(r"\d\d:\d\d:\d\d: ", "", line)
    p = re.sub(r" +", " ", line.rstrip("\n")).split(" ", 3)
    if len(p) < 3 or p[0] != "FP":
        continue
    if p[2] == "BUILD_FAILED":
        failed.add(p[1])
        continue
    m = re.search(r"EV=(\S+)", p[3] if len(p) > 3 else "")
    res[(p[1], p[2])] = float(m.group(1)) if m else None

labels = []
for lab, sha, _ in rows:
    date = subprocess.run(["git", "show", "-s", "--format=%cd", "--date=format:%d %b %Y", sha],
                          capture_output=True, text=True, check=True).stdout.strip()
    labels.append(f"{'429f003 (SDP)' if lab == '429f003-sdp' else lab}  {date}")

fig, axes = plt.subplots(len(DECKS), 1, figsize=(15, 2.7 * len(DECKS)), sharex=True)
for ax, (deck, title) in zip(axes, DECKS):
    ys = [None if lab in failed else res.get((lab, deck)) for lab, _, _ in rows]
    seg_x, seg_y = [], []                   # a stretch of builds with an eigenvalue, drawn as steps
    for i, y in enumerate(ys + [None]):
        if y is None:
            if seg_x:
                ax.step(seg_x, seg_y, where="post", color=BLUE, lw=2)
            seg_x, seg_y = [], []
            if i < len(ys):
                ax.axvspan(i - 0.5, i + 0.5, color=BAND, lw=0, zorder=0)
        else:
            seg_x.append(i)
            seg_y.append(y)
    xv = [i for i, y in enumerate(ys) if y is not None]
    yv = [ys[i] for i in xv]
    ax.plot(xv, yv, "o", ms=4.5, color=BLUE, mec="white", mew=1.2, zorder=3)

    lo, hi = min(yv), max(yv)
    span = (hi - lo) or abs(hi)
    prev, last, max_level = None, {True: (-99, -1), False: (-99, -1)}, 0
    for i, y in enumerate(ys):
        if y is None:
            prev = None
            continue
        moved = prev is not None and abs(y - prev) / abs(prev) > FP_ROUND
        if prev is None or moved:
            txt = f"{rows[i][1]}\n{y:+.4g}"
            if moved:
                rel = abs(y - prev) / abs(prev)
                txt += f"  ({rel:.1%})" if rel >= 1e-3 else f"  ({rel:.1e})"
            up = (y > prev) if moved else (y - lo < 0.5 * span)   # on the side the value moved to
            lx, lev = last[up]
            level = lev + 1 if i - lx < 4 else 0                  # stack labels of neighbouring builds
            last[up] = (i, level)
            max_level = max(max_level, level)
            off = 8 + 21 * level
            ax.annotate(txt, (i, y), xytext=(0, off if up else -off), textcoords="offset points",
                        ha="center", va="bottom" if up else "top", fontsize=7, color=INK2, linespacing=1.1,
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6, shrinkA=0, shrinkB=3)
                        if level else None)
        prev = y
    if lo < 0 < hi:
        ax.axhline(0, color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=1)
    pad = (0.45 + 0.3 * max_level) * span
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_title(f"{title}   [{deck}]", loc="left", fontsize=10, color=INK, pad=4)
    ax.grid(axis="x", alpha=0.15)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.ticklabel_format(axis="y", useOffset=False)

axes[-1].set_xticks(range(len(rows)))
axes[-1].set_xticklabels(labels, rotation=90, fontsize=7.5, color=INK2)
axes[-1].set_xlim(-0.6, len(rows) - 0.4)
fig.supylabel(r"Re $\omega^2 / \omega_A^2$ of the first eigenvalue  (negative = unstable)", fontsize=10)
fig.tight_layout()
save(fig, "fingerprint_series.png")
for deck, _ in DECKS:
    ys = [res.get((lab, deck)) for lab, _, _ in rows if lab not in failed]
    print(f"{deck:11s} builds with an eigenvalue: {sum(y is not None for y in ys):2d} of {len(ys)}")
