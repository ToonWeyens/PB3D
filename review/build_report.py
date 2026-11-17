#!/usr/bin/env python3
"""Build review/index.html, the review guide for the changes from tag 2.47 to
dev/vacuum-fixes, from the hand-written review/report.src.html.

The source page holds the prose. Markers in it are replaced by material
generated from git, so that no code change is typed by hand or left out:

  <!--@commit SHA PATH...-->  a card for commit SHA: message, what its own diff
                              to PATH holds (by the mechanical classifier), and
                              links to the blocks of the file diffs it changed
                              (after a section's last card: the blocks of the
                              2.47 -> tip diff its commits changed, each block
                              classified and listing the commits behind it)
  <!--@wholefile PATH-->      the whole 2.47 -> tip diff of a non-vacuum module
  <!--@group KEY-->           2.47 vs tip of one group of vacuum routines
  <!--@demo NAME-->           demo script, its output and its figure
  <!--@testsrc PATH-->        the source of a test file (new since 2.47)
  <!--@ledger-->              every (commit, Fortran file) pair and where it
                              is shown, plus the coverage checks

Run from the repository root of the worktree:
    python3 review/build_report.py            (runs the demos)
    python3 review/build_report.py --no-run   (uses review/demo_output/)
It exits with status 1 if any Fortran change is not shown somewhere.
"""
import base64
import html
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REV_A = "2.47"
REV_B = "6463fc2"                     # origin/dev/vacuum-fixes on 17 November 2026: the code from stage 5 on
REV_R = "9e9aba5"                     # the code that stages 0 to 4 were reviewed at
GUIDE_R = "df519a4"                   # the commit of this guide whose page holds stages 0 to 4 as reviewed
FROZEN_A = "<!-- ================================================================== STAGE 0 -->"
FROZEN_B = "<!-- ================================================================== STAGE 5 -->"
VAC_FILES = ["Modules/vac_ops.f90", "Modules/vac_utilities.f90", "Modules/vac_vars.f90"]
GH = "https://github.com/ToonWeyens/PB3D"   # public; both revisions are pushed there
PY = os.environ.get("PB3D_PY", os.path.join(ROOT, "..", "..", "..", ".venv", "bin", "python"))
if not os.path.exists(PY):
    PY = sys.executable


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def esc(s):
    return html.escape(s, quote=False)


# ---------------------------------------------------------------- diff parsing
def parse_diff(text):
    """Unified diff -> list of files {path, hunks: [{header, old0, new0, lines:[(tag, text)]}]}."""
    files, cur, hunk = [], None, None
    for line in text.splitlines():
        if line.startswith("diff --git"):
            m = re.match(r"diff --git a/(.*) b/(.*)", line)
            cur = {"path": m.group(2) if m else line, "hunks": [], "new": False}
            files.append(cur)
            hunk = None
        elif cur is None:
            continue
        elif line.startswith("new file mode"):
            cur["new"] = True
        elif line.startswith(("index ", "--- ", "+++ ", "similarity", "rename ", "deleted file")):
            continue
        elif line.startswith("@@"):
            m = re.match(r"@@ -(\d+)(?:,\d+)? \+(\d+)(?:,\d+)? @@(.*)", line)
            hunk = {"header": line, "old0": int(m.group(1)), "new0": int(m.group(2)),
                    "ctx": m.group(3).strip(), "lines": []}
            cur["hunks"].append(hunk)
        elif hunk is not None:
            if line.startswith("\\"):
                continue
            tag = line[:1] if line[:1] in "+- " else " "
            hunk["lines"].append((tag, line[1:]))
    return files


# ---------------------------------------------------------------- classifier
INTRINSICS = [("dsqrt", "sqrt"), ("dabs", "abs"), ("dlog", "log"), ("datan", "atan"),
              ("dexp", "exp"), ("dcos", "cos"), ("dsin", "sin")]
PETSC_TOKENS = [("petscbool", "logical"), ("vecgetarrayf90", "vecgetarray"),
                ("vecrestorearrayf90", "vecrestorearray"), ("petscviewersetformat", "petscviewerpushformat")]


def norm_ws(s):
    return re.sub(r"\s+", " ", s).strip()


def is_comment(s):
    t = s.strip()
    return t.startswith("!") and not t.startswith("!$")


def code_part(line):
    """Statement part of a Fortran line: the trailing comment (outside quotes) removed."""
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


def is_use_line(s):
    t = s.strip().lower()
    return t.startswith("use ") or re.fullmatch(r"&[\w\s,=>]*(,\s*&)?|&[\w\s,=>]*&?", t or "x") is not None


def classify_hunk(hunk, moved_ctx=None):
    rem = [t for g, t in hunk["lines"] if g == "-"]
    add = [t for g, t in hunk["lines"] if g == "+"]
    changed = rem + add
    nonblank = [t for t in changed if t.strip()]
    if not nonblank:
        return "whitespace"
    if Counter(norm_ws(t) for t in rem if t.strip()) == Counter(norm_ws(t) for t in add if t.strip()):
        return "whitespace"
    if all(is_comment(t) for t in nonblank):
        return "comment"
    if Counter(c for c in map(code_part, rem) if c) == Counter(c for c in map(code_part, add) if c):
        return "comment"                     # only comments (whole-line or trailing) changed
    pp = re.compile(r"^\s*#\s*(if|ifdef|ifndef|else|endif)\b")
    r2 = Counter(norm_ws(t) for t in rem if t.strip() and not pp.match(t))
    a2 = Counter(norm_ws(t) for t in add if t.strip() and not pp.match(t))
    if r2 == a2 and any(pp.match(t) for t in changed):
        return "ldebug-guard"
    lower = lambda L, pairs: Counter(_sub_all(norm_ws(t).lower(), pairs) for t in L if t.strip())  # noqa: E731
    if lower(rem, INTRINSICS) == lower(add, INTRINSICS):
        return "intrinsic"
    if lower(rem, PETSC_TOKENS) == lower(add, PETSC_TOKENS):
        return "petsc-api"
    code = [t for t in nonblank if not is_comment(t) and not pp.match(t)]
    if code and all(is_use_line(t) for t in code):
        return "use-list"
    if moved_ctx is not None and len(code) >= 5:
        crem, cadd = moved_ctx
        hit = sum(1 for g, t in hunk["lines"] if t.strip() and not is_comment(t) and not pp.match(t)
                  and ((g == "-" and norm_ws(t) in cadd) or (g == "+" and norm_ws(t) in crem)))
        if hit >= 0.9 * len(code):
            return "moved"
    # declarations moved under #if ldebug together with use lines
    if code and all(is_use_line(t) or _is_decl(t) for t in code) and any(pp.match(t) for t in changed) \
            and set(norm_ws(t) for t in rem if _is_decl(t)) <= set(norm_ws(t) for t in add if _is_decl(t)) | set():
        return "ldebug-guard"
    return "read"


def _is_decl(t):
    return re.match(r"\s*(integer|real|logical|character|complex|type\s*\(|petsc\w+)\b.*::", t, re.I) is not None


def _sub_all(s, pairs):
    for a, b in pairs:
        s = re.sub(r"\b" + a + r"\b", b, s)
    return s


BADGE = {   # hover texts; no double quotes (they go into a title attribute)
    "read": ("read", "read: none of the mechanical patterns fits, so the block changes statements. "
                     "A person has to look; it does not mean the change is important."),
    "use-list": ("auto", "use-list: every changed code line is a Fortran use statement or its continuation, "
                         "the line that imports names from another module. Only which names are imported, "
                         "or from which module, changed. A wrong one stops the build."),
    "ldebug-guard": ("auto", "ldebug-guard: the same lines, unchanged, now sit inside #if ldebug ... #endif, "
                             "so only debug builds compile them (mostly declarations of counters that only debug "
                             "output uses). A variable used outside the guard would stop the build."),
    "comment": ("auto", "comment: only comments changed (whole ! lines or the comment after a statement); "
                        "with the comments removed, the statements are identical."),
    "whitespace": ("auto", "whitespace: old and new lines are identical apart from spaces, indentation "
                           "and blank lines."),
    "intrinsic": ("auto", "intrinsic: after replacing each old double-precision function name by its generic "
                          "one (DSQRT by sqrt, DABS by abs, DLOG by log, DATAN by atan), old and new lines are "
                          "identical. For real(dp) arguments it is the same function."),
    "petsc-api": ("auto", "petsc-api: identical after replacing old PETSc names by their new ones (PetscBool by "
                          "logical, VecGetArrayF90 by VecGetArray, PetscViewerSetFormat by PetscViewerPushFormat)."),
    "moved": ("auto", "moved: 90 % or more of its code lines appear unchanged on the other side of the same "
                      "diff: code moved from one file to another (checked line by line by the tool of M4)."),
}


# ---------------------------------------------------------------- rendering
def _mapped(m, k):
    """File line number (and routine) of line k of a group text, or None."""
    if m is None:
        return (k, None)
    return m[k - 1] if 0 < k <= len(m) else (None, None)


def render_hunk(hunk, cls=None, path="", maps=None, hid=None, trail="", ghf=None):
    o, n = hunk["old0"], hunk["new0"]
    omap, nmap = (maps or (None, None))
    out = []
    first = None
    for g, t in hunk["lines"]:
        c = {"+": "a", "-": "d", " ": "c"}[g]
        cm = " cm" if (g != " " and is_comment(t)) else ""
        lo_n, _ = _mapped(omap, o)
        ln_n, rt = _mapped(nmap, n)
        if first is None and g != "-" and ln_n:
            first = (ln_n, rt)
        lo = str(lo_n) if (g != "+" and lo_n) else ""
        ln = str(ln_n) if (g != "-" and ln_n) else ""
        out.append(f'<span class="l {c}{cm}"><i>{lo}</i><i>{ln}</i>{esc(g + t) or " "}</span>')
        if g != "+":
            o += 1
        if g != "-":
            n += 1
    badge = ""
    if cls:
        kind, tip = BADGE[cls]
        badge = f'<span class="badge {kind}" title="{esc(tip)}">{"read" if kind == "read" else esc(cls)}</span>'
    if maps is not None:
        if first is None:
            lo_n, rt = _mapped(omap, hunk["old0"])
            label = f'{path}:{lo_n} at 2.47' + (f' ({rt})' if rt else '')
        else:
            label = f'{path}:{first[0]}' + (f' ({first[1]})' if first[1] else '')
        ctx = ""
    else:
        label = f'{path}:{first_changed_line(hunk)}'
        ctx = f'<span class="hctx">{esc(hunk["ctx"])}</span>' if hunk["ctx"] else ""
    mech = " mech" if (cls and cls != "read") else ""
    idattr = f' id="{hid}"' if hid else ""
    gh = f' data-f="{esc(ghf)}" data-c="on"' if ghf else ""
    tr = f'<span class="trailc">{trail}</span>' if trail else ""
    return (f'<div class="hunk{mech}"{idattr}><div class="hh">{badge}<code>{esc(label)}</code>{ctx}{tr}</div>'
            f'<pre class="diff"{gh}>{"".join(out)}</pre></div>')


def file_block(f, title_extra="", classify=True, open_=False, moved_ctx=None, maps=None):
    adds = sum(1 for h in f["hunks"] for g, _ in h["lines"] if g == "+")
    dels = sum(1 for h in f["hunks"] for g, _ in h["lines"] if g == "-")
    classes = [classify_hunk(h, moved_ctx) if classify else None for h in f["hunks"]]
    n_read = sum(1 for c in classes if c == "read")
    n_auto = sum(1 for c in classes if c and c != "read")
    summ = f'<code>{esc(f["path"])}</code> <span class="pm">+{adds} &minus;{dels}</span>'
    if classify:
        nh = len(classes)
        summ += (f' <span class="hc">{nh} hunk{"s" if nh != 1 else ""}: '
                 + (f'<b>{n_read} to read</b>' if n_read else '0 to read') + f', {n_auto} mechanical</span>')
    summ += title_extra
    body = "".join(render_hunk(h, c, f["path"], maps, ghf=f["path"] if maps is not None else None)
                   for h, c in zip(f["hunks"], classes))
    stats["hunks"].update(c for c in classes if c)
    kind = "has-read" if (not classify or n_read) else "all-mech"
    return (f'<details class="fdiff {kind}"{" open" if open_ else ""}><summary>{summ}</summary>{body}</details>')


stats = {"hunks": Counter()}
shown_pairs = defaultdict(set)        # (sha7, path) -> set of section ids


def fortran_paths_of_commit(sha):
    out = git("show", "--format=", "--name-only", sha)
    return [p for p in out.split() if re.search(r"\.(f90|F90|h)$", p) and p.startswith(("Modules/", "Include/"))]


def moved_context(files):
    crem = set(norm_ws(t) for f in files for h in f["hunks"] for g, t in h["lines"] if g == "-" and t.strip())
    cadd = set(norm_ws(t) for f in files for h in f["hunks"] for g, t in h["lines"] if g == "+" and t.strip())
    return crem, cadd


def class_summary(classes):
    n_read = sum(1 for c in classes if c == "read")
    mech = Counter(c for c in classes if c and c != "read")
    s = f'{len(classes)} hunk{"s" if len(classes) != 1 else ""}: {n_read} to read'
    if mech:
        s += ", mechanical " + ", ".join(f"{v} {k}" for k, v in mech.most_common())
    return s


def marker_commit(args, section):
    """A commit card: subject, message, what its own diff contains, and links to
    the blocks of the 2.47 -> tip diff whose lines it changed (shown below the
    cards of the same section)."""
    sha = args[0]
    paths = args[1:]
    if paths == ["ALL"] or not paths:
        paths = fortran_paths_of_commit(sha)
    elif paths == ["NONVAC"]:
        paths = [p for p in fortran_paths_of_commit(sha) if p not in VAC_FILES]
    s7 = sha[:7]
    full = git("rev-parse", sha).strip()
    subj = git("show", "-s", "--format=%ad %s", "--date=short", sha).strip()
    date, msg = subj.split(" ", 1)
    body = git("show", "-s", "--format=%b", sha).strip()
    files = parse_diff(git("show", "--format=", "-U3", sha, "--", *paths))
    for f in files:
        shown_pairs[(s7, f["path"])].add(section)
    ctx = moved_context(files)
    classes = [classify_hunk(h, ctx) for f in files for h in f["hunks"]]
    nfiles = len(files)
    own = f'{class_summary(classes)} in {nfiles} file{"s" if nfiles != 1 else ""}'
    chips = [block_chip(section, p, i) for p in paths for i in COMMIT_BLOCKS.get((s7, p), [])]
    gone = [p for p in paths if not COMMIT_BLOCKS.get((s7, p))]
    blocks = (f'<div class="cblocks"><span class="hc">its lines, 2.47 against {REV_B}, below:</span> '
              f'{" ".join(chips)}</div>') if chips else ""
    if gone:
        blocks += (f'<p class="note">Nothing of its change to {", ".join("<code>" + esc(p.split("/")[-1]) + "</code>" for p in gone)} '
                   f'is left at {REV_B}: later commits changed those lines again. Its own diff is on '
                   f'<a href="{GH}/commit/{full}" target="_blank" rel="noopener">GitHub</a>.</p>')
    vac = [p for p in fortran_paths_of_commit(sha) if p in VAC_FILES]
    if vac:
        blocks += ('<p class="note">Its change to the vacuum modules is in the routine diffs of '
                   '<a href="#stage5">stage 5</a>.</p>')
    msgbox = (f'<details class="cmsg"><summary>commit message</summary><pre class="src">{esc(body)}</pre></details>'
              if body else "")
    card = (f'<div class="commit" id="c-{s7}"><div class="ch">'
            f'<a class="sha" href="{GH}/commit/{full}" target="_blank" rel="noopener" title="the commit on GitHub">'
            f'<code>{esc(s7)}</code></a><span class="date">{esc(date)}</span> {esc(msg)}</div>'
            f'<div class="cown">its own diff: {esc(own)}</div>{blocks}{msgbox}</div>')
    if LAST_CARD.get(section) == s7:
        card += section_code(section)
    return card


# ---------------------------------------------------------------- non-vacuum modules: 2.47 -> tip blocks
NV = {}                               # path -> {"file", "classes", "trails", "labels"}
COMMIT_BLOCKS = defaultdict(list)     # (sha7, path) -> block indices (1-based)
CARD_SECTION = {}                     # sha7 -> section id of its commit card
LAST_CARD = {}                        # section id -> sha7 of its last commit card
SUBJECT = {}                          # sha7 -> "date subject"
SECTION_ORDER = {}                    # section id -> page order
SHOWN_BLOCKS = defaultdict(set)       # (path, i) -> sections showing the block


# ---------------------------------------------------------------- stages 0 to 4, as reviewed
# Toon reviewed stages 0 to 4 on the page of commit GUIDE_R, which shows the code at REV_R.
# That part of the page is inserted unchanged; the sections from stage 5 on show REV_B.
FROZEN = {"html": "", "sections": [], "cards": {}, "files": defaultdict(set), "old": defaultdict(set)}
AFTER = {}                            # (sha7, path) -> reviewed sections, for later changes shown only as a list


def load_frozen():
    page = git("show", f"{GUIDE_R}:review/index.html")
    part = page[page.index(FROZEN_A):page.index(FROZEN_B)]
    FROZEN["html"] = part
    sec = None
    for m in re.finditer(r'<section class="area" id="([^"]+)"|<div class="commit" id="c-([0-9a-f]{7})"'
                         r'|data-f="([^"]+)"', part):
        if m.group(1):
            sec = m.group(1)
            FROZEN["sections"].append(sec)
        elif m.group(2) and sec:
            FROZEN["cards"][m.group(2)] = sec
        elif m.group(3) and sec:
            FROZEN["files"][sec].add(m.group(3))
    # where the ledger of that page showed each (commit, file)
    led = page[page.index('<table class="ledger">'):]
    led = led[:led.index("</table>")]
    for row in re.finditer(r"<tr><td><code>([0-9a-f]+)</code></td><td>.*?</td><td><ul>(.*?)</ul></td></tr>", led, re.S):
        for li in re.finditer(r'<li class="(?:ok|miss)"><code>([^<]+)</code> &rarr; (.*?)</li>', row.group(2)):
            FROZEN["old"][(row.group(1)[:7], li.group(1))].update(re.findall(r'href="#([^"]+)"', li.group(2)))


def marker_frozen():
    sha = git("rev-parse", REV_R + "^{commit}").strip()
    return f'<div class="frozen" data-rev="{sha}" data-revname="{REV_R}">\n{FROZEN["html"]}</div>\n\n'


def after_review():
    """The commits after REV_R whose change to a non-vacuum module no section from stage 5
    on shows: their lines lie in code that stages 0 to 4 show as reviewed."""
    items = []
    for h in git("log", "--reverse", "--no-merges", "--format=%h", f"{REV_R}..{REV_B}").split():
        h7 = h[:7]
        paths = [p for p in fortran_paths_of_commit(h) if p not in VAC_FILES and not shown_pairs.get((h7, p))]
        if not paths:
            continue
        full = git("rev-parse", h).strip()
        date, msg = git("show", "-s", "--format=%ad %s", "--date=short", h).strip().split(" ", 1)
        renamed = {l.split("\t")[2]: l.split("\t")[1] for l in git("show", "-M", "--format=", "--name-status", h).splitlines()
                   if l.startswith("R")}
        lis = []
        for p in paths:
            names = {p, renamed.get(p, p)}
            secs = sorted((s for s, fs in FROZEN["files"].items() if fs & names), key=lambda s: SECTION_ORDER.get(s, 99))
            AFTER[(h7, p)] = secs
            where = (", ".join(f'<a href="#{s}">{s}</a>' for s in secs) if secs
                     else "<b>no reviewed section shows this file</b>")
            was = (f' (was <code>{esc(renamed[p].split("/")[-1])}</code>)' if p in renamed else "")
            lis.append(f'<li><code>{esc(p.split("/")[-1])}</code>{was}: in {where}</li>')
        items.append(f'<li><a class="sha" href="{GH}/commit/{full}" target="_blank" rel="noopener"><code>{h7}</code></a> '
                     f'<span class="date">{esc(date)}</span> {esc(msg)}<ul>{"".join(lis)}</ul></li>')
    return f'<ul class="after">{"".join(items)}</ul>' if items else '<p class="note">None.</p>'


def fslug(path):
    return re.sub(r"\W", "_", path.split("/")[-1].rsplit(".", 1)[0])


def block_id(section, path, i):
    return f"b-{section}-{fslug(path)}-{i}"


def block_sections(path, i):
    return sorted({CARD_SECTION[c] for c in NV[path]["trails"][i - 1] if c in CARD_SECTION},
                  key=lambda s: SECTION_ORDER.get(s, 99))


def block_chip(section, path, i):
    cls = NV[path]["classes"][i - 1]
    kind = "read" if cls == "read" else "auto"
    line = NV[path]["labels"][i - 1]
    return (f'<a class="bchip {kind}" href="#{block_id(section, path, i)}" title="{esc(BADGE[cls][1])}">'
            f'{esc(path.split("/")[-1])}:{line}</a>')


def tip_range(h, ntip):
    """Tip-side line range that holds every change of a 2.47 -> tip hunk; a deletion
    is represented by the tip lines on either side of it."""
    n, pos = h["new0"], []
    for g, _ in h["lines"]:
        if g == "+":
            pos.append(n)
            n += 1
        elif g == "-":
            pos += [n - 1, n]
        else:
            n += 1
    s, e = max(1, min(pos)), min(ntip, max(pos))
    return (s, e) if s <= e else (min(s, ntip), min(s, ntip))


def build_nonvac():
    """The 2.47 -> tip diff of every non-vacuum module, classified hunk by hunk, with
    the commits that changed the lines of each hunk (git log -L)."""
    paths = [p for p in git("diff", "--name-only", REV_A, REV_B, "--", "Modules", "Include").split()
             if p not in VAC_FILES and re.search(r"\.(f90|F90|h)$", p)]
    order = git("log", "--reverse", "--no-merges", "--format=%h", f"{REV_A}..{REV_B}").split()
    order = [h[:7] for h in order]
    rank = {h: k for k, h in enumerate(order)}
    for h in order:
        SUBJECT[h] = git("show", "-s", "--format=%ad %s", "--date=short", h).strip()
    files = parse_diff(git("diff", "-U3", REV_A, REV_B, "--", *paths))
    ctx = moved_context(files)
    for f in files:
        p = f["path"]
        ntip = len(git("show", f"{REV_B}:{p}").splitlines())
        classes, trails, labels = [], [], []
        for i, h in enumerate(f["hunks"], 1):
            classes.append(classify_hunk(h, ctx))
            s, e = tip_range(h, ntip)
            out = git("log", "--no-merges", "--format=%h", "--no-patch", f"-L{s},{e}:{p}", f"{REV_A}..{REV_B}")
            tr = sorted({x[:7] for x in re.findall(r"^[0-9a-f]{7,}$", out, re.M)} & set(rank), key=rank.get)
            if not tr:
                raise SystemExit(f"no commit found for hunk {i} of {p} (lines {s}-{e})")
            trails.append(tr)
            labels.append(first_changed_line(h))
            for c in tr:
                COMMIT_BLOCKS[(c, p)].append(i)
        NV[p] = {"file": f, "classes": classes, "trails": trails, "labels": labels}
        stats["hunks"].update(classes)


def first_changed_line(h):
    """Tip line number of the first changed line of a hunk (for a pure deletion: the
    tip line after it)."""
    n = h["new0"]
    for g, _ in h["lines"]:
        if g in "+-":
            return n
        n += 1
    return n


def trail_chips(section, path, i):
    out = []
    for c in NV[path]["trails"][i - 1]:
        sec = CARD_SECTION.get(c) or FROZEN["cards"].get(c, "")
        here = " here" if sec == section else ""
        out.append(f'<a class="cchip{here}" href="#c-{c}" title="{esc(SUBJECT.get(c, ""))}">'
                   f'<code>{c}</code>{(" " + esc(sec)) if sec else ""}</a>')
    secs = block_sections(path, i)
    word = "also in" if section in secs else "explained in"
    extra = "".join(f' <a class="also" href="#{block_id(s, path, i)}">{word} {esc(s)}</a>' for s in secs if s != section)
    return '<span class="hc">changed by</span> ' + " ".join(out) + extra


def code_file_block(section, p, idx):
    d = NV[p]
    hunks = [d["file"]["hunks"][i - 1] for i in idx]
    classes = [d["classes"][i - 1] for i in idx]
    adds = sum(1 for h in hunks for g, t in h["lines"] if g == "+")
    dels = sum(1 for h in hunks for g, t in h["lines"] if g == "-")
    n_read = sum(1 for c in classes if c == "read")
    nh = len(idx)
    whole = len(idx) == len(d["classes"])
    summ = (f'<code>{esc(p)}</code> <span class="pm">+{adds} &minus;{dels}</span> '
            f'<span class="hc">{nh} block{"s" if nh != 1 else ""}'
            + ("" if whole else f' of {len(d["classes"])}') + ': '
            + (f'<b>{n_read} to read</b>' if n_read else '0 to read') + f', {nh - n_read} mechanical</span>')
    body = ""
    for i, h, c in zip(idx, hunks, classes):
        SHOWN_BLOCKS[(p, i)].add(section)
        body += render_hunk(h, c, p, hid=block_id(section, p, i), trail=trail_chips(section, p, i), ghf=p)
    kind = "has-read" if n_read else "all-mech"
    return f'<details class="fdiff {kind}" id="f-{section}-{fslug(p)}"><summary>{summ}</summary>{body}</details>'


def section_code(section):
    """The blocks of the 2.47 -> tip diff that the commits carded in this section
    changed, file by file."""
    per = defaultdict(list)
    for p in NV:
        for i in range(1, len(NV[p]["classes"]) + 1):
            if section in block_sections(p, i):
                per[p].append(i)
    if not per:
        return ""
    files = "".join(code_file_block(section, p, per[p]) for p in per)
    nb = sum(len(v) for v in per.values())
    return (f'<h4 class="codeh">The code: {nb} block{"s" if nb != 1 else ""} of the diff 2.47 &rarr; {REV_B} '
            f'that these commits changed</h4>{files}')


def marker_wholefile(args, section):
    """The whole 2.47 -> tip diff of one non-vacuum module, shown again where a
    section is about that file."""
    p = args[0]
    return code_file_block(section, p, list(range(1, len(NV[p]["classes"]) + 1)))


# ---------------------------------------------------------------- vacuum groups
ROUT_RE = re.compile(r"^( {0,4})(?:(?:integer|real\(dp\)|logical|pure|elemental|recursive)\s+)*"
                     r"(function|subroutine)\s+(\w+)", re.I)


def routines(text):
    """Top-level procedures (indent <= 4) with their doc comment: name -> (start, end), 1-based inclusive."""
    lines = text.splitlines()
    out, i = {}, 0
    while i < len(lines):
        m = ROUT_RE.match(lines[i])
        if m and not re.match(r"\s*end\s", lines[i], re.I):
            name = m.group(3)
            s = i
            while s > 0 and (lines[s - 1].strip().startswith("!") and not lines[s - 1].strip().startswith("!---")):
                s -= 1
            e = i
            pat = re.compile(r"^\s*end\s+" + m.group(2) + r"\s+" + name + r"\b", re.I)
            while e < len(lines) and not pat.match(lines[e]):
                e += 1
            out[name.lower()] = (s + 1, e + 1)
            i = e + 1
        else:
            i += 1
    return out


GROUPS = {
    # key: (file, [2.47 routines], [tip routines], title)
    "ops-header": ("Modules/vac_ops.f90", ["@header"], ["@header"], "vac_ops: module header and declarations"),
    "store": ("Modules/vac_ops.f90", ["store_vac"], ["store_vac"], "store_vac (with store_vac_VMEC, store_vac_HEL and helpers)"),
    "calcgh": ("Modules/vac_ops.f90", ["calc_gh"], ["calc_gh"], "calc_GH"),
    "gh1": ("Modules/vac_ops.f90", ["calc_gh_1"], ["calc_gh_1"], "calc_GH_1 (field-line style)"),
    "gh2": ("Modules/vac_ops.f90", ["calc_gh_2"], ["calc_gh_2"], "calc_GH_2 (axisymmetric style)"),
    "diag2": ("Modules/vac_ops.f90", [], ["set_pot_2", "diag_indirect", "set_diag_vac"], "indirect diagonals of both styles (new)"),
    "res": ("Modules/vac_ops.f90", ["calc_vac_res"], ["calc_vac_res"], "calc_vac_res"),
    "check": ("Modules/vac_ops.f90", [], ["write_points_vac", "check_vac_exact", "check_vac_exact_1", "check_vac_exact_2"], "check_vac diagnostics (new)"),
    "plot": ("Modules/vac_ops.f90", ["vac_pot_plot"], ["eval_pot_ext", "glob2dis_cols", "sub_nearest", "grad_bnd_1", "check_vac_ext", "vac_pot_plot"], "POST vacuum potential plot"),
    "solve": ("Modules/vac_ops.f90", ["solve_phi_bem"], ["solve_phi_bem"], "solve_Phi_BEM"),
    "print": ("Modules/vac_ops.f90", ["print_output_vac"], ["print_output_vac"], "print_output_vac"),
    "u-header": ("Modules/vac_utilities.f90", ["@header"], ["@header"], "vac_utilities: module header"),
    "int1": ("Modules/vac_utilities.f90", ["calc_gh_int_1"], ["calc_gh_int_1"], "calc_GH_int_1 (field-line kernels)"),
    "int2": ("Modules/vac_utilities.f90", ["calc_gh_int_2"], ["calc_gh_int_2"], "calc_GH_int_2 (axisymmetric kernels)"),
    "dist": ("Modules/vac_utilities.f90", ["vec_dis2loc", "mat_dis2loc", "interlaced_vac_copy"],
             ["vec_dis2loc", "mat_dis2loc", "interlaced_vac_copy", "n_col1_loc", "vec_glob2dis", "dis_matvec", "dis_diag"],
             "distributed-matrix helpers"),
    "direct": ("Modules/vac_utilities.f90", [], ["vac_steps", "cross3", "solid_angle_node", "H_sing_node", "calc_diag_direct"],
               "direct diagonals of the field-line style (new)"),
    "fit": ("Modules/vac_utilities.f90", [], ["fit_diag"], "fit_diag: the least-squares fit of the indirect diagonals (new)"),
    "vars": ("Modules/vac_vars.f90", ["@all"], ["@all"], "vac_vars (whole module)"),
}
group_cover = defaultdict(lambda: defaultdict(set))   # (rev, file) -> line numbers covered


def group_text(rev, path, names):
    text = git("show", f"{rev}:{path}")
    lines = text.splitlines()
    rts = routines(text)
    chunks, covered, lmap = [], set(), []
    for nm in names:
        if nm == "@all":
            s, e = 1, len(lines)
        elif nm == "@header":
            first = min(s for s, _ in rts.values()) if rts else len(lines) + 1
            s, e = 1, first - 1
        else:
            if nm.lower() not in rts:
                raise SystemExit(f"routine {nm} not found in {rev}:{path}")
            s, e = rts[nm.lower()]
        chunks.append("\n".join(lines[s - 1:e]))
        covered.update(range(s, e + 1))
        if lmap:
            lmap += [(None, None)]          # the blank separator line between routines
        label = {"@all": None, "@header": "module header"}.get(nm, nm)
        lmap += [(k, label) for k in range(s, e + 1)]
    group_cover[(rev, path)]["lines"].update(covered)
    group_cover[(rev, path)]["n"] = {len(lines)}
    group_cover[(rev, path)]["blank_or_end"].update(
        i + 1 for i, t in enumerate(lines) if not t.strip() or re.match(r"\s*(contains|end module)", t, re.I))
    return "\n\n".join(chunks) + "\n", lmap


def marker_group(args, section):
    key = args[0]
    path, old_names, new_names, title = GROUPS[key]
    old, omap = group_text(REV_A, path, old_names) if old_names else ("", [])
    new, nmap = group_text(REV_B, path, new_names) if new_names else ("", [])
    with tempfile.TemporaryDirectory() as d:
        pa, pb = os.path.join(d, "a"), os.path.join(d, "b")
        open(pa, "w").write(old); open(pb, "w").write(new)
        r = subprocess.run(["git", "diff", "--no-index", "--histogram", "-U4", pa, pb],
                           capture_output=True, text=True)
    files = parse_diff(r.stdout)
    for f in files:
        f["path"] = path
    # commit trail of the tip ranges
    trail = []
    if new_names and new_names[0] not in ("@header", "@all"):
        rts = routines(git("show", f"{REV_B}:{path}"))
        in_range = set(git("log", "--format=%h", f"{REV_A}..{REV_B}").split())
        seen = []
        for nm in new_names:
            s, e = rts[nm.lower()]
            out = git("log", "--format=%h", "--no-patch", f"-L{s},{e}:{path}", REV_B)
            for h in re.findall(r"^[0-9a-f]{7,}$", out, re.M):
                h = h[:7]
                if h in in_range and h not in seen:
                    seen.append(h)
        for h in reversed(seen):
            trail.append(git("show", "-s", "--format=%h %ad %s", "--date=short", h).strip())
    shown_pairs[("group", key)].add(section)
    blocks = "".join(file_block(f, f' <span class="hc">{esc(title)}</span>', classify=False, maps=(omap, nmap))
                     for f in files)
    if not files:
        blocks = f'<p class="note">{esc(title)}: identical.</p>'
    blocks = _ranges_line(path, omap, nmap) + blocks + _full_text(path, new, nmap, f"now ({REV_B})") \
        + _full_text(path, old, omap, "at 2.47")
    tr = ""
    if trail:
        items = "".join(f'<li><a href="{GH}/commit/{t.split(" ", 1)[0]}" target="_blank" rel="noopener">'
                        f'<code>{esc(t.split(" ", 1)[0])}</code></a> {esc(t.split(" ", 1)[1])}</li>' for t in trail)
        tr = f'<details class="trail"><summary>commits that touched these lines ({len(trail)})</summary><ol>{items}</ol></details>'
    return f'<div class="group">{blocks}{tr}</div>'


def _ranges(m):
    """Contiguous file line ranges of a group map, with routine names."""
    out, cur = [], None
    for k, rt in m:
        if k is None:
            continue
        if cur and rt == cur[2] and k == cur[1] + 1:
            cur[1] = k
        else:
            cur = [k, k, rt]
            out.append(cur)
    return out


def _ranges_line(path, omap, nmap):
    fn = path.split("/")[-1]
    def fmt(m):
        return ", ".join(f"{a}&ndash;{b}" + (f" ({esc(rt)})" if rt and rt != "module header" else "")
                         for a, b, rt in _ranges(m)) or "absent"
    return (f'<p class="ranges"><code>{esc(fn)}</code> lines now: {fmt(nmap)}; at 2.47: {fmt(omap)}. '
            f'Line numbers in the diff are file line numbers (left: 2.47, right: now).</p>')


def _full_text(path, text, lmap, when):
    if not text.strip():
        return ""
    side = "n" if when.startswith("now") else "o"
    rows = []
    for (k, _), t in zip(lmap, text.splitlines()):
        cm = " cm" if is_comment(t) else ""
        rows.append(f'<span class="l c{cm}"><i>{k or ""}</i>{esc(t) or " "}</span>')
    n = sum(1 for k, _ in lmap if k)
    return (f'<details class="fulltext"><summary>full text {esc(when)} ({n} lines)</summary>'
            f'<pre class="diff" data-f="{esc(path)}" data-c="{side}">{"".join(rows)}</pre></details>')


# ---------------------------------------------------------------- demos
DEMOS = {
    "shear_odd_spline": ("shear_odd_spline.py", "shear_odd_spline.png"),
    "edge_weight_1d": ("edge_weight_1d.py", "edge_weight_1d.png"),
    "redistribute_owner": ("redistribute_owner.py", "redistribute_owner.png"),
    "nyquist_modes": ("nyquist_modes.py", None),
    "implicit_save": ("implicit_save.f90", None),
    "fourier_basis_check": ("fourier_basis_check.py", None),
    "halfcell_integral": ("halfcell_integral.py", "halfcell_integral.png"),
    "near_field_lines": ("near_field_lines.py", "near_field_lines.png"),
    "torus_bem_toy": ("torus_bem_toy.py", "torus_bem_toy.png"),
    "setup_modes_sim": ("../../tests/reference/setup_modes_sim.py", None),
    "pure_move_check": ("../tools/pure_move_check.py", None),
    "statement_diff_225aa72": ("../tools/statement_diff.py 225aa72 Modules/vac_ops.f90", None),
    "statement_diff_028735c": ("../tools/statement_diff.py 028735c Modules/vac_ops.f90 Modules/vac_utilities.f90 Modules/vac_vars.f90", None),
}


def run_demo(name, run):
    script, _ = DEMOS[name]
    out_dir = os.path.join(ROOT, "review", "demo_output")
    os.makedirs(out_dir, exist_ok=True)
    cache = os.path.join(out_dir, name + ".txt")
    if run or not os.path.exists(cache):
        demo_dir = os.path.join(ROOT, "review", "demos")
        if script.split()[0].endswith(".f90"):
            with tempfile.TemporaryDirectory() as d:
                exe = os.path.join(d, "a.out")
                fc = "/opt/homebrew/bin/gfortran" if os.path.exists("/opt/homebrew/bin/gfortran") else "gfortran"
                subprocess.run([fc, script, "-o", exe], cwd=demo_dir, check=True)
                out = subprocess.run([exe], capture_output=True, text=True).stdout
        else:
            parts = script.split()
            out = subprocess.run([PY, *parts], cwd=demo_dir, capture_output=True, text=True).stdout
        open(cache, "w").write(out)
    return open(cache).read()


def marker_demo(args, section, run):
    name = args[0]
    script, fig = DEMOS[name]
    sfile, *sargs = script.split()
    src = open(os.path.join(ROOT, "review", "demos", sfile)).read()
    out = run_demo(name, run)
    rel = os.path.normpath(os.path.join("review", "demos", sfile))
    cmd = (f"gfortran {rel} -o /tmp/a && /tmp/a" if sfile.endswith(".f90")
           else " ".join(["python3", rel, *sargs]))
    img = ""
    if fig:
        b64 = base64.b64encode(open(os.path.join(ROOT, "review", "figures", fig), "rb").read()).decode()
        img = f'<figure><img alt="{esc(name)}" src="data:image/png;base64,{b64}"></figure>'
    kind = "tool" if "/tools/" in "/" + rel else "demo"
    return (f'<div class="demo"><div class="dh"><span class="badge demo">{kind}</span> <code>{esc(" ".join([rel, *sargs]))}</code>'
            f' <button class="copy" data-copy="{esc(cmd)}" type="button">copy command</button></div>'
            f'{img}<pre class="out">{esc(out.strip())}</pre>'
            f'<details><summary>script</summary><pre class="src">{esc(src)}</pre></details></div>')


def marker_testsrc(args, section):
    path = args[0]
    src = git("show", f"{REV_B}:{path}")
    shown_pairs[("test", path)].add(section)
    n = src.count("\n")
    rows = "".join(f'<span class="l c"><i>{k}</i>{esc(t) or " "}</span>' for k, t in enumerate(src.splitlines(), 1))
    return (f'<details class="fdiff"><summary><code>{esc(path)}</code> <span class="pm">{n} lines, new</span></summary>'
            f'<pre class="diff" data-f="{esc(path)}" data-c="n">{rows}</pre></details>')


# ---------------------------------------------------------------- ledger
def ledger():
    commits = git("log", "--reverse", "--no-merges", "--format=%h", f"{REV_A}..{REV_B}").split()
    later = {h[:7] for h in git("log", "--no-merges", "--format=%h", f"{REV_R}..{REV_B}").split()}
    frozen = set(FROZEN["sections"])
    rows, missing = [], []
    for h in commits:
        paths = fortran_paths_of_commit(h)
        if not paths:
            continue
        h7 = h[:7]
        subj = git("show", "-s", "--format=%ad %s", "--date=short", h).strip()
        cells = []
        for p in paths:
            if p in VAC_FILES:
                where = "vacuum routine diffs (V-sections)"
                ok = True
            else:
                secs = set(shown_pairs.get((h7, p), set()))
                if h7 not in later:          # shown in stages 0 to 4 as reviewed
                    secs |= FROZEN["old"].get((h7, p.split("/")[-1]), set()) & frozen
                links = [f'<a href="#{s}">{s}</a>' for s in sorted(secs, key=lambda s: SECTION_ORDER.get(s, 99))]
                if AFTER.get((h7, p)):
                    links.append('<a href="#after">changed after the review</a>')
                ok = bool(links)
                where = ", ".join(links) if ok else "NOT SHOWN"
                if not ok:
                    missing.append((h, p))
            cells.append(f'<li class="{"ok" if ok else "miss"}"><code>{esc(p.split("/")[-1])}</code> &rarr; {where}</li>')
        rows.append(f'<tr><td><code>{esc(h)}</code></td><td>{esc(subj)}</td><td><ul>{"".join(cells)}</ul></td></tr>')
    # vacuum coverage
    vac_msgs = []
    for path in VAC_FILES:
        for rev in (REV_A, REV_B):
            info = group_cover.get((rev, path))
            n = len(git("show", f"{rev}:{path}").splitlines())
            covered = info["lines"] if info else set()
            ignorable = info["blank_or_end"] if info else set()
            unc = [i for i in range(1, n + 1) if i not in covered and i not in ignorable]
            vac_msgs.append((path, rev, n, unc))
            if unc:
                missing.append((f"{rev}", f"{path} lines {unc[:12]}{'...' if len(unc) > 12 else ''}"))
    vac_rows = "".join(
        f'<li class="{"miss" if unc else "ok"}"><code>{esc(p)}</code> at {esc(r)}: {n} lines, '
        f'{"all in a group" if not unc else str(len(unc)) + " lines outside every group: " + esc(str(unc[:20]))}</li>'
        for p, r, n, unc in vac_msgs)
    # non-vacuum modules: every block that a commit carded from stage 5 on changed is shown there
    # (the other blocks are in stages 0 to 4 as reviewed, whose ledger checked them against REV_R)
    for p, d in NV.items():
        for i in range(1, len(d["classes"]) + 1):
            if block_sections(p, i) and not SHOWN_BLOCKS.get((p, i)):
                missing.append(("block", f"{p} block {i} in no section"))
    for (h7, p), secs in AFTER.items():
        if not secs:
            missing.append((h7, f"{p}: changed after the review, in no section"))
    nv_lines = sum(1 for (p, i) in SHOWN_BLOCKS for g, _ in NV[p]["file"]["hunks"][i - 1]["lines"] if g in "+-")
    untraced = sorted({(h, p) for (h, p) in shown_pairs if p in NV and not COMMIT_BLOCKS.get((h, p))})
    # tests
    test_files = [p for p in git("diff", "--name-only", REV_A, REV_B, "--", "tests").split()
                  if p.endswith((".f90", ".F90")) and "vendor" not in p]
    test_missing = [p for p in test_files if ("test", p) not in shown_pairs]
    for p in test_missing:
        missing.append(("tests", p))
    hc = stats["hunks"]
    hunk_line = ", ".join(f"{k}: {v}" for k, v in sorted(hc.items(), key=lambda kv: -kv[1]))
    status = (f"<p class='okline'>Every Fortran change between 2.47 and {REV_B} is shown in this page: up to "
              f"{REV_R} in stages 0 to 4 as you reviewed them or in the sections from stage 5 on, after {REV_R} "
              f"in the sections from stage 5 on or in the list of changes after the review.</p>" if not missing
              else "<p class='missline'>Not shown: " + esc(str(missing)) + "</p>")
    nblocks = sum(len(d["classes"]) for d in NV.values())
    nshown = len({(p, i) for (p, i) in SHOWN_BLOCKS})
    unt = ("; ".join(f"<code>{h}</code> in <code>{esc(p.split('/')[-1])}</code>" for h, p in untraced)
           if untraced else "none")
    twice = sorted((p.split("/")[-1], i, sorted(v)) for (p, i), v in SHOWN_BLOCKS.items()
                   if len(v - {"I"}) > 1)
    tw = "; ".join(f"<code>{esc(p)}</code> block {i} in {', '.join(v)}" for p, i, v in twice) or "none"
    return (f'{status}<h4>Non-vacuum modules from stage 5 on: {nshown} of the {nblocks} blocks of the diff '
            f'2.47 &rarr; {REV_B} ({len(NV)} files), {nv_lines} changed lines</h4>'
            f'<p>The other blocks lie in code that stages 0 to 4 show at {REV_R}, as reviewed; the ledger of that '
            f'page (commit {GUIDE_R}) checked that they showed every change up to {REV_R}. '
            f'Each block is shown in the section of every commit that changed it; blocks in two sections: {tw}. '
            f'Blocks by classifier verdict: {esc(hunk_line)}. Every block names at least one commit '
            f'(<code>git log -L</code> over its lines). Commits whose change to a file is not in any block, '
            f'because later commits changed the same lines again: {unt}.</p>'
            f'<h4>Vacuum modules: every line of both versions lies in one routine group</h4><ul class="cov">{vac_rows}</ul>'
            f'<h4>New test sources shown: {len(test_files) - len(test_missing)} of {len(test_files)} '
            f'(the vendored test-drive framework is left out)</h4>'
            f'<h4>Every commit that touched a Fortran module, and where its change is shown</h4>'
            f'<div class="tablewrap"><table class="ledger"><thead><tr><th>commit</th><th>date and subject</th>'
            f'<th>files &rarr; section</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>'), missing


# ---------------------------------------------------------------- fingerprint
# (deck, equilibrium, boundary, what the deck exercises) -- review/fingerprint/decks
FP_DECKS = [
    ("hel_fixed", "cbm18a", "fixed", "HELENA path, no vacuum"),
    ("hel_free4", "cbm18a", "free, BC 4", "axisymmetric vacuum (style 2), explicit edge condition"),
    ("hel_free2", "cbm18a", "free, BC 2", "axisymmetric vacuum, variational edge condition"),
    ("c005_fixed", "circular &epsilon; = 0.05", "fixed", "HELENA path, no vacuum"),
    ("c005_free2", "circular &epsilon; = 0.05", "free, BC 2", "axisymmetric vacuum"),
    ("twin_fixed", "VMEC twin, &epsilon; = 0.33", "fixed", "VMEC path, no vacuum"),
    ("twin_free2", "VMEC twin, &epsilon; = 0.33", "free, BC 2", "field-line vacuum (style 1), 64 lines"),
    ("w7x_fixed", "W7-X", "fixed", "VMEC path with ntor &gt; 0"),
    ("w7x_free2", "W7-X", "free, BC 2", "field-line vacuum on a 3-D surface, 16 lines &times; 61 points"),
]
FP_ROUND = 1e-10        # relative change below which a difference counts as rounding


def fingerprint():
    base = os.path.join(ROOT, "review", "fingerprint")
    rows = [l.split() for l in open(os.path.join(base, "commits.txt")) if l.strip()]
    res, failed = {}, {}
    rpath = os.path.join(base, "results.txt")
    for line in (open(rpath) if os.path.exists(rpath) else []):
        line = re.sub(r"\x1b?\[[0-9;]*m", "", line)          # colour codes of PB3D's log
        line = re.sub(r"\d\d:\d\d:\d\d: ", "", line)             # its time stamps
        p = re.sub(r" +", " ", line.rstrip("\n")).split(" ", 3)
        if len(p) < 3 or p[0] != "FP":
            continue
        if p[2] == "BUILD_FAILED":
            failed[p[1]] = p[3] if len(p) > 3 else ""
            continue
        lab, deck, rest = p[1], p[2], p[3] if len(p) > 3 else ""
        m = re.search(r"rc=(-?\d+) t=(\d+)s EV=(\S+) IM=(\S+) PREC=(\S+)", rest)
        if m:
            res[(lab, deck)] = dict(ev=float(m.group(3)), evs=m.group(3), im=m.group(4), prec=m.group(5),
                                    t=int(m.group(2)), rc=int(m.group(1)))
        else:
            m = re.search(r"rc=(-?\d+) t=(\d+)s NO_EV ?(.*)", rest)
            why = (m.group(3) if m else rest).split(" >> ")[0].strip()      # the error, without the call chain
            res[(lab, deck)] = dict(ev=None, why=why, rc=int(m.group(1)) if m else None)
    if not res and not failed:
        return '<p class="gap">The fingerprint has not been run yet.</p>'

    def sections_of(sha):
        secs = set()
        for (h, p), s in shown_pairs.items():
            if h == sha[:7]:
                secs |= s
        for (h, b), s in FROZEN["old"].items():
            if h == sha[:7]:
                secs |= s & set(FROZEN["sections"])
        secs.discard("group")
        if any(p in VAC_FILES for p in fortran_paths_of_commit(sha)):
            secs.add("V")
        return sorted(secs, key=lambda s: SECTION_ORDER.get(s, SECTION_ORDER.get("V0", 0) + 0.5))

    head = "".join(f'<th title="{d[3]}">{d[1]}<br><span class="fpb">{d[2]}</span></th>' for d in FP_DECKS)
    prev = {d[0]: None for d in FP_DECKS}          # last result per deck (from the previous built row)
    body, moved = [], []
    for lab, sha, sdp in rows:
        subj = git("show", "-s", "--format=%ad %s", "--date=short", sha).strip()
        secs = sections_of(sha) if lab != "2.47" else []
        sec_html = " ".join(f'<a href="#{s}">{s}</a>' if s != "V" else
                            '<a href="#V0" title="vacuum routine diffs, sections V0 to V11">V</a>' for s in secs)
        tag = ' <span class="fpb">STRUMPACK</span>' if sdp == "sdp" else ""
        first = (f'<td class="fpc"><code>{esc(lab)}</code>{tag}<div class="fps">{esc(subj)}</div>'
                 f'<div class="fps">{sec_html}</div></td>')
        if lab in failed:
            why = failed[lab]
            short = ("written for the PETSc 3.9 Fortran interface, does not compile against PETSc 3.25 "
                     "(there InsertMode is a derived type, not an integer)") \
                if re.search(r"einsertmode|petsc", why, re.I) else why[:160]
            body.append(f'<tr class="fpfail">{first}<td colspan="{len(FP_DECKS)}">build failed: '
                        f'<span title="{esc(why)}">{esc(short)}</span></td></tr>')
            continue
        cells, row_moves = [], []
        for d in FP_DECKS:
            r = res.get((lab, d[0]))
            pr = prev[d[0]]
            if r is None:
                cells.append('<td class="fpn">not run</td>')
                continue
            if r["ev"] is None:
                why = r["why"]
                cls = "fpn" if (pr is None or pr["ev"] is None) else "fpn fplost"
                if pr is not None and pr["ev"] is not None:
                    row_moves.append((d, pr, r, None))
                cells.append(f'<td class="{cls}" title="{esc(why)}">no eigenvalue</td>')
            elif pr is None or pr["ev"] is None:
                cells.append(f'<td class="fpv fpnew" title="{esc(r["evs"])}; SLEPc relative precision {esc(r["prec"])}">'
                             f'{r["ev"]:+.6e}</td>')
                if pr is not None:
                    row_moves.append((d, pr, r, None))
            else:
                rel = abs(r["ev"] - pr["ev"]) / abs(pr["ev"])
                if r["evs"] == pr["evs"]:
                    cells.append(f'<td class="fpsame" title="{esc(r["evs"])}, bit for bit as above">=</td>')
                elif rel < FP_ROUND:
                    cells.append(f'<td class="fpsame" title="{esc(r["evs"])}; relative change {rel:.1e}">'
                                 f'&asymp;</td>')
                else:
                    cells.append(f'<td class="fpv fpmoved" title="{esc(r["evs"])}; was {esc(pr["evs"])}; '
                                 f'SLEPc relative precision {esc(r["prec"])}">'
                                 f'{r["ev"]:+.6e}<div class="fps">{rel:.1e}</div></td>')
                    row_moves.append((d, pr, r, rel))
            prev[d[0]] = r
        if row_moves:
            moved.append((lab, subj, secs, row_moves))
        body.append(f'<tr>{first}{"".join(cells)}</tr>')

    def mv(d, pr, r, rel):
        name = f'{d[1]} {d[2]}'
        if r["ev"] is None:
            return f'<li>{name}: {pr["ev"]:+.6e} &rarr; no eigenvalue <span class="fps">({esc(r["why"][:140])})</span></li>'
        if pr["ev"] is None:
            return f'<li>{name}: no eigenvalue &rarr; {r["ev"]:+.6e}</li>'
        return f'<li>{name}: {pr["ev"]:+.6e} &rarr; {r["ev"]:+.6e} (relative change {rel:.1e})</li>'
    summary = "".join(
        f'<li><code>{esc(lab)}</code> {esc(subj)} '
        f'{" ".join(f"<a href=#{s}>{s}</a>" if s != "V" else "<a href=#V0>V</a>" for s in secs)}'
        f'<ul>{"".join(mv(*x) for x in mvs)}</ul></li>' for lab, subj, secs, mvs in moved)
    n_built = sum(1 for lab, _, _ in rows if lab not in failed and any((lab, d[0]) in res for d in FP_DECKS))
    decks = "".join(f'<tr><td><code>{d[0]}</code></td><td>{d[1]}</td><td>{d[2]}</td><td>{d[3]}</td></tr>'
                    for d in FP_DECKS)
    return (f'<h4>The nine decks</h4><div class="tablewrap"><table class="fpdecks"><thead><tr><th>deck</th>'
            f'<th>equilibrium</th><th>boundary</th><th>what it runs through</th></tr></thead><tbody>{decks}'
            f'</tbody></table></div>'
            f'<h4>Commits that moved an eigenvalue by more than {FP_ROUND:.0e} (relative), or made one appear or '
            f'disappear</h4><ul class="fpmoves">{summary}</ul>'
            f'<h4>Eigenvalue per commit and deck ({n_built} builds)</h4>'
            f'<p class="note"><b>=</b> the same digits as the row above; <b>&asymp;</b> a relative change below '
            f'{FP_ROUND:.0e}; a number is shown where the eigenvalue first appears or moved, with its relative '
            f'change under it (a relative change above 1 means the sign changed). Hover a cell for the full value or the error message.</p>'
            f'<div class="tablewrap"><table class="fp"><thead><tr><th>commit</th>{head}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table></div>')


# ---------------------------------------------------------------- main
def main():
    run = "--no-run" not in sys.argv
    src = open(os.path.join(ROOT, "review", "report.src.html")).read()
    section = "top"
    out = []
    pos = 0
    pat = re.compile(r'<section[^>]*\bid="([^"]+)"|<!--@(\w+)\s*([^>]*?)-->')
    sec = "top"
    load_frozen()
    for m in pat.finditer(src):             # where each commit card is, before any rendering
        if m.group(1):
            sec = m.group(1)
            SECTION_ORDER.setdefault(sec, len(SECTION_ORDER))
        elif m.group(2) == "frozen":
            for s in FROZEN["sections"]:
                SECTION_ORDER.setdefault(s, len(SECTION_ORDER))
        elif m.group(2) == "commit":
            CARD_SECTION[m.group(3).split()[0][:7]] = sec
            LAST_CARD[sec] = m.group(3).split()[0][:7]
    build_nonvac()
    ledger_at = fp_at = after_at = None
    for m in pat.finditer(src):
        out.append(src[pos:m.start()])
        pos = m.end()
        if m.group(1):
            section = m.group(1)
            out.append(m.group(0))
            continue
        kind, args = m.group(2), m.group(3).split()
        if kind == "commit":
            out.append(marker_commit(args, section))
        elif kind == "wholefile":
            out.append(marker_wholefile(args, section))
        elif kind == "group":
            out.append(marker_group(args, section))
        elif kind == "demo":
            out.append(marker_demo(args, section, run))
        elif kind == "figure":
            b64 = base64.b64encode(open(os.path.join(ROOT, "review", "figures", args[0]), "rb").read()).decode()
            out.append(f'<figure><img alt="{esc(args[0])}" src="data:image/png;base64,{b64}"></figure>')
        elif kind == "testsrc":
            out.append(marker_testsrc(args, section))
        elif kind == "ledger":
            ledger_at = len(out)
            out.append("")
        elif kind == "fingerprint":
            fp_at = len(out)
            out.append("")
        elif kind == "frozen":
            out.append(marker_frozen())
        elif kind == "afterreview":
            after_at = len(out)
            out.append("")
        else:
            raise SystemExit(f"unknown marker {kind}")
    out.append(src[pos:])
    aft = after_review()
    if after_at is not None:
        out[after_at] = aft
    led, missing = ledger()
    if ledger_at is not None:
        out[ledger_at] = led
    if fp_at is not None:
        out[fp_at] = fingerprint()
    page = "".join(out)
    page = (page.replace("@@SHA_A@@", git("rev-parse", REV_A + "^{commit}").strip())
                .replace("@@SHA_B@@", git("rev-parse", REV_B + "^{commit}").strip())
                .replace("@@NAME_B@@", REV_B).replace("@@NAME_R@@", REV_R))
    left = sorted(set(re.findall(r"@@[A-Z_0-9]+@@", page)))
    if left:
        raise SystemExit(f"placeholders left in the page: {left}")
    path = os.path.join(ROOT, "review", "index.html")
    open(path, "w").write(page)
    print(f"wrote {os.path.relpath(path, ROOT)} ({len(page) / 1e6:.2f} MB)")
    print("block classes, non-vacuum file diffs:", dict(stats["hunks"]))
    if missing:
        print("NOT SHOWN:")
        for m in missing:
            print("  ", m)
        sys.exit(1)
    print("coverage: every Fortran change is shown")


if __name__ == "__main__":
    main()
