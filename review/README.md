# Review guide: PB3D 2.47 to dev/vacuum-fixes (9e9aba5)

A guide for reviewing by hand every Fortran change between tag 2.47 and
`dev/vacuum-fixes` at 9e9aba5 (6 November 2026), built on 14 November 2026.
Open `index.html` in a browser; the page has an index, a status and a notes
box per section, and the full diffs.

- `index.html`: the guide (also published privately as a claude.ai artifact).
  Generated; do not edit by hand.
- `report.src.html`: the prose, with markers that `build_report.py` fills
  from git (per-commit diffs, per-routine vacuum diffs, demos, test sources,
  the eigenvalue table of appendix A, ledger).
- `build_report.py`: `python3 review/build_report.py` rebuilds `index.html`
  and exits with status 1 if any Fortran change between 2.47 and 9e9aba5 is
  missing from the page (`--no-run` reuses `demo_output/`).
- `demos/`: one script per area; each prints its numbers and writes one
  figure to `figures/`. Run with `~/Code/PB3D/.venv/bin/python` or after
  `uv sync`. `demo_output/` holds their captured output.
- `tools/pure_move_check.py`: checks that the two module splits moved code
  unchanged.
- `tools/statement_diff.py SHA files...`: the Fortran statements a commit
  changed, with comments and whitespace ignored (for refactoring commits).
- `tools/timeline.py`: the commit timeline figure.
- `fingerprint/`: appendix A, the first eigenvalue of nine decks at every
  commit that changes the compiled code (48 commits), built and run on
  Snellius. `make_decks.sh` writes `decks/`, `submit.sh` unpacks the commits
  there and submits `fp.sbatch` as a job array, `fetch.sh` copies the job
  outputs to `logs/` and the result lines to `results.txt`, which
  `build_report.py` turns into the table.
