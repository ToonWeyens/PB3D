#!/bin/bash
# make_vmec_twin.sh -- the VMEC twin of a HELENA fixture (issue #4, step 6e).
#
# PB3D's --export_HEL writes a VMEC INDATA file from the HELENA mapping: the
# boundary as Fourier coefficients in the geometrical angle, pressure and iota
# as cubic splines against the normalised toroidal flux, with R_0 = B_0 = 1
# because HELENA is dimensionless (PB3D normalises either file the same way, so
# the SI scale does not matter). VMEC++ (the vmecpp pip package) then solves
# the fixed-boundary equilibrium and writes the wout file; vmec_twin.py sets
# the resolution (NS up to 319, FTOL 1e-15) and prints q on axis and at the
# edge for the comparison with HELENA.
#
# Usage: make_vmec_twin.sh <PB3D-executable> <HELENA-file> <python-with-vmecpp> [MPOL]
# Output in the current directory: input.<name> and wout_<name>.nc, <name> being
# the HELENA file's basename without .12. For circ_eps033_q104 (30 October 2026):
# q 1.04112 on axis and 1.03821 at the edge from both codes, VMEC++ in 3 s.
set -euo pipefail
PB3D_EXE=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
HEL=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
PY=$3
MPOL=${4:-12}
NAME=$(basename "$HEL" .12)
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=$(pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# the export needs the full normal range; the rest of the deck does not matter
sed -e 's/^\( *min_r_sol *= *\).*/\10.0/' -e 's/^\( *n_r_sol *= *\).*/\150/' \
    -e 's/^\( *min_n_par_X *= *\).*/\150/' "$HERE/../regression/input_circ_eps005" > "$WORK/input"
ln -s "$HEL" "$WORK/$NAME"
# answers to the interactive export: no vertical shift, default plot settings,
# prescribe iota (ncurr = 0), no perturbation, MPOL poloidal modes, the pause
printf 'no\nno\n0\nno\nyes\n%d\n\n\n\n\n' "$MPOL" > "$WORK/answers"
cd "$WORK"
if ! "$PB3D_EXE" input "$NAME" --export_HEL < answers > export.log 2>&1; then
    tail -20 export.log; exit 1
fi
cp "input.$NAME" "$OUT/"
cd "$OUT"
"$PY" -I "$HERE/vmec_twin.py" "input.$NAME" "wout_$NAME.nc"
