#!/bin/bash
# make_decks.sh -- write the nine decks of the per-commit fingerprint into
# review/fingerprint/decks/, from the regression decks of the tip (9e9aba5)
# and the W7-X deck of the 6 November bisect (check_vac line removed).
# Every variable in these decks is in the input namelist of every commit
# from 2.47 to 9e9aba5, so the same deck runs unchanged on all of them.
#
# Run from the repository root:  bash review/fingerprint/make_decks.sh <w7x-deck>
set -euo pipefail
W7X=${1:?path to input_cv_w7x_na16_p61 (from ~/pb3d/decks on Snellius)}
TIP=9e9aba5
D=review/fingerprint/decks
mkdir -p $D
reg() { git show $TIP:tests/regression/$1; }
bc() { sed "s/^\( *BC_style *= *\).*/\1$1/"; }
lines() {   # prepend alpha_style 2 and n_alpha to a deck
    { echo '&inputdata_PB3D'; echo '    alpha_style             = 2'
      echo "    n_alpha                 = $1"; sed '1d'; }
}
reg input_cbm18a_fixed_bnd | bc "1 1"            > $D/hel_fixed
reg input_cbm18a_free_bnd  | bc "1 4"            > $D/hel_free4
reg input_cbm18a_free_bnd  | bc "1 2"            > $D/hel_free2
reg input_circ_eps005      | bc "1 1"            > $D/c005_fixed
reg input_circ_eps005      | bc "1 2"            > $D/c005_free2
reg input_circ_eps033      | bc "1 1"            > $D/twin_fixed
reg input_circ_eps033      | bc "1 2" | lines 64 > $D/twin_free2
grep -v 'check_vac' "$W7X" | bc "1 1"            > $D/w7x_fixed
grep -v 'check_vac' "$W7X" | bc "1 2"            > $D/w7x_free2
ls -la $D
