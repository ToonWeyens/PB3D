#!/bin/bash
# Regression test: the input checks of the perturbation modes in PB3D, and
# POST on the output of a fixed-boundary run with three solutions.
#
# Usage: post_checks.sh <PB3D-executable> <POST-executable> <input-deck> <equilibrium-file>
#   <input-deck>  a fixed-boundary deck with the lines "n_sol_requested = 1"
#                 and "n_mod_X = ...", such as input_cbm18a_fixed_bnd
#
# Checks:
#  1. PB3D refuses a deck that sets max_sec_X without min_sec_X ("min_sec_X
#     not set"), and a deck that sets max_sec_X together with n_mod_X, each
#     with an exit status that is not zero.
#  2. PB3D runs the deck with n_sol_requested = 3 and retain_all_sol, so that
#     POST has three solutions to treat.
#  3. POST writes the energy decomposition of all three solutions to
#     POST_out_EN_R1.txt, and writes no file fort.<unit>.
#  4. POST with PB3D_rich_lvl = 0 does only the equilibrium output and
#     returns exit status 0.
#  5. POST with PB3D_rich_lvl = 5 stops, with an exit status that is not
#     zero, and names 1 as the highest Richardson level in the file.
set -u

PB3D_EXE=$1
POST_EXE=$2
INPUT=$3
EQ_FILE=$4

OK=1
fail() {
    echo "FAIL: $1"
    OK=0
}

# fresh scratch dir, as PB3D and POST write several output files
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"

# the namelist of the deck, without the lines after its closing "/"
sed '/^\//q' "$INPUT" > deck_base
if ! grep -q '^ *n_sol_requested *= *1 *$' deck_base || \
    ! grep -q '^ *n_mod_X *=' deck_base; then
    echo "FAIL: the deck $INPUT has no line 'n_sol_requested = 1' or no line 'n_mod_X = ...'"
    exit 1
fi

# 1. input checks of the perturbation modes (singleton MPI, see cbm18a_fixed_bnd.sh)
check_refused() {   # name, sed script, expected message
    mkdir "$1"
    sed "$2" deck_base > "$1/deck"
    (cd "$1" && "$PB3D_EXE" deck "$EQ_FILE" > run.log 2>&1)
    local status=$?
    if [ "$status" -eq 0 ]; then
        fail "PB3D accepted the deck '$1', which it has to refuse"
    elif ! grep -qF "$3" "$1/run.log"; then
        fail "PB3D refused the deck '$1' without the message '$3'; log tail:"
        tail -10 "$1/run.log"
    fi
}
check_refused only_max_sec_X \
    '/^ *n_mod_X *=/d; s/^\/$/    max_sec_X = 15\n\//' \
    'min_sec_X not set'
check_refused max_sec_X_and_n_mod_X \
    's/^\/$/    max_sec_X = 15\n\//' \
    'Cannot prescribe min_sec_X or max_sec_X together with n_mod_X'

# 2. a run with three solutions
mkdir run
sed 's/^\( *n_sol_requested *= *\)1 *$/\13/; s/^\/$/    retain_all_sol = .true.\n\//' \
    deck_base > run/deck
(cd run && "$PB3D_EXE" deck "$EQ_FILE" > run.log 2>&1)
STATUS=$?
if [ "$STATUS" -ne 0 ] || [ ! -f run/PB3D_out.h5 ]; then
    echo "FAIL: PB3D exited with status $STATUS or wrote no PB3D_out.h5; log tail:"
    tail -20 run/run.log
    exit 1
fi
N_SOL=$(grep -c '^ *#*[0-9]' run/PB3D_out_EV_R_1.txt)
if [ "$N_SOL" -ne 3 ]; then
    echo "FAIL: PB3D kept $N_SOL solutions, not 3:"
    cat run/PB3D_out_EV_R_1.txt
    exit 1
fi

# POST input with only the energy decomposition; the extra lines are added
# before the closing "/"
write_post_input() {   # file, extra lines
    cat > "$1" <<EOF
&inputdata_POST
    plot_resonance = .false.
    plot_flux_q = .false.
    plot_kappa = .false.
    plot_magn_grid = .false.
    plot_B = .false.
    plot_J = .false.
    plot_sol_xi = .false.
    plot_sol_Q = .false.
    plot_E_rec = .true.
    plot_vac_pot = .false.
    n_theta_plot = 51
    n_zeta_plot = 1
$2
/
EOF
}
run_post() {   # name, extra input lines; sets POST_STATUS
    mkdir "$1"
    ln -s ../run/PB3D_out.h5 "$1/PB3D_out.h5"
    write_post_input "$1/input_post" "$2"
    (cd "$1" && "$POST_EXE" input_post PB3D_out.h5 > post.log 2>&1)
    POST_STATUS=$?
}

# 3. energy decomposition of all solutions in one file
run_post post_energy ""
if [ "$POST_STATUS" -ne 0 ]; then
    fail "POST exited with status $POST_STATUS; log tail:"
    tail -20 post_energy/post.log
else
    N_EN=$(grep -c '^# Eigenvalue' post_energy/POST_out_EN_R1.txt 2>/dev/null)
    if [ "$N_EN" != 3 ]; then
        fail "POST_out_EN_R1.txt holds ${N_EN:-no} solutions, not 3"
    fi
    if ls post_energy/fort.* > /dev/null 2>&1; then
        fail "POST wrote $(ls post_energy/fort.* | xargs -n1 basename | tr '\n' ' ')"
    fi
fi

# 4. equilibrium output only
run_post post_lvl_0 "    PB3D_rich_lvl = 0"
if [ "$POST_STATUS" -ne 0 ]; then
    fail "POST with PB3D_rich_lvl = 0 exited with status $POST_STATUS; log tail:"
    tail -20 post_lvl_0/post.log
elif ! grep -qF 'so only equilibrium output will be done' post_lvl_0/post.log; then
    fail "POST with PB3D_rich_lvl = 0 did not say that it does only the equilibrium output"
fi

# 5. a level above the highest one in the file
run_post post_lvl_5 "    PB3D_rich_lvl = 5"
if [ "$POST_STATUS" -eq 0 ]; then
    fail "POST accepted PB3D_rich_lvl = 5 for a file with one Richardson level"
elif ! grep -qF 'PB3D_rich_lvl cannot be higher than 1' post_lvl_5/post.log; then
    fail "POST refused PB3D_rich_lvl = 5 without naming 1 as the highest level; log tail:"
    tail -10 post_lvl_5/post.log
fi

if [ "$OK" -ne 1 ]; then
    exit 1
fi
echo "PASS: mode input checks, POST energy file with 3 solutions, PB3D_rich_lvl 0 and 5"
