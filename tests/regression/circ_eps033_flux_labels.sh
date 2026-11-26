#!/bin/bash
# Physics regression test: the poloidal and the toroidal flux as normal
# coordinate must describe the same plasma (issue #22).
#
# Usage: circ_eps033_flux_labels.sh <PB3D-executable> <POST-executable> <input-deck> <wout-file> [python3]
#
# One Fourier mode, (m,n) = (11,10), on the circular VMEC twin
# wout_circ_eps033_q104.nc at fixed boundary: with the poloidal flux
# (use_pol_flux_F, prim_X = n = 10, secondary mode number m = 11) and with the
# toroidal flux (prim_X = m = 11, secondary mode number n = 10). An
# axisymmetric equilibrium does not couple different n, so one mode is the
# largest set that both labels have in common. The two runs differ in the
# field-aligned grid (lines of constant zeta - q theta followed in theta, or
# lines of constant iota zeta - theta followed in zeta) and in the expansion
# of the geodesic displacement U, in 1/n or in 1/m (U_style 3 in the deck).
#
# The test asks for the same eigenvalue to 1e-3 (measured 1.8e-4 on 25
# November 2026), and for each term of POST's energy decomposition, divided
# by the kinetic energy, to 1e-2 (measured 3.3e-3 or less): the normal and
# geodesic kinetic energies, and the line-bending, ballooning and kink terms.
# The normal ballooning term holds the normal curvature and the geodesic one
# the geodesic curvature. In the build of 2.55 the toroidal label had the
# normal curvature multiplied by VMEC's safety factor and a geodesic
# curvature built with the formula of the poloidal label: its normal
# ballooning term was -1.03 times and its geodesic one -2.0 times that of the
# poloidal label.
#
# Then the same mode at free boundary through the 3-D vacuum (alpha_style 2,
# 32 field lines), with tol_SLEPC = 1e-12:
#   - with the poloidal and with the toroidal flux. The eigenvalues must agree
#     to FREE_RTOL. 32 lines are 3.2 per toroidal period of the mode for
#     the poloidal flux and 2.9 per poloidal period for the toroidal flux, so
#     the vacuum response is some per cent off in both runs (it is within
#     0.5 % of the axisymmetric one at 64 lines); the test is that the vacuum
#     of the toroidal flux has the right size and sign, at a quarter of the
#     time of 64 lines.
#   - with the toroidal flux and prim_X = 11 and -11, each with the two
#     toroidal mode numbers closest to resonance (n_mod_X = 2). With -11 the
#     modes are the complex conjugates of those with 11, so the eigenvalue
#     must be the same to 1e-9.
#   - with the toroidal flux on the copy of the twin with the magnetic field
#     reversed (reverse_field.py, issue #19): the same eigenvalue to 1e-9.
#     This run needs a python3 with numpy and netCDF4 or scipy as fifth
#     argument, and is skipped with a note without it.
set -u

PB3D_EXE=$1
POST_EXE=$2
INPUT=$3
WOUT_FILE=$4

PYTHON=${5:-}
REVERSE=$(cd "$(dirname "$0")" && pwd)/reverse_field.py

EV_RTOL=1e-3
RTOL=1e-2
FREE_RTOL=2e-2
SAME_RTOL=1e-9
N_ALPHA=32

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT
cd "$WORKDIR"
# PB3D recognises a VMEC file by its .nc extension and limits path lengths
ln -s "$WOUT_FILE" wout_eq.nc

# one deck per label: fixed boundary, one mode, the namelist without the lines
# after its closing "/"
make_deck() {   # make_deck <deck> <use_pol_flux_F> <prim_X> <secondary mode number>
    sed -e '/^\//q' "$INPUT" | sed \
        -e 's/BC_style                = 1 2/BC_style                = 1 1/' \
        -e 's/tol_SLEPC               = 1.0E-6/tol_SLEPC               = 1.0E-12/' \
        -e "s/use_pol_flux_F          = .true./use_pol_flux_F          = $2/" \
        -e "s/prim_X                  = 10/prim_X                  = $3/" \
        -e "s/    n_mod_X                 = 10/    min_sec_X               = $4\\
    max_sec_X               = $4/" > "$1"
    if ! grep -q "min_sec_X *= $4" "$1" || ! grep -q "use_pol_flux_F *= $2" "$1"; then
        echo "FAIL: the deck $INPUT has not the lines that this test changes"
        exit 1
    fi
}
make_deck deck_pol .true. 10 11
make_deck deck_tor .false. 11 10

cat > input_post <<EOF_POST
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
/
EOF_POST

# run PB3D and POST for one label; the eigenvalue and the energy terms
# divided by the kinetic energy go to terms_<label>, one "name value" per line
run_label() {   # run_label <label>
    mkdir "$1"
    (cd "$1" && "$PB3D_EXE" ../deck_$1 ../wout_eq.nc > run.log 2>&1)
    STATUS=$?
    if [ "$STATUS" -ne 0 ] || [ ! -f "$1/PB3D_out.h5" ]; then
        echo "FAIL: PB3D exited with status $STATUS or wrote no PB3D_out.h5 ($1); log tail:"
        tail -20 "$1/run.log"
        exit 1
    fi
    (cd "$1" && "$POST_EXE" ../input_post PB3D_out.h5 > post.log 2>&1)
    STATUS=$?
    if [ "$STATUS" -ne 0 ] || [ ! -f "$1/POST_out_EN_R1.txt" ]; then
        echo "FAIL: POST exited with status $STATUS or wrote no POST_out_EN_R1.txt ($1); log tail:"
        tail -20 "$1/post.log"
        exit 1
    fi
    awk '/^# Eigenvalue 1$/ { k = NR }
        k && NR == k+1 { ek = $4; print "eigenvalue", $1 }
        k && NR == k+3 { print "E_kin_n", $1/ek; print "E_kin_g", $2/ek }
        k && NR == k+5 { print "line_bending_n", $1/ek; print "line_bending_g", $2/ek
            print "ballooning_n", $3/ek; print "ballooning_g", $4/ek
            print "kink_n", $5/ek; print "kink_g", $6/ek }' \
        "$1/POST_out_EN_R1.txt" > terms_$1
    rm -f "$1/PB3D_out.h5"
}
run_label pol
run_label tor

echo "term, poloidal flux, toroidal flux (energies divided by the kinetic energy):"
paste terms_pol terms_tor | awk -v evrtol="$EV_RTOL" -v rtol="$RTOL" '
    { tol = ($1 == "eigenvalue") ? evrtol : rtol
      rel = ($4 - $2) / $2; if (rel < 0) rel = -rel
      printf "%-16s %14.6e %14.6e   %8.1e\n", $1, $2, $4, rel
      n++
      if ($1 != $3 || rel > tol) { print "FAIL:", $1, "differs between the flux labels by", rel; bad = 1 }
      if (rel > worst) worst = rel }
    END { if (n != 9) { print "FAIL: " n " terms read instead of 9"; exit 1 }
          if (bad) exit 1
          print "fixed boundary: largest relative difference", worst }' || exit 1

# free boundary through the 3-D vacuum: the decks with the field lines added
make_free() {   # make_free <deck> <fixed-boundary deck> [sed expression]
    {
        echo '&inputdata_PB3D'
        echo '    alpha_style             = 2'
        echo "    n_alpha                 = $N_ALPHA"
        sed -e '1d; s/BC_style                = 1 1/BC_style                = 1 2/' -e "${3:-}" "$2"
    } > "$1"
}
run_free() {   # run_free <deck> <equilibrium>: the eigenvalue into EV
    rm -f PB3D_out*
    "$PB3D_EXE" "$1" "$2" > run.log 2>&1
    STATUS=$?
    EV=$(awk '$1 == 1 {print $2}' PB3D_out_EV_R_1.txt 2>/dev/null)
    if [ "$STATUS" -ne 0 ] || [ -z "$EV" ]; then
        echo "FAIL: PB3D exited with status $STATUS or gave no accepted eigenvalue ($1 on $2); log tail:"
        tail -20 run.log
        exit 1
    fi
}
make_free free_pol deck_pol
make_free free_tor deck_tor
TWO_MODES='/max_sec_X/d; s/    min_sec_X               = 10/    n_mod_X                 = 2/'
make_free free_tor_2 deck_tor "$TWO_MODES"
make_free free_tor_2_neg deck_tor "$TWO_MODES; s/prim_X                  = 11/prim_X                  = -11/"
if ! grep -q 'n_mod_X *= 2' free_tor_2 || grep -q '_sec_X' free_tor_2 || \
    ! grep -q 'prim_X *= -11' free_tor_2_neg; then
    echo "FAIL: the decks with two toroidal mode numbers were not made"
    exit 1
fi
run_free free_pol wout_eq.nc
EV_POL=$EV
run_free free_tor wout_eq.nc
EV_TOR=$EV
run_free free_tor_2 wout_eq.nc
EV_POS=$EV
run_free free_tor_2_neg wout_eq.nc
EV_NEG=$EV
EV_REV=
if [ -z "$PYTHON" ]; then
    echo "NOTE: the copy with the field reversed skipped: no python3"
elif ! "$PYTHON" "$REVERSE" wout_eq.nc wout_rev.nc > reverse.log 2>&1; then
    echo "NOTE: the copy with the field reversed skipped: reverse_field.py failed with $PYTHON"
else
    run_free free_tor wout_rev.nc
    EV_REV=$EV
fi

echo "free boundary, $N_ALPHA field lines: poloidal flux $EV_POL, toroidal flux $EV_TOR (rtol: $FREE_RTOL)"
echo "toroidal flux: m = 11 and -11 with two n: $EV_POS and $EV_NEG, field reversed: ${EV_REV:-skipped} (rtol: $SAME_RTOL)"
awk -v evp="$EV_POL" -v evt="$EV_TOR" -v ev2="$EV_POS" -v evn="$EV_NEG" -v evr="$EV_REV" \
    -v frtol="$FREE_RTOL" -v srtol="$SAME_RTOL" 'BEGIN {
    ok = 1
    rel = (evt - evp) / evp; if (rel < 0) rel = -rel
    reln = (evn - ev2) / ev2; if (reln < 0) reln = -reln
    if (evp >= 0 || evt >= 0) { print "FAIL: free boundary not unstable:", evp, evt; ok = 0 }
    if (rel > frtol) { print "FAIL: free boundary, the flux labels differ by", rel; ok = 0 }
    if (reln > srtol) { print "FAIL: toroidal flux, m = -11 and m = 11 differ by", reln; ok = 0 }
    relr = 0
    if (evr != "") {
        relr = (evr - evt) / evt; if (relr < 0) relr = -relr
        if (relr > srtol) { print "FAIL: toroidal flux, the copy with the field reversed differs by", relr; ok = 0 }
    }
    if (ok) { print "PASS: free boundary:", rel, "(flux labels)", reln, "(m = -11)", relr, "(field reversed)"; exit 0 }
    exit 1
}'
