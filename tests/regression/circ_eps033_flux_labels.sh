#!/bin/bash
# Physics regression test: the poloidal and the toroidal flux as normal
# coordinate must describe the same plasma (issue #22).
#
# Usage: circ_eps033_flux_labels.sh <PB3D-executable> <POST-executable> <input-deck> <wout-file>
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
set -u

PB3D_EXE=$1
POST_EXE=$2
INPUT=$3
WOUT_FILE=$4

EV_RTOL=1e-3
RTOL=1e-2

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
          print "PASS: largest relative difference", worst }'
