#!/bin/bash
# fp.sbatch -- one element of the per-commit fingerprint (review guide,
# appendix A): build PB3D at one commit and run the nine decks of
# review/fingerprint/decks at the same time, one process each (MPI singleton).
#
# Submitted by submit.sh as a job array; element i takes line i+1 of
# ~/pb3d/fp/commits.txt, "<label> <sha> <sdp|nosdp>": sdp builds against
# STRUMPACK-Dense 1.1.1 (needed to compile before 429f003), nosdp makes sure
# CMake does not find it, so the vacuum is solved with ScaLAPACK pdgesv.
# Source trees are in ~/pb3d/fp/src/<sha> (git archive, unpacked by submit.sh).
# Prints one line per deck:
#   FP <label> <deck> rc=<rc> t=<s>s EV=<real> IM=<imag> PREC=<rel. precision>
# or FP <label> <deck> rc=<rc> t=<s>s NO_EV <first ERROR line or log tail>,
# and FP <label> BUILD_FAILED <first error lines> if the build fails.
#SBATCH --job-name=fp
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=/home/%u/pb3d/fp/logs/fp-%A_%a.out
set -u
source "$HOME/pb3d/env.sh"
FP=$PB3D_ROOT/fp
FX=$PB3D_FIXTURE_DIR
read -r LABEL SHA SDP <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $FP/commits.txt)"
SRC=$FP/src/$SHA
B=$FP/build/$LABEL
echo "label=$LABEL sha=$SHA strumpack=$SDP node=$(hostname) start=$(date '+%F %T')"

if [ "$SDP" = sdp ]; then STR=$PB3D_ROOT/deps/strumpack-dense; else STR=$FP/no-strumpack; fi
rm -rf "$B"
t0=$(date +%s)
if ! cmake -S "$SRC" -B "$B" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
        -DPSPLINE_DIR="$PSPLINE_DIR" -DLIBSTELL_DIR="$LIBSTELL_DIR" -DSTRUMPACK_DIR="$STR" \
        -DPB3D_FIXTURE_DIR="$FX" > $FP/logs/cmake-$LABEL.log 2>&1; then
    echo "FP $LABEL BUILD_FAILED cmake: $(grep -m3 -i 'error' $FP/logs/cmake-$LABEL.log | tr '\n' '|' | cut -c1-300)"
    exit 0
fi
grep -i 'strumpack' $FP/logs/cmake-$LABEL.log | head -2
if ! cmake --build "$B" -j 16 --target PB3D > $FP/logs/build-$LABEL.log 2>&1; then
    echo "FP $LABEL BUILD_FAILED make: $(grep -m3 -i 'error' $FP/logs/build-$LABEL.log | tr '\n' '|' | cut -c1-300)"
    exit 0
fi
echo "build time $(( $(date +%s) - t0 )) s"
mkdir -p $FP/bin; cp "$B/PB3D" "$FP/bin/PB3D-$LABEL"; EXE=$FP/bin/PB3D-$LABEL

run() {   # run <deck> <equilibrium>
    local deck=$1 eq=$2 r=$FP/runs/$LABEL/$1 link=eq
    rm -rf "$r"; mkdir -p "$r"; cd "$r"
    cp $FP/decks/$deck input
    mkdir -p Data                                       # older commits write vacuum files into Data/
    case "$eq" in *.nc) link=$(basename "$eq");; esac    # PB3D detects VMEC by the .nc name
    ln -sf "$eq" "$link"
    local t=$(date +%s)
    GFORTRAN_UNBUFFERED_ALL=y "$EXE" input "$link" > run.log 2>&1; local rc=$?
    local dt=$(( $(date +%s) - t ))
    if [ -f PB3D_out_EV_R_1.txt ] && awk '$1 == 1 {f=1} END {exit !f}' PB3D_out_EV_R_1.txt; then
        echo "FP $LABEL $deck rc=$rc t=${dt}s $(awk '$1 == 1 {print "EV="$2, "IM="$3, "PREC="$4; exit}' PB3D_out_EV_R_1.txt)" > fp.line
    else
        local why=$(grep -m1 -A2 'ERROR' run.log | tr '\n' ' ' | tr -s ' ' | cut -c1-250)
        [ -z "$why" ] && why="log tail: $(tail -3 run.log | tr '\n' '|' | tr -s ' ' | cut -c1-250)"
        echo "FP $LABEL $deck rc=$rc t=${dt}s NO_EV $why" > fp.line
    fi
    rm -f PB3D_out.h5 fort.*
}
run hel_fixed  $FX/cbm18a_extended_vac.12 &
run hel_free4  $FX/cbm18a_extended_vac.12 &
run hel_free2  $FX/cbm18a_extended_vac.12 &
run c005_fixed $FX/circ_eps005_q104.12 &
run c005_free2 $FX/circ_eps005_q104.12 &
run twin_fixed $FX/wout_circ_eps033_q104.nc &
run twin_free2 $FX/wout_circ_eps033_q104.nc &
run w7x_fixed  $FX/wout_w7xflip.nc &
run w7x_free2  $FX/wout_w7xflip.nc &
wait
cat $FP/runs/$LABEL/*/fp.line | sed 's/\x1b\[[0-9;]*m//g'     # without the colour codes of the log
rm -rf "$B"
echo "end=$(date '+%F %T')"
