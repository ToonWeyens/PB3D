#!/bin/bash
# bisect.sbatch -- build each commit of claude/vac-one-convention (f407afa = base of the branch
# = the tree of the cv7 series up to comments) in ~/pb3d/bisect/build-<c> and run the W7-X
# check_vac deck at 16 lines x 61 points with one process, with and without the vac_diag_style
# line, to locate the commit that changes the 3-D vacuum response (issue #4, 6 November 2026).
#SBATCH --job-name=d8_bisect
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=01:00:00
#SBATCH --output=/home/%u/pb3d/logs/%x-%j.out
set -u
source "$HOME/pb3d/env.sh"
B=$PB3D_ROOT/bisect; D=$PB3D_ROOT/decks; FX=$PB3D_FIXTURE_DIR
for c in f407afa f3ebf11 e07f0fc e027949 8d64b87 49e2441 028735c; do
    echo "=== $c: $(date '+%T')"
    cmake -S $B/src-$c -B $B/build-$c -DPSPLINE_DIR=$PSPLINE_DIR -DLIBSTELL_DIR=$LIBSTELL_DIR -DPB3D_FIXTURE_DIR=$PB3D_FIXTURE_DIR -DCMAKE_BUILD_TYPE=Release > $B/cmake-$c.log 2>&1 \
        && cmake --build $B/build-$c -j 16 --target PB3D > $B/build-$c.log 2>&1 || { echo "BUILD FAILED $c: $(grep -m3 -i 'error' $B/build-$c.log | tr '\n' '|')"; continue; }
    for deck in input_cv_w7x_na16_p61 input_cv_w7x_na16_p61_d1 input_cv_w7x_na16_p61_d2; do
        r=$B/run-$c-$deck; rm -rf $r; mkdir -p $r; cd $r
        cp $D/$deck input; ln -sf $FX/wout_w7xflip.nc wout_w7xflip.nc
        $B/build-$c/PB3D input wout_w7xflip.nc > run.log 2>&1; rc=$?
        ev=$( [ -f PB3D_out_EV_R_1.txt ] && awk '$1 == 1 {print $2}' PB3D_out_EV_R_1.txt || echo NO_EV )
        echo "RESULT $c $deck rc=$rc EV=$ev"
        grep -h "CHECKVAC identity x \|CHECKVAC identity z \|CHECKVAC solve axis_ring \|CHECKVAC solve point\|CHECKVAC diagonal n_bnd\|ERROR" run.log | cut -c1-150 | sed "s/^/   /"
        rm -f PB3D_out.h5
    done
done
echo "=== done $(date '+%T')"
