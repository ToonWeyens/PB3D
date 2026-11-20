#!/bin/bash
# mem_run.sbatch -- one PB3D run under /usr/bin/time -v, to measure the peak
# resident memory and the wall time of one executable on one deck (step 9 of
# issue #4: trigonometric basis before/after).  Positional, like run_deck.sbatch:
#   sbatch mem_run.sbatch <exe> <deck> <eq-file>
#SBATCH --job-name=pb3d-mem
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=/home/%u/pb3d/logs/%x-%j.out
set -u
source "$HOME/pb3d/env.sh"
EXE=$1; DECK=$2; EQ_FILE=$3
rdir=$PB3D_ROOT/runs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}; mkdir -p "$rdir"; cd "$rdir"
cp "$DECK" input
EQ_LINK=eq; case "$EQ_FILE" in *.nc) EQ_LINK=$(basename "$EQ_FILE");; esac
ln -sf "$EQ_FILE" "$EQ_LINK"
echo "exe=$EXE deck=$DECK eq=$EQ_FILE node=$(hostname)"
/usr/bin/time -v "$EXE" input "$EQ_LINK" > run.log 2> time.log; rc=$?
ev=$(awk '$1 == 1 {print $2}' PB3D_out_EV_R_1.txt 2>/dev/null)
rss=$(awk -F: '/Maximum resident set size/ {gsub(/ /,"",$2); print $2}' time.log)
wall=$(awk -F'): ' '/Elapsed \(wall clock\)/ {print $2}' time.log)
echo "MEM exe=$(basename $(dirname $EXE)) deck=$(basename $DECK) rc=$rc wall=$wall maxrss_kB=$rss EV=$ev"
rm -f PB3D_out.h5 fort.*
