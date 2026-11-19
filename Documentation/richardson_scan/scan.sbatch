#!/bin/bash
# scan.sbatch <list-file> [parallel-runs]
# Runs every line "<name> <deck> <eq-file>" of the list as one PB3D run (one process each),
# several at a time, each in ~/pb3d/rich/runs/<name>. Prints "EV <name> rc time values...".
#SBATCH --job-name=rich_scan
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=01:00:00
#SBATCH --output=/home/%u/pb3d/rich/logs/%x-%j.out
set -u
source "$HOME/pb3d/env.sh"
export R=$PB3D_ROOT/rich
LIST=${1:?list file}
P=${2:-12}
one() {
    name=$1; deck=$2; eq=$3
    rdir=$R/runs/$name; rm -rf "$rdir"; mkdir -p "$rdir"; cd "$rdir" || exit 1
    cp "$deck" input
    link=eq; case "$eq" in *.nc) link=$(basename "$eq");; esac
    ln -sf "$eq" "$link"
    t0=$(date +%s)
    "$R/build/PB3D" input "$link" > run.log 2>&1; rc=$?
    dt=$(( $(date +%s) - t0 ))
    vals=""
    for f in $(ls PB3D_out_EV_R_*.txt 2>/dev/null | sort -t_ -k5 -n); do
        vals="$vals $(awk '$1 == 1 {print $2}' "$f")"
    done
    [ -n "$vals" ] || vals=" NO_EV ($(tail -2 run.log | tr '\n' '|' | cut -c1-200))"
    echo "EV $name rc=$rc t=${dt}s$vals" | tee result
    rm -f PB3D_out.h5 fort.*
}
export -f one
grep -v '^#' "$LIST" | xargs -P "$P" -L 1 bash -c 'one "$0" "$1" "$2"'
echo "=== done $(date '+%T')"
