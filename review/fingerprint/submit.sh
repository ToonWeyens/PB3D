#!/bin/bash
# submit.sh -- prepare and submit the per-commit fingerprint on Snellius
# (review guide, appendix A). Run from the repository root on the Mac, with
# an SSH control connection open at $TMPDIR/pb3d-ssh-ctl:
#   bash review/fingerprint/submit.sh [array-spec]
# array-spec defaults to all elements, e.g. "49" submits only the tip (a test).
# It writes ~/pb3d/fp/{commits.txt,decks/,src/<sha>/} on Snellius and submits
# fp.sbatch as one job array, at most 12 elements running at a time.
# Afterwards: bash review/fingerprint/fetch.sh, then `rm -rf ~/pb3d/fp` there.
set -euo pipefail
CTL="-o ControlPath=$TMPDIR/pb3d-ssh-ctl"
H=tweyens@snellius.surf.nl
REV_A=2.47; REV_B=9e9aba5
# commits that change what is compiled into PB3D, oldest first; tag 2.47 itself first
SHAS=( $(git rev-parse --short "$REV_A^{commit}")
       $(git log --reverse --format=%h "$REV_A..$REV_B" -- Modules PB3D.f90 CMakeLists.txt cmake include) )
# before 429f003 the vacuum module needs STRUMPACK-Dense to compile; 429f003 is built both ways
FIRST_OPTIONAL=429f003
list=review/fingerprint/commits.txt; : > $list
sdp=sdp
for s in "${SHAS[@]}"; do
    if [ "$s" = "$FIRST_OPTIONAL" ]; then echo "$s-sdp $s sdp" >> $list; sdp=nosdp; fi
    echo "$s $s $sdp" >> $list
done
sed -i '' "1s/^[^ ]* /2.47 /" $list          # label the tag by its name
N=$(wc -l < $list | tr -d ' ')
echo "$N builds:"; cat $list | tr '\n' ';'; echo

ssh $CTL $H 'mkdir -p ~/pb3d/fp/src ~/pb3d/fp/logs ~/pb3d/fp/decks ~/pb3d/fp/runs ~/pb3d/fp/no-strumpack'
scp -q $CTL review/fingerprint/decks/* $H:pb3d/fp/decks/
scp -q $CTL review/fingerprint/fp.sbatch $list $H:pb3d/fp/
for s in "${SHAS[@]}"; do
    git archive --format=tar.gz "$s" | ssh $CTL $H "rm -rf ~/pb3d/fp/src/$s; mkdir -p ~/pb3d/fp/src/$s; tar -xzf - -C ~/pb3d/fp/src/$s"
done
ARR=${1:-0-$((N - 1))%12}
ssh $CTL $H "du -sh ~/pb3d/fp/src; cd ~/pb3d/fp; sbatch --array=$ARR fp.sbatch"
