#!/bin/bash
# fetch.sh -- copy the fingerprint results back from Snellius:
# every FP line of the job array's output into review/fingerprint/results.txt,
# the job outputs into review/fingerprint/logs/.
set -euo pipefail
CTL="-o ControlPath=$TMPDIR/pb3d-ssh-ctl"
H=tweyens@snellius.surf.nl
mkdir -p review/fingerprint/logs
scp -q $CTL "$H:pb3d/fp/logs/fp-*.out" review/fingerprint/logs/
cat review/fingerprint/logs/fp-*.out | grep '^FP ' > review/fingerprint/results.txt
wc -l review/fingerprint/results.txt
