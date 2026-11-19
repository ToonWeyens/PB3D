#!/bin/bash
# usage: CLUSTER=<user>@<host> ./fetch.sh   (needs an open ssh control connection at $TMPDIR/pb3d-ssh-ctl)
# fetch.sh -- collect all result lines from Snellius into results.txt
ssh -S $TMPDIR/pb3d-ssh-ctl $CLUSTER 'cat ~/pb3d/rich/runs/*/result 2>/dev/null' > results.txt
wc -l results.txt; grep -c "NO_EV" results.txt
ssh -S $TMPDIR/pb3d-ssh-ctl $CLUSTER 'squeue --me -h -o "%j %T %M"'
