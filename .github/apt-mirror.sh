#!/usr/bin/env bash
# Make apt download from archive.ubuntu.com on the GitHub-hosted runners that
# run on AWS. Run it before "apt-get update".
#
# Those runners list us-east-2.ec2.archive.ubuntu.com first in
# /etc/apt/apt-mirrors.txt. On 18 November 2026 that mirror sent 60-150 kB/s,
# so the 192 MB of the full-stack jobs took 20-50 minutes to download. The
# runners that run on Azure use azure.archive.ubuntu.com, which sends several
# MB/s, and this script leaves them alone.
set -euo pipefail

files=$(grep -lE '[a-z0-9-]+\.ec2\.archive\.ubuntu\.com' \
  /etc/apt/apt-mirrors.txt /etc/apt/sources.list \
  /etc/apt/sources.list.d/* 2>/dev/null || true)

if [ -n "$files" ]; then
  sudo sed -i -E \
    's|https?://[a-z0-9-]+\.ec2\.archive\.ubuntu\.com/|http://archive.ubuntu.com/|g' \
    $files
  echo "apt mirror: EC2 mirror replaced by archive.ubuntu.com in" $files
else
  echo "apt mirror: no EC2 mirror listed, nothing changed"
fi

if [ -f /etc/apt/apt-mirrors.txt ]; then
  cat /etc/apt/apt-mirrors.txt
fi
