#!/usr/bin/env bash
# Container for PRECHECK_EVSHARE (4c6ca367fd). Same recipe as ../build_container.sh; new names only.
set -euo pipefail
here=/home/singjc/.jcode/scratch/evshare-20261126/container
name=opendia-evshare-20261126
mkdir -p "$here"
cd "$here"
rm -f payload.tar image_check.txt  # leftovers of the first attempt (help check used --help; the flag is advanced)
test ! -e "$name.sif"
repo=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment
build=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build
test -z "$(git -C "$repo" status --porcelain)"
git -C "$repo" rev-parse HEAD > source_commit.txt
grep -q '^4c6ca367fd' source_commit.txt
sha256sum "$build/bin/OpenDIA" > binary.sha256
cp /home/singjc/.jcode/scratch/dia2/ctr/pd/files.txt files.txt
cp /home/singjc/.jcode/scratch/dia2/ctr/pd/Dockerfile Dockerfile
cp -a /home/singjc/.jcode/scratch/dia2/ctr/pd/TOOLS TOOLS
printf 'FROM %s:base\nCOPY TOOLS /home/singjc/Documents/github/OpenMS-r1-sparse-experiment/share/OpenMS/TOOLS\n' "$name" > Dockerfile.tools
tar cf payload.tar -P --transform 's,^/,,' --files-from files.txt
docker build -t "$name:base" .
docker build -t "$name:tools" -f Dockerfile.tools .
docker run --rm "$name:tools" bash -c 'set -euo pipefail; sha256sum "$(which OpenDIA)"; OpenDIA --help > /tmp/h.txt 2>&1; head -3 /tmp/h.txt; OpenDIA --helphelp 2>&1 | grep -c -- -evidence_share' > image_check.txt
local_hash=$(cut -d' ' -f1 binary.sha256)
grep -q "^$local_hash " image_check.txt
singularity build "$name.sif" "docker-daemon://$name:tools" > sif_build.log 2>&1
sha256sum "$name.sif" > container.sha256
singularity exec --cleanenv "$name.sif" bash -c 'set -euo pipefail; sha256sum "$(which OpenDIA)"; OpenDIA --help > /tmp/h.txt 2>&1; head -3 /tmp/h.txt; OpenDIA --helphelp 2>&1 | grep -c -- -evidence_share' > sif_binary.sha256
grep -q "^$local_hash " sif_binary.sha256
rm -f payload.tar
cat binary.sha256 container.sha256 sif_binary.sha256
