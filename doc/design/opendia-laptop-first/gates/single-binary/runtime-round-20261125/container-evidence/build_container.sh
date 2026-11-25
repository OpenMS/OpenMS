#!/usr/bin/env bash
# Container for the runtime round (mem5 + rt fixes, bdb0159ff4). Same recipe as mem5-20261125/build_container.sh, new names only.
set -euo pipefail
here=/home/singjc/.jcode/scratch/rt2-20261125/container
name=opendia-rt2-20261125
mkdir -p "$here"
cd "$here"
test ! -e payload.tar
test ! -e "$name.sif"
repo=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment
build=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build
test -z "$(git -C "$repo" status --porcelain)"
git -C "$repo" rev-parse HEAD > source_commit.txt
sha256sum "$build/bin/OpenDIA" > binary.sha256
cp /home/singjc/.jcode/scratch/dia2/ctr/pd/files.txt files.txt
cp /home/singjc/.jcode/scratch/dia2/ctr/pd/Dockerfile Dockerfile
cp -a /home/singjc/.jcode/scratch/dia2/ctr/pd/TOOLS TOOLS
printf 'FROM %s:base\nCOPY TOOLS /home/singjc/Documents/github/OpenMS-r1-sparse-experiment/share/OpenMS/TOOLS\n' "$name" > Dockerfile.tools
tar cf payload.tar -P --transform 's,^/,,' --files-from files.txt
docker build -t "$name:base" .
docker build -t "$name:tools" -f Dockerfile.tools .
docker run --rm "$name:tools" bash -c 'set -euo pipefail; sha256sum "$(which OpenDIA)"; OpenDIA --help > /tmp/h.txt 2>&1; head -3 /tmp/h.txt' > image_check.txt
local_hash=$(cut -d' ' -f1 binary.sha256)
grep -q "^$local_hash " image_check.txt
singularity build "$name.sif" "docker-daemon://$name:tools" > sif_build.log 2>&1
sha256sum "$name.sif" > container.sha256
singularity exec --cleanenv "$name.sif" bash -c 'set -euo pipefail; sha256sum "$(which OpenDIA)"; OpenDIA --help > /tmp/h.txt 2>&1; head -3 /tmp/h.txt' > sif_binary.sha256
grep -q "^$local_hash " sif_binary.sha256
cat binary.sha256 container.sha256
