#!/usr/bin/env bash
set -euo pipefail
here=/home/singjc/.jcode/scratch/compact-20261124/container
mkdir -p "$here"
cd "$here"
test ! -e payload.tar
test ! -e opendia-compact-parquet-20261124.sif
repo=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment
build=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build
test -z "$(git -C "$repo" status --porcelain)"
git -C "$repo" rev-parse HEAD > source_commit.txt
sha256sum "$build/bin/OpenDIA" > binary.sha256
cp /home/singjc/.jcode/scratch/dia2/ctr/pd/files.txt files.txt
cp /home/singjc/.jcode/scratch/dia2/ctr/pd/Dockerfile Dockerfile
cp -a /home/singjc/.jcode/scratch/dia2/ctr/pd/TOOLS TOOLS
printf 'FROM opendia-compact-parquet:20261124\nCOPY TOOLS /home/singjc/Documents/github/OpenMS-r1-sparse-experiment/share/OpenMS/TOOLS\n' > Dockerfile.tools
tar cf payload.tar -P --transform 's,^/,,' --files-from files.txt
docker build -t opendia-compact-parquet:20261124 .
docker build -t opendia-compact-parquet:20261124-tools -f Dockerfile.tools .
docker run --rm opendia-compact-parquet:20261124-tools bash -c 'set -euo pipefail; sha256sum "$(which OpenDIA)"; OpenDIA --help > /tmp/opendia_help.txt 2>&1; head -14 /tmp/opendia_help.txt' > image_check.txt
local_hash=$(cut -d' ' -f1 binary.sha256)
grep -q "^$local_hash " image_check.txt
singularity build opendia-compact-parquet-20261124.sif docker-daemon://opendia-compact-parquet:20261124-tools > sif_build.log 2>&1
sha256sum opendia-compact-parquet-20261124.sif > container.sha256
singularity exec --cleanenv opendia-compact-parquet-20261124.sif bash -c 'set -euo pipefail; sha256sum "$(which OpenDIA)"; OpenDIA --help > /tmp/opendia_sif_help.txt 2>&1; head -14 /tmp/opendia_sif_help.txt' > sif_binary.sha256
grep -q "^$local_hash " sif_binary.sha256
cat binary.sha256 container.sha256
