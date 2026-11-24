#!/usr/bin/env bash
set -euo pipefail
work=/home/singjc/.jcode/scratch/compact-20261124
sif="$work/container/opendia-compact-parquet-20261124.sif"
bins=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/src/tests/class_tests/bin
fixture=/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/tests/topp/TargetedFileConverter_10_input.pqp
mkdir "$work/container-smoke"
sha256sum "$bins/OpenDIAWorkflow_test" "$bins/OpenDIALibraryParquet_test" "$fixture" > "$work/container-smoke/artifacts.sha256"
for name in OpenDIAWorkflow_test OpenDIALibraryParquet_test; do
 mkdir "$work/container-smoke/$name"
 singularity exec --cleanenv --no-home --bind "$bins:/validation-bin:ro" --bind "$fixture:/fixture.pqp:ro" --bind "$work/container-smoke/$name:/validation" --pwd /validation "$sif" "/validation-bin/$name" /validation /fixture.pqp > "$work/container-smoke/$name.log" 2>&1
 cat "$work/container-smoke/$name.log"
done
