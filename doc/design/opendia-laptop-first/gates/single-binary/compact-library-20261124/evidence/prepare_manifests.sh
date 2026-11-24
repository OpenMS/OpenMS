#!/usr/bin/env bash
# Root-owned preflight construction. Run only after the plan is frozen and pushed,
# the isolated SIF passes its binary hash checks, and raw01 local gates pass.
set -euo pipefail
repo=/home/singjc/Documents/github/OpenMS-wp00-design-docs
bundle="$repo/doc/design/opendia-laptop-first/gates/single-binary/astral-full-parquet"
work=/nfs/research/juan/singj/opendia/astral-parquet-20261124
container=/nfs/research/juan/singj/opendia/opendia-compact-parquet-20261124.sif
local_sif=/home/singjc/.jcode/scratch/compact-20261124/container/opendia-compact-parquet-20261124.sif
baseline=/nfs/research/juan/singj/opendia/astral-full-20261123
raw=/nfs/research/juan/singj/opendia/proteobench-astral-entrapment/transfer/raw
cd "$bundle"
test -s PREFLIGHT_APPROVED.md
test -s source_provenance.txt
test -s "$local_sif"
test ! -e preflight.sha256
test ! -e runtime.sha256
# Keep the runtime manifest small: executable scripts only, no PQP/raw/container scans.
for file in bootstrap.sbatch run.sbatch evaluate.sbatch submit.sh compare_resources.py; do
 sha256sum "$file"
done > runtime.sha256
cp runtime.sha256 preflight.sha256
for file in PLAN_PARQUET.md PREFLIGHT_APPROVED.md source_provenance.txt test_compare_resources.py runtime.sha256; do
 sha256sum "$file"
done >> preflight.sha256
hash=$(sha256sum "$local_sif" | cut -d' ' -f1)
printf '%s  %s\n' "$hash" "$container" >> preflight.sha256
printf '%s  %s\n' \
 200dea366ea358f9f97674715c070a6bd4ade8a839d2310fbf214fbf845847c7 "$baseline/library/search6.pqp" \
 95a19dec2d169d6b19f973b27eda9763df26fcad69d1a3847115cb76186c49a0 "$baseline/library/scoring12.pqp" \
 65d7c288949ac22290e93b6a4ad69d6eae5f944f8c0cbdc93dba75729f97a125 "$raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw" \
 d945298ce44cd49ef2337636f5e0cf4b7ec3a139504e09d5fcb5ec1d4a03215a "$raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw" \
 bee9061cb1528b2e39f5baa05fbe94806e7c3d28578ba0c625afe256f0566e10 "$raw/LFQ_Astral_DIA_15min_50ng_Human_03.raw" >> preflight.sha256
printf 'Prepared manifests for %s. This script does not authorize submission.\n' "$work"
