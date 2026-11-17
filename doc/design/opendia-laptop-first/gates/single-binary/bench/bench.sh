#!/bin/bash
# Fair speed/memory benchmark: OpenDIA vs DIA-NN 2.2 on the same laptop, same raw (.raw read directly by both),
# 8 threads each, run one after the other. DIA-NN gets the SDK-8 .NET via a read-only mount (no change to the image).
# DIA-NN settings as diann-pairlib-20261113 (its own decoys, same target library); OpenDIA default workflow.
B=/home/singjc/.jcode/scratch/e2e/bench
RAWD=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw
G=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
cd $B
for R in 01; do
  echo "== DIA-NN raw$R $(date -Is)"
  /usr/bin/time -v docker run --rm --cpus 8 -u $(id -u):$(id -g) -e HOME=/w/tmp -e DOTNET_ROOT=/dotnet \
    -e DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1 -v /home/singjc/.jcode/scratch/dotnet-sdk8:/dotnet:ro -v $RAWD:/in:ro -v $B:/w -v /usr/bin/time:/usr/bin/time:ro \
    --entrypoint /bin/sh diann:2.2.0 -c "export PATH=/dotnet:\$PATH; /usr/bin/time -v /usr/diann-2.2.0/diann-linux \
    --f /in/LFQ_Astral_DIA_15min_50ng_Human_$R.raw --lib /w/lib_targets.tsv.skyline.speclib --threads 8 --verbose 1 \
    --out /w/diann_raw$R.parquet --qvalue 1 --unimod4 --min-pr-charge 1 --max-pr-charge 5 --temp /w/tmp" > diann_raw$R.log 2>&1
  echo "== OpenDIA raw$R $(date -Is)"
  DOTNET_ROOT=/home/singjc/.dotnet/runtime8 /usr/bin/time -v /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIAWorkflow \
    -in $RAWD/LFQ_Astral_DIA_15min_50ng_Human_$R.raw -lib $G/gap-20261114/armE_raw01/library.pqp \
    -scoring_lib $G/lib12full/pairlib-view/lib12/library12_full.pqp -out_dir $B/opendia_raw$R -threads 8 -seed 0 > opendia_raw$R.log 2>&1
done
echo "== done $(date -Is)"
