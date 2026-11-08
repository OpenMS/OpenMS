#!/bin/bash
# Diagnostic re-score of the SPENT raw 02 with the frozen stage-2 design (014cb0dbf1, raw 01 RESULT.md).
# Not a gate: stage-2 features were chosen from the raw 02 post-mortem, so this is optimistic by construction.
# Reuses raw 02 calibration, seeds and windows read-only; only the -stage2 flags differ from raw 02 run.sh.
set -e
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
R=$B/astral-entrapment-raw02-20261108
O=$B/astral-raw02-stage2-check-20261108
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_02.raw
LIB=$B/astral-entrapment-20261107/lib/library.pqp
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
N=$(($(wc -l < $R/all_windows.tsv) - 1))
echo cnn > STATUS
/usr/bin/time -v env PATH=/usr/bin:/bin $BIN/OpenDIACNN -mode crossfit -in $IN -lib $LIB -windows $R/all_windows.tsv -seeds $R/seeds.tsv -model $O/crossfit/model -out $O/crossfit/crossfit.tsv -max_samples $N -max_tensor_gib 24 -threads 5 -folds 3 -seed 42 -epochs 40 -patience 6 -batch_size 32 -paired_tdc -pooled_tdc -rt_calibration $R/cal/cal.json -apply_calibration_coordinates -ss_iterations 0 -stage2 -stage2_features_out $O/crossfit/features.tsv > crossfit/stdout.log 2> crossfit/stderr.log || { echo FAILED > STATUS; exit 1; }
{ cmp <(cut -f1-6,9 crossfit/crossfit.tsv) <(cut -f1-6,9 $R/crossfit/crossfit.tsv) && echo "logit/fold/window columns IDENTICAL to raw 02 evaluation"; } > CMP.txt 2>&1 || true
~/.venvs/redeem-external-dev-uv-20261022/bin/python evaluate.py crossfit/crossfit.tsv $LIB > EVALUATION.txt 2>&1
echo done > STATUS
