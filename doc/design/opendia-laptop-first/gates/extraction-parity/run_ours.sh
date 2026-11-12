#!/bin/bash
# OpenDIA xic search on raw01 (same command as xic-peaks-20261112/run_look.sh) plus the native-scan trace dump of
# every precursor DIA-NN reported (7553 ids, q <= 0.5) for the trace-by-trace extraction parity check.
set -o pipefail
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
X=$B/extraction-parity-20261112
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
R=$B/rtlin-20261111/raw01
LIB=$B/astral-entrapment-20261107/lib/library.pqp
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $X; mkdir -p ours; sha256sum $BIN/OpenDIA > ours/binary.sha256
echo running > ours/STATUS
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration $R/cal/cal.json -out ours/search.tsv -features_out ours/search.features.tsv \
  -min_transitions 1 -peak_picking xic -xic_max_peaks 5 -xic_traces_out ours/traces.tsv -xic_traces_precursors dump_ids.txt \
  > ours/stdout.log 2> ours/stderr.log && echo done > ours/STATUS || echo failed > ours/STATUS
cmp ours/search.features.tsv $B/xic-peaks-20261112/raw01/xic/search.features.tsv > ours/CMP.txt 2>&1 && echo IDENTICAL >> ours/CMP.txt
