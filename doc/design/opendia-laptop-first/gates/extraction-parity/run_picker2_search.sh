#!/bin/bash
# Same calibration as picker2/ (k estimated 1.69); search variants: cap off (co-elution fix only) and cap 2.5 (reference).
set -o pipefail
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107; X=$B/extraction-parity-20261112
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin; LIB=$B/astral-entrapment-20261107/lib/library.pqp
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet; export PATH="$DOTNET_ROOT:$PATH"
run() { D=$X/$1; mkdir -p $D; cd $D; echo search > STATUS
  /usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration $X/picker2/cal/cal.json -out search.tsv -features_out search.features.tsv \
    -min_transitions 1 -peak_picking xic -xic_max_peaks 5 -xic_boundary_half_widths $2 -xic_traces_out traces.tsv -xic_traces_precursors $X/dump_ids.txt \
    > search.stdout.log 2> search.stderr.log && echo done > STATUS || echo failed > STATUS; }
run picker2 -1 &

wait
