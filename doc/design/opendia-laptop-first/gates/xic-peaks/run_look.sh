#!/bin/bash
# xic look, raw01: xic search (timed) then legacy byte-identity check. See PLAN_XIC_LOOK.md.
set -o pipefail
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
X=$B/xic-peaks-20261112
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
R=$B/rtlin-20261111/raw01
LIB=$B/astral-entrapment-20261107/lib/library.pqp
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
st() { echo "$1" > $X/STATUS; date -u +"%FT%TZ $1" >> $X/STATUS.log; }
cd $X; sha256sum $BIN/OpenDIA > binary.sha256
st xic01
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration $R/cal/cal.json -out raw01/xic/search.tsv -features_out raw01/xic/search.features.tsv \
  -min_transitions 1 -peak_picking xic -xic_max_peaks 5 > raw01/xic/stdout.log 2> raw01/xic/stderr.log || { st failed-xic01; exit 1; }
st legacy01
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration $R/cal/cal.json -out raw01/legacy/search.tsv -features_out raw01/legacy/search.features.tsv \
  -min_transitions 1 -local_peak_split -max_apexes_per_group 5 -min_apex_separation 10 > raw01/legacy/stdout.log 2> raw01/legacy/stderr.log || { st failed-legacy01; exit 1; }
cmp raw01/legacy/search.features.tsv $R/pipeline/search.features.tsv > raw01/legacy/CMP.txt 2>&1 && echo IDENTICAL >> raw01/legacy/CMP.txt
cmp raw01/legacy/search.tsv $R/pipeline/search.tsv >> raw01/legacy/CMP.txt 2>&1 && echo IDENTICAL_TSV >> raw01/legacy/CMP.txt
rm -f raw01/legacy/search.tsv
st done
