#!/bin/bash
# Recalibrate raw01 with the fitted fragment tolerance (OpenDIACalibration: strongest-centroid errors at apex +-1
# spectra, signal/background crossing), same options as rtlin-20261111/chain_data.sh, then the xic search with
# the trace dump for the 7553 DIA-NN precursors, then trace parity.
set -o pipefail
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107; X=$B/extraction-parity-20261112; D=$X/fitppm
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin; LIB=$B/astral-entrapment-20261107/lib/library.pqp
IN=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet; export PATH="$DOTNET_ROOT:$PATH"
mkdir -p $D/cal; cd $D; sha256sum $BIN/OpenDIA $BIN/OpenDIACalibrate > binary.sha256
echo cal > STATUS
/usr/bin/time -v $BIN/OpenDIACalibrate -in $IN -lib $LIB -panel_targets 40000 -window_panel_targets 300000 -out cal/cal.json \
  -out_anchors cal/anchors.tsv -out_diagnostics cal/diag.tsv > cal/stdout.log 2> cal/stderr.log || { echo failed-cal > STATUS; exit 1; }
echo search > STATUS
/usr/bin/time -v $BIN/OpenDIA -in $IN -lib $LIB -calibration cal/cal.json -out search.tsv -features_out search.features.tsv \
  -min_transitions 1 -peak_picking xic -xic_max_peaks 5 -xic_traces_out traces.tsv -xic_traces_precursors $X/dump_ids.txt \
  > search.stdout.log 2> search.stderr.log || { echo failed-search > STATUS; exit 1; }
echo parity > STATUS
cd $X && $JCODE_SCRATCH_DIR/torchenv/.venv/bin/python parity2.py fitppm/TRACE_PARITY_raw01.txt fitppm > fitppm/parity.log 2>&1 || { echo failed-parity > fitppm/STATUS; exit 1; }
echo done > $D/STATUS
