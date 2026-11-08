#!/bin/bash
# Diagnostic (step A): raw 01 calibration with the frozen settings plus the experimental anchor peak/trace outputs.
# Default outputs must be byte-identical to ../window300k. Not an evaluation.
B=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107
O=$B/astral-design-raw01-20261107/anchor-peaks
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
export PATH="$DOTNET_ROOT:$PATH"
cd $O
echo running > STATUS
/usr/bin/time -v /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIACalibrate \
  -in /home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw \
  -lib $B/astral-entrapment-20261107/lib/library.pqp \
  -panel_targets 40000 -window_panel_targets 300000 -out cal.json -out_anchors anchors.tsv -out_diagnostics diag.tsv \
  -out_anchor_peaks anchor_peaks.tsv -out_anchor_traces anchor_traces.tsv > stdout.log 2> stderr.log
echo $? > EXIT
{ cmp cal.json ../window300k/cal.json && cmp anchors.tsv ../window300k/anchors.tsv && cmp diag.tsv ../window300k/diag.tsv && echo "default calibration outputs IDENTICAL to window300k"; } > CMP.txt 2>&1
echo done > STATUS
