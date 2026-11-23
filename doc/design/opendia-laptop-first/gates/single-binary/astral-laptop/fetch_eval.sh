#!/bin/bash
# Fetch the 15 Astral laptop-profile runs from codon and evaluate them with the frozen eval_run.py.
set -e
H=/home/singjc/.jcode/scratch/dia2/gap7/astral
SSH="ssh -o BatchMode=yes -o ConnectTimeout=30 -o ControlMaster=no -o ControlPath=none codon"
W=/nfs/research/juan/singj/opendia/astral-laptop-20261123/runs
$SSH "become pst_prd bash -lc 'cd $W && tar cf - raw*_s*.log raw*_s*/*.precursors.tsv'" 2>/dev/null | tar xf - -C $H
EV=/home/singjc/Documents/github/OpenMS-wp00-design-docs/doc/design/opendia-laptop-first/gates/single-binary/workflow/eval_run.py
LIB=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp
python3 $EV $LIB $(ls $H/raw*_s*/LFQ*.precursors.tsv | sort) > $H/EVAL.txt 2>&1
for f in $(ls $H/raw*_s*.log | sort); do
  n=$(basename $f .log)
  w=$(grep -a "Elapsed (wall" $f | sed 's/.*: //'); m=$(grep -a "Maximum resident" $f | sed 's/.*: //')
  u=$(grep -a "User time" $f | sed 's/.*: //'); s=$(grep -a "System time" $f | sed 's/.*: //'); c=$(grep -a "Percent of CPU" $f | sed 's/.*: //')
  node=$(head -1 $f | cut -d' ' -f2); ex=$(grep -a "^exit" $f)
  echo "$n wall=$w rss_kb=$m user_s=$u sys_s=$s cpu=$c node=$node $ex"
done > $H/COST.txt
cat $H/EVAL.txt $H/COST.txt
