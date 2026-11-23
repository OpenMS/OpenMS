#!/bin/bash
# PRECHECK_ASTQ evaluation: fetch ladder outputs, build eval_pdec.py layout <P>_<arm>_s<seed>, run per raw.
set -e
H=/home/singjc/.jcode/scratch/dia2/gap7/astq; mkdir -p $H/fetch $H/lay
SSH="ssh -o BatchMode=yes -o ConnectTimeout=20 -o ServerAliveInterval=10 -o ControlMaster=no -o ControlPath=none codon"
$SSH "become pst_prd bash -lc 'cd /nfs/research/juan/singj/opendia/astq-20261123/runs && tar cf - *.log */LFQ*.precursors.tsv'" 2>/dev/null | tar xf - -C $H/fetch
L0=/home/singjc/.jcode/scratch/dia2/gap7/astral_gbmfix
for R in 01 02; do for s in 0 1 2 3 4; do
  mkdir -p $H/lay/raw${R}_L0_s$s; ln -sf $L0/raw${R}_s$s/LFQ_Astral_DIA_15min_50ng_Human_$R.precursors.tsv $H/lay/raw${R}_L0_s$s/
  for A in L1 L2 L3 L4 L5; do mkdir -p $H/lay/raw${R}_${A}_s$s; ln -sf $H/fetch/raw${R}_${A}_s$s/LFQ_Astral_DIA_15min_50ng_Human_$R.precursors.tsv $H/lay/raw${R}_${A}_s$s/; done
done; done
echo "== anchors"
cmp $H/fetch/raw01_A0_s0/LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv $L0/raw01_s0/LFQ_Astral_DIA_15min_50ng_Human_01.precursors.tsv && echo "A0 raw01_s0 SAME"
cmp $H/fetch/raw02_A0_s2/LFQ_Astral_DIA_15min_50ng_Human_02.precursors.tsv $L0/raw02_s2/LFQ_Astral_DIA_15min_50ng_Human_02.precursors.tsv && echo "A0 raw02_s2 SAME"
grep -H "^exit" $H/fetch/*.log | grep -v "exit 0" || echo "all exit 0"
P=/home/singjc/Documents/github/pyprophet/.venv/bin/python
EV=/home/singjc/Documents/github/OpenMS-wp00-design-docs/doc/design/opendia-laptop-first/gates/single-binary/pdec/eval_pdec.py
LIB=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp
for R in raw01 raw02; do echo "== $R"; $P $EV $LIB $H/lay $R; done
