#!/bin/bash
# PRECHECK_BLEND evaluation: fetch outputs, build eval_pdec.py layouts (<set>_<arm>_s<seed>), run the frozen evaluator.
# G0 = existing references if all 3 gates are SAME, else the fresh G0 runs (frozen rule).
set -euo pipefail
H=/home/singjc/.jcode/scratch/blend-20261125/eval; mkdir -p $H/fetch $H/lay
SSH="ssh -o BatchMode=yes -o ConnectTimeout=30 -o ServerAliveInterval=15 -o ControlMaster=no -o ControlPath=none codon"
C=/nfs/research/juan/singj/opendia
$SSH "become pst_prd bash -c 'cd $C/blend-20261125/runs && tar cf - *.log */*.precursors.tsv'" 2>/dev/null | tar xf - -C $H/fetch
grep -h "GATE" $H/fetch/*_G0_s0.log | sort | uniq -c
grep -H "^exit" $H/fetch/*.log | grep -v "exit 0" || echo "all exit 0"
ndiff=$(cat $H/fetch/raw01_G0_s0.log $H/fetch/D_G0_s0.log $H/fetch/H_G0_s0.log | grep -c "GATE DIFF" || true)
echo "gate diffs: $ndiff"
if [ "$ndiff" = 0 ]; then
  mkdir -p $H/ref
  $SSH "become pst_prd bash -c 'cd $C && tar cf - astral-gbmfix-20261123/runs/raw0[123]_s[0-4]/*.precursors.tsv entlib-gbt-20261122/runs/[DH]_gbt_s[0-4]/*.precursors.tsv'" 2>/dev/null | tar xf - -C $H/ref
  g0() { case $1 in raw*) echo $H/ref/astral-gbmfix-20261123/runs/${1}_s$2;; *) echo $H/ref/entlib-gbt-20261122/runs/${1}_gbt_s$2;; esac; }
else
  g0() { echo $H/fetch/${1}_G0_s$2; }
fi
for P in raw01 raw02 raw03 D H; do for s in 0 1 2 3 4; do
  for A in G0 B5; do
    d=$H/lay/${P}_${A}_s$s; mkdir -p $d
    src=$([ $A = G0 ] && g0 $P $s || echo $H/fetch/${P}_B5_s$s)
    for f in $src/*.precursors.tsv; do [ "$(basename $f)" = precursors.tsv ] || ln -sf $f $d/; done
  done
done; done
PY=/home/singjc/Documents/github/pyprophet/.venv/bin/python
EV=/home/singjc/Documents/github/OpenMS-wp00-design-docs/doc/design/opendia-laptop-first/gates/single-binary/pdec/eval_pdec.py
sha256sum $EV
LA=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp
LH=/home/singjc/.jcode/scratch/blend-20261125/eval/hela_r5.predicted_decoys.pqp
[ -s $LH ] || $SSH "become pst_prd bash -c 'cat $C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'" > $LH
for P in raw01 raw02 raw03; do echo "== $P"; $PY $EV $LA $H/lay $P; done
# All three raws pooled: a layout with prefix "astral", seed = 10 * raw + seed.
for R in 1 2 3; do for s in 0 1 2 3 4; do for A in G0 B5; do
  d=$H/lay/astral_${A}_s$((10 * R + s)); mkdir -p $d; ln -sf $H/lay/raw0${R}_${A}_s$s/*.precursors.tsv $d/
done; done; done
echo "== astral (15 runs)"; $PY $EV $LA $H/lay astral
for P in D H; do echo "== $P"; $PY $EV $LH $H/lay $P; done
