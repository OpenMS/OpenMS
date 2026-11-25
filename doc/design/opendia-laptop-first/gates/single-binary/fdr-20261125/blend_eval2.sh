#!/bin/bash
# PRECHECK_BLEND evaluation, stage 2 (no change to the frozen rule or evaluator): the G0 reference fetch in blend_eval.sh
# broke mid-stream (tar: Unexpected EOF). Gate result was 0 diffs, so G0 = the existing references. This refetches each
# reference file individually with a size check against the remote, then builds the same layouts and runs eval_pdec.py.
set -euo pipefail
H=/home/singjc/.jcode/scratch/blend-20261125/eval
SSH="ssh -o BatchMode=yes -o ConnectTimeout=30 -o ServerAliveInterval=15 -o ControlMaster=no -o ControlPath=none codon"
C=/nfs/research/juan/singj/opendia
list=$($SSH "become pst_prd bash -c 'cd $C && stat -c \"%s %n\" astral-gbmfix-20261123/runs/raw0[123]_s[0-4]/*.precursors.tsv entlib-gbt-20261122/runs/[DH]_gbt_s[0-4]/*.precursors.tsv'" 2>/dev/null)
echo "$list" | wc -l
echo "$list" | while read -r size rel; do
  dst=$H/ref/$rel
  if [ -f "$dst" ] && [ "$(stat -c %s "$dst")" = "$size" ]; then continue; fi
  mkdir -p "$(dirname "$dst")"
  for try in 1 2 3; do
    $SSH "become pst_prd bash -c 'cat $C/$rel'" 2>/dev/null > "$dst.part" || true
    if [ "$(stat -c %s "$dst.part")" = "$size" ]; then mv "$dst.part" "$dst"; break; fi
    echo "retry $try $rel"
  done
  [ -f "$dst" ] && [ "$(stat -c %s "$dst")" = "$size" ] || { echo "FAILED $rel"; exit 1; }
done
echo "refs complete: $(find $H/ref -name '*.precursors.tsv' ! -name precursors.tsv | wc -l) run tables"
g0() { case $1 in raw*) echo $H/ref/astral-gbmfix-20261123/runs/${1}_s$2;; *) echo $H/ref/entlib-gbt-20261122/runs/${1}_gbt_s$2;; esac; }
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
LH=$H/hela_r5.predicted_decoys.pqp
LS=$($SSH "become pst_prd bash -c 'stat -c %s $C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'" 2>/dev/null)
[ -f $LH ] && [ "$(stat -c %s $LH)" = "$LS" ] || $SSH "become pst_prd bash -c 'cat $C/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'" 2>/dev/null > $LH
[ "$(stat -c %s $LH)" = "$LS" ] || { echo "HeLa library fetch incomplete"; exit 1; }
for P in raw01 raw02 raw03; do echo "== $P"; $PY $EV $LA $H/lay $P; done
for R in 1 2 3; do for s in 0 1 2 3 4; do for A in G0 B5; do
  d=$H/lay/astral_${A}_s$((10 * R + s)); mkdir -p $d; ln -sf $H/lay/raw0${R}_${A}_s$s/*.precursors.tsv $d/
done; done; done
echo "== astral (15 runs)"; $PY $EV $LA $H/lay astral
for P in D H; do echo "== $P"; $PY $EV $LH $H/lay $P; done
