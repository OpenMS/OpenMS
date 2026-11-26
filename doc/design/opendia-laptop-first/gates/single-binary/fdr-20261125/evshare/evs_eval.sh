#!/bin/bash
# PRECHECK_EVSHARE evaluation (written before any output). Applies the frozen gate rule, fetches every run table
# with a sha256 check against codon, builds <set>_<arm>_s<seed> layouts and runs the frozen evaluators.
set -euo pipefail
H=/home/singjc/.jcode/scratch/evshare-20261126/eval
D=/home/singjc/Documents/github/OpenMS-wp00-design-docs/doc/design/opendia-laptop-first/gates/single-binary/fdr-20261125
SSH="ssh -n -o BatchMode=yes -o ConnectTimeout=30 -o ServerAliveInterval=15 -o ControlMaster=no -o ControlPath=none codon"
C=/nfs/research/juan/singj/opendia; W=$C/evshare-20261126
PY=/home/singjc/Documents/github/pyprophet/.venv/bin/python
mkdir -p $H
# Gate (frozen rule): all three no-flag cells SAME on every file -> existing G0 references, else fresh G0 runs.
$SSH "become pst_prd bash -c 'cd $W/runs && grep -H -E \"GATE|^exit\" raw03_G0_s0.log D_G0_s0.log H_G0_s0.log'" 2>/dev/null | tee $H/gate.txt
ndiff=$(grep -c "GATE DIFF" $H/gate.txt || true); nsame=$(grep -c "GATE SAME" $H/gate.txt || true)
[ "$nsame" -gt 0 ] || { echo "gate incomplete"; exit 1; }
echo "gate same $nsame diff $ndiff" | tee -a $H/gate.txt
if [ "$ndiff" = 0 ]; then
  g0() { case $1 in raw03) echo astral-gbmfix-20261123/runs/raw03_s$2;; D|H) echo entlib-gbt-20261122/runs/${1}_gbt_s$2;; *) echo evshare-20261126/runs/${1}_G0_s$2;; esac; }
else
  g0() { echo evshare-20261126/runs/${1}_G0_s$2; }
fi
# Fetch: every per-run table (A and the G0 reference), sha256 checked against the remote file.
want=$(for P in raw03 D H N1 N2; do for s in 0 1 2 3 4; do echo evshare-20261126/runs/${P}_A_s$s; g0 $P $s; done; done)
$SSH "become pst_prd bash -c 'cd $C && for d in $(echo $want); do for f in \$d/*.precursors.tsv; do [ \$(basename \$f) = precursors.tsv ] || sha256sum \$f; done; done'" 2>/dev/null > $H/remote.sha256
echo "remote tables: $(wc -l < $H/remote.sha256) (expect 50)"
while read -r sum rel; do
  dst=$H/fetch/$rel
  if [ -f "$dst" ] && [ "$(sha256sum < "$dst" | cut -d' ' -f1)" = "$sum" ]; then continue; fi
  mkdir -p "$(dirname "$dst")"
  for try in 1 2 3; do
    ssh -o BatchMode=yes -o ConnectTimeout=30 -o ServerAliveInterval=15 -o ControlMaster=no -o ControlPath=none codon \
      "become pst_prd bash -c 'cat $C/$rel'" < /dev/null 2>/dev/null > "$dst.part" || true
    if [ "$(sha256sum < "$dst.part" | cut -d' ' -f1)" = "$sum" ]; then mv "$dst.part" "$dst"; break; fi
    echo "retry $try $rel"
  done
  [ -f "$dst" ] || { echo "FAILED $rel"; exit 1; }
done < $H/remote.sha256
# Layout. A missing table (failed run) is left out of the layout and printed as MISSING. The evaluator then cannot see it,
# so RESULT_EVSHARE adds it by hand as N = 0 (N @ FDP_c 1 % = 0) in the median, as the frozen rule says.
rm -rf $H/lay; mkdir -p $H/lay
for P in raw03 D H N1 N2; do for s in 0 1 2 3 4; do for A in G0 A; do
  src=$H/fetch/$([ $A = G0 ] && g0 $P $s || echo evshare-20261126/runs/${P}_A_s$s)
  d=$H/lay/${P}_${A}_s$s
  ls $src/*.precursors.tsv > /dev/null 2>&1 || { echo "MISSING ${P}_${A}_s$s (N = 0)"; continue; }
  mkdir -p $d; for f in $src/*.precursors.tsv; do [ "$(basename $f)" = precursors.tsv ] || ln -sf $f $d/; done
done; done; done
LA=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp
LH=/home/singjc/.jcode/scratch/blend-20261125/eval/hela_r5.predicted_decoys.pqp
sha256sum $D/stepA/eval_pdec_relabel.py $D/stepA/eval_pdec_instr.py $D/stepA/label_audit_astral.tsv $D/stepA/label_audit_hela.tsv $LA $LH
{ echo "== PRIMARY (relabel)"
  echo "== raw03"; RELABEL=$D/stepA/label_audit_astral.tsv $PY $D/stepA/eval_pdec_relabel.py $LA $H/lay raw03
  for P in D H N1 N2; do echo "== $P"; RELABEL=$D/stepA/label_audit_hela.tsv $PY $D/stepA/eval_pdec_relabel.py $LH $H/lay $P; done
  echo "== SECONDARY (instr, no relabel)"
  echo "== raw03"; $PY $D/stepA/eval_pdec_instr.py $LA $H/lay raw03
  for P in D H N1 N2; do echo "== $P"; $PY $D/stepA/eval_pdec_instr.py $LH $H/lay $P; done; } | tee $H/EVAL_EVSHARE.txt
# Cost (descriptive).
$SSH "become pst_prd bash -c 'cd $W/runs && for f in *.log; do echo \$f \$(grep -E \"^node|Elapsed|Maximum resident|^exit\" \$f | sed \"s/.*: //\" | tr \"\\n\" \" \"); done'" 2>/dev/null > $H/COST_EVSHARE.txt
cat $H/COST_EVSHARE.txt
