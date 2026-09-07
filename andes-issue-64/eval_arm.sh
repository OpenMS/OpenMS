#!/usr/bin/env bash
# eval_arm.sh <name> — score one arm with the repo's benchmark scripts (quick tier, one fraction)
set -uo pipefail
NAME="$1"
B=/home/user/andes-bench
G=/home/user/andes/docs/benchmarks/glyco
cd "$B/$NAME"
PIN="$NAME.glyco.pin"; PSMS="$NAME.t.psms"
{
  echo "=== $NAME: eval_yield"
  python3 $G/eval_yield.py "$PIN" "$PSMS"
  echo; echo "=== $NAME: eval_entrap (1:1 shuffled-self database)"
  python3 $G/eval_entrap.py "$PIN" "$PSMS" 0.01 "$B/databases/mouse_entrap.fasta"
  echo; echo "=== $NAME: score_vs_truth pGlyco2"
  python3 $G/score_vs_truth.py --run MouseLiver-Z-T-1 --buckets buckets_pglyco2.tsv $G/truth/pglyco2_mouse_liver.tsv.gz "$PIN" "$PSMS"
  echo; echo "=== $NAME: score_vs_truth MSFragger"
  python3 $G/score_vs_truth.py --run MouseLiver-Z-T-1 --buckets buckets_msfragger.tsv $G/truth/msfragger_mouse_liver.tsv.gz "$PIN" "$PSMS"
  echo; echo "=== $NAME: agreement pGlyco2"
  python3 $G/agreement.py --run MouseLiver-Z-T-1 $G/truth/pglyco2_mouse_liver.tsv.gz "$PSMS"
  echo; echo "=== $NAME: agreement MSFragger"
  python3 $G/agreement.py --run MouseLiver-Z-T-1 $G/truth/msfragger_mouse_liver.tsv.gz "$PSMS"
} > eval.txt 2>&1
echo "wrote $B/$NAME/eval.txt"
