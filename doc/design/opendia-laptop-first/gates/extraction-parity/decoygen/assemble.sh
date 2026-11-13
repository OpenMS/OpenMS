#!/bin/bash
# PLAN_DECOYGEN: assemble each arm's run root (arm<X>/home) and its FROZEN.sha256. The PAIRLIB chain checks FROZEN.sha256
# in the run root before it starts. Scripts: the PAIRLIB scripts as frozen (parent.py with PAIRLIB amendment 2), with
# build_pqp12.py replaced and null_rates.py and pooled_rho.py added (decoygen-20261113/scripts).
set -e
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/decoygen-20261113
P=../pairlib-20261113
for a in A B C D; do
  H=arm$a/home
  [ "$(cat $H/lib/CONVERT_EXIT)" = 0 ] || { echo "arm $a: library not converted"; exit 1; }
  mkdir -p $H/scripts
  cp chain_raw.sh $P/chain_score.sh PLAN_DECOYGEN.md run_raw.sbatch run_score.sbatch $H/
  for s in build12.py build_voxels.py calib_tail.py deficit_features.py diagnose.py evaluate_v2.py features12.py make_population.py \
           make_w.py null_auc.py parent.py train_v4.py; do cp $P/scripts/$s $H/scripts/; done
  cp scripts/build_pqp12.py scripts/null_rates.py scripts/pooled_rho.py scripts/reuse_cal.py $H/scripts/
  (cd $H && sha256sum PLAN_DECOYGEN.md chain_raw.sh chain_score.sh run_raw.sbatch run_score.sbatch scripts/*.py \
     lib/library.pqp lib/pred12.parquet lib/decoy12.tsv raw01/cal/cal.json raw01/cal/anchors.tsv raw02/cal/cal.json raw02/cal/anchors.tsv > FROZEN.sha256)
  echo "arm $a: $(wc -l < $H/FROZEN.sha256) frozen files"
done
