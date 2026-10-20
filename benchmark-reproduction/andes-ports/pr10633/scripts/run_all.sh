#!/bin/bash
# Order of the benchmark of PR #10633 (PIP-ECHO donor_max_q) as it was run. Detection needs 10-13 GB and must run
# alone on a 15 GB machine; the combine arms need ~1.5 GB each and can run four at a time.
# Paths: /home/user/bench (BENCH in env.sh), OpenMS installs in /opt/openms-{idtools,pr,base,search}, Comet and
# Percolator in /opt/thirdparty (from OpenMS/THIRDPARTY fd8c89f, Linux/x86_64).
set -euo pipefail
S=$(dirname $0); B=/home/user/bench
mkdir -p $B/data $B/data_pb $B/db $B/run $B/results
# 1. inputs
$S/download_pb.sh $B/data_pb                       # ProteoBench mzML (as #10585) + HYE FASTA zip
curl -sS -o $B/data/ProteoBenchFASTA_MixedSpecies_HYE.fasta \
  https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/ProteoBenchFASTA_MixedSpecies_HYE.fasta
# 2. databases
#    plain HYE (ProteoBench arms, #10585 recipe): DecoyDatabase -method reverse -enzyme Trypsin  -> sha256 63410dbe...
source $S/env.sh; use_openms /opt/openms-search   # DecoyDatabase of 2f5b38d1 and of 786fc787 write the same file
DecoyDatabase -in $B/data/ProteoBenchFASTA_MixedSpecies_HYE.fasta -out $B/db/HYE_td.fasta \
  -decoy_string DECOY_ -decoy_string_position prefix -method reverse -enzyme Trypsin
#    entrapment (transfer FDP): frozen #10364 protocol: hye_td (efbf5ac2...) -> make_entrapment.py -> fb77e3df...
python3 -I $S/make_hye_td.py $B/data/ProteoBenchFASTA_MixedSpecies_HYE.fasta $B/db/hye_td_frozen.fasta
python3 -I $S/make_entrapment.py $B/db/hye_td_frozen.fasta $B/db/HYE_entrap_td.fasta
# 3. identification (Comet -> PSMFeatureExtractor -> Percolator -> q <= 0.10, PEP as score), both databases
$S/id_all.sh hye $B/db/HYE_td.fasta 4
$S/id_all.sh entrap $B/db/HYE_entrap_td.fasta 4
# 4. per database: PIP-ECHO feature detection once (PR build), then all arms from the same checkpoints
for db in hye entrap; do
  $S/detect_all.sh pr $db
  MZTAB=1 $S/lfq_arm.sh base $db base_off &      # merge-base build, same checkpoints (same revision label)
  MZTAB=1 $S/lfq_arm.sh pr $db pr_off &          # PR build, option at its default (off)
  $S/lfq_arm.sh pr $db q0.01  -PipEcho:donor_max_q 0.01 &
  $S/lfq_arm.sh pr $db q0.001 -PipEcho:donor_max_q 0.001 &
  wait
  $S/lfq_arm.sh pr $db pr_off_fdr0.01  -PipEcho:fdr 0.01 &
  $S/lfq_arm.sh pr $db q0.001_fdr0.01 -PipEcho:donor_max_q 0.001 -PipEcho:fdr 0.01 &
  [ $db = entrap ] && $S/lfq_arm.sh pr $db q0.005 -PipEcho:donor_max_q 0.005 &
  wait
  for a in $(ls $B/lfq/$db | grep -v -E "^(feat_pe|detect)$"); do $S/evaluate_arm.sh $db $a; done
done
# 5. identity check (base_off vs pr_off) and summary tables
for db in hye entrap; do for a in base_off pr_off; do for f in lfq.consensusXML lfq.mzTab msstats.csv rtlocal.msstats.csv; do
  echo -e "$db\t$a\t$f\t$($S/canon_hash.sh $B/lfq/$db/$a/$f)"; done; done; done > $B/results/identity.tsv
python3 -I $S/collect.py
