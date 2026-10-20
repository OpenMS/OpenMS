#!/bin/bash
# Second dataset (Agilent 6550, PXD069532): entrapment database, IDs, detection, three arms.
cd /home/user/bench
# data2/: the six mzML, the dataset FASTA (with rev_ decoys) and FragPipe-LFQ-MBR.zip from
# archive.openms.de/openms/testfiles/20260717-ProteomicsLFQ-Multi-Instrument/PXD069532-Agilent-6550/
awk '/^>/{keep = ($0 !~ /^>rev_/)} keep' data2/2026-07-17-decoys-contam-GCF_000695525.1_BOL_protein.fas > data2/BOL_targets.fasta
mkdir -p db2 && python3 -I scripts/make_entrapment.py data2/BOL_targets.fasta db2/BOL_entrap_td.fasta
for r in 202505616 202505617 202505618 202505622 202505623 202505624; do
  scripts/id_one_agilent.sh entrap /home/user/bench/db2/BOL_entrap_td.fasta $r 4 || echo "FAILED id $r"
done
echo "AGILENT_IDS_DONE"
export RUNS="202505616 202505617 202505618 202505622 202505623 202505624"
for r in $RUNS; do DETECT_RUN=$r scripts/lfq_arm_agilent.sh pr entrap detect || echo "FAILED detect $r"; done
echo "AGILENT_DETECT_DONE"
(scripts/lfq_arm_agilent.sh pr entrap pr_off; echo "ARM2_DONE pr_off $?") >> logs_arms2.txt 2>&1 &
(scripts/lfq_arm_agilent.sh pr entrap q0.01 -PipEcho:donor_max_q 0.01; echo "ARM2_DONE q0.01 $?") >> logs_arms2.txt 2>&1 &
(scripts/lfq_arm_agilent.sh pr entrap q0.001 -PipEcho:donor_max_q 0.001; echo "ARM2_DONE q0.001 $?") >> logs_arms2.txt 2>&1 &
wait
echo PHASE5_DONE
