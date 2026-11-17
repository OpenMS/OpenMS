#!/bin/bash
# Final HeLa diaPASEF benchmark (200 ng, Evosep 200SPD): OpenDIAWorkflow vs DIA-NN 2.2, same laptop, 8 threads each,
# run one after the other. Same library as the OpenSWATH/PyProphet references (openswath-hela-published):
# 200623_hela_evosep_fractions_library_ptypic_decoy.pqp (OpenDIA, its decoys) and the targets-only detecting
# transitions of the same library for DIA-NN (--gen-decoys, as in the Astral runs).
H=/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published
B=/home/singjc/.jcode/scratch/e2e/hela
D=$H/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d
cd $B
echo "== OpenDIA $(date -Is)"
/usr/bin/time -v /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin/OpenDIAWorkflow -in $D \
  -lib $H/200623_hela_evosep_fractions_library_ptypic_decoy.pqp -out_dir $B/opendia -threads 8 -seed 0 > $B/opendia.log 2>&1
echo "== DIA-NN $(date -Is)"
mkdir -p $B/tmp
/usr/bin/time -v docker run --rm --cpus 8 -u $(id -u):$(id -g) -v $H:/in:ro -v $B:/w -v /usr/bin/time:/usr/bin/time:ro \
  --entrypoint /usr/bin/time diann:2.2.0 -v /usr/diann-2.2.0/diann-linux \
  --f /in/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.d --lib /in/library_detecting_targets_only.tsv \
  --threads 8 --verbose 1 --out /w/diann_report.parquet --qvalue 1 --gen-decoys --unimod4 \
  --min-pr-charge 1 --max-pr-charge 5 --temp /w/tmp > $B/diann.log 2>&1
echo "== done $(date -Is)"
