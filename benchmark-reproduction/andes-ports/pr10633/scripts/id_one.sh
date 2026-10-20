#!/bin/bash
# Identification of one run as in #10585: Comet -> PSMFeatureExtractor -> Percolator (PEP), quantms-style filter
# at Percolator q <= 0.10 with PEP back as main score (ProteomicsLFQ applies -psmFDR 0.01 itself).
# usage: id_one.sh <db label: hye|entrap> <target-decoy fasta> <run> [threads]
set -euo pipefail
source /home/user/bench/scripts/env.sh
use_openms /opt/openms-idtools
db=$1; fasta=$2; r=$3; t=${4:-4}
out=$BENCH/ids/$db; tmp=$BENCH/tmp/$db; mkdir -p $out $tmp
log=$out/$r.log
[ -s $out/$r.idXML ] && exit 0
start=$(date +%s)
CometAdapter -in $DATA_PB/$r.mzML -out $tmp/$r.comet.idXML -database $fasta -comet_executable /opt/thirdparty/comet.exe \
  -precursor_mass_tolerance 10 -precursor_error_units ppm -isotope_error 0/1 \
  -fragment_mass_tolerance 0.02 -fragment_bin_offset 0.0 -instrument high_res \
  -enzyme Trypsin -missed_cleavages 2 -min_peptide_length 7 -max_peptide_length 40 \
  -fixed_modifications "Carbamidomethyl (C)" -variable_modifications "Oxidation (M)" "Acetyl (Protein N-term)" \
  -max_variable_mods_in_peptide 3 -precursor_charge 2:5 -num_hits 1 \
  -PeptideIndexing:decoy_string DECOY_ -PeptideIndexing:IL_equivalent -PeptideIndexing:unmatched_action warn \
  -threads $t -no_progress > $log 2>&1
PSMFeatureExtractor -in $tmp/$r.comet.idXML -out $tmp/$r.feat.idXML -threads $t -no_progress >> $log 2>&1
rm -f $tmp/$r.comet.idXML
PercolatorAdapter -in $tmp/$r.feat.idXML -out $tmp/$r.perc.idXML -percolator_executable /opt/thirdparty/percolator \
  -score_type pep -post_processing_tdc -decoy_pattern DECOY_ -subset_max_train 300000 -seed 42 -threads $t -no_progress >> $log 2>&1
rm -f $tmp/$r.feat.idXML
IDScoreSwitcher -in $tmp/$r.perc.idXML -out $tmp/$r.q.idXML -new_score MS:1001491 -new_score_orientation lower_better \
  -new_score_type q-value -old_score PEP -no_progress >> $log 2>&1
IDFilter -in $tmp/$r.q.idXML -out $tmp/$r.qf.idXML -score:psm 0.10 -no_progress >> $log 2>&1
IDScoreSwitcher -in $tmp/$r.qf.idXML -out $out/$r.idXML -new_score PEP -new_score_orientation lower_better \
  -new_score_type "Posterior Error Probability" -old_score q-value -no_progress >> $log 2>&1
rm -f $tmp/$r.q.idXML $tmp/$r.qf.idXML $tmp/$r.perc.idXML
echo "$db $r $(( $(date +%s) - start )) s" >> $out/timing.txt
