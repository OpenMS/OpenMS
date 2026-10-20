#!/bin/bash
# Pooled Percolator (in-process, PercolatorAdapter) over all runs per database and seed, I/L-equivalent
# re-indexing, removal of PSMs with Percolator q > 0.15 (a score threshold: TDC and knockoff q-values up to 0.1 are
# unchanged, both count from the top), then FalseDiscoveryRate arms at q <= 0.01/0.05/0.10 and evaluation against
# the donor genotypes. FDR outputs are deleted after evaluation (md5 of the baseline/off pair at q 0.01 is kept).
set -e
R=${R:-/home/user/bench/results}; DB=/home/user/data/db; S=/home/user/bench/scripts
PR=/opt/openms-pr; FIX=/opt/openms-prfix; BASE=/opt/openms-base; CLS=/opt/openms-classfdr
export OPENMS_TOOL_REGISTRY_PATH=$PR/registry
run() { local root=$1; shift; LD_LIBRARY_PATH=$root/lib $root/bin/"$@"; }
SEEDS=${SEEDS:-"1 42 137"}; DBS=${DBS:-"ko4 ko0"}
mkdir -p $R/perc $R/fdr $R/eval $R/logs
arm_job() {
  local db=$1 arm=$2 s=$3 c=$4 p=$5 common="-PSM true -protein false -FDR:PSM $4"
  local o=$R/fdr/${db}_${arm}_s${s}_q$c.idXML t0=$(date +%s.%N)
        case $arm in
          off) run $PR FalseDiscoveryRate -in $p -out $o $common > $o.log 2>&1 ;;
          ko_onetier) run $PR FalseDiscoveryRate -in $p -out $o $common -PSM_variant_database $DB/var_ko4_onetier.fasta > $o.log 2>&1 ;;
          ko_tiers) run $PR FalseDiscoveryRate -in $p -out $o $common -PSM_variant_database $DB/var_ko4_tiers.fasta > $o.log 2>&1 ;;
          ko_onetier_fix) run $FIX FalseDiscoveryRate -in $p -out $o $common -PSM_variant_database $DB/var_ko4_onetier.fasta > $o.log 2>&1 ;;
          ko_tiers_fix) run $FIX FalseDiscoveryRate -in $p -out $o $common -PSM_variant_database $DB/var_ko4_tiers.fasta > $o.log 2>&1 ;;
          class) python3 -I $S/annotate_variant_class.py $p $R/fdr/tmp_cls_s${s}_q$c.idXML 2> /dev/null
                 run $CLS FalseDiscoveryRate -in $R/fdr/tmp_cls_s${s}_q$c.idXML -out $o $common -PSM_class variant_class > $o.log 2>&1
                 rm -f $R/fdr/tmp_cls_s${s}_q$c.idXML ;;
        esac
        echo "$db $arm s$s q$c FDR $(echo "$(date +%s.%N) - $t0" | bc) s" >> $R/logs/timing.txt
        if [ $db = ko4 ] && [ $arm = off ] && [ $c = 0.01 ]; then
          run $BASE FalseDiscoveryRate -in $p -out $R/fdr/${db}_base_s${s}_q$c.idXML $common > $R/fdr/${db}_base_s${s}_q$c.idXML.log 2>&1
          (cd $R/fdr && md5sum ${db}_base_s${s}_q$c.idXML ${db}_off_s${s}_q$c.idXML >> $R/logs/baseline_md5.txt)
          rm -f $R/fdr/${db}_base_s${s}_q$c.idXML
        fi
        python3 -I $S/evaluate.py $DB/var_ko4_tiers.fasta /home/user/data/pxd022787/personal.fasta $DB/human_sprot.fasta \
          ${db}_${arm}_s${s}_q$c $R/eval/${db}_${arm}_s${s}_q$c.json $o > /dev/null
        rm -f $o
}
for db in $DBS; do
  for s in $SEEDS; do
    p=$R/perc/${db}_s$s.idXML
    if [ ! -e $p ]; then
      start=$(date +%s)
      run $PR PercolatorAdapter -in $R/$db/*.idXML -out $R/perc/tmp_${db}_s$s.idXML -seed $s -score_type svm \
        -post_processing_tdc -enzyme trypsinp -subset_max_train 500000 -threads 4 > $R/logs/perc_${db}_s$s.log 2>&1
      run $PR PeptideIndexer -in $R/perc/tmp_${db}_s$s.idXML -fasta $DB/var_${db}_trim_td.fasta -out $R/perc/tmp2_${db}_s$s.idXML \
        -IL_equivalent -enzyme:name Trypsin/P -enzyme:specificity full -decoy_string DECOY_ -decoy_string_position prefix \
        -missing_decoy_action warn -unmatched_action warn -threads 4 > $R/logs/pi_${db}_s$s.log 2>&1
      rm -f $R/perc/tmp_${db}_s$s.idXML
      python3 -I $S/filter_percolator_q.py $R/perc/tmp2_${db}_s$s.idXML $p 0.15 2> $R/logs/qfilter_${db}_s$s.log
      rm -f $R/perc/tmp2_${db}_s$s.idXML
      echo "$db s$s rescoring+indexing $(( $(date +%s) - start )) s" >> $R/logs/timing.txt
    fi
    for c in 0.01 0.05 0.10; do
      if [ $db = ko4 ]; then
        arms="off ko_onetier ko_tiers ko_onetier_fix ko_tiers_fix"
      else
        arms="off class"
      fi
      for arm in $arms; do arm_job $db $arm $s $c $p & done
      wait
    done
  done
done
echo RESCORE_DONE
