#!/bin/bash
# search the converted MS2 mzML files (convert.sh) one by one against the variant+knockoff database (ko4)
# and the variant-only database (ko0), both trimmed (trim_variant_entries.py); delete each mzML afterwards.
export LD_LIBRARY_PATH=/opt/openms-pr/lib OPENMS_TOOL_REGISTRY_PATH=/opt/openms-pr/registry
B=/opt/openms-pr/bin
D=/home/user/data/pxd022787; R=/home/user/bench/results; L=/home/user/bench/logs
source /home/user/bench/scripts/search_params.sh
ANNOT="-Search:annotate:PSM fragment_mz_error_median_ppm precursor_mz_error_ppm matched_prefix_ions_fraction matched_suffix_ions_fraction num_matched_peaks matched_prefix_ions matched_suffix_ions longest_peptide_ion_sequence matched_ion_current hyperscore_zscore ln_num_candidates matched_ion_current_fraction complementary_ions_fraction"
BATCH=${BATCH:-1}
search_batch() {
  local files=("$@") ; local tag=$(date +%H%M%S)
  for db in ko4 ko0; do
    local fasta=/home/user/data/db/var_${db}_trim_td.fasta
    local outs=(); for f in "${files[@]}"; do outs+=("$R/$db/$(basename ${f%.mzML}).idXML"); done
    eval $B/ProSE -in "${files[@]}" -database $fasta -out_idxml "${outs[@]}" -threads 4 $PROSE_PARAMS $ANNOT \
      -Search:database:chunk_size 100000 -no_progress > $L/prose_${db}_$tag.log 2>&1 || { echo "ProSE $db failed ($tag)"; return 1; }
  done
  for f in "${files[@]}"; do basename ${f%.mzML} >> $L/searched.txt; rm -f $f $f.ready; done
}
while true; do
  mz=($(for m in $D/mzml/*.mzML; do [ -e "$m.ready" ] && echo $m; done))
  fetch_done=$(grep -c FETCH_DONE /home/user/bench/fetch.log)
  pending=$(( $(ls $D/raw/*.raw 2>/dev/null | wc -l) + $(ls $D/mzml/tmp_* 2>/dev/null | wc -l) ))
  if [ ${#mz[@]} -ge $BATCH ] || { [ $fetch_done -gt 0 ] && [ $pending -eq 0 ] && [ ${#mz[@]} -gt 0 ]; }; then
    search_batch "${mz[@]:0:$BATCH}" || exit 1
  elif [ $fetch_done -gt 0 ] && [ $pending -eq 0 ]; then
    echo PROCESS_DONE; exit 0
  else
    sleep 20
  fi
done
