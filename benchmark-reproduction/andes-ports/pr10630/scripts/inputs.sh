set -e
F=/home/user/bench/frozen; B=$F/broad_benchmark; mkdir -p $B/data $B/raw
cd $B/data
[ -f hye_targets.fasta ] || curl -sSfL -o hye_targets.fasta https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/ProteoBenchFASTA_MixedSpecies_HYE.fasta
sha256sum hye_targets.fasta
python3 /home/user/bench/scripts/make_dbs.py $B
for d in $(python3 -c "import json;print(' '.join(json.load(open('$F/inputs.json'))))"); do
  [ -f $B/data/$d.mzML ] && { echo have $d; continue; }
  read url sha <<<$(python3 -c "import json;s=json.load(open('$F/inputs.json'))['$d']['source'];print(s['url'],s['sha256'])")
  ext=mzML; [[ $url == *.gz ]] && ext=mzML.gz
  curl -sSfL --retry 3 -o $B/raw/$d.$ext "$url"
  echo "$sha  $B/raw/$d.$ext" | sha256sum -c -
  python3 $B/prepare.py $d
  sel=$(python3 -c "import json;print(json.load(open('$F/inputs.json'))['$d']['selected_sha256'])")
  echo "$sel  $B/data/$d.mzML" | sha256sum -c -
  rm -f $B/raw/$d.$ext
done
echo INPUTS_DONE
