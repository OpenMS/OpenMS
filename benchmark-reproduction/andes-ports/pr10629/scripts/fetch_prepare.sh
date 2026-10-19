# download one dataset at a time, verify sha256, sample 8000 spectra with the frozen prepare.py, verify, delete raw
set -u
W=$1; shift
for id in "$@"; do
  [ -f $W/data/$id.mzML ] && { echo "have $id"; continue; }
  read url sha sel < <(python3 -I -c "import json,sys;m=json.load(open(sys.argv[1]))[sys.argv[2]];print(m['source']['url'],m['source']['sha256'],m['selected_sha256'])" $W/../inputs.json $id)
  ext=mzML; [[ $url == *.gz ]] && ext=mzML.gz
  for i in 1 2 3 4; do curl -sSfL -o $W/raw/$id.$ext "$url" && break; sleep $((2**i)); done
  echo "$sha  $W/raw/$id.$ext" | sha256sum -c - || { echo "RAW MISMATCH $id"; continue; }
  (cd $W/.. && python3 -I $W/prepare.py $id)
  echo "$sel  $W/data/$id.mzML" | sha256sum -c - || echo "SELECTED MISMATCH $id"
  rm -f $W/raw/$id.$ext
done
echo FETCH_DONE
