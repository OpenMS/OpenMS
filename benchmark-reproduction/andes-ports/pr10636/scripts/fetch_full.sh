#!/bin/sh
# Whole runs for the isofull setting; the .gz sha256 is checked against inputs.json (source of the frozen subsets).
P=/home/user/bench/proto/z/prose-andes-reproduction; O=/home/user/bench/dl/full
while read id url; do
  [ -f $O/$id.mzML ] && continue
  curl -sSfL --retry 3 "$url" -o $O/$id.mzML.gz
  want=$(python3 -c "import json;print(json.load(open('$P/inputs.json'))['$id']['source']['sha256'])")
  got=$(sha256sum $O/$id.mzML.gz | cut -c1-64)
  [ "$want" = "$got" ] || { echo "SHA MISMATCH $id"; continue; }
  gunzip $O/$id.mzML.gz && echo "OK $id $got"
done < "$(dirname "$0")/full_urls.txt"
echo FETCH_DONE
