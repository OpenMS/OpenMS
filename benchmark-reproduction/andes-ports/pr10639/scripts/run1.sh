#!/bin/bash
# run1.sh <prefix:/opt/openms-X> <ini> <outdir> <threads> <db> <in...> -- [extra ProSE args]
# One ProSE search, whole process under /usr/bin/time -v; writes time.txt, search.log, *.idXML, *.pin
set -e
P=$1; INI=$2; OUT=$3; T=$4; DB=$5; shift 5
IN=(); while [ $# -gt 0 ] && [ "$1" != "--" ]; do IN+=("$1"); shift; done; [ "$1" == "--" ] && shift
mkdir -p "$OUT"; IDX=(); PIN=()
for f in "${IN[@]}"; do b=$(basename "$f"); b=${b%%.*}; IDX+=("$OUT/$b.idXML"); PIN+=("$OUT/$b.pin"); done
export LD_LIBRARY_PATH=$P/lib OPENMS_DATA_PATH=$P/share/OpenMS OMP_NUM_THREADS=$T
/usr/bin/time -v -o "$OUT/time.txt" "$P/bin/ProSE" -ini "$INI" -in "${IN[@]}" -database "$DB" \
  -out_idxml "${IDX[@]}" -out_pin "${PIN[@]}" -threads "$T" -no_progress "$@" > "$OUT/search.log" 2>&1 || { echo "FAILED $OUT"; tail -5 "$OUT/search.log"; exit 1; }
w=$(grep "Elapsed (wall" "$OUT/time.txt" | awk '{print $NF}'); r=$(grep "Maximum resident" "$OUT/time.txt" | awk '{print $NF}')
echo "$OUT wall=$w maxrss_kb=$r"
