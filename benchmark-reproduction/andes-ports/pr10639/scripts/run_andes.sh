#!/bin/bash
# run_andes.sh <outdir> <slice_da> <db> <mzML> <mods> <max_mods> [threads]  - the frozen #10364 andes command (minus --candidate-index mmap and --mmap-window-cache-candidates, removed in andes HEAD, where the out-of-core path is the only one), with the given slice width
set -e
OUT=$1; SL=$2; DB=$3; MZ=$4; MODS=$5; MM=$6; T=${7:-4}; mkdir -p "$OUT"
/usr/bin/time -v -o "$OUT/time.txt" /home/user/bench/andes/target/release/andes --spectrum "$MZ" --database "$DB" \
  --output-pin "$OUT/native.pin" --output-tsv "$OUT/native.tsv" --threads $T \
  --fragment-index-slice-da $SL --decoy-strategy none --decoy-prefix DECOY \
  --precursor-tol 20ppm --charge 2..5 --top-n 10 --max-missed-cleavages 2 --min-length 7 --max-length 40 \
  --max-mods $MM --mods "$MODS" --precursor-cal off --isotope-error=-1..2 --protocol standard > "$OUT/search.log" 2>&1 || { echo "FAILED $OUT"; tail -3 "$OUT/search.log"; }
echo "$OUT wall=$(grep 'Elapsed (wall' $OUT/time.txt | awk '{print $NF}') maxrss_kb=$(grep 'Maximum resident' $OUT/time.txt | awk '{print $NF}')"
