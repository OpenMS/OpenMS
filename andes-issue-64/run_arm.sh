#!/usr/bin/env bash
# run_arm.sh <name> [extra andes flags...]  — one arm of the #64 A/B on MouseLiver-Z-T-1
set -uo pipefail
NAME="$1"; shift
B=/home/user/andes-bench
export DOTNET_ROOT=/home/user/mm/env/lib/dotnet
mkdir -p "$B/$NAME"; cd "$B/$NAME"
echo "binary: $(sha256sum /home/user/andes/target/release/andes)" > provenance.txt
echo "commit: $(git -C /home/user/andes rev-parse HEAD)" >> provenance.txt
echo "db: $(sha256sum $B/databases/mouse_entrap.fasta)" >> provenance.txt
echo "raw: 2f0142b7d79a73cb7fca84f62835af36402d24353be708aacbc0b62d66e5f071 MouseLiver-Z-T-1.raw" >> provenance.txt
echo "host: $(nproc) threads, $(free -g | awk '/Mem/{print $2}') GB, $(date -u +%FT%TZ)" >> provenance.txt
echo "flags: $*" >> provenance.txt
S=$(date +%s)
/home/user/andes/target/release/andes --spectrum "$B/spectra/MouseLiver-Z-T-1.raw" \
  --database "$B/databases/mouse_entrap.fasta" --glyco --decoy-strategy sequon-reverse \
  --threads 4 --rss-probe --output-pin "$NAME.pin" "$@" > andes.log 2>&1
echo "andes exit $? wall $(( $(date +%s) - S )) s" | tee -a provenance.txt
PIN="$NAME.glyco.pin"
if [ -s "$PIN" ]; then
  S=$(date +%s)
  /home/user/mm/env/bin/percolator --seed 42 -Y --only-psms=false \
     --results-psms "$NAME.t.psms" --decoy-results-psms "$NAME.d.psms" "$PIN" > percolator.log 2>&1
  echo "percolator exit $? wall $(( $(date +%s) - S )) s" | tee -a provenance.txt
fi
echo DONE >> provenance.txt
