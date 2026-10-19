#!/bin/bash
# usage: [REFINE=] run_andes.sh <variant: rule|norule> <dataset> [threads]
# Changes to the frozen command: no mmap candidate index (--refine needs the in-RAM index) and andes-generated
# reversed decoys on the target-only HYE FASTA (--refine refuses --decoy-strategy none)
# andes 3dded7b (rule, bigbio/andes#102) or its parent 72ee0a5 (norule), frozen #10364 andes command + --refine
set -euo pipefail
v=$1; ds=$2; th=${3:-4}; refine=${REFINE---refine}  # REFINE= (empty) for the no-refine control
D=$W/broad_benchmark/data
out=$S/runs/andes/$v${refine:+}; [ -z "$refine" ] && out=${out}_norefine; out=$out/$ds; mkdir -p $out
/usr/bin/time -v $S/ext/bin/andes-$v --spectrum $D/$ds.mzML --database $S/db/hye_targets_only.fasta --output-pin $out/native.pin \
  --output-tsv $out/native.tsv --threads $th \
  --decoy-strategy reverse --decoy-prefix DECOY_ --precursor-tol 20ppm --charge 2..5 \
  --top-n 10 --max-missed-cleavages 2 --min-length 7 --max-length 40 --max-mods 1 --mods $D/mods_none.txt \
  --precursor-cal off --isotope-error=-1..2 --protocol standard $refine > $out/andes.log 2>&1
# drop RT-prediction columns as in the frozen PIN normalization
python3 -I - $out/native.pin $out/input.pin <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
drop = {"DeltaRT", "AbsDeltaRT", "DeltaRTNorm"}
with open(src) as f, open(dst, "w") as g:
    h = f.readline().rstrip("\n").split("\t")
    keep = [i for i, n in enumerate(h) if n not in drop]
    g.write("\t".join(h[i] for i in keep) + "\n")
    for line in f:
        r = line.rstrip("\n").split("\t")
        g.write("\t".join(r[i] for i in keep[:-1]) + "\t" + "\t".join(r[keep[-1]:]) + "\n")
PY
for seed in 1 42 137; do mkdir -p $out/s$seed
  $PERC -Y -U --seed $seed --results-psms $out/s$seed/target.tsv --decoy-results-psms $out/s$seed/decoy.tsv $out/input.pin > $out/s$seed/percolator.log 2>&1
done
echo "done andes $v $ds"
