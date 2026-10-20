#!/bin/bash
# Collect all arms into results/*.tsv: ProteoBench metrics (raw = ProteomicsLFQ out_msstats.csv, rt_local = normalized)
# and transfer / entrapment counts. Usage: collect.sh <runs_dir> <out_dir>
R=$1; O=$2; S=$(dirname "$0"); mkdir -p $O
PY=${PB_PYTHON:-/home/user/pbenv/bin/python}
echo -e "db\tarm\tnormalization\tprecursors\tmedian_abs_eps\tecoli_abs_eps\tcv" > $O/proteobench.tsv
for d in $R/hye/*/; do a=$(basename $d); [ -s $d/norm.csv ] || continue
  for f in out_msstats:none norm:rt_local; do
    $PY -I $S/score_proteobench.py $d/${f%%:*}.csv 2>/dev/null | tail -1 | awk -v a=$a -v n=${f#*:} -F'\t' '{print "hye\t" a "\t" n "\t" $2 "\t" $3 "\t" $4 "\t" $5}' >> $O/proteobench.tsv
  done
done
head -1 $R/hye/off_q05/transfers.tsv | sed 's/^file/db\tarm/' > $O/transfers.tsv
for db in hye entrap; do for d in $R/$db/*/; do a=$(basename $d); [ -s $d/transfers.tsv ] || continue
  tail -n +2 $d/transfers.tsv | cut -f2- | sed "s/^/$db\t$a\t/" >> $O/transfers.tsv; done; done
for db in hye entrap; do for d in $R/$db/*/; do [ -s $d/sha256.txt ] && sed "s|  |  $db/$(basename $d)/|" $d/sha256.txt; done; done > $O/output_sha256.txt
