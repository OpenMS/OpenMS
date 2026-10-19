# andes docs/benchmarks/reproduce/run.sh (bigbio/andes fa0195a): PSMs at q <= 0.01 and entrapment hits, verbatim awk
for f in "$@"; do
  read -r n ent < <(awk -F'\t' '
    NR==1{for(i=1;i<=NF;i++) if($i=="q-value") q=i; p=NF; next}
    $q<=0.01{c++; all=1; for(i=p;i<=NF;i++) if($i !~ /ENTRAP_/){all=0; break}; e+=all}
    END{print c+0, e+0}' "$f")
  echo -e "$f\t$n\t$ent"
done
