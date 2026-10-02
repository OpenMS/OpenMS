#!/bin/sh
# probe2.sh <dataset> <binary> <name> <json-overrides>: search a dataset's no-candidate subset with its frozen parameters plus overrides
D=$1; R=../../../prose-andes-reproduction
export OMP_NUM_THREADS=2 OPENMS_DATA_PATH=$R/nightly/pyopenms/share/OpenMS
DB=$(python3 -c "import json;print([r['database'] for r in json.load(open('$R/broad_benchmark/manifest.json')) if r['id']=='$D'][0])")
python3 -c "import json,sys;p=json.load(open('../../results/$D/base/params.json'));p.update(json.loads(sys.argv[1]));json.dump(p,open('$D.$3.params.json','w'))" "$4"
$2 ${D}_nocand.mzML $R/broad_benchmark/data/${DB}_td.fasta $D.$3.params.json $D.$3.tsv $D.$3.idXML > $D.$3.log 2>&1
python3 - "$D" "$3" <<'PY'
import csv, json, re, sys, collections
d, name = sys.argv[1:]
want = json.load(open(f"{d}_nocand.json"))
seq = lambda s: re.sub('[^A-Z]', '', re.sub(r'\([^)]*\)|\[[^]]*\]', '', s)).replace('L', 'I')
hits = collections.defaultdict(list)
with open(f"{d}.{name}.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        hits[r["scan"].split("scan=")[-1]].append((float(r["score"]), seq(r["peptide"])))
top1 = sum(1 for s, q in want.items() if hits.get(s) and max(hits[s])[1] == q)
anyc = sum(1 for s, q in want.items() if q in [h[1] for h in hits.get(s, [])])
print(f"{d:16s} {name:10s} with candidates {sum(1 for s in want if hits.get(s)):4d}/{len(want)}  Sage peptide among candidates {anyc:4d}  rank 1 {top1:4d}")
PY
