import csv, glob, os, re, sys, collections, json
csv.field_size_limit(10**9)
root = sys.argv[1]
ds_cnt = collections.Counter(); run_cnt = collections.Counter(); multi = 0; n=0
pairs = collections.Counter()
for d in sorted(glob.glob(os.path.join(root, "datasets", "*"))):
    fs = glob.glob(os.path.join(d, "*.sdrf.tsv"))
    if not fs: continue
    orgs_ds = set()
    for f in fs:
        try:
            with open(f, newline="", encoding="utf-8", errors="replace") as fh:
                r = csv.reader(fh, delimiter="\t"); hdr = [h.strip().lower() for h in next(r)]
                idx = [i for i, h in enumerate(hdr) if h.startswith("characteristics[organism]")]
                fi = next((i for i, h in enumerate(hdr) if h == "comment[data file]"), None)
                files = collections.defaultdict(set)
                for row in r:
                    for i in idx:
                        if i < len(row):
                            v = re.sub(r"\s*\(.*\)$", "", row[i].strip()).strip()
                            if v and v.lower() not in ("not available", "not applicable", "na", ""):
                                v = " ".join(v.split()[:2]).capitalize()
                                files[row[fi] if fi is not None and fi < len(row) else len(files)].add(v)
                                orgs_ds.add(v)
                for fn, os_ in files.items():
                    for o in os_: run_cnt[o] += 1
        except StopIteration:
            continue
    if orgs_ds:
        n += 1
        for o in orgs_ds: ds_cnt[o] += 1
        if len(orgs_ds) > 1:
            multi += 1
            for a in orgs_ds:
                for b in orgs_ds:
                    if a < b: pairs[(a, b)] += 1
tot_ds = sum(ds_cnt.values()); tot_run = sum(run_cnt.values())
print("datasets", n, "multi-organism", multi, "distinct organisms", len(ds_cnt))
for label, C, tot in (("datasets", ds_cnt, tot_ds), ("raw files", run_cnt, tot_run)):
    cum = 0
    print(f"\n# by {label}")
    for k, (o, c) in enumerate(C.most_common(), 1):
        cum += c
        if k <= 60 or cum / tot < 0.99:
            print(f"{k}\t{o}\t{c}\t{cum/tot:.3f}")
        if cum / tot >= 0.99: break
print("\n# top co-occurring pairs"); [print(a, b, c) for (a, b), c in pairs.most_common(15)]
json.dump({"ds": ds_cnt, "run": run_cnt}, open("orgcounts.json", "w"))
