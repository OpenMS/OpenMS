"""Per-spectrum comparison of Sage v0.14.7 with ProSE arms on the 20-file suite.

A spectrum counts as accepted by an engine when Percolator accepts it (q <= 0.01) in at least two of the three seeds;
its peptide is the stripped sequence (I = L, modifications removed) of that majority. For every spectrum accepted by
one engine only, the other engine's candidate list (top 10 of its PIN) gives the reason: no_candidate (scan absent from
the PIN), absent_top10, rank1_not_accepted (the same sequence is the native rank-one candidate) or below_rank1.
Spectra both engines accept with different sequences are 'conflict'."""
import collections, csv, gzip, json, statistics, sys
from pathlib import Path
B = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(B / "prose-andes-reproduction/broad_benchmark"))
from analyze import sequence
SAGE = B / "sage-benchmark/results"
PROSE = B / "pr-untangled/results"
MANIFEST = json.loads((B / "prose-andes-reproduction/broad_benchmark/manifest.json").read_text())
sys.path.insert(0, str(B / "pr-untangled"))
from evaluate import group_of


def opener(p):
    p = Path(p)
    return p.open() if p.exists() else gzip.open(str(p) + ".gz", "rt")


def accepted(folder):
    votes = collections.defaultdict(list)
    for s in (1, 42, 137):
        with (folder / f"s{s}" / "target.tsv").open() as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if float(r["q-value"]) <= 0.01:
                    votes[int(r["PSMId"].rsplit("_", 2)[-2])].append(sequence(r["peptide"]))
    out = {}
    for scan, seqs in votes.items():
        seq, n = collections.Counter(seqs).most_common(1)[0]
        if n >= 2: out[scan] = seq
    return out


def candidates(pin, score_col):
    by_scan = collections.defaultdict(list)
    with opener(pin) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            by_scan[int(r["ScanNr"])].append((-float(r[score_col]), r["Label"] == "1", sequence(r["Peptide"])))
    return {k: [(t, s) for _, t, s in sorted(v)] for k, v in by_scan.items()}


def classify(scan, seq, other_acc, other_cands):
    if other_acc.get(scan) == seq: return "same"
    if scan in other_acc: return "conflict"
    c = other_cands.get(scan)
    if not c: return "no_candidate"
    seqs = [s for _, s in c]
    if seq not in seqs: return "absent_top10"
    return "rank1_not_accepted" if seqs[0] == seq else "below_rank1"


def main(ARMS):
    rows = []
    for r in MANIFEST:
        d = r["id"]
        if not (SAGE / d / "summary.json").exists(): continue
        sage_acc = accepted(SAGE / d)
        sage_c = candidates(SAGE / d / "input.pin", "ln(hyperscore)")
        for arm in ARMS:
            folder = PROSE / d / arm
            if not (folder / "summary.json").exists(): continue
            prose_acc = accepted(folder)
            prose_c = candidates(folder / "native.pin", "score")
            sage_only = collections.Counter(classify(s, q, prose_acc, prose_c) for s, q in sage_acc.items() if s not in prose_acc or prose_acc[s] != q)
            prose_only = collections.Counter(classify(s, q, sage_acc, sage_c) for s, q in prose_acc.items() if s not in sage_acc or sage_acc[s] != q)
            both = sum(1 for s, q in sage_acc.items() if prose_acc.get(s) == q)
            rows.append(dict(dataset=d, group=group_of(d), arm=arm, sage=len(sage_acc), prose=len(prose_acc), both_same=both,
                             conflict=sage_only["conflict"],
                             **{f"sage_only_{k}": sage_only[k] for k in ["no_candidate", "absent_top10", "below_rank1", "rank1_not_accepted"]},
                             **{f"prose_only_{k}": prose_only[k] for k in ["no_candidate", "absent_top10", "below_rank1", "rank1_not_accepted"]}))
    out = B / "pr-untangled/sage-compare/per_spectrum.tsv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    keys = [k for k in rows[0] if k not in ("dataset", "group", "arm")]
    for arm in ARMS:
        print(f"\n== Sage vs ProSE {arm} (majority-accepted spectra, summed over files)")
        print("group".ljust(22) + "".join(k.replace("sage_only_", "S:").replace("prose_only_", "P:")[:14].rjust(15) for k in keys))
        for g in dict.fromkeys(x["group"] for x in rows):
            sub = [x for x in rows if x["arm"] == arm and x["group"] == g]
            if sub: print(g.ljust(22) + "".join(str(sum(x[k] for x in sub)).rjust(15) for k in keys))


if __name__ == "__main__":
    main(sys.argv[1:] or ["base", "pr10335"])
