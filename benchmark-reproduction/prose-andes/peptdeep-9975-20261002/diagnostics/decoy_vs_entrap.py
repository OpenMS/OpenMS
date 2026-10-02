"""Do reversed decoys represent false targets for the PeptDeep MS2 features?

Rank-1 hits of a pd_on search, split into targets, decoys and entrapment (shuffled ENTRAP_ proteins, or for
Velos human non-UPS1 peptides absent from the sample). Within bins of the native score, compares the feature
medians of decoys and entrapment hits: both are false, so a decoy model of false targets needs them equal."""
import gzip, shutil, statistics, sys, tempfile
sys.path.insert(0, ".")
import pyopenms as p
import entrapment_velos as ev

def load(path):
    with tempfile.NamedTemporaryFile(suffix=".idXML") as tmp:
        with gzip.open(path, "rb") as f: shutil.copyfileobj(f, tmp)
        tmp.flush()
        prots, peps = [], p.PeptideIdentificationList()
        p.IdXMLFile().load(tmp.name, prots, peps)
    return peps

velos_cat = None
def category(hit, dataset):
    accs = [e.getProteinAccession() for e in hit.getPeptideEvidences()]
    if all(a.startswith("DECOY_") for a in accs): return "decoy"
    if dataset.startswith("velos"):
        global velos_cat
        if velos_cat is None:
            entries = ev.read_fasta(); db, _ = ev.classes(entries); velos_cat = db
        s = hit.getSequence().toUnmodifiedString().replace("I", "L")
        if s in velos_cat["sample"] or s in velos_cat["ups1"]: return "target"
        return "entrapment" if s in velos_cat["human"] else "target"
    t = [a for a in accs if not a.startswith("DECOY_")]
    return "entrapment" if t and all(a.startswith("ENTRAP_") for a in t) else "target"

FEATS = ["ms2_cosine", "ms2_frac_pred_found", "rt_abs_error"]
for dataset, arm in ([a.split(":") for a in sys.argv[1:]] if __name__ == "__main__" else []):
    peps = load(f"results/{dataset}/{arm}/search.idXML.gz")
    rows = []
    for pi in peps:
        hits = pi.getHits()
        if not hits: continue
        h = hits[0]
        rows.append((h.getScore(), category(h, dataset), [float(h.getMetaValue(f)) for f in FEATS]))
    rows.sort(key=lambda r: r[0])
    nb = 5
    print(f"{dataset} {arm}: rank-1 hits {len(rows)}; per native-score quintile, median of {FEATS} (n)")
    for b in range(nb):
        part = rows[b * len(rows) // nb:(b + 1) * len(rows) // nb]
        line = f"  q{b+1} score {part[0][0]:.1f}-{part[-1][0]:.1f}"
        for cat in ("target", "decoy", "entrapment"):
            sel = [r[2] for r in part if r[1] == cat]
            if len(sel) >= 10:
                med = [statistics.median(x[i] for x in sel) for i in range(len(FEATS))]
                line += f" | {cat[:5]} {med[0]:.3f} {med[1]:.3f} {med[2]:7.1f} ({len(sel)})"
        print(line)
