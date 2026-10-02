"""Accepted entrapment PSMs of a shuffled-entrapment arm: is the spectrum's other candidate the entrapment
peptide's own unshuffled target (same composition, same mass)? Such 'sibling' confusions share most fragments
and, by construction, the predicted retention time."""
import csv, gzip, shutil, sys, tempfile, collections
import pyopenms as p

def load(path):
    with tempfile.NamedTemporaryFile(suffix=".idXML") as tmp:
        with gzip.open(path, "rb") as f: shutil.copyfileobj(f, tmp)
        tmp.flush()
        prots, peps = [], p.PeptideIdentificationList()
        p.IdXMLFile().load(tmp.name, prots, peps)
    return peps

for spec in sys.argv[1:]:
    dataset, arm = spec.split(":")
    peps = load(f"results/{dataset}/{arm}/search.idXML.gz")
    by_scan = {}
    for pi in peps:
        ref = pi.getMetaValue("spectrum_reference") if pi.metaValueExists("spectrum_reference") else ""
        scan = str(ref).split("scan=")[-1]
        by_scan[scan] = pi
    stats = collections.Counter()
    for seed in (1, 42, 137):
        for row in csv.DictReader(open(f"results/{dataset}/{arm}/s{seed}/target.tsv"), delimiter="\t"):
            if float(row["q-value"]) > 0.01: continue
            prots = row["proteinIds"].split("\t") if "proteinIds" in row else []
            stats["accepted"] += 1
        with open(f"results/{dataset}/{arm}/s{seed}/target.tsv") as f:
            next(f)
            for line in f:
                x = line.rstrip("\n").split("\t")
                if float(x[2]) > 0.01: continue
                accs = [a for a in x[5:] if a]
                if not all(a.startswith("ENTRAP_") for a in accs): continue
                stats["entrapment"] += 1
                scan = x[0].rsplit("_", 2)[-2]
                pi = by_scan.get(scan)
                if pi is None: stats["no spectrum"] += 1; continue
                pep = x[4].split(".", 1)[1].rsplit(".", 1)[0] if "." in x[4] else x[4]
                ent = p.AASequence.fromString(pep).toUnmodifiedString()
                comp = sorted(ent)
                sib = None
                for r, h in enumerate(pi.getHits()):
                    s = h.getSequence().toUnmodifiedString()
                    if s != ent and sorted(s) == comp and s[-1] == ent[-1]:
                        sib = (r, h); break
                if sib is None: stats["no sibling among candidates"] += 1
                else:
                    stats[f"sibling target at rank {sib[0] + 1}"] += 1
    print(dataset, arm, dict(stats))
