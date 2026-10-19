"""Best hit per spectrum of a FalseDiscoveryRate output idXML -> TSV
(spectrum, q, decoy, entrapment, modification class, bench_mc, bench_z, bench_iso, bench_term, bench_len, class_q).

entrapment: every protein of the hit is ENTRAP_* (as EntrapmentEvaluator and evaluate_entrap_d.py count it).
modification class: as ClassFDR::classFromModifications - sorted names of non-fixed modifications, or "unmodified".
usage: extract.py in.idXML out.tsv fixed_mod_regex
"""
import re, sys, gzip
from lxml import etree as ET

fixed = re.compile(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] else None
out = gzip.open(sys.argv[2], "wt") if sys.argv[2].endswith(".gz") else open(sys.argv[2], "w")
out.write("spectrum\tq\tdecoy\tentrap\tmodclass\tmc\tz\tiso\tclass_q\tterm\tlen\n")
acc = {}
ctx = ET.iterparse(sys.argv[1], events=("end",), huge_tree=True)
for ev, el in ctx:
    tag = el.tag.rsplit("}", 1)[-1]
    if tag == "ProteinHit":
        acc[el.get("id")] = el.get("accession")
    elif tag == "PeptideIdentification":
        hits = [h for h in el if h.tag.rsplit("}", 1)[-1] == "PeptideHit"]
        if hits:
            hb = el.get("higher_score_better") == "true"
            h = max(hits, key=lambda x: float(x.get("score")) * (1 if hb else -1))
            up = {u.get("name"): u.get("value") for u in h if u.tag.rsplit("}", 1)[-1] == "UserParam"}
            accs = [acc[r] for r in (h.get("protein_refs") or "").split()]
            decoy = up.get("target_decoy") == "decoy"
            entrap = (not decoy) and bool(accs) and all(a.startswith("ENTRAP_") for a in accs)
            mods = sorted(set(m for m in re.findall(r"\(([^)]*)\)", h.get("sequence"))
                              if not (fixed and fixed.fullmatch(m))))
            mc = "+".join(mods) if mods else "unmodified"
            out.write("\t".join([el.get("spectrum_reference", ""), h.get("score"), str(int(decoy)), str(int(entrap)), mc,
                                 up.get("bench_mc", ""), up.get("bench_z", ""), up.get("bench_iso", ""),
                                 up.get("class_q-value", ""), up.get("bench_term", ""), up.get("bench_len", "")]) + "\n")
        el.clear()
out.close()
