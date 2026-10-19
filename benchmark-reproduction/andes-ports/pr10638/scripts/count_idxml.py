"""PSMs in a ProSE idXML (in-process Percolator, FDR:PSM 0.01): target PSMs, split by chimeric_rank, unique peptides,
entrapment PSMs (all proteins ENTRAP_) and the FDP 2*N_E/N.  usage: count_idxml.py <idXML>"""
import sys, json, xml.etree.ElementTree as ET
acc = {}; out = {}
for ev, el in ET.iterparse(sys.argv[1], events=("end",)):
    if el.tag == "ProteinHit":
        acc[el.get("id")] = el.get("accession")
    elif el.tag == "PeptideIdentification":
        for h in el.findall("PeptideHit"):
            uv = {u.get("name"): u.get("value") for u in h.findall("UserParam")}
            if uv.get("target_decoy") == "decoy": continue
            rank = uv.get("chimeric_rank", "1")
            prots = [acc.get(r, r) for r in (h.get("protein_refs") or "").split()]
            ent = bool(prots) and all(p.startswith("ENTRAP_") for p in prots)
            for k in ("all", "r" + rank):
                d = out.setdefault(k, {"psms": 0, "entrap": 0, "peps": set()})
                d["psms"] += 1; d["entrap"] += ent; d["peps"].add(h.get("sequence"))
        el.clear()
print(json.dumps({k: {"psms": v["psms"], "peptides": len(v["peps"]), "entrap": v["entrap"],
                      "fdp": round(2 * v["entrap"] / v["psms"], 5)} for k, v in out.items()}))
