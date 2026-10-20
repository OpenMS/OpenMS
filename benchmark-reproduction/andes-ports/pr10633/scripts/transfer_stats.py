"""Count PIP-ECHO transfers in a ProteomicsLFQ consensusXML and how many of them
transfer entrapment-only precursors.

A consensus feature is one precursor (sequence, charge). PIP-ECHO records its
accepted transfers in the meta values 'mbr_transfer_map_indices' and
'mbr_transfer_qvalues'. The precursor's class comes from the protein references
of its peptide hits: 'entrapment' if all referenced proteins are entrapment
proteins (accession prefix ENTRAP_, decoys ignored), 'target' if any is a target
protein (EntrapmentEvaluation's rule), 'decoy' if all are decoys.

Combined FDP estimate (Wen et al. 2025): FDP = N_E (1 + T/E) / (N_T + N_E),
T/E = 1 for a 1:1 entrapment database.
"""
import argparse, json, sys
import xml.etree.ElementTree as ET

ap = argparse.ArgumentParser()
ap.add_argument("cxml")
ap.add_argument("--label", default="")
ap.add_argument("--entrap", default="ENTRAP_")
ap.add_argument("--decoy", default="DECOY_")
ap.add_argument("--te", type=float, default=1.0, help="T/E ratio")
ap.add_argument("--qs", default="0.01,0.05")
ap.add_argument("--header", action="store_true")
ap.add_argument("--json", action="store_true")
ap.add_argument("--fasta", help="searched FASTA: an 'entrapment' peptide that also occurs in a target protein with I = L "
                "counts as target (as EntrapmentEvaluation does)")
a = ap.parse_args()
qs = [float(x) for x in a.qs.split(",")]

ph_acc = {}          # PH_n -> accession
n_maps = 0
rows = []            # (class, n_ident_runs, transfer_qs)

def parse_list(v):
    v = v.strip().strip("[]")
    return [x.strip() for x in v.split(",") if x.strip()]

def classify(accs):
    accs = [x for x in accs]
    if not accs:
        return "unknown"
    non_decoy = [x for x in accs if not x.startswith(a.decoy)]
    if not non_decoy:
        return "decoy"
    if all(x.startswith(a.entrap) for x in non_decoy):
        return "entrapment"
    return "target"

target_text = None
if a.fasta:
    seqs, cur, keep = [], [], False
    for line in open(a.fasta):
        if line.startswith(">"):
            if keep: seqs.append("".join(cur))
            ident = line[1:].split()[0] if len(line) > 1 else ""
            keep = not ident.startswith(a.decoy) and not ident.startswith(a.entrap)
            cur = []
        else:
            cur.append(line.strip())
    if keep: seqs.append("".join(cur))
    target_text = "#".join(seqs).replace("I", "L")

def unmodified(seq):
    import re
    return re.sub(r"\(.*?\)|\[.*?\]", "", seq).lstrip(".").rstrip(".")

for ev, el in ET.iterparse(a.cxml, events=("end",)):
    tag = el.tag
    if tag == "ProteinHit":
        ph_acc[el.get("id")] = el.get("accession")
    elif tag == "map":
        n_maps += 1
    elif tag == "consensusElement":
        maps = [int(e.get("map")) for e in el.iter("element")]
        accs, seqs = set(), set()
        for ph in el.iter("PeptideHit"):
            seqs.add((ph.get("sequence"), ph.get("charge")))
            for r in (ph.get("protein_refs") or "").split():
                accs.add(ph_acc.get(r, "?"))
        t_idx, t_q = [], []
        for up in el.findall("UserParam"):
            if up.get("name") == "mbr_transfer_map_indices":
                t_idx = [int(x) for x in parse_list(up.get("value"))]
            elif up.get("name") == "mbr_transfer_qvalues":
                t_q = [float(x) for x in parse_list(up.get("value"))]
        n_ident = len(set(maps) - set(t_idx))
        cls = classify(accs)
        if cls == "entrapment" and target_text is not None and \
                any(unmodified(sq).replace("I", "L") in target_text for sq, _ in seqs):
            cls = "target"
        rows.append((cls, n_ident, t_q, len(set(maps)), len(seqs)))
        el.clear()
    elif tag in ("UnassignedPeptideIdentification",):
        el.clear()

res = {"label": a.label, "maps": n_maps, "consensus_features": len(rows)}
for c in ("target", "entrapment", "decoy", "unknown"):
    res[f"precursors_{c}"] = sum(1 for r in rows if r[0] == c)
for q in qs:
    tag = f"q{q:g}"
    N = {c: 0 for c in ("target", "entrapment", "decoy", "unknown")}
    P = {c: 0 for c in N}
    for cls, n_ident, tq, nm, ns in rows:
        k = sum(1 for x in tq if x <= q)
        N[cls] += k
        P[cls] += 1 if k > 0 else 0
    nt, ne = N["target"], N["entrapment"]
    res[f"transfers_{tag}"] = sum(N.values())
    res[f"transfers_target_{tag}"] = nt
    res[f"transfers_entrap_{tag}"] = ne
    res[f"transfers_decoy_{tag}"] = N["decoy"]
    res[f"transfers_unknown_{tag}"] = N["unknown"]
    res[f"transfer_FDP_{tag}"] = (ne * (1 + a.te) / (nt + ne)) if nt + ne else 0.0
    res[f"precursors_with_transfer_{tag}"] = sum(P.values())
    res[f"precursors_with_transfer_entrap_{tag}"] = P["entrapment"]
# identified (donor) side
id_t = sum(1 for r in rows if r[1] > 0 and r[0] == "target")
id_e = sum(1 for r in rows if r[1] > 0 and r[0] == "entrapment")
res["identified_precursors_target"] = id_t
res["identified_precursors_entrap"] = id_e
res["identified_precursor_FDP"] = id_e * (1 + a.te) / (id_t + id_e) if id_t + id_e else 0.0
# quantified in >= 3 runs (identified or transferred at the default PIP-ECHO cut)
q3 = [r for r in rows if r[3] >= 3 and r[0] in ("target", "entrapment")]
res["precursors_min3"] = len(q3)
res["precursors_min3_entrap"] = sum(1 for r in q3 if r[0] == "entrapment")
res["precursors_min3_entrap_share"] = res["precursors_min3_entrap"] / len(q3) if q3 else 0.0

if a.json:
    print(json.dumps(res))
else:
    if a.header:
        print("\t".join(res.keys()))
    print("\t".join(f"{v:.5f}" if isinstance(v, float) else str(v) for v in res.values()))
