"""Properties of a target + entrapment (+ decoy) database: are the entrapment proteins/peptides mass-matched to their
targets, how many entrapment peptides collide with target peptides, and T/E.

  db_stats.py <db.fasta> [entrapment_prefix=ENTRAP_] [decoy_prefix=DECOY_]

Pairs an entrapment protein with its target by accession (ENTRAP_<acc> or ENTRAP_<acc>_<n>). Peptides: fully
tryptic, up to 2 missed cleavages, length 7-40, 'Trypsin' (no cleavage before P, as DecoyDatabase's default enzyme
and EntrapmentEvaluator's T/E) and 'Trypsin/P' (the search enzyme of the benchmark). Peptides are paired by position
in proteins whose cleavage sites stayed in place; 'entrapment_without_target_of_same_composition' counts
entrapment peptides for which no target peptide has the same amino acid composition (so none has the same mass)."""
import collections, sys

path = sys.argv[1]; EP = sys.argv[2] if len(sys.argv) > 2 else "ENTRAP_"; DP = sys.argv[3] if len(sys.argv) > 3 else "DECOY_"
entries, hdr, seq = [], None, []
for line in open(path):
    if line.startswith(">"):
        if hdr is not None: entries.append((hdr, "".join(seq)))
        hdr, seq = line[1:].split()[0] if line[1:].strip() else "", []
    else: seq.append(line.strip())
entries.append((hdr, "".join(seq)))
targets = {h: s for h, s in entries if not h.startswith(DP) and not h.startswith(EP)}
entraps = [(h, s) for h, s in entries if h.startswith(EP)]


def sites(s, p_rule):
    return [0] + [i + 1 for i in range(len(s) - 1) if s[i] in "KR" and (not p_rule or s[i + 1] != "P")] + [len(s)]


def digest_pos(s, p_rule):
    st = sites(s, p_rule); out = []
    for i in range(len(st) - 1):
        for j in range(i + 1, min(i + 4, len(st))):
            if 7 <= st[j] - st[i] <= 40: out.append((st[i], st[j]))
    return out


res = collections.OrderedDict(database=path, targets=len(targets), entrapment=len(entraps),
                              decoys=sum(h.startswith(DP) for h, _ in entries))
unpaired = same_len = same_comp = same_seq = 0
pairs = []
for h, s in entraps:
    acc = h[len(EP):]
    t = targets.get(acc) or targets.get(acc.rsplit("_", 1)[0])
    if t is None: unpaired += 1; continue
    pairs.append((t, s))
    same_len += len(t) == len(s); same_comp += sorted(t) == sorted(s); same_seq += t == s
res.update(unpaired=unpaired, protein_same_length=same_len, protein_same_composition=same_comp, protein_identical=same_seq)
for p_rule, name in ((True, "Trypsin"), (False, "Trypsin/P")):
    T = set(); TL = set()
    for s in targets.values():
        for a, b in digest_pos(s, p_rule): T.add(s[a:b])
    TL = {p.replace("I", "L") for p in T}
    E = set(); mut_pep = set(); site_mismatch = 0
    for h, s in entraps:
        for a, b in digest_pos(s, p_rule): E.add(s[a:b])
    for t, s in pairs:
        if sites(t, p_rule) != sites(s, p_rule):
            site_mismatch += 1; continue
        for a, b in digest_pos(s, p_rule):
            if sorted(s[a:b]) != sorted(t[a:b]): mut_pep.add(s[a:b])
    # composition (= mass) of every target peptide: an entrapment peptide without a target peptide of the same
    # composition has no mass-matched counterpart
    TC = {"".join(sorted(p)) for p in T}
    no_comp = sum("".join(sorted(p)) not in TC for p in E)
    E_only = E - T
    E_only_IL = {p for p in E if p.replace("I", "L") not in TL}
    res[name] = dict(target_peptides=len(T), entrapment_peptides=len(E), entrapment_equal_to_a_target=len(E) - len(E_only),
                     entrapment_equal_to_a_target_IL=len(E) - len(E_only_IL),
                     proteins_with_moved_cleavage_sites=site_mismatch,
                     paired_peptides_composition_changed=len(mut_pep),
                     entrapment_without_target_of_same_composition=no_comp, entrapment_without_target_of_same_composition_pct=round(100 * no_comp / max(len(E), 1), 3),
                     TE_exact=round(len(T) / max(len(E_only), 1), 6), TE_IL=round(len(TL) / max(len(E_only_IL), 1), 6))
import json; print(json.dumps(res, indent=1))
