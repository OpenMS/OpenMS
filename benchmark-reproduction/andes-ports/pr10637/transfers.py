"""Count PIP-ECHO transfers in a ProteomicsLFQ consensusXML (-out_cxml) and their entrapment share.

A consensus feature's transfers are its user params 'mbr_transfer_map_indices' / 'mbr_transfer_qvalues'.
A precursor (consensus feature) is 'entrapment' if every protein its peptide hits refer to starts with the
entrapment prefix (ENTRAP_). Transfer FDP for entrapment ratio r = 1: FDP = 2 * N_E / N (combined estimator).
Usage: transfers.py file.consensusXML [q thresholds ...]  -> one TSV row per threshold."""
import sys, re
import xml.etree.ElementTree as ET

def parse_list(v):
    v = v.strip().strip('[]')
    return [x.strip() for x in v.split(',') if x.strip()]

def analyse(path, thresholds, prefix='ENTRAP_', min_runs=3):
    acc = {}
    feats = []  # (n_identified_maps, [transfer q], is_entrapment, is_decoy)
    for ev, el in ET.iterparse(path, events=('end',)):
        tag = el.tag
        if tag == 'ProteinHit':
            acc[el.get('id')] = el.get('accession')
        elif tag == 'consensusElement':
            maps = [e.get('map') for e in el.iter('element')]
            up = {u.get('name'): u.get('value') for u in el.findall('UserParam')}
            tmaps = [int(x) for x in parse_list(up.get('mbr_transfer_map_indices', '[]'))]
            tq = [float(x) for x in parse_list(up.get('mbr_transfer_qvalues', '[]'))]
            refs = set()
            for ph in el.iter('PeptideHit'):
                refs.update((ph.get('protein_refs') or '').split())
            accs = [acc.get(r, r) for r in refs]
            entrap = bool(accs) and all(a.startswith(prefix) for a in accs)
            decoy = bool(accs) and all(a.startswith('DECOY_') for a in accs)
            n_ident = len(set(maps)) - len(set(tmaps))
            feats.append((n_ident, tq, entrap, decoy))
            el.clear()
    rows = []
    for t in thresholds:
        n = ne = prec = prec_e = prec_id = prec_id_e = 0
        for n_ident, tq, entrap, decoy in feats:
            k = sum(1 for q in tq if q <= t)
            n += k; ne += k if entrap else 0
            if n_ident + k >= min_runs:
                prec += 1; prec_e += entrap
            if n_ident >= min_runs:
                prec_id += 1; prec_id_e += entrap
        rows.append(dict(q=t, transfers=n, entrap_transfers=ne,
                         transfer_fdp=(2.0 * ne / n if n else 0.0),
                         consensus_features=len(feats),
                         prec_min3=prec, entrap_prec_min3=prec_e,
                         prec_min3_identified_only=prec_id, entrap_prec_min3_identified_only=prec_id_e,
                         entrap_features=sum(1 for f in feats if f[2])))
    return rows

if __name__ == '__main__':
    path = sys.argv[1]
    ts = [float(x) for x in sys.argv[2:]] or [0.01, 0.05]
    rows = analyse(path, ts)
    keys = list(rows[0].keys())
    print('file\t' + '\t'.join(keys))
    for r in rows:
        print(path + '\t' + '\t'.join(f'{r[k]:.4g}' if isinstance(r[k], float) else str(r[k]) for k in keys))
