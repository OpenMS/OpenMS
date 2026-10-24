"""Rough wrong-peak share of the real transfers: accepted transfers into reported precursor x run pairs at q (from the
consensusXML), the pairs without a detected feature of the precursor (status no candidate / candidate missing,
candidates.py), and the false-acceptance rates per pair from pseudo_eval.py (feature absent) plus the wrong-peak share
when the feature is present. Usage: fdp_estimate.py out.consensusXML.gz cand_pairs.tsv.gz"""
import sys, gzip, xml.etree.ElementTree as ET, pandas as pd
def pl(v): v = v.strip().strip('[]'); return [x.strip() for x in v.split(',') if x.strip()]
acc = {0.01: 0, 0.03: 0, 0.05: 0}
with gzip.open(sys.argv[1], 'rb') as fh:
    for ev, el in ET.iterparse(fh, events=('end',)):
        if el.tag == 'consensusElement':
            if el.find('.//PeptideHit') is not None:
                up = {u.get('name'): u.get('value') for u in el.findall('UserParam')}
                for q in [float(x) for x in pl(up.get('mbr_transfer_qvalues', '[]'))]:
                    for t in acc:
                        if q <= t: acc[t] += 1
            el.clear()
g = pd.read_csv(sys.argv[2], sep='\t')
absent = int(g.status.isin(['no_feature', 'outside_rt']).sum())
rates = {'random': {0.01: (0.0051, 0.0074), 0.03: (0.0118, 0.0078), 0.05: (0.0130, 0.0084)},
         'probe': {0.01: (0.0, 0.0), 0.03: (0.0014, 0.0064), 0.05: (0.0063, 0.0074)}}
print('pairs without the precursor\'s feature (no candidate at the predicted position):', absent)
for q, n in acc.items():
    fa, wp = rates['random'][q]
    wrong = fa * absent + wp * n
    print(f'q {q}: accepted transfers into reported pairs {n}; expected wrong ~{wrong:.0f} ({100 * wrong / max(n, 1):.1f} %)')
