"""Detail for the precursor x run pairs without a candidate (from candidates.py): how far the nearest unidentified
feature lies outside the RT window, whether an identified feature of another precursor occupies the spot, and the
species / condition of the missing run. Usage: candidates_detail.py link.tsv cand_pairs.tsv.gz out.consensusXML.gz"""
import sys, gzip, bisect, collections
import numpy as np, pandas as pd
import xml.etree.ElementTree as ET
link, pairs, cxml = sys.argv[1:4]
d = pd.read_csv(link, sep='\t')
g = pd.read_csv(pairs, sep='\t')
# identified features per run (donor rows repeat for every acceptor run; keep one)
donors = d.drop_duplicates(['donor_map', 'key'])[['donor_map', 'key', 'charge', 'mz', 'donor_rt']]
by_run = {}
for r, x in donors.groupby('donor_map'):
    x = x.sort_values('mz'); by_run[r] = (x.mz.values, x.donor_rt.values, x.charge.values, x.key.values)
def occupied(row):
    mz_, rt_, ch_, key_ = by_run[row.acceptor_map]
    tol = row.mz * 10e-6
    lo, hi = np.searchsorted(mz_, row.mz - tol), np.searchsorted(mz_, row.mz + tol, 'right')
    for i in range(lo, hi):
        if ch_[i] == row.charge and key_[i] != row.key and abs(rt_[i] - row.rt) <= row.gate:
            return True
    return False
g['other_id'] = g.apply(occupied, axis=1)
# species of reported precursors
def pl(v): v = v.strip().strip('[]'); return [x.strip() for x in v.split(',') if x.strip()]
acc = {}; species = {}
with gzip.open(cxml, 'rb') as fh:
    for ev, el in ET.iterparse(fh, events=('end',)):
        if el.tag == 'ProteinHit': acc[el.get('id')] = el.get('accession')
        elif el.tag == 'consensusElement':
            for ph in el.iter('PeptideHit'):
                sp = {a.split('_')[-1] for r in (ph.get('protein_refs') or '').split() for a in [acc.get(r, '')] if a and not a.startswith('DECOY_')}
                species[f"{ph.get('sequence')}_{el.get('charge')}"] = sp.pop() if len(sp) == 1 else 'mixed'
            el.clear()
g['species'] = g.key.map(species).fillna('?')
g['cond'] = np.where(g.acceptor_map < 3, 'A', 'B')   # map 0-2 = condition A, 3-5 = B (design order)
miss = g[g.status.isin(['outside_rt', 'no_feature'])]
print('pairs without candidate:', len(miss))
print('\noccupied by an identified feature of another precursor (same charge, 10 ppm, inside the RT window):')
print(miss.groupby('status').other_id.mean().mul(100).round(1).to_string())
o = miss[miss.status == 'outside_rt']
ratio = o.best_dt / o.gate
print('\noutside_rt: nearest unidentified feature at |dt| / half-window:')
print(pd.cut(ratio, [1, 1.5, 2, 4, 10, np.inf]).value_counts(normalize=True).sort_index().mul(100).round(1).to_string())
print('half-window (s): median %.1f, IQR %.1f-%.1f' % (o.gate.median(), o.gate.quantile(.25), o.gate.quantile(.75)))
acc_ = g[g.status == 'accepted']
print('\nfor comparison, half-window of accepted pairs: median %.1f' % acc_.gate.median())
print('\nshare without candidate by species and acceptor condition (A: yeast x2, E. coli x0.25 relative to B):')
t = g.assign(nocand=g.status.isin(['outside_rt', 'no_feature'])).groupby(['species', 'cond']).nocand.agg(['mean', 'size'])
t['mean'] = (t['mean'] * 100).round(1); print(t[t['size'] > 500].to_string())
print('\nshare without candidate by number of identified runs:')
print(g.assign(nocand=g.status.isin(['outside_rt', 'no_feature'])).groupby('n_ident').nocand.mean().mul(100).round(1).to_string())
