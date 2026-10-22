"""Why do precursors identified in few runs get no PIP-ECHO candidate transfer in the others?

Inputs: the PIPECHO_DIAG dump of link_donors_and_acceptors (one row per donor feature and acceptor run) and the
consensusXML of the same ProteomicsLFQ run (accepted transfers, PIP-ECHO fdr 0.05).
For every precursor (sequence + charge, PIP-ECHO's key) and every run in which it has no donor (identified feature):
  accepted       a transfer was accepted into that run
  candidate      at least one donor found a candidate acceptor, but no transfer was accepted (FDR)
  grid_miss      a compatible unidentified feature lies inside the window but the grid search did not return it
  outside_rt     unidentified features of that charge within the m/z tolerance exist, but none inside the RT window
  no_feature     no unidentified feature of that charge within the m/z tolerance anywhere in the run
Usage: candidates.py link.tsv out.consensusXML[.gz] out_prefix"""
import sys, gzip, collections
import xml.etree.ElementTree as ET
import pandas as pd

link_tsv, cxml, out = sys.argv[1:4]
d = pd.read_csv(link_tsv, sep='\t')
d = d[d.widen == d.widen.max()]          # last widening pass is the one used
runs = sorted(set(d.donor_map) | set(d.acceptor_map))
ident = d.groupby('key').donor_map.apply(lambda s: set(s)).to_dict()

def parse_list(v):
    v = v.strip().strip('[]'); return [x.strip() for x in v.split(',') if x.strip()]
accepted = collections.defaultdict(set)
reported = set()   # precursors ProteomicsLFQ reports (consensus features that keep peptide hits)
op = gzip.open if cxml.endswith('.gz') else open
with op(cxml, 'rb') as fh:
    for ev, el in ET.iterparse(fh, events=('end',)):
        if el.tag == 'consensusElement':
            up = {u.get('name'): u.get('value') for u in el.findall('UserParam')}
            idx = [int(x) for x in parse_list(up.get('mbr_transfer_map_indices', '[]'))]
            keys = {f"{ph.get('sequence')}_{el.get('charge')}" for ph in el.iter('PeptideHit')}
            reported.update(keys)
            for key in keys:
                accepted[key].update(idx)
            el.clear()

g = d.groupby(['key', 'acceptor_map']).agg(found=('found', 'max'), n_mz=('n_mz', 'max'), n_in=('n_in_window', 'max'),
                                            best_dt=('best_dt', lambda s: s.abs().min()), gate=('gate', 'median'),
                                            charge=('charge', 'first'), mz=('mz', 'first'), rt=('center', 'median')).reset_index()
def status(r):
    if r.acceptor_map in accepted.get(r.key, ()): return 'accepted'
    if r.found: return 'candidate'
    if r.n_in > 0: return 'grid_miss'
    if r.n_mz > 0: return 'outside_rt'
    return 'no_feature'
g = g[g.key.isin(reported)]
print('reported precursors with a donor:', g.key.nunique())
g['n_ident'] = g.key.map(lambda k: len(ident[k]))
g = g[[a not in ident[k] for k, a in zip(g.key, g.acceptor_map)]]
g['status'] = g.apply(status, axis=1)
g.to_csv(out + '_pairs.tsv.gz', sep='\t', index=False)
tab = pd.crosstab(g.n_ident, g.status, margins=True)
cols = [c for c in ['accepted', 'candidate', 'grid_miss', 'outside_rt', 'no_feature', 'All'] if c in tab.columns]
tab = tab[cols]
print('precursor x run pairs without an identification, by number of runs the precursor is identified in')
print(tab.to_string())
print((tab[cols[:-1]].div(tab['All'], axis=0) * 100).round(1).to_string())
tab.to_csv(out + '_status.tsv', sep='\t')
# precursors identified in 1-2 runs: how many could reach >= 3 runs
p = g[g.n_ident <= 2].groupby('key').agg(n_ident=('n_ident', 'first'),
        acc=('status', lambda s: (s == 'accepted').sum()), cand=('status', lambda s: s.isin(['accepted', 'candidate']).sum()))
p['reach_accepted'] = p.n_ident + p.acc >= 3
p['reach_if_all_candidates'] = p.n_ident + p.cand >= 3
p['reach_if_all_runs'] = True
print('\nprecursors identified in 1-2 runs:', len(p))
print(' reach >= 3 runs with accepted transfers:', int(p.reach_accepted.sum()))
print(' would reach >= 3 if every candidate were accepted:', int(p.reach_if_all_candidates.sum()))
