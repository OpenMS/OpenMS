"""PLAN_V11 descriptive diagnostic (labels only to count): accepted-set change at q <= 0.01, D11FUWT vs D11FWT, with the
unique-fragment evidence (natU) and the all-fragment evidence (natE) of gained / lost precursors."""
import sys, numpy as np, pandas as pd
sys.argv = ['x', 'D11FWT']; sys.path.insert(0, '.')
import importlib.util
spec = importlib.util.spec_from_file_location('ev', 'evaluate_v2.py')
src = open('evaluate_v2.py').read().split('rng = np.random.default_rng')[0]; ns = {}; exec(src, ns)
tdc, ent = ns['tdc'], ns['ent']
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{raw}/population.tsv', sep='\t'); XU = np.load(f'{raw}/natU.npy'); XE = np.load(f'{raw}/natE.npy')
    acc = {}
    for arm in ('D11FWT', 'D11FUWT'):
        s = np.load(f'{raw}/scores_{arm}.npy'); w, best = tdc(pop.assign(row=np.arange(len(pop))), s)
        acc[arm] = best[best.q <= 0.01].set_index('precursor_id')[['row', 'is_decoy']]
    a, b = acc['D11FWT'], acc['D11FUWT']
    for lab, ids, src_ in (('gained', b.index.difference(a.index), b), ('lost', a.index.difference(b.index), a)):
        d = src_.loc[ids].copy(); d['grp'] = np.where(d.is_decoy == 1, 'decoy', np.where(d.index.map(ent) == 1, 'entrapment', 'human'))
        d['U_apex'] = XU[d.row, 0]; d['E_apex'] = XE[d.row, 0]; d['n_unique'] = XU[d.row, 5]
        print(raw, lab, d.groupby('grp').agg(n=('row', 'size'), U_apex_med=('U_apex', 'median'), E_apex_med=('E_apex', 'median'), n_unique_med=('n_unique', 'median')).round(2).to_dict('index'))
    b2 = b.copy(); b2['grp'] = np.where(b2.is_decoy == 1, 'decoy', np.where(b2.index.map(ent) == 1, 'entrapment', 'human'))
    print(raw, 'D11FUWT accepted at q<=0.01 by group:', b2.grp.value_counts().to_dict())
