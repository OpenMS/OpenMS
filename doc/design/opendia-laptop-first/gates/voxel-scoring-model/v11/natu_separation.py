"""PLAN_V11 pre-plan look (raw 01 only; descriptive; labels only to measure). For natU (unique fragments) vs natE (all
fragments): (1) calibration on decoy windows, (2) score-matched AUC D vs DEC, ENT vs DEC, all-human vs DEC (D11FWT strata,
as NATE_SEPARATION_raw01.txt), (3) the key quantity for this plan: ENT vs DEC must move toward 0.5 while D vs DEC stays high."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, 'v4'); from dd_band import auc, ent
N = ['apex_logp', 'sum_logp', 'sig_scans', 'persist_logp', 'cos']
pop = pd.read_csv('raw01/population.tsv', sep='\t'); XE = np.load('raw01/natE.npy'); XU = np.load('raw01/natU.npy')
dec = pop.is_decoy.values == 1
for lab, A in (('E (all fragments)', XE), ('U (unique fragments)', XU)):
    print(f'{lab}: decoy windows P(apex_logp >= 2) {np.mean(A[dec, 0] >= 2):.3f}, >= 3 {np.mean(A[dec, 0] >= 3):.3f}; target {np.mean(A[~dec, 0] >= 2):.3f}, {np.mean(A[~dec, 0] >= 3):.3f}')
s = np.load('raw01/scores_D11FWT.npy').astype(float)
w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id'); w = w[w.s > -1e20].copy()
st = pd.read_csv('v3seq/misses_D11FWT_raw01.tsv', sep='\t').set_index('precursor_id').stage
w['grp'] = np.where(w.precursor_id.map(st) == 'D right window, below cut', 'D', np.where(w.is_decoy == 1, 'DEC', np.where(w.precursor_id.map(ent) == 1, 'ENT', 'HUM')))
for i, k in enumerate(N): w['E_' + k] = XE[w.row.values, i]; w['U_' + k] = XU[w.row.values, i]
w['U_n_unique'] = XU[w.row.values, 5]
top = w[w.s >= w[w.grp == 'DEC'].s.quantile(0.9)].copy()
print('median n_unique among winners above decoy q90:', top.groupby('grp').U_n_unique.median().to_dict())
for ref, lab in (('D', 'D'), ('HUM', 'all human')):
    g = top[top.grp.isin([ref, 'DEC', 'ENT'])].copy()
    edges = np.unique(np.r_[-np.inf, g[g.grp == ref].s.quantile(np.linspace(0.1, 0.9, 9)).values, np.inf]); g['st'] = pd.cut(g.s, edges, labels=False)
    print(f'== {lab} vs DEC | ENT vs DEC, score-matched')
    for k in [p + k for k in N for p in ('E_', 'U_')] + ['U_n_unique']:
        a1, a2, wt = [], [], []
        for _, t in g.groupby('st'):
            x, y, z = t[t.grp == ref][k], t[t.grp == 'DEC'][k], t[t.grp == 'ENT'][k]
            if len(x) >= 5 and len(y) >= 5 and len(z) >= 5: a1.append(auc(x, y)); a2.append(auc(z, y)); wt.append(len(x))
        print(f'  {k:16s} {np.average(a1, weights=wt):.3f} | ENT {np.average(a2, weights=wt):.3f}')
