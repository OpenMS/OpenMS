"""PRECHECK_DIANNVIEW (diagnostic). How DIA-NN itself rates the right-peak IDs we reject (R) vs those we accept (A):
DIA-NN Q.Value / PEP / Evidence / Empirical.Quality / FWHM / Precursor.Quantity rank, and how many R lie inside
DIA-NN's own N@FDP_c 1 % set (its true-1 % ranking cut) versus only in its looser q <= 0.01 tail."""
import sys, numpy as np, pandas as pd
exec(open('frag_precheck.py').read().split("for RAW, REF in")[0])
for RAW, REF in (('raw01', B + '/diann-pairlib-20261113/report_raw01.parquet'), ('raw02', B + '/diann-decoys-20261114/own/raw02/report.parquet')):
    D = f'{G}/armE_{RAW}'; pop = pd.read_csv(f'{D}/population.tsv', sep='\t'); pop['row'] = np.arange(len(pop))
    sc = np.load(f'{G}/cppscore/{RAW}_s0.npy').astype(float); w, best = tdc(pop, sc)
    tw = w[w.is_decoy == 0].sort_values('s', ascending=False); e = tw.precursor_id.map(ent).values
    fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(tw) + 1); n1 = int(np.where(fdp <= 0.01)[0].max() + 1)
    acc = set(tw.precursor_id.iloc[:n1]); W = w.set_index('precursor_id')
    ref = pd.read_parquet(REF); ref['key'] = ref['Modified.Sequence'] + '/' + ref['Precursor.Charge'].astype(str)
    ref = ref.sort_values('Q.Value').drop_duplicates('key'); ref = ref[ref.key.isin(tg.index)]
    ref['pid'] = tg.loc[ref.key, 'precursor_id'].values; ref['ent'] = ref.pid.map(ent).values
    # DIA-NN true-1 % cut by its own ranking (PEP then Q)
    ref = ref.sort_values(['Q.Value', 'PEP']); ef = np.cumsum(ref.ent.values) * (1 + 1 / R) / np.arange(1, len(ref) + 1)
    dn1 = int(np.where(ef <= 0.01)[0].max() + 1); ref['in_true1'] = np.arange(len(ref)) < dn1
    ref = ref[(ref['Q.Value'] <= 0.01) & (ref.ent == 0)].set_index('pid')
    ref['ours_rt'] = [pop.rt_apex.values[int(W.row[p])] if p in W.index else np.nan for p in ref.index]
    ref['grp'] = np.where(~ref.index.isin(W.index), 'unscored', np.where(np.abs(ref.ours_rt - ref.RT * 60) <= 6, np.where(ref.index.isin(acc), 'A', 'R'), np.where(ref.index.isin(acc), 'A_wrongRT', 'Wrong')))
    ref['lq'] = np.log10(ref['Precursor.Quantity'].clip(lower=1)); ref['qrank'] = ref.lq.rank(pct=True)
    print(f'== {RAW}: DIA-NN true-1 % cut at {dn1} of {len(ref)} non-entrapment q<=0.01')
    g = ref.groupby('grp')
    print(g.agg(n=('RT', 'size'), in_true1=('in_true1', 'mean'), q_med=('Q.Value', 'median'), pep_med=('PEP', 'median'), ev=('Evidence', 'median'),
                emp=('Empirical.Quality', 'median'), fwhm=('FWHM', 'median'), qty_pct=('qrank', 'median'), ms1corr=('Ms1.Profile.Corr', 'median')).round(4).to_string())
    r = ref[ref.grp == 'R']
    print('  R by DIA-NN Q.Value bin:', pd.cut(r['Q.Value'], [-1, 1e-4, 1e-3, 3e-3, 6e-3, 0.0101]).value_counts().sort_index().to_dict())
