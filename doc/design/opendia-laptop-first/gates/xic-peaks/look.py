"""XIC look read-outs R3-R7 (PLAN_XIC_LOOK.md). Labels (entrapment, DIA-NN) are used only here, for evaluation."""
import sqlite3, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/xic-peaks-20261112'
V = B + '/astral-voxel-pilot-20261110/v2'
R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
pre = pd.read_sql('select ID precursor_id, TRAML_ID traml, DECOY decoy from PRECURSOR', c)
pre['pair_key'] = pre.traml.str.replace('^DECOY_', '', regex=True)
pk = pre.set_index('precursor_id').pair_key

def load(path, cols):
    return pd.read_csv(path, sep='\t', usecols=cols)

xic = load(X + '/raw01/xic/search.features.tsv', ['precursor_id', 'is_decoy', 'rt_start', 'rt_apex', 'rt_end', 'transition_support',
                                                 'matched_scans', 'fwhm_seconds', 'coelution', 'peak_rank', 'peaks_in_window', 'apex_intensity'])
leg = load(B + '/rtlin-20261111/raw01/pipeline/search.features.tsv', ['precursor_id', 'is_decoy', 'rt_apex', 'transition_support', 'peak_width_seconds'])
print('R7 rows: xic', len(xic), 'legacy', len(leg), '| precursors with a row: xic', xic.precursor_id.nunique(), 'legacy', leg.precursor_id.nunique())
print('R7 rows per precursor xic', xic.groupby('precursor_id').size().describe().round(2).to_dict())
print('R7 peaks_in_window (rank 0) quantiles', xic[xic.peak_rank == 0].peaks_in_window.quantile([.1, .5, .9]).to_dict())
print('R7 xic peak width (rt_end-rt_start) s quantiles', (xic.rt_end - xic.rt_start).quantile([.1, .5, .9]).round(2).to_dict(),
      '| FWHM s', xic.fwhm_seconds.quantile([.1, .5, .9]).round(2).to_dict())

dn = pd.read_csv(B + '/diann-ourlib-20261112/diann_ourlib_raw01.tsv.gz', sep='\t')
dn = dn[(dn['Q.Value'] <= 0.01)].drop_duplicates('precursor_id').set_index('precursor_id')
dn['rt_s'] = dn.RT * 60; dn['fwhm_s'] = dn.FWHM * 60
dn = dn[dn.index.map(ent).fillna(1).values == 0]  # human targets only
W = pd.read_csv(V + '/v4/stageD_v16/score/reference_winners_raw01.tsv.gz', sep='\t', usecols=['precursor_id', 'full_pct', 'rt_apex', 'group'])
acc = set(W[(W.group == 'H') & (W.full_pct >= 0.999)].precursor_id)  # D11FWT accepted proxy, as check_width.py
print('DIA-NN human q<=0.01:', len(dn), '| also D11FWT-accepted:', len(set(dn.index) & acc))

def near(df, tol, rank0=False):
    d = df[df.precursor_id.isin(dn.index)]
    if rank0: d = d[d.peak_rank == 0]
    d = d.assign(err=(d.rt_apex - dn.rt_s.reindex(d.precursor_id).values).abs())
    best = d.groupby('precursor_id').err.min()
    return {t: round(float((best.reindex(dn.index) <= t).mean()), 4) for t in tol}
print('R3 any-candidate within 3/10 s: legacy', near(leg, [3, 10]), '| xic', near(xic, [3, 10]), '| xic rank0', near(xic, [3, 10], True))
for g, ids in [('DN', set(dn.index) - acc), ('BOTH', set(dn.index) & acc)]:
    sub = dn.loc[sorted(ids)]
    for name, df in [('legacy', leg), ('xic', xic)]:
        d = df[df.precursor_id.isin(sub.index)]; d = d.assign(err=(d.rt_apex - sub.rt_s.reindex(d.precursor_id).values).abs())
        b = d.groupby('precursor_id').err.min().reindex(sub.index)
        print(f'R3 {g} {name}: within 3 s {float((b <= 3).mean()):.3f}, within 10 s {float((b <= 10).mean()):.3f}')

d = xic[xic.precursor_id.isin(dn.index)]
d = d.assign(err=(d.rt_apex - dn.rt_s.reindex(d.precursor_id).values).abs()).sort_values('err').drop_duplicates('precursor_id')
d = d[d.err <= 3]
rho = spearmanr(d.fwhm_seconds, dn.fwhm_s.reindex(d.precursor_id)).statistic
print(f'R4 Spearman xic FWHM vs DIA-NN FWHM (n={len(d)}): {rho:.3f} | medians xic {d.fwhm_seconds.median():.2f} s, DIA-NN {dn.fwhm_s.reindex(d.precursor_id).median():.2f} s')
print('R4 nearest candidate rank among xic peaks (share rank 0/1/2+):', d.peak_rank.value_counts(normalize=True).sort_index().round(3).to_dict())
lb = leg[leg.precursor_id.isin(dn.index)].assign(err=lambda t: (t.rt_apex - dn.rt_s.reindex(t.precursor_id).values).abs()).sort_values('err').drop_duplicates('precursor_id')
lb = lb[lb.err <= 3]
print(f'R4 reference legacy width vs DIA-NN FWHM (n={len(lb)}): {spearmanr(lb.peak_width_seconds, dn.fwhm_s.reindex(lb.precursor_id)).statistic:.3f}')

def auc(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    r = pd.Series(np.concatenate([a, b])).rank().values
    return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))
x0 = xic[xic.peak_rank == 0].set_index('precursor_id')
lbest = leg.groupby('precursor_id').transition_support.max()
groups = {'DN': sorted(set(dn.index) - acc), 'BOTH': sorted(set(dn.index) & acc)}
tgt = x0[x0.is_decoy == 0].index; groups['ENT'] = [p for p in tgt if ent.get(p, 0) == 1]
dec_x = x0[x0.is_decoy == 1]; dec_l = lbest[lbest.index.isin(pre[pre.decoy == 1].precursor_id)]
for g, ids in groups.items():
    xs = x0.coelution.reindex(ids).fillna(0); ls = lbest.reindex(ids).fillna(0)
    print(f'R5 {g} (n={len(ids)}) AUC vs DEC: xic rank0 coelution {auc(xs, dec_x.coelution):.3f}, xic rank0 support {auc(x0.transition_support.reindex(ids).fillna(0), dec_x.transition_support):.3f}'
          f' | legacy best support {auc(ls, dec_l):.3f}')

# R6: single-score target-decoy competition on rank-0 co-elution (ties broken by apex intensity), pairs by pair_key.
s = x0.assign(score=x0.coelution + 1e-12 * np.log1p(x0.apex_intensity), pair=x0.index.map(pk))
best = s.sort_values(['score', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair').sort_values('score', ascending=False)
T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
q = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
ok = (q <= 0.01) & (best.is_decoy.values == 0)
e = best.index[ok].map(lambda p: ent.get(p, 0)).values
n = int(ok.sum()); ne = int(e.sum())
print(f'R6 untrained rank-0 co-elution TDC q<=0.01: N {n}, entrapment {ne}, FDP_c {ne * (1 + 1 / R) / max(n, 1):.4f}')
