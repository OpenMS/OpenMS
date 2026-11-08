"""Follow-ups to oswlike_score (raw 01, design diagnostic; labels target/decoy only):
 A. 'minus xcorr' collapsed to 0: rerun it with the SS loop initialised on lib_dotprod (initialisation artefact check).
 B. Candidate generation: same OSW-like subscores and the same SS scorer, but the candidate apexes are OUR streaming
    search's (search.features.tsv rt_apex, up to K=5 per precursor by apex intensity, snapped to the nearest XIC scan).
    If B ~ the XIC-maxima result, our candidate generator is not the gap; the gap is features + scorer.
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, '.')
import oswlike_feats as F
import oswlike_score as S

print('A. minus xcorr, SS initialised on lib_dotprod')
S.report('  LDA minus xcorr (init lib_dotprod)', S.crossfit([c for c in S.ALL if c not in S.GROUPS['xcorr']], main='lib_dotprod'))

# B. our search's candidates as apexes
s = pd.read_csv(F.D + '/pipeline/search.features.tsv', sep='\t', usecols=['precursor_id', 'rt_apex', 'apex_intensity'])
s = s.sort_values(['precursor_id', 'apex_intensity'], ascending=[True, False]).groupby('precursor_id').head(F.K)
row = pd.Series(np.arange(F.N), index=F.meta.precursor_id.values)
s = s[s.precursor_id.isin(row.index)]
ri = row.reindex(s.precursor_id).values
rtm = F.rt[ri]; n = F.npts[ri]
idx = np.nanargmin(np.abs(rtm - s.rt_apex.values[:, None]), axis=1)
s = s.assign(ri=ri, scan=idx, rank=s.groupby('precursor_id').cumcount())
given = np.full((F.N, F.K), -1, np.int64)
given[s.ri.values, s['rank'].values] = s.scan.values
okm = F.ok & (given[:, 0] >= 0)
res, ap = F.run(F.frag, F.ms1, F.lib, F.npts, F.rt, F.meta.cal_rt.fillna(-1e9).values, okm, given)
ii, kk = np.where(np.isfinite(res[:, :, 0]))
X = res[ii, kk].astype(np.float64); apex = ap[ii, kk]
m = F.meta.iloc[ii].reset_index(drop=True)
right = (m.decoy.values == 0) & m.dn_rt.notna().values & (np.abs(apex - m.dn_rt.values) <= 10)
cov = pd.Series(right).groupby(m.precursor_id.values).any().sum()
print(f'B. our search candidates: {len(X)} | confirmed precursors with the right peak among them {cov} of {int(F.meta.dn_rt.notna().sum())}', flush=True)
# swap the scorer's globals to this candidate set
S.X0 = X; S.prec = m.precursor_id.values; S.dec = m.decoy.values.astype(bool); S.fold = m.fold.values; S.pair = m.pair.values
S.entr = m.ent.fillna(0).values.astype(bool); S.right = right
S.report('  LDA all 20 on our search candidates', S.crossfit(S.ALL))
S.report('  MLP all 20 on our search candidates', S.crossfit(S.ALL, 'mlp'))
