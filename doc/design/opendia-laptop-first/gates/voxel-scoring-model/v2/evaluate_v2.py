"""Voxel v2 step 4: PLAN_V2.md reading. Only this script reads entrapment labels.
Per raw and arm: pooled paired TDC on the arm's logit over the v2 population (the tool's -pooled_tdc), precursor
target winners ranked by the standardised winner score, N at entrapment FDP_c 0.5/1/1.5/2 %, q <= 0.01 / 0.05 N and
FDP_c, entrapment-over-decoy winner AUC, and the paired bootstrap 90 % interval of N@1 % and N@2 % against BASE.
BASE is the frozen CNN + LDA stage2_score restricted to the same population (same metric, same TDC)."""
import sqlite3, sys, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110/v2'
R = 148983 / 151239
BASEF = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion/features.tsv',
         'raw02': B + '/astral-raw02-byion-design-20261109/byion/features.tsv'}
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
ARMS = sys.argv[1:] or ['V0', 'F0', 'S0', 'FS']

def tdc(pop, s):
    w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
    T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
    best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    return w, best

def n_at(s, e, t):
    o = np.argsort(-s, kind='stable'); e = e[o]
    f = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(f <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0

def read(raw, arm, pop):
    if arm == 'BASE':
        b = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'stage2_score']).set_index('sample_id').stage2_score
        return b.loc[pop.sample_id].values
    return np.load(f'{O}/{raw}/scores_{arm}.npy')

rng = np.random.default_rng(5)
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
    res = {}
    for arm in ['BASE'] + ARMS:
        try: s = read(raw, arm, pop)
        except FileNotFoundError: print(f'{raw} {arm}: missing'); continue
        w, best = tdc(pop, s)
        tw = w[w.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
        dw = w[w.is_decoy == 1]
        curve = [n_at(tw.s.values, tw.e.values, t) for t in (0.005, 0.01, 0.015, 0.02)]
        tb = best[best.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values)
        qn = {q: (int((tb.q <= q).sum()), tb[tb.q <= q].e.sum() * (1 + 1 / R) / max((tb.q <= q).sum(), 1)) for q in (0.01, 0.05)}
        ea = tw[tw.e == 1].s.values; auc = mannwhitneyu(ea, dw.s.values).statistic / (len(ea) * len(dw))
        res[arm] = tw.set_index('precursor_id')[['s', 'e']]
        line = (f'{raw} {arm:4s} N@FDP_c 0.5/1/1.5/2%: {curve} | q<=0.01 {qn[0.01][0]} FDP_c {qn[0.01][1]:.3f} | '
                f'q<=0.05 {qn[0.05][0]} FDP_c {qn[0.05][1]:.3f} | AUC(entrap>decoy) {auc:.4f}')
        if arm != 'BASE' and 'BASE' in res:
            m = res[arm].join(res['BASE'].s.rename('b'), how='inner')
            a_, b_, e_ = m.s.values, m.b.values, m.e.values.astype(int); d1, d2 = [], []
            for _ in range(1000):
                i = rng.integers(0, len(m), len(m))
                d1.append(n_at(a_[i], e_[i], 0.01) - n_at(b_[i], e_[i], 0.01)); d2.append(n_at(a_[i], e_[i], 0.02) - n_at(b_[i], e_[i], 0.02))
            line += f' | vs BASE 90%: @1% {np.percentile(d1, 5):+.0f}..{np.percentile(d1, 95):+.0f} @2% {np.percentile(d2, 5):+.0f}..{np.percentile(d2, 95):+.0f}'
        print(line, flush=True)
