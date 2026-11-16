"""PLAN_MS1 criterion 4 (evaluation; labels only here), seed 0.
Groups are fixed on the CUR run, as in PRECHECK_MS1: R = DIA-NN-supported (human, q <= 0.01, apex within 6 s) target
winners below CUR's N@FDP_c 1 % cut; E = entrapment winners with CUR score in [median CUR score of R, cut); ED = the 3
decoy winners nearest each E in CUR score. The AUC(E > ED) is then computed under the CUR score (about 0.5 by
construction) and under the NEW score of the same precursors (the two arms share search and population; only scoring
inputs differ). Criterion: AUC_NEW <= AUC_CUR + 0.02. Matching on the score under test would make the AUC 0.5 by
construction, so the groups must be fixed by one arm; CUR is the reference. Supplementary: AUC(entrapment winners >
decoy winners) over all winners below each run's own cut.
Usage: crit4.py <cur precursors.tsv> <new precursors.tsv> <raw01|raw02>"""
import sqlite3, sys, numpy as np, pandas as pd
CUR, NEW, RAW = sys.argv[1:4]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; G = B + '/gap-20261114'
sys.argv = ['x']
code = open(B + '/pairlib-20261113/scripts/evaluate_v2.py').read().split('rng = np.random.default_rng(5)')[0].replace("B + '/pairlib-20261113/lib/library.pqp", f"'{G}/armE_raw01/library.pqp")
exec(code)
from sklearn.metrics import roc_auc_score
c2 = sqlite3.connect(f'file:{G}/armE_raw01/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("select p.ID precursor_id,p.CHARGE z,p.DECOY decoy,pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID", c2)
pre['key'] = pre.ms + '/' + pre.z.astype(str); tg = pre[pre.decoy == 0].set_index('key')
REF = {'raw01': B + '/diann-pairlib-20261113/report_raw01.parquet', 'raw02': B + '/diann-decoys-20261114/own/raw02/report.parquet'}[RAW]
ref = pd.read_parquet(REF, columns=['Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value'])
ref['key'] = ref['Modified.Sequence'] + '/' + ref['Precursor.Charge'].astype(str); ref = ref[ref['Q.Value'] <= 0.01].drop_duplicates('key'); ref = ref[ref.key.isin(tg.index)]
ref['pid'] = tg.loc[ref.key, 'precursor_id'].values; ref = ref[ref.pid.map(ent) == 0].set_index('pid')
def load(p): return pd.read_csv(p, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score']).set_index('precursor_id')
def cut_of(t):
    T = t[t.is_decoy == 0].sort_values('score', ascending=False); e = T.index.map(ent).values
    fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(T) + 1); n1 = int(np.where(fdp <= 0.01)[0].max() + 1)
    return n1, T.score.iloc[n1 - 1], set(T.index[:n1])
cur, new = load(CUR), load(NEW)
n1, cut, acc = cut_of(cur)
Rj = [p for p in ref.index if p in cur.index and p not in acc and abs(cur.rt_apex[p] - ref.RT[p] * 60) <= 6]
lo = np.median(cur.score.loc[Rj])
T = cur[cur.is_decoy == 0]
E = T[(T.index.map(ent) == 1) & (T.score < cut) & (T.score >= lo)].index
D = cur[cur.is_decoy == 1].sort_values('score'); ds = D.score.values
ix = np.clip(np.searchsorted(ds, cur.score.loc[E].values), 1, len(ds) - 2)
ED = D.index[np.unique(np.r_[ix - 1, ix, ix + 1])]
common_E = [p for p in E if p in new.index]; common_D = [p for p in ED if p in new.index]
y = np.r_[np.ones(len(common_E)), np.zeros(len(common_D))]
a_cur = roc_auc_score(y, np.r_[cur.score.loc[common_E].values, cur.score.loc[common_D].values])
a_new = roc_auc_score(y, np.r_[new.score.loc[common_E].values, new.score.loc[common_D].values])
print(f'{RAW}: CUR N@1% {n1}; R {len(Rj)}; E {len(E)} ({len(common_E)} scored in NEW), ED {len(ED)} ({len(common_D)})')
print(f'{RAW}: AUC(E > ED) under CUR {a_cur:.4f}, under NEW {a_new:.4f}; criterion 4 (NEW <= CUR + 0.02): {"pass" if a_new <= a_cur + 0.02 else "FAIL"}')
for name, t in (('CUR', cur), ('NEW', new)):
    _, c_, _ = cut_of(t); b = t[t.score < c_]
    ee = b.index.map(ent).fillna(0).values == 1; m = (b.is_decoy.values == 1) | ee
    print(f'{RAW} {name}: supplementary AUC(entrapment winners > decoy winners, below own cut) {roc_auc_score(ee[m], b.score.values[m]):.4f}')
