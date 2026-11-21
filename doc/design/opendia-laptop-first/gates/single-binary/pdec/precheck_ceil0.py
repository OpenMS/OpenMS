"""PRECHECK_CEIL0 (diagnostic, 21 Nov 2026): is the r0 gap to DIA-NN a scoring problem or a candidate/feature problem?
Input: one -features_out run of the current default (opendia-pd10-rtrdef.sif) per run and level, seed 0.
DIA-NN labels and entrapment labels are used only for training the ceiling and for evaluation; nothing here is a method.

Arms (same method as precheck_scorer36.py: best window per precursor, paired target-decoy competition, 3-fold
cross-fitting by the library fold, HistGBM 300 rounds / 31 leaves, balanced class weights):
  A    tool score as written
  B    label-free GBM, 3 Percolator rounds from the tool score (positives = own q<=0.01 target winners)
  Bt   B, but negatives = windows of the hardest 2 % of decoy precursors by the current score (tail-focused)
  C    ceiling: positives = target windows within 6 s of DIA-NN RT for DIA-NN HeLa q<=0.01 IDs of the same library level
  Ct   C against hard decoys only
Reported: N at own q<=0.01 (all levels), N at combined FDP 1 % and bounds (r7 only), share of DIA-NN IDs.
Funnel of DIA-NN-only HeLa IDs (in DIA-NN q<=0.01, not in tool q<=0.01):
  A_low  a window within 6 s of DIA-NN RT exists and is the precursor's best window, but scores too low
  B2     a window within 6 s exists but another window of the precursor scores higher
  C      no window within 6 s
plus TPR at decoy FPR 1e-4 / 1e-3 / 1e-2 of the best right window, tool and ceiling, by DIA-NN quantity quartile.

Usage: precheck_ceil0.py <features.tsv> <lib.pqp> <diann report.parquet> <tag> [r]   (r = entrapment ratio; 0 for r0)
"""
import os, sys, sqlite3, time, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
SRC, L, DN, TAG = sys.argv[1:5]; R = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
OUT = os.path.dirname(os.path.abspath(SRC)); PQ = f'{OUT}/{TAG}_feat.parquet'; T0 = time.time()
if not os.path.exists(PQ):
    F = pd.read_csv(SRC, sep='\t', engine='c', low_memory=False, dtype={'candidate_row': np.int32, 'precursor_id': np.int64,
                    'is_decoy': np.int8, 'fold': np.int8, 'rt_apex': np.float64, 'q_value': np.float64})
    for c in F.columns:
        if F[c].dtype == np.float64 and c not in ('rt_apex', 'q_value', 'score'): F[c] = F[c].astype(np.float32)
    F.to_parquet(PQ)
F = pd.read_parquet(PQ)
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value']
X = [c for c in F.columns if c not in meta]
print(f'{TAG}: windows {len(F)}, precursors {F.precursor_id.nunique()}, features {len(X)}  [{time.time() - T0:.0f}s]', flush=True)
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.TRAML_ID t, p.DECOY d, p.CHARGE z, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
tg = pq[pq.d == 0]; k = 1 + 1 / R if R > 0 else np.nan
pair_id = pd.Series(pd.factorize(pq.key)[0], index=pq.pid); ent = pq.set_index('pid').e
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]
tgs = tg.assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
d['pid'] = tgs.pid.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64)
d = d[ent.reindex(d.pid).values == 0].sort_values('Q.Value').drop_duplicates('pid').set_index('pid')
print(f'DIA-NN HeLa q<=0.01 IDs mapped: {len(d)}', flush=True)
pid = F.precursor_id.values; dec = F.is_decoy.values.astype(bool); fold = F.fold.values; pair = pair_id.reindex(pid).values

def best_rows(score):
    o = np.lexsort((-score, pid)); first = np.r_[True, pid[o][1:] != pid[o][:-1]]; rows = o[first]
    b = pd.DataFrame({'row': rows, 'pid': pid[rows], 's': score[rows], 'dec': dec[rows], 'pair': pair[rows]})
    b = b.sort_values(['pair', 's', 'dec'], ascending=[True, False, False]).drop_duplicates('pair')
    b = b.sort_values(['s', 'dec'], ascending=[False, False])
    nd = np.cumsum(b.dec.values); nt = np.cumsum(~b.dec.values)
    b['q'] = np.minimum.accumulate((nd / np.maximum(nt, 1))[::-1])[::-1]
    return b

RES = {}
def evaluate(score, name):
    b = best_rows(score); t = b[~b.dec.values].copy(); t['e'] = ent.reindex(t.pid).values
    acc = t[t.q <= 0.01]; N, E = len(acc), int(acc.e.sum()); hela = acc[acc.e == 0]
    s = f'{name:40s} q<=0.01 N {N:6d} HeLa {len(hela):6d} (in DIA-NN {hela.pid.isin(d.index).mean():.3f}, DIA-NN covered {d.index.isin(hela.pid).mean():.3f})'
    n1 = np.nan
    if R > 0:
        fdp = np.cumsum(t.e.values) * k / np.arange(1, len(t) + 1); ok = np.where(fdp <= 0.01)[0]; n1 = int(ok.max() + 1) if len(ok) else 0
        s += f' lower {E / max(N, 1):.4f} combined {E * k / max(N, 1):.4f} N@FDP_c1% {n1}'
    print(s + f'  [{time.time() - T0:.0f}s]', flush=True)
    RES[name] = dict(hela=len(hela), n1=n1, ids=set(hela.pid)); return len(hela)

Xv = F[X].values.astype(np.float32)
def gbm(pos_mask, neg_mask, seed=0, max_neg=1_500_000):
    s = np.zeros(len(F)); rng = np.random.default_rng(seed)
    for f in range(3):
        tp = np.where(pos_mask & (fold != f))[0]; tn = np.where(neg_mask & (fold != f))[0]
        if len(tn) > max_neg: tn = rng.choice(tn, max_neg, replace=False)
        idx = np.r_[tp, tn]; y = np.r_[np.ones(len(tp)), np.zeros(len(tn))]
        w = np.r_[np.full(len(tp), len(tn) / max(len(tp), 1)), np.ones(len(tn))]
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(Xv[idx], y, sample_weight=w)
        te = np.where(fold == f)[0]; s[te] = m.decision_function(Xv[te])
    return s

def winners(score):
    b = best_rows(score); m = np.zeros(len(F), bool); m[b.row.values[(b.q <= 0.01) & ~b.dec.values]] = True; return m
def hard_decoys(score, frac=0.02):
    b = best_rows(score); bd = b[b.dec.values]; top = set(bd.pid.values[:max(1, int(frac * len(bd)))])
    return dec & np.isin(pid, list(top))

kept = F.score.values > -1e29
s_tool = np.where(kept, F.score.values, -1e30).astype(np.float64)
a = evaluate(s_tool, 'A tool score')
s = s_tool
for it in range(3):
    s = gbm(winners(s), dec, seed=it)
bn = evaluate(s, 'B label-free GBM (3 rounds)'); sB = s
s = s_tool
for it in range(3):
    s = gbm(winners(s), hard_decoys(s), seed=20 + it)
btn = evaluate(s, 'Bt label-free GBM, hard decoys (3 rounds)')
drt = np.abs(F.rt_apex.values - (d.RT * 60).reindex(pid).values)
right = (~dec) & (drt <= 6)
sC = gbm(right, dec, seed=9); cn = evaluate(sC, f'C DIA-NN-labelled ceiling ({right.sum()} pos)')
sCt = gbm(right, hard_decoys(s_tool), seed=19); ctn = evaluate(sCt, 'Ct DIA-NN-labelled, hard decoys')

# Funnel of DIA-NN-only IDs under the tool score
tool_ids = RES['A tool score']['ids']; miss = d.index[~d.index.isin(list(tool_ids))]
W = pd.DataFrame({'pid': pid, 's': s_tool, 'right': right, 'dec': dec})
W = W[~W.dec]
best = W.sort_values('s', ascending=False).drop_duplicates('pid').set_index('pid')
has_right = set(W.pid[W.right])
cls = pd.Series('C', index=miss)
cls[[p in has_right for p in miss]] = 'B2'
br = best.right.reindex(miss).fillna(False).values.astype(bool)
cls[br] = 'A_low'
not_scored = ~pd.Index(miss).isin(best.index)
print(f'\nfunnel of DIA-NN-only HeLa IDs ({len(miss)} of {len(d)}; {not_scored.sum()} have no window at all):',
      ' '.join(f'{c} {(cls == c).mean():.3f}' for c in ('A_low', 'B2', 'C')), flush=True)
q = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
print('  by DIA-NN quantity quartile (share of that quartile missed; A_low / B2 / C):')
for qq in ['Q1', 'Q2', 'Q3', 'Q4']:
    m_all = q.index[q == qq]; mm = cls.reindex(m_all).dropna()
    print(f'   {qq}: missed {len(mm) / len(m_all):.3f}  ' + ' / '.join(f'{(mm == c).mean():.3f}' for c in ('A_low', 'B2', 'C')))

# TPR of the best right window at fixed decoy FPR, tool vs ceiling
bd = best_rows(s_tool); Dd = bd[bd.dec.values]
for name, sc in (('tool', s_tool), ('B', sB), ('C', sC)):
    dsc = np.sort(sc[Dd.row.values])[::-1]
    Rw = pd.DataFrame({'pid': pid[right], 's': sc[right]}).sort_values('s', ascending=False).drop_duplicates('pid').set_index('pid')
    line = []
    for fpr in (1e-4, 1e-3, 1e-2):
        thr = dsc[max(0, int(fpr * len(dsc)) - 1)]
        line.append(f'FPR {fpr:g}: ' + ' '.join(f'{qq} {(Rw.s.reindex(q.index[q == qq]) > thr).mean():.3f}' for qq in ['Q1', 'Q2', 'Q3', 'Q4'])
                    + f' all {(Rw.s.reindex(d.index) > thr).mean():.3f}')
    print(f'TPR {name:4s} ' + ' | '.join(line), flush=True)
print(f'\nSUMMARY {TAG}: DIA-NN {len(d)}; HeLa at own q<=0.01: tool {a}, B {bn} ({bn / a - 1:+.1%}), Bt {btn} ({btn / a - 1:+.1%}), '
      f'C {cn} ({cn / a - 1:+.1%}), Ct {ctn} ({ctn / a - 1:+.1%}); tool/DIA-NN {a / len(d):.3f}, best ceiling/DIA-NN {max(cn, ctn) / len(d):.3f}'
      + (f"; N@FDP_c 1 %: tool {RES['A tool score']['n1']}, C {RES[next(n for n in RES if n.startswith('C '))]['n1']}" if R > 0 else ''))
