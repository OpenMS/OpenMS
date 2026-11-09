"""Fragment commonness (follow-up to c2st.py, declared before running). Hypothesis: the top-6 library fragments of
entrapment precursors are more 'common' (shared with many real tryptic peptides in the same isolation window, e.g.
low-mass y ions) than those of decoys, so their XICs pick up more real background. That would explain wider peaks,
higher S/N and library correlation in entrapment winners (c2st_marg02.txt), and it is computable without entrapment
labels.
commonness(precursor) = mean over its top-6 fragments of log1p(number of human-target top-6 fragments of OTHER
precursors in the same 2 Th window within +-17 ppm). Also reported: the share of top-6 fragments below 300 Th, and the
ion-type/ordinal mix from the ANNOTATION column.
Reading: the hypothesis is supported on a run if entrapment and decoy winners differ in mean commonness (AUC >= 0.53),
AND reweighting decoy winners to the entrapment commonness distribution (10 quantile bins) removes at least half of the
5% excess (MLP seeds 0/1/2, >=2 of 3), on both runs. Family O, full set."""
import os, sys, sqlite3, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
c = sqlite3.connect('file:merged.pqp?mode=ro', uri=True)
tr = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ f, t.LIBRARY_INTENSITY li, t.ID t, t.ANNOTATION ann from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
tr = tr.sort_values(['pid', 'li', 't'], ascending=[True, False, True]).groupby('pid').head(6)
pr = pd.read_sql('select ID pid, PRECURSOR_MZ mz, DECOY d, TRAML_ID tid from PRECURSOR', c)
acc = pd.read_sql('select pp.PRECURSOR_ID pid, pr.PROTEIN_ACCESSION a from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID', c)
ent = acc.assign(e=acc.a.str.contains('_p_target')).groupby('pid').e.all()
pr['ent'] = pr.pid.map(ent).fillna(False).astype(bool)
pr['win'] = np.floor((pr.mz - 380) / 2).astype(int)
tr = tr.merge(pr[['pid', 'win', 'd', 'ent']], on='pid')
hum = tr[(tr.d == 0) & (~tr.ent)]
ref = {w: (np.sort(g.f.values), g.pid.values[np.argsort(g.f.values)]) for w, g in hum.groupby('win')}
def count(fs, pid, w):
    a, ids = ref.get(w, (np.array([]), np.array([])))
    out = []
    for f in fs:
        lo, hi = np.searchsorted(a, f * (1 - 17e-6)), np.searchsorted(a, f * (1 + 17e-6))
        out.append(int((ids[lo:hi] != pid).sum()))
    return np.log1p(np.array(out)).mean()
grp = tr.groupby('pid')
COM = {pid: count(g.f.values, pid, g.win.values[0]) for pid, g in grp}
LOW = grp.f.apply(lambda f: (f < 300).mean()).to_dict()
tr['ion'] = tr.ann.str.extract(r'^([abcxyz])(\d+)')[0].fillna('?') + tr.ann.str.extract(r'^[abcxyz](\d+)')[0].fillna('').map(lambda s: s if s in ('1', '2', '3') else ('4+' if s else ''))
cls = np.where(tr.d == 1, 'decoy', np.where(tr.ent, 'entrapment', 'human'))
print('ion mix of top-6 (share):'); print(pd.crosstab(cls, tr.ion, normalize='index').round(3).to_string())
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
def band(ke, ne, sd):
    return f'{ke / ne / sd:.2f} [{max(0, ke - 1.645 * np.sqrt(ke)) / ne / sd:.2f},{(ke + 1.645 * np.sqrt(ke + 1)) / ne / sd:.2f}]'
for seed in (0, 1, 2):
    sc = S.crossfit(S.ALL, kind='mlp', seed=seed * 10)
    w = S.winners(np.arange(len(sc)), sc); kept, q = S.tdc_q(w, sc); t = sc[kept[(~S.dec[kept]) & (q <= 0.05)]].min()
    dw = w[S.dec[w]]; ew = w[(~S.dec[w]) & S.entr[w]]
    ce = np.array([COM.get(p, 0) for p in S.prec[ew]]); cd = np.array([COM.get(p, 0) for p in S.prec[dw]])
    le = np.array([LOW.get(p, 0) for p in S.prec[ew]]); ld = np.array([LOW.get(p, 0) for p in S.prec[dw]])
    auc = roc_auc_score(np.r_[np.ones(len(ce)), np.zeros(len(cd))], np.r_[ce, cd])
    edges = np.unique(np.quantile(np.r_[ce, cd], np.linspace(0, 1, 11))); be = np.clip(np.searchsorted(edges, ce, 'right') - 1, 0, len(edges) - 2); bd = np.clip(np.searchsorted(edges, cd, 'right') - 1, 0, len(edges) - 2)
    pe = np.bincount(be, minlength=len(edges) - 1) / len(be); pdd = np.bincount(bd, minlength=len(edges) - 1) / len(bd)
    wts = (pe / np.maximum(pdd, 1e-12))[bd]
    se = sc[ew] >= t; sd = sc[dw] >= t; r = se.mean() / sd.mean(); sdw = (wts * sd).sum() / wts.sum(); rw = se.mean() / sdw
    print(f'seed{seed}: commonness mean e {ce.mean():.3f} d {cd.mean():.3f} AUC {auc:.3f} | frac<300Th e {le.mean():.3f} d {ld.mean():.3f}'
          f' | 5% S_e/S_d {band(int(se.sum()), len(ew), sd.mean())} -> reweighted {band(int(se.sum()), len(ew), sdw)} | excess removed {(r - rw) / max(1e-9, r - 1):.2f}', flush=True)
